import json
from Bio import Entrez
from src.llm.client import LLMClient
from src.sra.search import search_sra, fetch_summary, fetch_details, safe_esearch
from src.sra.search import search_sra, fetch_summary, fetch_details, safe_esearch
from src.utils.xml_parser import parse_sra_xml, parse_sra_full_xml

class Pipeline:
    def __init__(self, ollama_threads=None):
        self.llm = LLMClient(num_threads=ollama_threads)
        self.bioprojects_seen = set()
        self.results = []

    def generate_query(self, topic):
        """
        Core 1: Generate SRA search query using LLM.
        """
        print(f"Generating search query for topic: {topic}...")
        
        keywords = topic
        
        # Expand query with synonyms
        from src.query.expander import QueryExpander
        expander = QueryExpander()
        expanded_concepts = expander.expand_query_string(keywords)
        
        expansion_text = ""
        if expanded_concepts:
            expansion_text = "\nSYNONYM SUGGESTIONS (Use if applicable):\n" + "\n".join(expanded_concepts) + "\n"
        
        prompt = f"""
Eres un asistente experto en bioinformática especializado en minería de datos de NCBI SRA.
Tu objetivo es generar una estrategia de búsqueda de **ALTO RECALL** (alta recuperación) para el tema: "{keywords}".

{expansion_text}
Queremos encontrar la MAYOR cantidad posible de estudios relevantes, incluso si eso implica algo de ruido.
Es preferible traer 10,000 resultados y filtrar después, que traer 0 por ser demasiado específico.

IMPORTANTE: Debes devolver la respuesta en formato JSON VÁLIDO.
Ten mucho cuidado con las comillas dentro de la query. Si usas comillas dobles dentro de la string de la query, DEBES escaparlas con backslash (\").

Reglas para ALTO RECALL:
1. **Usa [All Fields] generosamente**: Si un término no es estrictamente un organismo o estrategia, úsalo en [All Fields].
2. **Evita ANDs innecesarios**: No conectes todos los conceptos con AND. Si el usuario pide "sequía en tomate y pimiento", usa OR para los organismos.
3. **Expande sinónimos**: Usa OR para variantes (ej: "drought" OR "water stress" OR "water deficit").
4. **Organismos**: Identifica el nombre científico pero incluye también el nombre común en [All Fields] (ej: "Solanum lycopersicum"[Organism] OR "tomato"[All Fields]).
5. **Estrategia**: Si el usuario menciona una técnica (ej: RNA-Seq), úsala, pero incluye variantes (ej: "RNA-Seq"[Strategy] OR "transcriptome"[All Fields]).

Ejemplo de salida válida:
{{
  "natural_query": "Búsqueda amplia de RNA-Seq en tomate o pimiento bajo estrés hídrico",
  "esearch_query": "(\\\"Solanum lycopersicum\\\"[Organism] OR \\\"tomato\\\"[All Fields]) AND (\\\"drought\\\"[All Fields] OR \\\"water stress\\\"[All Fields])"
}}

Devuelve SOLO el JSON.
"""
        
        # Disable strict json mode to handle complex escapes manually if needed
        response = self.llm.generate_response(prompt, json_mode=True)
        
        if not response:
            return None

        try:
            # Attempt 1: Direct JSON load
            data = json.loads(response)
        except json.JSONDecodeError:
            # Attempt 2: Try to fix common escape issues or find JSON block
            import re
            try:
                # Find first { and last }
                match = re.search(r'\{.*\}', response, re.DOTALL)
                if match:
                    json_str = match.group(0)
                    data = json.loads(json_str)
                else:
                    raise ValueError("No JSON block found")
            except Exception as e:
                print(f"Error parsing query JSON: {e}")
                print(f"Raw response: {response}")
                return None

        print(f"Natural Query: {data.get('natural_query')}")
        return data.get('esearch_query')

    def enrich_metadata(self, metadata):
        """
        Fetches additional counts (SRA Experiments, BioSamples) for the BioProject.
        """
        bioproject = metadata.get('bioproject')
        if not bioproject:
            return metadata
            
        try:
            # Count SRA Experiments
            # We search SRA for the BioProject accession
            result = safe_esearch(db="sra", term=bioproject, retmax=0)
            metadata['sra_experiment_count'] = int(result['Count'])
            
            # Count BioSamples
            # We search BioSample for the BioProject accession
            result = safe_esearch(db="biosample", term=bioproject, retmax=0)
            metadata['biosample_count'] = int(result['Count'])
            
        except Exception as e:
            print(f"Error fetching counts for {bioproject}: {e}")
            metadata['sra_experiment_count'] = 0
            metadata['biosample_count'] = 0
            
        return metadata

    def analyze_bioproject(self, metadata):
        """
        Core 3: Format and contextualize BioProject using LLM.
        """
        from src.llm.prompts import get_analysis_prompt
        
        # Enrich with counts first
        metadata = self.enrich_metadata(metadata)
        
        # Build text representation for LLM
        # Use full_text if available (from fetch_details), otherwise fall back
        if 'full_text' in metadata:
            text_data = metadata['full_text']
        else:
            text_data = "\n".join([f"{k}: {v}" for k, v in metadata.items()])
        
        system_prompt, user_prompt = get_analysis_prompt(text_data)
        
        response = self.llm.generate_response(user_prompt, system_prompt=system_prompt, json_mode=True)
        
        try:
            return json.loads(response)
        except:
            print(f"Failed to parse LLM JSON response for {metadata.get('bioproject')}")
            return {}

    def run(self, topic, max_workers=15):
        from config.settings import TARGET_PROJECTS, BATCH_SIZE
        import concurrent.futures
        
        # 1. Generate Query
        query = self.generate_query(topic)
        if not query:
            print("Query generation failed.")
            return
        
        print(f"Generated Query: {query}")
        
        # Pagination Loop
        retstart = 0
        total_found = 0
        
        while len(self.bioprojects_seen) < TARGET_PROJECTS:
            print(f"\nFetching batch starting at {retstart} (Target Unique Projects: {TARGET_PROJECTS})...")
            
            # 2. Search SRA (Batch)
            search_result = search_sra(query, retstart=retstart, retmax=BATCH_SIZE)
            
            if not search_result or 'IdList' not in search_result:
                print("Search failed or no more results.")
                break
                
            id_list = search_result['IdList']
            count = int(search_result['Count'])
            total_found = count
            
            if not id_list:
                print("No more IDs returned in this batch.")
                break
                
            print(f"Batch contains {len(id_list)} records. Total available: {count}")

            # 3. Process and Filter in Sub-Batches
            # Fetch details (full XML) in chunks of 50 (efetch is heavier than esummary)
            chunk_size = 50
            new_unique_records = []
            
            for i in range(0, len(id_list), chunk_size):
                chunk_ids = id_list[i:i+chunk_size]
                print(f"   Fetching details for records {i}-{i+len(chunk_ids)}...")
                
                try:
                    xml_content = fetch_details(chunk_ids)
                except Exception as e:
                    print(f"   Error fetching details batch: {e}")
                    continue

                if not xml_content:
                    continue
                
                # Parse XML string to list of dicts
                details_list = parse_sra_full_xml(xml_content)

                for metadata in details_list:
                    if not metadata:
                        continue

                    bioproject = metadata.get('bioproject')
                    
                    # Filter: Unique BioProjects
                    if bioproject and bioproject in self.bioprojects_seen:
                        continue
                    
                    if bioproject:
                        self.bioprojects_seen.add(bioproject)
                    
                    new_unique_records.append(metadata)
                    
                    if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                        break
                
                if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                    break

            print(f"   Found {len(new_unique_records)} new unique BioProjects in this batch.")

            # 4. Parallel LLM Analysis
            if new_unique_records:
                print(f"   Analyzing {len(new_unique_records)} projects in parallel with {max_workers} threads...")
                
                with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                    future_to_meta = {executor.submit(self.analyze_bioproject, meta): meta for meta in new_unique_records}
                    
                    for future in concurrent.futures.as_completed(future_to_meta):
                        meta = future_to_meta[future]
                        try:
                            analysis = future.result()
                            full_record = {**meta, **analysis}
                            self.results.append(full_record)
                            print(f"      Analyzed: {meta.get('bioproject')}")
                        except Exception as exc:
                            print(f"      Error analyzing {meta.get('bioproject')}: {exc}")

            # Check if we processed all available results
            if retstart + len(id_list) >= count:
                print("Reached end of search results.")
                break
                
            # Check if target reached
            if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                print(f"Target of {TARGET_PROJECTS} unique BioProjects reached!")
                break
                
            # Move to next batch
            retstart += len(id_list)
            
        return self.results
