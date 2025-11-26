import json
from src.llm.client import LLMClient
from src.sra.search import search_sra, fetch_summary
from src.utils.xml_parser import parse_sra_xml

class Pipeline:
    def __init__(self):
        self.llm = LLMClient()
        self.bioprojects_seen = set()
        self.results = []

    def generate_query(self, topic):
        """
        Núcleo 1: Generar consulta de búsqueda SRA usando Mistral.
        """
        print(f"🧠 Generando consulta de búsqueda para el tema: {topic}...")
        
        keywords = topic
        
        prompt = f"""
Eres un asistente experto en bioinformática especializado en minería de datos de NCBI SRA.
Tu objetivo es generar una estrategia de búsqueda de **ALTO RECALL** (alta recuperación) para el tema: "{keywords}".

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
        
        # Desactivamos el modo json estricto en la llamada al cliente para evitar posibles problemas con la gramática forzada de Ollama
        # si tiene problemas con el escape complejo. Lo analizaremos manualmente.
        response = self.llm.generate_response(prompt, json_mode=True)
        
        if not response:
            return None

        try:
            # Intento 1: Carga directa de JSON
            data = json.loads(response)
        except json.JSONDecodeError:
            # Intento 2: Intentar arreglar problemas comunes de escape o encontrar el bloque JSON
            import re
            try:
                # Encontrar el primer { y el último }
                match = re.search(r'\{.*\}', response, re.DOTALL)
                if match:
                    json_str = match.group(0)
                    # A veces los modelos olvidan escapar comillas dentro del valor
                    # Este es un problema difícil de arreglar perfectamente con regex, pero podemos intentarlo
                    data = json.loads(json_str)
                else:
                    raise ValueError("No se encontró bloque JSON")
            except Exception as e:
                print(f"❌ Error analizando JSON de consulta: {e}")
                print(f"Respuesta cruda: {response}")
                return None

        print(f"📄 Consulta Natural: {data.get('natural_query')}")
        return data.get('esearch_query')

    def analyze_bioproject(self, metadata):
        """
        Núcleo 3: Formatear y contextualizar BioProject usando Mistral.
        """
        # Construir representación de texto para el LLM
        text_data = "\n".join([f"{k}: {v}" for k, v in metadata.items()])
        
        system_prompt = """Eres un asistente de bioinformática. 
        Analiza los metadatos del estudio SRA proporcionados.
        Devuelve un objeto JSON con las siguientes claves:
        - experimental_conditions: lista de cadenas (tratamientos/condiciones)
        - is_time_series: booleano
        - tissues_studied: lista de cadenas
        - summary: un resumen breve de 1 oración del contexto del estudio
        """
        
        prompt = f"""Analiza este estudio:\n\n{text_data}"""
        
        response = self.llm.generate_response(prompt, system_prompt=system_prompt, json_mode=True)
        
        try:
            return json.loads(response)
        except:
            print(f"⚠️ Falló el análisis de la respuesta JSON del LLM para {metadata.get('bioproject')}")
            return {}

    def run(self, topic):
        from config.settings import TARGET_PROJECTS, BATCH_SIZE
        import concurrent.futures
        
        # 1. Generar Consulta
        query = self.generate_query(topic)
        if not query:
            print("❌ Falló la generación de la consulta.")
            return
        
        print(f"📝 Consulta Generada: {query}")
        
        # Bucle de Paginación
        retstart = 0
        total_found = 0
        
        while len(self.bioprojects_seen) < TARGET_PROJECTS:
            print(f"\n🔄 Obteniendo lote comenzando en {retstart} (Objetivo Proyectos Únicos: {TARGET_PROJECTS})...")
            
            # 2. Buscar en SRA (Lote)
            search_result = search_sra(query, retstart=retstart, retmax=BATCH_SIZE)
            
            if not search_result or 'IdList' not in search_result:
                print("❌ La búsqueda falló o no hay más resultados.")
                break
                
            id_list = search_result['IdList']
            count = int(search_result['Count'])
            total_found = count
            
            if not id_list:
                print("⚠️ No se devolvieron más IDs en este lote.")
                break
                
            print(f"📊 El lote contiene {len(id_list)} registros. Total disponible: {count}")

            # 3. Procesar y Filtrar en Sub-Lotes
            # Obtener resúmenes en trozos de 200 para evitar problemas de longitud de URL o tiempos de espera
            chunk_size = 200
            new_unique_records = []
            
            for i in range(0, len(id_list), chunk_size):
                chunk_ids = id_list[i:i+chunk_size]
                print(f"   📥 Obteniendo resúmenes para registros {i}-{i+len(chunk_ids)}...")
                
                summary_list = fetch_summary(chunk_ids)
                if not summary_list:
                    continue
                
                for study_summary in summary_list:
                    if "ExpXml" not in study_summary:
                        continue

                    # Analizar XML
                    metadata = parse_sra_xml(study_summary["ExpXml"])
                    if not metadata:
                        continue

                    bioproject = metadata.get('bioproject')
                    
                    # Filtro: BioProjects Únicos
                    if bioproject in self.bioprojects_seen:
                        continue
                    
                    self.bioprojects_seen.add(bioproject)
                    new_unique_records.append(metadata)
                    
                    if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                        break
                
                if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                    break

            print(f"   ✨ Encontrados {len(new_unique_records)} nuevos BioProjects únicos en este lote.")

            # 4. Análisis LLM Paralelo
            if new_unique_records:
                print(f"   🤖 Analizando {len(new_unique_records)} proyectos en paralelo...")
                
                with concurrent.futures.ThreadPoolExecutor(max_workers=15) as executor:
                    future_to_meta = {executor.submit(self.analyze_bioproject, meta): meta for meta in new_unique_records}
                    
                    for future in concurrent.futures.as_completed(future_to_meta):
                        meta = future_to_meta[future]
                        try:
                            analysis = future.result()
                            full_record = {**meta, **analysis}
                            self.results.append(full_record)
                            print(f"      ✅ Analizado: {meta.get('bioproject')}")
                        except Exception as exc:
                            print(f"      ❌ Error analizando {meta.get('bioproject')}: {exc}")

            # Verificar si procesamos todos los resultados disponibles
            if retstart + len(id_list) >= count:
                print("🏁 Se alcanzó el final de los resultados de búsqueda.")
                break
                
            # Verificar si alcanzamos el objetivo
            if len(self.bioprojects_seen) >= TARGET_PROJECTS:
                print(f"🎉 ¡Objetivo de {TARGET_PROJECTS} BioProjects únicos alcanzado!")
                break
                
            # Mover al siguiente lote
            retstart += len(id_list)
            
        return self.results
