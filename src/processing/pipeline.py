"""
Main processing pipeline for MAInr.

This module orchestrates the entire workflow: query generation, SRA search,
metadata extraction, and LLM-based analysis.
"""

import json
import re
from typing import List, Dict, Optional, Any
import concurrent.futures

from Bio import Entrez
from tqdm import tqdm

from src.llm.client import LLMClient
from src.sra.search import search_sra, fetch_details, safe_esearch
from src.utils.xml_parser import parse_sra_full_xml
from src.utils.logger import get_logger
from config.settings import TARGET_PROJECTS, BATCH_SIZE
from src.constants import PipelineDefaults

logger = get_logger(__name__)


class Pipeline:
    """
    Main processing pipeline for SRA data mining.

    Handles the complete workflow from query generation to final analysis.
    """

    def __init__(
        self,
        ollama_threads: Optional[int] = None,
        target_projects: Optional[int] = None
    ) -> None:
        """
        Initialize the pipeline.

        Args:
            ollama_threads: Number of threads for Ollama inference (optional)
            target_projects: Target number of unique projects to retrieve (optional)
        """
        self.llm = LLMClient(num_threads=ollama_threads)
        self.bioprojects_seen: set = set()
        self.results: List[Dict[str, Any]] = []

        # Use provided target or default
        self.target_projects = target_projects or TARGET_PROJECTS

        logger.info(f"Pipeline initialized (target_projects={self.target_projects})")

    def generate_query(self, topic: str) -> Optional[str]:
        """
        Generate SRA search query using LLM.

        Args:
            topic: Research topic in natural language

        Returns:
            NCBI E-search compatible query string, or None if generation failed
        """
        logger.info(f"Generating search query for topic: {topic}")

        # Expand query with synonyms
        from src.query.expander import QueryExpander
        expander = QueryExpander()
        expanded_concepts = expander.expand_query_string(topic)

        expansion_text = ""
        if expanded_concepts:
            expansion_text = (
                "\nSYNONYM SUGGESTIONS (Use if applicable):\n"
                + "\n".join(expanded_concepts) + "\n"
            )
            logger.debug(f"Query expanded with {len(expanded_concepts)} synonym groups")

        prompt = f"""
You are an expert bioinformatics assistant specialized in NCBI SRA data mining.
Your goal is to generate a **HIGH RECALL** search strategy for the topic: "{topic}".

{expansion_text}
We want to find the LARGEST possible amount of relevant studies, even if it implies some noise.
It is preferable to retrieve 10,000 results and filter later, than to retrieve 0 by being too specific.

IMPORTANT: You must return the response in VALID JSON format.
Be very careful with quotes inside the query. If you use double quotes inside the query string, you MUST escape them with a backslash (\\").

Rules for HIGH RECALL:
1. **Use [All Fields] generously**: If a term is not strictly an organism or strategy, use it in [All Fields].
2. **Avoid unnecessary ANDs**: Do not connect all concepts with AND. If the user asks for "drought in tomato and pepper", use OR for the organisms.
3. **Expand synonyms**: Use OR for variants (e.g., "drought" OR "water stress" OR "water deficit").
4. **Organisms**: Identify the scientific name but also include the common name in [All Fields] (e.g., "Solanum lycopersicum"[Organism] OR "tomato"[All Fields]).
5. **Strategy**: If the user mentions a technique (e.g., RNA-Seq), use it, but include variants (e.g., "RNA-Seq"[Strategy] OR "transcriptome"[All Fields]).

Example of valid output:
{{
  "natural_query": "Broad search for RNA-Seq in tomato or pepper under water stress",
  "esearch_query": "(\\\"Solanum lycopersicum\\\"[Organism] OR \\\"tomato\\\"[All Fields]) AND (\\\"drought\\\"[All Fields] OR \\\"water stress\\\"[All Fields])"
}}

Return ONLY the JSON.
"""

        response = self.llm.generate_response(prompt, json_mode=True)

        if not response:
            logger.error("LLM failed to generate query")
            return None

        # Parse JSON response
        try:
            data = json.loads(response)
        except json.JSONDecodeError:
            # Try to extract JSON block from response
            logger.warning("Failed to parse JSON, attempting extraction")
            try:
                match = re.search(r'\{.*\}', response, re.DOTALL)
                if match:
                    json_str = match.group(0)
                    data = json.loads(json_str)
                else:
                    raise ValueError("No JSON block found in response")
            except Exception as e:
                logger.error(f"Error parsing query JSON: {e}")
                logger.debug(f"Raw response: {response}")
                return None

        natural_query = data.get('natural_query', '')
        esearch_query = data.get('esearch_query', '')

        logger.info(f"Natural Query: {natural_query}")
        logger.info(f"E-search Query: {esearch_query}")

        return esearch_query

    def enrich_metadata(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """
        Enrich metadata with additional counts from NCBI.

        Fetches:
        - Number of SRA experiments
        - Number of BioSamples

        Args:
            metadata: Base metadata dictionary

        Returns:
            Enriched metadata dictionary
        """
        bioproject = metadata.get('bioproject')
        if not bioproject:
            return metadata

        try:
            # Count SRA Experiments
            result = safe_esearch(db="sra", term=bioproject, retmax=0)
            metadata['sra_experiment_count'] = int(result['Count'])

            # Count BioSamples
            result = safe_esearch(db="biosample", term=bioproject, retmax=0)
            metadata['biosample_count'] = int(result['Count'])

            logger.debug(
                f"Enriched {bioproject}: "
                f"{metadata['sra_experiment_count']} experiments, "
                f"{metadata['biosample_count']} biosamples"
            )

        except Exception as e:
            logger.warning(f"Error fetching counts for {bioproject}: {e}")
            metadata['sra_experiment_count'] = 0
            metadata['biosample_count'] = 0

        return metadata

    def analyze_bioproject(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze and extract structured information from BioProject metadata using LLM.

        Args:
            metadata: BioProject metadata dictionary

        Returns:
            Dictionary with extracted analysis (tissues, cultivars, design, etc.)
        """
        from src.llm.prompts import get_analysis_prompt

        # Enrich with counts first
        metadata = self.enrich_metadata(metadata)

        # Build text representation for LLM
        if 'full_text' in metadata:
            text_data = metadata['full_text']
        else:
            text_data = "\n".join([f"{k}: {v}" for k, v in metadata.items()])

        system_prompt, user_prompt = get_analysis_prompt(text_data)

        response = self.llm.generate_response(
            user_prompt,
            system_prompt=system_prompt,
            json_mode=True
        )

        if not response:
            logger.warning(f"LLM analysis failed for {metadata.get('bioproject')}")
            return {}

        try:
            analysis = json.loads(response)
            return analysis
        except json.JSONDecodeError as e:
            logger.error(
                f"Failed to parse LLM JSON response for {metadata.get('bioproject')}: {e}"
            )
            logger.debug(f"Raw response: {response}")
            return {}

    def run(self, topic: str, max_workers: int = 15) -> List[Dict[str, Any]]:
        """
        Run the complete pipeline.

        Args:
            topic: Research topic in natural language
            max_workers: Number of parallel workers for LLM analysis

        Returns:
            List of analyzed BioProject records
        """
        # 1. Generate Query
        query = self.generate_query(topic)
        if not query:
            logger.error("Query generation failed, cannot proceed")
            return []

        logger.info("")
        logger.info("--- Starting SRA Search ---")
        logger.info(f"Query: {query}")

        # Pagination loop
        retstart = 0
        total_found = 0

        with tqdm(
            total=self.target_projects,
            desc="Collecting unique BioProjects",
            unit="project"
        ) as pbar:
            while len(self.bioprojects_seen) < self.target_projects:
                logger.debug(
                    f"Fetching batch starting at {retstart} "
                    f"({len(self.bioprojects_seen)}/{self.target_projects} collected)"
                )

                # 2. Search SRA (Batch)
                try:
                    search_result = search_sra(query, retstart=retstart, retmax=BATCH_SIZE)
                except Exception as e:
                    logger.error(f"Search failed: {e}")
                    break

                if not search_result or 'IdList' not in search_result:
                    logger.warning("Search failed or returned no results")
                    break

                id_list = search_result['IdList']
                count = int(search_result['Count'])
                total_found = count

                if not id_list:
                    logger.info("No more IDs returned in this batch")
                    break

                logger.info(
                    f"Batch: {len(id_list)} records (total available: {count})"
                )

                # 3. Process in sub-batches
                chunk_size = PipelineDefaults.FETCH_CHUNK_SIZE
                new_unique_records = []

                for i in range(0, len(id_list), chunk_size):
                    chunk_ids = id_list[i:i+chunk_size]

                    try:
                        xml_content = fetch_details(chunk_ids)
                    except Exception as e:
                        logger.error(f"Error fetching details for chunk: {e}")
                        continue

                    if not xml_content:
                        continue

                    # Parse XML
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
                            pbar.update(1)

                        new_unique_records.append(metadata)

                        if len(self.bioprojects_seen) >= self.target_projects:
                            break

                    if len(self.bioprojects_seen) >= self.target_projects:
                        break

                logger.info(f"Found {len(new_unique_records)} new unique BioProjects in this batch")

                # 4. Parallel LLM Analysis
                if new_unique_records:
                    logger.info(f"Analyzing {len(new_unique_records)} projects with {max_workers} workers...")

                    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                        # Submit all tasks
                        future_to_meta = {
                            executor.submit(self.analyze_bioproject, meta): meta
                            for meta in new_unique_records
                        }

                        # Process completed tasks with progress bar
                        with tqdm(
                            total=len(future_to_meta),
                            desc="Analyzing projects",
                            unit="project",
                            leave=False
                        ) as analysis_pbar:
                            for future in concurrent.futures.as_completed(future_to_meta):
                                meta = future_to_meta[future]
                                bioproject = meta.get('bioproject', 'Unknown')

                                try:
                                    analysis = future.result()
                                    full_record = {**meta, **analysis}
                                    self.results.append(full_record)
                                    logger.debug(f"Analyzed: {bioproject}")

                                except Exception as exc:
                                    logger.error(f"Error analyzing {bioproject}: {exc}")

                                finally:
                                    analysis_pbar.update(1)

                # Check if we've processed all available results
                if retstart + len(id_list) >= count:
                    logger.info("Reached end of search results")
                    break

                # Check if target reached
                if len(self.bioprojects_seen) >= self.target_projects:
                    logger.info(f"Target of {self.target_projects} unique BioProjects reached!")
                    break

                # Move to next batch
                retstart += len(id_list)

        logger.info("")
        logger.info(f"Pipeline completed: {len(self.results)} projects analyzed")
        return self.results
