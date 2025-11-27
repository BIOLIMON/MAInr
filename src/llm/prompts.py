def get_analysis_prompt(text_data):
    """
    Generates the prompt for analyzing BioProject metadata.
    """
    system_prompt = """You are an expert bioinformatics data curator specializing in plant biology.
Your task is to extract structured and accurate metadata from the description of an SRA study.
You must be rigorous and base your extraction ONLY on the provided text.

Return a VALID JSON object with the following keys:
- "tissues": List of plant tissues studied (e.g., ["root", "leaf", "fruit"]). Use standardized English terms. If not mentioned, use [].
- "cultivars": List of cultivars, genotypes, or lines mentioned (e.g., ["M82", "Heinz 1706", "Col-0"]). If none, use [].
- "time_course_days": MAXIMUM duration of the study in DAYS (float). 
    - If it is a time series (e.g., "0h, 6h, 24h"), convert the maximum time to days (24h = 1.0).
    - If it is "7 days", it is 7.0.
    - If it is NOT a time series, return 0.0.
- "library_count_stated": Integer number of libraries/samples that the text SAYS were sequenced (e.g., "we sequenced 12 libraries"). If not explicitly mentioned, return null.
- "experiment_count_stated": Integer number of distinct conditions or experiments mentioned.
- "summary": A concise 1-sentence summary in English about the study's objective.
"""

    user_prompt = f"""Analyze the following study and extract the required metadata in JSON:

{text_data}"""

    return system_prompt, user_prompt
