def get_analysis_prompt(text_data):
    """
    Generates the prompt for analyzing BioProject metadata.
    """
    system_prompt = """Eres un experto curador de datos bioinformáticos especializado en biología vegetal.
Tu tarea es extraer metadatos estructurados y precisos de la descripción de un estudio SRA.
Debes ser riguroso y basarte ÚNICAMENTE en el texto proporcionado.

Devuelve un objeto JSON VÁLIDO con las siguientes claves:
- "tissues": Lista de tejidos vegetales estudiados (ej: ["root", "leaf", "fruit"]). Usa términos en inglés estandarizados. Si no se menciona, usa [].
- "cultivars": Lista de cultivares, genotipos o líneas mencionadas (ej: ["M82", "Heinz 1706", "Col-0"]). Si no hay, usa [].
- "time_course_days": Duración MÁXIMA del estudio en DÍAS (float). 
    - Si es una serie temporal (ej: "0h, 6h, 24h"), convierte el tiempo máximo a días (24h = 1.0).
    - Si es "7 days", es 7.0.
    - Si NO es una serie temporal, devuelve 0.0.
- "library_count_stated": Número entero de librerías/muestras que el texto DICE que se secuenciaron (ej: "we sequenced 12 libraries"). Si no se menciona explícitamente un número total, devuelve null.
- "experiment_count_stated": Número entero de condiciones o experimentos distintos mencionados.
- "summary": Un resumen conciso de 1 oración en español sobre el objetivo del estudio.
"""

    user_prompt = f"""Analiza el siguiente estudio y extrae los metadatos requeridos en JSON:

{text_data}"""

    return system_prompt, user_prompt
