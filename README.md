# MAInr - Agente de Minería de Datos SRA

MAInr (Mining Agent for SRA) es una herramienta avanzada diseñada para buscar, extraer y analizar metadatos de la base de datos Sequence Read Archive (SRA) del NCBI. Utiliza Modelos de Lenguaje Grande (LLMs) para generar estrategias de búsqueda inteligentes y analizar el contexto de los estudios encontrados.

## Características

*   **Búsqueda Inteligente**: Genera consultas de búsqueda complejas y de alto recall utilizando modelos LLM (Mistral/Qwen) para asegurar que no se pierdan estudios relevantes.
*   **Minería Masiva**: Capaz de procesar miles de registros de SRA, manejando la paginación y la recuperación de datos de manera eficiente.
*   **Análisis Contextual**: Analiza los metadatos de cada BioProject utilizando LLMs para extraer información clave como condiciones experimentales, tejidos estudiados y si es una serie temporal.
*   **Procesamiento Paralelo**: Utiliza hilos para analizar múltiples proyectos simultáneamente, acelerando el proceso de minería.
*   **Salida Estructurada**: Genera archivos CSV con toda la información recopilada y analizada lista para su uso.

## Requisitos

*   Python 3.8+
*   Acceso a internet (para conectar con NCBI Entrez)
*   Una instancia de Ollama corriendo localmente (o accesible por red) con el modelo configurado (por defecto `qwen2.5:14b`).

## Instalación

1.  Clona este repositorio:
    ```bash
    git clone <url-del-repositorio>
    cd MAInr
    ```

2.  Instala las dependencias:
    ```bash
    pip install -r requirements.txt
    ```

3.  Configura Ollama:
    Asegúrate de tener Ollama instalado y ejecutándose.
    ```bash
    ollama serve
    ```
    Y descarga el modelo necesario:
    ```bash
    ollama pull qwen2.5:14b
    ```

## Configuración

El archivo `config/settings.py` contiene las configuraciones principales:

*   **ENTREZ_EMAIL**: Tu correo electrónico para el uso de la API de NCBI (Requerido por NCBI).
*   **ENTREZ_API_KEY**: Tu clave API de NCBI (Opcional pero recomendada para mayores límites de velocidad).
*   **OLLAMA_MODEL**: El modelo de Ollama a utilizar.
*   **OLLAMA_HOSTS**: Lista de hosts de Ollama. Soporta múltiples hosts para balanceo de carga si tienes varias GPUs.
*   **TARGET_PROJECTS**: Número objetivo de proyectos únicos a recuperar.
*   **BATCH_SIZE**: Tamaño del lote para las búsquedas en SRA.

## Uso
 
 Ejecuta el script principal desde la terminal pasando los argumentos necesarios:
 
 ```bash
 python3 mainr.py "tema de investigación" [opciones]
 ```
 
 ### Argumentos
 
 *   `topic`: El tema de investigación (ej. "drought stress in tomato"). Si no se proporciona, el sistema lo pedirá interactivamente.
 *   `-O`, `--output-dir`: Directorio donde se guardarán los resultados (por defecto: directorio actual).
 *   `-n`, `--num-workers`: Número de hilos para el procesamiento paralelo de análisis (por defecto: 15).
 *   `-t`, `--ollama-threads`: Número de hilos que utilizará Ollama para la inferencia (opcional).
 
 ### Ejemplos
 
 **Uso básico:**
 ```bash
 python3 mainr.py "drought stress in tomato"
 ```
 
 **Especificando directorio de salida y más hilos de trabajo:**
 ```bash
 python3 mainr.py "heat stress in arabidopsis" -O ./resultados -n 20
 ```
 
 **Controlando hilos de Ollama:**
 ```bash
 python3 mainr.py "cold stress in maize" -t 4
 ```
 
 ## Estructura del Proyecto
 
 *   `mainr.py`: Punto de entrada de la aplicación (CLI).
*   `config/`: Archivos de configuración.
*   `src/`: Código fuente.
    *   `llm/`: Cliente para interactuar con Ollama.
    *   `processing/`: Lógica principal del pipeline de procesamiento.
    *   `sra/`: Funciones para buscar y recuperar datos de NCBI SRA.
    *   `utils/`: Utilidades varias, incluyendo el parser de XML.

## Contribución

¡Las contribuciones son bienvenidas! Por favor, abre un issue o envía un pull request para mejoras o correcciones.

## Licencia

[MIT](LICENSE)
