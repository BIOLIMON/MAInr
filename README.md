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
 
### Credenciales NCBI Entrez
 
Para usar la API de NCBI, necesitas proporcionar un correo electrónico (requerido) y opcionalmente una API Key. Puedes hacerlo de dos formas:
 
1.  **Variables de Entorno**:
    Establece las variables `ENTREZ_EMAIL` y `ENTREZ_API_KEY` en tu sistema.
    ```bash
    export ENTREZ_EMAIL="tu_email@ejemplo.com"
    export ENTREZ_API_KEY="tu_api_key"
    ```
 
2.  **Argumentos de Línea de Comandos**:
    Pásalos al ejecutar el script:
    ```bash
    python3 mainr.py --email "tu_email@ejemplo.com" --api-key "tu_api_key" ...
    ```
 
### Otras Configuraciones
 
El archivo `config/settings.py` contiene otras configuraciones por defecto que pueden ser sobreescritas por variables de entorno o modificando el archivo:
 
*   **OLLAMA_MODEL**: El modelo de Ollama a utilizar (default: `qwen2.5:14b`).
*   **OLLAMA_URL**: URL de la instancia de Ollama (default: `http://localhost:11434`).
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
