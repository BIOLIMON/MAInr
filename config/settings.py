import os

# Configuración de Entrez
ENTREZ_EMAIL = ""
ENTREZ_API_KEY = ""

# Configuración del LLM
OLLAMA_MODEL = "qwen2.5:14b"
# Lista de hosts de Ollama para balanceo de carga.
# Si tienes múltiples GPUs, ejecuta múltiples instancias de Ollama en diferentes puertos.
OLLAMA_HOSTS = [
    "http://localhost:11434",
    # "http://localhost:11435", 
    # "http://localhost:11436"
]

# Configuración de Búsqueda
TARGET_PROJECTS = 100000  # Meta: efectivamente ilimitado
BATCH_SIZE = 10000        # Obtener grandes bloques a la vez
