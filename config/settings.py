import os

# Configuración de Entrez
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY", "")

# Configuración del LLM
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_URL = "http://localhost:11434"

# Configuración de Búsqueda
TARGET_PROJECTS = 100000  # Meta: efectivamente ilimitado
BATCH_SIZE = 10000        # Obtener grandes bloques a la vez
