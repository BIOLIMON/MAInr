import os

# Entrez Configuration
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY", "")

# LLM Configuration
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_URL = "http://localhost:11434"

# Search Configuration
TARGET_PROJECTS = 100000  # Goal: effectively unlimited
BATCH_SIZE = 10000        # Fetch large blocks at once
