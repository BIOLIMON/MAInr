"""
Configuration settings for MAInr.

This module provides default configuration values that can be overridden
via environment variables or command-line arguments.
"""

import os
from src.constants import PipelineDefaults, OllamaDefaults, NCBILimits

# NCBI Entrez Configuration
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY", "")

# LLM Configuration
OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", OllamaDefaults.DEFAULT_MODEL)
# Support comma-separated URLs for multi-GPU setups
OLLAMA_URL = os.getenv("OLLAMA_URL", OllamaDefaults.DEFAULT_URL)

# Search Configuration
# Using more reasonable defaults - can be overridden
TARGET_PROJECTS = int(os.getenv("TARGET_PROJECTS", PipelineDefaults.DEFAULT_TARGET_PROJECTS))
BATCH_SIZE = int(os.getenv("BATCH_SIZE", NCBILimits.RECOMMENDED_BATCH_SIZE))
