"""
Constants and default values for MAInr.

This module centralizes all constant values used throughout the application
to improve maintainability and prevent magic numbers in the code.
"""

from typing import Final


class NCBILimits:
    """NCBI API limits and constraints."""

    # NCBI recommends no more than 3 requests per second without API key
    # With API key: 10 requests per second
    RATE_LIMIT_NO_KEY: Final[int] = 3
    RATE_LIMIT_WITH_KEY: Final[int] = 10
    RATE_LIMIT_PERIOD: Final[int] = 1  # seconds

    # NCBI efetch may have issues with very large batch sizes
    MAX_BATCH_SIZE: Final[int] = 1000
    DEFAULT_BATCH_SIZE: Final[int] = 500
    RECOMMENDED_BATCH_SIZE: Final[int] = 100

    # Retry configuration
    MAX_RETRIES: Final[int] = 5
    RETRY_MIN_WAIT: Final[int] = 4  # seconds
    RETRY_MAX_WAIT: Final[int] = 60  # seconds


class PipelineDefaults:
    """Default values for pipeline configuration."""

    # Reasonable default for initial testing
    DEFAULT_TARGET_PROJECTS: Final[int] = 100

    # Maximum recommended for production without memory issues
    MAX_TARGET_PROJECTS: Final[int] = 10000

    # Chunk size for processing batches
    FETCH_CHUNK_SIZE: Final[int] = 50

    # Default number of worker threads
    DEFAULT_NUM_WORKERS: Final[int] = 15
    MAX_NUM_WORKERS: Final[int] = 100
    MIN_NUM_WORKERS: Final[int] = 1


class OllamaDefaults:
    """Default Ollama configuration."""

    DEFAULT_MODEL: Final[str] = "qwen2.5:14b"
    DEFAULT_URL: Final[str] = "http://localhost:11434"
    DEFAULT_TEMPERATURE: Final[float] = 0.1

    # Base port for multi-GPU cluster
    CLUSTER_BASE_PORT: Final[int] = 11500

    # Timeout for health checks
    HEALTH_CHECK_TIMEOUT: Final[int] = 30  # seconds
    HEALTH_CHECK_INTERVAL: Final[float] = 0.5  # seconds

    # Startup wait time for cluster (will be replaced by health checks)
    CLUSTER_STARTUP_WAIT: Final[int] = 5  # seconds


class GPUDefaults:
    """GPU and hardware detection defaults."""

    # Estimated VRAM per concurrent stream (GB)
    VRAM_PER_WORKER: Final[float] = 2.5

    # Model base size (GB) - for qwen2.5:14b
    MODEL_BASE_SIZE: Final[float] = 10.0

    # Safety margin for VRAM estimation (percentage)
    VRAM_SAFETY_MARGIN: Final[float] = 0.15  # 15%


class FileDefaults:
    """File and directory defaults."""

    DEFAULT_OUTPUT_DIR: Final[str] = "."
    CACHE_DIR: Final[str] = ".cache/mainr"
    LOG_DIR: Final[str] = "logs"
    CONFIG_DIR: Final[str] = "config"

    # Log file patterns
    MAIN_LOG_FILE: Final[str] = "mainr.log"
    OLLAMA_LOG_PATTERN: Final[str] = "ollama_gpu{}.log"

    # Cache file
    RESULTS_CACHE_FILE: Final[str] = "results_cache.pkl"


class ValidationLimits:
    """Validation limits for user inputs."""

    MIN_WORKERS: Final[int] = 1
    MAX_WORKERS: Final[int] = 100

    MIN_OLLAMA_THREADS: Final[int] = 1
    MAX_OLLAMA_THREADS: Final[int] = 64

    MIN_TARGET_PROJECTS: Final[int] = 1
    MAX_TARGET_PROJECTS: Final[int] = 100000

    MIN_BATCH_SIZE: Final[int] = 1
    MAX_BATCH_SIZE: Final[int] = 1000
