"""
Hardware detection utilities.

Provides functions for detecting and analyzing available hardware resources,
particularly GPU memory for estimating optimal concurrency settings.
"""

import subprocess
import shutil
from typing import List, Dict

from src.utils.logger import get_logger
from src.constants import GPUDefaults

logger = get_logger(__name__)


def get_gpu_info() -> List[Dict[str, any]]:
    """
    Detect available GPUs using nvidia-smi.

    Returns:
        List of dictionaries with GPU information:
        - name: GPU model name
        - total_memory_mb: Total VRAM in MB
        - free_memory_mb: Free VRAM in MB

        Returns empty list if nvidia-smi is not available or fails.
    """
    if not shutil.which("nvidia-smi"):
        logger.debug("nvidia-smi not found, cannot detect GPU info")
        return []

    try:
        # Query GPU details
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,memory.total,memory.free",
                "--format=csv,noheader,nounits"
            ],
            capture_output=True,
            text=True,
            check=True
        )

        gpus = []
        lines = result.stdout.strip().split('\n')

        for i, line in enumerate(lines):
            if not line:
                continue

            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 3:
                try:
                    gpu_info = {
                        "name": parts[0],
                        "total_memory_mb": int(parts[1]),
                        "free_memory_mb": int(parts[2])
                    }
                    gpus.append(gpu_info)
                    logger.debug(
                        f"GPU {i}: {gpu_info['name']} - "
                        f"{gpu_info['free_memory_mb']}MB free / {gpu_info['total_memory_mb']}MB total"
                    )
                except ValueError as e:
                    logger.warning(f"Could not parse GPU {i} info: {e}")

        logger.info(f"Detected {len(gpus)} GPU(s) with detailed info")
        return gpus

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running nvidia-smi: {e}")
        return []
    except Exception as e:
        logger.error(f"Unexpected error detecting GPUs: {type(e).__name__}: {e}")
        return []


def estimate_max_workers(
    gpus: List[Dict[str, any]],
    model_size_gb: float = GPUDefaults.VRAM_PER_WORKER,
    safety_margin: float = GPUDefaults.VRAM_SAFETY_MARGIN
) -> int:
    """
    Estimate the maximum number of concurrent workers based on free VRAM.

    This is a heuristic calculation assuming the model is already loaded
    and each concurrent stream requires additional VRAM for KV cache.

    Args:
        gpus: List of GPU info dictionaries from get_gpu_info()
        model_size_gb: Estimated VRAM per concurrent worker in GB (default: 2.5)
        safety_margin: Safety margin as a percentage (default: 0.15 = 15%)

    Returns:
        Estimated maximum workers (minimum 1)
    """
    if not gpus:
        logger.warning("No GPUs detected, defaulting to 1 worker (CPU mode)")
        return 1

    # Validate inputs
    if model_size_gb <= 0:
        logger.error(f"Invalid model_size_gb: {model_size_gb}, using default")
        model_size_gb = GPUDefaults.VRAM_PER_WORKER

    # Calculate total free VRAM
    total_free_vram_mb = sum(gpu['free_memory_mb'] for gpu in gpus)
    total_free_vram_gb = total_free_vram_mb / 1024

    # Apply safety margin
    usable_vram_gb = total_free_vram_gb * (1 - safety_margin)

    # Estimate workers
    estimated_workers = int(usable_vram_gb // model_size_gb)

    # Ensure at least 1 worker
    estimated_workers = max(1, estimated_workers)

    logger.info(
        f"VRAM estimation: {total_free_vram_gb:.2f}GB free, "
        f"{usable_vram_gb:.2f}GB usable (after {safety_margin*100:.0f}% margin), "
        f"~{estimated_workers} workers @ {model_size_gb}GB each"
    )

    return estimated_workers
