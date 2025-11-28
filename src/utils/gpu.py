"""
GPU detection and Ollama cluster management.

This module handles GPU detection and automatic launching of Ollama instances
for multi-GPU setups, including health checks and proper cleanup.
"""

import os
import shutil
import subprocess
import time
import signal
from pathlib import Path
from typing import List, Optional, Tuple

import requests

from src.utils.logger import get_logger
from src.constants import OllamaDefaults, FileDefaults

logger = get_logger(__name__)


def get_gpu_count() -> int:
    """
    Detect the number of available NVIDIA GPUs using nvidia-smi.

    Returns:
        Number of GPUs detected, or 0 if nvidia-smi is not available
    """
    if not shutil.which("nvidia-smi"):
        logger.debug("nvidia-smi not found, no GPUs detected")
        return 0

    try:
        # List GPU UUIDs and count them
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=uuid", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            check=True
        )

        lines = result.stdout.strip().split('\n')
        gpu_count = len([l for l in lines if l.strip()])

        logger.info(f"Detected {gpu_count} NVIDIA GPU(s)")
        return gpu_count

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running nvidia-smi: {e}")
        return 0
    except Exception as e:
        logger.error(f"Unexpected error detecting GPUs: {e}")
        return 0


def wait_for_ollama(url: str, timeout: int = OllamaDefaults.HEALTH_CHECK_TIMEOUT) -> bool:
    """
    Wait for an Ollama instance to become ready.

    Args:
        url: Ollama server URL
        timeout: Maximum time to wait in seconds

    Returns:
        True if server is ready, False if timeout
    """
    start = time.time()
    interval = OllamaDefaults.HEALTH_CHECK_INTERVAL

    logger.debug(f"Waiting for Ollama at {url} (timeout: {timeout}s)")

    while time.time() - start < timeout:
        try:
            response = requests.get(f"{url}/api/tags", timeout=2)
            if response.status_code == 200:
                logger.info(f"Ollama at {url} is ready")
                return True
        except (requests.ConnectionError, requests.Timeout):
            pass
        except Exception as e:
            logger.warning(f"Unexpected error checking Ollama health: {e}")

        time.sleep(interval)

    logger.error(f"Timeout waiting for Ollama at {url}")
    return False


class OllamaCluster:
    """
    Manages a multi-GPU Ollama cluster.

    Automatically launches one Ollama instance per GPU with isolated
    environments and proper resource management.
    """

    def __init__(self, base_port: int = OllamaDefaults.CLUSTER_BASE_PORT) -> None:
        """
        Initialize the cluster manager.

        Args:
            base_port: Base port number for the first GPU instance
        """
        self.base_port = base_port
        self.processes: List[Tuple[subprocess.Popen, any]] = []
        self.urls: List[str] = []
        self.gpu_count = get_gpu_count()

        logger.info(f"OllamaCluster initialized (base_port={base_port}, gpus={self.gpu_count})")

    def start(self) -> List[str]:
        """
        Start an Ollama instance for each detected GPU.

        Returns:
            List of URLs for the started instances, empty list on failure
        """
        if self.gpu_count <= 1:
            logger.info("Single GPU or no GPU detected. Skipping cluster launch.")
            return []

        logger.info(f"Starting Ollama cluster for {self.gpu_count} GPUs...")

        # Locate models directory
        home = Path.home()
        user_models = home / ".ollama" / "models"
        system_models = Path("/usr/share/ollama/.ollama/models")

        # Determine which models directory to use
        orig_models = user_models
        if not (user_models / "manifests").exists() or not any((user_models / "manifests").iterdir()):
            if system_models.exists():
                logger.info(f"Using system-wide models from {system_models}")
                orig_models = system_models
            else:
                logger.warning(f"No models found in {user_models} or {system_models}")

        # Launch instances
        for i in range(self.gpu_count):
            port = self.base_port + i

            # Setup isolated environment
            iso_home = home / ".ollama_cluster" / f"gpu{i}"
            iso_models_dir = iso_home / ".ollama" / "models"

            try:
                # Create directory structure
                iso_models_dir.parent.mkdir(parents=True, exist_ok=True)

                # Symlink models to save space
                if orig_models.exists():
                    if iso_models_dir.exists():
                        if iso_models_dir.is_symlink():
                            iso_models_dir.unlink()
                        elif iso_models_dir.is_dir():
                            shutil.rmtree(iso_models_dir)

                    try:
                        iso_models_dir.symlink_to(orig_models)
                    except OSError as e:
                        logger.warning(f"Could not symlink models for GPU {i}: {e}")
                        # On Windows or without permissions, try copying
                        logger.info(f"Attempting to copy models instead...")
                        shutil.copytree(orig_models, iso_models_dir, dirs_exist_ok=True)

                # Prepare environment variables
                env = os.environ.copy()
                env["CUDA_VISIBLE_DEVICES"] = str(i)
                env["OLLAMA_HOST"] = f"0.0.0.0:{port}"
                env["HOME"] = str(iso_home)

                # Launch process
                logger.info(f"  Launching Ollama on GPU {i} (Port {port})...")

                log_file_path = FileDefaults.OLLAMA_LOG_PATTERN.format(i)
                log_file = open(log_file_path, "w")

                proc = subprocess.Popen(
                    ["ollama", "serve"],
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    preexec_fn=os.setsid  # Create new process group for clean kill
                )

                self.processes.append((proc, log_file))
                url = f"http://localhost:{port}"
                self.urls.append(url)

                logger.debug(f"  Process started for GPU {i}, PID: {proc.pid}")

            except FileNotFoundError:
                logger.error(f"  'ollama' command not found for GPU {i}")
            except Exception as e:
                logger.error(f"  Failed to launch on GPU {i}: {type(e).__name__}: {e}")

        if not self.urls:
            logger.error("Failed to start any Ollama instances")
            return []

        # Wait for all instances to be ready
        logger.info("Waiting for cluster to initialize...")
        ready_urls = []

        for url in self.urls:
            if wait_for_ollama(url):
                ready_urls.append(url)
            else:
                logger.warning(f"Instance at {url} failed to become ready")

        if len(ready_urls) < len(self.urls):
            logger.warning(f"Only {len(ready_urls)}/{len(self.urls)} instances became ready")

        logger.info(f"Cluster started with {len(ready_urls)} ready instance(s)")
        return ready_urls if ready_urls else self.urls  # Return all URLs even if health checks failed

    def stop(self) -> None:
        """
        Stop all started Ollama processes.
        """
        if not self.processes:
            logger.debug("No processes to stop")
            return

        logger.info("Stopping Ollama cluster...")

        for i, (proc, log_file) in enumerate(self.processes):
            try:
                # Try graceful shutdown first
                logger.debug(f"  Stopping process {i} (PID: {proc.pid})...")
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                proc.wait(timeout=5)
                logger.debug(f"  Process {i} stopped gracefully")

            except subprocess.TimeoutExpired:
                # Force kill if needed
                logger.warning(f"  Process {i} didn't stop gracefully, forcing...")
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except Exception as e:
                    logger.error(f"  Error force-killing process {i}: {e}")

            except Exception as e:
                logger.error(f"  Error stopping process {i}: {e}")

            finally:
                if log_file:
                    try:
                        log_file.close()
                    except:
                        pass

        self.processes = []
        self.urls = []
        logger.info("Cluster stopped successfully")


def get_ollama_launch_commands(
    gpu_count: int,
    base_port: int = OllamaDefaults.CLUSTER_BASE_PORT
) -> List[str]:
    """
    Generate shell commands to launch Ollama instances manually.

    Useful for documentation or manual setup.

    Args:
        gpu_count: Number of GPUs
        base_port: Base port number

    Returns:
        List of shell commands
    """
    commands = []
    for i in range(gpu_count):
        port = base_port + i
        cmd = f"CUDA_VISIBLE_DEVICES={i} OLLAMA_HOST=0.0.0.0:{port} ollama serve"
        commands.append(cmd)
    return commands
