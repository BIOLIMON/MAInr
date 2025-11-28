"""
LLM Client for interacting with Ollama.

This module provides a robust client for communicating with Ollama instances,
supporting multi-host load balancing and automatic failover.
"""

import os
import random
from typing import Optional, List, Dict, Any

import ollama

from config.settings import OLLAMA_MODEL, OLLAMA_URL
from src.utils.logger import get_logger
from src.constants import OllamaDefaults

logger = get_logger(__name__)


class LLMClient:
    """
    Client for interacting with Ollama LLM instances.

    Supports multiple hosts for load balancing and automatic failover.
    Uses random selection for simple but effective load distribution.

    Attributes:
        model: Name of the Ollama model to use
        num_threads: Number of threads for Ollama inference (optional)
        hosts: List of Ollama host URLs
    """

    def __init__(
        self,
        model: str = OLLAMA_MODEL,
        num_threads: Optional[int] = None
    ) -> None:
        """
        Initialize the LLM client.

        Args:
            model: Ollama model name (default from config)
            num_threads: Number of threads for inference (optional)
        """
        self.model = model
        self.num_threads = num_threads

        # Parse URLs (comma-separated for multi-GPU setups)
        # Re-read from env var to support runtime updates
        env_urls = os.getenv("OLLAMA_URL", OLLAMA_URL)
        self.hosts = [url.strip() for url in env_urls.split(',') if url.strip()]

        if not self.hosts:
            logger.warning("No Ollama URLs configured, using default")
            self.hosts = [OllamaDefaults.DEFAULT_URL]

        logger.info(f"LLM Client initialized with {len(self.hosts)} host(s): {self.hosts}")
        logger.info(f"Using model: {self.model}")

    def generate_response(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        json_mode: bool = False
    ) -> Optional[str]:
        """
        Generate a response from the LLM.

        Tries multiple hosts if configured, using random selection for
        load balancing. Returns None if all hosts fail.

        Args:
            prompt: User prompt to send to the LLM
            system_prompt: Optional system prompt for context
            json_mode: Whether to request JSON-formatted output

        Returns:
            LLM response as string, or None if all hosts failed
        """
        # Build messages
        messages: List[Dict[str, str]] = []
        if system_prompt:
            messages.append({'role': 'system', 'content': system_prompt})
        messages.append({'role': 'user', 'content': prompt})

        # Build options
        options: Dict[str, Any] = {'temperature': OllamaDefaults.DEFAULT_TEMPERATURE}
        if self.num_threads:
            options['num_threads'] = self.num_threads

        # Build request kwargs
        kwargs: Dict[str, Any] = {
            'model': self.model,
            'messages': messages,
            'options': options
        }

        if json_mode:
            kwargs['format'] = 'json'

        # Try each host once (random order for load balancing)
        hosts_to_try = list(self.hosts)
        random.shuffle(hosts_to_try)
        
        last_error: Optional[Exception] = None

        for host in hosts_to_try:

            try:
                client = ollama.Client(host=host)
                response = client.chat(**kwargs)

                # Success
                logger.debug(f"Successfully got response from {host}")
                return response['message']['content']

            except ollama.ResponseError as e:
                logger.warning(f"Ollama response error from {host}: {e}")
                last_error = e

            except ConnectionError as e:
                logger.warning(f"Connection error to {host}: {e}")
                last_error = e

            except Exception as e:
                logger.error(f"Unexpected error from {host}: {type(e).__name__}: {e}")
                last_error = e

        # All hosts failed
        logger.error(
            f"All {len(self.hosts)} Ollama host(s) failed. "
            f"Last error: {type(last_error).__name__}: {last_error}"
        )
        return None
