import os
import ollama
from config.settings import OLLAMA_MODEL, OLLAMA_URL
import itertools
import random

class LLMClient:
    def __init__(self, model=OLLAMA_MODEL, num_threads=None):
        self.model = model
        self.num_threads = num_threads
        
        # Parse URLs (comma-separated)
        # Re-read from env var to support runtime updates (e.g. from mainr.py automation)
        env_urls = os.getenv("OLLAMA_URL", OLLAMA_URL)
        self.hosts = [url.strip() for url in env_urls.split(',') if url.strip()]
        if not self.hosts:
            self.hosts = ["http://localhost:11434"]
            
        # Round-robin iterator
        self.host_cycle = itertools.cycle(self.hosts)

    def get_next_client(self):
        """
        Returns an Ollama client for the next host in the cycle.
        """
        host = next(self.host_cycle)
        return ollama.Client(host=host), host

    def generate_response(self, prompt, system_prompt=None, json_mode=False):
        """
        Generates a response from the LLM, trying multiple hosts if necessary.
        """
        messages = []
        if system_prompt:
            messages.append({'role': 'system', 'content': system_prompt})
        
        messages.append({'role': 'user', 'content': prompt})

        options = {'temperature': 0.1}
        if self.num_threads:
            options['num_thread'] = self.num_threads
        
        kwargs = {
            'model': self.model,
            'messages': messages,
            'options': options
        }
        
        if json_mode:
            kwargs['format'] = 'json'

        # Try up to len(self.hosts) times (once per host)
        attempts = len(self.hosts)
        last_error = None
        
        # Start from a random position to distribute load if multiple clients are created
        # Actually itertools.cycle is per instance. If we create one LLMClient per thread, we want them desynchronized.
        # But Pipeline creates ONE LLMClient. So round-robin is fine for sequential calls.
        # Wait, Pipeline.run uses ThreadPoolExecutor. 
        # If Pipeline shares one LLMClient instance across threads, next(cycle) is not thread-safe?
        # itertools.cycle IS thread-safe in CPython (GIL), but let's be careful.
        # Better: Pick a random host for each request? Or just rely on cycle.
        # Random is statistically fine for load balancing.
        
        # Let's use a simple retry loop over the hosts
        # We'll shuffle the list for this request to avoid thundering herd if one is down
        current_hosts = list(self.hosts)
        # random.shuffle(current_hosts) # Optional, but good for distribution
        
        # Actually, for true round robin across threads, we need a lock.
        # But random choice is easier and effective for load balancing.
        
        for _ in range(attempts):
            # Simple Random Load Balancing
            host = random.choice(self.hosts)
            
            try:
                client = ollama.Client(host=host)
                response = client.chat(**kwargs)
                return response['message']['content']
            except Exception as e:
                # print(f"Error communicating with LLM at {host}: {e}")
                last_error = e
                # Try next host
                continue
        
        print(f"All LLM hosts failed. Last error: {last_error}")
        return None
