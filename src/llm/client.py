import ollama
from config.settings import OLLAMA_MODEL, OLLAMA_URL

class LLMClient:
    def __init__(self, model=OLLAMA_MODEL, num_threads=None):
        self.model = model
        self.host = OLLAMA_URL
        self.num_threads = num_threads

    def get_client(self):
        return ollama.Client(host=self.host)

    def generate_response(self, prompt, system_prompt=None, json_mode=False):
        """
        Generates a response from the LLM.
        """
        messages = []
        if system_prompt:
            messages.append({'role': 'system', 'content': system_prompt})
        
        messages.append({'role': 'user', 'content': prompt})

        options = {'temperature': 0.1}
        if self.num_threads:
            options['num_thread'] = self.num_threads
        
        # Note: The Ollama python client might handle format='json' differently depending on the version.
        # For now we rely on the prompt to request JSON if needed,
        # but we can pass format='json' if the library supports it.
        kwargs = {
            'model': self.model,
            'messages': messages,
            'options': options
        }
        
        if json_mode:
            kwargs['format'] = 'json'

        try:
            client = self.get_client()
            response = client.chat(**kwargs)
            return response['message']['content']
        except Exception as e:
            print(f"Error communicating with LLM: {e}")
            return None
