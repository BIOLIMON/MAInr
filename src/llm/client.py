import ollama
from config.settings import OLLAMA_MODEL, OLLAMA_URL

class LLMClient:
    def __init__(self, model=OLLAMA_MODEL):
        self.model = model
        self.host = OLLAMA_URL

    def get_client(self):
        return ollama.Client(host=self.host)

    def generate_response(self, prompt, system_prompt=None, json_mode=False):
        """
        Genera una respuesta del LLM.
        """
        messages = []
        if system_prompt:
            messages.append({'role': 'system', 'content': system_prompt})
        
        messages.append({'role': 'user', 'content': prompt})

        options = {'temperature': 0.1}
        
        # Nota: El cliente python de Ollama podría manejar format='json' de manera diferente dependiendo de la versión.
        # Por ahora confiaremos en el prompt para solicitar JSON si es necesario,
        # pero podemos pasar format='json' si la librería lo soporta.
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
            print(f"❌ Error comunicándose con el LLM: {e}")
            return None
