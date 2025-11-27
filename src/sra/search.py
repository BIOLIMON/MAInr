from Bio import Entrez
from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY, BATCH_SIZE
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import urllib.error

if ENTREZ_EMAIL:
    Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((urllib.error.URLError, urllib.error.HTTPError, RuntimeError))
)
def search_sra(query, retstart=0, retmax=BATCH_SIZE):
    """
    Busca en la base de datos SRA y devuelve la lista de IDs.
    """
    print(f"Buscando en SRA: {query} (Lote: {retstart}-{retstart+retmax})")
    try:
        handle = Entrez.esearch(
            db="sra",
            term=query,
            retstart=retstart,
            retmax=retmax
        )
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        print(f"Error buscando en SRA (intento fallido): {e}")
        raise # Re-raise para que tenacity lo capture

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((urllib.error.URLError, urllib.error.HTTPError, RuntimeError))
)
def fetch_summary(id_list):
    """
    Obtiene resúmenes para una lista de IDs (o un solo ID).
    Nota: A menudo es mejor obtener uno por uno si necesitamos procesarlos individualmente
    y manejar errores con gracia, o en lotes.
    """
    # Para esta implementación, obtendremos uno por uno en el bucle principal para coincidir con el flujo de la lógica original
    # pero proporcionamos una función para obtener un solo resumen aquí.
    try:
        handle = Entrez.esummary(db="sra", id=id_list, rettype="docsum")
        summary = Entrez.read(handle)
        handle.close()
        return summary
    except Exception as e:
        print(f"Error obteniendo resumen para ID {id_list} (intento fallido): {e}")
        raise # Re-raise para que tenacity lo capture
