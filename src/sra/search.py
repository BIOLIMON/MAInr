from Bio import Entrez
from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY, BATCH_SIZE

Entrez.email = ENTREZ_EMAIL
Entrez.api_key = ENTREZ_API_KEY

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
        print(f"Error buscando en SRA: {e}")
        return None

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
        print(f"Error obteniendo resumen para ID {id_list}: {e}")
        return None
