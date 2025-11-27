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
def safe_esearch(db, term, retstart=0, retmax=20):
    """
    Generic safe search with retries.
    """
    try:
        handle = Entrez.esearch(
            db=db,
            term=term,
            retstart=retstart,
            retmax=retmax
        )
        results = Entrez.read(handle)
        handle.close()
        return results
    except Exception as e:
        # print(f"Error searching {db} (failed attempt): {e}") # Too verbose for high concurrency
        raise

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((urllib.error.URLError, urllib.error.HTTPError, RuntimeError))
)
def search_sra(query, retstart=0, retmax=BATCH_SIZE):
    """
    Searches the SRA database and returns the list of IDs.
    """
    print(f"Searching SRA: {query} (Batch: {retstart}-{retstart+retmax})")
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
        print(f"Error searching SRA (failed attempt): {e}")
        raise # Re-raise for tenacity to catch

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((urllib.error.URLError, urllib.error.HTTPError, RuntimeError))
)
def fetch_summary(id_list):
    """
    Fetches summaries for a list of IDs.
    """
    try:
        handle = Entrez.esummary(db="sra", id=id_list, rettype="docsum")
        summary = Entrez.read(handle)
        handle.close()
        return summary
    except Exception as e:
        print(f"Error fetching summary for IDs (failed attempt): {e}")
        raise # Re-raise for tenacity to catch

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=4, max=60),
    retry=retry_if_exception_type((urllib.error.URLError, urllib.error.HTTPError, RuntimeError))
)
def fetch_details(id_list):
    """
    Fetches full details for a list of IDs using efetch.
    Provides more information than fetch_summary.
    """
    try:
        # rettype="full" and retmode="xml" retrieves the complete record
        handle = Entrez.efetch(db="sra", id=id_list, rettype="full", retmode="xml")
        results = handle.read() # Return raw XML string
        handle.close()
        return results
    except Exception as e:
        print(f"Error fetching details for IDs (failed attempt): {e}")
        raise

