"""
NCBI SRA search functions.

This module provides functions for searching and fetching data from the NCBI SRA
database using Entrez utilities, with automatic retry logic.
"""

import urllib.error
from typing import List, Dict, Any

from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from config.settings import ENTREZ_EMAIL, ENTREZ_API_KEY, BATCH_SIZE
from src.utils.logger import get_logger
from src.constants import NCBILimits

logger = get_logger(__name__)

# Configure Entrez
if ENTREZ_EMAIL:
    Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY


@retry(
    stop=stop_after_attempt(NCBILimits.MAX_RETRIES),
    wait=wait_exponential(
        multiplier=1,
        min=NCBILimits.RETRY_MIN_WAIT,
        max=NCBILimits.RETRY_MAX_WAIT
    ),
    retry=retry_if_exception_type((
        urllib.error.URLError,
        urllib.error.HTTPError,
        RuntimeError
    ))
)
def safe_esearch(
    db: str,
    term: str,
    retstart: int = 0,
    retmax: int = 20
) -> Dict[str, Any]:
    """
    Generic safe search with automatic retries.

    Args:
        db: NCBI database name (e.g., 'sra', 'biosample')
        term: Search term/query
        retstart: Starting index for pagination (default: 0)
        retmax: Maximum results to return (default: 20)

    Returns:
        Parsed Entrez result dictionary

    Raises:
        urllib.error.URLError: Network-related errors (will retry)
        urllib.error.HTTPError: HTTP errors (will retry)
        RuntimeError: Other runtime errors (will retry)
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

        logger.debug(
            f"esearch({db}): term='{term[:50]}...', "
            f"found={results.get('Count', 0)} total"
        )

        return results

    except Exception as e:
        logger.debug(f"esearch retry: {type(e).__name__}: {e}")
        raise  # Re-raise for tenacity


@retry(
    stop=stop_after_attempt(NCBILimits.MAX_RETRIES),
    wait=wait_exponential(
        multiplier=1,
        min=NCBILimits.RETRY_MIN_WAIT,
        max=NCBILimits.RETRY_MAX_WAIT
    ),
    retry=retry_if_exception_type((
        urllib.error.URLError,
        urllib.error.HTTPError,
        RuntimeError
    ))
)
def search_sra(
    query: str,
    retstart: int = 0,
    retmax: int = BATCH_SIZE
) -> Dict[str, Any]:
    """
    Search the SRA database.

    Args:
        query: E-search compatible query string
        retstart: Starting index for pagination (default: 0)
        retmax: Maximum results to return (default: from config)

    Returns:
        Parsed Entrez result dictionary with 'IdList' and 'Count'

    Raises:
        urllib.error.URLError: Network-related errors (will retry)
        urllib.error.HTTPError: HTTP errors (will retry)
        RuntimeError: Other runtime errors (will retry)
    """
    logger.debug(f"Searching SRA: retstart={retstart}, retmax={retmax}")

    try:
        handle = Entrez.esearch(
            db="sra",
            term=query,
            retstart=retstart,
            retmax=retmax
        )
        results = Entrez.read(handle)
        handle.close()

        logger.debug(
            f"SRA search returned {len(results.get('IdList', []))} IDs "
            f"(total available: {results.get('Count', 0)})"
        )

        return results

    except Exception as e:
        logger.warning(f"SRA search retry: {type(e).__name__}: {e}")
        raise  # Re-raise for tenacity


@retry(
    stop=stop_after_attempt(NCBILimits.MAX_RETRIES),
    wait=wait_exponential(
        multiplier=1,
        min=NCBILimits.RETRY_MIN_WAIT,
        max=NCBILimits.RETRY_MAX_WAIT
    ),
    retry=retry_if_exception_type((
        urllib.error.URLError,
        urllib.error.HTTPError,
        RuntimeError
    ))
)
def fetch_summary(id_list: List[str]) -> List[Dict[str, Any]]:
    """
    Fetch summaries for a list of SRA IDs.

    Args:
        id_list: List of SRA ID strings

    Returns:
        List of summary dictionaries

    Raises:
        urllib.error.URLError: Network-related errors (will retry)
        urllib.error.HTTPError: HTTP errors (will retry)
        RuntimeError: Other runtime errors (will retry)
    """
    logger.debug(f"Fetching summaries for {len(id_list)} IDs")

    try:
        handle = Entrez.esummary(db="sra", id=id_list, rettype="docsum")
        summary = Entrez.read(handle)
        handle.close()

        logger.debug(f"Successfully fetched {len(summary)} summaries")
        return summary

    except Exception as e:
        logger.warning(f"Summary fetch retry: {type(e).__name__}: {e}")
        raise  # Re-raise for tenacity


@retry(
    stop=stop_after_attempt(NCBILimits.MAX_RETRIES),
    wait=wait_exponential(
        multiplier=1,
        min=NCBILimits.RETRY_MIN_WAIT,
        max=NCBILimits.RETRY_MAX_WAIT
    ),
    retry=retry_if_exception_type((
        urllib.error.URLError,
        urllib.error.HTTPError,
        RuntimeError
    ))
)
def fetch_details(id_list: List[str]) -> str:
    """
    Fetch full details for a list of SRA IDs using efetch.

    Provides more complete information than fetch_summary, including
    abstracts and experimental designs.

    Args:
        id_list: List of SRA ID strings

    Returns:
        Raw XML string containing detailed records

    Raises:
        urllib.error.URLError: Network-related errors (will retry)
        urllib.error.HTTPError: HTTP errors (will retry)
        RuntimeError: Other runtime errors (will retry)
    """
    logger.debug(f"Fetching details for {len(id_list)} IDs")

    try:
        # rettype="full" and retmode="xml" retrieves complete records
        handle = Entrez.efetch(
            db="sra",
            id=id_list,
            rettype="full",
            retmode="xml"
        )
        results = handle.read()  # Return raw XML string
        handle.close()

        logger.debug(f"Successfully fetched details ({len(results)} bytes XML)")
        return results

    except Exception as e:
        logger.warning(f"Details fetch retry: {type(e).__name__}: {e}")
        raise  # Re-raise for tenacity
