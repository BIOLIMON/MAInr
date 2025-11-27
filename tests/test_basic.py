import sys
import os
import unittest

# Add the project root to the path so we can import src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.sra.search import search_sra, fetch_summary
from Bio import Entrez

class TestSRASearch(unittest.TestCase):
    def setUp(self):
        # Set a default email for testing if not already set
        if not Entrez.email:
            Entrez.email = "n.mulleraguirre@gmail.com"

    def test_search_sra_basic(self):
        """Test a simple SRA search returns results."""
        query = '("Solanum lycopersicum"[Organism]) AND "RNA-Seq"[Strategy]'
        print(f"\nTesting search with query: {query}")
        result = search_sra(query, retmax=5)
        self.assertIsNotNone(result)
        self.assertIn('IdList', result)
        self.assertGreater(len(result['IdList']), 0)
        print(f"Found {len(result['IdList'])} results.")
        return result['IdList']

    def test_fetch_summary(self):
        """Test fetching summary for a known ID."""
        # Use a known ID or fetch one first
        query = '("Solanum lycopersicum"[Organism]) AND "RNA-Seq"[Strategy]'
        search_result = search_sra(query, retmax=1)
        if not search_result or not search_result['IdList']:
            self.skipTest("No IDs found to test fetch_summary")
        
        id_to_fetch = search_result['IdList'][0]
        print(f"\nTesting fetch_summary for ID: {id_to_fetch}")
        summary = fetch_summary([id_to_fetch])
        self.assertIsNotNone(summary)
        self.assertGreater(len(summary), 0)
        self.assertIn('Id', summary[0])
        self.assertEqual(str(summary[0]['Id']), str(id_to_fetch))
        print("Summary fetched successfully.")

if __name__ == '__main__':
    unittest.main()
