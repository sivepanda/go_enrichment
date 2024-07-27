import unittest
from go_enrichment.enrichment import GOEnrichment

class TestGOEnrichment(unittest.TestCase):
    def test_enrichment(self):
        enrichment = GOEnrichment()
        results = enrichment.perform_enrichment()
        self.assertFalse(results.empty, "The results should not be empty.")
        self.assertIn(Adjusted_P_Value, results.columns, "Results should contain Adjusted_P_Value column.")

if __name__ == '__main__':
    unittest.main()
