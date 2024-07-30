import unittest
from go_enrichment.enrichment import enrichment_analysis
from go_enrichment.utils import parse_go_data

class TestEnrichment(unittest.TestCase):
    def test_enrichment(self):
        go_data = parse_go_data('test_data.json')
        results = enrichment_analysis(go_data)
        self.assertIsInstance(results, pd.DataFrame)
        self.assertTrue('GO_category' in results.columns)
        self.assertTrue('chromosome' in results.columns)
        self.assertTrue('odds_ratio' in results.columns)
        self.assertTrue('p_value' in results.columns)

if __name__ == '__main__':
    unittest.main()

