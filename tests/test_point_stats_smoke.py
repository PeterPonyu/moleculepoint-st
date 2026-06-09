import unittest
from moleculepoint_st import build_fixture, colocalization_fraction, run_synthetic_smoke, summarize_gene

class SmokeTests(unittest.TestCase):
    def test_colocalization_exceeds_null(self):
        points = build_fixture()
        self.assertGreater(colocalization_fraction(points, "GENE_A", "GENE_B"), colocalization_fraction(points, "GENE_A", "GENE_C"))
    def test_smoke_reports_low_count_instability(self):
        report = run_synthetic_smoke()
        self.assertEqual(report.claim_status.value, "locked")
        self.assertGreater(report.metrics["separation"], 0.5)
        self.assertEqual(summarize_gene(build_fixture(), "GENE_LOW").unstable_low_count, True)

if __name__ == "__main__":
    unittest.main()
