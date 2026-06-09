from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from moleculepoint_st.contracts import ClaimGateEvidence
from moleculepoint_st.parity import Gate2ParityResult
from moleculepoint_st.real_smoke import RealSmokeConfig
from moleculepoint_st.result_wiring import emit_gate2_results


class ResultsContractWiringTests(unittest.TestCase):
    def test_vendored_contract_matches_recorded_sha256(self):
        package_dir = Path(__file__).resolve().parents[1] / "src" / "moleculepoint_st"
        contract_path = package_dir / "results_contract.py"
        expected = (package_dir / "results_contract.sha256").read_text(encoding="utf-8").strip().split()[0]
        observed = hashlib.sha256(contract_path.read_bytes()).hexdigest()
        self.assertEqual(observed, expected)

    def test_gate2_emits_contract_files(self):
        report = Gate2ParityResult(
            rows=[
                {
                    "method_id": "local_segmentation_light_csr",
                    "provenance": "RAN",
                    "same_molecules": True,
                    "segmentation_required": False,
                    "shared_metric_gene_count": 135,
                    "mean_nn_csr_ratio": 0.535397,
                },
                {
                    "method_id": "segmented_reference",
                    "provenance": "RAN",
                    "segmentation_required": True,
                    "shared_metric_gene_count": 135,
                },
            ],
            differentiator={
                "shared_metric_agreement": {
                    "spearman_rho": 0.813604,
                    "pearson_r": 0.793427,
                    "gene_count_finite": 135,
                }
            },
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        config = RealSmokeConfig(molecules_path=Path("data/processed/example/molecules.parquet"))
        with tempfile.TemporaryDirectory() as tmp:
            paths = emit_gate2_results(report, config, results_dir=Path(tmp))
            metrics = json.loads(paths["metrics"].read_text(encoding="utf-8"))
            self.assertEqual(metrics["project"], "moleculepoint-st")
            self.assertEqual(metrics["metrics"]["gate2.shared_metric_agreement.spearman_rho"], 0.813604)
            self.assertEqual(metrics["metrics"]["gate2.local_segmentation_light_csr.mean_nn_csr_ratio"], 0.535397)
            self.assertTrue(paths["run_metadata"].is_file())


if __name__ == "__main__":
    unittest.main()
