import tempfile
import unittest
from pathlib import Path

from moleculepoint_st.contracts import ClaimGateEvidence
from moleculepoint_st.data_paths import find_repo_root, processed_data_path
from moleculepoint_st.real_smoke import RealSmokeConfig, RealSmokeResult, _mean_nearest_neighbor


class RealDataSmokeUnitTests(unittest.TestCase):
    def test_result_keeps_claim_locked(self):
        result = RealSmokeResult({"mean_nn_csr_ratio": 0.5}, ClaimGateEvidence(public_data_smoke=True))
        self.assertEqual(result.claim_status.value, "locked")
        self.assertEqual(
            result.to_jsonable()["missing_claim_evidence"],
            ["baseline_comparison", "ablation", "failure_modes", "license_review"],
        )

    def test_config_rejects_empty_gene_count(self):
        with self.assertRaises(ValueError):
            RealSmokeConfig(molecules_path=Path("molecules.parquet"), max_genes=0).validate()

    def test_repo_root_path_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data" / "processed").mkdir(parents=True)
            anchor = root / "src" / "package" / "module.py"
            anchor.parent.mkdir(parents=True)
            anchor.touch()
            self.assertEqual(find_repo_root(anchor), root)
            self.assertEqual(
                processed_data_path("fixture_card", anchor=anchor),
                root / "data" / "processed" / "fixture_card",
            )

    def test_nearest_neighbor_distance_on_square(self):
        import numpy as np

        xy = np.array([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
        self.assertAlmostEqual(_mean_nearest_neighbor(xy), 4.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
