from __future__ import annotations

import unittest

import numpy as np

from moleculepoint_st.contracts import ClaimGateEvidence
from moleculepoint_st.parity import Gate2ParityResult, _neighbor_count_rate, _ripley_h_l_max


class Gate2ParityUnitTests(unittest.TestCase):
    def test_neighbor_count_rate_uses_radius(self):
        xy = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 0.0]])
        self.assertAlmostEqual(_neighbor_count_rate(xy, radius=1.5), 2.0 / 3.0)

    def test_shared_metric_is_finite_on_clustered_points(self):
        xy = np.array([[0.0, 0.0], [0.2, 0.0], [0.0, 0.2], [8.0, 8.0]])
        score = _ripley_h_l_max(xy, (0.0, 10.0, 0.0, 10.0), max_radii=4)
        self.assertTrue(np.isfinite(score))

    def test_full_gate_evidence_stays_locked_without_license(self):
        result = Gate2ParityResult(
            rows=[{"method_id": "local", "provenance": "RAN"}],
            differentiator={"same_data_ran_baseline_present": True},
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        payload = result.to_jsonable()
        self.assertEqual(payload["claim_status"], "locked")
        self.assertEqual(payload["missing_claim_evidence"], ["license_review"])


if __name__ == "__main__":
    unittest.main()
