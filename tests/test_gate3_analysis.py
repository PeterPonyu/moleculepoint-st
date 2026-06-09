from __future__ import annotations

import unittest

from moleculepoint_st.contracts import ClaimGateEvidence
from moleculepoint_st.gate3 import Gate3Result


class Gate3AnalysisUnitTests(unittest.TestCase):
    def test_gate3_result_reports_license_as_only_missing_item(self):
        result = Gate3Result(
            ablation={"removed_component": "null_correction"},
            failure_mode={"floor_statement": "low-count floor"},
            evidence=ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True),
        )
        self.assertEqual(result.claim_status.value, "locked")
        self.assertEqual(result.to_jsonable()["missing_claim_evidence"], ["license_review"])


if __name__ == "__main__":
    unittest.main()
