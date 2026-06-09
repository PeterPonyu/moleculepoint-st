from __future__ import annotations

import unittest

from moleculepoint_st.claim_status import (
    graduation_claim_status,
    graduation_claim_status_line,
    graduation_evidence,
    missing_evidence,
    signed_recovery_review,
    validation_marker_approved,
)
from moleculepoint_st.contracts import ClaimStatus


class ClaimStatusIntegrityTests(unittest.TestCase):
    def test_public_status_is_evidence_derived_preliminary(self):
        self.assertFalse(validation_marker_approved())
        self.assertFalse(signed_recovery_review())
        evidence = graduation_evidence()
        self.assertEqual(evidence.missing(), ())
        self.assertEqual(graduation_claim_status(), ClaimStatus.PRELIMINARY)
        self.assertTrue(missing_evidence())
        self.assertTrue(graduation_claim_status_line().startswith("preliminary missing="))

    def test_claim_status_does_not_require_governance_docs(self):
        self.assertEqual(graduation_claim_status().value, "preliminary")


if __name__ == "__main__":
    unittest.main()
