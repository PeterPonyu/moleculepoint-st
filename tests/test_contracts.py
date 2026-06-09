import unittest
from moleculepoint_st import ClaimGateEvidence, ClaimStatus, MoleculePoint, evaluate_claim_gate

class ContractTests(unittest.TestCase):
    def test_claim_lock(self):
        self.assertEqual(evaluate_claim_gate(ClaimGateEvidence(failure_modes=True)), ClaimStatus.LOCKED)
    def test_point_contract(self):
        p = MoleculePoint("m1", "G", 1.0, 2.0)
        self.assertEqual(p.sample_id, "synthetic")

if __name__ == "__main__":
    unittest.main()
