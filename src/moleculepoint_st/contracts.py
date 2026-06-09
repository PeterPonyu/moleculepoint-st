from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class ClaimStatus(str, Enum):
    LOCKED = "locked"
    PRELIMINARY = "preliminary"
    REVIEW_READY = "review_ready"
    VALIDATED = "validated"


@dataclass(frozen=True)
class ClaimGateEvidence:
    public_data_smoke: bool = False
    baseline_comparison: bool = False
    ablation: bool = False
    failure_modes: bool = False
    license_review: bool = False

    def missing(self) -> tuple[str, ...]:
        required = {
            "public_data_smoke": self.public_data_smoke,
            "baseline_comparison": self.baseline_comparison,
            "ablation": self.ablation,
            "failure_modes": self.failure_modes,
            "license_review": self.license_review,
        }
        return tuple(name for name, present in required.items() if not present)


def evaluate_claim_gate(evidence: ClaimGateEvidence, *, human_signed: bool = False) -> ClaimStatus:
    if evidence.missing():
        return ClaimStatus.LOCKED
    return ClaimStatus.VALIDATED if human_signed else ClaimStatus.REVIEW_READY



@dataclass(frozen=True)
class MoleculePoint:
    molecule_id: str
    gene: str
    x: float
    y: float
    sample_id: str = "synthetic"


@dataclass(frozen=True)
class PointPatternSummary:
    gene: str
    count: int
    mean_nearest_neighbor: float
    unstable_low_count: bool


@dataclass(frozen=True)
class PatternStabilityReport:
    metrics: dict[str, float]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)
