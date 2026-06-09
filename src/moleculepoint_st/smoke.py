from __future__ import annotations

import random
from .contracts import ClaimGateEvidence, MoleculePoint, PatternStabilityReport
from .point_stats import colocalization_fraction, summarize_gene


def build_fixture(seed: int = 3) -> list[MoleculePoint]:
    rng = random.Random(seed)
    points: list[MoleculePoint] = []
    idx = 0
    for center in [(1.0, 1.0), (4.0, 4.0), (7.0, 2.0)]:
        for _ in range(10):
            x = center[0] + rng.gauss(0, 0.08)
            y = center[1] + rng.gauss(0, 0.08)
            points.append(MoleculePoint(f"m{idx}", "GENE_A", x, y))
            idx += 1
            points.append(MoleculePoint(f"m{idx}", "GENE_B", x + rng.gauss(0, 0.05), y + rng.gauss(0, 0.05)))
            idx += 1
    for _ in range(30):
        points.append(MoleculePoint(f"m{idx}", "GENE_C", rng.uniform(0, 8), rng.uniform(0, 8)))
        idx += 1
    for _ in range(2):
        points.append(MoleculePoint(f"m{idx}", "GENE_LOW", rng.uniform(0, 8), rng.uniform(0, 8)))
        idx += 1
    return points


def run_synthetic_smoke() -> PatternStabilityReport:
    points = build_fixture()
    clustered = colocalization_fraction(points, "GENE_A", "GENE_B")
    null = colocalization_fraction(points, "GENE_A", "GENE_C")
    low = summarize_gene(points, "GENE_LOW")
    metrics = {
        "colocalized_fraction": round(clustered, 4),
        "null_fraction": round(null, 4),
        "separation": round(clustered - null, 4),
        "low_count_unstable": 1.0 if low.unstable_low_count else 0.0,
    }
    return PatternStabilityReport(metrics, ClaimGateEvidence(ablation=True, failure_modes=True))
