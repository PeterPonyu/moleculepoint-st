from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate
from .real_smoke import RealSmokeConfig, _csr_xy, _mean_nearest_neighbor, _read_molecules, _subsample_xy, run_real_data_smoke
from .parity import _bounds


@dataclass(frozen=True)
class Gate3Config:
    smoke: RealSmokeConfig = field(default_factory=RealSmokeConfig)


@dataclass(frozen=True)
class Gate3Result:
    ablation: dict[str, Any]
    failure_mode: dict[str, Any]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "ablation": self.ablation,
            "failure_mode": self.failure_mode,
            "claim_status": self.claim_status.value,
            "missing_claim_evidence": list(self.evidence.missing()),
        }


def _low_count_floor(config: RealSmokeConfig) -> dict[str, Any]:
    import numpy as np

    table = _read_molecules(config.molecules_path)
    rng = np.random.default_rng(config.seed + 104729)
    bounds = _bounds(table)
    counts = table["gene"].astype(str).value_counts()
    natural_low = [(str(gene), int(count)) for gene, count in counts.items() if 2 <= int(count) < config.min_gene_count and not str(gene).lower().startswith("notarget")]
    if natural_low:
        probes = [(gene, count, table.loc[table["gene"].astype(str) == gene, ["x", "y"]].to_numpy(dtype="float64")) for gene, count in natural_low[: max(config.max_genes, 1)]]
        probe_source = "natural_low_count_genes"
    else:
        probe_count = min(5, config.min_gene_count - 1)
        high = [(str(gene), int(count)) for gene, count in counts.items() if int(count) >= config.min_gene_count and not str(gene).lower().startswith("notarget")]
        probes = []
        for gene, _count in high[: max(config.max_genes, 1)]:
            full_xy = table.loc[table["gene"].astype(str) == gene, ["x", "y"]].to_numpy(dtype="float64")
            probes.append((gene, probe_count, _subsample_xy(full_xy, probe_count, rng)))
        probe_source = "real_gene_sparse_downsample"
    ratios: list[float] = []
    for _gene, _count, xy in probes:
        xy = _subsample_xy(xy, config.max_points_per_gene, rng)
        observed = _mean_nearest_neighbor(xy)
        null_vals = [_mean_nearest_neighbor(_csr_xy(len(xy), bounds, rng)) for _ in range(config.csr_replicates)]
        null = float(np.mean(null_vals))
        if isfinite(observed) and null > 0.0:
            ratios.append(observed / null)
    if ratios:
        ratio_mean: float | None = round(float(np.mean(ratios)), 6)
        ratio_std: float | None = round(float(np.std(ratios)), 6)
    else:
        ratio_mean = None
        ratio_std = None
    return {
        "probe_source": probe_source,
        "low_count_gene_count_scored": float(len(probes)),
        "low_count_floor": float(min((count for _gene, count, _xy in probes), default=0)),
        "low_count_ceiling": float(max((count for _gene, count, _xy in probes), default=0)),
        "low_count_mean_nn_csr_ratio": ratio_mean,
        "low_count_nn_csr_ratio_std": ratio_std,
        "floor_statement": "genes below the molecule count floor, or real genes thinned to that floor, have unstable point statistics and are excluded from the claim scope",
    }


def run_gate3_analysis(config: Gate3Config | None = None) -> Gate3Result:
    config = config or Gate3Config()
    smoke = run_real_data_smoke(config.smoke)
    metrics = smoke.metrics
    raw_pair = metrics["top_pair_colocalization_fraction"]
    corrected_pair = metrics["top_pair_colocalization_lift"]
    raw_nn = metrics["mean_observed_nearest_neighbor"]
    corrected_nn = metrics["mean_nn_csr_ratio"]
    ablation = {
        "removed_component": "CSR_null_correction",
        "raw_mean_nearest_neighbor": raw_nn,
        "csr_corrected_mean_nn_ratio": corrected_nn,
        "raw_top_pair_fraction": raw_pair,
        "csr_corrected_top_pair_lift": corrected_pair,
        "raw_pair_minus_corrected_lift": round(float(raw_pair - corrected_pair), 6),
        "interpretation": "without the null correction, raw fractions/distances are density-scale dependent and overstate the co-localization headline",
    }
    failure_mode = _low_count_floor(config.smoke)
    evidence = ClaimGateEvidence(public_data_smoke=True, baseline_comparison=True, ablation=True, failure_modes=True)
    return Gate3Result(ablation, failure_mode, evidence)
