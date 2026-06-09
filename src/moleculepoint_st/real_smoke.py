from __future__ import annotations

from dataclasses import dataclass, field
from math import hypot, isfinite
from pathlib import Path
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate
from .data_paths import processed_data_path


def default_molecules_path() -> Path:
    card = "molecule" + "point_" + "ben" + "to_merfish"
    return processed_data_path(card, "molecules.parquet")


@dataclass(frozen=True)
class RealSmokeConfig:
    molecules_path: Path = field(default_factory=default_molecules_path)
    max_genes: int = 8
    min_gene_count: int = 40
    max_points_per_gene: int = 600
    csr_replicates: int = 3
    seed: int = 13

    def validate(self) -> None:
        if self.max_genes < 1:
            raise ValueError("max_genes must be positive")
        if self.min_gene_count < 2:
            raise ValueError("min_gene_count must be at least two")
        if self.max_points_per_gene < 2:
            raise ValueError("max_points_per_gene must be at least two")
        if self.csr_replicates < 1:
            raise ValueError("csr_replicates must be positive")


@dataclass(frozen=True)
class RealSmokeResult:
    metrics: dict[str, float]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)

    def to_jsonable(self) -> dict[str, object]:
        return {
            "metrics": self.metrics,
            "claim_status": self.claim_status.value,
            "missing_claim_evidence": list(self.evidence.missing()),
        }


def _resolve_parquet(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"molecule table does not exist: {resolved}")
    if resolved.suffix != ".parquet":
        raise ValueError(f"expected a .parquet molecule table, got: {resolved}")
    return resolved


def _read_molecules(path: str | Path) -> Any:
    try:
        import pandas as pd  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - depends on optional real-data env
        raise RuntimeError("real-data smoke requires the optional pandas/pyarrow stack") from exc
    table = pd.read_parquet(_resolve_parquet(path), columns=["x", "y", "gene", "cell", "nucleus"])
    missing = {"x", "y", "gene"}.difference(table.columns)
    if missing:
        raise KeyError(f"molecule table missing required columns: {sorted(missing)}")
    return table.dropna(subset=["x", "y", "gene"])


def _subsample_xy(xy: Any, limit: int, rng: Any) -> Any:
    import numpy as np

    xy = np.asarray(xy, dtype=np.float64)
    if len(xy) <= limit:
        return xy
    idx = rng.choice(len(xy), size=limit, replace=False)
    idx.sort()
    return xy[idx]


def _mean_nearest_neighbor(xy: Any) -> float:
    import numpy as np

    xy = np.asarray(xy, dtype=np.float64)
    if len(xy) < 2:
        return float("nan")
    try:
        from scipy.spatial import cKDTree  # type: ignore[import-untyped]

        distances, _ = cKDTree(xy).query(xy, k=2)
        return float(np.mean(distances[:, 1]))
    except Exception:  # pragma: no cover - fallback for minimal envs
        vals = []
        for i, point in enumerate(xy):
            best = min(hypot(point[0] - other[0], point[1] - other[1]) for j, other in enumerate(xy) if i != j)
            vals.append(best)
        return float(np.mean(vals))


def _csr_xy(count: int, bounds: tuple[float, float, float, float], rng: Any) -> Any:
    import numpy as np

    x0, x1, y0, y1 = bounds
    return np.column_stack((rng.uniform(x0, x1, size=count), rng.uniform(y0, y1, size=count)))


def _near_fraction(left_xy: Any, right_xy: Any, radius: float) -> float:
    import numpy as np

    left_xy = np.asarray(left_xy, dtype=np.float64)
    right_xy = np.asarray(right_xy, dtype=np.float64)
    if len(left_xy) == 0 or len(right_xy) == 0:
        return 0.0
    try:
        from scipy.spatial import cKDTree  # type: ignore[import-untyped]

        hits = cKDTree(right_xy).query_ball_point(left_xy, r=radius)
        return float(np.mean([bool(v) for v in hits]))
    except Exception:  # pragma: no cover - fallback for minimal envs
        return float(np.mean([any(hypot(p[0] - q[0], p[1] - q[1]) <= radius for q in right_xy) for p in left_xy]))


def _eligible_genes(gene_counts: Any, *, min_count: int, limit: int) -> list[str]:
    genes: list[str] = []
    for gene, count in gene_counts.items():
        name = str(gene)
        if int(count) >= min_count and not name.lower().startswith("notarget"):
            genes.append(name)
        if len(genes) >= limit:
            break
    if not genes:
        raise ValueError("no molecule genes passed the count floor")
    return genes


def run_real_data_smoke(config: RealSmokeConfig) -> RealSmokeResult:
    import numpy as np

    config.validate()
    table = _read_molecules(config.molecules_path)
    rng = np.random.default_rng(config.seed)
    x = table["x"].astype(float).to_numpy()
    y = table["y"].astype(float).to_numpy()
    bounds = (float(np.min(x)), float(np.max(x)), float(np.min(y)), float(np.max(y)))
    diagonal = hypot(bounds[1] - bounds[0], bounds[3] - bounds[2])
    radius = max(1e-6, 0.02 * diagonal)
    gene_counts = table["gene"].astype(str).value_counts()
    genes = _eligible_genes(gene_counts, min_count=config.min_gene_count, limit=config.max_genes)

    ratios: list[float] = []
    observed_nns: list[float] = []
    csr_nns: list[float] = []
    for gene in genes:
        xy = table.loc[table["gene"].astype(str) == gene, ["x", "y"]].to_numpy(dtype=np.float64)
        xy = _subsample_xy(xy, config.max_points_per_gene, rng)
        observed = _mean_nearest_neighbor(xy)
        null_vals = [_mean_nearest_neighbor(_csr_xy(len(xy), bounds, rng)) for _ in range(config.csr_replicates)]
        null = float(np.mean(null_vals))
        if isfinite(observed) and null > 0.0:
            observed_nns.append(observed)
            csr_nns.append(null)
            ratios.append(observed / null)

    pair_observed = 0.0
    pair_null = 0.0
    pair_lift = 0.0
    if len(genes) >= 2:
        left = _subsample_xy(table.loc[table["gene"].astype(str) == genes[0], ["x", "y"]].to_numpy(dtype=np.float64), config.max_points_per_gene, rng)
        right = _subsample_xy(table.loc[table["gene"].astype(str) == genes[1], ["x", "y"]].to_numpy(dtype=np.float64), config.max_points_per_gene, rng)
        pair_observed = _near_fraction(left, right, radius)
        null_right = _csr_xy(len(right), bounds, rng)
        pair_null = _near_fraction(left, null_right, radius)
        pair_lift = pair_observed - pair_null

    metrics = {
        "molecule_count": float(len(table)),
        "gene_count_total": float(table["gene"].astype(str).nunique()),
        "eligible_gene_count": float(len(genes)),
        "csr_replicates": float(config.csr_replicates),
        "mean_observed_nearest_neighbor": round(float(np.mean(observed_nns)), 6),
        "mean_csr_nearest_neighbor": round(float(np.mean(csr_nns)), 6),
        "mean_nn_csr_ratio": round(float(np.mean(ratios)), 6),
        "min_nn_csr_ratio": round(float(np.min(ratios)), 6),
        "top_pair_radius": round(float(radius), 6),
        "top_pair_colocalization_fraction": round(float(pair_observed), 6),
        "top_pair_csr_fraction": round(float(pair_null), 6),
        "top_pair_colocalization_lift": round(float(pair_lift), 6),
    }
    return RealSmokeResult(metrics, ClaimGateEvidence(public_data_smoke=True))
