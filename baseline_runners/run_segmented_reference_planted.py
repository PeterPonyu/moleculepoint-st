#!/usr/bin/env python
"""Run the frozen Bento checkout on the planted point-pattern benchmark.

This runner is intentionally outside src/tests/scripts: upstream imports and names stay isolated,
and only generic output files are consumed by track code or META review.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _import_bento(checkout: Path) -> tuple[Any, str]:
    sys.path.insert(0, str(checkout.resolve()))
    import dask

    original_set = dask.config.set

    def guarded_set(arg=None, *args, **kwargs):
        if isinstance(arg, dict) and arg.get("dataframe.query-planning") is False:
            arg = {k: v for k, v in arg.items() if k != "dataframe.query-planning"}
        return original_set(arg or {}, *args, **kwargs)

    dask.config.set = guarded_set
    try:
        import bento as bt  # type: ignore[import-not-found]
    finally:
        dask.config.set = original_set
    features = bt.tl.list_point_features()
    names = set(features.keys()) if isinstance(features, dict) else set(features.index)
    if "ripley" not in names:
        raise RuntimeError("Bento imported, but the ripley point feature is unavailable")
    return bt, getattr(bt, "__version__", "unknown")


def _git_head(checkout: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(checkout), "rev-parse", "HEAD"], check=True, capture_output=True, text=True)
    except Exception:
        return "unknown"
    return out.stdout.strip() or "unknown"


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    finite = np.isfinite(scores)
    labels = labels[finite]
    scores = scores[finite]
    positives = int(labels.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    ranked = labels[order]
    hits = np.cumsum(ranked)
    precision = hits / (np.arange(len(ranked)) + 1)
    return float(np.sum(precision[ranked]) / positives)


def _precision_recall_at_k(labels: np.ndarray, scores: np.ndarray, k: int) -> tuple[float, float]:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    finite = np.isfinite(scores)
    labels = labels[finite]
    scores = scores[finite]
    positives = int(labels.sum())
    if positives == 0 or k <= 0:
        return float("nan"), float("nan")
    k = min(k, len(labels))
    order = np.argsort(-scores, kind="mergesort")[:k]
    tp = int(labels[order].sum())
    return float(tp / k), float(tp / positives)


def _score_groups(bt: Any, molecules: pd.DataFrame, cells: pd.DataFrame, min_points: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    import geopandas as gpd
    from shapely.geometry import Point

    feature = bt.tl.RipleyStats("cell_boundaries")
    cell_meta = cells.set_index("cell").to_dict(orient="index")
    rows: list[dict[str, Any]] = []
    groups_scored = 0
    groups_skipped = 0
    for (cell_id, gene), df in molecules.groupby(["cell", "gene"], sort=True):
        count = int(len(df))
        if count < min_points:
            groups_skipped += 1
            continue
        meta = cell_meta[str(cell_id)]
        work = gpd.GeoDataFrame(
            df[["x", "y"]].copy(),
            geometry=[Point(float(x), float(y)) for x, y in zip(df["x"], df["y"], strict=True)],
        )
        work["cell_boundaries_span"] = float(meta["span"])
        work["cell_boundaries_minx"] = float(meta["minx"])
        work["cell_boundaries_miny"] = float(meta["miny"])
        work["cell_boundaries_maxx"] = float(meta["maxx"])
        work["cell_boundaries_maxy"] = float(meta["maxy"])
        work["cell_boundaries_area"] = float(meta["area"])
        stats = feature.extract(work)
        score = float(stats["l_max"])
        if np.isfinite(score):
            rows.append({"cell": str(cell_id), "gene": str(gene), "l_max": score, "molecule_count": count})
            groups_scored += 1
        else:
            groups_skipped += 1
    per_cell = pd.DataFrame(rows)
    summary = {"cell_gene_groups_scored": int(groups_scored), "cell_gene_groups_skipped": int(groups_skipped)}
    if per_cell.empty:
        return pd.DataFrame(columns=["gene", "segmented_score", "scored_molecule_count", "scored_cell_count"]), summary
    gene_rows = []
    for gene, df in per_cell.groupby("gene", sort=True):
        weights = df["molecule_count"].to_numpy(dtype=float)
        scores = df["l_max"].to_numpy(dtype=float)
        gene_rows.append(
            {
                "gene": str(gene),
                "segmented_score": float(np.average(scores, weights=weights)),
                "scored_molecule_count": int(weights.sum()),
                "scored_cell_count": int(len(df)),
            }
        )
    return pd.DataFrame(gene_rows), summary


def run(args: argparse.Namespace) -> dict[str, Any]:
    start = time.perf_counter()
    bt, version = _import_bento(args.checkout)
    molecules = pd.read_parquet(args.molecules_path)
    cells = pd.read_parquet(args.cell_boundaries_path)
    nuclei = pd.read_parquet(args.nucleus_boundaries_path)
    truth = pd.read_csv(args.truth_genes_path)
    gene_scores, scoring_summary = _score_groups(bt, molecules, cells, args.min_points_per_cell_gene)
    merged = truth.merge(gene_scores, on="gene", how="left")
    merged["segmented_score"] = merged["segmented_score"].astype(float)
    labels = merged["is_clustered"].to_numpy(dtype=bool)
    scores = merged["segmented_score"].to_numpy(dtype=float)
    positive_count = int(labels.sum())
    precision_at_pos, recall_at_pos = _precision_recall_at_k(labels, scores, positive_count)
    auprc = _average_precision(labels, scores)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scores_path = args.out_dir / "segmented_reference_gene_scores.csv"
    metrics_path = args.out_dir / "segmented_reference_recovery.json"
    provenance_path = args.out_dir / "provenance.json"
    merged.to_csv(scores_path, index=False)
    runtime = time.perf_counter() - start
    metrics = {
        "method_id": "Bento",
        "provenance": "RAN",
        "metric": "RipleyStats_l_max_weighted_by_cell_gene_molecule_count",
        "target": "planted_clustered_gene_recovery",
        "auprc": round(float(auprc), 6),
        "precision_at_planted_positive_count": round(float(precision_at_pos), 6),
        "recall_at_planted_positive_count": round(float(recall_at_pos), 6),
        "n_genes": int(len(merged)),
        "n_planted_clustered_genes": positive_count,
        "n_molecules": int(len(molecules)),
        "n_cells": int(len(cells)),
        "n_nucleus_boundaries": int(len(nuclei)),
        "genes_with_finite_scores": int(np.isfinite(scores).sum()),
        "min_points_per_cell_gene": int(args.min_points_per_cell_gene),
        "runtime_seconds": round(float(runtime), 3),
        "bento_version": version,
        "checkout_head": _git_head(args.checkout),
        **scoring_summary,
    }
    provenance = {
        "baseline": "Bento",
        "status": "ran",
        "source_checkout": str(args.checkout),
        "source_license": "BSD-2-Clause",
        "uses_cell_boundaries": True,
        "uses_nucleus_boundaries": True,
        "outputs": {"gene_scores": str(scores_path), "metrics": str(metrics_path)},
        **metrics,
    }
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"segmented_reference": metrics, "outputs": provenance["outputs"]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkout", type=Path, default=Path("../baselines/bento-tools-original"))
    parser.add_argument("--molecules-path", type=Path, required=True)
    parser.add_argument("--cell-boundaries-path", type=Path, required=True)
    parser.add_argument("--nucleus-boundaries-path", type=Path, required=True)
    parser.add_argument("--truth-genes-path", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-points-per-cell-gene", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        payload = run(args)
    except Exception as exc:
        print(f"Bento planted-GT reference run failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
