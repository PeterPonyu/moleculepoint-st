from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .contracts import ClaimGateEvidence
from .gate3 import Gate3Result
from .parity import Gate2ParityResult
from .real_smoke import RealSmokeConfig, RealSmokeResult
from .results_contract import dataset_card_id, write_results

PROJECT_ID = "moleculepoint-st"


def _numeric_items(payload: Mapping[str, Any], prefix: str = "") -> dict[str, float | None]:
    metrics: dict[str, float | None] = {}
    for key, value in payload.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if value is None:
            metrics[name] = None
        elif isinstance(value, bool):
            continue
        elif isinstance(value, (int, float)):
            metrics[name] = float(value)
        elif isinstance(value, Mapping):
            metrics.update(_numeric_items(value, name))
    return metrics


def _method_id(value: Any) -> str:
    return str(value or "method").replace(".", "_").replace(" ", "_").replace("/", "_")


def _config_paths(config: RealSmokeConfig) -> list[str]:
    return [str(config.molecules_path)]


def _metadata(config: RealSmokeConfig, report: RealSmokeResult | None, *, notes: str) -> dict[str, Any]:
    metrics = report.metrics if report is not None else {}
    return {
        "dataset_paths": _config_paths(config),
        "n_obs": metrics.get("molecule_count"),
        "n_vars": metrics.get("gene_count_total"),
        "seed": config.seed,
        "deterministic": True,
        "num_threads": 1,
        "reproducibility_level": "seeded",
        "normalization": {"applied": False, "method": "raw molecule coordinates"},
        "interpretability": {
            "claim_scope": "shared_metric_recovery_not_best_accuracy",
            "cell_segmentation_required": False,
        },
        "notes": notes,
        "provenance": {
            "max_genes": config.max_genes,
            "min_gene_count": config.min_gene_count,
            "max_points_per_gene": config.max_points_per_gene,
            "csr_replicates": config.csr_replicates,
        },
    }


def emit_real_smoke_results(
    report: RealSmokeResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        report.metrics,
        outputs=outputs,
        run_metadata=_metadata(config, report, notes="real-data point-pattern smoke emitted through the vendored results contract"),
        results_dir=results_dir,
    )


def emit_gate2_results(
    report: Gate2ParityResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics: dict[str, float | None] = {}
    metrics.update(_numeric_items(report.differentiator, "gate2"))
    for row in report.rows:
        method = _method_id(row.get("method_id"))
        metrics.update(_numeric_items(row, f"gate2.{method}"))
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(
            config,
            None,
            notes="gate-2 same-data shared-metric parity emitted through the vendored results contract",
        ),
        results_dir=results_dir,
    )


def emit_gate3_results(
    report: Gate3Result,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics: dict[str, float | None] = {}
    metrics.update(_numeric_items(report.ablation, "gate3.ablation"))
    metrics.update(_numeric_items(report.failure_mode, "gate3.failure_mode"))
    synthetic_report = RealSmokeResult(dict(metrics), evidence=ClaimGateEvidence(public_data_smoke=True))
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(
            config,
            synthetic_report,
            notes="gate-3 null-correction ablation and failure-mode floor emitted through the vendored results contract",
        ),
        results_dir=results_dir,
    )
