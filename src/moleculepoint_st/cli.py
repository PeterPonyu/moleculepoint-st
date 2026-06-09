from __future__ import annotations

import argparse
import json
from pathlib import Path
from .claim_status import graduation_claim_status
from .contracts import ClaimGateEvidence, evaluate_claim_gate
from .real_smoke import RealSmokeConfig, default_molecules_path
from .smoke import run_synthetic_smoke


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="moleculepoint-st")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("smoke-synthetic")
    real = sub.add_parser("smoke-real")
    real.add_argument("--molecules-path", "--st-path", dest="molecules_path", type=Path, default=None)
    real.add_argument("--max-genes", type=int, default=8)
    real.add_argument("--min-gene-count", type=int, default=40)
    real.add_argument("--max-points-per-gene", type=int, default=600)
    real.add_argument("--csr-replicates", type=int, default=3)
    real.add_argument("--seed", type=int, default=13)
    real.add_argument("--results-dir", type=Path, default=None)
    gate2 = sub.add_parser("gate2-parity")
    gate2.add_argument("--molecules-path", "--st-path", dest="molecules_path", type=Path, default=None)
    gate2.add_argument("--reference-path", type=Path, default=None)
    gate2.add_argument("--shared-metric-path", type=Path, default=None)
    gate2.add_argument("--out-path", type=Path, default=None)
    gate2.add_argument("--max-genes", type=int, default=8)
    gate2.add_argument("--min-gene-count", type=int, default=40)
    gate2.add_argument("--max-points-per-gene", type=int, default=600)
    gate2.add_argument("--csr-replicates", type=int, default=3)
    gate2.add_argument("--seed", type=int, default=13)
    gate2.add_argument("--results-dir", type=Path, default=None)
    gate3 = sub.add_parser("gate3-analysis")
    gate3.add_argument("--molecules-path", "--st-path", dest="molecules_path", type=Path, default=None)
    gate3.add_argument("--out-path", type=Path, default=None)
    gate3.add_argument("--max-genes", type=int, default=8)
    gate3.add_argument("--min-gene-count", type=int, default=40)
    gate3.add_argument("--max-points-per-gene", type=int, default=600)
    gate3.add_argument("--csr-replicates", type=int, default=3)
    gate3.add_argument("--seed", type=int, default=13)
    gate3.add_argument("--results-dir", type=Path, default=None)
    sub.add_parser("claim-status")
    args = parser.parse_args(argv)
    if args.command == "smoke-synthetic":
        report = run_synthetic_smoke()
        print(json.dumps({"metrics": report.metrics, "claim_status": report.claim_status.value}, sort_keys=True))
        return 0
    if args.command == "smoke-real":
        from .real_smoke import run_real_data_smoke
        from .result_wiring import emit_real_smoke_results

        config = _real_config_from_args(args)
        report = run_real_data_smoke(config)
        payload = report.to_jsonable()
        payload["contract_results"] = _paths_payload(emit_real_smoke_results(report, config, results_dir=args.results_dir))
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate2-parity":
        from .parity import Gate2ParityConfig, run_gate2_parity
        from .result_wiring import emit_gate2_results

        config = Gate2ParityConfig(smoke=_real_config_from_args(args), reference_path=args.reference_path, shared_metric_path=args.shared_metric_path)
        report = run_gate2_parity(config)
        payload = report.to_jsonable()
        outputs = {}
        if args.out_path is not None:
            outputs["parity_table"] = args.out_path
        payload["contract_results"] = _paths_payload(emit_gate2_results(report, config.smoke, results_dir=args.results_dir, outputs=outputs))
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate3-analysis":
        from .gate3 import Gate3Config, run_gate3_analysis
        from .result_wiring import emit_gate3_results

        config = Gate3Config(smoke=_real_config_from_args(args))
        report = run_gate3_analysis(config)
        payload = report.to_jsonable()
        outputs = {}
        if args.out_path is not None:
            outputs["gate3_table"] = args.out_path
        payload["contract_results"] = _paths_payload(emit_gate3_results(report, config.smoke, results_dir=args.results_dir, outputs=outputs))
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "claim-status":
        print(graduation_claim_status().value)
        return 0
    print(evaluate_claim_gate(ClaimGateEvidence()).value)
    return 0


def _real_config_from_args(args: argparse.Namespace) -> RealSmokeConfig:
    return RealSmokeConfig(
        molecules_path=args.molecules_path or default_molecules_path(),
        max_genes=args.max_genes,
        min_gene_count=args.min_gene_count,
        max_points_per_gene=args.max_points_per_gene,
        csr_replicates=args.csr_replicates,
        seed=args.seed,
    )


def _write_json_if_requested(path: Path | None, payload: dict[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _paths_payload(paths: dict[str, Path]) -> dict[str, str]:
    return {key: str(path) for key, path in paths.items()}


if __name__ == "__main__":
    raise SystemExit(main())
