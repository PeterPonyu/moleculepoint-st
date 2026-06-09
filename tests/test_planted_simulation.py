from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


def _load_simulator():
    script = Path(__file__).resolve().parents[1] / "scripts" / "data" / "simulate_point_patterns.py"
    spec = importlib.util.spec_from_file_location("simulate_point_patterns", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_planted_simulation_writes_truth_and_boundaries(tmp_path: Path):
    sim = _load_simulator()
    manifest = sim.simulate(
        sim.SimulationConfig(
            out_dir=tmp_path,
            seed=7,
            grid_size=2,
            n_genes=8,
            clustered_genes=4,
            coloc_pairs=2,
            molecules_per_gene=64,
        )
    )
    assert manifest["n_molecules"] == 512
    assert manifest["n_clustered_genes"] == 4
    assert manifest["n_colocalized_pairs"] == 2
    molecules = pd.read_parquet(tmp_path / "molecules.parquet")
    truth = pd.read_csv(tmp_path / "truth_genes.csv")
    cells = pd.read_parquet(tmp_path / "cell_boundaries.parquet")
    assert {"x", "y", "gene", "cell", "nucleus"}.issubset(molecules.columns)
    assert truth["is_clustered"].sum() == 4
    assert len(cells) == 4
