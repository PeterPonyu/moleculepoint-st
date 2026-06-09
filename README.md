# MoleculePoint-ST

MoleculePoint-ST computes segmentation-light single-molecule point-pattern summaries and shared-metric parity diagnostics.

This repository is a conservative public code surface: method implementation, command-line entry points, tests, and the byte-locked results schema. Background citations are listed in `BASELINE_REFERENCES.md`.

## Install

```bash
python -m pip install -e .
```

The lightweight unit tests run without bundled datasets. Real-data commands expect local spatial-omics inputs and expose `--help` for path overrides.

## Command-line usage

```bash
python -m moleculepoint_st.cli smoke-synthetic
python -m moleculepoint_st.cli smoke-real --help
python -m moleculepoint_st.cli gate2-parity --help
python -m moleculepoint_st.cli gate3-analysis --help
python -m moleculepoint_st.cli claim-status
```

Commands emit JSON to stdout. Gate commands also write uniform contract outputs under `results/<project>/` via the vendored `results_contract.py`.

## Evidence-derived status

`python -m moleculepoint_st.cli claim-status` reads `evidence/summary.json`, derives the public claim label from the visible evidence and stated validation bar, and currently prints `preliminary` plus the missing-evidence list. It does not rely on private governance documents or a hardcoded validation constant.

The current public-safe scope is preliminary: same-data shared-metric recovery is visible, but this is not yet an improvement claim. The planted-GT robustness upgrade remains locked until META verifies the checkpoint and the segmentation-corruption sweep is run.

## Citations and references

See `BASELINE_REFERENCES.md` for papers, code references, and citation context.
