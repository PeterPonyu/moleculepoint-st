# MoleculePoint-ST baseline references

Verification date: 2026-06-08

## Baseline decision summary

| Role | Baseline | Decision |
|---|---|---|
| Primary | smoppix | Use as the first open-code/public-artifact starting point because it directly matches this track's input-output problem. |
| Secondary | See table below | Use only for comparison, adapter design, and ablation inspiration; do not copy implementation. Bento has a same-data isolated shared-metric run artifact (segmentation-based Ripley vector); smoppix remains reference-reported pending an R/Bioconductor environment. |

## Primary baseline

- Paper title: smoppix: unified nonparametric analysis of single-molecule spatial omics data using probabilistic indices
- Venue/date: Genome Biology, 2026
- Article URL: https://link.springer.com/article/10.1186/s13059-026-03976-5
- Code/artifact URL: https://github.com/sthawinke/smoppix and https://git.bioconductor.org/packages/smoppix
- Verification date: 2026-06-08
- Default branch/artifact: main
- Observed HEAD SHA or DOI: 2a318bceacb32138e8c9e80fa05c0480be2b353a
- Local audit checkout/artifact: `baselines/smoppix-original`
- License note: GPL-2 in DESCRIPTION; Bioconductor release also reports GPL-2
- Local use: Primary public-code reference for aggregation, colocalization, gradients, and vicinity tests on single-molecule spatial omics data. Same-data run status: `REFERENCE_REPORTED`; a dedicated R/Bioconductor environment is needed for a true run.
- Fallback: If this public code/artifact becomes unavailable, mark this track `deferred-unverified` until a comparable open-code baseline is found.
- Verification command/evidence:
  - `git ls-remote --symref <repo> HEAD` for GitHub/Git/Bioconductor repositories.
  - `zenodo.org/api/records/<record>` plus local file checks for Zenodo artifacts.
  - Same-data R run not attempted in this gate because the package is GPL-2 and requires a dedicated R/Bioconductor stack.

## Secondary verified references

| Baseline | Role | Code/artifact URL | Local audit checkout | Branch/artifact | Observed SHA/DOI | License note |
|---|---|---|---|---|---:|---|
| SmoppixPaper | smoppix paper analysis repository | https://github.com/maerelab/SmoppixPaper | `not cloned yet` | `main` | `876642aab205` | GPL-3.0 reported by GitHub API |
| Bento | Subcellular spatial transcriptomics analysis toolkit | https://github.com/ckmah/bento-tools | `baselines/bento-tools-original` | `master` | `3e82209c19f8` | BSD-2-Clause license file observed locally; same-data isolated Python shared-metric artifacts written to `experiments/gate2_baseline_comparison/reference_rows.json` and `experiments/gate2_baseline_comparison/shared_metric_reference_vector.json` |
| starfish | Image-based transcriptomics pipelines | https://github.com/spacetx/starfish | `not cloned yet` | `master` | `5c86bdaa024c` | MIT reported by GitHub API |
| spicyR | Spatial interaction analysis reference | https://github.com/SydneyBioX/spicyR | `not cloned yet` | `devel` | `64675c3eadc6` | license requires re-check |

## Brand independence note

Reference names in this file are provenance labels only. Local package names, CLI commands, figure labels, and manuscript novelty claims must use `MoleculePoint-ST` terminology and the independent refinements in `README.md`, not upstream branding.
