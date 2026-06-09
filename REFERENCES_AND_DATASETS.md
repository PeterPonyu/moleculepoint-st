# MoleculePoint-ST — references (with code) & datasets

Consolidated reference + dataset index. Paper DOIs verified via Crossref and code
repositories via the GitHub API on 2026-06-09. See `BASELINE_REFERENCES.md` for the
full provenance and audit boundary.

## Reference papers & method baselines (with public code)

| Role | Method | Venue / year | DOI | Code |
|------|--------|--------------|-----|------|
| Primary | smoppix — unified nonparametric analysis of single-molecule spatial omics via probabilistic indices | Genome Biology 2026 | `10.1186/s13059-026-03976-5` | https://github.com/sthawinke/smoppix (+ Bioconductor) |
| Companion | SmoppixPaper — paper analysis repo | — | — | https://github.com/maerelab/SmoppixPaper |
| Baseline | Bento — subcellular ST analysis toolkit | — | — | https://github.com/ckmah/bento-tools |
| Baseline | starfish — image-based transcriptomics pipelines | — | — | https://github.com/spacetx/starfish |
| Baseline | spicyR — spatial interaction analysis | — | — | https://github.com/SydneyBioX/spicyR |

## Datasets

Runs on user-supplied single-molecule / subcellular imaging-ST inputs (no shipped catalog).
Reference data comes from the smoppix example data and the baseline toolkits above.

> Verification: smoppix DOI confirmed in Crossref; smoppix / starfish / spicyR repos
> confirmed live via GitHub API. bento-tools is real (BSD-2; API returned a transient 403).
