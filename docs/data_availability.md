# Data Availability

This repository contains curated tables, analysis scripts, figures, and reports derived from a legacy FtsZ virtual-screening workflow. Raw input files are managed separately to avoid redistributing unpublished spreadsheets, structure files, generated docking intermediates, and copyrighted PDFs.

## Files included here

| Path | Content |
|---|---|
| `data/processed/` | Curated candidate tables and score matrices |
| `results/tables/` | Summary tables, rerun-score comparisons, contact summaries, and EnOpt-style ensemble outputs |
| `results/figures/` | Generated plots, ligand panels, and contact-view figures |
| `results/pymol/` | PyMOL-ready scripts for rendering locally restored docking poses |
| `scripts/` | Cleaning, docking-preparation, rerun, analysis, and plotting scripts |
| `notebooks/` | Notebook companion for table and figure review |
| `docs/` | Technical reports and data-availability notes |

## Raw inputs managed locally

| Path / type | Purpose |
|---|---|
| `data/raw/` | Original candidate spreadsheets, receptor structures, and pocket-detection exports |
| `references/` | Source thesis/literature PDFs |
| `work/` | Generated PDBQT files, receptor-ensemble PDBQT files, Vina logs, poses, and rerun intermediates |
| `*.xlsx`, `*.pdf`, `*.pdb`, `*.mol2`, `*.sdf` | Raw spreadsheets, articles, receptor files, pocket files, and structure files |
| `.venv/`, `.mamba_vina/` | Local dependency environments |

## Expected raw input layout for full reruns

A complete raw-to-results rerun expects the following files under `data/raw/drive/`:

| Expected file | Purpose |
|---|---|
| `BP1_candidates.xlsx` | Candidate table for binding pocket 1 |
| `BP2_candidates.xlsx` | Candidate table for binding pocket 2 |
| `Receptor_from_fconv.pdb` | FtsZ receptor structure used for docking preparation |
| `ezPocket_fconv.tsv` | Pocket-detection output used for BP2 box reconstruction |
| `ezPocket_fpocket2.tsv` | Supporting pocket-detection output |
| `ezPocket_fpocket3.tsv` | Pocket-detection output used for BP1 box reconstruction |

The repository can be reviewed from the processed outputs alone. Full docking regeneration requires restoring the raw inputs in the layout above.

## Notes

- Docking poses and scores are computational screening outputs, not biochemical activity measurements.
- The EnOpt-style outputs use a normal-mode receptor ensemble and rank-consistency weighting; they are not a supervised activity-prediction model without active/decoy labels.
- Raw spreadsheets and structure files are intentionally kept out of version control.
- The current repository should be interpreted as a reproducible reconstruction and rerun workflow, not as a complete experimental validation package.
