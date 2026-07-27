# Project Index

This index summarizes where to find the main materials in the FtsZ coumarin virtual-screening reconstruction repository.

## Recommended reading order

| Step | File | Purpose |
|---:|---|---|
| 1 | `README.md` | Project overview, workflow, key figures, and reproducibility commands |
| 2 | `docs/technical_report_cn.md` | Chinese technical report with methodology, results, and limitations |
| 3 | `docs/technical_report.md` | English technical report with detailed tables and reproducibility notes |
| 4 | `notebooks/ftsZ_reconstruction.ipynb` | Notebook companion for reviewing processed tables and figures |
| 5 | `docs/data_availability.md` | Raw-input layout for complete local reruns |

## Main results

| Output | Path |
|---|---|
| Workflow summary figure | `results/figures/workflow_summary.png` |
| Score distribution by pocket | `results/figures/score_distribution_by_pocket.png` |
| Top candidate structure panel | `results/figures/top_ligand_structures.png` |
| Legacy-vs-rerun Vina score comparison | `results/figures/vina_rerun_vs_original_scores.png` |
| BP1/BP2 overlap plot | `results/figures/bp_score_scatter_overlap.png` |
| Binding-mode contact views | `results/figures/binding_mode_top_hits.png` |
| Vina rerun table | `results/tables/vina_rerun_top_structures.csv` |
| Residue-contact summary | `results/tables/binding_contacts_top_hits.csv` |
| Ensemble-ready score matrix | `data/processed/enopt_score_matrix.csv` |

## Folder map

| Folder | Contents |
|---|---|
| `data/processed/` | Curated candidate tables, pocket boxes, and score matrices |
| `docs/` | Technical reports, data audit, and data-availability notes |
| `notebooks/` | Review notebook for tables and figures |
| `results/figures/` | Generated figures used in README and reports |
| `results/tables/` | CSV outputs from cleaning, rerun, correlation, and contact analyses |
| `results/pymol/` | PyMOL-ready scripts for locally restored receptor/pose files |
| `scripts/` | Reproducible Python scripts for cleaning, docking preparation, reruns, contacts, and figures |

## Reproducibility notes

The repository contains processed outputs and generated figures. A complete raw-to-results rerun requires restoring the raw files listed in `docs/data_availability.md` under `data/raw/drive/`.

Quick figure regeneration:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/draw_top_ligands.py
python scripts/make_summary_figure.py
```

Full local rerun after restoring raw inputs:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
micromamba create -y -p .mamba_vina -f environment-vina.yml

python scripts/clean_candidates.py
python scripts/derive_pocket_boxes.py
python scripts/prepare_docking_inputs.py --top-n 10 --drug-like-only --max-heavy-atoms 35
python scripts/run_vina_batch.py --exhaustiveness 4 --cpu 4
python scripts/analyze_binding_contacts.py
python scripts/draw_top_ligands.py
python scripts/draw_binding_mode_views.py
python scripts/make_summary_figure.py
```
