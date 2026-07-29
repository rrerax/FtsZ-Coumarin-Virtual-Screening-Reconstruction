# FtsZ Coumarin Virtual Screening Reconstruction

![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![RDKit](https://img.shields.io/badge/RDKit-SMILES%20standardization-green)
![AutoDock Vina](https://img.shields.io/badge/AutoDock%20Vina-1.2.7-orange)
![Workflow](https://img.shields.io/badge/workflow-reconstruction%20and%20rerun-purple)

A reproducible reconstruction of a legacy high-throughput virtual screening workflow for **coumarin-like compounds against FtsZ binding pockets**. The repository converts spreadsheet-based screening outputs into curated molecular tables, reconstructed pocket boxes, rerun docking results, residue-contact summaries, receptor-ensemble docking outputs, and EnOpt-style reranking figures.

The current analysis is limited to computational docking reconstruction and candidate prioritization. Docking scores and poses are not interpreted as experimental activity measurements.

For file-by-file navigation guides, see `PROJECT_INDEX.md` and `PROJECT_INDEX_CN.md`.

## Summary

| Output | Result |
|---|---:|
| Method-level records parsed | 204 |
| Valid SMILES records | 204 / 204 |
| Unique pocket-compound records | 127 |
| Unique desalted structures | 90 |
| Ensemble-ready matrix rows | 88 |
| Unique-structure Vina reruns | 90 / 90 completed |
| Receptor conformations in ensemble | 5 |
| Multi-conformation docking jobs | 450 / 450 completed |

![Workflow summary](results/figures/workflow_summary.png)

## Methods

- **Data reconstruction**: parses legacy BP1/BP2 screening workbooks into normalized candidate tables.
- **Ligand standardization**: validates SMILES, removes salts, canonicalizes structures, merges duplicate entries, and computes RDKit descriptors.
- **Pocket reconstruction**: maps BP1 and BP2 to source pocket-detection outputs using pocket centers, volumes, and residue anchors.
- **Docking rerun**: rebuilds receptor/ligand inputs and reruns AutoDock Vina for all 90 cleaned unique structures across BP1/BP2.
- **Contact analysis**: summarizes receptor residues within 4 Å of the top rerun pose for each pocket.
- **EnOpt-style extension**: generates a small normal-mode receptor ensemble, runs multi-conformation docking, and exports conformation-weighted reranking outputs.

## Workflow

```mermaid
flowchart LR
    A["Legacy screening files"] --> B["Parse BP1/BP2 tables"]
    B --> C["Validate and canonicalize SMILES"]
    C --> D["Deduplicate structures"]
    D --> E["Reconstruct docking boxes"]
    E --> F["Prepare receptor and ligands"]
    F --> G["AutoDock Vina rerun"]
    G --> H["Score comparison and contact analysis"]
    F --> I["Build receptor ensemble"]
    I --> J["Multi-conformation docking"]
    J --> K["EnOpt-style weighted reranking"]
    D --> L["BP1/BP2 baseline score matrix"]
```

## Key Figures

| Score distribution | Top candidate structures |
|---|---|
| ![Score distribution](results/figures/score_distribution_by_pocket.png) | ![Top ligand structures](results/figures/top_ligand_structures.png) |

| Rerun score comparison | Pocket overlap |
|---|---|
| ![Rerun vs original](results/figures/vina_rerun_vs_original_scores.png) | ![BP overlap](results/figures/bp_score_scatter_overlap.png) |

| Binding-mode contact views |
|---|
| ![Binding-mode contact views](results/figures/binding_mode_top_hits.png) |

| EnOpt-style reranking | Conformation score matrix |
|---|---|
| ![EnOpt-style reranking](results/figures/enopt_style_reranking_vs_legacy.png) | ![Conformation score matrix](results/figures/ensemble_score_matrix_heatmap.png) |

## Repository Contents

| Path | Description |
|---|---|
| `PROJECT_INDEX.md` | File-by-file navigation guide |
| `PROJECT_INDEX_CN.md` | Chinese file-by-file navigation guide |
| `docs/technical_report.md` | Technical report with methodology, results, caveats, and rerun commands |
| `docs/technical_report_cn.md` | Chinese technical report covering the same reconstruction and limitations |
| `docs/enopt_extension.md` | EnOpt-style ensemble extension note |
| `docs/data_availability.md` | Data-availability notes and raw-input layout for full reruns |
| `notebooks/ftsZ_reconstruction.ipynb` | Notebook companion for reviewing tables and figures |
| `data/processed/bp_unique_structures.csv` | Structure-deduplicated BP1/BP2 candidate table |
| `data/processed/enopt_score_matrix.csv` | Two-pocket ensemble-ready score matrix |
| `results/tables/vina_rerun_unique_structures.csv` | Legacy-vs-rerun docking score comparison for all cleaned unique structures |
| `results/tables/binding_contacts_top_hits.csv` | 4 Å residue-contact summary for top rerun hits |
| `results/tables/ensemble_vina_scores.csv` | Multi-conformation docking score table |
| `results/tables/enopt_style_score_matrix.csv` | Conformation-level EnOpt-style score matrix |
| `results/tables/enopt_style_weights.csv` | Conformation weights used for weighted reranking |
| `results/tables/enopt_style_top_hits.csv` | Top hits by ensemble and EnOpt-style metrics |
| `results/figures/binding_mode_top_hits.png` | BP1/BP2 binding-mode contact visualization |
| `results/figures/enopt_style_reranking_vs_legacy.png` | Weighted ensemble score comparison |
| `results/figures/ensemble_score_matrix_heatmap.png` | Compound × receptor-conformation score heatmap |
| `results/pymol/` | PyMOL-ready rendering scripts for locally restored docking poses |

## Reproducibility

The repository includes processed tables, final figures, scripts, and reports. Original spreadsheets, receptor structures, pocket-detection exports, article PDFs, and docking intermediates are managed outside version control; see `docs/data_availability.md` for details.

Review the processed outputs and regenerate figures:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m compileall -q scripts
python scripts/draw_top_ligands.py
python scripts/analyze_enopt_style.py
python scripts/make_summary_figure.py
```

Run the complete workflow after restoring raw inputs under `data/raw/drive/`:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
micromamba create -y -p .mamba_vina -f environment-vina.yml

python scripts/clean_candidates.py
python scripts/derive_pocket_boxes.py
python scripts/prepare_docking_inputs.py --top-n 0 --max-heavy-atoms 100
python scripts/run_vina_batch.py --exhaustiveness 4 --cpu 4
python scripts/analyze_binding_contacts.py
python scripts/build_receptor_ensemble.py --target-ca-rmsd 0.65 --modes 2
python scripts/run_ensemble_vina.py --exhaustiveness 3 --cpu 4
python scripts/analyze_enopt_style.py
python scripts/draw_top_ligands.py
python scripts/draw_binding_mode_views.py
python scripts/make_summary_figure.py
```

AutoDock Vina and OpenBabel are isolated in `.mamba_vina` via `environment-vina.yml`; Python analysis dependencies are listed in `requirements.txt`.

## Project Structure

```text
.
├── data/processed/              # Curated derived tables
├── docs/                        # Technical reports and data-availability notes
├── notebooks/                   # Notebook companion
├── results/figures/             # Generated analysis figures
├── results/pymol/               # PyMOL-ready local rendering scripts
├── results/tables/              # Summary and rerun-result tables
├── scripts/                     # Cleaning, docking, analysis, and plotting scripts
├── environment-vina.yml         # Vina/OpenBabel conda environment
└── requirements.txt             # Python analysis dependencies
```

## Methodological Notes

- BP1 and BP2 are treated as distinct binding pockets, not conformers of one binding site.
- Docking boxes were reconstructed from pocket-detection outputs and residue anchors because the original docking configuration was not available.
- The current EnOpt-style extension uses conformation-weighted consensus reranking without experimental active/decoy labels; it is not a supervised EnOpt model.
- Docking scores and contact views are interpreted only as computational prioritization evidence rather than biochemical potency measurements.

## Future Work

- Recover the original docking configuration if available and compare against the reconstructed boxes.
- Replace the normal-mode pilot ensemble with experimentally resolved conformers or MD snapshots.
- Add active/decoy labels or a documented benchmark set before supervised EnOpt training.
- Refine selected binding modes with PyMOL ray-traced figures or 2D interaction diagrams.
