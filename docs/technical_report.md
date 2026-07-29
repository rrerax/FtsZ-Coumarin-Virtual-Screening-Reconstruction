# Technical Report: FtsZ Coumarin HTVS Reconstruction and Rerun

## Executive summary

This stage reconstructs a legacy FtsZ/coumarin virtual-screening workflow into a reproducible project structure. The BP1/BP2 candidate spreadsheets were cleaned, SMILES strings were validated and canonicalized, two FtsZ pocket boxes were reconstructed from source ezPocket outputs and residue anchors, and an AutoDock Vina rerun was completed for all 90 cleaned unique structures across BP1/BP2. A small normal-mode receptor ensemble was then added to test EnOpt-style conformation-weighted reranking mechanics.

Key numbers:

| pocket   |   method_records |   valid_method_records |   unique_ids |   min_score |   median_score |   max_score |
|:---------|-----------------:|-----------------------:|-------------:|------------:|---------------:|------------:|
| BP1      |              101 |                    101 |           55 |      -7.396 |         -6.98  |      -6.318 |
| BP2      |              103 |                    103 |           72 |      -9.116 |         -7.429 |      -5.969 |

Additional structure-level deduplication produced **90 unique pocket-structure records**, an **88-row BP1/BP2 score matrix**, a **5-conformation receptor ensemble**, and **450/450 completed multi-conformation docking jobs** for the full cleaned unique-structure rerun set.

## Source inputs used

| Input | Local path | Role |
|---|---|---|
| Thesis/project notes | managed locally | Original project narrative, methods, BP1/BP2 interpretation |
| EnOpt reference paper | managed locally | Methodological reference for ensemble/ML reranking |
| BP1 workbook | `data/raw/drive/BP1_candidates.xlsx` | Source table for pocket 1 |
| BP2 workbook | `data/raw/drive/BP2_candidates.xlsx` | Source table for pocket 2 |
| Receptor PDB | `data/raw/drive/Receptor_from_fconv.pdb` | Receptor source used for local docking preparation |
| Pocket outputs | `data/raw/drive/ezPocket_*.tsv` | Pocket-center and volume source files |

## Data cleaning

The old workbooks preserve legacy column labels, so the first spreadsheet block was interpreted by content rather than by header text:

- column 1: compound/source ID
- column 2: original docking score, kcal/mol
- column 3: Tanimoto similarity to the coumarin reference query
- column 4: SMILES
- column 5: fingerprint method

RDKit was used to validate SMILES, canonicalize structures, remove salts by keeping the largest fragment, and calculate light drug-likeness descriptors.

Primary cleaned outputs:

- `data/processed/bp_method_records.csv`
- `data/processed/bp_unique_compounds.csv`
- `data/processed/bp_unique_structures.csv`
- `data/processed/enopt_score_matrix.csv`

## Pocket reconstruction

The two binding pockets were mapped to source pocket-detection outputs using the residue anchors described in the original project notes.

| pocket   | source                          |   center_x |   center_y |   center_z |   size_x |   size_y |   size_z | thesis_residue_anchor                                         |
|:---------|:--------------------------------|-----------:|-----------:|-----------:|---------:|---------:|---------:|:--------------------------------------------------------------|
| BP1      | ezPocket_fpocket3.tsv / Pocket1 |   -50.6447 |    35.0742 |   15.2025  |     22.5 |     22.5 |     22.5 | Ala39B, Asn41B, Leu47B, Met49B, Ser50B, Lys55B                |
| BP2      | ezPocket_fconv.tsv / Pocket0    |   -63.4691 |    20.8178 |    1.73489 |     22.5 |     22.5 |     22.5 | Met169B, Glu185B, Asn189B, Ile225B, Gly226B, Ser227B, Arg304B |

The Vina box sizes are estimated from source pocket centers/volumes and residue anchors. If the original ezSMDock/Vina box configuration is recovered later, these values should be replaced and the rerun can be repeated without changing the rest of the pipeline.

## Original candidate ranking after structure deduplication

### BP1 top structures

| compound_ids            | source_databases          |   best_docking_score_kcal_mol |   max_tanimoto_to_reference |   mol_weight |   logp | passes_lipinski_light   |
|:------------------------|:--------------------------|------------------------------:|----------------------------:|-------------:|-------:|:------------------------|
| 251157;ZINC000103969950 | ZINC15;numeric library ID |                        -7.396 |                       0.88  |      226.231 |  2.96  | True                    |
| CHEMBL1929523           | ChEMBL23                  |                        -7.388 |                       0.6   |      236.226 |  3.692 | True                    |
| CHEMBL3306562           | ChEMBL23                  |                        -7.338 |                       0.75  |      237.214 |  2.523 | True                    |
| DB03623                 | DrugBank                  |                        -7.328 |                       0.316 |      272.307 |  4.156 | True                    |
| MolPort-002-367-171     | MolPort                   |                        -7.327 |                       0.846 |      190.198 |  2.115 | True                    |
| CHEMBL329500            | ChEMBL23                  |                        -7.27  |                       0.913 |      288.302 |  4.319 | True                    |
| PEY                     | PDB ligand                |                        -7.242 |                       0.659 |      178.234 |  3.993 | True                    |
| HMDB34138               | HMDB                      |                        -7.241 |                       0.786 |      284.223 |  2.809 | True                    |
| 19824                   | numeric library ID        |                        -7.229 |                       0.282 |      275.351 |  4.226 | True                    |
| CHEMBL3409185           | ChEMBL23                  |                        -7.206 |                       0.571 |      321.126 |  3.979 | True                    |

### BP2 top structures

| compound_ids                           | source_databases                   |   best_docking_score_kcal_mol |   max_tanimoto_to_reference |   mol_weight |   logp | passes_lipinski_light   |
|:---------------------------------------|:-----------------------------------|------------------------------:|----------------------------:|-------------:|-------:|:------------------------|
| 746361                                 | numeric library ID                 |                        -9.116 |                       0.441 |      344.326 |  4.714 | True                    |
| 3543                                   | numeric library ID                 |                        -9.093 |                       0.265 |      380.451 |  7.993 | False                   |
| 118810                                 | numeric library ID                 |                        -8.495 |                       0.694 |      282.295 |  3.386 | True                    |
| 50651                                  | numeric library ID                 |                        -8.307 |                       0.241 |      336.778 |  4.932 | True                    |
| MolPort-002-321-924                    | MolPort                            |                        -8.303 |                       0.562 |      408.322 |  1.385 | True                    |
| 686393;CHEMBL392451                    | ChEMBL23;numeric library ID        |                        -8.211 |                       0.633 |      746.68  |  7.295 | False                   |
| CHEMBL1364708;MolPort-003-378-434      | ChEMBL23;MolPort                   |                        -8.175 |                       0.475 |      278.263 |  3.912 | True                    |
| 4368550;CHEMBL1801014;ZINC000000057885 | ChEMBL23;ZINC15;numeric library ID |                        -8.141 |                       0.913 |      288.302 |  4.319 | True                    |
| 7524                                   | numeric library ID                 |                        -8.124 |                       0.315 |      673.8   |  1.356 | False                   |
| 656562                                 | numeric library ID                 |                        -8.08  |                       0.818 |      296.325 |  5.253 | False                   |

## AutoDock Vina rerun

A full rerun was completed for all 90 cleaned unique structures from BP1 and BP2. Ligand 3D coordinates were regenerated from the cleaned SMILES strings, receptor and ligand PDBQT files were prepared with an OpenBabel fallback, and Vina 1.2.7 was run with exhaustiveness 4 and seed 20260724. One large BP2 ligand required a loose RDKit 3D-coordinate fallback, but the full cleaned set completed successfully.

Top rerun hits by pocket:

| pocket   |   rerun_count |   successful_jobs | top_rerun_compound   |   legacy_score_kcal_mol |   rerun_score_kcal_mol |   score_delta |
|:---------|--------------:|------------------:|:---------------------|------------------------:|-----------------------:|--------------:|
| BP1      |            37 |                37 | DB00776              |                  -6.883 |                  -8.02 |        -1.137 |
| BP2      |            53 |                53 | 686393;CHEMBL392451  |                  -8.211 |                 -10.76 |        -2.549 |

Correlation between legacy worksheet scores and rerun Vina scores:

| scope   |   n |   pearson_r |   pearson_p |   spearman_r |   spearman_p |
|:--------|----:|------------:|------------:|-------------:|-------------:|
| BP1     |  37 |       0.528 |      0.0008 |        0.543 |       0.0005 |
| BP2     |  53 |       0.743 |      0      |        0.771 |       0      |
| overall |  90 |       0.737 |      0      |        0.737 |       0      |

Figures:

![Score distribution](../results/figures/score_distribution_by_pocket.png)

![Top 10 scores](../results/figures/top10_scores_by_pocket.png)

![Rerun vs original](../results/figures/vina_rerun_vs_original_scores.png)


## Binding-pose contact summary

For the best full-rerun hit from each pocket, receptor residues within 4 Å of the top Vina pose were summarized as a lightweight binding-mode check. This does not replace manual PyMOL/Discovery Studio inspection, but it confirms that rerun poses are placed around residues discussed in the thesis.

| pocket   | compound_ids        | residue   |   min_distance_a |   contact_count_within_4a | has_polar_candidate_within_3_5a   |
|:---------|:--------------------|:----------|-----------------:|--------------------------:|:----------------------------------|
| BP1      | DB00776             | LEU47B    |            3.197 |                        14 | False                             |
| BP1      | DB00776             | LYS55B    |            3.228 |                        13 | True                              |
| BP1      | DB00776             | ASP57B    |            3.249 |                         9 | True                              |
| BP1      | DB00776             | SER50B    |            3.285 |                         4 | True                              |
| BP1      | DB00776             | LEU56B    |            3.383 |                         7 | True                              |
| BP1      | DB00776             | MET49B    |            3.435 |                         4 | False                             |
| BP1      | DB00776             | ALA39B    |            3.58  |                         3 | False                             |
| BP1      | DB00776             | ASN41B    |            3.607 |                         3 | False                             |
| BP1      | DB00776             | LEU48B    |            3.612 |                         2 | False                             |
| BP2      | 686393;CHEMBL392451 | GLU185B   |            2.881 |                        11 | True                              |
| BP2      | 686393;CHEMBL392451 | SER182B   |            3.083 |                        15 | True                              |
| BP2      | 686393;CHEMBL392451 | ARG304B   |            3.103 |                        11 | True                              |
| BP2      | 686393;CHEMBL392451 | MET169B   |            3.141 |                        10 | False                             |
| BP2      | 686393;CHEMBL392451 | VAL186B   |            3.264 |                         7 | True                              |
| BP2      | 686393;CHEMBL392451 | LEU188B   |            3.361 |                         1 | False                             |
| BP2      | 686393;CHEMBL392451 | SER227B   |            3.463 |                         8 | False                             |
| BP2      | 686393;CHEMBL392451 | ASN189B   |            3.474 |                         6 | False                             |
| BP2      | 686393;CHEMBL392451 | ILE225B   |            3.483 |                         3 | False                             |
| BP2      | 686393;CHEMBL392451 | PRO245B   |            3.491 |                         1 | False                             |
| BP2      | 686393;CHEMBL392451 | GLY170B   |            3.591 |                        11 | False                             |
| BP2      | 686393;CHEMBL392451 | ILE240B   |            3.591 |                         6 | False                             |
| BP2      | 686393;CHEMBL392451 | GLN192B   |            3.69  |                         1 | False                             |
| BP2      | 686393;CHEMBL392451 | GLY226B   |            3.788 |                         1 | False                             |
| BP2      | 686393;CHEMBL392451 | VAL174B   |            3.824 |                         1 | False                             |
| BP2      | 686393;CHEMBL392451 | SER244B   |            3.986 |                         1 | False                             |

Additional output:

- `results/tables/binding_contacts_top_hits.csv`
- `results/figures/top_ligand_structures.png`
- `results/figures/binding_mode_top_hits.png`
- `results/pymol/` with PyMOL-ready rendering scripts for locally restored docking poses

![Top ligand structures](../results/figures/top_ligand_structures.png)

![Binding-mode contact views](../results/figures/binding_mode_top_hits.png)

## EnOpt-style ensemble update

The cleaned BP1/BP2 matrix remains useful as a pocket-level screening baseline, while the added EnOpt-style module now tests conformation-level reranking on the same FtsZ system. Five receptor conformations were generated from the reconstructed FtsZ receptor using ANM normal-mode perturbations: the original receptor plus positive and negative perturbations along the first two non-trivial modes. All 90 cleaned unique structures were redocked across these five conformations, producing 450 completed Vina jobs.

Core outputs:

- `data/processed/enopt_score_matrix.csv` for the original BP1/BP2 pocket-level baseline
- `results/tables/receptor_ensemble_manifest.csv` for generated receptor conformations
- `results/tables/ensemble_vina_scores.csv` for conformation-level Vina scores
- `results/tables/enopt_style_score_matrix.csv` for compound × conformation features
- `results/tables/enopt_style_weights.csv` for rank-consistency conformation weights
- `results/tables/enopt_style_top_hits.csv` for top candidates under ensemble-best, ensemble-mean, and weighted metrics

Top EnOpt-style weighted results from the pilot ensemble:

| pocket | top weighted compound | weighted score, kcal/mol | legacy score, kcal/mol |
|---|---:|---:|---:|
| BP1 | CHEMBL329500 | -7.691 | -7.270 |
| BP2 | 85Z | -9.702 | -7.594 |

Methodological note for CADD review: this is an **EnOpt-style conformation-weighted consensus reranking**, not a strict supervised EnOpt model. The source files do not contain experimental active/decoy labels, so the current weighting is an exploratory rank-consistency procedure rather than activity-trained machine learning.

![EnOpt-style reranking](../results/figures/enopt_style_reranking_vs_legacy.png)

![Conformation-level score matrix](../results/figures/ensemble_score_matrix_heatmap.png)

## Reproducibility commands

The repository includes processed tables and final figures. A full raw-to-docking rerun requires restoring the raw inputs listed in `docs/data_availability.md`.

Public review / figure regeneration:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m compileall -q scripts
python scripts/draw_top_ligands.py
python scripts/analyze_enopt_style.py
python scripts/make_summary_figure.py
```

Full local rerun with raw source files:

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
python scripts/draw_top_ligands.py
python scripts/draw_binding_mode_views.py
python scripts/make_summary_figure.py
python scripts/build_receptor_ensemble.py --target-ca-rmsd 0.65 --modes 2
python scripts/run_ensemble_vina.py --exhaustiveness 3 --cpu 4
python scripts/analyze_enopt_style.py
```

Vina and OpenBabel were installed in `.mamba_vina` from conda-forge for the local rerun.

## Limitations and next steps

- The exact original docking box file was not available, so box sizes are reconstructed estimates.
- One large BP2 candidate required a loose RDKit 3D-coordinate fallback during ligand preparation; the full cleaned set nevertheless completed successfully.
- The receptor ensemble is a normal-mode pilot around one reconstructed structure, not a curated experimental-structure or MD-snapshot ensemble.
- The source files lack experimental active/inactive labels; strict supervised EnOpt model training is outside the scope of this stage.
- Docking score changes should be interpreted as rerun consistency evidence, not biological validation.
- The next rigorous extension is to replace the pilot ensemble with curated FtsZ conformers or MD snapshots and evaluate the reranking against active/decoy benchmarks.
