# Data Audit — FtsZ Coumarin HTVS Reconstruction

## Inputs
- `BP1_candidates.xlsx`: thesis spreadsheet labeled as pocket 1 / BP1.
- `BP2_candidates.xlsx`: thesis spreadsheet labeled as pocket 2 / BP2.
- `Receptor_from_fconv.pdb`: FtsZ receptor copied from the ezPocket result folder.
- `ezPocket_fconv.tsv`: two pocket centers and volumes from the source project folder.

## Cleaning rule
The first spreadsheet block is interpreted as: compound ID, docking score, Tanimoto similarity, SMILES, and fingerprint method. The spreadsheet column labels are inherited from the old workbook and are therefore not used literally.

## Counts
| pocket   |   method_records |   valid_smiles |   unique_compounds |
|:---------|-----------------:|---------------:|-------------------:|
| BP1      |              101 |            101 |                 55 |
| BP2      |              103 |            103 |                 72 |

## EnOpt readiness note
The cleaned score matrix is EnOpt-compatible at the data-shape level, but the two score columns currently represent two binding pockets, not two conformations of the same pocket. For a methodologically strict EnOpt run, add multiple receptor conformations for the same binding site and experimental active/inactive labels or a documented surrogate label.

## Key outputs
- `data/processed/bp_method_records.csv`
- `data/processed/bp_unique_compounds.csv`
- `data/processed/bp_unique_structures.csv`
- `data/processed/enopt_score_matrix.csv`
- `results/tables/top10_by_pocket.csv`
- `results/tables/top10_unique_structures_by_pocket.csv`
- `results/figures/score_distribution_by_pocket.png`
- `results/figures/top10_scores_by_pocket.png`
- `results/figures/bp_score_scatter_overlap.png`
- `results/figures/top_ligand_structures.png`
- `results/figures/binding_mode_top_hits.png`
- `results/pymol/`: PyMOL-ready binding-view scripts
