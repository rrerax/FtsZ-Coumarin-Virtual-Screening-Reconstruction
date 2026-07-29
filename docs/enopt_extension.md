# EnOpt-style Ensemble Extension

This extension adds a small receptor-conformation ensemble to the FtsZ HTVS reconstruction. The current implementation uses normal-mode perturbations around the reconstructed FtsZ receptor and redocks all 90 cleaned BP1/BP2 unique structures across the ensemble.

The output is described as EnOpt-style because no experimental active/decoy labels are available in the source material. The weighted score is therefore an exploratory conformation-weighted consensus score, not a supervised EnOpt model.

## Outputs

- `results/tables/receptor_ensemble_manifest.csv`
- `results/tables/ensemble_vina_scores.csv`
- `results/tables/enopt_style_score_matrix.csv`
- `results/tables/enopt_style_weights.csv`
- `results/tables/enopt_style_top_hits.csv`
- `results/tables/enopt_style_correlations.csv`
- `results/figures/enopt_style_reranking_vs_legacy.png`
- `results/figures/ensemble_score_matrix_heatmap.png`
- `results/figures/enopt_style_conformation_weights.png`

## Top hits

| pocket   | metric                              | top_compound_ids    | top_ligand_name   |   top_score_kcal_mol |   legacy_score_kcal_mol |
|:---------|:------------------------------------|:--------------------|:------------------|---------------------:|------------------------:|
| BP1      | ensemble_best_score_kcal_mol        | DB00776             | BP1_DB00776       |               -7.797 |                  -6.883 |
| BP1      | ensemble_mean_score_kcal_mol        | CHEMBL329500        | BP1_CHEMBL329500  |               -7.691 |                  -7.27  |
| BP1      | enopt_style_weighted_score_kcal_mol | CHEMBL329500        | BP1_CHEMBL329500  |               -7.691 |                  -7.27  |
| BP2      | ensemble_best_score_kcal_mol        | 686393;CHEMBL392451 | BP2_686393        |              -10.46  |                  -8.211 |
| BP2      | ensemble_mean_score_kcal_mol        | 85Z                 | BP2_85Z           |               -9.691 |                  -7.594 |
| BP2      | enopt_style_weighted_score_kcal_mol | 85Z                 | BP2_85Z           |               -9.702 |                  -7.594 |

## Correlation with legacy worksheet scores

| pocket   | metric                              |   n |   pearson_r |   pearson_p |   spearman_r |   spearman_p |
|:---------|:------------------------------------|----:|------------:|------------:|-------------:|-------------:|
| BP1      | ensemble_best_score_kcal_mol        |  37 |       0.519 |      0.001  |        0.492 |       0.002  |
| BP1      | ensemble_mean_score_kcal_mol        |  37 |       0.53  |      0.0007 |        0.527 |       0.0008 |
| BP1      | enopt_style_weighted_score_kcal_mol |  37 |       0.531 |      0.0007 |        0.527 |       0.0008 |
| BP2      | ensemble_best_score_kcal_mol        |  53 |       0.752 |      0      |        0.802 |       0      |
| BP2      | ensemble_mean_score_kcal_mol        |  53 |       0.747 |      0      |        0.802 |       0      |
| BP2      | enopt_style_weighted_score_kcal_mol |  53 |       0.752 |      0      |        0.804 |       0      |

## Method boundary

The current BP1/BP2 candidate set and normal-mode ensemble are suitable for testing the mechanics of ensemble reranking. A strict supervised EnOpt analysis would require experimental active/decoy labels or a benchmark set, plus a receptor ensemble designed for the same binding site.
