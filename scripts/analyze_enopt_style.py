#!/usr/bin/env python3
"""Create ensemble and EnOpt-style rankings from multi-conformation docking scores."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import softmax
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"


def minmax(series: pd.Series) -> pd.Series:
    if series.max() == series.min():
        return pd.Series(0.0, index=series.index)
    return (series - series.min()) / (series.max() - series.min())


def rank_inverse_weights(scores: pd.DataFrame) -> pd.Series:
    columns = [column for column in scores.columns if column.startswith("conf")]
    weight_rows = []
    for pocket, group in scores.groupby("pocket"):
        old = group["old_best_docking_score_kcal_mol"]
        correlations = {}
        for column in columns:
            valid = group[[column, "old_best_docking_score_kcal_mol"]].dropna()
            if len(valid) < 3 or valid[column].nunique() < 2:
                correlations[column] = 0.0
                continue
            corr = spearmanr(valid[column], valid["old_best_docking_score_kcal_mol"]).statistic
            correlations[column] = 0.0 if np.isnan(corr) else max(float(corr), 0.0)
        raw = pd.Series(correlations).fillna(0.0)
        weights = pd.Series(softmax(raw / 0.35), index=columns)
        for column, weight in weights.items():
            weight_rows.append({"pocket": pocket, "conformation_id": column, "weight": round(float(weight), 4)})
    return pd.DataFrame(weight_rows)


def build_score_matrix(raw: pd.DataFrame) -> pd.DataFrame:
    pivot = raw.pivot_table(
        index=["pocket", "ligand_name", "compound_ids", "source_databases", "desalted_smiles", "old_best_docking_score_kcal_mol"],
        columns="conformation_id",
        values="ensemble_vina_score_kcal_mol",
        aggfunc="min",
    ).reset_index()
    pivot.columns.name = None
    conf_columns = [column for column in pivot.columns if str(column).startswith("conf")]
    pivot["ensemble_best_score_kcal_mol"] = pivot[conf_columns].min(axis=1)
    pivot["ensemble_mean_score_kcal_mol"] = pivot[conf_columns].mean(axis=1)
    pivot["ensemble_std_score_kcal_mol"] = pivot[conf_columns].std(axis=1)
    pivot["ensemble_range_score_kcal_mol"] = pivot[conf_columns].max(axis=1) - pivot[conf_columns].min(axis=1)
    weights = rank_inverse_weights(pivot)
    weighted_rows = []
    for _, row in pivot.iterrows():
        pocket_weights = weights[weights["pocket"] == row["pocket"]].set_index("conformation_id")["weight"]
        available = [column for column in conf_columns if column in pocket_weights.index and pd.notna(row[column])]
        if not available:
            weighted_rows.append(np.nan)
            continue
        normalized = pocket_weights.loc[available] / pocket_weights.loc[available].sum()
        weighted_rows.append(float(sum(row[column] * normalized[column] for column in available)))
    pivot["enopt_style_weighted_score_kcal_mol"] = weighted_rows
    for metric in ["ensemble_best_score_kcal_mol", "ensemble_mean_score_kcal_mol", "enopt_style_weighted_score_kcal_mol"]:
        pivot[f"rank_{metric.replace('_score_kcal_mol', '')}"] = pivot.groupby("pocket")[metric].rank(method="min", ascending=True).astype(int)
    return pivot, weights


def summarize_rankings(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pocket, group in matrix.groupby("pocket"):
        for metric in ["ensemble_best_score_kcal_mol", "ensemble_mean_score_kcal_mol", "enopt_style_weighted_score_kcal_mol"]:
            top = group.sort_values(metric).head(1).iloc[0]
            rows.append(
                {
                    "pocket": pocket,
                    "metric": metric,
                    "top_compound_ids": top["compound_ids"],
                    "top_ligand_name": top["ligand_name"],
                    "top_score_kcal_mol": round(float(top[metric]), 3),
                    "legacy_score_kcal_mol": round(float(top["old_best_docking_score_kcal_mol"]), 3),
                }
            )
    return pd.DataFrame(rows)


def summarize_correlations(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["ensemble_best_score_kcal_mol", "ensemble_mean_score_kcal_mol", "enopt_style_weighted_score_kcal_mol"]
    for pocket, group in matrix.groupby("pocket"):
        for metric in metrics:
            valid = group[[metric, "old_best_docking_score_kcal_mol"]].dropna()
            if len(valid) < 3 or valid[metric].nunique() < 2:
                continue
            pearson = pearsonr(valid[metric], valid["old_best_docking_score_kcal_mol"])
            spearman = spearmanr(valid[metric], valid["old_best_docking_score_kcal_mol"])
            rows.append(
                {
                    "pocket": pocket,
                    "metric": metric,
                    "n": len(valid),
                    "pearson_r": round(float(pearson.statistic), 3),
                    "pearson_p": round(float(pearson.pvalue), 4),
                    "spearman_r": round(float(spearman.statistic), 3),
                    "spearman_p": round(float(spearman.pvalue), 4),
                }
            )
    return pd.DataFrame(rows)


def write_figures(matrix: pd.DataFrame, weights: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(8, 6))
    plot_data = matrix.copy()
    plot_data["label"] = plot_data["compound_ids"].str.split(";").str[0]
    sns.scatterplot(
        data=plot_data,
        x="old_best_docking_score_kcal_mol",
        y="enopt_style_weighted_score_kcal_mol",
        hue="pocket",
        s=90,
    )
    plt.xlabel("Legacy worksheet score (kcal/mol)")
    plt.ylabel("EnOpt-style weighted ensemble score (kcal/mol)")
    plt.title("Multi-conformation ensemble reranking")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "enopt_style_reranking_vs_legacy.png", dpi=220)
    plt.close()

    heatmap_data = matrix.set_index(["pocket", "compound_ids"])[[column for column in matrix.columns if str(column).startswith("conf")]]
    plt.figure(figsize=(10, max(5, len(heatmap_data) * 0.28)))
    sns.heatmap(heatmap_data, cmap="viridis_r", cbar_kws={"label": "Vina score (kcal/mol)"})
    plt.title("Conformation-level docking score matrix")
    plt.xlabel("Receptor conformation")
    plt.ylabel("Pocket / compound")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "ensemble_score_matrix_heatmap.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    sns.barplot(data=weights, x="conformation_id", y="weight", hue="pocket")
    plt.xticks(rotation=35, ha="right")
    plt.xlabel("Receptor conformation")
    plt.ylabel("Rank-consistency weight")
    plt.title("EnOpt-style conformation weights")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "enopt_style_conformation_weights.png", dpi=220)
    plt.close()


def write_method_note(matrix: pd.DataFrame, summary: pd.DataFrame, correlations: pd.DataFrame) -> None:
    lines = [
        "# EnOpt-style Ensemble Extension",
        "",
        "This extension adds a small receptor-conformation ensemble to the FtsZ HTVS reconstruction. The current implementation uses normal-mode perturbations around the reconstructed FtsZ receptor and redocks all 90 cleaned BP1/BP2 unique structures across the ensemble.",
        "",
        "The output is described as EnOpt-style because no experimental active/decoy labels are available in the source material. The weighted score is therefore an exploratory conformation-weighted consensus score, not a supervised EnOpt model.",
        "",
        "## Outputs",
        "",
        "- `results/tables/receptor_ensemble_manifest.csv`",
        "- `results/tables/ensemble_vina_scores.csv`",
        "- `results/tables/enopt_style_score_matrix.csv`",
        "- `results/tables/enopt_style_weights.csv`",
        "- `results/tables/enopt_style_top_hits.csv`",
        "- `results/tables/enopt_style_correlations.csv`",
        "- `results/figures/enopt_style_reranking_vs_legacy.png`",
        "- `results/figures/ensemble_score_matrix_heatmap.png`",
        "- `results/figures/enopt_style_conformation_weights.png`",
        "",
        "## Top hits",
        "",
        summary.to_markdown(index=False),
        "",
        "## Correlation with legacy worksheet scores",
        "",
        correlations.to_markdown(index=False) if not correlations.empty else "No correlation summary available.",
        "",
        "## Method boundary",
        "",
        "The current BP1/BP2 candidate set and normal-mode ensemble are suitable for testing the mechanics of ensemble reranking. A strict supervised EnOpt analysis would require experimental active/decoy labels or a benchmark set, plus a receptor ensemble designed for the same binding site.",
    ]
    (ROOT / "docs" / "enopt_extension.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    raw = pd.read_csv(TABLE_DIR / "ensemble_vina_scores.csv")
    raw = raw[raw["returncode"] == 0].copy()
    matrix, weights = build_score_matrix(raw)
    summary = summarize_rankings(matrix)
    correlations = summarize_correlations(matrix)

    matrix.to_csv(TABLE_DIR / "enopt_style_score_matrix.csv", index=False)
    weights.to_csv(TABLE_DIR / "enopt_style_weights.csv", index=False)
    summary.to_csv(TABLE_DIR / "enopt_style_top_hits.csv", index=False)
    correlations.to_csv(TABLE_DIR / "enopt_style_correlations.csv", index=False)
    write_figures(matrix, weights)
    write_method_note(matrix, summary, correlations)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
