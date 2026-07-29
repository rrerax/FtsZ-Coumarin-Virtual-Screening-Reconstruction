#!/usr/bin/env python3
"""Create a compact workflow summary figure for the README."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR = ROOT / "results" / "tables"
PROCESSED_DIR = ROOT / "data" / "processed"


def metric_card(ax, title: str, value: str, subtitle: str, color: str) -> None:
    ax.axis("off")
    ax.add_patch(
        plt.Rectangle(
            (0, 0),
            1,
            1,
            transform=ax.transAxes,
            facecolor=color,
            alpha=0.12,
            edgecolor=color,
            linewidth=2,
            clip_on=False,
        )
    )
    ax.text(0.05, 0.72, title, fontsize=12, weight="bold", color="#243447", transform=ax.transAxes)
    ax.text(0.05, 0.38, value, fontsize=26, weight="bold", color=color, transform=ax.transAxes)
    ax.text(0.05, 0.14, subtitle, fontsize=9.5, color="#5f6b7a", transform=ax.transAxes)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(TABLE_DIR / "data_quality_summary.csv")
    unique_structures = pd.read_csv(PROCESSED_DIR / "bp_unique_structures.csv")
    rerun = pd.read_csv(TABLE_DIR / "vina_rerun_unique_structures.csv")
    contacts = pd.read_csv(TABLE_DIR / "binding_contacts_top_hits.csv")

    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    fig = plt.figure(figsize=(19, 11.5), facecolor="white")
    grid = fig.add_gridspec(3, 4, height_ratios=[0.82, 1.55, 1.55], hspace=0.58, wspace=0.78)

    fig.suptitle("FtsZ Coumarin HTVS Reconstruction", fontsize=24, weight="bold", x=0.03, ha="left", y=0.98)
    fig.text(
        0.03,
        0.935,
        "Legacy HTVS records curated into a reproducible two-pocket rerun and EnOpt-style conformation ensemble",
        fontsize=12.5,
        color="#5f6b7a",
    )

    method_records = int(summary["method_records"].sum())
    unique_count = int(len(unique_structures))
    rerun_success = int((rerun["returncode"] == 0).sum())
    ensemble = pd.read_csv(TABLE_DIR / "ensemble_vina_scores.csv")
    ensemble_success = int((ensemble["returncode"] == 0).sum())
    ensemble_total = int(len(ensemble))

    metric_card(fig.add_subplot(grid[0, 0]), "Cleaned records", f"{method_records}", "BP1 + BP2 method-level rows", "#2563eb")
    metric_card(fig.add_subplot(grid[0, 1]), "Unique structures", f"{unique_count}", "After salt stripping and deduplication", "#16a34a")
    metric_card(fig.add_subplot(grid[0, 2]), "Vina rerun", f"{rerun_success}/{len(rerun)}", "Cleaned unique structures completed", "#ea580c")
    metric_card(fig.add_subplot(grid[0, 3]), "EnOpt-style jobs", f"{ensemble_success}/{ensemble_total}", "Five-conformation rerun completed", "#7c3aed")

    ax1 = fig.add_subplot(grid[1, :2])
    sns.histplot(
        data=unique_structures,
        x="best_docking_score_kcal_mol",
        hue="pocket",
        bins=14,
        kde=True,
        element="step",
        ax=ax1,
    )
    ax1.set_title("Original candidate score distribution", fontsize=14, weight="bold")
    ax1.set_xlabel("Best worksheet docking score (kcal/mol; lower is better)")
    ax1.set_ylabel("Unique structures")

    ax2 = fig.add_subplot(grid[1, 2:])
    sns.scatterplot(
        data=rerun,
        x="old_best_docking_score_kcal_mol",
        y="new_vina_score_kcal_mol",
        hue="pocket",
        s=90,
        ax=ax2,
    )
    low = min(rerun["old_best_docking_score_kcal_mol"].min(), rerun["new_vina_score_kcal_mol"].min()) - 0.3
    high = max(rerun["old_best_docking_score_kcal_mol"].max(), rerun["new_vina_score_kcal_mol"].max()) + 0.3
    ax2.plot([low, high], [low, high], linestyle="--", color="#8a8a8a", linewidth=1.2)
    ax2.set_title("Docking rerun vs legacy worksheet scores", fontsize=14, weight="bold")
    ax2.set_xlabel("Original score (kcal/mol)")
    ax2.set_ylabel("Rerun Vina score (kcal/mol)")

    ax3 = fig.add_subplot(grid[2, :2])
    top_plot = unique_structures.groupby("pocket", group_keys=False).head(6).copy()
    top_plot["label"] = top_plot["compound_ids"].str.split(";").str[0] + " (" + top_plot["pocket"] + ")"
    sns.barplot(
        data=top_plot.sort_values("best_docking_score_kcal_mol"),
        x="best_docking_score_kcal_mol",
        y="label",
        hue="pocket",
        dodge=False,
        ax=ax3,
    )
    ax3.set_title("Top structures by original score", fontsize=14, weight="bold")
    ax3.set_xlabel("Best docking score (kcal/mol)")
    ax3.set_ylabel("")
    ax3.tick_params(axis="y", labelsize=9)

    ax4 = fig.add_subplot(grid[2, 2:])
    contact_plot = contacts.sort_values(["pocket", "min_distance_a"]).groupby("pocket", group_keys=False).head(7).copy()
    contact_plot["label"] = contact_plot["residue"] + " (" + contact_plot["pocket"] + ")"
    sns.barplot(
        data=contact_plot.sort_values("min_distance_a"),
        x="min_distance_a",
        y="label",
        hue="pocket",
        dodge=False,
        ax=ax4,
    )
    ax4.axvline(3.5, color="#8a8a8a", linestyle="--", linewidth=1)
    ax4.set_title("Nearest residues for best rerun hits", fontsize=14, weight="bold")
    ax4.set_xlabel("Minimum ligand-residue distance (Å)")
    ax4.set_ylabel("")
    ax4.tick_params(axis="y", labelsize=9)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.spines[["top", "right"]].set_visible(False)

    fig.text(
        0.03,
        0.012,
        "Note: this reconstruction reports computational screening consistency only; docking scores are not experimental activity evidence.",
        fontsize=10.5,
        color="#6b7280",
    )
    out = FIGURE_DIR / "workflow_summary.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
