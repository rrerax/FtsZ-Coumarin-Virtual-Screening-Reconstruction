#!/usr/bin/env python3
"""Run AutoDock Vina for prepared FtsZ ligand PDBQT files."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
WORK_DIR = ROOT / "work"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
VINA = ROOT / ".mamba_vina" / "bin" / "vina"


VINA_ROW = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")


def parse_best_affinity(stdout: str) -> float | None:
    for line in stdout.splitlines():
        match = VINA_ROW.match(line)
        if match:
            return float(match.group(1))
    return None


def run_vina(row: pd.Series, exhaustiveness: int, cpu: int, seed: int) -> dict[str, object]:
    output_dir = WORK_DIR / "docking" / row["pocket"]
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{row['ligand_name']}_vina_out.pdbqt"
    log_path = output_dir / f"{row['ligand_name']}_vina.log"
    command = [
        str(VINA),
        "--receptor",
        str(ROOT / row["receptor_pdbqt"]),
        "--ligand",
        str(ROOT / row["pdbqt_path"]),
        "--center_x",
        str(row["box_center_x"]),
        "--center_y",
        str(row["box_center_y"]),
        "--center_z",
        str(row["box_center_z"]),
        "--size_x",
        str(row["box_size_x"]),
        "--size_y",
        str(row["box_size_y"]),
        "--size_z",
        str(row["box_size_z"]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        "5",
        "--cpu",
        str(cpu),
        "--seed",
        str(seed),
        "--out",
        str(out_path),
    ]
    completed = subprocess.run(command, text=True, capture_output=True)
    log_path.write_text(completed.stdout + "\n" + completed.stderr)
    return {
        "pocket": row["pocket"],
        "compound_ids": row["compound_ids"],
        "source_databases": row["source_databases"],
        "desalted_smiles": row["desalted_smiles"],
        "old_best_docking_score_kcal_mol": row["old_best_docking_score_kcal_mol"],
        "new_vina_score_kcal_mol": parse_best_affinity(completed.stdout),
        "score_delta_new_minus_old": None,
        "ligand_name": row["ligand_name"],
        "out_pdbqt": str(out_path.relative_to(ROOT)) if out_path.exists() else "",
        "log_path": str(log_path.relative_to(ROOT)),
        "returncode": completed.returncode,
    }


def write_comparison_figure(results: pd.DataFrame) -> None:
    valid = results.dropna(subset=["new_vina_score_kcal_mol"]).copy()
    if valid.empty:
        return
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=valid,
        x="old_best_docking_score_kcal_mol",
        y="new_vina_score_kcal_mol",
        hue="pocket",
        s=80,
    )
    low = min(valid["old_best_docking_score_kcal_mol"].min(), valid["new_vina_score_kcal_mol"].min()) - 0.3
    high = max(valid["old_best_docking_score_kcal_mol"].max(), valid["new_vina_score_kcal_mol"].max()) + 0.3
    plt.plot([low, high], [low, high], linestyle="--", color="gray", linewidth=1)
    plt.xlabel("Original thesis/worksheet score (kcal/mol)")
    plt.ylabel("Rerun Vina score (kcal/mol)")
    plt.title("FtsZ top-candidate docking rerun vs legacy worksheet scores")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "vina_rerun_vs_original_scores.png", dpi=220)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=str(WORK_DIR / "docking_manifest.csv"))
    parser.add_argument("--exhaustiveness", type=int, default=4)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument("--limit", type=int, default=0, help="Limit number of prepared rows for smoke testing. 0 means all.")
    args = parser.parse_args()

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest["prepared"] == True].copy()
    if args.limit > 0:
        manifest = manifest.head(args.limit)

    rows = [run_vina(row, args.exhaustiveness, args.cpu, args.seed) for _, row in manifest.iterrows()]
    results = pd.DataFrame(rows)
    if not results.empty:
        results["score_delta_new_minus_old"] = results["new_vina_score_kcal_mol"] - results["old_best_docking_score_kcal_mol"]
    out_csv = TABLE_DIR / "vina_rerun_top_structures.csv"
    results.to_csv(out_csv, index=False)
    write_comparison_figure(results)
    print(results[["pocket", "compound_ids", "old_best_docking_score_kcal_mol", "new_vina_score_kcal_mol", "returncode"]].to_string(index=False))


if __name__ == "__main__":
    main()
