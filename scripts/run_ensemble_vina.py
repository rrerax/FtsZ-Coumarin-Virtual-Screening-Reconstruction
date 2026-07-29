#!/usr/bin/env python3
"""Run Vina docking across the receptor ensemble."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK_DIR = ROOT / "work"
TABLE_DIR = ROOT / "results" / "tables"
VINA = ROOT / ".mamba_vina" / "bin" / "vina"
VINA_ROW = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+")


def parse_best_affinity(stdout: str) -> float | None:
    for line in stdout.splitlines():
        match = VINA_ROW.match(line)
        if match:
            return float(match.group(1))
    return None


def run_vina(row: pd.Series, conformation: pd.Series, exhaustiveness: int, cpu: int, seed: int) -> dict[str, object]:
    output_dir = WORK_DIR / "ensemble_docking" / conformation["conformation_id"] / row["pocket"]
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{row['ligand_name']}_vina_out.pdbqt"
    log_path = output_dir / f"{row['ligand_name']}_vina.log"
    command = [
        str(VINA),
        "--receptor",
        str(ROOT / conformation["pdbqt_path"]),
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
    score = parse_best_affinity(completed.stdout)
    return {
        "pocket": row["pocket"],
        "compound_ids": row["compound_ids"],
        "source_databases": row["source_databases"],
        "desalted_smiles": row["desalted_smiles"],
        "ligand_name": row["ligand_name"],
        "conformation_id": conformation["conformation_id"],
        "conformation_mode": conformation["mode"],
        "conformation_direction": conformation["direction"],
        "old_best_docking_score_kcal_mol": row["old_best_docking_score_kcal_mol"],
        "ensemble_vina_score_kcal_mol": score,
        "out_pdbqt": str(out_path.relative_to(ROOT)) if out_path.exists() else "",
        "log_path": str(log_path.relative_to(ROOT)),
        "returncode": completed.returncode,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=str(WORK_DIR / "docking_manifest.csv"))
    parser.add_argument("--ensemble", default=str(TABLE_DIR / "receptor_ensemble_manifest.csv"))
    parser.add_argument("--exhaustiveness", type=int, default=3)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    if not VINA.exists():
        raise FileNotFoundError(f"Vina not found: {VINA}")
    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest["prepared"] == True].copy()
    ensemble = pd.read_csv(args.ensemble)
    ensemble = ensemble[ensemble["prepared"] == True].copy()
    if args.limit > 0:
        manifest = manifest.head(args.limit)

    rows: list[dict[str, object]] = []
    total = len(manifest) * len(ensemble)
    counter = 0
    for _, conformation in ensemble.iterrows():
        for _, row in manifest.iterrows():
            counter += 1
            print(f"[{counter}/{total}] {conformation['conformation_id']} {row['ligand_name']}")
            rows.append(run_vina(row, conformation, args.exhaustiveness, args.cpu, args.seed))

    results = pd.DataFrame(rows)
    out_csv = TABLE_DIR / "ensemble_vina_scores.csv"
    results.to_csv(out_csv, index=False)
    success = int((results["returncode"] == 0).sum()) if not results.empty else 0
    print(f"Wrote {out_csv.relative_to(ROOT)} ({success}/{len(results)} successful jobs)")


if __name__ == "__main__":
    main()
