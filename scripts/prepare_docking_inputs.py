#!/usr/bin/env python3
"""Prepare FtsZ receptor and ligand PDBQT files for Vina reruns."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ROOT / "data" / "processed"
WORK_DIR = ROOT / "work"
RECEPTOR_PDBQT = WORK_DIR / "receptors" / "FtsZ_receptor_obabel_gasteiger.pdbqt"
OBABEL = ROOT / ".mamba_vina" / "bin" / "obabel"
MEEKO_LIGAND = ROOT / ".venv" / "bin" / "mk_prepare_ligand.py"


def sanitize_message(message: str) -> str:
    return message.replace(str(ROOT) + "/", "").replace(str(ROOT), ".")


def safe_name(value: str, max_length: int = 80) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value[:max_length].strip("_") or "ligand"


def embed_to_sdf(smiles: str, out_path: Path, name: str, properties: dict[str, object]) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260724
    params.useRandomCoords = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=20260724, useRandomCoords=True, maxAttempts=1000)
    if status != 0:
        return False
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    mol.SetProp("_Name", name)
    for key, value in properties.items():
        mol.SetProp(str(key), str(value))
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()
    return True


def embed_to_sdf_loose(smiles: str, out_path: Path, name: str, properties: dict[str, object]) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDG()
    params.randomSeed = 20260724
    params.useRandomCoords = True
    params.ignoreSmoothingFailures = True
    params.enforceChirality = False
    params.maxIterations = 200
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=50)
    except Exception:
        pass
    mol.SetProp("_Name", name)
    for key, value in properties.items():
        mol.SetProp(str(key), str(value))
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()
    return True


def ensure_receptor() -> None:
    if RECEPTOR_PDBQT.exists() and RECEPTOR_PDBQT.stat().st_size > 0:
        return
    RECEPTOR_PDBQT.parent.mkdir(parents=True, exist_ok=True)
    source = ROOT / "data" / "raw" / "drive" / "Receptor_from_fconv.pdb"
    command = [
        str(OBABEL),
        str(source),
        "-xr",
        "-h",
        "--partialcharge",
        "gasteiger",
        "-O",
        str(RECEPTOR_PDBQT),
    ]
    subprocess.run(command, check=True)


def prepare_ligand_pdbqt(sdf_path: Path, pdbqt_path: Path) -> tuple[bool, str]:
    command = [str(MEEKO_LIGAND), "-i", str(sdf_path), "-o", str(pdbqt_path)]
    try:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=45)
        if completed.returncode == 0:
            return True, sanitize_message(completed.stderr.strip() or completed.stdout.strip() or "prepared with Meeko")
        meeko_message = sanitize_message(completed.stderr.strip() or completed.stdout.strip())
    except subprocess.TimeoutExpired:
        meeko_message = "Meeko timed out after 45 seconds"

    fallback_command = [
        str(OBABEL),
        str(sdf_path),
        "--partialcharge",
        "gasteiger",
        "-O",
        str(pdbqt_path),
    ]
    fallback = subprocess.run(fallback_command, text=True, capture_output=True, timeout=45)
    if fallback.returncode == 0 and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
        return True, f"prepared with OpenBabel fallback; Meeko issue: {meeko_message[:180]}"
    return False, f"Meeko issue: {meeko_message[:180]}; OpenBabel issue: {sanitize_message(fallback.stderr or fallback.stdout)[:180]}"


def prepare_large_ligand_with_obabel(smiles: str, sdf_path: Path, pdbqt_path: Path, name: str) -> tuple[bool, str]:
    if embed_to_sdf_loose(smiles, sdf_path, name, {}):
        prepared, message = prepare_ligand_pdbqt(sdf_path, pdbqt_path)
        return prepared, f"loose RDKit 3D fallback; {message}"
    smi_path = sdf_path.with_suffix(".smi")
    smi_path.write_text(f"{smiles}\t{name}\n")
    command = [
        str(OBABEL),
        "-ismi",
        str(smi_path),
        "-osdf",
        "-O",
        str(sdf_path),
        "--gen3d",
        "--partialcharge",
        "gasteiger",
    ]
    try:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=120)
    except subprocess.TimeoutExpired:
        return False, "OpenBabel SDF generation timed out after 120 seconds"
    if completed.returncode != 0 or not sdf_path.exists() or sdf_path.stat().st_size == 0:
        return False, f"OpenBabel SDF issue: {sanitize_message(completed.stderr or completed.stdout)[:180]}"
    prepared, message = prepare_ligand_pdbqt(sdf_path, pdbqt_path)
    return prepared, f"OpenBabel 3D fallback; {message}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=10, help="Top structures per pocket to prepare. Use 0 for all.")
    parser.add_argument("--drug-like-only", action="store_true", help="Filter to rows passing the light Lipinski screen before top-N selection.")
    parser.add_argument("--max-heavy-atoms", type=int, default=45, help="Skip very large structures before 3D embedding.")
    parser.add_argument("--rdkit-heavy-atom-limit", type=int, default=45, help="Use OpenBabel 3D fallback above this heavy-atom count.")
    args = parser.parse_args()

    ensure_receptor()
    boxes = json.loads((PROCESSED_DIR / "pocket_boxes.json").read_text())
    structures = pd.read_csv(PROCESSED_DIR / "bp_unique_structures.csv")
    if args.drug_like_only:
        structures = structures[structures["passes_lipinski_light"] == True].copy()
    structures = structures[structures["heavy_atoms"] <= args.max_heavy_atoms].copy()
    if args.top_n > 0:
        structures = structures.groupby("pocket", group_keys=False).head(args.top_n)

    manifest_rows = []
    for _, row in structures.iterrows():
        pocket = row["pocket"]
        ligand_name = safe_name(f"{pocket}_{row['compound_ids'].split(';')[0]}")
        print(f"Preparing {ligand_name}...", flush=True)
        ligand_dir = WORK_DIR / "ligands" / pocket
        ligand_dir.mkdir(parents=True, exist_ok=True)
        sdf_path = ligand_dir / f"{ligand_name}.sdf"
        pdbqt_path = ligand_dir / f"{ligand_name}.pdbqt"
        properties = {
            "pocket": pocket,
            "compound_ids": row["compound_ids"],
            "source_databases": row["source_databases"],
            "old_best_docking_score_kcal_mol": row["best_docking_score_kcal_mol"],
            "max_tanimoto_to_reference": row["max_tanimoto_to_reference"],
        }
        prepared = False
        prep_message = "embedding failed"
        if int(row["heavy_atoms"]) > args.rdkit_heavy_atom_limit:
            prepared, prep_message = prepare_large_ligand_with_obabel(row["desalted_smiles"], sdf_path, pdbqt_path, ligand_name)
        else:
            embedded = embed_to_sdf(row["desalted_smiles"], sdf_path, ligand_name, properties)
            if not embedded:
                prepared, prep_message = prepare_large_ligand_with_obabel(row["desalted_smiles"], sdf_path, pdbqt_path, ligand_name)
            else:
                prepared, prep_message = prepare_ligand_pdbqt(sdf_path, pdbqt_path)
        manifest_rows.append(
            {
                "pocket": pocket,
                "compound_ids": row["compound_ids"],
                "source_databases": row["source_databases"],
                "desalted_smiles": row["desalted_smiles"],
                "old_best_docking_score_kcal_mol": row["best_docking_score_kcal_mol"],
                "max_tanimoto_to_reference": row["max_tanimoto_to_reference"],
                "ligand_name": ligand_name,
                "sdf_path": str(sdf_path.relative_to(ROOT)),
                "pdbqt_path": str(pdbqt_path.relative_to(ROOT)) if prepared else "",
                "receptor_pdbqt": str(RECEPTOR_PDBQT.relative_to(ROOT)),
                "box_center_x": boxes[pocket]["center"][0],
                "box_center_y": boxes[pocket]["center"][1],
                "box_center_z": boxes[pocket]["center"][2],
                "box_size_x": boxes[pocket]["size"][0],
                "box_size_y": boxes[pocket]["size"][1],
                "box_size_z": boxes[pocket]["size"][2],
                "prepared": prepared,
                "prep_message": prep_message,
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(WORK_DIR / "docking_manifest.csv", index=False)
    print(manifest[["pocket", "compound_ids", "prepared", "prep_message"]].to_string(index=False))


if __name__ == "__main__":
    main()
