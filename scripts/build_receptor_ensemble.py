#!/usr/bin/env python3
"""Build a small normal-mode FtsZ receptor ensemble for docking reruns."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import subprocess
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="pkg_resources is deprecated.*")
from prody import ANM, parsePDB


ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ROOT / "data" / "processed"
WORK_DIR = ROOT / "work"
ENSEMBLE_DIR = WORK_DIR / "receptors" / "ensemble"
TABLE_DIR = ROOT / "results" / "tables"
OBABEL = ROOT / ".mamba_vina" / "bin" / "obabel"
DEFAULT_RECEPTOR = PROCESSED_DIR / "FtsZ_receptor_noH.pdb"
FALLBACK_RECEPTOR = ROOT / "data" / "raw" / "drive" / "Receptor_from_fconv.pdb"


@dataclass(frozen=True)
class AtomLine:
    line: str
    key: tuple[str, str, str]
    coord: np.ndarray


def parse_pdb_atom_lines(path: Path) -> tuple[list[str], list[AtomLine]]:
    lines = path.read_text(errors="ignore").splitlines()
    atom_lines: list[AtomLine] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21].strip(), line[22:26].strip(), line[27].strip())
        coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
        atom_lines.append(AtomLine(line=line, key=key, coord=coord))
    return lines, atom_lines


def write_perturbed_pdb(path: Path, lines: list[str], atom_vectors: dict[tuple[str, str, str], np.ndarray], scale: float) -> None:
    output: list[str] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            output.append(line)
            continue
        key = (line[21].strip(), line[22:26].strip(), line[27].strip())
        shift = atom_vectors.get(key)
        if shift is None:
            output.append(line)
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float) + shift * scale
        output.append(f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}")
    path.write_text("\n".join(output) + "\n")


def residue_key_from_atom(atom) -> tuple[str, str, str]:
    return (atom.getChid().strip(), str(atom.getResnum()).strip(), atom.getIcode().strip())


def build_mode_vectors(receptor_path: Path, mode_count: int) -> list[tuple[str, dict[tuple[str, str, str], np.ndarray], float]]:
    structure = parsePDB(str(receptor_path))
    ca_atoms = structure.select("protein and name CA")
    if ca_atoms is None or len(ca_atoms) < 20:
        raise ValueError("Not enough CA atoms for normal-mode ensemble generation")
    anm = ANM("FtsZ_ANM")
    anm.buildHessian(ca_atoms)
    anm.calcModes(max(mode_count + 2, 6))
    rows = []
    for mode_index in range(mode_count):
        mode = anm[mode_index]
        array = mode.getArrayNx3()
        vectors: dict[tuple[str, str, str], np.ndarray] = {}
        for atom, vector in zip(ca_atoms, array):
            vectors[residue_key_from_atom(atom)] = np.array(vector, dtype=float)
        rms = float(np.sqrt(np.mean(np.sum(array**2, axis=1))))
        rows.append((f"mode{mode_index + 1}", vectors, rms))
    return rows


def prepare_receptor_pdbqt(pdb_path: Path, pdbqt_path: Path) -> tuple[bool, str]:
    command = [str(OBABEL), "-ipdb", str(pdb_path), "-opdbqt", "-O", str(pdbqt_path), "-xr", "--partialcharge", "gasteiger"]
    completed = subprocess.run(command, text=True, capture_output=True)
    message = (completed.stdout + "\n" + completed.stderr).strip().replace("\n", " | ")
    message = message.replace(str(ROOT) + "/", "").replace(str(ROOT), ".")
    return completed.returncode == 0 and pdbqt_path.exists(), message


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", default="", help="Optional receptor PDB path. Defaults to processed FtsZ receptor.")
    parser.add_argument("--target-ca-rmsd", type=float, default=0.65)
    parser.add_argument("--modes", type=int, default=2)
    args = parser.parse_args()

    receptor_path = Path(args.receptor) if args.receptor else (DEFAULT_RECEPTOR if DEFAULT_RECEPTOR.exists() else FALLBACK_RECEPTOR)
    if not receptor_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_path}")
    if not OBABEL.exists():
        raise FileNotFoundError(f"OpenBabel not found: {OBABEL}")

    ENSEMBLE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    lines, atom_lines = parse_pdb_atom_lines(receptor_path)
    base_pdb = ENSEMBLE_DIR / "FtsZ_conf00_base.pdb"
    base_pdb.write_text("\n".join(lines) + "\n")

    rows = []
    base_pdbqt = ENSEMBLE_DIR / "FtsZ_conf00_base.pdbqt"
    prepared, message = prepare_receptor_pdbqt(base_pdb, base_pdbqt)
    rows.append(
        {
            "conformation_id": "conf00_base",
            "mode": "base",
            "direction": 0,
            "target_ca_rmsd_a": 0.0,
            "pdb_path": str(base_pdb.relative_to(ROOT)),
            "pdbqt_path": str(base_pdbqt.relative_to(ROOT)),
            "prepared": prepared,
            "prep_message": message,
        }
    )

    mode_vectors = build_mode_vectors(receptor_path, args.modes)
    conformation_number = 1
    for mode_name, vectors, vector_rms in mode_vectors:
        if vector_rms <= 0:
            continue
        scale = args.target_ca_rmsd / vector_rms
        for direction in (-1, 1):
            conformation_id = f"conf{conformation_number:02d}_{mode_name}_{'minus' if direction < 0 else 'plus'}"
            pdb_path = ENSEMBLE_DIR / f"FtsZ_{conformation_id}.pdb"
            pdbqt_path = ENSEMBLE_DIR / f"FtsZ_{conformation_id}.pdbqt"
            write_perturbed_pdb(pdb_path, lines, vectors, scale * direction)
            prepared, message = prepare_receptor_pdbqt(pdb_path, pdbqt_path)
            rows.append(
                {
                    "conformation_id": conformation_id,
                    "mode": mode_name,
                    "direction": direction,
                    "target_ca_rmsd_a": args.target_ca_rmsd,
                    "pdb_path": str(pdb_path.relative_to(ROOT)),
                    "pdbqt_path": str(pdbqt_path.relative_to(ROOT)),
                    "prepared": prepared,
                    "prep_message": message,
                }
            )
            conformation_number += 1

    ensemble = pd.DataFrame(rows)
    ensemble.to_csv(TABLE_DIR / "receptor_ensemble_manifest.csv", index=False)
    print(ensemble[["conformation_id", "mode", "direction", "target_ca_rmsd_a", "prepared"]].to_string(index=False))


if __name__ == "__main__":
    main()
