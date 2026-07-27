#!/usr/bin/env python3
"""Summarize receptor residues near rerun Vina poses."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RECEPTOR = ROOT / "work" / "receptors" / "FtsZ_receptor_obabel_gasteiger.pdbqt"
TABLE_DIR = ROOT / "results" / "tables"


def parse_atoms(path: Path, first_model_only: bool = False) -> list[dict[str, object]]:
    atoms = []
    in_first_model = not first_model_only
    for line in path.read_text(errors="ignore").splitlines():
        if first_model_only and line.startswith("MODEL"):
            in_first_model = line.split()[-1] == "1"
            continue
        if first_model_only and line.startswith("ENDMDL") and in_first_model:
            break
        if not in_first_model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        element = atom_name[0].upper()
        if element == "H":
            continue
        atoms.append(
            {
                "atom_name": atom_name,
                "residue_name": line[17:20].strip(),
                "chain": line[21].strip(),
                "residue_number": line[22:26].strip(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": element,
            }
        )
    return atoms


def distance(a: dict[str, object], b: dict[str, object]) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def main() -> None:
    receptor_atoms = parse_atoms(RECEPTOR)
    rerun = pd.read_csv(TABLE_DIR / "vina_rerun_top_structures.csv")
    best_by_pocket = rerun.sort_values("new_vina_score_kcal_mol").groupby("pocket", group_keys=False).head(1)
    rows = []
    for hit in best_by_pocket.itertuples(index=False):
        ligand_atoms = parse_atoms(ROOT / hit.out_pdbqt, first_model_only=True)
        residue_contacts: dict[str, list[tuple[float, str, str, str]]] = defaultdict(list)
        for ligand_atom in ligand_atoms:
            for rec_atom in receptor_atoms:
                dist = distance(ligand_atom, rec_atom)
                if dist <= 4.0:
                    residue = f"{rec_atom['residue_name']}{rec_atom['residue_number']}{rec_atom['chain']}"
                    contact_type = "polar_candidate" if dist <= 3.5 and ligand_atom["element"] in {"N", "O", "S"} and rec_atom["element"] in {"N", "O", "S"} else "close_contact"
                    residue_contacts[residue].append((dist, ligand_atom["atom_name"], rec_atom["atom_name"], contact_type))
        for residue, contacts in residue_contacts.items():
            contacts = sorted(contacts, key=lambda item: item[0])
            rows.append(
                {
                    "pocket": hit.pocket,
                    "compound_ids": hit.compound_ids,
                    "new_vina_score_kcal_mol": hit.new_vina_score_kcal_mol,
                    "residue": residue,
                    "min_distance_a": round(contacts[0][0], 3),
                    "contact_count_within_4a": len(contacts),
                    "has_polar_candidate_within_3_5a": any(item[3] == "polar_candidate" for item in contacts),
                    "closest_ligand_atom": contacts[0][1],
                    "closest_receptor_atom": contacts[0][2],
                }
            )
    out = pd.DataFrame(rows).sort_values(["pocket", "min_distance_a"])
    out.to_csv(TABLE_DIR / "binding_contacts_top_hits.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
