#!/usr/bin/env python3
"""Draw PyMOL-style binding-mode views for top rerun hits."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"
PYMOL_DIR = ROOT / "results" / "pymol"
RECEPTOR_PDB = ROOT / "data" / "processed" / "FtsZ_receptor_noH.pdb"
FALLBACK_RECEPTOR_PDB = ROOT / "data" / "raw" / "drive" / "Receptor_from_fconv.pdb"


ELEMENT_COLORS = {
    "C": "#3f4652",
    "N": "#2f6de1",
    "O": "#d83f31",
    "S": "#d5a11e",
    "P": "#cc7a00",
    "F": "#27ae60",
    "CL": "#27ae60",
    "BR": "#8b4513",
    "I": "#6b3fa0",
}

COVALENT_RADII = {
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "S": 1.05,
    "P": 1.07,
    "F": 0.57,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
}


@dataclass(frozen=True)
class Atom:
    atom_name: str
    residue_name: str
    chain: str
    residue_number: str
    x: float
    y: float
    z: float
    element: str

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def residue_label(self) -> str:
        return f"{self.residue_name}{self.residue_number}{self.chain}"


def infer_element(line: str, atom_name: str) -> str:
    field_element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if field_element:
        return field_element
    cleaned = "".join(character for character in atom_name.upper() if character.isalpha())
    if cleaned.startswith(("CL", "BR")):
        return cleaned[:2]
    return cleaned[:1] or "C"


def parse_atoms(path: Path, first_model_only: bool = False) -> list[Atom]:
    atoms: list[Atom] = []
    in_first_model = not first_model_only
    for line in path.read_text(errors="ignore").splitlines():
        if first_model_only and line.startswith("MODEL"):
            in_first_model = line.split()[-1] == "1"
            continue
        if first_model_only and line.startswith("ENDMDL") and in_first_model:
            break
        if not in_first_model or not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        element = infer_element(line, atom_name)
        if element == "H":
            continue
        atoms.append(
            Atom(
                atom_name=atom_name,
                residue_name=line[17:20].strip() or "UNL",
                chain=line[21].strip(),
                residue_number=line[22:26].strip() or "1",
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=element,
            )
        )
    return atoms


def distance(first_atom: Atom, second_atom: Atom) -> float:
    return float(np.linalg.norm(first_atom.coord - second_atom.coord))


def likely_bonded(first_atom: Atom, second_atom: Atom) -> bool:
    radius_sum = COVALENT_RADII.get(first_atom.element, 0.76) + COVALENT_RADII.get(second_atom.element, 0.76)
    cutoff = min(2.05, radius_sum + 0.45)
    separation = distance(first_atom, second_atom)
    return 0.35 < separation <= cutoff


def atom_color(atom: Atom) -> str:
    return ELEMENT_COLORS.get(atom.element, "#7f8c8d")


def draw_bonds(axis, atoms: list[Atom], color: str, linewidth: float, alpha: float) -> None:
    for first_index, first_atom in enumerate(atoms):
        for second_atom in atoms[first_index + 1 :]:
            if likely_bonded(first_atom, second_atom):
                xs = [first_atom.x, second_atom.x]
                ys = [first_atom.y, second_atom.y]
                zs = [first_atom.z, second_atom.z]
                axis.plot(xs, ys, zs, color=color, linewidth=linewidth, alpha=alpha)


def residue_centroid(atoms: list[Atom]) -> np.ndarray:
    return np.mean([atom.coord for atom in atoms], axis=0)


def nearest_atom_pair(ligand_atoms: list[Atom], residue_atoms: list[Atom]) -> tuple[Atom, Atom, float]:
    best_ligand = ligand_atoms[0]
    best_receptor = residue_atoms[0]
    best_distance = math.inf
    for ligand_atom in ligand_atoms:
        for receptor_atom in residue_atoms:
            separation = distance(ligand_atom, receptor_atom)
            if separation < best_distance:
                best_ligand = ligand_atom
                best_receptor = receptor_atom
                best_distance = separation
    return best_ligand, best_receptor, best_distance


def plot_binding_mode(axis, pocket: str, hit: pd.Series, receptor_atoms: list[Atom], contact_rows: pd.DataFrame) -> None:
    ligand_atoms = parse_atoms(ROOT / str(hit["out_pdbqt"]), first_model_only=True)
    top_contacts = contact_rows.sort_values("min_distance_a").head(9)
    residue_labels = set(top_contacts["residue"])
    residue_atoms_by_label: dict[str, list[Atom]] = {
        residue_label: [atom for atom in receptor_atoms if atom.residue_label == residue_label]
        for residue_label in residue_labels
    }
    residue_atoms_by_label = {label: atoms for label, atoms in residue_atoms_by_label.items() if atoms}

    ligand_coords = np.array([atom.coord for atom in ligand_atoms])
    ligand_center = ligand_coords.mean(axis=0)

    for residue_label, atoms in residue_atoms_by_label.items():
        residue_coords = np.array([atom.coord for atom in atoms])
        sidechain_color = "#7b8794"
        axis.scatter(
            residue_coords[:, 0],
            residue_coords[:, 1],
            residue_coords[:, 2],
            s=26,
            c=sidechain_color,
            alpha=0.55,
            depthshade=True,
        )
        draw_bonds(axis, atoms, sidechain_color, linewidth=1.0, alpha=0.45)
        centroid = residue_centroid(atoms)
        text_vector = centroid - ligand_center
        vector_length = float(np.linalg.norm(text_vector)) or 1.0
        text_position = centroid + (text_vector / vector_length) * 0.55
        axis.text(
            text_position[0],
            text_position[1],
            text_position[2],
            residue_label,
            fontsize=7,
            color="#2d3748",
            ha="center",
        )

    for row in top_contacts.head(6).itertuples(index=False):
        residue_atoms = residue_atoms_by_label.get(row.residue, [])
        if not residue_atoms:
            continue
        ligand_atom, receptor_atom, separation = nearest_atom_pair(ligand_atoms, residue_atoms)
        contact_color = "#d9480f" if bool(row.has_polar_candidate_within_3_5a) else "#6c757d"
        axis.plot(
            [ligand_atom.x, receptor_atom.x],
            [ligand_atom.y, receptor_atom.y],
            [ligand_atom.z, receptor_atom.z],
            linestyle="--",
            color=contact_color,
            linewidth=1.2,
            alpha=0.75,
        )
        midpoint = (ligand_atom.coord + receptor_atom.coord) / 2
        axis.text(midpoint[0], midpoint[1], midpoint[2], f"{separation:.1f} Å", fontsize=6, color=contact_color)

    for element in sorted({atom.element for atom in ligand_atoms}):
        element_atoms = [atom for atom in ligand_atoms if atom.element == element]
        coords = np.array([atom.coord for atom in element_atoms])
        axis.scatter(
            coords[:, 0],
            coords[:, 1],
            coords[:, 2],
            s=78,
            c=atom_color(element_atoms[0]),
            edgecolor="white",
            linewidth=0.55,
            label=f"Ligand {element}",
            depthshade=True,
        )
    draw_bonds(axis, ligand_atoms, color="#111827", linewidth=2.0, alpha=0.95)

    all_coords = [atom.coord for atom in ligand_atoms]
    for residue_atoms in residue_atoms_by_label.values():
        all_coords.extend(atom.coord for atom in residue_atoms)
    coords = np.array(all_coords)
    center = coords.mean(axis=0)
    span = max(float(np.ptp(coords[:, 0])), float(np.ptp(coords[:, 1])), float(np.ptp(coords[:, 2])), 6.0)
    axis.set_xlim(center[0] - span * 0.48, center[0] + span * 0.48)
    axis.set_ylim(center[1] - span * 0.48, center[1] + span * 0.48)
    axis.set_zlim(center[2] - span * 0.48, center[2] + span * 0.48)
    axis.set_box_aspect((1, 1, 0.82))
    axis.view_init(elev=23, azim=-50 if pocket == "BP1" else -36)
    axis.set_axis_off()
    axis.set_title(
        f"{pocket}: {hit['compound_ids']}\nVina rerun score {hit['new_vina_score_kcal_mol']:.3f} kcal/mol",
        fontsize=12,
        fontweight="bold",
        pad=8,
    )


def write_pymol_script(pocket: str, hit: pd.Series, contact_rows: pd.DataFrame) -> Path:
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)
    ligand_name = str(hit["ligand_name"])
    script_path = PYMOL_DIR / f"{ligand_name}_binding_view.pml"
    residues = "+".join(str(residue).replace(pocket, "") for residue in contact_rows["residue"].head(9))
    residue_selection = " or ".join(
        f"(chain {str(residue)[-1]} and resi {str(residue)[3:-1]})" for residue in contact_rows["residue"].head(9)
    )
    script_path.write_text(
        "\n".join(
            [
                "reinitialize",
                f"load {RECEPTOR_PDB.relative_to(ROOT).as_posix()}, receptor",
                f"load {str(hit['out_pdbqt'])}, ligand",
                "hide everything",
                "show cartoon, receptor",
                "color slate, receptor",
                "show sticks, ligand",
                "color yelloworange, ligand",
                f"select contact_residues, {residue_selection}",
                "show sticks, contact_residues",
                "color marine, contact_residues",
                "show surface, byres contact_residues around 4 of ligand",
                "set transparency, 0.55",
                "distance contacts, ligand, contact_residues, 4.0",
                "set dash_color, orange",
                "set label_size, 18",
                "set ray_opaque_background, off",
                "bg_color white",
                "orient ligand or contact_residues",
                "zoom ligand or contact_residues, 4",
                f"png results/figures/{ligand_name}_pymol_binding_view.png, width=1800, height=1400, dpi=250, ray=1",
                f"# Contact residue labels: {residues}",
            ]
        )
        + "\n"
    )
    return script_path


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    receptor_path = RECEPTOR_PDB if RECEPTOR_PDB.exists() else FALLBACK_RECEPTOR_PDB
    receptor_atoms = parse_atoms(receptor_path)
    rerun = pd.read_csv(TABLE_DIR / "vina_rerun_unique_structures.csv")
    contacts = pd.read_csv(TABLE_DIR / "binding_contacts_top_hits.csv")
    best_hits = rerun.sort_values("new_vina_score_kcal_mol").groupby("pocket", group_keys=False).head(1)

    fig = plt.figure(figsize=(15, 8.5), facecolor="white")
    for subplot_index, hit in enumerate(best_hits.sort_values("pocket").to_dict("records"), start=1):
        pocket = hit["pocket"]
        axis = fig.add_subplot(1, 2, subplot_index, projection="3d")
        pocket_contacts = contacts[contacts["pocket"] == pocket]
        plot_binding_mode(axis, pocket, pd.Series(hit), receptor_atoms, pocket_contacts)
        write_pymol_script(pocket, pd.Series(hit), pocket_contacts)

    fig.suptitle("FtsZ Top Rerun Hits: Binding-Mode Contact Views", fontsize=19, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.035,
        "Ligands are shown as colored sticks/spheres; nearby receptor residues are grey sticks; dashed lines mark closest 4 Å contacts. Computational docking poses only.",
        ha="center",
        fontsize=10,
        color="#6b7280",
    )
    fig.tight_layout(rect=(0.02, 0.06, 0.98, 0.94))
    combined_path = FIGURE_DIR / "binding_mode_top_hits.png"
    fig.savefig(combined_path, dpi=220, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)

    for hit in best_hits.sort_values("pocket").to_dict("records"):
        pocket = hit["pocket"]
        single_fig = plt.figure(figsize=(8, 7), facecolor="white")
        axis = single_fig.add_subplot(1, 1, 1, projection="3d")
        pocket_contacts = contacts[contacts["pocket"] == pocket]
        plot_binding_mode(axis, pocket, pd.Series(hit), receptor_atoms, pocket_contacts)
        single_fig.tight_layout()
        single_path = FIGURE_DIR / f"binding_mode_{pocket}_{hit['ligand_name'].split('_', 1)[1]}.png"
        single_fig.savefig(single_path, dpi=240, bbox_inches="tight", pad_inches=0.12)
        plt.close(single_fig)
        print(single_path.relative_to(ROOT))

    print(combined_path.relative_to(ROOT))
    print(f"PyMOL-ready scripts written to {PYMOL_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
