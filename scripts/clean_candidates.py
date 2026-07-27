#!/usr/bin/env python3
"""Clean thesis BP1/BP2 candidate spreadsheets into reproducible tables."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw" / "drive"
PROCESSED_DIR = ROOT / "data" / "processed"
TABLE_DIR = ROOT / "results" / "tables"
FIGURE_DIR = ROOT / "results" / "figures"


@dataclass(frozen=True)
class PocketInput:
    pocket: str
    file_name: str


INPUTS = [
    PocketInput("BP1", "BP1_candidates.xlsx"),
    PocketInput("BP2", "BP2_candidates.xlsx"),
]


def infer_source_database(compound_id: str) -> str:
    identifier = str(compound_id).strip()
    upper = identifier.upper()
    if upper.startswith("CHEMBL"):
        return "ChEMBL23"
    if upper.startswith("SCHEMBL"):
        return "SureChEMBL"
    if upper.startswith("ZINC"):
        return "ZINC15"
    if upper.startswith("DB"):
        return "DrugBank"
    if upper.startswith("HMDB"):
        return "HMDB"
    if upper.startswith("LSM"):
        return "LINCS"
    if upper.startswith("MOLPORT"):
        return "MolPort"
    if re.fullmatch(r"[A-Z0-9]{3}", upper):
        return "PDB ligand"
    if re.fullmatch(r"\d+", upper):
        return "numeric library ID"
    return "unknown"


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return mol
    return max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())


def canonicalize_smiles(smiles: object) -> dict[str, object]:
    if pd.isna(smiles):
        return {"valid_smiles": False, "canonical_smiles": None, "desalted_smiles": None}

    raw_smiles = str(smiles).strip()
    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
        return {"valid_smiles": False, "canonical_smiles": None, "desalted_smiles": None}

    desalted = largest_fragment(mol)
    return {
        "valid_smiles": True,
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "desalted_smiles": Chem.MolToSmiles(desalted, canonical=True),
        "mol_weight": round(Descriptors.MolWt(desalted), 3),
        "logp": round(Crippen.MolLogP(desalted), 3),
        "hbd": Lipinski.NumHDonors(desalted),
        "hba": Lipinski.NumHAcceptors(desalted),
        "tpsa": round(rdMolDescriptors.CalcTPSA(desalted), 3),
        "rotatable_bonds": Lipinski.NumRotatableBonds(desalted),
        "heavy_atoms": Descriptors.HeavyAtomCount(desalted),
    }


def read_main_candidate_block(input_spec: PocketInput) -> pd.DataFrame:
    path = RAW_DIR / input_spec.file_name
    df = pd.read_excel(path)

    main_columns = list(df.columns[:5])
    records = df[main_columns].copy()
    records.columns = [
        "compound_id",
        "docking_score_kcal_mol",
        "tanimoto_to_reference",
        "smiles",
        "fingerprint_method",
    ]
    records["pocket"] = input_spec.pocket
    records["source_file"] = input_spec.file_name

    records = records.dropna(subset=["compound_id", "smiles"]).copy()
    records["compound_id"] = records["compound_id"].astype(str).str.strip()
    records["smiles"] = records["smiles"].astype(str).str.strip()
    records["fingerprint_method"] = records["fingerprint_method"].astype(str).str.strip()
    records["docking_score_kcal_mol"] = pd.to_numeric(records["docking_score_kcal_mol"], errors="coerce")
    records["tanimoto_to_reference"] = pd.to_numeric(records["tanimoto_to_reference"], errors="coerce")
    records = records.dropna(subset=["docking_score_kcal_mol", "tanimoto_to_reference"])

    descriptor_rows = records["smiles"].apply(canonicalize_smiles).apply(pd.Series)
    records = pd.concat([records.reset_index(drop=True), descriptor_rows.reset_index(drop=True)], axis=1)
    records["source_database"] = records["compound_id"].apply(infer_source_database)
    records["passes_lipinski_light"] = (
        (records["mol_weight"] <= 500)
        & (records["logp"] <= 5)
        & (records["hbd"] <= 5)
        & (records["hba"] <= 10)
    )
    return records


def build_unique_table(method_records: pd.DataFrame) -> pd.DataFrame:
    valid = method_records[method_records["valid_smiles"]].copy()
    grouped = (
        valid.groupby(["pocket", "compound_id", "desalted_smiles"], dropna=False)
        .agg(
            source_database=("source_database", "first"),
            best_docking_score_kcal_mol=("docking_score_kcal_mol", "min"),
            median_docking_score_kcal_mol=("docking_score_kcal_mol", "median"),
            max_tanimoto_to_reference=("tanimoto_to_reference", "max"),
            fingerprint_methods=("fingerprint_method", lambda values: ";".join(sorted(set(map(str, values))))),
            observation_count=("fingerprint_method", "size"),
            original_smiles=("smiles", "first"),
            canonical_smiles=("canonical_smiles", "first"),
            mol_weight=("mol_weight", "first"),
            logp=("logp", "first"),
            hbd=("hbd", "first"),
            hba=("hba", "first"),
            tpsa=("tpsa", "first"),
            rotatable_bonds=("rotatable_bonds", "first"),
            heavy_atoms=("heavy_atoms", "first"),
            passes_lipinski_light=("passes_lipinski_light", "first"),
        )
        .reset_index()
    )
    return grouped.sort_values(["pocket", "best_docking_score_kcal_mol", "compound_id"]).reset_index(drop=True)


def build_score_matrix(unique_records: pd.DataFrame) -> pd.DataFrame:
    matrix = unique_records.pivot_table(
        index="desalted_smiles",
        columns="pocket",
        values="best_docking_score_kcal_mol",
        aggfunc="min",
    ).reset_index()

    metadata = (
        unique_records.sort_values("best_docking_score_kcal_mol")
        .groupby("desalted_smiles", as_index=False)
        .agg(
            compound_ids=("compound_id", lambda values: ";".join(sorted(set(map(str, values))))),
            source_databases=("source_database", lambda values: ";".join(sorted(set(map(str, values))))),
            max_tanimoto_to_reference=("max_tanimoto_to_reference", "max"),
            mol_weight=("mol_weight", "first"),
            logp=("logp", "first"),
            hbd=("hbd", "first"),
            hba=("hba", "first"),
            tpsa=("tpsa", "first"),
            rotatable_bonds=("rotatable_bonds", "first"),
            heavy_atoms=("heavy_atoms", "first"),
            passes_lipinski_light=("passes_lipinski_light", "first"),
        )
    )
    matrix = matrix.merge(metadata, on="desalted_smiles", how="left")
    pocket_columns = [column for column in ["BP1", "BP2"] if column in matrix.columns]
    matrix["ensemble_best_score_kcal_mol"] = matrix[pocket_columns].min(axis=1, skipna=True)
    matrix["ensemble_mean_score_kcal_mol"] = matrix[pocket_columns].mean(axis=1, skipna=True)
    matrix["observed_pocket_count"] = matrix[pocket_columns].notna().sum(axis=1)
    matrix["rank_ensemble_best"] = matrix["ensemble_best_score_kcal_mol"].rank(method="min", ascending=True).astype(int)
    matrix["rank_ensemble_mean"] = matrix["ensemble_mean_score_kcal_mol"].rank(method="min", ascending=True).astype(int)
    ordered = [
        "desalted_smiles",
        "compound_ids",
        "source_databases",
        "BP1",
        "BP2",
        "ensemble_best_score_kcal_mol",
        "ensemble_mean_score_kcal_mol",
        "observed_pocket_count",
        "rank_ensemble_best",
        "rank_ensemble_mean",
        "max_tanimoto_to_reference",
        "mol_weight",
        "logp",
        "hbd",
        "hba",
        "tpsa",
        "rotatable_bonds",
        "heavy_atoms",
        "passes_lipinski_light",
    ]
    return matrix[ordered].sort_values("rank_ensemble_best").reset_index(drop=True)


def build_unique_structure_table(unique_records: pd.DataFrame) -> pd.DataFrame:
    structure_records = (
        unique_records.groupby(["pocket", "desalted_smiles"], dropna=False)
        .agg(
            compound_ids=("compound_id", lambda values: ";".join(sorted(set(map(str, values))))),
            source_databases=("source_database", lambda values: ";".join(sorted(set(map(str, values))))),
            best_docking_score_kcal_mol=("best_docking_score_kcal_mol", "min"),
            median_docking_score_kcal_mol=("median_docking_score_kcal_mol", "median"),
            max_tanimoto_to_reference=("max_tanimoto_to_reference", "max"),
            fingerprint_methods=("fingerprint_methods", lambda values: ";".join(sorted(set(";".join(map(str, values)).split(";"))))),
            source_id_count=("compound_id", "nunique"),
            mol_weight=("mol_weight", "first"),
            logp=("logp", "first"),
            hbd=("hbd", "first"),
            hba=("hba", "first"),
            tpsa=("tpsa", "first"),
            rotatable_bonds=("rotatable_bonds", "first"),
            heavy_atoms=("heavy_atoms", "first"),
            passes_lipinski_light=("passes_lipinski_light", "first"),
        )
        .reset_index()
    )
    return structure_records.sort_values(["pocket", "best_docking_score_kcal_mol", "compound_ids"]).reset_index(drop=True)


def write_figures(unique_records: pd.DataFrame, score_matrix: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(8, 5))
    sns.histplot(
        data=unique_records,
        x="best_docking_score_kcal_mol",
        hue="pocket",
        bins=14,
        kde=True,
        element="step",
    )
    plt.xlabel("Best docking score (kcal/mol; lower is better)")
    plt.ylabel("Unique compounds")
    plt.title("FtsZ candidate score distribution by pocket")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "score_distribution_by_pocket.png", dpi=220)
    plt.close()

    top_records = unique_records.groupby("pocket", group_keys=False).head(10).copy()
    top_records["compound_label"] = top_records["compound_id"] + " (" + top_records["pocket"] + ")"
    plt.figure(figsize=(10, 7))
    sns.barplot(
        data=top_records.sort_values("best_docking_score_kcal_mol"),
        y="compound_label",
        x="best_docking_score_kcal_mol",
        hue="pocket",
        dodge=False,
    )
    plt.xlabel("Best docking score (kcal/mol; lower is better)")
    plt.ylabel("Compound")
    plt.title("Top 10 candidate compounds per FtsZ pocket")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "top10_scores_by_pocket.png", dpi=220)
    plt.close()

    overlap = score_matrix[score_matrix["BP1"].notna() & score_matrix["BP2"].notna()].copy()
    if not overlap.empty:
        plt.figure(figsize=(6, 6))
        sns.scatterplot(
            data=overlap,
            x="BP1",
            y="BP2",
            size="max_tanimoto_to_reference",
            sizes=(40, 160),
            legend=False,
        )
        plt.xlabel("BP1 score (kcal/mol)")
        plt.ylabel("BP2 score (kcal/mol)")
        plt.title("Compounds observed in both pocket candidate sets")
        plt.tight_layout()
        plt.savefig(FIGURE_DIR / "bp_score_scatter_overlap.png", dpi=220)
        plt.close()


def write_markdown_summary(
    method_records: pd.DataFrame,
    unique_records: pd.DataFrame,
    score_matrix: pd.DataFrame,
) -> None:
    lines = [
        "# Data Audit — FtsZ Coumarin HTVS Reconstruction",
        "",
        "## Inputs",
        "- `BP1_candidates.xlsx`: thesis spreadsheet labeled as pocket 1 / BP1.",
        "- `BP2_candidates.xlsx`: thesis spreadsheet labeled as pocket 2 / BP2.",
        "- `Receptor_from_fconv.pdb`: FtsZ receptor copied from the ezPocket result folder.",
        "- `ezPocket_fconv.tsv`: two pocket centers and volumes from the source project folder.",
        "",
        "## Cleaning rule",
        "The first spreadsheet block is interpreted as: compound ID, docking score, Tanimoto similarity, SMILES, and fingerprint method. The spreadsheet column labels are inherited from the old workbook and are therefore not used literally.",
        "",
        "## Counts",
    ]
    raw_summary = method_records.groupby("pocket").agg(
        method_records=("compound_id", "size"),
        valid_smiles=("valid_smiles", "sum"),
    )
    unique_summary = unique_records.groupby("pocket").agg(unique_compounds=("compound_id", "size"))
    summary = raw_summary.join(unique_summary)
    lines.append(summary.to_markdown())
    lines += [
        "",
        "## EnOpt readiness note",
        "The cleaned score matrix is EnOpt-compatible at the data-shape level, but the two score columns currently represent two binding pockets, not two conformations of the same pocket. For a methodologically strict EnOpt run, add multiple receptor conformations for the same binding site and experimental active/inactive labels or a documented surrogate label.",
        "",
        "## Key outputs",
        "- `data/processed/bp_method_records.csv`",
        "- `data/processed/bp_unique_compounds.csv`",
        "- `data/processed/bp_unique_structures.csv`",
        "- `data/processed/enopt_score_matrix.csv`",
        "- `results/tables/top10_by_pocket.csv`",
        "- `results/tables/top10_unique_structures_by_pocket.csv`",
        "- `results/figures/score_distribution_by_pocket.png`",
        "- `results/figures/top10_scores_by_pocket.png`",
    ]
    if (FIGURE_DIR / "bp_score_scatter_overlap.png").exists():
        lines.append("- `results/figures/bp_score_scatter_overlap.png`")
    if (FIGURE_DIR / "top_ligand_structures.png").exists():
        lines.append("- `results/figures/top_ligand_structures.png`")
    if (FIGURE_DIR / "binding_mode_top_hits.png").exists():
        lines.append("- `results/figures/binding_mode_top_hits.png`")
    if (ROOT / "results" / "pymol").exists():
        lines.append("- `results/pymol/`: PyMOL-ready binding-view scripts")
    (ROOT / "docs" / "data_audit.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    method_records = pd.concat([read_main_candidate_block(item) for item in INPUTS], ignore_index=True)
    unique_records = build_unique_table(method_records)
    unique_structures = build_unique_structure_table(unique_records)
    score_matrix = build_score_matrix(unique_records)

    method_records.to_csv(PROCESSED_DIR / "bp_method_records.csv", index=False)
    unique_records.to_csv(PROCESSED_DIR / "bp_unique_compounds.csv", index=False)
    unique_structures.to_csv(PROCESSED_DIR / "bp_unique_structures.csv", index=False)
    score_matrix.to_csv(PROCESSED_DIR / "enopt_score_matrix.csv", index=False)

    top10 = unique_records.groupby("pocket", group_keys=False).head(10)
    top10.to_csv(TABLE_DIR / "top10_by_pocket.csv", index=False)
    top10_structures = unique_structures.groupby("pocket", group_keys=False).head(10)
    top10_structures.to_csv(TABLE_DIR / "top10_unique_structures_by_pocket.csv", index=False)

    quality_summary = (
        method_records.groupby("pocket")
        .agg(
            method_records=("compound_id", "size"),
            valid_method_records=("valid_smiles", "sum"),
            unique_ids=("compound_id", "nunique"),
            min_score=("docking_score_kcal_mol", "min"),
            median_score=("docking_score_kcal_mol", "median"),
            max_score=("docking_score_kcal_mol", "max"),
        )
        .reset_index()
    )
    quality_summary.to_csv(TABLE_DIR / "data_quality_summary.csv", index=False)

    write_figures(unique_records, score_matrix)
    write_markdown_summary(method_records, unique_records, score_matrix)

    print("Method records:", len(method_records))
    print("Unique pocket-compound records:", len(unique_records))
    print("Unique pocket-structure records:", len(unique_structures))
    print("EnOpt score matrix rows:", len(score_matrix))
    print(quality_summary.to_string(index=False))


if __name__ == "__main__":
    main()
