#!/usr/bin/env python3
"""Draw 2D structures for top FtsZ candidate ligands."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "results" / "figures"


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    top = pd.read_csv(ROOT / "results" / "tables" / "top10_unique_structures_by_pocket.csv")
    selected = top.groupby("pocket", group_keys=False).head(6).copy()
    mols = [Chem.MolFromSmiles(smiles) for smiles in selected["desalted_smiles"]]
    legends = [
        f"{row.pocket}: {row.compound_ids}\nscore={row.best_docking_score_kcal_mol:.3f}"
        for row in selected.itertuples(index=False)
    ]
    image = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        subImgSize=(340, 240),
        legends=legends,
        useSVG=False,
    )
    image.save(FIGURE_DIR / "top_ligand_structures.png")
    print(FIGURE_DIR / "top_ligand_structures.png")


if __name__ == "__main__":
    main()
