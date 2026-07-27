#!/usr/bin/env python3
"""Derive reproducible Vina boxes for the two thesis FtsZ pockets."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw" / "drive"
PROCESSED_DIR = ROOT / "data" / "processed"


POCKET_SOURCES = [
    {
        "pocket": "BP1",
        "source": "ezPocket_fpocket3.tsv / Pocket1",
        "center_x": -50.6447,
        "center_y": 35.0742,
        "center_z": 15.2025,
        "volume_a3": 669.1080,
        "thesis_residue_anchor": "Ala39B, Asn41B, Leu47B, Met49B, Ser50B, Lys55B",
        "selection_reason": "Nearest source pocket center to the BP1 interaction residues described in the original project notes.",
    },
    {
        "pocket": "BP2",
        "source": "ezPocket_fconv.tsv / Pocket0",
        "center_x": -63.4691,
        "center_y": 20.8178,
        "center_z": 1.73489,
        "volume_a3": 879.9980,
        "thesis_residue_anchor": "Met169B, Glu185B, Asn189B, Ile225B, Gly226B, Ser227B, Arg304B",
        "selection_reason": "Nearest source pocket center to the BP2 interaction residues described in the original project notes.",
    },
]


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for row in POCKET_SOURCES:
        cube_side_from_volume = float(row["volume_a3"]) ** (1 / 3)
        box_size = max(22.5, cube_side_from_volume + 12.0)
        rows.append(
            {
                **row,
                "cube_side_from_volume_a": round(cube_side_from_volume, 3),
                "size_x": round(box_size, 3),
                "size_y": round(box_size, 3),
                "size_z": round(box_size, 3),
                "box_note": "Estimated from source pocket center and volume; replace with original ezSMDock/Vina box if recovered.",
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(PROCESSED_DIR / "pocket_boxes.csv", index=False)

    config = {
        row["pocket"]: {
            "center": [row["center_x"], row["center_y"], row["center_z"]],
            "size": [row["size_x"], row["size_y"], row["size_z"]],
            "source": row["source"],
            "note": row["box_note"],
        }
        for row in rows
    }
    (PROCESSED_DIR / "pocket_boxes.json").write_text(json.dumps(config, indent=2) + "\n")

    print(table[["pocket", "source", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]].to_string(index=False))


if __name__ == "__main__":
    main()
