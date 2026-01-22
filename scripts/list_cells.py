#!/usr/bin/env python3
"""List notebook cells with indices, IDs, and first lines.

Usage:
    python scripts/list_cells.py notebooks/2_integration.ipynb
    python scripts/list_cells.py notebooks/*.ipynb --add-ids  # Add IDs to old notebooks
"""
import json
import sys
import uuid
from pathlib import Path


def list_cells(nb_path: str, max_line_len: int = 50) -> None:
    """List all cells in a notebook with their indices and first lines."""
    with open(nb_path) as f:
        nb = json.load(f)

    print(f"\nNotebook: {nb_path}")
    print(f"nbformat: {nb.get('nbformat', '?')}.{nb.get('nbformat_minor', '?')}")
    print(f"Total cells: {len(nb['cells'])}")
    print()
    print(f"{'Idx':<4} | {'Cell ID':<12} | {'Type':<8} | First Line")
    print("-" * 90)

    for i, cell in enumerate(nb["cells"]):
        cell_id = cell.get("id", "-")[:12]
        cell_type = cell["cell_type"]
        source = cell["source"]

        if isinstance(source, list):
            first_line = source[0] if source else "(empty)"
        else:
            first_line = source.split("\n")[0] if source else "(empty)"

        first_line = first_line.replace("\n", "").strip()[:max_line_len]
        print(f"{i:<4} | {cell_id:<12} | {cell_type:<8} | {first_line}")


def add_cell_ids(nb_path: str) -> None:
    """Add cell IDs to notebooks that don't have them (upgrades to nbformat 4.5)."""
    with open(nb_path) as f:
        nb = json.load(f)

    modified = False
    for cell in nb["cells"]:
        if "id" not in cell:
            cell["id"] = str(uuid.uuid4())[:8]
            modified = True

    if modified:
        nb["nbformat_minor"] = 5
        with open(nb_path, "w") as f:
            json.dump(nb, f, indent=1)
        print(f"Added cell IDs to {nb_path}")
    else:
        print(f"All cells already have IDs in {nb_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    add_ids = "--add-ids" in sys.argv
    paths = [p for p in sys.argv[1:] if not p.startswith("--")]

    for path in paths:
        for nb_file in Path(".").glob(path) if "*" in path else [Path(path)]:
            if add_ids:
                add_cell_ids(str(nb_file))
            else:
                list_cells(str(nb_file))
