# scripts_and_artifacts/morisette_stuff/xlsx_to_csv.py
"""
Convert Morisette xlsx files to CSV and remove the originals.

Usage:
    python scripts_and_artifacts/morisette_stuff/xlsx_to_csv.py
"""

from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[2] / "data"

# (xlsx glob pattern in data dir, output csv name)
CONVERSIONS = [
    ("beta_tubulin/*mutations*.xlsx",     "beta_tubulin/beta_tubulin_mutations.csv"),
    ("beta_tubulin/*modifications*.xlsx", "beta_tubulin/beta_tubulin_modifications.csv"),
    ("alpha_tubulin/*mutations*.xlsx",    "alpha_tubulin/alpha_tubulin_mutations.csv"),
    ("alpha_tubulin/*modifications*.xlsx","alpha_tubulin/alpha_tubulin_modifications.csv"),
]


def convert():
    for glob_pat, csv_name in CONVERSIONS:
        xlsx_files = list(DATA_DIR.glob(glob_pat))
        if not xlsx_files:
            continue

        xlsx_path = xlsx_files[0]
        csv_path = DATA_DIR / csv_name

        if csv_path.exists():
            print(f"  CSV already exists, skipping: {csv_path.name}")
            continue

        print(f"  Converting: {xlsx_path.name} -> {csv_path.name}")
        df = pd.read_excel(xlsx_path, header=None)
        df.to_csv(csv_path, index=False, header=False)

        # Remove the xlsx
        xlsx_path.unlink()
        print(f"  Deleted: {xlsx_path.name}")


if __name__ == "__main__":
    print("Converting Morisette xlsx files to CSV...")
    convert()
    print("Done.")
