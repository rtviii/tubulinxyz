import pandas as pd
import sys
from pathlib import Path

def main():
    csv_path = Path("tubulin_ligand_diversity.csv")
    if not csv_path.exists():
        print(f"Error: {csv_path} not found. Please run analyze_ligands.py first.")
        return

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # --- 1. Overview by Category ---
    print("\n" + "="*60)
    print("  TUBULIN LIGAND CENSUS: CATEGORY BREAKDOWN")
    print("="*60)

    # Group by Category and sum structure counts
    cat_stats = df.groupby('Category').agg({
        'Chemical_ID': 'count',
        'Structure_Count': 'sum'
    }).rename(columns={'Chemical_ID': 'Unique Ligands', 'Structure_Count': 'Total Occurrences'})

    print(cat_stats.sort_values('Total Occurrences', ascending=False).to_string())

    # --- 2. Top Interesting Ligands Visualization ---
    print("\n" + "="*80)
    print("  TOP 25 BIO-RELEVANT LIGANDS (Excluding Ions/Buffers)")
    print("  Values represent number of PDB structures containing the ligand.")
    print("="*80)

    # Filter out common buffers/ions
    mask = ~df['Category'].str.contains("Common/Ion/Buffer", na=False)
    # Also optionally filter out GTP/GDP if user wants purely exogenous things,
    # but they are 'Nucleotide Analog' or 'Common' in the CSV usually.
    # Let's trust the Category assignment from the previous script.

    interesting = df[mask].sort_values(by='Structure_Count', ascending=False).head(25)

    if interesting.empty:
        print("No interesting ligands found.")
        return

    max_val = interesting['Structure_Count'].max()
    bar_width = 40

    # Table Header
    # ID  | Count | Name                           | Bar
    header = f"{'ID':<5} | {'Count':<5} | {'Type':<15} | {'Name (Truncated)':<30} | Distribution"
    print(header)
    print("-" * len(header))

    for _, row in interesting.iterrows():
        chem_id = str(row['Chemical_ID'])
        count = int(row['Structure_Count'])
        name = str(row['Name'])
        cat = str(row['Category'])

        # Truncate name
        if len(name) > 28:
            name = name[:25] + "..."

        # Simplify Category for display
        if "Drug" in cat:
            cat_disp = "Drug/Lead"
        elif "Nucleotide" in cat:
            cat_disp = "Nucleotide"
        else:
            cat_disp = "Other"

        # Draw Bar
        bar_len = int((count / max_val) * bar_width)
        bar_str = "█" * bar_len
        if count > 0 and bar_len == 0:
            bar_str = "▏" # Tiny bar for small non-zero counts

        print(f"{chem_id:<5} | {count:<5} | {cat_disp:<15} | {name:<30} | {bar_str} {count}")

    print("-" * len(header))

    # --- 3. DrugBank Stats ---
    total_relevant = len(df[mask])
    with_db = len(df[mask & (df['Is_DrugBank'] == True)])
    print(f"\nMetadata Insight: {with_db} out of {total_relevant} bio-relevant ligands have DrugBank IDs linked.")

if __name__ == "__main__":
    main()
