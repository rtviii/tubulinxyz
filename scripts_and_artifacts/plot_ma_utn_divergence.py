#!/usr/bin/env python3
"""
Visualize where MA positions diverge from UTN consensus positions.

Produces two figures:
1. Global alignment overview: full MA with UTN-mapped vs orphan positions highlighted
2. Zoomed CTT region: the actual divergent positions with residue labels

Output: scripts_and_artifacts/ma_utn_divergence_alpha.png
        scripts_and_artifacts/ma_utn_divergence_beta.png
"""

import sys
sys.path.insert(0, ".")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from Bio import AlignIO

from lib.etl.utn_mapper import get_mapper
from lib.etl.sequence_alignment import ConsensusCalculator
from api.config import settings


def get_msa_data(family: str):
    """Get MSA sequences, consensus, and orphan positions."""
    mapper = get_mapper(family)

    if family == "tubulin_alpha":
        msa_path = settings.PROJECT_ROOT / "data" / "alpha_tubulin" / "alpha_tubulin.afasta"
    else:
        msa_path = settings.PROJECT_ROOT / "data" / "beta_tubulin" / "beta_tubulin.afasta"

    alignment = AlignIO.read(str(msa_path), "fasta")
    consensus = ConsensusCalculator(msa_path)
    ma_length = consensus.length

    # Classify each MA position
    orphan_positions = []
    mapped_positions = []
    for ma_pos in range(1, ma_length + 1):
        utn = mapper.master_to_utn_map.get(ma_pos)
        if utn is None:
            orphan_positions.append(ma_pos)
        else:
            mapped_positions.append(ma_pos)

    return {
        "family": family,
        "alignment": alignment,
        "consensus": consensus,
        "ma_length": ma_length,
        "utn_length": mapper._stats["utn_length"],
        "orphan_positions": orphan_positions,
        "mapped_positions": mapped_positions,
        "mapper": mapper,
    }


def plot_family(data: dict, output_path: str):
    """Create the two-panel figure for one family."""
    family = data["family"]
    alignment = data["alignment"]
    consensus = data["consensus"]
    ma_length = data["ma_length"]
    orphans = data["orphan_positions"]
    mapper = data["mapper"]

    family_label = family.replace("tubulin_", "").capitalize()
    n_seqs = len(alignment)
    seq_names = [r.id.split("|")[0] for r in alignment]

    fig, axes = plt.subplots(2, 1, figsize=(18, 10), gridspec_kw={"height_ratios": [3, 4]})
    fig.suptitle(
        f"{family_label}-tubulin: MA positions vs UTN consensus mapping\n"
        f"MA length = {ma_length}, UTN length = {data['utn_length']}, "
        f"{len(orphans)} MA-only positions (no UTN equivalent)",
        fontsize=13, fontweight="bold", y=0.98,
    )

    # =========================================================================
    # Panel 1: Global overview - horizontal bar showing mapped vs orphan
    # =========================================================================
    ax1 = axes[0]

    # Draw the full MA as a colored bar
    colors = []
    for pos in range(1, ma_length + 1):
        if pos in orphans:
            colors.append("#e74c3c")  # red for orphan
        else:
            colors.append("#2ecc71")  # green for mapped

    ax1.barh(0, width=1, left=range(ma_length), height=0.6, color=colors, edgecolor="none")

    # Mark orphan positions with vertical lines and labels
    for pos in orphans:
        ax1.axvline(x=pos - 0.5, color="#e74c3c", linewidth=0.8, alpha=0.7, linestyle="--")

    # Add tick marks every 50 positions
    tick_positions = list(range(0, ma_length, 50)) + [ma_length - 1]
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels([str(p + 1) for p in tick_positions], fontsize=8)
    ax1.set_yticks([])
    ax1.set_xlim(-1, ma_length + 1)
    ax1.set_xlabel("Master Alignment position", fontsize=10)
    ax1.set_title("Global overview: green = mapped to UTN, red = MA-only (no UTN equivalent)", fontsize=10)

    # Add region annotations
    ax1.annotate("N-terminus", xy=(10, 0.5), fontsize=8, color="gray", ha="left")
    if orphans:
        first_orphan = min(orphans)
        ax1.annotate(
            f"CTT region\n(orphans start at MA {first_orphan})",
            xy=(first_orphan - 0.5, 0.5),
            xytext=(first_orphan - 60, 0.9),
            fontsize=9, color="#c0392b", fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.5),
        )

    legend_mapped = mpatches.Patch(color="#2ecc71", label=f"Mapped to UTN ({len(data['mapped_positions'])} positions)")
    legend_orphan = mpatches.Patch(color="#e74c3c", label=f"MA-only, no UTN ({len(orphans)} positions)")
    ax1.legend(handles=[legend_mapped, legend_orphan], loc="upper left", fontsize=9)

    # =========================================================================
    # Panel 2: Zoomed CTT - show actual residues per isotype at orphan positions
    # =========================================================================
    ax2 = axes[1]

    # Define the zoom window: from ~20 positions before first orphan to end
    if orphans:
        zoom_start = max(1, min(orphans) - 20)
    else:
        zoom_start = max(1, ma_length - 30)
    zoom_end = ma_length
    zoom_range = list(range(zoom_start, zoom_end + 1))
    n_cols = len(zoom_range)

    # Build a residue grid: rows = sequences, columns = zoom positions
    residue_grid = []
    for record in alignment:
        row = []
        for ma_pos in zoom_range:
            idx = ma_pos - 1  # 0-based
            if idx < len(record.seq):
                row.append(str(record.seq[idx]))
            else:
                row.append("-")
        residue_grid.append(row)

    # Also add the consensus row
    consensus_row = []
    for ma_pos in zoom_range:
        consensus_row.append(consensus.get_residue(ma_pos - 1))

    # Also add the UTN consensus row (via mapper)
    utn_row = []
    from lib.etl.utn_mapper import _FAMILY_CONFIG
    utn_consensus_str = _FAMILY_CONFIG[family][0]
    for ma_pos in zoom_range:
        utn_pos = mapper.master_to_utn_map.get(ma_pos)
        if utn_pos is not None and 1 <= utn_pos <= len(utn_consensus_str):
            utn_row.append(utn_consensus_str[utn_pos - 1])
        else:
            utn_row.append("-")

    all_rows = residue_grid + [consensus_row, utn_row]
    all_labels = seq_names + ["Our consensus", "UTN consensus"]
    n_rows = len(all_rows)

    # Color coding
    for row_idx, (label, row) in enumerate(zip(all_labels, all_rows)):
        for col_idx, (residue, ma_pos) in enumerate(zip(row, zoom_range)):
            is_orphan = ma_pos in orphans

            if residue in ("-", "."):
                color = "#f5f5f5" if not is_orphan else "#fde8e8"
                text_color = "#cccccc"
            elif is_orphan:
                color = "#fadbd8"
                text_color = "#c0392b"
            elif row_idx >= n_seqs:  # consensus rows
                color = "#d5f5e3"
                text_color = "#1a5632"
            else:
                color = "#eafaf1"
                text_color = "#2c3e50"

            rect = plt.Rectangle(
                (col_idx, n_rows - 1 - row_idx), 1, 1,
                facecolor=color, edgecolor="#dddddd", linewidth=0.3,
            )
            ax2.add_patch(rect)
            ax2.text(
                col_idx + 0.5, n_rows - 1 - row_idx + 0.5, residue,
                ha="center", va="center", fontsize=7, fontfamily="monospace",
                color=text_color, fontweight="bold" if row_idx >= n_seqs else "normal",
            )

    # Axis setup
    ax2.set_xlim(0, n_cols)
    ax2.set_ylim(0, n_rows)
    ax2.set_xticks([i + 0.5 for i in range(n_cols)])
    ax2.set_xticklabels([str(p) for p in zoom_range], fontsize=6, rotation=90)
    ax2.set_yticks([i + 0.5 for i in range(n_rows)])
    ax2.set_yticklabels(list(reversed(all_labels)), fontsize=8)

    # Highlight orphan columns with red border
    for ma_pos in orphans:
        if ma_pos in zoom_range:
            col_idx = zoom_range.index(ma_pos)
            rect = plt.Rectangle(
                (col_idx, 0), 1, n_rows,
                facecolor="none", edgecolor="#e74c3c", linewidth=2,
            )
            ax2.add_patch(rect)

    # Draw separator line between sequences and consensus rows
    ax2.axhline(y=2, color="black", linewidth=1.5)

    ax2.set_title(
        f"Zoomed view: MA positions {zoom_start}-{zoom_end} (C-terminal region)",
        fontsize=10,
    )
    ax2.set_xlabel("Master Alignment position", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output_path}")


def main():
    out_dir = Path("scripts_and_artifacts")

    for family in ["tubulin_alpha", "tubulin_beta"]:
        print(f"\nProcessing {family}...")
        data = get_msa_data(family)
        label = family.replace("tubulin_", "")
        output = str(out_dir / f"ma_utn_divergence_{label}.png")
        plot_family(data, output)

        # Print summary
        print(f"  Orphan MA positions: {data['orphan_positions']}")
        for pos in data["orphan_positions"]:
            utn = data["mapper"].master_to_utn_map.get(pos)
            cons = data["consensus"].get_residue(pos - 1)
            print(f"    MA {pos}: consensus={cons}, UTN={utn}")


if __name__ == "__main__":
    main()
