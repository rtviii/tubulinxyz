import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Dict, List, Any
from mpl_toolkits.mplot3d import Axes3D
import json
from pathlib import Path


class VisualizationUtils:
    """Handles all visualization tasks"""

    def __init__(self):
        self.ALPHA_COLOR = "#1f77b4"
        self.BETA_COLOR = "#ff7f0e"
        self.N_TERMINUS_COLOR = "green"
        self.debug_dir = Path("debug_output")
        self.debug_dir.mkdir(exist_ok=True)

    def create_pymol_script(self, pdb_id: str):
        """Create PyMOL script to visualize N-terminus analysis"""
        nterm_file = self.debug_dir / "nterm_debug.json"
        if not nterm_file.exists():
            print("Warning: No N-terminus debug data found")
            return

        with open(nterm_file) as f:
            data = json.load(f)

        script_content = f"""# PyMOL script for {pdb_id.upper()} N-terminus analysis
# Load structure (download first: wget https://files.rcsb.org/download/{pdb_id.upper()}.cif)
load {pdb_id.upper()}.cif
hide everything
show cartoon
color grey, all

# Color N-terminus residues in green
"""

        # Add N-terminus residues
        for chain_id, chain_data in data["n_terminus_data"].items():
            resnum = chain_data["n_terminus_resnum"]
            script_content += (
                f"select nterm_{chain_id}, chain {chain_id} and resi {resnum}\n"
            )
            script_content += f"color green, nterm_{chain_id}\n"
            script_content += f"show spheres, nterm_{chain_id}\n"

        script_content += f"\n# Show connections as distances\n"
        for source, target in data["connections"].items():
            script_content += (
                f"distance conn_{source}_{target}, nterm_{source}, nterm_{target}\n"
            )

        script_content += f"""
# Color connected pairs in red
"""

        for source, target in data["connections"].items():
            script_content += f"color red, chain {source} or chain {target}\n"

        script_content += f"""
# Summary: {len(data['connections'])} connections from {data['total_chains']} chains
# Success rate: {len(data['connections'])/data['total_chains']*100:.1f}%

zoom all
"""

        script_file = self.debug_dir / f"{pdb_id}_visualize.pml"
        with open(script_file, "w") as f:
            f.write(script_content)

        # Create batch script
        batch_script = f"""#!/bin/bash
# Download and visualize {pdb_id.upper()}
echo "Downloading {pdb_id.upper()}.cif..."
wget -O {pdb_id.upper()}.cif https://files.rcsb.org/download/{pdb_id.upper()}.cif
echo "Starting PyMOL..."
pymol {pdb_id.upper()}.cif -d "@{pdb_id}_visualize.pml"
"""

        batch_file = self.debug_dir / f"view_{pdb_id}.sh"
        with open(batch_file, "w") as f:
            f.write(batch_script)
        batch_file.chmod(0o755)

    def create_grid_visualization(self, grid_data, pdb_id: str):
        """Create 2D grid visualization"""
        if not grid_data.subunits:
            return

        max_pf = max((s.protofilament for s in grid_data.subunits), default=0)
        max_su = max((s.subunitIndex for s in grid_data.subunits), default=0)

        fig_width = max(8, (max_pf + 1) * 1.0)
        fig_height = max(6, (max_su + 1) * 1.0)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Plot subunits
        for s in grid_data.subunits:
            x, y = s.protofilament, s.subunitIndex
            color = self.ALPHA_COLOR if s.monomerType == "alpha" else self.BETA_COLOR
            ax.add_patch(patches.Circle((x, y), 0.45, fc=color, ec="black", lw=0.5))
            ax.text(
                x,
                y,
                s.auth_asym_id,
                ha="center",
                va="center",
                fontsize=6,
                weight="normal",
                color="white",
            )

        # Set up the plot
        ax.set_aspect("equal")
        ax.set_xlim(-0.7, max_pf + 0.7)
        ax.set_ylim(-0.7, max_su + 0.7)
        ax.set_xticks(range(max_pf + 1))
        ax.set_xticklabels([f"PF{i}" for i in range(max_pf + 1)])
        ax.set_yticks(range(max_su + 1))
        ax.set_xlabel("Protofilament Index", fontsize=12)
        ax.set_ylabel("Subunit Index", fontsize=12)
        ax.invert_yaxis()
        ax.set_title(
            f"Grid Layout - {pdb_id.upper()} ({grid_data.structure_type})",
            fontsize=14,
            weight="bold",
        )

        # Add legend
        alpha_patch = patches.Patch(color=self.ALPHA_COLOR, label="α-Tubulin")
        beta_patch = patches.Patch(color=self.BETA_COLOR, label="β-Tubulin")
        ax.legend(
            handles=[alpha_patch, beta_patch],
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )

        # Add grid
        ax.grid(True, linestyle=":", alpha=0.6)

        # plt.tight_layout(rect=[0, 0, 0.9, 1])
        # grid_viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        # plt.savefig(grid_viz_filename, dpi=150)
        # plt.close()

    def create_protofilament_tracing_plot(
        self,
        positions: Dict[str, np.ndarray],
        protofilaments: List[List[str]],
        tubulin_chains: Dict[str, str],
        pdb_id: str,
    ):
        """Create visualization showing traced protofilaments"""
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(111, projection="3d")

        # Color palette for different protofilaments
        if len(protofilaments) > 0:
            colors = plt.cm.Set3(np.linspace(0, 1, len(protofilaments)))
        else:
            colors = []

        # Plot each protofilament
        for pf_idx, pf_chains in enumerate(protofilaments):
            pf_color = colors[pf_idx] if pf_idx < len(colors) else "gray"

            # Get positions for this protofilament
            pf_positions = []
            for chain_id in pf_chains:
                if chain_id in positions:
                    pf_positions.append(positions[chain_id])

            if len(pf_positions) > 1:
                pf_positions = np.array(pf_positions)

                # Plot the protofilament path
                ax.plot(
                    pf_positions[:, 0],
                    pf_positions[:, 1],
                    pf_positions[:, 2],
                    color=pf_color,
                    linewidth=3,
                    alpha=0.7,
                    label=f"PF{pf_idx} ({len(pf_chains)})",
                )

                # Plot individual chains
                for i, chain_id in enumerate(pf_chains):
                    if chain_id in positions:
                        pos = positions[chain_id]
                        chain_type = tubulin_chains.get(chain_id, "unknown")
                        marker = "o" if chain_type == "alpha" else "s"
                        ax.scatter(
                            pos[0],
                            pos[1],
                            pos[2],
                            color=pf_color,
                            s=100,
                            marker=marker,
                            alpha=0.8,
                            edgecolors="black",
                        )
                        ax.text(pos[0], pos[1], pos[2], chain_id, fontsize=6)

        # Set labels and title
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.set_title(f"Protofilament Tracing - {pdb_id.upper()}")

        # Add legend if not too many protofilaments
        if len(protofilaments) <= 15:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        # plt.tight_layout()
        # pf_viz_filename = f"protofilaments_{pdb_id}.png"
        # plt.savefig(pf_viz_filename, dpi=150, bbox_inches="tight")
        # plt.close()

        # Create PyMOL script after creating the plot
        self.create_pymol_script(pdb_id)