import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from typing import Dict, List, Any
from mpl_toolkits.mplot3d import Axes3D


class VisualizationUtils:
    """Handles all visualization tasks"""

    def __init__(self):
        self.ALPHA_COLOR = "#1f77b4"
        self.BETA_COLOR = "#ff7f0e"
        self.INTERFACE_ALPHA_COLOR = "lightblue"
        self.INTERFACE_BETA_COLOR = "lightcoral"
        self.N_TERMINUS_COLOR = "green"

    def create_interface_visualization(
        self,
        positions: Dict[str, np.ndarray],
        interfaces: Dict[str, Dict[str, List]],
        chain_directions: Dict[str, Any],
        tubulin_chains: Dict[str, str],
        pdb_id: str,
    ):
        """Create 3D visualization showing positions, interfaces, and N-terminus"""
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(111, projection="3d")

        # Plot chain positions
        for chain_id, pos in positions.items():
            if chain_id in tubulin_chains:
                color = (
                    self.ALPHA_COLOR
                    if tubulin_chains[chain_id] == "alpha"
                    else self.BETA_COLOR
                )
                ax.scatter(pos[0], pos[1], pos[2], c=color, s=100, alpha=0.8)
                ax.text(pos[0], pos[1], pos[2], chain_id, fontsize=8)

        # Plot interface residues
        for chain_id, chain_interfaces in interfaces.items():
            if chain_id not in positions:
                continue

            chain_pos = positions[chain_id]
            chain_type = tubulin_chains.get(chain_id, "unknown")
            interface_color = (
                self.INTERFACE_ALPHA_COLOR
                if chain_type == "alpha"
                else self.INTERFACE_BETA_COLOR
            )

            for neighbor_id, interface_residues in chain_interfaces.items():
                if interface_residues:
                    # Calculate interface center
                    interface_coords = []
                    for residue in interface_residues:
                        for atom in residue:
                            interface_coords.append(atom.get_coord())

                    if interface_coords:
                        interface_center = np.mean(interface_coords, axis=0)
                        # Plot small circles for interface residues
                        ax.scatter(
                            interface_center[0],
                            interface_center[1],
                            interface_center[2],
                            c=interface_color,
                            s=30,
                            alpha=0.6,
                            marker="o",
                        )

                        # Draw line between chain center and interface center
                        ax.plot(
                            [chain_pos[0], interface_center[0]],
                            [chain_pos[1], interface_center[1]],
                            [chain_pos[2], interface_center[2]],
                            color=interface_color,
                            alpha=0.3,
                            linewidth=1,
                        )

        # Plot N-terminus residues as green triangles
        for chain_id, direction_data in chain_directions.items():
            if chain_id in positions and "n_terminus_residue" in direction_data:
                n_terminus = direction_data["n_terminus_residue"]
                if n_terminus:
                    # Get N-terminus position
                    n_term_coords = []
                    for atom in n_terminus:
                        n_term_coords.append(atom.get_coord())

                    if n_term_coords:
                        n_term_pos = np.mean(n_term_coords, axis=0)
                        ax.scatter(
                            n_term_pos[0],
                            n_term_pos[1],
                            n_term_pos[2],
                            c=self.N_TERMINUS_COLOR,
                            s=50,
                            alpha=0.9,
                            marker="^",
                        )

        # Set labels and title
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.set_title(
            f"Interface Analysis - {pdb_id.upper()}\nChain Positions, Interfaces, and N-terminus"
        )

        # Add legend
        legend_elements = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.ALPHA_COLOR,
                markersize=10,
                label="α-Tubulin",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.BETA_COLOR,
                markersize=10,
                label="β-Tubulin",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.INTERFACE_ALPHA_COLOR,
                markersize=6,
                label="α Interface",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.INTERFACE_BETA_COLOR,
                markersize=6,
                label="β Interface",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor=self.N_TERMINUS_COLOR,
                markersize=8,
                label="N-terminus",
            ),
        ]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        interface_viz_filename = f"interface_analysis_{pdb_id}.png"
        plt.savefig(interface_viz_filename, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"🎨 Saved interface visualization: {interface_viz_filename}")

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

        plt.tight_layout(rect=[0, 0, 0.9, 1])
        grid_viz_filename = f"grid_{pdb_id}_{grid_data.structure_type}.png"
        plt.savefig(grid_viz_filename, dpi=150)
        plt.close()
        print(f"🎨 Saved grid visualization: {grid_viz_filename}")

    def create_orientation_plot(
        self,
        positions: Dict[str, np.ndarray],
        orientations: Dict[str, np.ndarray],
        tubulin_chains: Dict[str, str],
        interfaces: Dict[str, Dict[str, List]],
        chain_directions: Dict[str, Any],
        pdb_id: str,
    ):
        """Create 3D orientation plot with interface and N-terminus information"""
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(111, projection="3d")

        # Plot chain positions and orientations
        for chain_id, pos in positions.items():
            if chain_id in tubulin_chains and chain_id in orientations:
                color = (
                    self.ALPHA_COLOR
                    if tubulin_chains[chain_id] == "alpha"
                    else self.BETA_COLOR
                )
                orientation = orientations[chain_id]

                # Plot chain position
                ax.scatter(pos[0], pos[1], pos[2], c=color, s=100, alpha=0.8)

                # Plot orientation vector
                end_pos = pos + orientation * 20  # Scale for visibility
                ax.quiver(
                    pos[0],
                    pos[1],
                    pos[2],
                    orientation[0],
                    orientation[1],
                    orientation[2],
                    color=color,
                    alpha=0.7,
                    arrow_length_ratio=0.1,
                    length=20,
                )

                # Add chain label
                ax.text(pos[0], pos[1], pos[2], chain_id, fontsize=8)

        # Plot interface residues as small circles
        for chain_id, chain_interfaces in interfaces.items():
            if chain_id not in positions:
                continue

            chain_type = tubulin_chains.get(chain_id, "unknown")
            interface_color = (
                self.INTERFACE_ALPHA_COLOR
                if chain_type == "alpha"
                else self.INTERFACE_BETA_COLOR
            )

            for neighbor_id, interface_residues in chain_interfaces.items():
                if interface_residues:
                    # Calculate interface center
                    interface_coords = []
                    for residue in interface_residues:
                        for atom in residue:
                            interface_coords.append(atom.get_coord())

                    if interface_coords:
                        interface_center = np.mean(interface_coords, axis=0)
                        ax.scatter(
                            interface_center[0],
                            interface_center[1],
                            interface_center[2],
                            c=interface_color,
                            s=20,
                            alpha=0.6,
                            marker="o",
                        )

        # Plot N-terminus residues as green triangles
        for chain_id, direction_data in chain_directions.items():
            if chain_id in positions and "n_terminus_residue" in direction_data:
                n_terminus = direction_data["n_terminus_residue"]
                if n_terminus:
                    # Get N-terminus position
                    n_term_coords = []
                    for atom in n_terminus:
                        n_term_coords.append(atom.get_coord())

                    if n_term_coords:
                        n_term_pos = np.mean(n_term_coords, axis=0)
                        ax.scatter(
                            n_term_pos[0],
                            n_term_pos[1],
                            n_term_pos[2],
                            c=self.N_TERMINUS_COLOR,
                            s=50,
                            alpha=0.9,
                            marker="^",
                        )

        # Set labels and title
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.set_title(f"Subunit Orientations with Interfaces - {pdb_id.upper()}")

        # Add legend
        legend_elements = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.ALPHA_COLOR,
                markersize=10,
                label="α-Tubulin",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.BETA_COLOR,
                markersize=10,
                label="β-Tubulin",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.INTERFACE_ALPHA_COLOR,
                markersize=6,
                label="α Interface",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=self.INTERFACE_BETA_COLOR,
                markersize=6,
                label="β Interface",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="^",
                color="w",
                markerfacecolor=self.N_TERMINUS_COLOR,
                markersize=8,
                label="N-terminus",
            ),
        ]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        orientation_viz_filename = f"orientations_{pdb_id}.png"
        plt.savefig(orientation_viz_filename, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"🎨 Saved orientation visualization: {orientation_viz_filename}")

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
            pf_color = colors[pf_idx]

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
                    label=f"PF{pf_idx}",
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

        # Add legend
        if len(protofilaments) <= 10:  # Only show legend if not too many protofilaments
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        pf_viz_filename = f"protofilaments_{pdb_id}.png"
        plt.savefig(pf_viz_filename, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"🎨 Saved protofilament tracing visualization: {pf_viz_filename}")

    def create_debug_visualization(self, data: Dict[str, Any], pdb_id: str, title: str):
        """Create generic debug visualization"""
        fig, ax = plt.subplots(figsize=(12, 8))

        # Simple text-based debug info
        debug_text = f"Debug Info for {pdb_id.upper()}\n\n"
        for key, value in data.items():
            if isinstance(value, (int, float, str)):
                debug_text += f"{key}: {value}\n"
            elif isinstance(value, (list, dict)):
                debug_text += f"{key}: {len(value)} items\n"

        ax.text(
            0.1,
            0.9,
            debug_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            fontfamily="monospace",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(title)
        ax.axis("off")

        debug_filename = f"debug_{pdb_id}_{title.lower().replace(' ', '_')}.png"
        plt.savefig(debug_filename, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"🎨 Saved debug visualization: {debug_filename}")
