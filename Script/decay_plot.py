import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ast
from typing import List

def plot_static_and_spin_decay(
    static_file: str,
    spin_file: str,
    transform_matrix: List[List[float]],
    scale_diagonal: List[float],
    output_static: str,
    output_xzz: str,
    output_comparison: str
):
    """
    Generates three plots:
    1. χ⁰ (static bare susceptibility) vs transformed distance.
    2. χᶻᶻ (longitudinal spin susceptibility) vs transformed distance.
    3. Comparison of χ⁰ and χᶻᶻ on the same plot.

    Args:
        static_file: CSV containing dx, dy, dz, χ⁰↑, χ⁰↓ values (e.g., formated_data.csv or merged .dat file).
        spin_file: CSV containing 'j-coordinate' and computed χᶻᶻ values.
        transform_matrix: 3x3 matrix to transform lattice directions to Cartesian space.
        scale_diagonal: Real-space scaling vector [a₁, a₂, a₃] (e.g., lattice constants).
        output_static: Path to save plot of χ⁰ decay.
        output_xzz: Path to save plot of χᶻᶻ decay.
        output_comparison: Path to save comparison plot.
    """

    # === Load static bare susceptibility χ⁰ ===
    static_df = pd.read_csv(static_file, sep=';')
    static_df['coord'] = static_df.apply(lambda row: (float(row['dx']), float(row['dy']), float(row['dz'])), axis=1)
    static_df['chi0'] = 0.5 * (static_df['χ⁰↑'] + static_df['χ⁰↓'])  # Average up/down

    coords_static = np.array(static_df['coord'].tolist())
    chi0_vals = static_df['chi0'].values

    # === Load longitudinal spin susceptibility χᶻᶻ ===
    spin_df = pd.read_csv(spin_file, sep=';')
    # Convert 'j-coordinate' string to float tuple
    spin_df[['x', 'y', 'z']] = spin_df['j-coordinate'].apply(lambda s: pd.Series([float(v) for v in ast.literal_eval(s)]))
    coords_spin = spin_df[['x', 'y', 'z']].values
    xzz_vals = spin_df['Xzz'].values

    # === Apply transformation: lattice → real space ===
    T = np.array(transform_matrix)
    S = np.diag(scale_diagonal)

    coords_static_transformed = coords_static @ S @ T
    coords_spin_transformed = coords_spin @ S @ T

    # === Compute Euclidean distances from origin ===
    dists_static = np.linalg.norm(coords_static_transformed, axis=1)
    dists_spin = np.linalg.norm(coords_spin_transformed, axis=1)

    # === Plot 1: χ⁰ decay ===
    sort_idx_static = np.argsort(dists_static)
    plt.figure()
    plt.plot(
        dists_static[sort_idx_static],
        chi0_vals[sort_idx_static], # type: ignore
        'o-', color='blue', linewidth=1, markersize=2
    )
    plt.title("χ⁰ vs Transformed Distance")
    plt.xlabel("Distance (Transformed)")
    plt.ylabel("χ⁰ (avg)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_static, dpi=100)
    plt.close()
    print(f"[plot_static_and_decay] Static decay plot saved to: {output_static}")

    # === Plot 2: χᶻᶻ decay ===
    sort_idx_spin = np.argsort(dists_spin)
    plt.figure()
    plt.plot(
        dists_spin[sort_idx_spin],
        xzz_vals[sort_idx_spin], # type: ignore
        'o-', color='darkgreen', linewidth=1, markersize=2
    )
    plt.title("χᶻᶻ vs Transformed Distance")
    plt.xlabel("Distance (Transformed)")
    plt.ylabel("χᶻᶻ (longitudinal)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_xzz, dpi=100)
    plt.show()
    print(f"[plot_static_and_decay] LSS decay plot saved to: {output_xzz}")

    # === Plot 3: Comparison χ⁰ vs χᶻᶻ ===
    plt.figure()
    plt.plot(
        dists_static[sort_idx_static],
        chi0_vals[sort_idx_static], # type: ignore
        'o-', label='χ⁰ (avg)', color='blue', linewidth=1, markersize=2
    )
    plt.plot(
        dists_spin[sort_idx_spin],
        xzz_vals[sort_idx_spin], # type: ignore
        'o-', label='χᶻᶻ', color='darkred', linewidth=1, markersize=2
    )
    plt.title("Comparison of χ⁰ and χᶻᶻ vs Transformed Distance")
    plt.xlabel("Distance (Transformed)")
    plt.ylabel("Response Value")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_comparison, dpi=100)
    plt.show()
    print(f"[plot_static_and_decay] Comparison decay plot saved to: {output_comparison}")
