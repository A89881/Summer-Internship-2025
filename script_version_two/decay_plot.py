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
    output_static: str = r"data\static_decay_plot.png",
    output_xzz: str = r"data\xzz_decay_plot.png",
    output_comparison: str = r"data\comparison_decay_plot.png"
):
    """
    Generates two plots:
    1. Decay of average static bare response χ⁰ vs distance.
    2. Comparison of χ⁰ and χzz vs distance.
    
    Args:
        static_file: CSV with dx, dy, dz, χ⁰↑, χ⁰↓ (e.g., formated_data.csv).
        spin_file: CSV with j-coordinate and Xzz (e.g., xzz_output_iter.csv).
        transform_matrix: 3x3 transformation matrix for unit vectors.
        scale_diagonal: List of 3 scale factors for scaling coordinates.
        output_static: Output path for static-only plot.
        output_comparison: Output path for static vs χzz comparison plot.
    """
    # --- Load and prepare static file ---
    static_df = pd.read_csv(static_file, sep=';')
    static_df['coord'] = static_df.apply(lambda row: (row['dx'], row['dy'], row['dz']), axis=1)
    static_df['chi0'] = 0.5 * (static_df['χ⁰↑'] + static_df['χ⁰↓'])

    coords_static = np.array(static_df['coord'].tolist())
    chi0_vals = static_df['chi0'].values

    # --- Load and prepare spin file ---
    spin_df = pd.read_csv(spin_file)
    spin_df[['x', 'y', 'z']] = spin_df['j-coordinate'].apply(lambda s: pd.Series(ast.literal_eval(s)))
    coords_spin = spin_df[['x', 'y', 'z']].values
    xzz_vals = spin_df['Xzz'].values

    # --- Transformation matrices ---
    T = np.array(transform_matrix)
    S = np.diag(scale_diagonal)

    # Apply transformation
    coords_static_transformed = coords_static @ S @ T
    coords_spin_transformed = coords_spin @ S @ T

    # Compute radial distances
    dists_static = np.linalg.norm(coords_static_transformed, axis=1)
    dists_spin = np.linalg.norm(coords_spin_transformed, axis=1)

    # --- Plot 1: Static bare susceptibility decay ---
    sort_idx_static = np.argsort(dists_static)
    plt.figure()
    plt.plot(dists_static[sort_idx_static], chi0_vals[sort_idx_static], 'o-', color='blue', linewidth=1, markersize=2) # type: ignore
    plt.title("χ⁰ vs. Transformed Distance")
    plt.xlabel("Distance from Origin (Transformed)")
    plt.ylabel("χ⁰ (static avg)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_static, dpi=300)
    plt.close()
    print(f"Static decay plot saved to: {output_static}")
    
    # === Plot 2: χzz-only decay plot ===
    sort_idx_spin = np.argsort(dists_spin)

    plt.figure()
    plt.plot(
        dists_spin[sort_idx_spin],
        xzz_vals[sort_idx_spin], # type: ignore
        'o-', color='darkgreen', linewidth=1, markersize=2
    )
    plt.title("χᶻᶻ vs. Transformed Distance")
    plt.xlabel("Distance from Origin (Transformed)")
    plt.ylabel("χᶻᶻ (longitudinal)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_xzz, dpi=300)
    plt.show()
    print(f"χzz-only decay plot saved to: {output_xzz}")


    # --- Plot 3: Comparison plot ---
    sort_idx_spin = np.argsort(dists_spin)
    plt.figure()
    plt.plot(dists_static[sort_idx_static], chi0_vals[sort_idx_static], 'o-', label='χ⁰ (static avg)', color='blue', markersize=2, linewidth=1) # type: ignore
    plt.plot(dists_spin[sort_idx_spin], xzz_vals[sort_idx_spin], 'o-', label='χᶻᶻ (spin response)', color='darkred', markersize=2, linewidth=1) # type: ignore
    plt.title("Comparison of χ⁰ and χᶻᶻ vs. Transformed Distance")
    plt.xlabel("Distance from Origin (Transformed)")
    plt.ylabel("Response Value")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_comparison, dpi=300)
    plt.show()
    print(f"Comparison decay plot saved to: {output_comparison}")


plot_static_and_spin_decay(
    static_file="data/formatted_data.csv",
    spin_file="data/xzz_output_iter.csv",
    transform_matrix=[
        [-0.5, 0.5, 0.5],  # x-axis
        [0.5, -0.5, 0.5],  # y-axis
        [0.5, 0.5, -0.5]   # z-axis
    ],
    scale_diagonal=[5.42, 5.42, 5.42],  # Example lattice scaling in a₀
)
