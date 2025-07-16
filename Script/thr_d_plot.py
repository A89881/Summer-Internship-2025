import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap, BoundaryNorm
import ast
import numpy as np
import os

# === Function: Visualize χᶻᶻ from a CSV file ===
def phys_plot(input_file: str):
    """
    Generates:
    1. 2D pixel heatmaps (Xzz per (x, y) layer at each z-slice)
    2. 3D scatter plot of all χᶻᶻ points
    Also checks symmetry under inversion (x, y, z) → (-x, -y, -z)

    Args:
        input_file: Path to the CSV containing columns:
            - 'j-coordinate' (as stringified 3-tuples)
            - 'Xzz' (spin susceptibility at that site)
    """

    # === Load and parse input CSV ===
    df = pd.read_csv(input_file)
    # Convert 'j-coordinate' strings to columns: x, y, z
    df[['x', 'y', 'z']] = df['j-coordinate'].apply(lambda s: pd.Series(ast.literal_eval(s)))

    # === Inversion symmetry check ===
    # For AFM Cr, Xzz(i) ≈ Xzz(–i) → check max difference
    mirrored = df.copy()
    mirrored[['x', 'y', 'z']] *= -1
    merged = pd.merge(df, mirrored, on=['x', 'y', 'z'], suffixes=('', '_mirrored'))
    diff = (merged['Xzz'] - merged['Xzz_mirrored']).abs()
    print("Max diff under inversion symmetry:", diff.max())

    # === Setup: Output folder for 2D z-layer plots ===
    output_dir = os.path.join(os.path.dirname(input_file), "z_layer_plots")
    os.makedirs(output_dir, exist_ok=True)

    # === Setup: Discrete colormap centered at 0 ===
    N_BINS = 11
    vmax = df['Xzz'].abs().max()
    vmin = -vmax  # Symmetric range about 0
    bounds = np.linspace(vmin, vmax, N_BINS + 1)
    norm = BoundaryNorm(boundaries=bounds, ncolors=N_BINS)
    cmap = plt.get_cmap('coolwarm', N_BINS)  # Perceptually clear diverging colormap

    # === 2D Heatmap per z-layer ===
    z_layers = sorted(df['z'].unique())

    for z_val in z_layers:
        layer_df = df[df['z'] == z_val]

        # Define spatial range
        x_vals = sorted(layer_df['x'].unique())
        y_vals = sorted(layer_df['y'].unique())

        # Grid indices
        x_index = {v: i for i, v in enumerate(x_vals)}
        y_index = {v: i for i, v in enumerate(y_vals)}
        grid = np.full((len(y_vals), len(x_vals)), np.nan)  # y rows × x cols

        # Fill grid with Xzz values
        for _, row in layer_df.iterrows():
            i = y_index[row['y']]
            j = x_index[row['x']]
            grid[i, j] = row['Xzz']

        # Create heatmap
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(
            grid,
            origin='lower',
            cmap=cmap,
            norm=norm,
            extent=[
                min(x_vals) - 0.5, max(x_vals) + 0.5,
                min(y_vals) - 0.5, max(y_vals) + 0.5
            ], # type: ignore
            interpolation='none',
            aspect='equal'
        )

        ax.set_title(f'Xzz Heatmap at z = {z_val}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, pad=0.01, ticks=bounds)
        cbar.set_label('Xzz (discretized)')
        cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])

        # Save each 2D plot
        plt.tight_layout()
        output_path = os.path.join(output_dir, f'Xzz_z{z_val:.1f}.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"Saved: {output_path}")

    # === Final 3D Scatter Plot ===
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(
        df['x'], df['y'], df['z'],
        c=df['Xzz'],
        cmap=cmap,
        norm=norm,
        s=30, # type: ignore
        edgecolor='black', linewidth=0.3
    )

    # Add 3D colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, ticks=bounds)
    cbar.set_label('Xzz (discretized)')
    cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z') # type: ignore
    ax.set_title('3D Scatter Plot of Xzz (Discrete Colors)')

    plt.tight_layout()
    plt.show()
