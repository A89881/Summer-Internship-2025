import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import BoundaryNorm
import ast
import numpy as np
import os

def phys_plot_by_i(input_file: str):
    """
    For each unique i:
        - Show interactive 3D scatter plot of Xzz for all sites of that i
        - Save 2D heatmaps per z-layer in a folder specific to that i
    """

    # === Load CSV with proper delimiter ===
    df = pd.read_csv(input_file, sep=';')

    # Ensure 'j-coordinate' column exists
    if 'j-coordinate' in df.columns:
        df[['x', 'y', 'z']] = df['j-coordinate'].apply(lambda s: pd.Series(ast.literal_eval(s)))
    else:
        raise ValueError("CSV must contain 'j-coordinate' column (semicolon-separated).")
    
    df_inv = df.copy()
    df_inv[['x', 'y', 'z']] *= -1

    merged = pd.merge(
        df, df_inv,
        on=['x', 'y', 'z'],
        suffixes=('', '_inv')
    )

    # Compute difference in Xzz
    merged['diff'] = np.abs(merged['Xzz'] - merged['Xzz_inv'])

    # Print summary
    print(f"Number of sites checked: {len(merged)}")
    print(f"Max difference under inversion symmetry: {merged['diff'].max():.6e}")
    print(f"Mean difference: {merged['diff'].mean():.6e}")

    # Ensure 'i' column exists
    if 'i' not in df.columns:
        df['i'] = 1
    df['i'] = df['i'].astype(int)

    # Discrete colormap
    N_BINS = 11
    vmax = df['Xzz'].abs().max()
    vmin = -vmax
    bounds = np.linspace(vmin, vmax, N_BINS + 1)
    norm = BoundaryNorm(boundaries=bounds, ncolors=N_BINS)
    cmap = plt.get_cmap('coolwarm', N_BINS)

    # Loop over unique i
    for i_val in sorted(df['i'].unique()):
        df_i = df[df['i'] == i_val]

        # === 3D scatter plot (interactive) ===
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        sc = ax.scatter(
            df_i['x'], df_i['y'], df_i['z'],
            c=df_i['Xzz'], cmap=cmap, norm=norm,
            s=30, edgecolor='black', linewidth=0.3 # type: ignore
        )
        cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, ticks=bounds)
        cbar.set_label('Xzz')
        cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])
        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z') # type: ignore
        ax.set_title(f'3D Scatter Plot of Xzz (i={i_val})')
        plt.show()  # interactive plot

        # === 2D heatmaps per z-layer ===
        output_dir = os.path.join(os.path.dirname(input_file), f"plots_i{i_val}")
        os.makedirs(output_dir, exist_ok=True)

        z_layers = sorted(df_i['z'].unique())
        for z_val in z_layers:
            layer_df = df_i[df_i['z'] == z_val]

            x_vals = sorted(layer_df['x'].unique())
            y_vals = sorted(layer_df['y'].unique())
            x_index = {v: idx for idx, v in enumerate(x_vals)}
            y_index = {v: idx for idx, v in enumerate(y_vals)}

            grid = np.full((len(y_vals), len(x_vals)), np.nan)
            for _, row in layer_df.iterrows():
                grid[y_index[row['y']], x_index[row['x']]] = row['Xzz']

            fig, ax = plt.subplots(figsize=(6, 5))
            im = ax.imshow(
                grid, origin='lower', cmap=cmap, norm=norm,
                extent=[min(x_vals)-0.5, max(x_vals)+0.5, min(y_vals)-0.5, max(y_vals)+0.5], # type: ignore
                interpolation='none', aspect='equal'
            )
            ax.set_title(f'Xzz Heatmap: i={i_val}, z={z_val}')
            ax.set_xlabel('x'); ax.set_ylabel('y')
            ax.grid(True, linestyle='--', alpha=0.6)
            cbar = plt.colorbar(im, ax=ax, pad=0.01, ticks=bounds)
            cbar.set_label('Xzz')
            cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])

            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'i{i_val}_z{z_val:.1f}.png'), dpi=300)
            plt.close()
            print(f"Saved 2D heatmap for i={i_val}, z={z_val}")
