import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap, BoundaryNorm
import ast
import numpy as np
import os

# === Load CSV data ===
file_path = r'data\xzz_output.csv'
df = pd.read_csv(file_path)

# === Parse j-coordinate strings to (x, y, z)
df[['x', 'y', 'z']] = df['j-coordinate'].apply(lambda s: pd.Series(ast.literal_eval(s)))

# === Symmetry check ===
mirrored = df.copy()
mirrored[['x', 'y', 'z']] *= -1
merged = pd.merge(df, mirrored, on=['x', 'y', 'z'], suffixes=('', '_mirrored'))
diff = (merged['Xzz'] - merged['Xzz_mirrored']).abs()
print("Max diff under inversion symmetry:", diff.max())

# === Output directory for 2D layer plots
output_dir = 'z_layer_plots'
os.makedirs(output_dir, exist_ok=True)

# === Discrete color setup
N_BINS = 9  # Number of discrete color steps
vmin = df['Xzz'].min()
vmax = df['Xzz'].max()
bounds = np.linspace(vmin, vmax, N_BINS + 1)
norm = BoundaryNorm(boundaries=bounds, ncolors=N_BINS)
cmap = plt.get_cmap('seismic', N_BINS)  # discrete seismic map

# === Plot 2D scatter plots layer by layer ===
z_layers = sorted(df['z'].unique())

for z_val in z_layers:
    layer_df = df[df['z'] == z_val]

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        layer_df['x'], layer_df['y'],
        c=layer_df['Xzz'],
        cmap=cmap,
        norm=norm,
        s=120, edgecolor='black', linewidth=0.5
    )

    ax.set_title(f'Xzz Scatter at z = {z_val}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)

    cbar = plt.colorbar(sc, ax=ax, pad=0.01, ticks=bounds)
    cbar.set_label('Xzz (discretized)')
    cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])  # scientific notation

    plt.tight_layout()
    output_path = os.path.join(output_dir, f'Xzz_z{z_val:.1f}.png')
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Saved: {output_path}")

# === 3D scatter with discrete colormap ===
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

cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, ticks=bounds)
cbar.set_label('Xzz (discretized)')
cbar.ax.set_yticklabels([f"{b:.1e}" for b in bounds])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z') # type: ignore
ax.set_title('3D Scatter Plot of Xzz (Discrete Colors)')
plt.tight_layout()
plt.show()
