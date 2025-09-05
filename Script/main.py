# === Main driver script for computing longitudinal spin susceptibility χ^zz ===
# Supports AFM-Cr, Bcc-Fe, NM-Cr depending on chosen input.
# Includes data formatting, k-site determination, susceptibility calculation, and plotting.

from format import *
from determine_k_site import *
from solving_equation_mult import *
from thr_d_plot import *
from decay_plot import *
import time as t
import os

# === MATERIAL INPUT SETUP ===
# Choose dataset by changing the `url` path to point to desired material input
# Crucial to have this format as the header: 
# (Please copy and paste it into the top and align it correctly)
# i    j    dx   dy   dz   Jij         χ⁰↑             χ⁰↓

# url = r"Bcc-Fe\chfile-1.dat"         # Ferromagnetic Bcc Iron
# url = r"AFM-Cr\AFM-chfile.dat"       # Antiferromagnetic Chromium (2 sublattices)
url = r"NM-Cr\NM-chfile-1.dat"       # Nonmagnetic Chromium
# url = r"Multi_AFM-Cr\AFM-chfile.dat" # Multisite dependant Ferromagnetic Bcc Iron (Testing only)
# url = r"NA-BCC-NM\chfile-1.dat"      # Nonmagnetic BCC Sodium
# url = r"Cu-FCC\chfile-1.dat"         # FCC Copper
# url = r"TaSe\chfile-TaSe.dat"         

base_folder = os.path.dirname(url)

# === PHYSICAL PARAMETERS ===
radius = 5  # Cutoff radius for valid k-site neighbors
min = -10   # Grid range minimum
max = 10    # Grid range maximum

# === SITE SHIFT RULES FOR SUBLATTICES ===
# For AFM-Cr, sublattices are at (0,0,0) and (0.5,0.5,0.5) → shift between site types
# Defaults to, zero-shift if empty dict, or if site shift not given
# Input format: (i, j): (dx, dy, dz) shifting vector'

shift_rules_null = {

}

shift_rules_AFM = {
    (1, 2): (-0.5, -0.5, -0.5),
    (2, 1): (0.5, 0.5, 0.5),
}

shift_rules_TaSe = {
    (1, 2): (0.0, 0.0, 0.5),
    (2, 1): (0.0, 0.0, -0.5)
}

# === TRANSFORMATION MATRIX (depends on lattice symmetry and geometry) ===
base_change_matrix_NM = [
    [1.0, 0.0, 0.0],
    [-0.333333333333333, 0.942809033637852, 0.0],
    [-0.333333333333333, -0.471404564484194, 0.816496546528280]
]
base_change_matrix_AFM = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
]
base_change_matrix_Bcc = [
    [-0.5, 0.5, 0.5],
    [0.5, -0.5, 0.5],
    [0.5, 0.5, -0.5]
]
base_change_matrix_FCC = [
    [-0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [-0.5, 0.5, 0.0]
]

base_change_matrix_TaSe = [
    [0.866025403784439, -0.5, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 3.85960152170524]
]

# === STEP 1: Data Formatting with Optional Sublattice Shifting ===
time_start = t.time()
f_url = format_data(
    url, 
    output_file=os.path.join(base_folder, "formated_data.csv"),
    shift_map=shift_rules_null
)

# === STEP 2: Determine Valid K-sites in Radius ===
k_pot_url = det_K_pot(
    min, 
    max, 
    radius, 
    output_file=os.path.join(base_folder, "k_pot_coords.csv")
)

# === STEP 3: Find K-sites Valid for Each J-site Using Shifted Coordinates + Lattice Transform ===
k_suit_url = det_K_suit(
    f_url, 
    k_pot_url, 
    radius, 
    base_change_matrix_NM, 
    output_file=os.path.join(base_folder, "k_pot_coords.csv")
)

# === STEP 4: Group K by J into JSON (for input to Dyson solver) or CSV (debug view) ===
# Use one (recommended JSON for Reliability) of these depending on what your solver expects:
# k_match_url = det_K_match_csv(...)   # Easier to understand data format (Recommended for debugging)
k_match_url = det_K_match_json(f_url, k_suit_url, output_file=os.path.join(base_folder, "neighbouring_k_to_j.json"))

# === STEP 5: Define Kernel Parameters and Solve for χ^zz ===
# === ADVANCED: Example for future site-dependent U implementations (Multi-U KERNEL systems) ===
X_file = compute_Xzz_all_site_dependent(
    xij_file=f_url,
    kfile=k_match_url,
    site_map_file=f_url,
    U_params=[
        (2.0, 2.0, 0.0, 0.0),  # U matrix for site type 1  # U↑↑, U↓↓, U↑↓, U↓↑
        # (2.5, 2.5, 0.0, 0.0),  # U matrix for site type 2  # U↑↑, U↓↓, U↑↓, U↓↑
    ],
    output_file=os.path.join(base_folder, "xzz_output_iter.csv"),
    temp_txt=os.path.join(base_folder, "debug_iter.txt")
)

# === STEP 6: Merge Raw + Computed Data (useful for plotting/validation/export) ===
merge_format_and_xzz(f_url, X_file, output_file=os.path.join(base_folder, "formated_output.dat"))

time_end = t.time()
print(f"Data preparation and χ^zz computation done. Runtime: {time_end - time_start:.2f} s")

# === === === === === === === === === === === === === === === === ===
# === Optional: Enable Visualisation for Decay and Comparison Plots ===
# === === === === === === === === === === === === === === === === ===

# # AFM-Cr Setup (enabled by default)
# Length_scale = 5.32988627
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix=base_change_matrix_AFM,
#     scale_diagonal=[Length_scale] * 3,
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# Bcc-Fe Setup (enabled by default)
# Length_scale = 5.42 # Bcc - Fe Length scaling
# Length_scale = 7.98 # NA - BCC Length scaling
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix= base_change_matrix_Bcc,
#     scale_diagonal=[Length_scale] * 3,
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# NM-Cr Setup (enabled by default)
Length_scale = 4.640997226251
plot_static_and_spin_decay(
    static_file=f_url,
    spin_file=X_file,
    transform_matrix=base_change_matrix_NM,
    scale_diagonal=[Length_scale] * 3,
    output_static=os.path.join(base_folder, "static_decay_plot.png"),
    output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
    output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
)

# # Cu-FCC Setup (enabled by default)
# Length_scale = 6.760364108039488
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix= base_change_matrix_FCC,
#     scale_diagonal=[Length_scale] * 3,
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# # TaSe2 Setup (enabled by default)
# Length_scale = 6.760364108039488
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix= base_change_matrix_TaSe,
#     scale_diagonal=[Length_scale] * 3,
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# === OPTIONAL 3D - Plot ===
phys_plot(X_file)
time_end = t.time()
print(f"Plotting finished. Runtime: {time_end - time_start:.2f} s")


