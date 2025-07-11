from format import *
from determine_k_site import *
from solving_equation import *
from thr_d_plot import *
from decay_plot import *
import time as t
import os

#Any file type as long as column seperated by space and not any other divisor
url = r"AFM-Cr\AFM-chfile.dat"
base_folder = os.path.dirname(url)
radius = 5
min = -10
max = 10
time_start = t.time()
f_url = format_data(url, output_file=os.path.join(base_folder, "formated_data.csv"))
k_pot_url = det_K_pot(min, max, radius, output_file=os.path.join(base_folder, "k_pot_coords.csv"))
k_suit_url = det_K_suit(f_url, k_pot_url, radius, output_file=os.path.join(base_folder, "k_pot_coords.csv"))
k_match_url = det_K_match(f_url, k_suit_url, output_file=os.path.join(base_folder, "neighbouring_k_to_j.csv"))
U_kernel = (2.0, 2.0, 0.0, 0.0)  # U↑↑, U↓↓, U↑↓, U↓↑
X_file = compute_Xzz_all(f_url, k_match_url, U_kernel, output_file=os.path.join(base_folder, "xzz_output_iter.csv"), temp_txt=os.path.join(base_folder, "debug_iter.txt"))
merge_format_and_xzz(f_url, X_file, output_file=os.path.join(base_folder, "formated_output.dat"))
time_end = t.time()
print(f"run time {time_end-time_start}") 
time_start = t.time()
phys_plot(X_file)
Length_scale = 5.32988627
plot_static_and_spin_decay(
    static_file=f_url,
    spin_file=X_file,
    transform_matrix=[
        [1.0, 0.0, 0.0],  # x-axis (row vector of x-unit vector)
        [0.0, 1.0, 0.0],  # y-axis
        [0.0, 0.0, 1.0]   # z-axis
    ],
    scale_diagonal=[Length_scale, Length_scale, Length_scale],  # Example lattice scaling in a₀
    output_static=os.path.join(base_folder, "static_decay_plot.png"),
    output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
    output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
)
time_end = t.time()
print(f"run time {time_end-time_start}") 
