from format import *
from determine_k_site import *
from solving_equation import *
from thr_d_plot import *
from decay_plot import *
import time as t
import os

#Any file type as long as column seperated by space and not any other divisor
url = r"AFM-Cr\AFM-chfile.dat"
# url = r"Bcc-Fe\chfile-1.dat"
# url = r"NM-Cr\NM-chfile-1.dat"
base_folder = os.path.dirname(url)
radius = 5
min = -10
max = 10
shift_rules = {
    (2, 1): (0.5, 0.5, 0.5),
    (1, 2): (-0.5, -0.5, -0.5),
    # (2, 1): (0.0, 0.0, 0.0),
    # (1, 2): (0.0, 0.0, 0.0),
    (1, 1): (0.0, 0.0, 0.0),
    (2, 2): (0.0, 0.0, 0.0)
}
time_start = t.time()
f_url = format_data(url, output_file=os.path.join(base_folder, "formated_data.csv"))
k_pot_url = det_K_pot(min, max, radius, output_file=os.path.join(base_folder, "k_pot_coords.csv"))
k_suit_url = det_K_suit(f_url, k_pot_url, radius, output_file=os.path.join(base_folder, "k_pot_coords.csv"))

"""
det_k_match_json for calculation, but det_k_match for csv output (easier to look interpret)
"""
# k_match_url = det_K_match_csv(f_url, k_suit_url, output_file=os.path.join(base_folder, "neighbouring_k_to_j.csv"))  
k_match_url = det_K_match_json(f_url, k_suit_url, output_file=os.path.join(base_folder, "neighbouring_k_to_j.json"))  

U_kernel = (2.0, 2.0, 0.0, 0.0)  # U↑↑, U↓↓, U↑↓, U↓↑
X_file = compute_Xzz_all(f_url, k_match_url, U_kernel, output_file=os.path.join(base_folder, "xzz_output_iter.csv"), temp_txt=os.path.join(base_folder, "debug_iter.txt"))
merge_format_and_xzz(f_url, X_file, output_file=os.path.join(base_folder, "formated_output.dat"))
time_end = t.time()
print(f"run time {time_end-time_start}") 

#Graphical tools and visualisation
time_start = t.time()
# phys_plot(X_file)

"""AFM-Cr - PLOT DATA"""
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

"""Bcc-Fe - PLOT DATA"""
# Length_scale = 5.42
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix=[
#         [-0.5, 0.5, 0.5],  # x-axis (row vector of x-unit vector)
#         [0.5, -0.5, 0.5],  # y-axis
#         [0.5, 0.5, -0.5]   # z-axis
#     ],
#     scale_diagonal=[Length_scale, Length_scale, Length_scale],  # Example lattice scaling in a₀
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# time_end = t.time()
# print(f"run time {time_end-time_start}") 

"""NM-Cr - PLOT DATA"""
# Length_scale = 4.640997226251
# plot_static_and_spin_decay(
#     static_file=f_url,
#     spin_file=X_file,
#     transform_matrix=[
#         [1.0, 0.0, 0.0],  # x-axis (row vector of x-unit vector)
#         [-0.333333333333333, 0.942809033637852, 0.0],  # y-axis
#         [-0.333333333333333, -0.471404564484194, 0.816496546528280]   # z-axis
#     ],
#     scale_diagonal=[Length_scale, Length_scale, Length_scale],  # Example lattice scaling in a₀
#     output_static=os.path.join(base_folder, "static_decay_plot.png"),
#     output_xzz=os.path.join(base_folder, "xzz_decay_plot.png"),
#     output_comparison=os.path.join(base_folder, "comparison_decay_plot.png")
# )

# time_end = t.time()
# print(f"run time {time_end-time_start}") 

"""Generallized script run"""
# compute_Xzz_all_site_dependent(
#     xij_file='data/formated_data.csv',
#     kfile='data/k_match.csv',
#     U_params=[
#         (2.0, 2.0, 0.0, 0.0),  # U for site type 1
#         (1.5, 1.5, 0.0, 0.0)   # U for site type 2
#     ],
#     output_file='data/xzz_output.csv',
#     temp_txt='data/debug_log.txt'
# )
