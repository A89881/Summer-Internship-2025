from format import *
from determine_k_site import *
from solving_equation import *
import time as t

url = r"Data\chfile-1.csv"
radius = 5
min = -10
max = 10
time_start = t.time()
f_url = format_data(url)
k_pot_url = det_K_pot(min, max, radius)
k_suit_url = det_K_suit(f_url, k_pot_url, radius)
k_match_url = det_K_match(f_url, k_suit_url)
u_kernel_list = [2.0, 0.0, 0.0, 2.0]
results = compute_all_chi_zz(k_match_url, f_url, u_kernel_list)
time_end = t.time()
print(f"run time {time_end-time_start}") 