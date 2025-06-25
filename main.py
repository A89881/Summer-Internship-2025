import csv
import pandas as pd
import numpy as np
from format import *
from determine_k_site import *

url = r"Data\raw-data-set.csv"
radius = 5
min = -10
max = 10
f_url = format_data(url)
k_url = det_K_pot(min, max, radius)
k2_url = det_K_suit(f_url, k_url, radius)
det_K_match(f_url, k2_url)





