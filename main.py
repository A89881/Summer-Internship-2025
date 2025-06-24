import csv
import pandas as pd
import numpy as np

df = pd.read_csv(r"Data\formatted-data.csv", sep=";", index_col=False)

# Determine the Range in Which i and j are to each other 
# Assume if i is origin also --> rij = (dx-0, dy-0, dz-0)
# for col in df.columns[2:5]:
#     print(min(df[col]))
#     print(max(df[col]))

# Determine possible K that are smaller than 25
def determine_K1(max, R):
    k_x_vals = []
    k_y_vals = []
    k_z_vals = []
    for i in range(0, max+1):
        for j in range(0, max+1):
            for k in range(0, max+1):
                if 0 < i**2+j**2+k**2 <= R**2:
                    k_x_vals.append(i)
                    k_y_vals.append(j)
                    k_z_vals.append(k)
    # return [k_x_vals, k_y_vals, k_z_vals]
    return k_x_vals, k_y_vals, k_z_vals
lst1, lst2, lst3 = determine_K1(10, 5)

df2 = pd.DataFrame(list(zip(*[lst1, lst2, lst3])))
df2.columns = ["k_dx", "k_dy", "k_dz"]
df2.to_csv(r"Data\k_pot_coord.csv", sep=";", index=False)



