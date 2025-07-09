import numpy as np
import pandas as pd
import ast
from collections import defaultdict

"""
Version 0: Without txt output
"""

def parse_xij_file(xij_path):
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (row['dx'], row['dy'], row['dz'])
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

def parse_k_contrib_file(kfile_path):
    contrib_map = {}
    with open(kfile_path, 'r') as f:
        next(f)
        for line in f:
            j_str, k_str_list = line.strip().split(';')
            j_coord = ast.literal_eval(j_str)
            k_coords = ast.literal_eval(k_str_list)
            contrib_map[j_coord] = k_coords
    return contrib_map

def construct_A_int(x_map_up, x_map_down, k_coords, U):
    N = len(k_coords)
    A = np.zeros((2*N, 2*N))  # [↑↑, ↑↓, ↓↑, ↓↓] block representation
    for i, kq in enumerate(k_coords):
        for j, kq_prime in enumerate(k_coords):
            x_up = x_map_up.get(tuple(np.subtract(kq, kq_prime)), 0)
            x_down = x_map_down.get(tuple(np.subtract(kq, kq_prime)), 0)

            A[2*i, 2*j]     = x_up * U[0]     # ↑↑
            A[2*i, 2*j+1]   = x_up * U[2]     # ↑↓
            A[2*i+1, 2*j]   = x_down * U[3]   # ↓↑
            A[2*i+1, 2*j+1] = x_down * U[1]   # ↓↓

    return A

def construct_X_int(x_map_up, x_map_down, k_coords):
    X = np.zeros((2*len(k_coords), 1))
    for i, kq in enumerate(k_coords):
        X[2*i]   = x_map_up.get((0, 0, 0), 0)   # Diagonal X^{↑}
        X[2*i+1] = x_map_down.get((0, 0, 0), 0) # Diagonal X^{↓}
    return X

def compute_Kij(x_map_up, x_map_down, j_coord, k_coords, K_int, U):
    xij_up   = x_map_up.get(j_coord, 0)
    xij_down = x_map_down.get(j_coord, 0)

    K_up = xij_up
    K_down = xij_down

    for q, kq in enumerate(k_coords):
        xik_up = x_map_up.get(tuple(np.subtract(j_coord, kq)), 0)
        xik_down = x_map_down.get(tuple(np.subtract(j_coord, kq)), 0)

        K_up += xik_up * (U[0]*K_int[2*q][0] + U[2]*K_int[2*q+1][0])
        K_down += xik_down * (U[3]*K_int[2*q][0] + U[1]*K_int[2*q+1][0])

    return K_up + K_down

def compute_Kzz_all(xij_file, kfile, U, output_file='Kzz_output.csv'):
    x_map_up, x_map_down = parse_xij_file(xij_file)
    contrib_map = parse_k_contrib_file(kfile)

    output_rows = []

    for j_coord, k_coords in contrib_map.items():
        if not k_coords:
            continue

        A_int = construct_A_int(x_map_up, x_map_down, k_coords, U)
        X_int = construct_X_int(x_map_up, x_map_down, k_coords)

        try:
            K_int = np.linalg.solve(np.eye(A_int.shape[0]) - A_int, X_int)
        except np.linalg.LinAlgError:
            print(f"Warning: Singular matrix for j = {j_coord}. Skipping.")
            continue

        Kzz = compute_Kij(x_map_up, x_map_down, j_coord, k_coords, K_int, U)
        output_rows.append({'j-coordinate': j_coord, 'Kzz': Kzz})

    pd.DataFrame(output_rows).to_csv(output_file, index=False)
    print(f"Saved Kzz results to: {output_file}")