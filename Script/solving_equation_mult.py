import numpy as np
import pandas as pd
import ast
import json
from typing import Dict, List, Tuple

# === 1. Parse bare response Xij file ===
def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float],
                                           Dict[Tuple[float, float, float], float]]:
    """
    Parses spin-resolved bare response (χ⁰↑, χ⁰↓) into dictionaries keyed by relative position (dx, dy, dz).
    """
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

# === 2. Parse j→k map from JSON or CSV ===
def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    """
    Parses j → [k₁, k₂, ...] contribution map from a CSV or JSON file.
    """
    def to_float_tuple(t) -> Tuple[float, float, float]:
        return tuple(float(x) for x in t) # type: ignore

    if kfile_path.endswith(".json"):
        with open(kfile_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
        return {
            to_float_tuple(ast.literal_eval(j_str)): [to_float_tuple(k) for k in k_list]
            for j_str, k_list in json_data.items()
        }
    else:
        contrib_map = {}
        with open(kfile_path, 'r') as f:
            next(f)  # skip header
            for line in f:
                j_str, k_str_list = line.strip().split(';')
                j_coord = to_float_tuple(ast.literal_eval(j_str))
                k_coords = [to_float_tuple(k) for k in ast.literal_eval(k_str_list)]
                contrib_map[j_coord] = k_coords
        return contrib_map

# === 3. Get spin-resolved block for a relative position ===
def get_x_block(rel: Tuple[float, float, float],
                x_map_up: Dict,
                x_map_down: Dict) -> Tuple[float, float]:
    """
    Gets χ⁰↑, χ⁰↓ at a relative vector `rel`. Returns (0.0, 0.0) if not present.
    """
    if rel == (0.0, 0.0, 0.0):
        return 0.0, 0.0
    return x_map_up.get(rel, 0.0), x_map_down.get(rel, 0.0)

# === 4. Direct Dyson solver ===
def compute_Kij_direct(Xij: np.ndarray,
                       Xik_list: List[np.ndarray],
                       U: np.ndarray) -> np.ndarray:
    """
    Solves Dyson-like equation: Kij = (I - S·U)⁻¹ · Xij where S = Σ Xik.
    """
    S = sum(Xik_list)
    identity = np.eye(2)
    A = identity - S @ U
    try:
        return np.linalg.solve(A, Xij)
    except np.linalg.LinAlgError:
        raise RuntimeError("Matrix inversion failed. Possibly singular or ill-conditioned.")

# === 5. Main driver for site-dependent kernel ===
def compute_Xzz_all_site_dependent(xij_file: str,
                                   kfile: str,
                                   site_map_file: str,
                                   U_params: List[Tuple[float, float, float, float]],
                                   output_file: str,
                                   temp_txt: str) -> str:
    """
    Computes χ^zz for each j-site with site-dependent U kernels.

    Args:
        xij_file: Formatted file with dx,dy,dz,χ⁰↑,χ⁰↓.
        kfile: CSV or JSON mapping j → [k₁, k₂, ...].
        site_map_file: CSV file with i,j,dx,dy,dz (to get site index of j).
        U_params: List of 4-tuple kernel matrices, one per site type.
        output_file: Output CSV file.
        temp_txt: Debug output path.
    """
    x_map_up, x_map_down = parse_xij_file(xij_file)
    contrib_map = parse_k_contrib_file(kfile)

    # Load site index for each j
    site_df = pd.read_csv(site_map_file, sep=';')
    j_site_type = {
        (float(row['dx']), float(row['dy']), float(row['dz'])): int(row['j'])
        for _, row in site_df.iterrows()
    }

    results = []
    with open(temp_txt, 'w', encoding='utf-8') as debug_out:
        for j_idx, (j_coord, k_coords) in enumerate(contrib_map.items()):
            if not k_coords:
                continue

            try:
                site_index = j_site_type.get(j_coord, 0)
                U_vals = U_params[site_index - 1]  # 1-based site index

                U = np.array([
                    [U_vals[0], U_vals[2]],
                    [U_vals[3], U_vals[1]]
                ])

                xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
                Xij = np.diag([xij_up, xij_down])

                Xik_list = []
                for kq in k_coords:
                    rel = tuple(np.subtract(kq, (0.0, 0.0, 0.0)))
                    x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
                    Xik_list.append(np.diag([x_up, x_down]))

                Kij = compute_Kij_direct(Xij, Xik_list, U)
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})

                # Write first few debug entries
                if j_idx < 3:
                    debug_out.write(f"=== j = {j_coord} ===\n")
                    debug_out.write(f"Xij:\n{Xij}\n")
                    debug_out.write(f"Sum(Xik):\n{sum(Xik_list)}\n")
                    debug_out.write(f"Kij:\n{Kij}\n")
                    debug_out.write(f"χzz = {chi_zz}\n\n")

            except Exception as e:
                print(f"[ERROR] j = {j_coord}: {e}")
                continue

    pd.DataFrame(results).to_csv(output_file, index=False)
    print(f"Done: χzz data written to {output_file}")
    print(f"Debug output written to {temp_txt}")
    return output_file
