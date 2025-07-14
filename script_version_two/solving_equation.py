import numpy as np
import pandas as pd
import ast
import json
from typing import Dict, List, Tuple

def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float],
                                           Dict[Tuple[float, float, float], float]]:
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    """
    Parses j → [k₁, k₂, ...] mapping from JSON or CSV and ensures float-typed keys consistently.
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
            next(f)
            for line in f:
                j_str, k_str_list = line.strip().split(';')
                j_coord = to_float_tuple(ast.literal_eval(j_str))
                k_coords = [to_float_tuple(k) for k in ast.literal_eval(k_str_list)]
                contrib_map[j_coord] = k_coords
        return contrib_map

def get_x_block(rel: Tuple[float, float, float],
                x_map_up: Dict,
                x_map_down: Dict) -> Tuple[float, float]:
    if rel == (0.0, 0.0, 0.0):
        return 0.0, 0.0
    x_up = x_map_up.get(rel, 0.0)
    x_down = x_map_down.get(rel, 0.0)
    return x_up, x_down

def compute_Kij_iterative(Xij: np.ndarray,
                          Xik_list: List[np.ndarray],
                          U: np.ndarray,
                          max_iter: int = 1000,
                          tol: float = 1e-8) -> np.ndarray:
    S = sum(Xik_list)
    Kij = Xij.copy()
    for _ in range(max_iter):
        Kij_new = Xij + S @ U @ Kij
        if np.linalg.norm(Kij_new - Kij, ord='fro') < tol:
            return Kij_new
        Kij = Kij_new
    raise RuntimeError("Kij iteration did not converge.")

def compute_Xzz_all(xij_file: str,
                    kfile: str,
                    U_params: Tuple[float, float, float, float],
                    output_file: str,
                    temp_txt: str):
    U = np.array([
        [U_params[0], U_params[2]],
        [U_params[3], U_params[1]]
    ])

    x_map_up, x_map_down = parse_xij_file(xij_file)
    contrib_map = parse_k_contrib_file(kfile)
    results = []

    with open(temp_txt, 'w', encoding='utf-8') as debug_out:
        for j_idx, (j_coord, k_coords) in enumerate(contrib_map.items()):
            if not k_coords:
                continue
            try:
                xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
                Xij = np.diag([xij_up, xij_down])

                Xik_list = []
                
                 # Ensure deterministic k-order and safe float coordinates
                # k_coords = sorted([kq for kq in k_coords])
                for kq in k_coords:
                    rel = tuple(np.subtract(kq, (0.0, 0.0, 0.0)))
                    x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
                    Xik_list.append(np.diag([x_up, x_down]))

                Kij = compute_Kij_iterative(Xij, Xik_list, U)
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})

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
    print(f"Done: The string url is: {output_file} (Result)")
    print(f"Done: The string url is: {temp_txt} (Debug)")
    return output_file
