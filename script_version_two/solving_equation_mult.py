from typing import Dict, List, Tuple, Union
import numpy as np
import pandas as pd
import ast


def parse_xij_file_with_sites(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float],
                                                      Dict[Tuple[float, float, float], float],
                                                      Dict[Tuple[float, float, float], int]]:
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    j_site_map = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = float(row['χ⁰↑'])
        x_map_down[coord] = float(row['χ⁰↓'])
        j_site_map[coord] = int(row['j']) if 'j' in df.columns else 1
    return x_map_up, x_map_down, j_site_map


def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    contrib_map = {}
    with open(kfile_path, 'r') as f:
        next(f)
        for line in f:
            j_str, k_str_list = line.strip().split(';')
            j_coord = ast.literal_eval(j_str)
            k_coords = ast.literal_eval(k_str_list)
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


def compute_Kij_site_dependent(Xij: np.ndarray,
                               XikU_blocks: List[np.ndarray],
                               max_iter: int = 100,
                               tol: float = 1e-8) -> np.ndarray:
    Kij = Xij.copy()
    for _ in range(max_iter):
        Kij_new = Xij + sum([XU @ Kij for XU in XikU_blocks])
        if np.linalg.norm(Kij_new - Kij, ord='fro') < tol:
            return Kij_new
        Kij = Kij_new
    raise RuntimeError("Kij iteration did not converge")


def compute_Xzz_all_site_dependent(xij_file: str,
                                   kfile: str,
                                   U_params: Union[Tuple[float, float, float, float],
                                                   List[Tuple[float, float, float, float]]],
                                   output_file: str,
                                   temp_txt: str):
    if isinstance(U_params[0], float) or isinstance(U_params[0], int):
        U_params = [U_params]  # treat as one U kernel

    U_matrices = [np.array([[u[0], u[2]], [u[3], u[1]]]) for u in U_params]

    x_map_up, x_map_down, site_map = parse_xij_file_with_sites(xij_file)
    contrib_map = parse_k_contrib_file(kfile)
    results = []

    with open(temp_txt, 'w', encoding='utf-8') as debug_out:
        for j_idx, (j_coord, k_coords) in enumerate(contrib_map.items()):
            if not k_coords:
                continue

            try:
                xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
                Xij = np.diag([xij_up, xij_down])

                XikU_blocks = []
                for kq in k_coords:
                    rel = tuple(np.subtract(kq, (0.0, 0.0, 0.0)))
                    x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
                    site_type = site_map.get(rel, 1)
                    U_site = U_matrices[site_type - 1 if site_type - 1 < len(U_matrices) else 0]
                    Xik = np.diag([x_up, x_down])
                    XikU_blocks.append(Xik @ U_site)

                Kij = compute_Kij_site_dependent(Xij, XikU_blocks)
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})

                if j_idx < 3:
                    debug_out.write(f"=== j = {j_coord} ===\n")
                    debug_out.write(f"Xij:\n{Xij}\n")
                    debug_out.write(f"Kij:\n{Kij}\n")
                    debug_out.write(f"χzz = {chi_zz}\n\n")

            except Exception as e:
                print(f"[ERROR] j = {j_coord}: {e}")
                continue

    pd.DataFrame(results).to_csv(output_file, index=False)
    print(f"Done: The string url is: {output_file} (Result)")
    print(f"Done: The string url is: {temp_txt} (Debug)")
    return output_file
