import numpy as np
import pandas as pd
import ast
import json
from typing import Dict, List, Tuple

# === 1. Parse bare response Xij file ===
def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float],
                                           Dict[Tuple[float, float, float], float]]:
    """
    Parses the formatted bare response file (e.g. formated_data.csv), where each row contains:
    - A relative coordinate (dx, dy, dz)
    - The corresponding spin-up (χ⁰↑) and spin-down (χ⁰↓) components

    Returns:
        Two dictionaries mapping (dx, dy, dz) to χ⁰↑ and χ⁰↓ respectively.
    """
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

# === 2. Parse j-to-k site mapping file (from JSON or CSV) ===
def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    """
    Loads the mapping from each j-site to its contributing k-sites.

    Returns:
        Dictionary of the form:
        (j_dx, j_dy, j_dz) → [(k1_dx, k1_dy, k1_dz), (k2_dx, ...), ...]
    """
    def to_float_tuple(t) -> Tuple[float, float, float]:
        return tuple(float(x) for x in t)  # type: ignore # converts stringified tuples to float tuples

    if kfile_path.endswith(".json"):
        with open(kfile_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
        return {
            to_float_tuple(ast.literal_eval(j_str)): [to_float_tuple(k) for k in k_list]
            for j_str, k_list in json_data.items()
        }
    else:  # assume CSV
        contrib_map = {}
        with open(kfile_path, 'r') as f:
            next(f)  # skip header
            for line in f:
                j_str, k_str_list = line.strip().split(';')
                j_coord = to_float_tuple(ast.literal_eval(j_str))
                k_coords = [to_float_tuple(k) for k in ast.literal_eval(k_str_list)]
                contrib_map[j_coord] = k_coords
        return contrib_map

# === 3. Lookup χ⁰↑ and χ⁰↓ block for a given relative displacement ===
def get_x_block(rel: Tuple[float, float, float],
                x_map_up: Dict,
                x_map_down: Dict) -> Tuple[float, float]:
    """
    Returns the spin-resolved χ⁰ values for a relative position.
    If rel is (0,0,0) or not found in the data, it returns (0.0, 0.0).
    """
    if rel == (0.0, 0.0, 0.0):
        return 0.0, 0.0
    return x_map_up.get(rel, 0.0), x_map_down.get(rel, 0.0)

# === 4. Dyson-like equation solver ===
def compute_Kij_direct(Xij: np.ndarray,
                       Xik_list: List[np.ndarray],
                       U: np.ndarray) -> np.ndarray:
    """
    Solves the Dyson-like equation:
        Kij = (I - Σ Xik · U)^(-1) · Xij

    Args:
        Xij: 2x2 matrix of static response between i and j.
        Xik_list: List of 2x2 matrices (for each k-contribution).
        U: 2x2 exchange-correlation kernel for site type of j.

    Returns:
        2x2 matrix representing Kij (full response).
    """
    S = sum(Xik_list)
    identity = np.eye(2)
    A = identity - S @ U
    try:
        return np.linalg.solve(A, Xij)
    except np.linalg.LinAlgError:
        raise RuntimeError("Matrix inversion failed. Possibly singular or ill-conditioned.")

# === 5. Full χ^zz evaluation with site-dependent U kernels ===
def compute_Xzz_all_site_dependent(xij_file: str,
                                   kfile: str,
                                   site_map_file: str,
                                   U_params: List[Tuple[float, float, float, float]],
                                   output_file: str,
                                   temp_txt: str) -> str:
    """
    Computes longitudinal spin susceptibility χ^zz for each j-site, using a different
    U kernel per site-type as indicated in site_map_file.

    Args:
        xij_file: Path to file containing χ⁰↑, χ⁰↓ and displacement vectors.
        kfile: JSON or CSV that maps each j-site to its contributing k-sites.
        site_map_file: CSV with dx, dy, dz and j (site type) for every j-site.
        U_params: List of U kernels (in tuple form), one for each site type.
        output_file: CSV output of j-coordinates and their computed χ^zz.
        temp_txt: Debug log path for inspecting internal matrix components.

    Returns:
        Path to the CSV output file.
    """
    # Load bare response functions
    x_map_up, x_map_down = parse_xij_file(xij_file)
    # Load j-to-k contribution map
    contrib_map = parse_k_contrib_file(kfile)
    # Load j-site types from format_file
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
                # === Identify site type and corresponding U matrix ===
                site_index = j_site_type.get(j_coord, 1)  # fallback to site type 1
                site_index = min(site_index, len(U_params))  # prevent index out-of-bounds
                U_vals = U_params[site_index - 1]  # convert 1-based to 0-based index
                U = np.array([
                    [U_vals[0], U_vals[2]],
                    [U_vals[3], U_vals[1]]
                ])

                # === Get static Xij for current j ===
                xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
                Xij = np.diag([xij_up, xij_down])

                # === Accumulate contributions from each k-site ===
                Xik_list = []
                for kq in k_coords:
                    rel = tuple(np.subtract(kq, (0.0, 0.0, 0.0)))  # shift is already applied before
                    x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
                    Xik_list.append(np.diag([x_up, x_down]))

                # === Dyson equation solver ===
                Kij = compute_Kij_direct(Xij, Xik_list, U)

                # === χ^zz = (K↑↑ + K↓↓ - K↑↓ - K↓↑)/4 ===
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})

                # === Optional debug info for first few sites ===
                if j_idx < 3:
                    debug_out.write(f"=== j = {j_coord} ===\n")
                    debug_out.write(f"Xij:\n{Xij}\n")
                    debug_out.write(f"Sum(Xik):\n{sum(Xik_list)}\n")
                    debug_out.write(f"Kij:\n{Kij}\n")
                    debug_out.write(f"χzz = {chi_zz}\n\n")

            except Exception as e:
                print(f"[ERROR] j = {j_coord}: {e}")
                continue

    # === Write results to output file ===
    pd.DataFrame(results).to_csv(output_file, index=False)
    print(f"Done: χzz data written to {output_file}")
    print(f"Debug output written to {temp_txt}")
    return output_file
