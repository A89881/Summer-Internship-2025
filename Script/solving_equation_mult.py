import numpy as np
import pandas as pd
import ast
import json
from typing import Dict, List, Tuple

# === 1. Parse bare response Xij file with (i,j) included ===
def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[int,int,Tuple[float,float,float]], float],
                                           Dict[Tuple[int,int,Tuple[float,float,float]], float]]:
    """
    Parses bare response file where each row has:
    i ; j ; dx ; dy ; dz ; χ⁰↑ ; χ⁰↓
    Returns two dictionaries with keys (i,j,(dx,dy,dz)) -> χ⁰↑ / χ⁰↓
    """
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        key = (int(row['i']), int(row['j']), (float(row['dx']), float(row['dy']), float(row['dz'])))
        x_map_up[key] = row['χ⁰↑']
        x_map_down[key] = row['χ⁰↓']
    return x_map_up, x_map_down

# === 2. Parse j-to-k site mapping file (JSON only) ===
def parse_k_contrib_file(kfile_path: str):
    """
    Loads the mapping from each j-site (i,j,(dx,dy,dz)) to its contributing k-sites.
    Returns:
        Dictionary of the form:
        (i, j, (dx,dy,dz)) -> [(j_for_k1, (dx,dy,dz)), ...]
    """
    def to_float_tuple(t):
        return tuple(float(x) for x in t)

    if not kfile_path.endswith(".json"):
        raise ValueError("Currently only JSON format supported for kfile.")

    with open(kfile_path, 'r', encoding='utf-8') as f:
        json_data = json.load(f)

    contrib_map = {}
    for j_str, k_list in json_data.items():
        # j_str is "(i,j,dx,dy,dz)"
        j_tuple = ast.literal_eval(j_str)
        i_val, j_val, dx, dy, dz = j_tuple
        j_key = (int(i_val), int(j_val), (float(dx), float(dy), float(dz)))

        k_entries = []
        for k_entry in k_list:
            j_for_k, kx, ky, kz = k_entry
            k_entries.append((int(j_for_k), (float(kx), float(ky), float(kz))))

        contrib_map[j_key] = k_entries
    return contrib_map

# === 3. Lookup χ⁰↑ and χ⁰↓ for given (i,j,(dx,dy,dz)) ===
def get_x_block(rel_key: Tuple[int,int,Tuple[float,float,float]],
                x_map_up: Dict,
                x_map_down: Dict) -> Tuple[float, float]:
    """
    Returns χ⁰↑ and χ⁰↓ for a given (i,j,(dx,dy,dz)) key.
    Returns (0,0) if not found.
    """
    return x_map_up.get(rel_key, 0.0), x_map_down.get(rel_key, 0.0)

# === 4. Dyson-like equation solver ===
def compute_Kij_direct(Xij: np.ndarray,
                       Xik_list: List[np.ndarray],
                       U: np.ndarray) -> np.ndarray:
    """
    Solves the Dyson-like equation:
        Kij = (I - Σ Xik · U)^(-1) · Xij
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
    Computes longitudinal spin susceptibility χ^zz for each j-site with site-dependent U kernels.
    Output CSV columns: i ; j ; j-coordinate ; N_k ; Xzz
    """
    # Load bare response functions
    x_map_up, x_map_down = parse_xij_file(xij_file)

    # Load j-to-k contribution map
    contrib_map = parse_k_contrib_file(kfile)

    # Load j-site types from site_map_file (dx,dy,dz -> site type)
    site_df = pd.read_csv(site_map_file, sep=';')
    j_site_type = {
        (float(row['dx']), float(row['dy']), float(row['dz'])): int(row['j'])
        for _, row in site_df.iterrows()
    }

    results = []
    with open(temp_txt, 'w', encoding='utf-8') as debug_out:
        for j_idx, (j_key, k_entries) in enumerate(contrib_map.items()):
            try:
                i_val, j_val, j_coord = j_key

                if not k_entries:
                    continue

                # Determine site type and corresponding U matrix
                site_index = j_site_type.get(j_coord, 1)  # fallback
                site_index = min(site_index, len(U_params))
                U_vals = U_params[site_index - 1]
                U = np.array([[U_vals[0], U_vals[2]],
                              [U_vals[3], U_vals[1]]])

                # Get static Xij for current j
                Xij_up, Xij_down = get_x_block((i_val, j_val, j_coord), x_map_up, x_map_down)
                Xij = np.diag([Xij_up, Xij_down])

                # Accumulate contributions from each k-site
                Xik_list = []
                for j_for_k, k_coord in k_entries:
                    key_k = (i_val, j_for_k, k_coord)
                    x_up, x_down = get_x_block(key_k, x_map_up, x_map_down)
                    Xik_list.append(np.diag([x_up, x_down]))

                # Solve Dyson-like equation
                Kij = compute_Kij_direct(Xij, Xik_list, U)

                # Compute χ^zz
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0

                # Append results
                results.append({
                    'i': i_val,
                    'j': j_val,
                    'j-coordinate': str(j_coord),
                    'N_k': len(k_entries),
                    'Xzz': chi_zz
                })

                # Debug output for first few sites
                if j_idx < 3:
                    debug_out.write(f"=== j = {j_coord} (i={i_val}, j={j_val}) ===\n")
                    debug_out.write(f"Xij:\n{Xij}\n")
                    debug_out.write(f"Sum(Xik):\n{sum(Xik_list)}\n")
                    debug_out.write(f"Kij:\n{Kij}\n")
                    debug_out.write(f"χzz = {chi_zz}\n\n")

            except Exception as e:
                print(f"[ERROR] j = {j_key}: {e}")
                continue

    # Write output CSV
    pd.DataFrame(results).to_csv(output_file, sep=';', index=False)
    print(f"[compute_Xzz] χzz data written to: {output_file}")
    print(f"[compute_Xzz] Debug output written to: {temp_txt}")
    return output_file
