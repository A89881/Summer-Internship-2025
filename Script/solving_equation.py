import numpy as np
import pandas as pd
import ast
import json
from typing import Dict, List, Tuple

# === 1. Parse bare response Xij file into spin-resolved maps ===
def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float],
                                           Dict[Tuple[float, float, float], float]]:
    """
    Parses the formatted Xij CSV file and returns two dictionaries:
    one for spin-up (χ⁰↑) and one for spin-down (χ⁰↓) values indexed by (dx,dy,dz).

    Args:
        xij_path: Path to the CSV containing dx, dy, dz, χ⁰↑, χ⁰↓.

    Returns:
        Tuple of dictionaries: x_map_up and x_map_down
    """
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

# === 2. Parse K-contribution map from JSON or CSV (j → [k₁, k₂, ...]) ===
def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    """
    Reads the j → k-site contribution mapping from JSON or CSV format.

    Returns:
        Dictionary mapping j-coordinates to a list of k-coordinates.
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
            next(f)  # Skip header
            for line in f:
                j_str, k_str_list = line.strip().split(';')
                j_coord = to_float_tuple(ast.literal_eval(j_str))
                k_coords = [to_float_tuple(k) for k in ast.literal_eval(k_str_list)]
                contrib_map[j_coord] = k_coords
        return contrib_map

# === 3. Get χ⁰↑ and χ⁰↓ for a given relative vector (dx,dy,dz) ===
def get_x_block(rel: Tuple[float, float, float],
                x_map_up: Dict,
                x_map_down: Dict) -> Tuple[float, float]:
    """
    Retrieves χ⁰↑ and χ⁰↓ values for a relative vector.

    Returns:
        Tuple (χ⁰↑, χ⁰↓), or (0.0, 0.0) if not found.
    """
    if rel == (0.0, 0.0, 0.0):
        return 0.0, 0.0
    x_up = x_map_up.get(rel, 0.0)
    x_down = x_map_down.get(rel, 0.0)
    return x_up, x_down

# === 4. Dyson-like solver: compute full response using matrix solve ===
def compute_Kij_direct(Xij: np.ndarray,
                       Xik_list: List[np.ndarray],
                       U: np.ndarray) -> np.ndarray:
    """
    Solves the Dyson-like equation using a direct matrix inversion:
        Kij = (I - S·U)^(-1) · Xij
    where S = sum of diagonal Xik matrices.

    Args:
        Xij: Diagonal 2x2 matrix from j-site.
        Xik_list: List of diagonal Xik matrices from k-sites.
        U: 2x2 exchange-correlation kernel.

    Returns:
        Full 2x2 Kij matrix.
    """
    S = sum(Xik_list)
    identity = np.eye(2)

    try:
        A = identity - S @ U
        Kij = np.linalg.solve(A, Xij)
        return Kij
    except np.linalg.LinAlgError:
        raise RuntimeError("Direct solver failed: matrix may be singular or ill-conditioned.")

# === 5. Compute χ^zz across all j-sites ===
def compute_Xzz_all(xij_file: str,
                    kfile: str,
                    U_params: Tuple[float, float, float, float],
                    output_file: str,
                    temp_txt: str) -> str:
    """
    Computes the longitudinal spin susceptibility χ^zz for each j-site.

    Args:
        xij_file: Formatted file with χ⁰↑, χ⁰↓ for each relative coordinate.
        kfile: JSON/CSV with contributing k-sites for each j-site.
        U_params: Tuple (U↑↑, U↓↓, U↑↓, U↓↑) in row-major format.
        output_file: CSV file to store χ^zz values.
        temp_txt: Debug file to store intermediate Kij, Xij info.

    Returns:
        Path to the output CSV with χ^zz values.
    """
    # Assemble exchange-correlation kernel matrix U
    U = np.array([
        [U_params[0], U_params[2]],
        [U_params[3], U_params[1]]
    ])

    x_map_up, x_map_down = parse_xij_file(xij_file)
    contrib_map = parse_k_contrib_file(kfile_path=kfile)
    results = []

    with open(temp_txt, 'w', encoding='utf-8') as debug_out:
        for j_idx, (j_coord, k_coords) in enumerate(contrib_map.items()):
            if not k_coords:
                continue
            try:
                # Get static Xij for the j-site
                xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
                Xij = np.diag([xij_up, xij_down])

                # Construct sum over k contributions
                Xik_list = []
                for kq in k_coords:
                    rel = tuple(np.subtract(kq, (0.0, 0.0, 0.0)))  # Absolute in this implementation
                    x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
                    Xik_list.append(np.diag([x_up, x_down]))

                # Solve Dyson-like equation
                Kij = compute_Kij_direct(Xij, Xik_list, U)

                # Compute χ^zz = (K↑↑ + K↓↓ - K↑↓ - K↓↑)/4
                chi_zz = (Kij[0, 0] + Kij[1, 1] - Kij[0, 1] - Kij[1, 0]) / 4.0
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})

                # Write debug information for the first 3 j-sites
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
