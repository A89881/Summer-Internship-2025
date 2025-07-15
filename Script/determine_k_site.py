import pandas as pd
import ast  # For safely converting string representations to Python tuples
import json  # For saving results to JSON format
from typing import Tuple
import numpy as np

# === Step 1: Generate potential k-site displacements within a radius ===
def det_K_pot(min: int, max: int, R: int, output_file: str) -> str:
    """
    Computes all possible integer vectors (k_dx, k_dy, k_dz) such that
    their squared Euclidean norm is >0 and <= R^2.

    Args:
        min: Minimum coordinate value for each axis.
        max: Maximum coordinate value for each axis.
        R: Radius cutoff.
        output_file: Path to save the potential k-site vectors.
    Returns:
        Path to saved CSV.
    """
    k_x_vals, k_y_vals, k_z_vals = [], [], []
    for i in range(min, max + 1):
        for j in range(min, max + 1):
            for k in range(min, max + 1):
                if 0 < i**2 + j**2 + k**2 <= R**2:
                    k_x_vals.append(i)
                    k_y_vals.append(j)
                    k_z_vals.append(k)

    df = pd.DataFrame(list(zip(k_x_vals, k_y_vals, k_z_vals)), columns=["k_dx", "k_dy", "k_dz"])
    df.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file} (Result)")
    return output_file

# === Utility: Convert numpy values to standard Python tuple types ===
def to_tuple(x):
    return tuple(float(i) if isinstance(i, np.floating) else int(i) for i in x)

# === Utility: Format tuples cleanly for output as strings ===
def clean_tuple_str(t: Tuple[float, float, float]) -> str:
    """
    Converts a tuple to a string representation without trailing .0 for integers.
    E.g., (2.0, -1.0, 0.5) -> "(2, -1, 0.5)"
    """
    return "(" + ", ".join(str(int(x)) if x == int(x) else str(x) for x in t) + ")"

# === Step 2: Determine valid k-site neighbors for each j-site using base transformation ===
def det_K_suit(f_url: str,
               k_url: str,
               R: int,
               base_change: list,
               output_file: str) -> str:
    """
    Given a list of j-coordinates and potential k-vectors, determine which
    (j,k) pairs are valid based on a radius R in a transformed basis.

    Args:
        f_url: CSV path containing j-site dx,dy,dz.
        k_url: CSV path containing potential k displacements.
        R: Radius cutoff in transformed space.
        base_change: 3x3 matrix to transform relative (k-j) vectors.
        output_file: Output CSV with matching (j,k) pairs.

    Returns:
        Path to saved CSV file with columns ['j-coordinate', 'k-coordinate'].
    """
    df_j = pd.read_csv(f_url, sep=";", index_col=False)
    df_k = pd.read_csv(k_url, sep=";", index_col=False)

    j_coords = list(zip(df_j["dx"], df_j["dy"], df_j["dz"]))
    k_coords = list(zip(df_k["k_dx"], df_k["k_dy"], df_k["k_dz"]))
    T = np.array(base_change)  # Transformation matrix

    matched_j, matched_k = [], []

    for j in j_coords:
        j_vec = np.array(j, dtype=float)
        for k in k_coords:
            k_vec = np.array(k, dtype=float)
            rel_vec = k_vec - j_vec
            rel_transformed = rel_vec @ T
            dist2 = np.dot(rel_transformed, rel_transformed)
            if 0 < dist2 <= R**2:
                matched_j.append(to_tuple(j_vec))
                matched_k.append(to_tuple(k_vec))

    df_out = pd.DataFrame({
        "j-coordinate": [clean_tuple_str(j) for j in matched_j],
        "k-coordinate": [clean_tuple_str(k) for k in matched_k]
    })

    df_out.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file} (Clean tuple output)")
    return output_file

# === Step 3A: Create JSON mapping of j-sites → list of k-sites ===
def det_K_match_json(f_url: str, k_url: str, output_file: str) -> str:
    """
    From matched (j,k) CSV, generate a JSON file mapping each j-site to a list of contributing k-sites.

    Args:
        f_url: Formatted CSV with j-site dx,dy,dz columns.
        k_url: Matched CSV with 'j-coordinate' and 'k-coordinate' columns.
        output_file: Path to save JSON file.

    Returns:
        Path to output JSON.
    """
    df_formatted = pd.read_csv(f_url, sep=";", index_col=False)
    df_matches = pd.read_csv(k_url, sep=";", index_col=False)

    # Parse string tuples into actual tuple objects
    df_matches['j-tuple'] = df_matches['j-coordinate'].apply(lambda s: tuple(map(float, ast.literal_eval(s))))
    df_matches['k-tuple'] = df_matches['k-coordinate'].apply(lambda s: tuple(map(float, ast.literal_eval(s))))

    # Ensure k-tuples are present in formatted j-coordinates (validity check)
    known_j_coords = set(tuple(row) for row in df_formatted[['dx', 'dy', 'dz']].values)
    df_valid = df_matches[df_matches['k-tuple'].isin(known_j_coords)]

    # Build mapping
    grouped_dict = {}
    for _, row in df_valid.iterrows():
        j = row['j-tuple']
        k = row['k-tuple']
        grouped_dict.setdefault(j, []).append(k)

    # Convert keys and values to JSON-safe formats
    json_ready = {str(j): [list(k) for k in k_list] for j, k_list in grouped_dict.items()}

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(json_ready, f, indent=2)

    print(f"Done: The string url is: {output_file} (Result)")
    return output_file

# === Step 3B: Create CSV mapping of j-sites → list of k-sites (same logic as JSON but saved differently) ===
def det_K_match_csv(f_url: str, k_url: str, output_file: str) -> str:
    """
    Generate a CSV file grouping each j-site with its matched k-site list.

    Args:
        f_url: Formatted CSV with dx,dy,dz columns.
        k_url: Matched (j,k) CSV.
        output_file: Path to save grouped CSV.

    Returns:
        Path to saved CSV.
    """
    df_j = pd.read_csv(f_url, sep=";", index_col=False)
    df_k = pd.read_csv(k_url, sep=";", index_col=False)

    j_set = set(tuple(row) for row in df_j[['dx', 'dy', 'dz']].values)

    df_k['k-tuple'] = df_k['k-coordinate'].apply(lambda x: tuple(map(float, ast.literal_eval(x))))
    df_k['j-tuple'] = df_k['j-coordinate'].apply(lambda x: tuple(map(float, ast.literal_eval(x))))

    df_k_filtered = df_k[df_k['k-tuple'].isin(j_set)]

    grouped = df_k_filtered.groupby('j-tuple')['k-tuple'].apply(list).reset_index()

    grouped.to_csv(output_file, sep=";", index=False, header=["j-coordinate", "k-coordinates"])
    print(f"Done: The string url is: {output_file} (Result)")
    return output_file
