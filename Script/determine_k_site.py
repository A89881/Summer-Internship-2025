import pandas as pd
import ast
import json
import numpy as np
from typing import Tuple, List

# ================================================================
# Utility functions
# ================================================================

def clean_tuple_str(t: Tuple[float, float, float]) -> str:
    """
    Convert a numeric 3-tuple into a clean string.
    Integers are printed without decimal, floats kept as floats.
    Example:
      (0.0, 1.5, 2.0) -> "(0, 1.5, 2)"
    """
    return "(" + ", ".join(str(int(x)) if float(x).is_integer() else str(float(x)) for x in t) + ")"

def _detect_column(df: pd.DataFrame, candidates: List[str], name: str) -> str:
    """
    Detect which column in a DataFrame matches a given logical field.

    Args:
        df (pd.DataFrame): input dataframe
        candidates (List[str]): possible column names
        name (str): logical label (e.g. 'i' or 'j') for error message

    Returns:
        str: column name found in df
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Could not find column for {name}. Tried {candidates}. Got {list(df.columns)}")


# ================================================================
# Main Step 1: Determine valid (j,k) pairs
# ================================================================

def det_K_suit(f_url: str,
               R: float,
               base_change: List[List[float]],
               output_file: str) -> str:
    """
    From the formatted site CSV (f_url), determine valid (j,k) pairs using:

    Conditions:
      - |k| <= R (norm in transformed basis)
      - |j - k| <= R (relative distance in transformed basis)
      - Same i-site index between j and k

    Args:
        f_url (str): input CSV path (formatted data file)
        R (float): cutoff radius
        base_change (List[List[float]]): 3x3 matrix mapping lattice basis to Cartesian
        output_file (str): output CSV path

    Returns:
        str: output CSV path

    Output CSV columns:
        i; j; j-coordinate; k_j; k-coordinate
    """
    # Read input CSV
    df = pd.read_csv(f_url, sep=";", index_col=False)

    # Detect required columns (robust to name variations)
    i_col = _detect_column(df, ["i", "site_i"], "i")
    j_col = _detect_column(df, ["j", "site_j"], "j")
    dx_col = _detect_column(df, ["dx"], "dx")
    dy_col = _detect_column(df, ["dy"], "dy")
    dz_col = _detect_column(df, ["dz"], "dz")

    # Lattice transformation matrix
    T = np.array(base_change, dtype=float)
    if T.shape != (3, 3):
        raise ValueError("base_change must be a 3x3 matrix")

    # Collect all sites as (i, j, coord)
    coords = []
    for _, row in df.iterrows():
        coords.append({
            "i": int(row[i_col]),
            "j": int(row[j_col]),
            "coord": np.array([float(row[dx_col]), float(row[dy_col]), float(row[dz_col])])
        })

    # Prepare results
    out_rows = []
    R2 = float(R) ** 2  # cutoff squared for efficiency

    # Loop over all (j, k) pairs
    for row_j in coords:
        j_vec = row_j["coord"]

        for row_k in coords:
            # Must belong to same i-site
            if row_j["i"] != row_k["i"]:
                continue

            k_vec = row_k["coord"]

            # Condition 1: |k| <= R (in transformed basis)
            if np.dot(k_vec @ T, k_vec @ T) > R2:
                continue

            # Condition 2: |j - k| <= R (relative distance in transformed basis)
            rel = j_vec - k_vec
            if np.dot(rel @ T, rel @ T) > R2:
                continue

            # Store valid (j, k) pair
            out_rows.append({
                "i": row_j["i"],
                "j": row_j["j"],
                "j-coordinate": clean_tuple_str(tuple(j_vec)),
                "k_j": row_k["j"],
                "k-coordinate": clean_tuple_str(tuple(k_vec))
            })

    # Save output CSV
    df_out = pd.DataFrame(out_rows, columns=["i", "j", "j-coordinate", "k_j", "k-coordinate"])
    df_out.to_csv(output_file, sep=";", index=False)
    print(f"[det_K_suit] Wrote {len(df_out)} rows written to: {output_file}")
    return output_file


# ================================================================
# Main Step 2: Convert CSV pairs to JSON mapping
# ================================================================

def det_K_match_json(det_k_csv: str, json_out: str) -> str:
    """
    Convert CSV from det_K_suit to JSON mapping of j-site → k-sites.

    JSON structure:
    {
      "(i, j, j_dx, j_dy, j_dz)": [
        [k_j, k_dx, k_dy, k_dz], ...
      ]
    }

    Args:
        det_k_csv (str): input CSV file (output of det_K_suit)
        json_out (str): output JSON path

    Returns:
        str: output JSON path
    """
    df = pd.read_csv(det_k_csv, sep=";")

    def parse_tuple(s: str) -> List[float]:
        """Parse tuple string '(x, y, z)' into list of floats"""
        return [float(x) for x in ast.literal_eval(s)]

    mapping = {}
    for _, row in df.iterrows():
        i_val = int(row["i"])
        j_val = int(row["j"])
        j_coord = parse_tuple(row["j-coordinate"])
        k_j = int(row["k_j"])
        k_coord = parse_tuple(row["k-coordinate"])

        # Unique key per (i, j, j-coord)
        key = f"({i_val}, {j_val}, {j_coord[0]}, {j_coord[1]}, {j_coord[2]})"

        # Append k-site data
        mapping.setdefault(key, []).append([k_j, k_coord[0], k_coord[1], k_coord[2]])

    # Save JSON
    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(mapping, f, indent=2)

    print(f"[det_K_match_json] Wrote {len(mapping)} keys written to: {json_out}")
    return json_out
