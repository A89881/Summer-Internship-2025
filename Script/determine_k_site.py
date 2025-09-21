import pandas as pd
import ast
import json
import numpy as np
from typing import Tuple, List

# --- utilities ---
def clean_tuple_str(t: Tuple[float, float, float]) -> str:
    return "(" + ", ".join(str(int(x)) if float(x).is_integer() else str(float(x)) for x in t) + ")"

def _detect_column(df: pd.DataFrame, candidates: List[str], name: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Could not find column for {name}. Tried {candidates}. Got {list(df.columns)}")

# --- main step 1 ---
def det_K_suit(f_url: str,
               R: float,
               base_change: List[List[float]],
               output_file: str) -> str:
    """
    From the formatted site CSV (f_url), determine valid (j,k) pairs using:
      - |k| <= R (after base_change)
      - |j - k| <= R (after base_change)
      - same i-site between j and k

    Output CSV columns: i; j; j-coordinate; k_j; k-coordinate
    """
    df = pd.read_csv(f_url, sep=";", index_col=False)

    i_col = _detect_column(df, ["i", "site_i"], "i")
    j_col = _detect_column(df, ["j", "site_j"], "j")
    dx_col = _detect_column(df, ["dx"], "dx")
    dy_col = _detect_column(df, ["dy"], "dy")
    dz_col = _detect_column(df, ["dz"], "dz")

    T = np.array(base_change, dtype=float)
    if T.shape != (3, 3):
        raise ValueError("base_change must be a 3x3 matrix")

    coords = []
    for _, row in df.iterrows():
        coords.append({
            "i": int(row[i_col]),
            "j": int(row[j_col]),
            "coord": np.array([float(row[dx_col]), float(row[dy_col]), float(row[dz_col])])
        })

    out_rows = []
    R2 = float(R) ** 2

    for row_j in coords:
        j_vec = row_j["coord"]
        for row_k in coords:
            if row_j["i"] != row_k["i"]:
                continue

            k_vec = row_k["coord"]

            # check |k| <= R
            if np.dot(k_vec @ T, k_vec @ T) > R2:
                continue

            # check |j - k| <= R
            rel = j_vec - k_vec
            if np.dot(rel @ T, rel @ T) > R2:
                continue

            out_rows.append({
                "i": row_j["i"],
                "j": row_j["j"],
                "j-coordinate": clean_tuple_str(tuple(j_vec)),
                "k_j": row_k["j"],
                "k-coordinate": clean_tuple_str(tuple(k_vec))
            })

    df_out = pd.DataFrame(out_rows, columns=["i", "j", "j-coordinate", "k_j", "k-coordinate"])
    df_out.to_csv(output_file, sep=";", index=False)
    print(f"det_K_suit wrote {len(df_out)} rows → {output_file}")
    return output_file

# --- main step 2 ---
def det_K_match_json(det_k_csv: str, json_out: str) -> str:
    """
    Convert CSV from det_K_suit to JSON mapping:
    {
      "(i, j, j_dx, j_dy, j_dz)": [
        [k_j, k_dx, k_dy, k_dz], ...
      ]
    }
    """
    df = pd.read_csv(det_k_csv, sep=";")

    def parse_tuple(s: str) -> List[float]:
        return [float(x) for x in ast.literal_eval(s)]

    mapping = {}
    for _, row in df.iterrows():
        i_val = int(row["i"])
        j_val = int(row["j"])
        j_coord = parse_tuple(row["j-coordinate"])
        k_j = int(row["k_j"])
        k_coord = parse_tuple(row["k-coordinate"])

        key = f"({i_val}, {j_val}, {j_coord[0]}, {j_coord[1]}, {j_coord[2]})"
        mapping.setdefault(key, []).append([k_j, k_coord[0], k_coord[1], k_coord[2]])

    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(mapping, f, indent=2)

    print(f"det_K_match_json wrote {len(mapping)} keys → {json_out}")
    return json_out
