import pandas as pd
import ast  
import json 
from typing import Optional, Tuple
import numpy as np

# Determine the range in which i and j are to each other 
# (most likely symmetrical and known rangde)
# Assume if i is origin also --> rij = (dx-0, dy-0, dz-0)
# Determine possible K that are smaller than R^2

# Determine Potential K's smaller than R in size
def det_K_pot(min:int, max:int, R:int, output_file):
    k_x_vals = []
    k_y_vals = []
    k_z_vals = []
    for i in range(min, max+1):
        for j in range(min, max+1):
            for k in range(min, max+1):
                if 0 < i**2+j**2+k**2 <= R**2:
                    k_x_vals.append(i)
                    k_y_vals.append(j)
                    k_z_vals.append(k)
    df = pd.DataFrame(list(zip(k_x_vals, k_y_vals, k_z_vals)))
    df.columns = ["k_dx", "k_dy", "k_dz"]
    output_path = output_file
    df.to_csv(output_path, sep=";", index=False)
    print(f"Done: The string url is: {output_path} (Result)") # type: ignore
    return output_path

def to_tuple(x):
    return tuple(float(i) if isinstance(i, np.floating) else int(i) for i in x)

def clean_tuple_str(t: Tuple[float, float, float]) -> str:
    """Formats tuple as a string, removing .0 when the value is an integer."""
    return "(" + ", ".join(
        str(int(x)) if x == int(x) else str(x)
        for x in t
    ) + ")"

def det_K_suit(f_url: str,
               k_url: str,
               R: int,
               base_change: list,
               output_file: str) -> str:
    import numpy as np
    import pandas as pd

    df_j = pd.read_csv(f_url, sep=";", index_col=False)
    df_k = pd.read_csv(k_url, sep=";", index_col=False)

    j_coords = list(zip(df_j["dx"], df_j["dy"], df_j["dz"]))
    k_coords = list(zip(df_k["k_dx"], df_k["k_dy"], df_k["k_dz"]))
    T = np.array(base_change)

    matched_j = []
    matched_k = []

    for j in j_coords:
        j_vec = np.array(j, dtype=float)
        for k in k_coords:
            k_vec = np.array(k, dtype=float)
            rel_vec = k_vec - j_vec
            rel_transformed = rel_vec @ T
            dist2 = np.dot(rel_transformed, rel_transformed)

            if 0 < dist2 <= R**2:
                matched_j.append(to_tuple(j_vec))   # Convert to clean Python tuple
                matched_k.append(to_tuple(k_vec))   # Convert to clean Python tuple

    df_out = pd.DataFrame({
        "j-coordinate": [clean_tuple_str(j) for j in matched_j],
        "k-coordinate": [clean_tuple_str(k) for k in matched_k]
    })

    df_out.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file} (Clean tuple output)")
    return output_file


# === Step 2: Determine suitable K-sites near each j ===
def det_K_suit_proto(f_url: str,
               k_url: str,
               R: int,
               output_file: str) -> str:
    """
    Determines valid (j, k) pairs where k lies within radius R of j.
    Now assumes dx,dy,dz have already been shifted in `format_data`.

    Args:
        f_url: Path to formatted file with shifted dx,dy,dz.
        k_url: Path to file with k_dx, k_dy, k_dz.
        R: Cutoff radius for neighborhood.
        output_file: Destination path for matched j-k pairs.

    Returns:
        Path to output CSV with 'j-coordinate', 'k-coordinate' columns.
    """
    df_j = pd.read_csv(f_url, sep=";", index_col=False)
    df_k = pd.read_csv(k_url, sep=";", index_col=False)

    j_coords = list(zip(df_j["dx"], df_j["dy"], df_j["dz"]))
    k_coords = list(zip(df_k["k_dx"], df_k["k_dy"], df_k["k_dz"]))

    matched_j = []
    matched_k = []

    for j in j_coords:
        for k in k_coords:
            dist2 = sum((j_i - k_i)**2 for j_i, k_i in zip(j, k))
            if 0 < dist2 <= R**2:
                matched_j.append(j)
                matched_k.append(k)

    df_out = pd.DataFrame({
        "j-coordinate": [str(j) for j in matched_j],
        "k-coordinate": [str(k) for k in matched_k]
    })
    df_out.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file} (Result)")
    return output_file

# === Step 3: Group suitable K by J (match step) ===
def det_K_match_json(f_url: str, k_url: str, output_file: str) -> str:
    """
    Generate a JSON file mapping each j-site to its contributing k-sites.
    Input CSVs must contain columns: 'dx','dy','dz' and 'j-coordinate','k-coordinate'.
    """
    # === Load files ===
    df_formatted = pd.read_csv(f_url, sep=";", index_col=False)
    df_matches = pd.read_csv(k_url, sep=";", index_col=False)

    # Convert coordinate strings to tuples
    df_matches['j-tuple'] = df_matches['j-coordinate'].apply(lambda s: tuple(map(float, ast.literal_eval(s))))
    df_matches['k-tuple'] = df_matches['k-coordinate'].apply(lambda s: tuple(map(float, ast.literal_eval(s))))

    # Ensure j-tuple exists in the formatted data
    known_j_coords = set(tuple(row) for row in df_formatted[['dx', 'dy', 'dz']].values)

    # Filter valid matches where the k-site exists in the system
    df_valid = df_matches[df_matches['k-tuple'].isin(known_j_coords)]

    # Group k-sites by j-tuple
    grouped_dict = {}
    for _, row in df_valid.iterrows():
        j = row['j-tuple']
        k = row['k-tuple']
        grouped_dict.setdefault(j, []).append(k)

    # Convert to JSON-safe format: tuple → str keys, tuple → list values
    json_ready = {str(j): [list(k) for k in k_list] for j, k_list in grouped_dict.items()}

    # Save to JSON
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(json_ready, f, indent=2)

    print(f"Done: The string url is: {output_file} (Result)") # type: ignore
    return output_file

# Determine K's that match and hence group subsequent J-coordinate with neighbouring
def det_K_match_csv(f_url, k_url, output_file):
    i_df_j = pd.read_csv(f_url, sep=";", index_col=False)  # File with dx, dy, dz
    i_df_k = pd.read_csv(k_url, sep=";", index_col=False)  # File with j and k coordinates
    j_set = set(list(zip(*[i_df_j[col] for col in i_df_j.columns[2:5]])))
    
    # Convert k-coordinate strings to tuples (e.g., "(-5, 0, 0)" → (-5.0, 0.0, 0.0))
    k_list_k = i_df_k['k-coordinate'].apply(lambda x: tuple(map(float, ast.literal_eval(x))))
    i_df_k['k-tuple'] = k_list_k  # Store as a new column
    df_k_filtered = i_df_k[i_df_k['k-tuple'].isin(j_set)]
    
    # Group matching k-coordinates by j-coordinate
    # Convert j-coordinate strings to tuples (same as k-coordinate)
    j_tuples = i_df_k['j-coordinate'].apply(lambda x: tuple(map(float, ast.literal_eval(x))))
    df_k_filtered['j-tuple'] = j_tuples
    
    grouped = df_k_filtered.groupby('j-tuple')['k-tuple'].apply(list).reset_index()
    output_path = output_file
    grouped.to_csv(output_path, sep=";", index=False, header=["j-coordinate", "k-coordinates"])
    print(f"Done: The string url is: {output_path} (Result)") # type: ignore
    # return grouped
    return output_path
