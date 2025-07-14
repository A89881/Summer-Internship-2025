import pandas as pd # Datafram imports
import ast  # For safely converting string tuples to actual tuples
import json

# Determine the range in which i and j are to each other 
# (most likely symmetrical and known rangde)
# Assume if i is origin also --> rij = (dx-0, dy-0, dz-0)
# Determine possible K that are smaller than R^2

# Determine Potential K's smaller than 5 in size
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

# === Step 2: Determine suitable K-sites near each j ===
# def det_K_suit(f_url: str,k_url: str,R: int,output_file: str,
#     shift_map: dict = None) -> str: # type: ignore
#     shift_map = shift_map or {}

#     """
#     shift_map: Optional dictionary like {(1,2): (-0.5, -0.5, -0.5), (2,1): (0.5, 0.5, 0.5)}
#                if not provided, defaults to 0-vector for all pairs.
#     """
#     df_j = pd.read_csv(f_url, sep=";", index_col=False)
#     df_k = pd.read_csv(k_url, sep=";", index_col=False)

#     j_tuples = list(zip(df_j["i"], df_j["j"], df_j["dx"], df_j["dy"], df_j["dz"]))
#     k_tuples = list(zip(df_k["k_dx"], df_k["k_dy"], df_k["k_dz"]))

#     j_coords = []
#     k_coords = []

#     for (i_site, j_site, j_dx, j_dy, j_dz) in j_tuples:
#         shift = shift_map.get((i_site, j_site), (0.0, 0.0, 0.0))

#         for (k_dx, k_dy, k_dz) in k_tuples:
#             kx_s, ky_s, kz_s = k_dx + shift[0], k_dy + shift[1], k_dz + shift[2]
#             dist2 = (j_dx - kx_s)**2 + (j_dy - ky_s)**2 + (j_dz - kz_s)**2
#             if 0 < dist2 <= R**2:
#                 j_coords.append((j_dx, j_dy, j_dz))
#                 k_coords.append((k_dx, k_dy, k_dz))

#     df = pd.DataFrame({
#         "j-coordinate": [str(j) for j in j_coords],
#         "k-coordinate": [str(k) for k in k_coords]
#     })
#     df.to_csv(output_file, sep=";", index=False)
#     print(f"Done: The string url is: {output_file} (Has been rewritten)")
#     return output_file

def det_K_suit(f_url: str, k_url: str, R: int, output_file: str, shift_map: dict = None, default_shift_val: float = 0.0):
    """
    Determine suitable k-sites for each j-site based on optional sublattice shift.
    shift_map: {(i_site, j_site): (dx, dy, dz)}
    default_shift_val: if a rule exists (e.g. Cr AFM), apply ±shift
    """
    if shift_map is None:
        shift_map = {}

    i_df_j = pd.read_csv(f_url, sep=";", index_col=False)
    i_df_k = pd.read_csv(k_url, sep=";", index_col=False)

    j_list = list(zip(i_df_j["i"], i_df_j["j"],
                      i_df_j["dx"], i_df_j["dy"], i_df_j["dz"]))
    k_list = list(zip(i_df_k["k_dx"], i_df_k["k_dy"], i_df_k["k_dz"]))

    suit_k, suit_j = [], []

    for (i_site, j_site, j_dx, j_dy, j_dz) in j_list:
        shift = shift_map.get((i_site, j_site), (0.0, 0.0, 0.0))

        for (k_dx, k_dy, k_dz) in k_list:
            # Apply site-dependent shift
            kx_shifted = k_dx + shift[0]
            ky_shifted = k_dy + shift[1]
            kz_shifted = k_dz + shift[2]

            val = (j_dx - kx_shifted)**2 + (j_dy - ky_shifted)**2 + (j_dz - kz_shifted)**2
            if 0 < val <= R**2:
                suit_k.append((k_dx, k_dy, k_dz))
                suit_j.append((j_dx, j_dy, j_dz))

    df = pd.DataFrame(list(zip(suit_j, suit_k)), columns=["j-coordinate", "k-coordinate"])
    df.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file} (Shift-aware output)")
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
def det_K_match(f_url, k_url, output_file):
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
    print(f"Done: The string url is: {output_path} (Debugg)") # type: ignore
    # return grouped
    return output_path
