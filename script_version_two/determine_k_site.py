import pandas as pd # Datafram imports
import ast  # For safely converting string tuples to actual tuples

# Determine the range in which i and j are to each other 
# (most likely symmetrical and known rangde)
# Assume if i is origin also --> rij = (dx-0, dy-0, dz-0)
# Determine possible K that are smaller than R^2

# Determine Potential K's smaller than 5 in size
def det_K_pot(min:int, max:int, R:int):
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
    output_path = r"data\k_pot_coord.csv" # type: ignore
    df.to_csv(output_path, sep=";", index=False)
    print(f"Done: The string url is: {output_path} (Result)") # type: ignore
    return output_path


# Determine K's in the proximity of R to j-coordinates
def det_K_suit(f_url:str, k_url: str, R: int):
    i_df_j = pd.read_csv(f_url, sep=";", index_col=False)
    i_df_k = pd.read_csv(k_url, sep=";", index_col=False)
    j_list = list(zip(*[i_df_j[col] for col in i_df_j.columns[2:5]]))
    k_list = list(zip(*[i_df_k[col] for col in i_df_k.columns]))
    suit_k = []
    suit_j = []
    
    for _, (j_dx, j_dy, j_dz) in enumerate(j_list):
        for _, (k_dx, k_dy, k_dz) in enumerate(k_list):
            val = j_dx**2+j_dy**2+j_dz**2+k_dx**2+k_dy**2+k_dz**2-2*(j_dx*k_dx+j_dy*k_dy+j_dz*k_dz) 
            if 0 < val <= 25:
                suit_k.append((k_dx,k_dy,k_dz))
                suit_j.append((j_dx,j_dy,j_dz))
    df = pd.DataFrame(list(zip(suit_j,suit_k)))
    df.columns = ["j-coordinate", "k-coordinate"]
    output_path = r"data\k_pot_coord.csv" # type: ignore
    df.to_csv(output_path, sep=";", index=False)
    print(f"Done: The string url is: {output_path} (Has been rewritten)") # type: ignore
    return output_path

# Determine K's that match and hence group subsequent J-coordinate with neighbouring
def det_K_match(f_url, k_url):
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
    output_path = r"data\neighbouring_k_to_j.csv"
    grouped.to_csv(output_path, sep=";", index=False, header=["j-coordinate", "k-coordinates"])
    print(f"Done: The string url is: {output_path} (Result)") # type: ignore
    # return grouped
    return output_path