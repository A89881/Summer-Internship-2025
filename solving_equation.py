import numpy as np
import pandas as pd
import ast
from typing import Tuple, List, Dict

def parse_j_k_mapping(file_path: str) -> Dict[Tuple[float, float, float], List[Tuple[float, float, float]]]:
    mapping = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.strip():  # Skip empty lines
                continue
            try:
                j_str, k_str = line.strip().split(';')
                j_site = ast.literal_eval(j_str.strip())
                k_sites = ast.literal_eval(k_str.strip())
                if isinstance(j_site, tuple) and isinstance(k_sites, list):
                    mapping[j_site] = k_sites
            except Exception as e:
                print(f"Error parsing line: {line.strip()}\n{e}")
    return mapping

def load_bare_response(file_path: str) -> pd.DataFrame:
    return pd.read_csv(file_path)

def collect_bare_responses(k_sites: List[Tuple[float, float, float]], bare_df: pd.DataFrame) -> Dict[Tuple, Tuple[float, float]]:
    response_map = {}
    for k in k_sites:
        match = bare_df[(bare_df['dx'] == k[0]) & (bare_df['dy'] == k[1]) & (bare_df['dz'] == k[2])]
        if not match.empty:
            row = match.iloc[0]
            response_map[k] = (row['X_up'], row['X_down'])
    return response_map

def solve_dyson_recursive(X_up_map: Dict[Tuple, float], X_dn_map: Dict[Tuple, float], 
                          U: np.ndarray, tol=1e-8, max_iter=1000) -> Tuple[float, float, float, float]:
    K_upup = K_dndn = K_updn = K_dnup = 0.0

    for _ in range(max_iter):
        K_upup_new = sum(X_up_map[k] + X_up_map[k] * (U[0, 0] * K_upup + U[0, 1] * K_dnup) for k in X_up_map)
        K_dndn_new = sum(X_dn_map[k] + X_dn_map[k] * (U[1, 1] * K_dndn + U[1, 0] * K_updn) for k in X_dn_map)
        K_updn_new = sum(X_up_map[k] * (U[0, 1] * K_dndn + U[0, 0] * K_updn) for k in X_up_map)
        K_dnup_new = sum(X_dn_map[k] * (U[1, 0] * K_upup + U[1, 1] * K_dnup) for k in X_dn_map)

        if (abs(K_upup - K_upup_new) < tol and abs(K_dndn - K_dndn_new) < tol and
            abs(K_updn - K_updn_new) < tol and abs(K_dnup - K_dnup_new) < tol):
            break

        K_upup, K_dndn, K_updn, K_dnup = K_upup_new, K_dndn_new, K_updn_new, K_dnup_new

    return K_upup, K_dndn, K_updn, K_dnup

def compute_chi_zz(K_upup: float, K_dndn: float, K_updn: float, K_dnup: float) -> float:
    return 0.25 * (K_upup + K_dndn - K_updn - K_dnup)

def compute_all_chi_zz(j_k_file: str, bare_response_file: str, U_list: List[float], output_csv: str = "chi_zz_results.csv") -> pd.DataFrame:
    j_k_map = parse_j_k_mapping(j_k_file)
    bare_df = load_bare_response(bare_response_file)
    U = np.array([[U_list[0], U_list[1]], [U_list[2], U_list[3]]])

    results = []

    for j_site, k_sites in j_k_map.items():
        responses = collect_bare_responses(k_sites, bare_df)
        X_up_map = {k: val[0] for k, val in responses.items()}
        X_dn_map = {k: val[1] for k, val in responses.items()}

        if not X_up_map or not X_dn_map:
            print(f"Warning: No valid k-sites found for j-site {j_site}")
            continue

        K_upup, K_dndn, K_updn, K_dnup = solve_dyson_recursive(X_up_map, X_dn_map, U)
        chi_zz = compute_chi_zz(K_upup, K_dndn, K_updn, K_dnup)

        results.append({
            "j_dx": j_site[0],
            "j_dy": j_site[1],
            "j_dz": j_site[2],
            "chi_zz": chi_zz
        })

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_csv, index=False)
    print(f"Saved χ^zz results to '{output_csv}'")
    return df_out