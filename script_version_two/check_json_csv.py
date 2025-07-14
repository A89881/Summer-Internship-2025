# import pandas as pd
# import json
# import ast

# def load_csv_k_map(csv_path: str) -> dict:
#     df = pd.read_csv(csv_path, sep=";")
#     csv_map = {}

#     for _, row in df.iterrows():
#         j_key = tuple(map(float, ast.literal_eval(row['j-coordinate'])))

#         try:
#             k_list_raw = ast.literal_eval(row['k-coordinates'])
#             # Make sure each k is a tuple of floats
#             k_list = [tuple(map(float, k)) for k in k_list_raw]
#         except Exception as e:
#             print(f"[WARN] Failed to parse k-coordinates for j = {row['j-coordinate']}: {e}")
#             continue

#         csv_map[j_key] = sorted(k_list)

#     return csv_map


# def load_json_k_map(json_path: str) -> dict:
#     with open(json_path, 'r', encoding='utf-8') as f:
#         data = json.load(f)
#     return {
#         tuple(map(float, ast.literal_eval(j_str))): sorted([tuple(k) for k in k_list])
#         for j_str, k_list in data.items()
#     }

# def compare_k_maps(csv_map: dict, json_map: dict):
#     all_keys = set(csv_map.keys()) | set(json_map.keys())
#     discrepancies = []

#     for j in sorted(all_keys):
#         csv_k = csv_map.get(j)
#         json_k = json_map.get(j)
#         if csv_k != json_k:
#             discrepancies.append((j, csv_k, json_k))

#     if not discrepancies:
#         print("CSV and JSON K-maps match perfectly.")
#     else:
#         print(f"Found {len(discrepancies)} mismatched J-sites:")
#         for j, csv_k, json_k in discrepancies:
#             print(f"  j = {j}")
#             print(f"    CSV:  {csv_k}")
#             print(f"    JSON: {json_k}")
#             print()

# === Run the comparison ===
csv_file = r"AFM-Cr\neighbouring_k_to_j.csv"
json_file = r"AFM-Cr\neighbouring_k_to_j.json"

# csv_map = load_csv_k_map(csv_file)
# json_map = load_json_k_map(json_file)
# compare_k_maps(csv_map, json_map)

import pandas as pd

df1 = pd.read_csv(csv_file)
df2 = pd.read_csv(json_file)

merged = pd.merge(df1, df2, on="j-coordinate", suffixes=('_csv', '_json'))
merged["diff"] = abs(merged["Xzz_csv"] - merged["Xzz_json"])

print("Max diff:", merged["diff"].max())
print("Rows with large diff:")
print(merged[merged["diff"] > 1e-6])
