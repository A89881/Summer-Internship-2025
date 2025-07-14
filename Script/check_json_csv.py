import pandas as pd
import json
import ast

def load_csv_k_map(csv_path: str) -> dict:
    df = pd.read_csv(csv_path, sep=";")
    csv_map = {}

    for _, row in df.iterrows():
        j_key = tuple(map(float, ast.literal_eval(row['j-coordinate'])))

        try:
            k_list_raw = ast.literal_eval(row['k-coordinates'])
            # Make sure each k is a tuple of floats
            k_list = [tuple(map(float, k)) for k in k_list_raw]
        except Exception as e:
            print(f"[WARN] Failed to parse k-coordinates for j = {row['j-coordinate']}: {e}")
            continue

        csv_map[j_key] = sorted(k_list)

    return csv_map


def load_json_k_map(json_path: str) -> dict:
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return {
        tuple(map(float, ast.literal_eval(j_str))): sorted([tuple(map(float, k)) for k in k_list])
        for j_str, k_list in data.items()
    }

def compare_k_maps(csv_map: dict, json_map: dict):
    all_keys = set(csv_map.keys()) | set(json_map.keys())
    discrepancies = []

    for j in sorted(all_keys):
        csv_k = csv_map.get(j)
        json_k = json_map.get(j)
        if csv_k != json_k:
            discrepancies.append((j, csv_k, json_k))

    if not discrepancies:
        print("CSV and JSON K-maps match perfectly.")
    else:
        print(f"Found {len(discrepancies)} mismatched J-sites:")
        for j, csv_k, json_k in discrepancies:
            print(f"  j = {j}")
            print(f"    CSV:  {csv_k}")
            print(f"    JSON: {json_k}")
            print()

def compare_numeric_data(csv_path: str, json_path: str):
    # Load CSV data
    csv_df = pd.read_csv(csv_path, sep=";")
    
    # Load JSON data by first converting it to a DataFrame
    with open(json_path, 'r', encoding='utf-8') as f:
        json_data = json.load(f)
    json_df = pd.DataFrame.from_dict(json_data, orient='index')
    json_df.reset_index(inplace=True)
    json_df.rename(columns={'index': 'j-coordinate'}, inplace=True)
    
    # Merge data on j-coordinate
    merged = pd.merge(csv_df, json_df, on="j-coordinate", suffixes=('_csv', '_json'))
    
    # Compare numeric columns (assuming Xzz is present)
    if 'Xzz_csv' in merged.columns and 'Xzz_json' in merged.columns:
        merged["diff"] = abs(merged["Xzz_csv"] - merged["Xzz_json"])
        print("Numeric comparison results:")
        print("Max diff:", merged["diff"].max())
        print("Rows with large diff (>1e-6):")
        print(merged[merged["diff"] > 1e-6])
    else:
        print("No Xzz column found for numeric comparison")

# === Run the comparison ===
csv_file = r"Bcc-Fe\neighbouring_k_to_j.csv"
json_file = r"Bcc-Fe\neighbouring_k_to_j.json"



# Compare k-coordinates mapping
print("Comparing k-coordinates mapping...")
csv_map = load_csv_k_map(csv_file)
json_map = load_json_k_map(json_file)
compare_k_maps(csv_map, json_map)

# # Compare numeric data (like Xzz values)
# print("\nComparing numeric data...")
# compare_numeric_data(csv_file, json_file)