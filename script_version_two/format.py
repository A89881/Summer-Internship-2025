import pandas as pd
# Reads data-set and formatts the columnns correctly
# Converts the units for distance vectors to lattice intergers
import csv
from typing import Union, Tuple
import chardet

def detect_encoding(file_path: str, num_bytes: int = 1024) -> str:
    with open(file_path, 'rb') as f:
        raw = f.read(num_bytes)
    return chardet.detect(raw)['encoding'] # type: ignore

def format_data(url: str, output_file: str):
    # Try automatic encoding fallback
    encoding = detect_encoding(url)
    i_df = pd.read_csv(url, sep=r"\s+", index_col=False, encoding=encoding)
    
    # Fill NaNs and fix columns
    if 'Jij' in i_df.columns:
        i_df['Jij'] = i_df['Jij'].fillna(0.0)
    else:
        i_df['Jij'] = 0.0

    if 'χ⁰↑' not in i_df.columns and 'χ⁰↓' in i_df.columns:
        i_df['χ⁰↑'] = i_df['χ⁰↓'].copy()  # copy down to up

    # If χ⁰↑ exists but has NaNs (non-magnetic), copy from χ⁰↓
    i_df['χ⁰↑'] = i_df['χ⁰↑'].fillna(i_df['χ⁰↓'])

    # Make sure order and types are correct
    i_df = i_df[['i', 'j', 'dx', 'dy', 'dz', 'Jij', 'χ⁰↑', 'χ⁰↓']]
    i_df = i_df.astype({
        'i': int, 'j': int,
        'dx': int, 'dy': int, 'dz': int,
        'Jij': float, 'χ⁰↑': float, 'χ⁰↓': float
    })

    # Save
    i_df.to_csv(output_file, sep=";", index=False)
    print(f"Done: The string url is: {output_file}")
    return output_file

def merge_format_and_xzz(
    format_file: str,
    xzz_file: str,
    output_file) -> None:
    """
    Merge Xzz values into the format file and output a .dat file.
    
    Parameters:
        format_file: path to the original format file (.csv or .txt)
        xzz_file: path to the Xzz file (CSV with j-coordinate, Xzz)
        output_file: path to the output .dat file
    """
    # Parse Xzz file
    xzz_map = {}
    with open(xzz_file, "r", encoding="utf-8") as xfile:
        reader = csv.reader(xfile)
        next(reader)  # Skip header
        for row in reader:
            coord_str, xzz_val = row
            coord = tuple(int(float(x)) for x in eval(coord_str))  # Convert to (dx, dy, dz)
            xzz_map[coord] = float(xzz_val)

    # Read format file and write output
    with open(format_file, "r", encoding="utf-8") as ffile, open(output_file, "w", encoding="utf-8") as outfile:
        reader = csv.reader(ffile, delimiter=';')
        header = next(reader)  # Skip original header

        # Write new header with fixed-width spacing
        outfile.write(f"{'i':<4}{'j':<4}{'dx':>5}{'dy':>5}{'dz':>5}{'Jij':>11}"
                      f"{'χ⁰↑':>18}{'χ⁰↓':>18}{'Xzz':>15}\n")

        for row in reader:
            i, j, dx, dy, dz, Jij, chi_up, chi_down = row
            dx_i, dy_i, dz_i = int(dx), int(dy), int(dz)
            coord = (dx_i, dy_i, dz_i)
            xzz_val = xzz_map.get(coord, 0.0)

            outfile.write(
                f"{int(i):<4}{int(j):<4}{dx_i:5}{dy_i:5}{dz_i:5}"
                f"{float(Jij):11.6f}{float(chi_up):18.6E}{float(chi_down):18.6E}{xzz_val:15.6E}\n"
            )
    print(f"Done: The string url is: {output_file} (Result)")
