import pandas as pd
import csv
from typing import Dict, Tuple

def format_data(url: str) -> str:
    """Read dataset and format columns correctly with optimized pandas operations."""
    output_path = r"data\formatted-data.csv"
    # Read CSV with optimized parameters
    pd.read_csv(url, sep=r"\s+", index_col=False).to_csv(output_path, sep=";", index=False)
    print(f"Done: The string url is: {output_path}")
    return output_path

def merge_format_and_xzz(
    format_file: str,
    xzz_file: str,
    output_file: str = r"data\formatted_output.dat"
) -> None:
    """
    Optimized merge of Xzz values into the format file with faster parsing and writing.
    
    Parameters:
        format_file: path to the original format file (.csv or .txt)
        xzz_file: path to the Xzz file (CSV with j-coordinate, Xzz)
        output_file: path to the output .dat file
    """
    # Optimized Xzz file parsing using dictionary comprehension
    with open(xzz_file, "r", encoding="utf-8") as xfile:
        reader = csv.reader(xfile)
        next(reader)  # Skip header
        xzz_map = {
            tuple(int(float(x)) for x in eval(coord_str)): float(xzz_val)
            for coord_str, xzz_val in reader
        }

    # Precompute header string
    header = (
        f"{'i':<4}{'j':<4}{'dx':>5}{'dy':>5}{'dz':>5}{'Jij':>11}"
        f"{'χ⁰↑':>18}{'χ⁰↓':>18}{'Xzz':>15}\n"
    )
    
    # Process format file with buffered writing
    with (
        open(format_file, "r", encoding="utf-8") as ffile,
        open(output_file, "w", encoding="utf-8", buffering=8192) as outfile
    ):
        reader = csv.reader(ffile, delimiter=';')
        next(reader)  # Skip original header
        outfile.write(header)
        
        # Process rows in bulk
        for row in reader:
            i, j, dx, dy, dz, Jij, chi_up, chi_down = row
            coord = (int(dx), int(dy), int(dz))
            xzz_val = xzz_map.get(coord, 0.0)
            
            outfile.write(
                f"{int(i):<4}{int(j):<4}{coord[0]:5}{coord[1]:5}{coord[2]:5}"
                f"{float(Jij):11.6f}{float(chi_up):18.6E}{float(chi_down):18.6E}"
                f"{xzz_val:15.6E}\n"
            )
    
    print(f"Done: The string url is: {output_file}")