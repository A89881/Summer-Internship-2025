import pandas as pd
# Reads data-set and formatts the columnns correctly
# Converts the units for distance vectors to lattice intergers
import csv
from typing import Union, Tuple

def format_data(url: str):
    i_df = pd.read_csv(url, sep=r"\s+", index_col=False)
    # for col in i_df.columns[2:5]:
    #     i_df[col] = round(i_df[col] / 2.71)
    # i_df = i_df.drop(["Jij"], axis=1)
    output_path = r"Data\formatted-data.csv"
    i_df.to_csv(output_path, sep=";", index=False)
    print(f"Done: The string url is: {output_path}")
    return output_path

import csv
from typing import Tuple

def merge_format_and_xzz(
    format_file: str,
    xzz_file: str,
    output_file=r"Data\formatted_output.dat") -> None:
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
    print(f"Done: The string url is: {output_file}")
