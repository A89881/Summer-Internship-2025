import pandas as pd
import csv
from typing import Optional, Tuple, Dict
import chardet

# Script to detect encoding regarding the inputed data file
def detect_encoding(file_path: str, num_bytes: int = 1024) -> str:
    with open(file_path, 'rb') as f:
        raw = f.read(num_bytes)
    return chardet.detect(raw)['encoding'] # type: ignore

def format_data(url: str,
                output_file: str,
                shift_map: Optional[Dict[Tuple[int, int], Tuple[float, float, float]]] = None,
                scale_response: float = -1.0) -> str:
    """
    Load, clean, shift, and optionally scale response data.
    Args:
        url: Path to raw input data (whitespace-separated).
        output_file: Destination CSV path.
        shift_map: Optional site-wise coordinate shift {(i,j): (dx,dy,dz)}.
        scale_response: Multiplier for χ⁰↑ and χ⁰↓ (e.g. -1 to invert signs).
    Returns:
        Path to the saved formatted CSV.
    """
    if shift_map is None:
        shift_map = {}

    encoding = detect_encoding(url)
    df = pd.read_csv(url, sep=r"\s+", index_col=False, encoding=encoding)

    # Add missing or fallback columns
    df['Jij'] = df.get('Jij', 0.0).fillna(0.0) # type: ignore
    if 'χ⁰↑' not in df.columns and 'χ⁰↓' in df.columns:
        df['χ⁰↑'] = df['χ⁰↓'].copy()
    df['χ⁰↑'] = df['χ⁰↑'].fillna(df['χ⁰↓'])

    # Enforce column order and types
    df = df[['i', 'j', 'dx', 'dy', 'dz', 'Jij', 'χ⁰↑', 'χ⁰↓']].astype({
        'i': int, 'j': int,
        'dx': float, 'dy': float, 'dz': float,
        'Jij': float, 'χ⁰↑': float, 'χ⁰↓': float
    })

    # Apply site-dependent coordinate shift
    shifted = []
    for _, row in df.iterrows():
        i_site, j_site = row['i'], row['j']
        dx, dy, dz = row['dx'], row['dy'], row['dz']
        shift = shift_map.get((i_site, j_site), (0.0, 0.0, 0.0))
        shifted.append((dx + shift[0], dy + shift[1], dz + shift[2]))
    df[['dx', 'dy', 'dz']] = pd.DataFrame(shifted, index=df.index)

    # Apply scaling to static response columns
    df['χ⁰↑'] *= scale_response
    df['χ⁰↓'] *= scale_response

    # Save and return
    df.to_csv(output_file, sep=';', index=False)
    print(f"Done: The string url is: {output_file} (Shift-augmented, Scaled={scale_response})")
    return output_file

import csv

def merge_format_and_xzz(format_file: str, xzz_file: str, output_file: str) -> None:
    """
    Merge Xzz values into the original format file and output a .dat file for analysis.

    Parameters:
        format_file: path to the original format file (.csv or .txt)
        xzz_file: path to the Xzz CSV file (must have columns: i;j;j-coordinate;N_k;Xzz)
        output_file: path to the output .dat file
    """
    # === 1. Parse Xzz file into a mapping of coordinates to Xzz values ===
    xzz_map = {}
    with open(xzz_file, "r", encoding="utf-8") as xf:
        reader = csv.DictReader(xf, delimiter=';')
        for row in reader:
            # Extract j-coordinate string and convert to tuple of floats
            coord_str = row['j-coordinate']
            coord = tuple(float(x) for x in eval(coord_str))
            xzz_val = float(row['Xzz'])
            xzz_map[coord] = xzz_val

    # === 2. Read format file and write merged .dat output ===
    with open(format_file, "r", encoding="utf-8") as ffile, \
         open(output_file, "w", encoding="utf-8") as outfile:
        
        reader = csv.reader(ffile, delimiter=';')
        header = next(reader)  # skip original header

        # Write new header with fixed-width spacing
        outfile.write(f"{'i':<4}{'j':<4}{'dx':>7}{'dy':>7}{'dz':>7}"
                      f"{'Jij':>12}{'χ⁰↑':>15}{'χ⁰↓':>15}{'Xzz':>15}\n")

        for row in reader:
            i, j, dx, dy, dz, Jij, chi_up, chi_down = row
            dx_f, dy_f, dz_f = float(dx), float(dy), float(dz)
            coord = (dx_f, dy_f, dz_f)

            xzz_val = xzz_map.get(coord, 0.0)  # default to 0.0 if not found

            outfile.write(
                f"{int(i):<4}{int(j):<4}{dx_f:7.3f}{dy_f:7.3f}{dz_f:7.3f}"
                f"{float(Jij):12.6f}{float(chi_up):15.6E}{float(chi_down):15.6E}{xzz_val:15.6E}\n"
            )

    print(f"Done: merged .dat file written to {output_file}")

