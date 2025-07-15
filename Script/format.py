import pandas as pd
import csv
from typing import Optional, Tuple, Dict
import chardet

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
            # dx_i, dy_i, dz_i = int(dx), int(dy), int(dz)
            dx_i, dy_i, dz_i = float(dx), float(dy), float(dz)
            coord = (dx_i, dy_i, dz_i)
            xzz_val = xzz_map.get(coord, 0.0)

            outfile.write(
                f"{int(i):<4}{int(j):<4}{dx_i:5}{dy_i:5}{dz_i:5}"
                f"{float(Jij):11.6f}{float(chi_up):18.6E}{float(chi_down):18.6E}{xzz_val:15.6E}\n"
            )
    print(f"Done: The string url is: {output_file} (Result)")
