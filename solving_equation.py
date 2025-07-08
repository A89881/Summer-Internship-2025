import numpy as np
import pandas as pd
import ast
from typing import Dict, List, Tuple

def parse_xij_file(xij_path: str) -> Tuple[Dict[Tuple[float, float, float], float], 
                                         Dict[Tuple[float, float, float], float]]:
    """Parse Xij file with χ⁰↑ and χ⁰↓ values for relative coordinates."""
    df = pd.read_csv(xij_path, sep=';')
    x_map_up = {}
    x_map_down = {}
    for _, row in df.iterrows():
        coord = (float(row['dx']), float(row['dy']), float(row['dz']))
        x_map_up[coord] = row['χ⁰↑']
        x_map_down[coord] = row['χ⁰↓']
    return x_map_up, x_map_down

def parse_k_contrib_file(kfile_path: str) -> Dict[Tuple[float, float, float], 
                                                List[Tuple[float, float, float]]]:
    """Parse j-to-k mapping file."""
    contrib_map = {}
    with open(kfile_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            j_str, k_str_list = line.strip().split(';')
            j_coord = ast.literal_eval(j_str)
            k_coords = ast.literal_eval(k_str_list)
            contrib_map[j_coord] = k_coords
    return contrib_map

def get_x_block(rel: Tuple[float, float, float], 
               x_map_up: Dict, 
               x_map_down: Dict) -> Tuple[float, float]:
    """Get X matrix block for relative coordinate with dipole fallback."""
    x_up = x_map_up.get(rel, 0.0)
    x_down = x_map_down.get(rel, 0.0)
    if x_up == 0.0 or x_down == 0.0:  # Only apply dipole if exactly zero
        r = max(np.linalg.norm(rel), 1e-6)  # type: ignore # Avoid division by zero
        dipole = 1.0 / (r**3)
        x_up = x_up if x_up != 0.0 else dipole
        x_down = x_down if x_down != 0.0 else dipole
    return x_up, x_down

def construct_A_int(x_map_up: Dict, 
                   x_map_down: Dict, 
                   k_coords: List[Tuple[float, float, float]], 
                   U: np.ndarray) -> np.ndarray:
    """Build coupling matrix A_int (2Nx2N) from X_kqkq' blocks."""
    N = len(k_coords)
    A = np.zeros((2*N, 2*N))
    for q, kq in enumerate(k_coords):
        for qp, kqp in enumerate(k_coords):
            rel = tuple(np.subtract(kq, kqp))  # Vector kq - kqp
            x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
            # Fill 2x2 block: X_kqkqp * U
            A[2*q:2*q+2, 2*qp:2*qp+2] = np.diag([x_up, x_down]) @ U
    return A

def construct_X_int(x_map_up: Dict, 
                   x_map_down: Dict, 
                   k_coords: List[Tuple[float, float, float]]) -> np.ndarray:
    """Build X_int vector (2Nx2) from X_kqj blocks."""
    X = np.zeros((2*len(k_coords), 2))
    origin = (0.0, 0.0, 0.0)
    for q, kq in enumerate(k_coords):
        rel = tuple(np.subtract(kq, origin))  # Vector kq - i (i=origin)
        x_up, x_down = get_x_block(rel, x_map_up, x_map_down)
        X[2*q:2*q+2, :] = np.diag([x_up, x_down])
    return X

def compute_chi_zz(x_map_up: Dict, 
                  x_map_down: Dict, 
                  j_coord: Tuple[float, float, float], 
                  k_coords: List[Tuple[float, float, float]], 
                  K_int: np.ndarray, 
                  U: np.ndarray) -> float:
    """Compute Xzz = 1/4(K↑↑ + K↓↓ - K↑↓ - K↓↑)."""
    # Get Xij (2x2 diagonal)
    xij_up, xij_down = get_x_block(j_coord, x_map_up, x_map_down)
    
    # Sum over intermediates: Σ Xik_q * U * Kk_qj
    K_sum = np.zeros((2,2))
    origin = (0.0, 0.0, 0.0)
    for q, kq in enumerate(k_coords):
        rel = tuple(np.subtract(kq, origin))  # Vector kq - i
        xik_up, xik_down = get_x_block(rel, x_map_up, x_map_down)
        K_block = K_int[2*q:2*q+2, :]  # 2x2 block of Kk_qj
        K_sum += np.diag([xik_up, xik_down]) @ U @ K_block
    
    # Full Kij = Xij + Σ terms
    Kij = np.diag([xij_up, xij_down]) + K_sum
    
    # χᶻᶻ = ¼(K↑↑ + K↓↓ - K↑↓ - K↓↑)
    return (Kij[0,0] + Kij[1,1] - Kij[0,1] - Kij[1,0]) / 4.0

def compute_Xzz_all(xij_file: str, 
                   kfile: str, 
                   U_params: Tuple[float, float, float, float], 
                   output_file: str = r'Data\Xzz_output.csv', 
                   temp_txt: str = r'Data\debugg.txt'):
    """Main computation workflow."""
    # Initialize interaction matrix (2x2)
    U = np.array([
        [U_params[0], U_params[2]],  # U↑↑, U↑↓
        [U_params[3], U_params[1]]   # U↓↑, U↓↓
    ])
    
    # Load data
    x_map_up, x_map_down = parse_xij_file(xij_file)
    contrib_map = parse_k_contrib_file(kfile)
    
    # Process each j-site
    results = []
    with open(temp_txt, 'w') as debug_out:
        for j_idx, (j_coord, k_coords) in enumerate(contrib_map.items()):
            if not k_coords:
                continue
            
            try:
                # Construct and solve block matrix equation
                A_int = construct_A_int(x_map_up, x_map_down, k_coords, U)
                X_int = construct_X_int(x_map_up, x_map_down, k_coords)
                K_int = np.linalg.solve(np.eye(A_int.shape[0]) - A_int, X_int)
                
                # Compute χᶻᶻ
                chi_zz = compute_chi_zz(x_map_up, x_map_down, j_coord, k_coords, K_int, U)
                results.append({'j-coordinate': str(j_coord), 'Xzz': chi_zz})
                
                # Save first 3 intermediates for debugging
                if j_idx < 3:
                    debug_out.write(f"===== J-site #{j_idx + 1}: j = {j_coord} =====\n\n")
                    debug_out.write("A_int Matrix (2Nx2N):\n")
                    np.savetxt(debug_out, A_int, fmt="%.6e", delimiter='\t')
                    debug_out.write("\n\nX_int Matrix (2Nx2):\n")
                    np.savetxt(debug_out, X_int, fmt="%.6e", delimiter='\t')
                    debug_out.write("\n\nK_int Solution (2Nx2):\n")
                    np.savetxt(debug_out, K_int, fmt="%.6e", delimiter='\t')
                    debug_out.write("\n\n")
                    
            except np.linalg.LinAlgError:
                print(f"[WARN] Singular matrix for j = {j_coord}. Skipping.")
                continue
    
    # Save final results
    pd.DataFrame(results).to_csv(output_file, index=False)
    print(f"[SUCCESS] Results saved to {output_file}")
    print(f"[DEBUG] First 3 intermediates saved to {temp_txt}")