import numpy as np

# Define TDDFT kernels (diagonal: up/down)
X_up = 0.8
X_down = 0.9
X = np.diag([X_up, X_down])  # 2x2 matrix

# Define XC response matrix U
U = np.array([
    [2.0, 0.3],
    [0.3, 1.8]
])

# Compute K-matrix via Dyson equation: K = X (I - U X)^(-1)
I = np.eye(2)
K = X @ np.linalg.inv(I - U @ X)

# Compute longitudinal susceptibility
chi_zz = 0.25 * (K[0, 0] + K[1, 1] - K[0, 1] - K[1, 0])
print(f"Longitudinal susceptibility χ^zz = {chi_zz:.4f}")
