import numpy as np

def solve_dyson_equation(X_up, X_down, U):
    """
    Solves the Dyson-like matrix equations for the longitudinal K-functions:
    K↑↑ = X↑ + X↑U↑↑K↑↑ + X↑U↑↓K↓↑
    K↓↓ = X↓ + X↓U↓↓K↓↓ + X↓U↓↑K↑↓
    K↑↓ = X↑U↑↓K↓↓ + X↑U↑↑K↑↓
    K↓↑ = X↓U↓↑K↑↑ + X↓U↓↓K↓↑

    Parameters:
    - X_up, X_down: susceptibility-like matrices (numpy.ndarray)
    - U: interaction kernel, a 2x2 block of 2D numpy arrays:

    Returns:
    - K: dictionary with keys 'uu', 'dd', 'ud', 'du' representing the 4 components
    """
    # Shorten names
    Xu, Xd = X_up, X_down
    Uuu, Uud = U[0][0], U[0][1]
    Udu, Udd = U[1][0], U[1][1]

    # Size check
    assert Xu.shape == Xd.shape, "X_up and X_down must have same shape"
    n = Xu.shape[0]

    # Initialize all K matrices
    Kuu = np.zeros((n, n))
    Kdd = np.zeros((n, n))
    Kud = np.zeros((n, n))
    Kdu = np.zeros((n, n))

    # Flatten everything for solving a big linear system
    I = np.eye(n * n)

    def kron_diag(A, B):
        # Kronecker product but suitable for AXB in vectorized form
        return np.kron(A, B)

    # Dyson system (linearized): AX = B
    # Vectorize all matrices with .reshape(-1, 1)
    A11 = I - kron_diag(Xu, Uuu)
    A12 = -kron_diag(Xu, Uud)
    A21 = -kron_diag(Xd, Udu)
    A22 = I - kron_diag(Xd, Udd)

    # Solve for Kuu and Kdu first
    A_top = np.hstack((A11, A12))
    A_bot = np.hstack((A21, A22))
    A_full = np.vstack((A_top, A_bot))

    B1 = Xu.reshape(-1, 1)
    B2 = Xd.reshape(-1, 1)
    B_full = np.vstack((B1, B2))

    Kvec = np.linalg.solve(A_full, B_full)
    Kuu = Kvec[:n*n].reshape(n, n)
    Kdu = Kvec[n*n:].reshape(n, n)

    # Then solve for Kud and Kdd
    A11_ = I - kron_diag(Xu, Uuu)
    A12_ = -kron_diag(Xu, Uud)
    A21_ = -kron_diag(Xd, Udu)
    A22_ = I - kron_diag(Xd, Udd)

    A_top2 = np.hstack((A12_, A11_))  # reorder to Kud, Kdd
    A_bot2 = np.hstack((A22_, A21_))
    A_full2 = np.vstack((A_top2, A_bot2))

    # Right-hand side for second system
    B1_ = np.zeros((n*n, 1))  # K↑↓ has no inhomogeneous term
    B2_ = np.zeros((n*n, 1))  # K↓↓ has no inhomogeneous term

    Kvec2 = np.linalg.solve(A_full2, np.vstack((B1_, B2_)))
    Kud = Kvec2[:n*n].reshape(n, n)
    Kdd = Kvec2[n*n:].reshape(n, n)

    return {
        'uu': Kuu,
        'dd': Kdd,
        'ud': Kud,
        'du': Kdu
    }


def compute_chi_zz(K):
    """
    Computes longitudinal spin susceptibility from K matrices:
    χ^zz = 1/4 (K↑↑ + K↓↓ - K↑↓ - K↓↑)
    """
    return 0.25 * (K['uu'] + K['dd'] - K['ud'] - K['du'])
