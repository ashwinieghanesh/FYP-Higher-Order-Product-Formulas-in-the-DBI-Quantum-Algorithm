from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
from scipy.linalg import expm


def diag_pinch(H: np.ndarray) -> np.ndarray:
    """
    Return the diagonal matrix corresponding to the diagonal of H.

    Parameters:
        H (np.ndarray): Input square matrix.

    Returns:
        np.ndarray: Diagonal matrix with the same diagonal as H.
    """
    return np.diag(H.diagonal())


def sigma(H: np.ndarray) -> np.ndarray:
    """
    Return the off-diagonal part of a Hamiltonian.

    Parameters:
        H (np.ndarray): Input Hamiltonian.

    Returns:
        np.ndarray: Off-diagonal part of H.
    """
    return H - diag_pinch(H)


def norm_sigma(A: np.ndarray) -> float:
    """
    Return the Frobenius norm of the off-diagonal part of a matrix.

    Parameters:
        A (np.ndarray): Input matrix.

    Returns:
        float: Norm of the off-diagonal part of A.
    """
    return np.linalg.norm(sigma(A))


def get_basic_bhmm_d_operator(N: int, power: int = 1) -> np.ndarray:
    """
    Return a diagonal BHMM operator.

    Parameters:
        N (int): Size of the matrix.
        power (int): Exponent applied to the linear diagonal values.

    Returns:
        np.ndarray: Diagonal BHMM matrix.
    """
    return np.diag([x ** power for x in np.linspace(0, 1, N)])


def exact_commutator_exponential(A: np.ndarray, B: np.ndarray, x: float = 1) -> np.ndarray:
    """
    Compute the exact commutator exponential e^(x[A,B]).

    Parameters:
        A (np.ndarray): Matrix A.
        B (np.ndarray): Matrix B.
        x (float): Scaling factor for the exponent.

    Returns:
        np.ndarray: Exponential of the commutator.
    """
    comm = A @ B - B @ A
    return expm(x * comm)


def exact_sum_hamiltonian_exponential(A: np.ndarray, B: np.ndarray, x: float = 1) -> np.ndarray:
    """
    Compute the exponential of the sum of two Hamiltonians.

    Parameters:
        A (np.ndarray): Matrix A.
        B (np.ndarray): Matrix B.
        x (float): Scaling factor for the exponent.

    Returns:
        np.ndarray: Exponential of the sum A + B.
    """
    return expm(1j * (A + B) * x)


def gc(A: np.ndarray, B: np.ndarray, x: float = 1, order: int = 2) -> np.ndarray:
    """
    Compute the product formula (GC) for given order.

    Parameters:
        A (np.ndarray): Matrix A.
        B (np.ndarray): Matrix B.
        x (float): Scaling factor for the exponent.
        order (int): Order of the GC product formula (0,1,2,3).

    Returns:
        np.ndarray: Resulting unitary from the GC formula.
    """
    if order == 0:
        return exact_commutator_exponential(A, B, x)
    elif order == 1:
        return expm(1j * A * x) @ expm(1j * B * x)
    elif order == 2:
        V_A = expm(1j * A * x)
        V_B = expm(1j * B * x)
        V_A_dag = expm(-1j * A * x)
        V_B_dag = expm(-1j * B * x)
        return V_A @ V_B @ V_A_dag @ V_B_dag
    elif order == 3:
        a = (np.sqrt(5) - 1) / 2
        b = (np.sqrt(5) + 1) / 2
        c = (3 - np.sqrt(5)) / 2
        V_A1 = expm(1j * A * x * a)
        V_B1 = expm(1j * B * x * a)
        V_A2 = expm(-1j * A * x)
        V_B2 = expm(-1j * B * x * b)
        V_A3 = expm(1j * A * x * c)
        V_B3 = expm(1j * B * x)
        return V_A1 @ V_B1 @ V_A2 @ V_B2 @ V_A3 @ V_B3
    else:
        raise ValueError("Invalid order. Must be 0,1,2, or 3.")


def db_r_on_matrix(H: np.ndarray, x: float, order: int = 2) -> np.ndarray:
    """
    Apply a single DBR rotation to a Hamiltonian.

    Parameters:
        H (np.ndarray): Input Hamiltonian.
        x (float): Rotation duration.
        order (int): Order of the GC product formula.

    Returns:
        np.ndarray: Updated Hamiltonian after rotation.
    """
    D = diag_pinch(H)
    V = gc(H, D, x, order=order)
    return V.T.conj() @ H @ V


def run_dbr_on_grid(H_0: np.ndarray, s_min: float = 0, s_max: float = 0.5,
                    nmb_grid_points: int = 50, order: int = 2):
    """
    Evaluate the norm of sigma(H) over a grid of s values.

    Returns:
        tuple: (s_list, norms_sigma_H_s)
    """
    s_list = np.linspace(s_min, s_max, nmb_grid_points)
    norms_sigma_H_s = []
    for s in s_list:
        H_s = db_r_on_matrix(H_0, x=s, order=order)
        norms_sigma_H_s.append(norm_sigma(H_s))
    return s_list, norms_sigma_H_s


def multiple_rotations(H: np.ndarray, s_min: float = 0, s_max: float = 0.5,
                       nmb_grid_points: int = 50, rotation_step: int = 1,
                       order: int = 2):
    """
    Compute norms after repeated GC rotations for different orders.

    Returns:
        tuple: (s_list, norms_sigma_H_s)
    """
    s_list = np.linspace(s_min, s_max, nmb_grid_points)
    norms_sigma_H_s = []
    D = diag_pinch(H)

    for s in s_list:
        V = None
        if order in [2, 3]:
            step_size = s / rotation_step
            V = np.eye(H.shape[0], dtype=complex)
            for _ in range(rotation_step):
                V_step = gc(D, H, np.sqrt(step_size), order=order)
                V = V @ V_step
        elif order == 0:
            if rotation_step == 1:
                V = gc(H, D, s, order=order)
        else:
            raise ValueError("Unsupported order or rotation_step.")

        if V is not None:
            H_rot = V.conj().T @ H @ V
            norms_sigma_H_s.append(norm_sigma(H_rot))

    return s_list, norms_sigma_H_s


def higher_order_pf(A: np.ndarray, s_min: float, s_max: float, points: int,
                    order: int = 2):
    """
    Evaluate GC performance against exact commutator exponential.

    Returns:
        tuple: (s_list, norm_perf)
    """
    B = diag_pinch(A)
    norm_perf = []
    s_list = np.linspace(s_min, s_max, points)

    for s in s_list:
        V_pf = gc(B, A, s, order=order)
        if order == 1:
            V_exact = exact_sum_hamiltonian_exponential(B, A, x=s)
            norm_perf.append(np.linalg.norm(V_pf - V_exact))
        else:
            V_exact = exact_commutator_exponential(A, B, x=s**2)
            norm_perf.append(np.linalg.norm(V_pf - V_exact))

    return s_list, norm_perf


def perform_k_iterations(H, s_min, s_max, grid_points, rep_no, k, order,
                         s_tot=None, s_min_accumulated=None, norm_sigma_min=None):
    """
    Recursively perform k DBR iterations.

    Returns:
        tuple: Updated H, total s, accumulated s_min, and norms.
    """
    if k == 0:
        return H, s_tot, s_min_accumulated, norm_sigma_min

    if s_tot is None:
        s_tot = 0
    if s_min_accumulated is None:
        s_min_accumulated = [0]
    if norm_sigma_min is None:
        norm_sigma_min = [norm_sigma(H)]

    s_list, all_norms_sigma_dbr = multiple_rotations(H, s_min, s_max,
                                                      grid_points, rep_no, order)
    s_opt = s_list[np.argmin(all_norms_sigma_dbr)]
    s_tot += s_opt
    s_min_accumulated.append(s_tot)
    norm_sigma_min.append(np.min(all_norms_sigma_dbr))
    H = db_r_on_matrix(H, s_opt, order)

    return perform_k_iterations(H, s_min, s_max, grid_points, rep_no,
                                k - 1, order, s_tot, s_min_accumulated,
                                norm_sigma_min)




