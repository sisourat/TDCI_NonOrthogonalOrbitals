"""
libdyn.py

Numba-optimized solver for the time-dependent Schrödinger equation:
    i * hbar * d/dt |psi(t)> = H(t) |psi(t)>

Usage:
- If H is a constant numpy array, exact evolution via diagonalization is used.
- If H is a callable H(t) -> ndarray, the Hamiltonian is sampled on the time
  grid (and midpoints) and a numba-jitted RK4 solver is used.

Requires:
  numpy, numba, scipy (optional if you want matrix exponential instead of diag),
  matplotlib for the demo plot.
"""

import numpy as np
from numpy.linalg import eigh, norm
import math
from numba import njit, types
from numba.typed import List

# ---------------------------
# Numba-compiled helpers
# ---------------------------

@njit
def _matvec_mul(H, v, N):
    """Compute H @ v with explicit loops (H: (N,N), v: (N,)) -> out (N,)"""
    out = np.zeros(N, dtype=np.complex128)
    for i in range(N):
        s = 0+0j
        row = H[i]
        for j in range(N):
            s += row[j] * v[j]
        out[i] = s
    return out

@njit
def _compute_rhs(H, psi, hbar, N):
    """Return -1j/hbar * H @ psi"""
    tmp = _matvec_mul(H, psi, N)
    coeff = -1j / hbar
    for i in range(N):
        tmp[i] = coeff * tmp[i]
    return tmp

@njit
def rk4_step_precomputed(Hn, Hhalf, Hn1, psi_n, dt, hbar, N):
    """
    Single RK4 step using precomputed H at t_n, t_n+dt/2, t_n+dt:
      Hn:   Hamiltonian at t_n (N,N)
      Hhalf: Hamiltonian at t_n + dt/2 (N,N)
      Hn1:  Hamiltonian at t_n + dt (N,N)
    psi_n: state at t_n (N,)
    returns psi_{n+1} (N,)
    """
    k1 = _compute_rhs(Hn, psi_n, hbar, N)

    psi_tmp = np.empty(N, dtype=np.complex128)
    for i in range(N):
        psi_tmp[i] = psi_n[i] + 0.5 * dt * k1[i]
    k2 = _compute_rhs(Hhalf, psi_tmp, hbar, N)

    for i in range(N):
        psi_tmp[i] = psi_n[i] + 0.5 * dt * k2[i]
    k3 = _compute_rhs(Hhalf, psi_tmp, hbar, N)

    for i in range(N):
        psi_tmp[i] = psi_n[i] + dt * k3[i]
    k4 = _compute_rhs(Hn1, psi_tmp, hbar, N)

    psi_np1 = np.empty(N, dtype=np.complex128)
    for i in range(N):
        psi_np1[i] = psi_n[i] + dt/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])

    return psi_np1

@njit
def rk4_integrate_numba(Hn_arr, Hhalf_arr, H1_arr, psi0, dt_array, normalize, hbar):
    """
    Integrate for all steps using RK4 with pre-sampled Hamiltonians.
    Hn_arr, Hhalf_arr, H1_arr: arrays of shape (M, N, N) where M = nt-1
    psi0: initial state (N,)
    dt_array: array of length M with step sizes
    Returns psi_t: array (M+1, N)
    """
    M, N, _ = Hn_arr.shape[0], psi0.size, psi0.size
    psi_t = np.zeros((M+1, N), dtype=np.complex128)
    psi_t[0] = psi0.copy()

    for step in range(M):
        Hn = Hn_arr[step]
        Hhalf = Hhalf_arr[step]
        H1 = H1_arr[step]
        dt = dt_array[step]

        psi_n = psi_t[step]
        psi_np1 = rk4_step_precomputed(Hn, Hhalf, H1, psi_n, dt, hbar, N)

        if normalize:
            # compute norm
            s = 0.0
            for i in range(N):
                val = psi_np1[i]
                s += (val.real*val.real + val.imag*val.imag)
            if s > 0:
                inv = 1.0 / math.sqrt(s)
                for i in range(N):
                    psi_np1[i] = psi_np1[i] * inv

        psi_t[step+1] = psi_np1

    return psi_t

# ---------------------------
# High-level solver
# ---------------------------

def _is_hermitian(mat, atol=1e-12):
    return np.allclose(mat, mat.conj().T, atol=atol)

def solve_tdse_numba(H, psi0, t_grid, hbar=1.0, method='auto', normalize=True):
    """
    Solve the TDSE using either exact diagonalization (time-independent H)
    or a numba-optimized RK4 (time-dependent H).

    Parameters
    ----------
    H : ndarray (N,N) or callable H(t) -> ndarray (N,N)
    psi0 : ndarray (N,) complex
    t_grid : 1D array of times (increasing)
    hbar : float
    method : 'auto' | 'expm' | 'numba_rk4'
        'expm' uses exact diagonalization (only for time-independent H).
        'numba_rk4' uses JIT RK4 and requires either time-dependent H or
        pre-sampled H.
        'auto' selects 'expm' for constant H else 'numba_rk4'.
    normalize : bool
        Whether to renormalize after each step.
    Returns
    -------
    psi_t : ndarray shape (len(t_grid), N)
    """
    t_grid = np.asarray(t_grid)
    nt = len(t_grid)
    if nt < 2:
        raise ValueError("t_grid must contain at least two times")

    psi0 = np.asarray(psi0, dtype=np.complex128).reshape((-1,))
    N = psi0.size

    # If H is a numpy array -> time-independent
    if isinstance(H, np.ndarray):
        H_const = H.astype(np.complex128)
        if method == 'auto':
            method = 'expm'
        if method == 'numba_rk4':
            # we can still run RK4 with constant H by sampling
            pass

        if method == 'expm':
            # Use exact diagonalization: U(t) = exp(-i H t / hbar)
            # Diagonalize H once
            if not _is_hermitian(H_const):
                raise ValueError("Hamiltonian should be Hermitian for physical evolution (time-independent path).")
            evals, evecs = eigh(H_const)
            # psi(t) = V * diag(exp(-i*e*t/hbar)) * V^dagger * psi0
            coeffs = evecs.conj().T.dot(psi0)
            psi_t = np.zeros((nt, N), dtype=np.complex128)
            for i, t in enumerate(t_grid):
                ph = np.exp(-1j * evals * t / hbar)
                psi_t[i] = evecs.dot(ph * coeffs)
                if normalize:
                    psi_t[i] = psi_t[i] / norm(psi_t[i])
            return psi_t

    # Fall back to numba RK4 (time-dependent or requested)
    if method == 'auto':
        method = 'numba_rk4'

    if method != 'numba_rk4':
        raise ValueError("Unknown method requested.")

    # -------------------------
    # Pre-sample Hamiltonian arrays H(t_n), H(t_n + dt/2), H(t_n+dt)
    # -------------------------
    # This allows the numba integrator to be pure numeric (no Python callables).
    dt_array = np.diff(t_grid)
    M = nt - 1
    Hn_arr = np.zeros((M, N, N), dtype=np.complex128)
    Hhalf_arr = np.zeros((M, N, N), dtype=np.complex128)
    H1_arr = np.zeros((M, N, N), dtype=np.complex128)

    # Helper to get H(t) for both array or callable
    def _get_H_at(t):
        if isinstance(H, np.ndarray):
            return H.astype(np.complex128)
        else:
            Ht = H(t)
            return np.asarray(Ht, dtype=np.complex128)

    for i in range(M):
        t_n = t_grid[i]
        dt = dt_array[i]
        Hn_arr[i] = _get_H_at(t_n)
        Hhalf_arr[i] = _get_H_at(t_n + 0.5 * dt)
        H1_arr[i] = _get_H_at(t_n + dt)

    # Call numba integrator
    psi_t = rk4_integrate_numba(Hn_arr, Hhalf_arr, H1_arr, psi0, dt_array, normalize, hbar)
    return psi_t

# ---------------------------
# Demo / Example: driven two-level system
# ---------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=np.complex128)

    # parameters
    Omega = 1.0
    Delta = 0.5
    omega_drive = 1.2

    def H_driven(t):
        # H(t) = -0.5 * Delta * sigma_z + 0.5 * Omega * cos(omega_drive * t) * sigma_x
        return -0.5 * Delta * sigma_z + 0.5 * Omega * math.cos(omega_drive * t) * sigma_x

    psi0 = np.array([1+0j, 0+0j])
    t_max = 50.0
    nt = 2001
    t = np.linspace(0.0, t_max, nt)

    print("Compiling and running numba-integrated RK4 (first run may JIT compile)...")
    psi_t = solve_tdse_numba(H_driven, psi0, t, hbar=1.0, method='numba_rk4', normalize=True)

    pop0 = np.abs(psi_t[:,0])**2
    pop1 = np.abs(psi_t[:,1])**2

    plt.figure(figsize=(8,4.5))
    plt.plot(t, pop0)
    plt.plot(t, pop1)
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title("Two-level driven system populations (Numba RK4)")
    plt.legend(["|0> population", "|1> population"])
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("Final state at t =", t[-1])
    print(psi_t[-1])
    print("Normalization:", np.vdot(psi_t[-1], psi_t[-1]).real)
