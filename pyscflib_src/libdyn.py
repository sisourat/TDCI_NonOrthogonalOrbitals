import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Define the RHS function with hmat_interp as an argument
def rhs(t, psi, hmat_interp):
     """Right-hand side of the TDSE with hmat_interp as an argument."""
     hmat = hmat_interp(t)  # Interpolate H(t) at time t
     #print(t,hmat[0,0].real)
     return -1j * np.dot(hmat, psi)

def solve_tdse(hmat_interp, psi0, t_grid):

 sol = solve_ivp(
    lambda t, psi: rhs(t, psi, hmat_interp),  # Pass hmat_interp as an argument
    (t_grid[0], t_grid[-1]),
    psi0,
    t_eval=t_grid,
    method='DOP853',
    rtol=1e-6,
    atol=1e-8,
    vectorized=False  # Ensure the function is not treated as vectorized
    )

 # Extract the solution
 return sol.y.T

