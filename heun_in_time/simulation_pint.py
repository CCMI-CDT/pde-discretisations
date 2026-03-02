from time import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from simulation import plotsol, solve, update_u_n

# Advection Speed
a = 1.0

# Space Discretisation
M = 100
L = 100
delta_x = L / M

# Time Discretisation
delta_t = 0.5
tmax = 10
N = int(tmax / delta_t)  # Number of time steps

# Courant Number
c = (a * delta_x) / delta_t

# Initial Condition
u_0 = np.ones(M)
u_0[: (M // 2)] = 0
# Full Initial Condition in time and space
x0 = np.zeros((N, M))
x0[0, :] = u_0


def matvec(x):
    """
    Creates a Parallel in time Operator.
    Identity in block diagonal, -P=-update_u_n in the first sub-diagonal.
    [I, 0, 0, ..., 0, 0]
    [-P, I, 0, ..., 0, 0]
    [0, -P, I, ..., 0, 0]
    ...
    [0, 0, 0, ..., -P, I]
    """
    # Reshape x to (N, M)
    x_reshaped = x.reshape((N, M))

    # Initialize the output vector
    y = np.empty_like(x_reshaped)

    # First block row remains unchanged
    y[0, :] = x_reshaped[0, :]

    # Vectorized computation for the rest of the rows.
    y[1:, :] = x_reshaped[1:, :] - update_u_n(x_reshaped[:-1, :], c)

    return y.flatten()


# TODO: Implement PinT preconditioner for gmres.

# Create Linear Operator
A = LinearOperator(shape=(N * M, N * M), matvec=matvec)

# Solve Ax = b using GMRES
b = x0.flatten()
tol = 1e-12
t0 = time()
x, info = gmres(A, b, rtol=tol, maxiter=5)
t1 = time()
print(f"Time taken for GMRES solution: {t1 - t0:.4f} seconds")

# Select the final solution at the last time step
u_final = x.reshape((N, M))[-1, :]

# Compute Solution without gmres for comparison
t0 = time()
u_serial = solve(u_0, c, N - 1, delta_t)
t1 = time()
print(f"Time taken for serial solution: {t1 - t0:.4f} seconds")

# Compute error between the two solutions
error = np.linalg.norm(u_final - u_serial)
print(f"Error between GMRES solution and serial solution: {error}")


# Plot the final solution
plotsol(u_final, N - 1, delta_t)
# Plot the final solution
plotsol(u_serial, N - 1, delta_t)
# Plot the error
plotsol(np.abs(u_final - u_serial), N - 1, delta_t)
