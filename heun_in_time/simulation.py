import time

import matplotlib.pyplot as plt
import numpy as np


# Helper Functions
def plotsol_animate(u_n, i):
    plt.plot(range(len(u_n)), u_n)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(f"Solution of the Advection Equation, t = {i}")
    plt.savefig("Advection.png")
    time.sleep(0.1)
    plt.close()


def plotsol(u_n, i):
    plt.plot(range(len(u_n)), u_n)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(f"Solution of the Advection Equation, t = {i}")
    plt.show()
    plt.close()

def compute_u_hat(u_n: np.ndarray, c: float) -> np.ndarray:
    """
    Computes the intermediate variable u_hat.
    """
    return -(c / 2) * u_n + (1 + c / 2) * np.roll(u_n, 1)


def update_u_hat(u_hat: np.ndarray, c: float) -> np.ndarray:
    """
    Computes the update to u_hat.
    """
    return (c / 4) * (np.roll(u_hat, 1) - u_hat)


def update_u_n(u_n: np.ndarray, c: float) -> np.ndarray:
    """
    Computes the next time step.
    """
    # 1. Compute prediction
    u_hat = compute_u_hat(u_n, c)

    # 2. Update u_hat
    u_hat = update_u_hat(u_hat, c)

    # 3. Final integration step
    return -(c / 4) * u_n + (1 + c / 4) * np.roll(u_n, 1) + u_hat


def solve(u_n: np.ndarray, c: float, N: int, animate: bool = False) -> np.ndarray:
    """
    Solves the advection equation for N time steps.
    """
    for i in range(N):
        u_n = update_u_n(u_n, c)
        if animate:
            plotsol_animate(u_n, i)
    return u_n


if __name__ == "__main__":
    # Simulation Parameters
    delta_t = 0.5
    M = 10
    a = 1
    L = 100
    delta_x = L / M
    N = int(L / delta_t)
    c = (a * delta_x) / delta_t

    # Initial Condition
    u_n = np.ones(M)
    u_n[: (M // 2)] = 0

    # Solve the PDE
    u_n = solve(u_n, c, N, animate=True)

    # Plot the final solution
    plotsol(u_n, N)
