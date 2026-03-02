import os
import time

import matplotlib.pyplot as plt
import numpy as np


# -------------------------
# Results Directory
# -------------------------
def create_results_dir():
    """
    Creates a time stamped results directory.
    """
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    results_dir = f"results/{timestamp}"
    os.makedirs(results_dir, exist_ok=True)
    return results_dir


# Create Results Directory
RESULTS_DIR = create_results_dir()


# -------------------------
# Helper Functions
# -------------------------


def plotsol_animate(u_n, i, delta_t, sleep_time=0.1, results_dir=RESULTS_DIR):
    """
    Plots the solution at each time step for animation.
    Parameters:
    - u_n: The solution at the current time step.
    - i: The current time step index.
    - delta_t: The time step size.
    """
    if results_dir is None:
        results_dir = create_results_dir()
    # Plot the solution
    plt.plot(range(len(u_n)), u_n)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(f"Solution of the Advection Equation, t = {i * delta_t}")
    plt.savefig(f"{results_dir}/Advection.png")
    time.sleep(sleep_time)
    plt.close()


def plotsol(u_n, i, delta_t):
    """
    Plots the solution at a time step.
    Parameters:
    - u_n: The solution at the current time step.
    - i: The current time step index.
    - delta_t: The time step size.
    """
    plt.plot(range(len(u_n)), u_n)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(f"Solution of the Advection Equation, t = {i * delta_t}")
    plt.show()
    plt.close()


# ----------------------------
# Heun's Method Implementation
# ----------------------------
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


def solve(
    u_n: np.ndarray,
    c: float,
    N: int,
    delta_t: float,
    animate: bool = False,
    sleep_time: float = 0.1,
) -> np.ndarray:
    """
    Solves the advection equation for N time steps.
    """
    for i in range(N):
        u_n = update_u_n(u_n, c)
        if animate:
            plotsol_animate(u_n, i, delta_t, sleep_time=sleep_time)
    return u_n


# -------------------
# Main Simulation
# ------------------

if __name__ == "__main__":
    # -------------------------
    # Simulation Parameters
    # -------------------------

    # Advection Speed
    a = 3.0

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
    u_n = np.ones(M)
    u_n[: (M // 2)] = 0

    # Solve the PDE
    u_n = solve(u_n, c, N, delta_t, animate=True, sleep_time=0.1)

    # Dump Configuration in txt file
    with open(f"{RESULTS_DIR}/configuration.txt", "w") as f:
        if a is not None:
            f.write(f"Advection Speed: {a}\n")
        if M is not None and L is not None and delta_x is not None:
            f.write(f"Space Discretisation: M = {M}, L = {L}, delta_x = {delta_x}\n")
        if delta_t is not None and tmax is not None and N is not None:
            f.write(
                f"Time Discretisation: delta_t = {delta_t}, tmax = {tmax}, N = {N}\n"
            )
        f.write(f"Courant Number: c = {c}\n")
