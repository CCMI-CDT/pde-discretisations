import time

import matplotlib.pyplot as plt
import numpy as np

delta_t = 0.5
M = 10
a = 1
L = 100
delta_x = L / M
N = int(L / delta_t)
c = (a * delta_x) / delta_t

# Define Intermediate Variable
u_hat = np.zeros(M)

# Initial Condition
u_n = np.ones(M)
u_n[: (M // 2)] = 0


def plotsol(u_n, i):
    plt.plot(range(len(u_n)), u_n)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(f"Solution of the Advection Equation, t = {i}")
    plt.savefig("Advection.png")
    time.sleep(0.1)
    plt.close()


for i in range(N):
    # Compute u_hat
    u_hat = -(c / 2) * u_n + (1 + c / 2) * np.roll(u_n, 1)

    # Update U
    u_hat = (c / 4) * (np.roll(u_hat, 1) - u_hat)
    u_n = -(c / 4) * u_n + (1 + c / 4) * np.roll(u_n, 1) + u_hat
    # Plot Solution
    plotsol(u_n, i)

# plotsol(u_n,N)
