import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-L",
        "--Length",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "-n",
        "--time_steps",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-c",
        "--Courant",
        default=1.1,
    )
    parser.add_argument(
        "-a",
        "--speed",
        default=1.0,
    )

    args = parser.parse_args()
    L = args.Length
    n = args.time_steps
    dx = L / n
    x = np.arange(0.0, L, dx)

    un_0 = 0.0 * x
    un_0[np.where(x < 0.5)] = 1

    a = args.speed
    c = args.Courant
    dt = dx * (c / a)

    un_1 = un_0 - c * (un_0 - np.roll(un_0, 1))

    t = 0
    tmax = 2.0

    # --- set up figure ---
    fig, ax = plt.subplots()
    (line,) = ax.plot(x, un_0)
    ax.set_xlim(0, L)
    ax.set_ylim(-0.2, 1.2)
    ax.set_xlabel("x")
    ax.set_ylabel("u")

    writer = PillowWriter(fps=30)

    norms = [
        np.linalg.norm(un_0),
        np.linalg.norm(un_1),
    ]  # for Von Neumann Stability condition check

    with writer.saving(fig, "./group_4/visual/leapfrog.gif", dpi=100):
        # write initial frame
        ax.set_title(f"Leapfrog at t={t:.2f}")
        writer.grab_frame()

        while t < tmax - dt / 2:
            t += dt
            un_2 = un_0 - c * (np.roll(un_1, -1) - np.roll(un_1, 1))
            un_0 = un_1
            un_1 = un_2
            norms.append(
                np.linalg.norm(un_2)
            )  # for Von Neumann Stability condition check

            # update plot
            line.set_ydata(un_1)
            ax.set_title(f"Leapfrog at t={t:.2f}")
            writer.grab_frame()

    plt.close(fig)
    plt.plot(norms)
    plt.title(r"$c=$" + str(c))
    plt.show()
