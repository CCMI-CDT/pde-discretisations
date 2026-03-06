import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

def phi(k, dx):
    return k * dx

def A(k, c, dx):
    return -1j * c * np.sin(phi(k, dx)) + np.emath.sqrt(1 - c**2 *(np.sin(phi(k, dx)))**2)

def norm_A(k, c, dx):
    return np.abs(A(k, c, dx))

def arg_A(k, c, dx):
    A_num = A(k, c, dx)
    A_exact = np.exp(-1j * c * phi(k, dx))
    return np.angle(A_num / A_exact)

def diffusion_error(c, dx, k):
    return norm_A(k, c, dx)

def dispersion_error(c, dx, k):
    return arg_A(k, c, dx)

dx = 0.1
k = np.linspace(-np.pi / dx, np.pi / dx, 100)

all_diffusion_errors = {}
all_dispersion_errors = {}

c_values = [0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

for c in c_values:
    all_diffusion_errors[c] = diffusion_error(c, dx, k)
    all_dispersion_errors[c] = dispersion_error(c, dx, k)

all_diffusion_errors_df = pd.DataFrame(all_diffusion_errors, index=phi(k, dx))
all_dispersion_errors_df = pd.DataFrame(all_dispersion_errors, index=phi(k, dx))

x = all_diffusion_errors_df.index.values
cs = all_diffusion_errors_df.columns.values

norm = mpl.colors.Normalize(vmin=cs.min(), vmax=cs.max())
cmap = plt.cm.viridis


fig = plt.figure(figsize=(20, 5))
gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], height_ratios=[1], figure=fig, wspace=0.4)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharex=ax1)
cax = fig.add_subplot(gs[0, 2])  # colorbar axis spanning both rows

# --- plotting loops stay the same ---
for c in cs:
    c_norm = norm(c)
    alpha = 1.0 - 0.6 * c_norm

    ax1.plot(x, all_diffusion_errors_df[c], color=cmap(c_norm), alpha=alpha, lw=2)
    ax2.plot(x, all_dispersion_errors_df[c], color=cmap(c_norm), alpha=alpha, lw=2)

ax1.set_ylabel("Diffusion error")
ax2.set_ylabel("Dispersion error")
ax2.set_xlabel(r"$\phi$")

# Properly positioned colorbar
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
fig.colorbar(sm, cax=cax, label="Courant number $c$")
plt.savefig("./group_4/visual/leapfrog_error.png", dpi=300)
