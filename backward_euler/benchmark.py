import time
import numpy as np
from scipy.linalg import circulant, lu_factor, lu_solve

# Setup (same as simulate with c=1.0, a=1.0, M=1000)
M = 1000
L = 1.0
c = 1.0
a = 1.0
t_max = 1.0
dx = L / M
dt = dx * c / a
x = np.arange(0, L, dx)
Nt = int(t_max / dt)

vals = np.zeros(M)
vals[0] = 1.0
vals[1] = c / 2.0
vals[-1] = -c / 2.0
A = circulant(vals)

# --- Direct (LU every step) ---
un = np.zeros_like(x)
un[np.where(x > 0.5)] = 1.0
t = 0.0
start = time.perf_counter()
while t < t_max - dt / 2:
    t += dt
    un = np.linalg.solve(A, un)
direct_time = time.perf_counter() - start

# --- Pre-factored LU ---
un = np.zeros_like(x)
un[np.where(x > 0.5)] = 1.0
t = 0.0
start = time.perf_counter()
lu_piv = lu_factor(A)
while t < t_max - dt / 2:
    t += dt
    un = lu_solve(lu_piv, un)
lu_time = time.perf_counter() - start

# --- FFT ---
un = np.zeros_like(x)
un[np.where(x > 0.5)] = 1.0
t = 0.0
start = time.perf_counter()
eigenvalues = np.fft.fft(A[:, 0])
while t < t_max - dt / 2:
    t += dt
    un = np.real(np.fft.ifft(np.fft.fft(un) / eigenvalues))
fft_time = time.perf_counter() - start

print(f"M = {M}, Nt = {Nt} time steps\n")
print(f"{'Method':<20} {'Time (s)':>10} {'Speedup vs Direct':>20}")
print("-" * 52)
print(f"{'Direct (LU/step)':<20} {direct_time:>10.4f} {'1.0x':>20}")
print(f"{'Pre-factored LU':<20} {lu_time:>10.4f} {direct_time/lu_time:>19.1f}x")
print(f"{'FFT':<20} {fft_time:>10.4f} {direct_time/fft_time:>19.1f}x")
