import numpy as np
import matplotlib.pyplot as pp
from scipy.linalg import circulant
import os
import argparse
'''
L = 1.0
M = 100
dx = L / M
x = np.arange(0, L, dx)
un = np.zeros_like(x)
un[np.where(x>0.5)] = 1.0
c = 1.0
a = 2.0
dt = dx*c/a
t=0.0
t_max = 1.0
tdump = 0.2
dumpt = 0.0
dcount = 0
'''

def stability_analysis(A):
    '''
    Check if the magnitude of eigenvalues of the 
    linear transformation for the backward Euler
    centered scheme is less than/equal to 1
    '''
    Ainv = np.linalg.inv(A)
    eigAinv = np.linalg.eigvals(Ainv)
    if max(np.abs(eigAinv)) > 1.0 + 1e-8:
        print("Unstable")
    else:
        print("Stable")

def build_matrix_centered(c, M):
    '''
    Build a circulant matrix for the backward Euler
    centered scheme applied to the 1-D advection equation
    '''
    vals = np.zeros(M)
    vals[0] = 1.0
    vals[1] = c/2
    vals[-1] = -c/2
    A = circulant(vals)   
    return A

def build_matrix_upwind(c, M):
    '''
    Build a circulant matrix for the backward Euler
    upwind scheme applied to the 1-D advection equation
    '''
    vals = np.zeros(M)
    vals[0] = 1.0 + c
    vals[-1] = -c
    A = circulant(vals)
    return A

def backward_euler(A, un, solver = 'fourier'):
    """
    If selected solver == 'direct'
        Solve Ax = b using direct inversion
        Computationally scales as O(M^3)
    If selected solver == 'fourier'
        Solve Ax = b using FFT
        Where A is our circulant matrix for the backward Euler 
        centered scheme for the 1D advection equation.
        computationally scales as O(M log M)
    """
    if solver == 'direct':
        un = np.linalg.solve(A, un)
        return un
    elif solver == 'fourier':
        # compute p = Fb (Fourier transform of b, a vertical stack of u_n)
        p = np.fft.fft(un) 
        # s = Fa_1
        a_1 = A[:,0] # First column of A
        s = np.fft.fft(a_1)
        # compute q = diag(s)^-1 p
        q = p / s
        # Solve for x = A^-1 b = F diag(s)^-1 F b
        x = np.fft.ifft(q).real
        return x 
    else:
        raise ValueError(f"Unknown solver: {solver}")

def dump(x, un, t, dcount, c, solver, results_dir='results'):
    # Create results directory with c subfolder if it doesn't exist
    c_dir = os.path.join(results_dir, f'c_{c}')
    os.makedirs(c_dir, exist_ok=True)
    pp.plot(x, un)
    pp.xlabel('x')
    pp.ylabel('y')
    pp.title(f'Solution at time {t}')
    filepath = os.path.join(c_dir, f'beup_solver_{solver}_c{c}_dump{dcount}.jpg')
    pp.savefig(filepath)
    pp.close()

def simulate(c, a, t_max=1.0, tdump=0.2, L=1.0, M=100, results_dir='results', method='upwind', solver = 'fourier'):
    dx = L / M
    x = np.arange(0, L, dx)
    dt = dx*c/a
    un = np.zeros_like(x)
    un = np.sin(2*np.pi*x)
    if method == 'upwind':
        A = build_matrix_upwind(c, M)
    elif method == 'centered':
        A = build_matrix_centered(c, M)
    else:
        raise ValueError
    t = 0.0
    dcount = 0
    dumpt = 0
    while t < t_max - dt/2:
        t += dt
        un = backward_euler(A, un, solver = solver)

        dumpt += dt
        if dumpt > tdump - dt/2:
            dump(x, un, t, dcount, c, solver, results_dir)
            dumpt -= tdump
            dcount += 1

def plot_diffusion_dispersion_error(c_values, M=100, method='upwind'):
    """
    Plot diffusion and dispersion errors for multiple CFL numbers
    for Fourier modes phi in [-pi, pi] using the correct phase factor.

    Diffusion error = |A / A_exact|
    Dispersion error = arg(A / A_exact)
    """
    # Fourier modes: l = -M/2 ... M/2-1
    l = np.arange(-M//2, M//2)
    phi = 2 * np.pi * l / M  # Standard 2π normalization, φ ∈ [-π, π]

    pp.figure(figsize=(12,5))

    # Diffusion error plot
    pp.subplot(1,2,1)
    for c in c_values:
        if method == 'upwind':
            # Upwind Backward Euler amplification
            A = 1 / (1 + c - c * np.exp(-1j * phi))
        elif method == 'centered':
            # Centered Backward Euler amplification
            A = 1 / (1 - 1j * c * np.sin(phi))
        else:
            raise ValueError("Unknown method")

        # Exact PDE amplification (pure advection)
        A_exact = np.exp(-1j * c * phi)

        diffusion_error = np.abs(A / A_exact)
        pp.plot(phi, diffusion_error, label=f'c={c}')

    pp.xlabel(r'$\phi$ (rad)')
    pp.ylabel(r'Diffusion error $|A/A_{exact}|$')
    pp.title(f'Diffusion error ({method} scheme)')
    pp.grid(True)
    pp.legend()
    pp.xlim([-np.pi, np.pi])

    # Dispersion error plot
    pp.subplot(1,2,2)
    for c in c_values:
        if method == 'upwind':
            A = 1 / (1 + c - c * np.exp(-1j * phi))
        elif method == 'centered':
            A = 1 / (1 - 1j * c * np.sin(phi))

        A_exact = np.exp(-1j * c * phi)
        dispersion_error = np.angle(A / A_exact)
        pp.plot(phi, dispersion_error, label=f'c={c}')

    pp.xlabel(r'$\phi$ (rad)')
    pp.ylabel('Dispersion error (rad)')
    pp.title(f'Dispersion error ({method} scheme)')
    pp.grid(True)
    pp.legend()
    pp.xlim([-np.pi, np.pi])

    pp.tight_layout()
    pp.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Backward Euler solver with Upwind and Centered schemes")
    parser.add_argument('-c', type=float, default=0.25, help="Courant number")
    parser.add_argument('-a', type=float, default=1.0, help="Advection velocity")
    parser.add_argument('--method', type=str, default="upwind", help="Space discretization method")
    parser.add_argument('--tmax', type=float, default=1.0)
    parser.add_argument('-m', type=int, default=100, help="Number of discretization points")
    parser.add_argument('--solver', type=str, default='fourier', help="Solver method for backward Euler ('fourier' or 'direct')") 
       
    args = parser.parse_args()
    
    simulate(c=args.c, a=args.a, t_max=args.tmax, M=args.m, method=args.method)