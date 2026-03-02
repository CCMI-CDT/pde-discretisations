import numpy as np
import matplotlib.pyplot as pp
from scipy.linalg import circulant
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
    eigAinv = np.linalg.eig(Ainv)
    print("Eigenvalues of A^-1:", np.abs(eigAinv))
    if max(np.abs(eigAinv)) > 1.0:
        print("Unstable")
    else:
        print("Stable")

def build_matrix(c, M):
    '''
    Build a circulant matrix for the backward Euler
    centered scheme applied to the 1-D advection equation
    '''
    vals = np.zeros(M)
    vals[0] = 1.0
    vals[1] = c/2.0
    vals[-1] = -c/2.0
    A = circulant(vals)
    stability_analysis(A)
    return A

def backward_euler(A, un):
    un = np.linalg.solve(A, un)
    return un



def dump(x, un, t, dcount, c):
    pp.plot(x, un)
    pp.xlabel('x')
    pp.ylabel('y')
    pp.title(f'Solution at time {t}')
    pp.savefig(f'beup_c{c}_dump{dcount}.jpg')
    pp.close()

def simulate(c, a, t_max=1.0, tdump=0.2, L=1.0, M=100):
    dx = L / M
    x = np.arange(0, L, dx)
    dt = dx*c/a
    un = np.zeros_like(x)
    un[np.where(x>0.5)] = 1.0
    A = build_matrix(c, M)
    t = 0.0
    dcount = 0
    dumpt = 0
    while t < t_max - dt/2:
        t += dt
        un = backward_euler(A, un)

        dumpt += dt
        if dumpt > tdump - dt/2:
            dump(x, un, t, dcount, c)
            dumpt -= tdump
            dcount += 1

for c in [0.25, 0.5, 0.75, 1.0]:
    simulate(c, 1.0)