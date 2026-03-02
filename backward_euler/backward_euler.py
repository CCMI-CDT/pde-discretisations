import numpy as np
import matplotlib.pyplot as pp
from scipy.linalg import circulant

L = 1.0
M = 100
dx = L / M
x = np.arange(0, L, dx)
un = np.zeros_like(x)
un[np.where(x>0.5)] = 1.0
c = 1.0
a = 1.0
dt = dx*c/a
t=0.0
t_max = 1.0
tdump = 0.2
dumpt = 0.0
dcount = 0

def build_matrix(c, M):

    vals = np.zeros(M)
    vals[0] = 1.0
    vals[1] = c/2
    vals[-1] = -c/2
    A = circulant(vals)
    
    return A

A = build_matrix(c, M)
def backward_euler(un):
    un = np.linalg.solve(A, un)
    return un

def dump(t):
    pp.plot(x, un)
    pp.xlabel('x')
    pp.ylabel('y')
    pp.title(f'Solution at time {t}')
    pp.savefig(f'beup_c{c}_dump{dcount}.jpg')
    pp.close()

while t < t_max - dt/2:
    t += dt
    un = backward_euler(un)

    dumpt += dt
    if dumpt > tdump - dt/2:
        dump(t)
        dumpt -= tdump
        dcount += 1