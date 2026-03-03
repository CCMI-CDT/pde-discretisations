import numpy as np
import matplotlib.pyplot as plt
import scipy


def dump(x, u, time, d_count):
    plt.figure()
    plt.plot(x, u)
    plt.ylim(-0.1, 1.1)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(f"Advection solution at t = {time:.2f}")
    plt.savefig(f"solution_{d_count:03d}.png")
    plt.close()
    

def main():

    L = 1
    n = 100
    dx = L/n
    x = np.arange(0. ,L, dx)
    un = np.where(x<0.5,1,0)
    c = 1
    a = 1
    dt = dx*c/a
    t = 0.
    t_max = 2
    t_dump = 0.2
    dump_t = 0.
    d_count = 0

    col = np.zeros(n)
    col[0] = 1
    col[1] = -c/4
    col[-1] = c/4

    Lmat = scipy.linalg.circulant(col)
    Rmat = Lmat.T

    while t < (t_max - dt / 2):
        # Lu_{n+1} = Ru_{n}
        t += dt
        # TODO: use roll instead of matrix multiplication
        rhs = Rmat @ un

        # TODO: FFT the matrix 
        # TODO: Compare with upwind CN
        un = scipy.linalg.solve(Lmat, rhs)
        
        if t >= dump_t - 1e-12:
            dump(x, un, t, d_count)
            dump_t += t_dump
            d_count += 1

if __name__=='__main__':
    main()
