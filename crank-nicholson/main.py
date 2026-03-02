
import numpy as np
import matplotlib.pyplot as plt
import scipy


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


def dump(t, u, d_count):
    plt.figure()
    plt.plot(x, u)
    plt.ylim(-0.1, 1.1)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(f"Advection solution at t = {t:.2f}")
    plt.savefig(f"solution_{d_count:03d}.png")
    plt.close()
    


def main():

    col = np.zeros(n)
    col[0] = 1
    col[1] = -c/4
    col[-1] = c/4

    A = scipy.linalg.circulant(col)

    B = A.T

    




    








if __name__=='__main__':
    main()
