from firedrake import *
import matplotlib.pyplot as plt
import os
import argparse
# from firedrake.pyplot import tricontour

def dudt(u, u_, dt):
    return ((u - u_)/dt)

def plot_results(sol_history_r, sol_history_i, mesh, dt, t_array, omega, results_dir):
        # Create results directory with c subfolder if it doesn't exist
    results_dir = os.path.join(results_dir, f'timestep_{dt}')
    os.makedirs(results_dir, exist_ok=True)
    filepath = os.path.join(results_dir, f'schrodinger_solution.jpg')

    # convert history to numpy array shape (ntime, npoints)
    sol_history_r = np.array(sol_history_r)
    sol_history_i = np.array(sol_history_i)

    # get spatial coordinates for plotting
    coords = mesh.coordinates.dat.data[0:-1:2]

    # create contour map: t along vertical axis (rows), x along horizontal
    T, X = np.meshgrid(np.array(t_array), coords, indexing='ij')

    fig, ax = plt.subplots(3, 3, figsize=(8,8))
    # plot real solution
    contours_1 = ax[0][0].contourf(X, T, sol_history_r)
    ax[0][0].set_title("Numerical Solution real(ψ)")
    fig.colorbar(contours_1, ax=ax[0][0], shrink=0.625)

    # plot imaginary solution
    contours_2 = ax[1][0].contourf(X, T, sol_history_i)
    ax[1][0].set_title("Numerical Solution imag(ψ)")
    fig.colorbar(contours_2, ax=ax[1][0], shrink=0.625)

    # plot magnitude
    contours_3 = ax[2][0].contourf(X, T, np.sqrt(np.pow(sol_history_i,2) + np.sqrt(np.pow(sol_history_r,2))))
    ax[2][0].set_title("Numerical Solution mag(ψ)")
    fig.colorbar(contours_3, ax=ax[2][0], shrink=0.625)

    contours_4 = ax[0][1].contourf(X, T, np.cos(k*X-omega*T))
    ax[0][1].set_title("Analytical Solution real(ψ)")
    fig.colorbar(contours_4, ax=ax[0][1], shrink=0.625)

    contours_5 = ax[1][1].contourf(X, T, np.sin(k*X-omega*T))
    ax[1][1].set_title("Analytical Solution imag(ψ)")
    fig.colorbar(contours_5, ax=ax[1][1], shrink=0.625)

    contours_6 = ax[2][1].contourf(X, T, np.sqrt(np.pow(np.sin(k*X-omega*T),2) +np.pow(np.cos(k*X-omega*T),2)))
    ax[2][1].set_title("Analytical Solution mag(ψ)")
    fig.colorbar(contours_6, ax=ax[2][1], shrink=0.625)

    contours_7 =ax[0][2].contourf(X, T, np.abs(sol_history_r - np.cos(k*X-omega*T)))
    ax[0][2].set_title("Point-wise error real(ψ)")
    fig.colorbar(contours_7, ax=ax[0][2], shrink=0.625)

    contours_8 = ax[1][2].contourf(X, T, np.abs(sol_history_i - np.sin(k*X-omega*T)))
    ax[1][2].set_title("Point-wise error imag(ψ)")
    fig.colorbar(contours_8, ax=ax[1][2], shrink=0.625)

    contours_9 = ax[2][2].contourf(X, T, np.abs(np.sqrt(np.pow(sol_history_i,2) + np.sqrt(np.pow(sol_history_r,2))) - np.sqrt(np.pow(np.sin(k*X-omega*T),2) +np.pow(np.cos(k*X-omega*T),2))))
    ax[2][2].set_title("Point-wise error mag(ψ)")
    fig.colorbar(contours_9, ax=ax[2][2], shrink=0.625)

    plt.tight_layout()
    plt.savefig(filepath)

if __name__ == "__main__":


    opts = PETSc.Options()
    num_cells = opts.getInt("--cells", default=100) # Number of cells
    t_max = opts.getReal("--tmax", default=0.1) # Time period
    dt = opts.getReal("--timestep", default = 0.0010) # Timestep
    k = opts.getReal("-k", default = 2.0*pi)
    discretization = opts.getString("-d", default="imp") #Type of discretization. imp for Implicit midpoint and be for Backward Euler
    omega = k ** 2

    # Create a periodic unit interval mesh 
    mesh = PeriodicIntervalMesh(num_cells, 1.0)
    x = SpatialCoordinate(mesh)

    # Create a vector function space
    V = VectorFunctionSpace(mesh, "CG", 1, dim=2)

    # Create test and trial functions
    v = TestFunction(V)
    psi_0 = Function(V) 
    psi_0.interpolate(as_vector([cos(k*x[0]), sin(k*x[0])]))
    psi_r0, psi_i0 = split(psi_0)

    psi = Function(V)
    psi.assign(psi_0)
    psi_r, psi_i = split(psi)

    if discretization == "be":
        # Residual in real and complex variables
        F1 = (inner(dudt(psi_i, psi_i0, dt), v[0]) + inner(grad(psi_r), grad(v[0]))) * dx
        F2 = (inner(dudt(psi_r, psi_r0, dt), v[1]) - inner(grad(psi_i), grad(v[1]))) * dx

    elif discretization == "imp":
        # Implicit midpoint scheme 
        psi_rh = (psi_r0 + psi_r) / 2
        psi_ih = (psi_i0 + psi_i) / 2

        # Residual in real and complex variables
        F1 = (inner(dudt(psi_i, psi_i0, dt), v[0]) + inner(grad(psi_rh), grad(v[0]))) * dx
        F2 = (inner(dudt(psi_r, psi_r0, dt), v[1]) - inner(grad(psi_ih), grad(v[1]))) * dx
    
    F = F1 + F2

    sol_history_r = []
    sol_history_i = []
    t_array = []
    t = 0
    while t < t_max - dt/2:
        t_array.append(t)
        t += dt
        solve(F == 0, psi)
        psi_0.assign(psi)
        # extract real component at each degree of freedom
        sol_r = [psi.dat.data[i][0] for i in range(len(psi.dat.data))]
        sol_history_r.append(sol_r)
        sol_i = [psi.dat.data[i][1] for i in range(len(psi.dat.data))]
        sol_history_i.append(sol_i)

    plot_results(sol_history_r, sol_history_i, mesh, dt, t_array, omega, results_dir="schrodinger_results/")

        

