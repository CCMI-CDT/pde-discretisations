from firedrake import *
import matplotlib.pyplot as plt
import os

# Time derivative
def dudt(u, u_, dt):
    return ((u - u_)/dt)

t_max = 0.1
dt = 1/200 # Time discretisation large enough to have negligible effect on error

k = 2 * pi # Frequency 
omega = k**2

cells = np.arange(5, 50) # Increasing set of basis functions
final_errors = []

for num_cells in cells: 
    mesh = PeriodicIntervalMesh(num_cells, 1.0) # Periodic domain (1D)
    x = SpatialCoordinate(mesh)

    V = VectorFunctionSpace(mesh, "CG", 1, dim=2) # Complex plane as \R^2

    v = TestFunction(V)

    psi_0 = Function(V)
    psi_0.interpolate(as_vector([cos(k*x[0]), sin(k*x[0])]))
    psi_r0, psi_i0 = split(psi_0)

    psi = Function(V)
    psi.assign(psi_0)
    psi_r, psi_i = split(psi) # Split into real and complex parts

    # Implicit midpoint scheme schrodinger/schrodinger_errors.py
    psi_rh = (psi_r0 + psi_r) / 2
    psi_ih = (psi_i0 + psi_i) / 2

    # Residual in real and complex variables
    F1 = (inner(dudt(psi_i, psi_i0, dt), v[0]) + inner(grad(psi_rh), grad(v[0]))) * dx
    F2 = (inner(dudt(psi_r, psi_r0, dt), v[1]) - inner(grad(psi_ih), grad(v[1]))) * dx
    F = F1 + F2

    t = 0.0
    t_ufl = Constant(0.0)                 
    psi_true = Function(V, name="psi_true")

    times = []
    l2_errors = []

    # step through time and record solution
    while t < t_max - dt/2:
        t += dt
        solve(F == 0, psi)
        psi_0.assign(psi)

        t_ufl.assign(t)

        #True solution 
        psi_true.interpolate(as_vector([cos(k * x[0] -  omega * t_ufl), sin(k * x[0] -  omega *  t_ufl)]))

        err_l2 = norms.errornorm(psi, psi_true) # L2 error 

        times.append(t)
        l2_errors.append(err_l2)
    
    # Only return final error
    final_errors.append(l2_errors[-1])

final_errors = np.array(final_errors)

# Reference line in log-space with fixed slope -2 through (xd, yd)
log_x = np.log(cells)
log_y = np.log(final_errors)
xd, yd = log_x[0], log_y[0]
slope = -2.0

# log(y_ref) = yd + slope * (log(x) - xd)
log_y_ref = yd + slope * (log_x - xd)
y_ref = np.exp(log_y_ref)

plt.figure()
plt.loglog(cells, final_errors, '-', label='L2 error')
plt.loglog(cells, y_ref, '--', label='Reference slope -2 through (xd, yd)')
plt.xlabel('Number of cells')
plt.ylabel('Final time L2 error')
plt.legend()

results_dir = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(results_dir, exist_ok=True)

plot_path = os.path.join(results_dir, f"L2_error_{dt}.png")
plt.savefig(plot_path)
plt.show()