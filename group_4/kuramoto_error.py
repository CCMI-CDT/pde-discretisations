import matplotlib.pyplot as plt
import numpy as np
from firedrake import *

mesh = PeriodicIntervalMesh(128, 3*pi)
x = SpatialCoordinate(mesh)

V_1 = FunctionSpace(mesh, "CG", 1)
V = V_1 * V_1

T = 1.0
dt_list = [0.1, 0.05, 0.025, 0.0125, 0.00625]

x_dof = Function(V_1).interpolate(x[0])
x_coords_all = x_dof.dat.data_ro.copy()


def run_for_dt(dt_value):
    timestep = Constant(dt_value)
    t_variable = Constant(0.0)

    v1, v2 = TestFunctions(V)
    sol = Function(V)
    u, w = split(sol)

    sol_prev = Function(V)
    u_prev, w_prev = split(sol_prev)

    # initial conditions
    sol_prev.sub(0).interpolate(Constant(0.0))
    sol_prev.sub(1).interpolate(Constant(0.0))

    F = (
        (u - u_prev) / timestep * v1 * dx
        + u * u.dx(0) * v1 * dx
        - u.dx(0) * v1.dx(0) * dx
        - w.dx(0) * v1.dx(0) * dx
        + w * v2 * dx
        + u.dx(0) * v2.dx(0) * dx
        - cos(t_variable) * sin((2.0/3.0) * x[0]) * v1 * dx
        - (2.0/3.0) * sin(t_variable)**2 * sin((2.0/3.0) * x[0]) * cos((2.0/3.0) * x[0]) * v1 * dx
        + (20.0/81.0) * sin(t_variable) * sin((2.0/3.0) * x[0]) * v1 * dx
    )

    nsteps = int(round(T / dt_value))
    t = 0.0

    for _ in range(nsteps):
        t += dt_value
        t_variable.assign(t)

        sol.assign(sol_prev)

        solve(
            F == 0, sol,
            solver_parameters={
                "snes_type": "newtonls",
                "ksp_type": "preonly",
                "pc_type": "lu"
            }
        )

        sol_prev.assign(sol)

    u_exact = Function(V_1).interpolate(sin(Constant(T)) * sin((2.0/3.0) * x[0]))
    u_num = sol.sub(0)

    err_L2 = errornorm(u_exact, u_num, norm_type="L2")

    return err_L2

errors_L2 = []

for i, dt in enumerate(dt_list):
    make_plot = (i == len(dt_list) - 1)
    err_L2 = run_for_dt(dt)
    errors_L2.append(err_L2)

plt.figure(figsize=(8, 5))
plt.loglog(dt_list, errors_L2, marker='o', label='L2 error at T=1')
plt.xlabel('dt')
plt.ylabel('Error')
plt.title('Temporal error analysis')
plt.grid(True, which='both')
plt.legend()
plt.show()