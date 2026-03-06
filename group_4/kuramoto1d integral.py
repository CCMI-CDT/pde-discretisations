import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from firedrake import *

#mesh and spaces
mesh = PeriodicIntervalMesh(128, 50.0)
x = SpatialCoordinate(mesh)
V_1 = FunctionSpace(mesh, "CG", 1)
V = V_1 * V_1

v1, v2 = TestFunctions(V)
sol = Function(V)
u, w = split(sol)
sol_prev = Function(V)
u_, w_ = split(sol_prev)

#initial conds
ic_u = sin(2 * pi * x[0] / 50.0)
ic_w = -(2 * pi / 50.0)**2 * sin(2 * pi * x[0] / 50.0)
sol_prev.sub(0).interpolate(ic_u)
sol_prev.sub(1).interpolate(ic_w)

#weak form
timestep = Constant(0.1)
F = (
    (u - u_) / timestep * v1 * dx
    + 0.5 * (u.dx(0))**2 * v1 * dx
    - u.dx(0) * v1.dx(0) * dx
    - w.dx(0) * v1.dx(0) * dx
    + w * v2 * dx
    + u.dx(0) * v2.dx(0) * dx
)

#timeand data storage
t = 0.0
T = 100.0
dt = float(timestep)

#need to hold data for plotting
u_data = [sol_prev.sub(0).dat.data_ro.copy()]
time_steps = [t]

while t < T:
    t += dt
    solve(F == 0, sol)
    sol_prev.assign(sol)
    u_field = sol.sub(0)
    mean_val = assemble(u_field * dx) / assemble(Constant(1.0) * dx(mesh))
    u_field.assign(u_field - mean_val)

    u_data.append(sol.sub(0).dat.data_ro.copy())
    time_steps.append(t)

#plotting, here we need to use matplot
u_array = np.array(u_data)
x = np.linspace(0, 50.0, 128)

""" plt.figure(figsize=(10, 6))
plt.pcolormesh(x_coords, time_steps, u_array, shading='gouraud', cmap='viridis')
plt.colorbar(label='u(x,t)')
plt.title('Kuramoto-Sivashinsky Equation')
plt.xlabel('x')
plt.ylabel('Time (t)')
plt.show() """


fig, ax = plt.subplots()

# Initial plot
line, = ax.plot(x, u_array[0], lw=2)

ax.set_xlim(x.min(), x.max())
ax.set_ylim(u_array.min(), u_array.max())
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")
ax.set_title("Kuramoto–Sivashinsky evolution")

def update(frame):
    line.set_ydata(u_array[frame])
    ax.set_title(f"Kuramoto–Sivashinsky evolution (t step = {frame})")
    return line,

anim = FuncAnimation(
    fig,
    update,
    frames=len(time_steps),
    interval=25,   # milliseconds between frames
    blit=True
)

# Save GIF
anim.save("ks_integral_evolution.gif", writer=PillowWriter(fps=30))

plt.close(fig)