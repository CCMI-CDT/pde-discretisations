##### COLAB 
import matplotlib.pyplot as plt
import numpy as np
from firedrake import *

#Setup Mesh and Spaces
mesh = PeriodicRectangleMesh(256, 256, 32*pi, 32*pi) 
x = SpatialCoordinate(mesh)
V_1 = FunctionSpace(mesh, "CG", 1)
V = V_1 * V_1

v1, v2 = TestFunctions(V)
sol = Function(V)
u, w = split(sol)
sol_prev = Function(V)
u_, w_ = split(sol_prev)

#Initial Conditions
ic_u = cos(x[0]/16.0) * (1 + sin(x[0]/16.0)) #sin(2 * pi * x[0] / 50.0)  
ic_w = -(cos(x[0]/16.0) * (1 + 4*sin(x[0]/16.0))) / 128.0
sol_prev.sub(0).interpolate(ic_u)
sol_prev.sub(1).interpolate(ic_w)

#Weak Form
timestep = Constant(0.1)
F = (
    (u - u_) / timestep * v1 * dx
    + u * u.dx(0) * v1 * dx
    - u.dx(0) * v1.dx(0) * dx
    - w.dx(0) * v1.dx(0) * dx
    + w * v2 * dx
    + u.dx(0) * v2.dx(0) * dx
)

#Time-Stepping and Data Storage
t = 0.0
T = 100.0  # Extended time slightly to see chaos develop
dt = float(timestep)

#Lists to hold data for plotting
u_data = [sol_prev.sub(0).dat.data_ro.copy()]
time_steps = [t]

while t < T:
    t += dt
    solve(F == 0, sol)
    sol_prev.assign(sol)
    
    u_data.append(sol.sub(0).dat.data_ro.copy())
    time_steps.append(t)

#Plotting
u_array = np.array(u_data)
x_coords = np.linspace(0, 32*pi, 256)

plt.figure(figsize=(10, 6))
plt.pcolormesh(x_coords, time_steps, u_array, shading='gouraud', cmap='viridis')
plt.colorbar(label='u(x,t)')
plt.title('Kuramoto-Sivashinsky Equation')
plt.xlabel('x')
plt.ylabel('Time (t)')
plt.show()