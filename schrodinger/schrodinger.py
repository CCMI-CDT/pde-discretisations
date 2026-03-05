from firedrake import *
import matplotlib.pyplot as plt
# from firedrake.pyplot import tricontour

def dudt(u, u_, dt):
    return ((u - u_)/dt)

num_cells = 100 # Number of cells in each direction for the mesh
checkpointing = False # Checkpointing flag
t_max = 0.1
dt = 1/1000
k = 2 * pi
omega = k**2

# Create a unit square mesh
mesh = PeriodicIntervalMesh(num_cells, 1.0)
x = SpatialCoordinate(mesh)

V = VectorFunctionSpace(mesh, "CG", 1, dim=2)


v = TestFunction(V)


# phi_rn, phi_in = split(Function(V))
psi_0 = Function(V)
psi_0.interpolate(as_vector([cos(k*x[0]), sin(k*x[0])]))
psi_r0, psi_i0 = split(psi_0)

psi = Function(V)
psi.assign(psi_0)
psi_r, psi_i = split(psi)

F1 = (inner(dudt(psi_i, psi_i0, dt), v[0]) + inner(grad(psi_r), grad(v[0]))) * dx
F2 = (inner(dudt(psi_r, psi_r0, dt), v[1]) - inner(grad(psi_i), grad(v[1]))) * dx
F = F1 + F2
# step through time and record solution
# use a proper scalar for t rather than overwriting with an array
num_steps = int(np.floor(t_max/dt))
sol_history_r = []
sol_history_i = []
t_values = np.linspace(0, num_steps*dt, num_steps+1)
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



plt.savefig("sol.png")


        

