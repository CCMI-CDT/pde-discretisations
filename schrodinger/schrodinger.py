from firedrake import *
import matplotlib.pyplot as plt
# from firedrake.pyplot import tricontour

def dudt(u, u_, dt):
    return ((u - u_)/dt)

num_cells = 100 # Number of cells in each direction for the mesh
checkpointing = False # Checkpointing flag
t_max = 1.0
dt = 1/num_cells

# Create a unit square mesh
mesh = PeriodicIntervalMesh(num_cells, 1.0)
x = SpatialCoordinate(mesh)

V = VectorFunctionSpace(mesh, "CG", 1, dim=2)


v = TestFunction(V)


# phi_rn, phi_in = split(Function(V))
psi_0 = Function(V)
psi_0.interpolate(as_vector([sin(4*pi*x[0]), sin(4*pi*x[0])]))
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
sol_history = []
t_values = np.linspace(0, num_steps*dt, num_steps+1)
t_array = []
t = 0
while t < t_max - dt/2:
    t_array.append(t)
    t += dt
    solve(F == 0, psi)
    psi_0.assign(psi)
    # extract real component at each degree of freedom
    sol = [psi.dat.data[i][0] for i in range(len(psi.dat.data))]
    sol_history.append(sol)

# convert history to numpy array shape (ntime, npoints)
sol_history = np.array(sol_history)

# get spatial coordinates for plotting
coords = mesh.coordinates.dat.data[0:-1:2]

print(len(coords))
print(len(t_array))

# create contour map: t along vertical axis (rows), x along horizontal
T, X = np.meshgrid(np.array(t_array), coords, indexing='ij')
print(np.shape(X), np.shape(T), np.shape(sol_history))
plt.figure()
cs = plt.contourf(X, T, sol_history, cmap='viridis')
plt.colorbar(cs, label='real(ψ)')
plt.xlabel('x')
plt.ylabel('t')
plt.title('Solution contour in x–t plane')
plt.show()


        

