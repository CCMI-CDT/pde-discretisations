from firedrake import *

from tqdm import tqdm

def heat_equation():
    mesh = UnitSquareMesh(30, 30)
    V = FunctionSpace(mesh, "CG", 1)

    dt = 0.0001

    u = TrialFunction(V)
    v = TestFunction(V)
    init = Function(V)
    output_function = Function(V)

    sigma = 0.01
    means = [[0.25, 0.75], [0.475, 0.475], [0.475, 0.525], [0.525, 0.475], [0.525, 0.525]]
    x, y = SpatialCoordinate(V.mesh())
    init.interpolate(
        (1 / (2 * pi * sigma**2)) * exp(-((x - means[0][0])**2 + (y - means[0][1])**2) / (2 * sigma**2)) +
        (1 / (2 * pi * sigma**2)) * exp(-((x - means[1][0])**2 + (y - means[1][1])**2) / (2 * sigma**2)) -
        (1 / (2 * pi * sigma**2)) * exp(-((x - means[2][0])**2 + (y - means[2][1])**2) / (2 * sigma**2)) -
        (1 / (2 * pi * sigma**2)) * exp(-((x - means[3][0])**2 + (y - means[3][1])**2) / (2 * sigma**2)) +
        (1 / (2 * pi * sigma**2)) * exp(-((x - means[4][0])**2 + (y - means[4][1])**2) / (2 * sigma**2))
    )

    a = (dt * inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(init, v) * dx

    problem = LinearVariationalProblem(a, L, output_function)
    solver = LinearVariationalSolver(problem, solver_parameters={'ksp_type': 'preonly', 'pc_type': 'lu'})

    outfile = VTKFile("heat_equation_sol.pvd")

    for i in tqdm(range(10)):
        solver.solve()

        outfile.write(output_function, time=i*dt)

        init.assign(output_function)


if __name__=='__main__':
    heat_equation()
