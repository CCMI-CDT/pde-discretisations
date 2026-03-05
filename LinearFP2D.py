from firedrake import * 
import matplotlib.pyplot as plt
from firedrake.pyplot import tripcolor, tricontour
import time 


n = 100
Mesh = UnitSquareMesh(n,n) #2 spatial dimension
Ex_Mesh = ExtrudedMesh(Mesh,n,(1/30),extrusion_type='uniform') #2 spatial, 1 time dimension

##1D spatial mesh
V = FunctionSpace(Mesh, "CG",2) #for test & trial functions
Vvec = VectorFunctionSpace(Mesh, "CG", 2) #For drift 

#2D spatial mesh
#horiz_elt = FiniteElement("CG",triangle,1)
#vert_elt = FiniteElement("CG",interval,1)
#elt = TensorProductElement(horiz_elt,vert_elt)
#V_2 = FunctionSpace(Mesh,elt)
#Vvec_2 = VectorFunctionSpace(Mesh,elt)


#Can't be expressed as a(u,v) = L[v]
u_ = Function(V, name="Velocity")
u = Function(V, name="VelocityNext")
v = TestFunction(V)

x = SpatialCoordinate(Mesh)
b = Function(Vvec)
b.interpolate(as_vector((x[0],x[1])))
D = Constant(0)
chi = interpolate(x[0]* (1.0 - x[0]) * x[1] * (1.0 - x[1]),V) #cutoff
#Inital condition
b_tilde = Function(Vvec)
b_tilde.project(chi*b)
D_tilde = Function(V)
D_tilde.project(chi * D)

u.assign(project(x[1]/100000, V)) #Setting initial conditions 
u_.assign(u)

timestep=1.0/n

F = (inner(b_tilde*u,nabla_grad(v)) + inner(D_tilde*nabla_grad(u),nabla_grad(v)))*dx - inner(v,(u - u_)/timestep)*dx


outfile = VTKFile("Fokker-Planck.pvd")
outfile.write(u)


t = 0.0
end = 5
while (t <= end):
    solve(F==0,u)
    u_.assign(u)
    t += timestep
    outfile.write(u)
    fig, axes = plt.subplots()
    colors = tripcolor(u, axes=axes)
    fig.colorbar(colors)
    plt.savefig("Linear_FP.png", dpi=200)
    print("Saved plot to Linear_FP.png")
    time.sleep(0.025)
