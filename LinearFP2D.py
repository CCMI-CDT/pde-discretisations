from firedrake import * 
import matplotlib.pyplot as plt
from firedrake.pyplot import tripcolor, tricontour
import time 


n = 50
Mesh = UnitSquareMesh(n,n) #2 spatial dimension
Ex_Mesh = ExtrudedMesh(Mesh,n,(1/30),extrusion_type='uniform') #2 spatial, 1 time dimension

##1D spatial mesh
V = FunctionSpace(Mesh, "CG",2) #for test & trial functions
Vvec = VectorFunctionSpace(Mesh, "CG", 2) #For drift 

#Messing around with extruding mesh, can ignore 
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
b.interpolate(as_vector((0.5-x[0],0.5-x[1])))
D = Constant(0.5)
chi = interpolate(x[0]* (1.0 - x[0]) * x[1] * (1.0 - x[1]),V) #cutoff
#Inital condition
b_tilde = Function(Vvec)
b_tilde.project(chi*b)
D_tilde = Function(V)
D_tilde.project(chi * D)

u.assign(project(
    exp(-40*((x[0]-0.5)**2 + (x[1]-0.5)**2)) ,
    V
))
mass = assemble(u*dx)
u.assign(u/mass)
u_.assign(u)

timestep=1.0/n

F = (-inner(b*u,nabla_grad(v)) + inner(D*nabla_grad(u),nabla_grad(v)))*dx - inner(v,(u_ - u)/timestep)*dx


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
    print(assemble(u*dx))
    time.sleep(1)
    plt.close()
