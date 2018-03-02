from fenics import *
import mshr

krylov_method = "gmres"

#
# Solve a stokes flow problem
#
L = 4.0
h = 1.0
Dp = 1.0
mesh = RectangleMesh(Point(-L/2.0, -h/2.0),Point(L/2.0,h/2.0) ,40,40, "right/left")
boundfunc = FacetFunction("uint", mesh)
boundfunc.set_all(0)
top = CompiledSubDomain(" x[1]>h/2.0-eps && on_boundary",h=h,L=L,eps=1.0e-12)
top.mark(boundfunc,1)
bot = CompiledSubDomain("-x[1]>h/2.0-eps && on_boundary",h=h,L=L,eps=1.0e-12)
bot.mark(boundfunc,2)
inlet = CompiledSubDomain("-x[0]>L/2.0-eps && on_boundary",h=h,L=L,eps=1.0e-12)
inlet.mark(boundfunc,3)
outlet = CompiledSubDomain("x[0]>L/2.0-eps && on_boundary",h=h,L=L,eps=1.0e-12)
outlet.mark(boundfunc,4)
# The mixed element
V2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
S1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V1 = VectorElement("Lagrange", mesh.ufl_cell(), 1)
S0 = FiniteElement("DG", mesh.ufl_cell(), 0)
TH = V2 * S1
W = FunctionSpace(mesh, TH)
# Set up the boundary conditions
noslip = Constant((0.0, 0.0))
bcs = [
    DirichletBC(W.sub(0),noslip,boundfunc,1),
    DirichletBC(W.sub(0),noslip,boundfunc,2),
    DirichletBC(W.sub(1),0,boundfunc,4),
]
# Define variational problem
(v,   p) = TrialFunctions(W)
(tv, tp) = TestFunctions(W)
n = FacetNormal(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundfunc)
f = Constant((0.0, 0.0))
a = inner(grad(tv), grad(v))*dx + div(tv)*p*dx + tp*div(v)*dx
L = inner(tv, f)*dx \
      - inner(tv, Constant(Dp)*n)*ds(3)
b = inner(grad(tv), grad(v))*dx + tp*p*dx
# Solve it
U = Function(W)
A, bb = assemble_system(a, L, bcs)
P, btmp = assemble_system(b, L, bcs)
solver = KrylovSolver(krylov_method, "amg")
solver.set_operators(A, P)
solver.solve(U.vector(),bb)
#solve(a==L,U,bcs)
v, p = U.split()

File('fuild_v.pvd') << v

#
# Now let's track some particles through it!
#
import fpl
P = fpl.Particles(1,2)

NT = 100
for t in xrange(NT):
    P.euler_step(v, 0.02)
    print P.x.get_local()
