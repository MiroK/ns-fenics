from dolfin import *
mesh = UnitSquareMesh(2,2)
V = FunctionSpace(mesh,"CG",1)
u0 = interpolate(Constant(1), V)
u1 = Function(V)

U = u0.vector()
print "Same:", U.id()==u0.vector().id()
print "Same:", all((U.array()-u0.vector().array()) < DOLFIN_EPS)

u0.assign(u1)
print "Same:", U.id()==u0.vector().id()
print "Same:", all((U.array()-u0.vector().array()) < DOLFIN_EPS)
