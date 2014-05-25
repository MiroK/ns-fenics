from dolfin import *
from numpy import random

mesh = UnitIntervalMesh(5)
V = FunctionSpace(mesh, 'DG', 0)
V = FunctionSpace(mesh, 'DG', 0)
V = MixedFunctionSpace([V, V])

u = Function(V)
U = u.vector()

for i in range(10):
    x = random.random(V.dim())
    u.vector()[:] = x

    print id(u.vector()), id(U)

    print u.vector().array()
    print U.array()
    print

    x = random.random(V.dim())
    U[:] = x
    print u.vector().array()
    print U.array()
    print
    print





