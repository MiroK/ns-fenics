from dolfin import *

mesh = UnitIntervalMesh(5)
V = FunctionSpace(mesh, 'DG', 0)

f = Expression('i', i=0)


def foo():
    u = Function(V)
    U = u.vector()

    for i in range(3):
        f.i = i
        ut = interpolate(f, V)
        u.assign(ut)

        print u.vector().norm('l2')
        print ut.vector().norm('l2')
        print U.norm('l2')
        print


def bar():
    u = Function(V)
    U = u.vector()

    for i in range(3):
        f.i = i
        ut = interpolate(f, V)
        U.zero()
        U.axpy(1, ut.vector())

        print u.vector().norm('l2')
        print ut.vector().norm('l2')
        print U.norm('l2')
        print

foo()
bar()
