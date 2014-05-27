'''
Current approach in mixed solver is to use the MixedFunctionSpace M with
M.dim() = m and matrices m x m, that is
Am = [[A, 0], [0, 0]], Btm = [[0, Bt], [0, 0]], Bm = [[0, 0], [B, 0]]. Alter-
native is to work with block matrices, that is operate on A, Bt, B (smaller
matrices) and use the system one when necessary. This should reduce memory
usage, but setting up the problem for block matrix in PETSc ... .
This is a comparison of speed between
    [[A, 0], [0, 0]]*[[f], [0]]
and
    [A]*[F].
'''

from dolfin import *

mesh = UnitSquareMesh(100, 100)

# Non-block
V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
M = MixedFunctionSpace([V, Q])

u, p = TrialFunctions(M)
v, q = TestFunctions(M)

a = inner(grad(u), grad(v))*dx
L = inner(Constant((0., 0.)), v)*dx

A, b = assemble_system(a, L)
x = Vector(b)

times = []
for j in range(2):
    timer = Timer('Non-block')
    timer.start()
    i = 0
    while i < 1000:
        i += 1
        A.mult(b, x)
    times.append(timer.stop())
nb_time = sum(times)/len(times)
nb_size = A.size(0)

# Block
u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = inner(Constant((0., 0.)), v)*dx

A, b = assemble_system(a, L)
x = Vector(b)

times = []
for j in range(2):
    timer = Timer('Block')
    timer.start()
    i = 0
    while i < 1000:
        i += 1
        A.mult(b, x)
    times.append(timer.stop())
b_time = sum(times)/len(times)
b_size = A.size(0)

print nb_time, nb_size
print b_time, b_size
print nb_time/b_time, float(nb_size)/b_size

# Non-block multiplication seems be aware of `zero blocks` quite well
# the speed up of block version is NOT WORTH the effert right now!
