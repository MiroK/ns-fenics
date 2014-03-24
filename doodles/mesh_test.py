from dolfin import *
from mpi4py import MPI as mpi

mesh = UnitSquareMesh(20, 20)

V = VectorFunctionSpace(mesh, 'CG', 1)
u = interpolate(Expression(('x[0]', '1')), V)

value = None
found = 0

try:
  print 'Process', MPI.process_number()
  value = u(1, 0)
  found = 1
except:
  pass

comm = mpi.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

found = comm.allgather(found)

assert len(found) == MPI.num_processes()

root = -1
for i in range(len(found)):
  if found:
    root = i
    break

value = comm.bcast(value, root=root)
print value

print errornorm(u, u)
