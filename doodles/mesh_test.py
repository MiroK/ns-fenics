from dolfin import *
from mpi4py import MPI as mpi
from commands import getoutput

mesh = UnitSquareMesh(100, 100)

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
  if found[i]:
    root = i
    break

print found, root
value = comm.bcast(value, root=root)
print value

#print errornorm(u, u)

import os 

def getMyMemoryUsage():                                                    
  mypid = os.getpid()                                                          
  mymemory = float(getoutput('ps -o rss %s' % mypid).split()[1])
  total_memory = MPI.sum(mymemory)
  return mymemory, total_memory  

print MPI.process_number(), getMyMemoryUsage()

def localize(f):
  if MPI.process_number() == 0:
    return f
  else:
    def void(*args, **kwargs):
      pass
    return void

@localize
def myprint(a):
  print a

myprint(77)
