from dolfin import *
import mpi4py.MPI as mpi

comm = mpi.COMM_WORLD
rank = comm.Get_rank()


def pprint(arg):
    if rank == 0:
        print arg

pprint('Hello')

if rank == 0:
    n = raw_input('Set n:')
else:
    n = None
n = comm.bcast(n, root=0)

print rank, n


timer = Timer('WHAAAT')
timer.start()
for i in range(100000*(rank+1)):
    j = i
time = timer.stop()

print rank, time

time = timing('WHAAAT')
time = MPI.sum(time)


print rank, time
