from dolfin import *

rank = MPI.process_number()

min_rank = MPI.min(rank)

print rank, min_rank
