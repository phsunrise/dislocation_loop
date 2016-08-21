from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

print "hello from %d/%d" % (rank, nprocs)
with open("out%d.txt"%rank, 'w') as f:
    f.write("rank=%d"%rank)
