from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

print rank, nprocs
with open("test/%d.out"%rank, 'w') as f:
    f.write(str(rank)+' '+str(nprocs)+'\n')
