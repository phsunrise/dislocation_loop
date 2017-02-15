import os

nprocs = 6 

## first generate the filelist
val = os.system("python core_2d.py -d")
if not val == 0:
    print "Aborting submission..."
    sys.exit(1) 

for rank in xrange(nprocs):
    with open("c2d_%03d.sbatch"%rank, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("\n")
        f.write("#SBATCH --job-name=c2d_%03d\n"%rank)
        f.write("#SBATCH --output=c2d_%03d.out\n"%rank)
        f.write("#SBATCH --error=c2d_%03d.err\n"%rank)
        f.write("#SBATCH --time=0:10:00\n")
        f.write("#SBATCH --qos=normal\n")
        f.write("#SBATCH --partition=iric\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks-per-node=1\n")
        if rank == nprocs-1:
            f.write("#SBATCH --mail-type=END\n")
            f.write("#SBATCH --mail-user=phsun@stanford.edu\n")
        f.write("\n")
        f.write("python core_2d.py -n %d -p %d\n"%(nprocs, rank))

    os.system("sbatch c2d_%03d.sbatch" % rank)
