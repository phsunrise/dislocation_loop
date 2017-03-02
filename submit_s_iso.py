import os

nprocs = 16 

## first generate the filelist
val = os.system("python s_iso.py -d")
if not val == 0:
    print "Aborting submission..."
    sys.exit(1)

for rank in xrange(nprocs):
    with open("s_%03d.sbatch"%rank, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("\n")
        f.write("#SBATCH --job-name=s_%03d\n"%rank)
        f.write("#SBATCH --output=s_%03d.out\n"%rank)
        f.write("#SBATCH --error=s_%03d.err\n"%rank)
        f.write("#SBATCH --time=1:00:00\n")
        f.write("#SBATCH --qos=normal\n")
        f.write("#SBATCH --partition=iric\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks-per-node=1\n")
        if rank == nprocs-1:
            f.write("#SBATCH --mail-type=END\n")
            f.write("#SBATCH --mail-user=phsun@stanford.edu\n")
        f.write("\n")
        f.write("python s_iso.py -n %d -p %d\n"%(nprocs, rank))

    os.system("sbatch s_%03d.sbatch" % rank)
