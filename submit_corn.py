import os

nproc = 16
for i in xrange(nproc):
    with open("submit_%d.txt"%i, 'w') as f:
        f.write("#$ -N job%d\n"%i)
        f.write("#$ -cwd\n")
        f.write("#$ -j y\n")
        f.write("python s_parallel.py %d %d" % (nproc, i))
    os.system("qsub submit_%d.txt"%i)
    os.system("rm submit_%d.txt"%i)
