#BSUB -W 0:15
#BSUB -n 1024 
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q bulletmpi-large 

mpiexec -np 1024 python s.py
