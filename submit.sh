#BSUB -W 3:00
#BSUB -n 16
#BSUB -e %J.err
#BSUB -o %J.out

mpiexec -np 16 python A3.py
