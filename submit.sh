#BSUB -W 3:00
#BSUB -n 11
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[hosts=1]"
#BSUB -q medium

mpiexec -np 11 python A3.py
