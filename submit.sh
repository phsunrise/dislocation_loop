#BSUB -W 3:00
#BSUB -n 16
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[hosts=1]"
#BSUB -q medium

mpiexec -np 16 python A3.py
