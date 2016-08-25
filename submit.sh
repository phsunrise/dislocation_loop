#BSUB -W 4:00
#BSUB -n 16 
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[hosts=1]"
#BSUB -q medium

mpiexec -np 16 python s.py
