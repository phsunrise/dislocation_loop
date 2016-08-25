#BSUB -n 8
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[hosts=1]"
#BSUB -q express 

mpiexec -np 8 python test.py
