#BSUB -W 24:00 
#BSUB -n 750 
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q bulletmpi-large 

mpiexec -np 750 python atoms_amplitude_poly.py
