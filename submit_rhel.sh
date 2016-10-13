#BSUB -W 1:00 
#BSUB -n 743 
#BSUB -q bulletmpi-large 

#/opt/lsf-openmpi/1.8.1/bin//mpirun -np 743 python atoms_amplitude_poly.py
#/opt/lsf-openmpi/1.8.1/bin//mpirun -np 743 python test.py
mpirun -np 743 python test.py
