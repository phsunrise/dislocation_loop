#BSUB -W 3:00
#BSUB -n 128 
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q psanaq 

source /reg/g/psdm/etc/ana_env.sh
#mpirun -np 63 python atoms_s.py
mpirun -np 128 python atoms_amplitude.py
#mpirun -np 16 python test.py
