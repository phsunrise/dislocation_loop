bsub -e Al_A3.err -o Al_A3.out -W 2:00 -q long -n 16 mpiexec -np 16 python A3.py
