bsub -W 2:00 -q long -o Al_s.out -e Al_s.err mpiexec -n 16 python s_parallel.py
