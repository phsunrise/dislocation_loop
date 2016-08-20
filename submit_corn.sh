#!/bin/bash

#$ -N s_calc
#$ -cwd
#$ -M phsun@stanford.edu
#$ -m besa
#$ -j y

mpiexec -n python s_parallel.py
