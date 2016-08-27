#!/bin/bash

#$ -N test
#$ -cwd
#$ -M phsun@stanford.edu
#$ -m bea
#$ -pe orte 8

## this joins normal output and error output into one file
#$ -j y

mpirun -np 8 python test.py
