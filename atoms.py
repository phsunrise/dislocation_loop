'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os

sample = sys.argv[1]
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

## in loop coordinates, there should be 
for qR in np.linspace(-5., 5., 51):
    if qR == 0:
        continue
    q = qR/R * eq

    ## first calculate the laue term

