import numpy as np
import matplotlib.pyplot as plt

xx = np.logspace(-2.,2.,20)
y1 = np.zeros(len(xx))
for i, x in enumerate(xx):
    if x < 1.: 
        y1[i] = 1./x**2
    else:
        y1[i] = 1./x**4
y2 = y1*xx/5.

plt.loglog(xx, y1)
plt.loglog(xx, y2)
plt.show()
