import numpy as np
import matplotlib.pyplot as plt

output = []
for i in xrange(20):
    for j in xrange(20):
        phi = 2.*np.pi*np.random.random_sample()
        th = np.arccos(2.*np.random.random_sample()-1)
        output.append([th, phi])

output = np.array(output)
np.save("uniformsphere.npy", output)

for th, phi in output:
    plt.scatter(np.sin(th)*np.cos(phi),\
                np.sin(th)*np.sin(phi))
plt.show()
