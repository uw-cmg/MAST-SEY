import os
import numpy as np

energies = [50, 100, 150, 200, 250, 500, 700, 900]
tey = []
sey = []

for energy in energies:
    print(energy)
    os.system("mast_sey -e {} -m 1000 -ins -ph > ttt".format(energy))
    res = np.loadtxt('ttt', comments='#')
    tey.append(res[1])
    sey.append(res[2])

data = np.column_stack([energies, tey, sey])
np.savetxt('results', data, fmt=['%d','%.2f','%.2f'])
