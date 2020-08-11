import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2, 2, 100)
eigv = np.loadtxt('energies.dat')
pot = np.loadtxt('potential.dat')

fig1, ax = plt.subplots()

plt.plot(pot[:, 0], pot[:, 1], color = 'black')

ax.set_prop_cycle('color', ['blue', 'red'])
for val in eigv:
    plt.plot(x, val * x + val)
