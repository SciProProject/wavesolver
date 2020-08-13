import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sl


x = np.linspace(-2, 2, 100)
eigv = np.loadtxt('output/energies.dat')
pot = np.loadtxt('output/potential.dat')
wavef = np.loadtxt('output/wavefuncs.dat')

eigval = np.empty((len(pot[:,0]), len(eigv)+1))
eigval[:, 0] = pot[:, 0]
for i, val in enumerate(eigv):
    eigval[:, i+1] = val


fig1, ax = plt.subplots()

plt.plot(pot[:, 0], pot[:, 1], color = 'black')

ax.set_prop_cycle('color', ['blue', 'red'])
for i in range(0,len(eigv)):
    plt.plot(eigval[:, 0], eigval[:,i+1], color = 'grey', alpha = 0.5)
    #plt.plot(wavef[:, 0], wavef[:, val])

#/home/fabian/scipro/input/schrodinger.inp