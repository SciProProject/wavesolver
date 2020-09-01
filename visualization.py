"""Visualization of the data retrived by the solvers module."""
import numpy as np
import matplotlib.pyplot as plt


def run():
    visualization()


def visualization():
    eigv = np.loadtxt('output/energies.dat')
    pot = np.loadtxt('output/potential.dat')
    wavef = np.loadtxt('output/wavefuncs.dat')
    for i, val in enumerate(eigv):
        b = wavef[:, i+1] + val
        wavef[:, i+1] = b


    eigval = np.empty((len(pot[:,0]), len(eigv)+1))
    eigval[:, 0] = pot[:, 0]
    for i, val in enumerate(eigv):
        eigval[:, i+1] = val


    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize = (8,6)) 

    ax1.plot(pot[:, 0], pot[:, 1], color = 'black')

    ax1.set_prop_cycle('color', ['blue', 'red'])
    for i in range(0,len(eigv)):
        ax1.plot(eigval[:, 0], eigval[:,i+1], color = 'grey', alpha = 0.5)
        ax1.plot(wavef[:, 0], wavef[:, i+1])
    
    ax1.set_xlabel('x [Bohr]')
    ax1.set_ylabel('Energy [Hartree]')
    ax1.set_title(r'Potential, eigenstates, $\langle x \rangle$')
    
    plt.show()
    
#/home/fabian/scipro/input/infpot.inp
#/home/fabian/scipro/input/finitepot.inp
#/home/fabian/scipro/input/cubspli.inp