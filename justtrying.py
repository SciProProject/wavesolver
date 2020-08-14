import numpy as np
import matplotlib.pyplot as plt


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


fig1, ax = plt.subplots(figsize = (4,6)) 

plt.plot(pot[:, 0], pot[:, 1], color = 'black')

ax.set_prop_cycle('color', ['blue', 'red'])
for i in range(0,len(eigv)):
    plt.plot(eigval[:, 0], eigval[:,i+1], color = 'grey', alpha = 0.5)
    plt.plot(wavef[:, 0], wavef[:, i+1])
    
plt.xlabel('x [Bohr]')
plt.ylabel('Energy [Hartree]')
plt.title(r'Potential, eigenstates, $\langle x \rangle$')

#/home/fabian/scipro/input/infpot.inp
#/home/fabian/scipro/input/finitepot.inp
#/home/fabian/scipro/input/cubspli.inp