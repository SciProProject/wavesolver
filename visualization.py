"""Visualization of the data retrived by the solvers module."""
import numpy as np
import matplotlib.pyplot as plt
import os.path as op


def run(path):
    """
    Interface for user to define plot ranges and scaling and visualizes
    data.

    Args:
        path  directory path form whitch to use data
    """
    scaling = input("Do you want to scale your axes manually? [y/n] ")
    if scaling == 'y':
        scalingy = input("Do you want to scale the y-axis? [y/n] ")
        if scalingy == 'y':
            y_min = input("What do you want as lower end for the y-axis? ")
            y_max = input("What do you want as top end for the y-axis? ")
            try:
                y_min = float(y_min)
                y_max = float(y_max)
                pary = [y_min, y_max]
            except:
                raise TypeError("""Your inputs {}, {} could not
be converted into float""".format(y_min, y_max))
        else:
            pary = None

        scalingx = input("Do you want to scale the x-axis? [y/n] ")
        if scalingx == 'y':
            x_min = input("What do you want as lower end for the x-axis? ")
            x_max = input("What do you want as top end for the x-axis? ")
            try:
                x_min = float(x_min)
                x_max = float(x_max)
                parx = [x_min, x_max]
            except:
                raise TypeError("""Your inputs {}, {} could not
be converted into float""".format(x_min, x_max))
        else:
            parx = None

        visualization(path, pary, parx)

    if scaling == 'n':
        visualization(path)


def visualization(path, pary=None, parx=None):
    """Visualizes the wavefunctions and eigenenergies for a given potential
       as well as the expected x values in a graph as well as the diviations
       of the x values. The solvers output files are being used.

       Args:
           path  directory path form whitch to use data
           pary  tuple disclosing the range of the y-axis (y_min, y_max)
           parx  tuple disclosing the range of the x-axis (x_min, x_max)
    """
    name_list = ["energies.dat", "potential.dat", "wavefuncs.dat",
                 "expvalues.dat"]
    try:
        data = []
        for names in enumerate(name_list):
            directory = op.join(path, names[1])
            data.append(np.loadtxt(directory))
        eigv, pot, wavef, expv = data
    except:
        raise FileNotFoundError("""Data was not found.""")

    # Scaling of the wavefunctions

    scal = 1 / (np.amax(abs(eigv)) + 2)
    for i, y_axis in enumerate(eigv):
        transfwavef = wavef[:, i+1] * scal + y_axis
        wavef[:, i+1] = transfwavef

    # Defining eigenvalues
    eigval = np.empty((len(pot[:, 0]), len(eigv) + 1))
    eigval[:, 0] = pot[:, 0]
    for i, val in enumerate(eigv):
        eigval[:, i+1] = val

    # Defining two subplots in one figure
    (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 6), sharey=True)[1]

    # Ploting potential
    ax1.plot(pot[:, 0], pot[:, 1], color='black')

    # Ploting eigenvalues, wavefunctions and expected vlaues
    ax1.set_prop_cycle('color', ['blue', 'red'])
    for val in enumerate(eigv):
        i = val[0]
        ax1.plot(eigval[:, 0], eigval[:, i+1], color='grey', alpha=0.5)
        ax1.plot(wavef[:, 0], wavef[:, i+1])
        ax1.plot(expv[i, 0], eigv[i], color='green', marker='x',
                 markersize=12, markeredgewidth=1.5)

    # Setting limits
    if pary:
        ax1.set_ylim(pary[0], pary[1])
    if parx:
        ax1.set_xlim(parx[0], parx[1])

    # Labeling first subplot
    ax1.set_xlabel('x [Bohr]')
    ax1.set_ylabel('Energy [Hartree]')
    ax1.set_title(r'Potential, eigenstates, $\langle x \rangle$')

    # Ploting standart diviation of expected values for respectiv eigenvlaues
    for i in range(0, len(eigv)):
        ax2.plot(eigval[:, 0], eigval[:, i+1], color='grey', alpha=0.5)
        ax2.plot(expv[i, 1], eigv[i], color='magenta',
                 marker='+', markersize=15, markeredgewidth=2)

    # Labeling seconde subplot
    sigmax = np.max(expv[:, 1]) + 0.05 + 0.01 * np.max(expv[:, 1])
    ax2.set_xlim(0, sigmax)
    ax2.set_xlabel('x [Bohr]')
    ax2.set_title(r'$\sigma_{x}$')

    plt.show()
