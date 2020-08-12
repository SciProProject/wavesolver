"""Module used to solve onedimensional timeindependent Schroedinger equations
   using a set of data from a given file.
"""
import scipy.linalg as sl
import scipy.interpolate as si
import numpy as np


def run(filedir):
    _input_reader(filedir)


def _input_reader(filedir):
    """Collects data from the file the user iputs whitch is used to solve
       the wavefunction in the solvers module.

       Args:
           filedir: directory of the inputfile

       Returns:
           mass, interpdata, methode, x_inp, pot, eigmin, eigmax and
           nump used in solving the schroedinger equation.
    """
    with open(filedir, "r") as fp:
        initlist = fp.readlines()
        datalist = []
        for tup in enumerate(initlist):
            part1 = tup[1].partition('#')
            part2 = part1[0].partition('\t')
            part3 = part2[0].partition('\n')
            datalist.append(part3[0])

    mass = float(datalist[0])
    interpdata = datalist[1].split(' ')
    interpdata[:2] = [float(item) for item in interpdata[:2]]
    interpdata[2] = int(interpdata[2])
    eigmin, eigmax = datalist[2].split(' ')
    eigmin = float(eigmin)
    eigmax = float(eigmax)
    methode = datalist[3]
    nump = float(datalist[4])
    x_inp = []
    pot = []
    for i in range(5, len(datalist)):
        xval, potval = datalist[i].split(' ')
        x_inp.append(float(xval))
        pot.append(float(potval))
    return mass, interpdata, methode, x_inp, pot, eigmin, eigmax, nump


def potential_interpolate(datalist):
    """interpolates the potentail for a new range of x values from a
       given set of potentails.

    Args:
        datalist = list created by the communicator modul function
        input_reader.

    Returns:
        the interpolated potentials and their x values as ndarray.
    """
    x_min, x_max, npoint = _input_reader(filedir)[1]
    methode, x_inp, pot = _input_reader(filedir)[2:5]
    if methode == "linear":
        func = si.interp1d(x_inp, pot)
        x_val = np.linspace(x_min, x_max, npoint)
        potx = func(x_val)
    if methode == "cspline":
        func = si.CubicSpline(x_inp, pot)
        x_val = np.linspace(x_min, x_max, npoint)
        potx = func(x_val)
    if methode == "polynomial":
        func = si.KroghInterpolator(x_inp, pot)
        x_val = np.linspace(x_min, x_max, npoint)
        potx = func(x_val)
    return x_val, potx


def schroedinger_equation_solver(datalist):
    """Solves the schroedingerequation from a given set of potentials.

    Args:
        datalist = list created by the communicator modul function
        input_reader.

    Returns:
        Eigenvalues and their normalized wavefunctions as ndarray.
    """
    mass = _input_reader(filedir)[0]
    x_min, x_max, npoint = _input_reader(filedir)[1]
    step = (x_max - x_min) / npoint
    potx = potential_interpolate(datalist)[1]
    coeff = 1 / (mass * step**2)
    # Calculation of the eigenvalues and wavefunctions
    maindia = np.array([])
    for i, val in enumerate(potx):
        elem = coeff + val
        maindia = np.append(maindia, elem)
    offdia = np.full(len(potx) - 1, -coeff / 2)
    eigval, wavef = sl.eigh_tridiagonal(maindia, offdia)
    # Normalization of the wavefunctions
    normw = np.empty((len(wavef), len(wavef)))
    for i in range(len(wavef[0])):
        norm = step * sum(wavef[:, i]**2)
        normw[:, i] = wavef[:, i] / norm**0.5
    return eigval, normw
