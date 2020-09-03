"""Module used to solve onedimensional timeindependent Schroedinger equations
   using a set of data from a given file.
"""
import scipy.linalg as sl
import scipy.interpolate as si
import numpy as np


def run(filedir):
    mass, interpdata, methode, potar, eigmin, eigmax, nump = _input_reader(filedir)[:]
    x_val, potx = _potential_interpolate(interpdata, methode, potar)
    eigval, normwf = _schroedinger_equation_solver(mass, interpdata,
                                                   potx, eigmin, eigmax)
    mat = np.empty((len(x_val), 2))
    mat[:, 0] = x_val
    mat[:, 1] = potx
    np.savetxt("output/potential.dat", mat)
    mat2 = np.empty((len(x_val), len(eigval) + 1))
    mat2[:, 0] = x_val
    for i, val in enumerate(np.transpose(normwf)):
        mat2[:, i + 1] = val
    np.savetxt("output/wavefuncs.dat", mat2)
    mat3 = np.empty((len(eigval), 1))#
    for i, val in enumerate(eigval):
        mat3[i] = val
    np.savetxt("output/energies.dat", mat3)
    mat4 = np.empty((len(eigval), 2))
    expx, sigma = expected_values(interpdata, normwf, x_val)
    mat4[:, 0] = expx
    mat4[:, 1] = sigma
    np.savetxt("output/expvalues.dat", mat4)


def _input_reader(filedir):
    """Collects data from the file the user iputs whitch is used to solve
       the wavefunction in the solvers module.

       Args:
           filedir: directory of the inputfile

       Returns:
           mass, interpdata, methode, potar, eigmin, eigmax and
           nump used in solving the schroedinger equation.
    """
    data = np.genfromtxt(filedir, dtype=str, delimiter="/n")
    mass = data[0].astype(np.float)
    x_min, x_max, npoint = data[1].split(' ')
    x_min = float(x_min)
    x_max = float(x_max)
    npoint = int(npoint)
    interpdata = [x_min, x_max, npoint]
    eigmin, eigmax = data[2].split(' ')
    eigmin = int(eigmin)
    eigmax = int(eigmax)
    methode = data[3]
    nump = data[4]
    nump = int(nump)
    potar = np.empty((nump, 2), dtype=float)
    for i in range(0, nump):
        x_inp, pot = data[5+i].split(' ')
        potar[i, 0] = x_inp
        potar[i, 1] = pot

    return mass, interpdata, methode, potar, eigmin, eigmax, nump


def _potential_interpolate(interpdata, methode, potar):
    """interpolates the potentail for a new range of x values from a
       given set of potentails.

    Args:
        datalist = list created by the communicator modul function
        input_reader.

    Returns:
        the interpolated potentials and their x values as ndarray.
    """
    if methode == "linear":
        func = si.interp1d(potar[:, 0], potar[:, 1])
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    if methode == "cspline":
        func = si.CubicSpline(potar[:, 0], potar[:, 1])
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    if methode == "polynomial":
        func = si.KroghInterpolator(potar[:, 0], potar[:, 1])
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    return x_val, potx


def _schroedinger_equation_solver(mass, interpdata, potx, eigmin, eigmax):
    """Solves the schroedingerequation from a given set of potentials.

    Args:
        datalist = list created by the communicator modul function
        input_reader.

    Returns:
        Eigenvalues and their normalized wavefunctions as ndarray.
    """
    x_min, x_max, npoint = interpdata
    step = (x_max - x_min) / npoint
    coeff = 1 / (mass * step**2)
    # Calculation of the eigenvalues and wavefunctions
    maindia = np.array([])
    for val in enumerate(potx):
        elem = coeff + val[1]
        maindia = np.append(maindia, elem)
    offdia = np.full(len(potx) - 1, -coeff / 2)
    eigval, wavef = sl.eigh_tridiagonal(maindia, offdia, False,
                                        'i', (eigmin - 1, eigmax - 1),
                                        True, 0.0, 'stebz')
    # Normalization of the wavefunctions
    normwf = np.empty((len(wavef), len(eigval)))
    for i in range(len(wavef[0])):
        norm = step * sum(wavef[:, i]**2)
        normwf[:, i] = wavef[:, i] / norm**0.5
    return eigval, normwf

def expected_values(interpdata, normwf, x_val):
    x_min, x_max, npoint = interpdata
    step = (x_max - x_min) / npoint
    expxlist = []
    sigma = []
    for i, val in enumerate(normwf[0, :]):
        expx = step * sum(normwf[:, i] * x_val * normwf[:, i])
        expx2 = step * sum(normwf[:, i] * x_val**2 * normwf[:, i])
        expxlist.append(expx)
        sigma.append((expx2 - expx**2)**0.5)
    return expxlist, sigma

