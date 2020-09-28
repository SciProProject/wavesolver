"""Module used to solve onedimensional timeindependent Schroedinger equations
   using a set of data from a given file.
"""
import scipy.linalg as sl
import scipy.interpolate as si
import numpy as np


def run(filedir):
    """Calculates the eigenvalues, wavefunctions, expected values and
       the interpolated potentials then returns the result in respectiv files.

       Args:
           filedir: directory of the inputfile

       Returns:
           four files into the output folder
    """
    mass, interpdata, methode,\
        potar, eigmin, eigmax = input_reader(filedir)[:]
    potential = potential_interpolate(interpdata, methode, potar)
    eigval, normwf = schroedinger_equation_solver(mass, interpdata,
                                                  potential, eigmin, eigmax)
    expvalues = expected_values(interpdata, normwf[:, 1:],
                                potential[:, 0], eigval)
    np.savetxt("output/potential.dat", potential)
    np.savetxt("output/wavefuncs.dat", normwf)
    np.savetxt("output/energies.dat", np.array(eigval))
    np.savetxt("output/expvalues.dat", expvalues)


def input_reader(filedir):
    """Collects data from the file the user iputs whitch is used to solve
       the wavefunction in the solvers module.

       Args:
           filedir: directory of the inputfile

       Returns:
           mass, interpdata, methode, potar, eigmin, eigmax
           used in solving the schroedinger equation.
    """
    try:
        data = np.genfromtxt(filedir, dtype=str, delimiter="/n")
    except:
        raise FileNotFoundError("""Check the the given file directory
{}
and make sure it leads to the right file""".format(filedir))
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

    return mass, interpdata, methode, potar, eigmin, eigmax


def potential_interpolate(interpdata, methode, potar):
    """interpolates the potentail for a new range of x values from a
       given set of potentails.

    Args:
        interpdata -- list with range and number of points to interpolate
        methode -- methode of interpolation (linear, cspline, polynomial)
        potar -- 2darray of x and potential values

    Returns:
        the interpolated potentials and their x values as ndarray.
    """
    if methode == "linear":
        func = si.interp1d(potar[:, 0], potar[:, 1])
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    elif methode == "cspline":
        func = si.CubicSpline(potar[:, 0], potar[:, 1], bc_type='natural',
                              extrapolate=True)
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    elif methode == "polynomial":
        func = si.KroghInterpolator(potar[:, 0], potar[:, 1])
        x_val = np.linspace(interpdata[0], interpdata[1], interpdata[2])
        potx = func(x_val)
    else:
        raise TypeError("""Your method {} was not understood,
choose from [linear, cspline, polynomial]""".format(methode))

    potential = np.empty((len(x_val), 2))
    potential[:, 0] = x_val
    potential[:, 1] = potx
    return potential


def schroedinger_equation_solver(mass, interpdata, potential, eigmin, eigmax):
    """Solves the schroedingerequation from a given set of potentials.

    Args:
        mass -- Mass of the partical
        interpdata -- list with range and number of points to interpolate
        potential -- 2darray of x and potential values after interpolation
        eigmin -- minimal amount of calculated Eigenvalues
        eigmax -- maximum amount of calculatzed Eigenvalues

    Returns:
        Eigenvalues and their normalized wavefunctions as ndarray.
    """
    x_min, x_max, npoint = interpdata
    step = (x_max - x_min) / npoint
    potx = potential[:, 1]
    coeff = 1 / (mass * step**2)
    # Calculation of the eigenvalues and wavefunctions
    maindia = np.array([])
    for val in enumerate(potx):
        elem = coeff + val[1]
        maindia = np.append(maindia, elem)
    offdia = np.full(len(potx) - 1, -coeff / 2)
    eigval, wavef = sl.eigh_tridiagonal(maindia, offdia, False,
                                        'i', (eigmin - 1, eigmax - 1),
                                        True, 0.0, 'auto')
    # Normalization of the wavefunctions
    normw = np.empty((len(wavef), len(eigval)))
    for i in range(len(wavef[0])):
        norm = step * sum(wavef[:, i]**2)
        normw[:, i] = wavef[:, i] / norm**0.5
    normwf = np.empty((len(potential[:, 0]), len(eigval) + 1))
    normwf[:, 0] = potential[:, 0]
    for i, val in enumerate(np.transpose(normw)):
        normwf[:, i + 1] = val
    return eigval, normwf


def expected_values(interpdata, normwf, x_val, eigval):
    """Calculates the expected values for the given
       wavefunctions, eigenenergies and potentials.

    Args:
        interpdata -- list with range and number of points to interpolate
        normwf -- ndarray of wavefunctions
        x_val -- list of x values
        eigval -- list of eigenenergies

    Returns:
        expected x values and its standrt diviation for the respectiv
        wavefunctions.
    """
    x_min, x_max, npoint = interpdata
    step = (x_max - x_min) / npoint
    expxlist = []
    sigma = []
    for val in enumerate(normwf[0, :]):
        pos = val[0]
        expx = step * sum(normwf[:, pos] * x_val * normwf[:, pos])
        expx2 = step * sum(normwf[:, pos] * x_val**2 * normwf[:, pos])
        expxlist.append(expx)
        sigma.append((expx2 - expx**2)**0.5)
    expvalues = np.empty((len(eigval), 2))
    expvalues[:, 0] = expxlist
    expvalues[:, 1] = sigma
    return expvalues
