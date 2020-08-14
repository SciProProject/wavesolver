"""Contains tests for solvers module"""
import pytest
import numpy as np
import solvers as sol

FILES = ["test_infpot.txt", "test_harmonic.txt", "test_pot.txt",
         "test_dualpot_lin.txt", "test_dualpot_scpline.txt",
         "test_asympot.txt"]
ERROR = 1e-2


@pytest.mark.parametrize("file", FILES)
def test_eigval(file):
    filedir = "testfiles/" + file
    h = 1
    m = sol._input_reader(filedir)[0]
    x_min, x_max = sol._input_reader(filedir)[1][:2]
    L = x_max - x_min
    w = 1
    eigmin, eigmax = sol._input_reader(filedir)[5:7]
    eigvallist = []
    if file == "test_infpot.txt":
        for n in range(eigmin, eigmax + 1):
            eigval = (4 * np.pi**2 * h**2) / (8 * m * L**2) * n**2
            eigvallist.append(eigval)
    elif file == "test_harmonic.txt":
        for n in range(eigmin - 1, eigmax):
            eigval = 1 / 2 * h * w * (n + 1 / 2)
            eigvallist.append(eigval)
#    elif file == "test_pot.txt":
#        eigvallist = np.ones((1, eigmax - eigmin + 1))
#    elif file == "test_dualpot_lin.txt":
#        eigvallist = np.ones((1, eigmax - eigmin + 1))
#    elif file == "test_dualpot_scpline.txt":
#        eigvallist = np.ones((1, eigmax - eigmin + 1))
#    elif file == "test_asympot.txt":
#        eigvallist = np.ones((1, eigmax - eigmin + 1))
    else:
        eigvallist = np.ones((1, eigmax - eigmin + 1))
    eigvalarray = np.array(eigvallist)
    sol.run(filedir)
    testeigarray = np.loadtxt("output/energies.dat")
    assert np.all(np.abs(eigvalarray - testeigarray) < ERROR)


@pytest.mark.parametrize("file", FILES)
def test_wavefunc(file):
    filedir = "testfiles/" + file
    x_min, x_max = sol._input_reader(filedir)[1][:2]
    L = x_max - x_min
    eigmin, eigmax = sol._input_reader(filedir)[5:7]
    x_vals = sol._potential_interpolate(filedir)[0]
    wavefuncsmat = np.empty((len(x_vals), eigmax - eigmin + 1))
    if file == "test_infpot.txt":
        for n in range(eigmin, eigmax + 1):
            wavefuncs = []
            for val in enumerate(x_vals):
                wavefunc = np.sqrt(2 / L) * np.sin(n * np.pi / L * val[0])
                wavefuncs.append(wavefunc)
            wavefuncsmat[:, n - 1] = wavefuncs
#    elif file == "test_harmonic.txt":
#
#    elif file == "test_pot.txt":
#
#    elif file == "test_dualpot_lin.txt":
#
#    elif file == "test_dualpot_scpline.txt":
#
#    elif file == "test_asympot.txt":
    else:
        wavefuncsmat = np.ones((len(x_vals), eigmax - eigmin + 1))
    sol.run(filedir)
    testwavefuncs = np.loadtxt("output/wavefuncs.dat")[:, 1:]
    assert np.all(np.abs(wavefuncsmat - testwavefuncs) < ERROR)
