"""Contains tests for solvers module"""
import pytest
import numpy as np
import solvers as sol

FILES = ["test_infpot.txt", "test_pot.txt", "test_harmonic.txt",
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
    if file == "test_pot.txt":
        for n in range(eigmin, eigmax + 1):
            if n % 2 == 0:
                eigval = (4 * np.pi**2 * h**2) / (8 * m * L**2) * (2 * n - 1)**2
            if n % 2 == 1:
                eigval = (4 * np.pi**2 * h**2) / (8 * m * L**2) * (2 * n)**2
            eigvallist.append(eigval)
    if file == "test_harmonic.txt":
        for n in range(eigmin - 1, eigmax):
            eigval = 1 / 2 * h * w * (n + 1 / 2)
            eigvallist.append(eigval)
    if file == "test_dualpot_lin.txt":
        eigvallist = np.ones((1, eigmax - eigmin + 1))
    if file == "test_dualpot_scpline.txt":
        eigvallist = np.ones((1, eigmax - eigmin + 1))
    if file == "test_asympot.txt":
        eigvallist = np.ones((1, eigmax - eigmin + 1))
    eigvalarray = np.array(eigvallist)
    sol.run(filedir)
    testeigarray = np.loadtxt("output/energies.dat")
    assert np.all(np.abs(eigvalarray - testeigarray) < ERROR)
