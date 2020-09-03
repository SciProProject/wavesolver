"""Contains tests for solvers module"""
import pytest
import numpy as np
import solvers as sol
import visualization as vis

FILES = ["test_infpot.txt", "test_harmonic.txt", "test_pot.txt",
         "test_dualpot_lin.txt", "test_dualpot_scpline.txt",
         "test_asympot.txt"]
ERROR = 1e-2


@pytest.mark.parametrize("file", FILES)
def test_eigval(file):
    filedir = "testfiles/" + file
    m = sol._input_reader(filedir)[0]
    x_min, x_max = sol._input_reader(filedir)[1][:2]
    L = x_max - x_min
    eigmin, eigmax = sol._input_reader(filedir)[4:6]
    eigvallist = []
    if file == "test_infpot.txt":
        for n in range(eigmin, eigmax + 1):
            eigval = (4 * np.pi**2) / (8 * m * L**2) * n**2
            eigvallist.append(eigval)
    elif file == "test_harmonic.txt":
        for n in range(eigmin - 1, eigmax):
            eigval = 1 / 2 * (n + 1 / 2)
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
    interpdata = sol._input_reader(filedir)[1]
    x_min, x_max = interpdata[:2]
    L = x_max - x_min
    methode, potar = sol._input_reader(filedir)[2:4]
    eigmin, eigmax = sol._input_reader(filedir)[4:6]
    x_vals = sol._potential_interpolate(interpdata, methode, potar)[0]
    wavefuncsmat = np.empty((len(x_vals), eigmax - eigmin + 1))
    if file == "test_infpot.txt":
        for n in range(eigmin, eigmax + 1):
            wavefuncs = []
            if n % 2 == 0:
                for val in enumerate(x_vals):
                    wavefunc = np.sqrt(2 / L) * np.sin(n * np.pi / L * val[1])
                    wavefuncs.append(wavefunc)
            else:
                for val in enumerate(x_vals):
                    wavefunc = np.sqrt(2 / L) * np.cos(n * np.pi / L * val[1])
                    wavefuncs.append(wavefunc)
            wavefuncsmat[:, n - 1] = wavefunc
    elif file == "test_harmonic.txt":
        wavefuncsmat = np.loadtxt("__unittestfiles__/test_harmonic.txt")[:, 1:]
    elif file == "test_pot.txt":
        wavefuncsmat = np.loadtxt("__unittestfiles__/test_pot.txt")[:, 1:]
    elif file == "test_dualpot_lin.txt":
        wavefuncsmat = np.loadtxt("__unittestfiles__/test_dualpot_lin.txt")[:, 1:]
    elif file == "test_dualpot_scpline.txt":
        wavefuncsmat = np.loadtxt("__unittestfiles__/test_dualpot_scpline.txt")[:, 1:]
    elif file == "test_asympot.txt":
        wavefuncsmat = np.loadtxt("__unittestfiles__/test_asympot.txt")[:, 1:]
    else:
        wavefuncsmat = np.ones((len(x_vals), eigmax - eigmin + 1))
    sol.run(filedir)
    testwavefuncs = np.loadtxt("output/wavefuncs.dat")[:, 1:]
    assert np.all(np.abs(np.abs(wavefuncsmat)**2 - np.abs(testwavefuncs)**2) < ERROR)


def _creat_testdata(file):
    filedir = "testfiles/" + file
    sol.run(filedir)
    vis.visualizer()
    check = input("Should this data be used as testdata (input Y or N).")
    if check == "Y":
        newdirwf = "__unittestfiles__/" + file.split(".")[0] + "_wf.dat"
        newdirenergy = "__unittestfiles__/" + file.split(".")[0] + "_energy.dat"
        testdatawf = np.loadtxt("output/wavefuncs.dat")
        np.savetxt(newdirwf, testdatawf)
        testdataenergy = np.loadtxt("output/energies.dat")
        np.savetxt(newdirenergy, testdataenergy)
    if check == "N":
        return None
