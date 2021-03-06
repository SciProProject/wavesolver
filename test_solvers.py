#!/usr/bin/env python3
"""Contains tests for solvers module"""
import pytest
import numpy as np
import matplotlib.pyplot as plt
import solvers as sol
import visualization as vis


FILES = ["test_infpot.txt", "test_harmonic.txt", "test_pot.txt",
         "test_dualpot_lin.txt", "test_dualpot_cspline.txt",
         "test_asympot.txt"]
ERROR = 5e-3


@pytest.mark.parametrize("file", FILES)
def test_eigval(file):
    """test if the eigenenergies created by the files in the testfiles folder
       are realistic.

       Args:
           file    filename from the testfiles folder
    """
    filedir = "__testfiles__/" + file
    mass = sol.input_reader(filedir)[0]
    x_min, x_max = sol.input_reader(filedir)[1][:2]
    length = x_max - x_min
    eigmin, eigmax = sol.input_reader(filedir)[4:6]
    # getting the eigenvalues for the specific problem from data or equations
    eigvallist = []
    if file == "test_infpot.txt":
        for nn in range(eigmin, eigmax + 1):
            eigval = (4 * np.pi**2) / (8 * mass * length**2) * nn**2
            eigvallist.append(eigval)
    elif file == "test_harmonic.txt":
        for nn in range(eigmin - 1, eigmax):
            eigval = 1 / 2 * (nn + 1 / 2)
            eigvallist.append(eigval)
    elif file == "test_pot.txt":
        eigvallist = np.loadtxt("__unittestfiles__/test_pot_energy.dat")
    elif file == "test_dualpot_lin.txt":
        eigvallist = np.loadtxt\
            ("__unittestfiles__/test_dualpot_lin_energy.dat")
    elif file == "test_dualpot_cspline.txt":
        eigvallist = np.loadtxt\
            ("__unittestfiles__/test_dualpot_cspline_energy.dat")
    elif file == "test_asympot.txt":
        eigvallist = np.loadtxt("__unittestfiles__/test_asympot_energy.dat")
    else:
        eigvallist = np.ones((1, eigmax - eigmin + 1))
    eigvalarray = np.array(eigvallist)
    sol.run(filedir, "__output__")
    testeigarray = np.loadtxt("__output__/energies.dat")
    assert np.all(np.abs(eigvalarray - testeigarray) < ERROR)


@pytest.mark.parametrize("file", FILES)
def test_wavefunc(file):
    """test if the wavefunctions created by the files in the testfiles folder
       are realistic.

       Args:
           file    filename from the testfiles folder
    """
    filedir = "__testfiles__/" + file
    interpdata = sol.input_reader(filedir)[1]
    x_min, x_max = interpdata[:2]
    length = x_max - x_min
    methode, potar = sol.input_reader(filedir)[2:4]
    eigmin, eigmax = sol.input_reader(filedir)[4:6]
    x_vals = sol.potential_interpolate(interpdata, methode, potar)[:, 0]
    # getting the wavefunctions for the specific problem from data or equations
    wavefuncsmat = np.empty((len(x_vals), eigmax - eigmin + 1))
    if file == "test_infpot.txt":
        for nn in range(eigmin, eigmax + 1):
            wavefuncs = []
            if nn % 2 == 0:
                for val in enumerate(x_vals):
                    wavefunc = np.sqrt(2 / length) * np.sin(nn * np.pi /
                                                            length * val[1])
                    wavefuncs.append(wavefunc)
            else:
                for val in enumerate(x_vals):
                    wavefunc = np.sqrt(2 / length) * np.cos(nn * np.pi /
                                                            length * val[1])
                    wavefuncs.append(wavefunc)
            wavefuncsmat[:, nn - 1] = wavefuncs
    elif file == "test_harmonic.txt":
        wavefuncsmat = np.loadtxt\
            ("__unittestfiles__/test_harmonic_wf.dat")[:, 1:]
    elif file == "test_pot.txt":
        wavefuncsmat = np.loadtxt\
            ("__unittestfiles__/test_pot_wf.dat")[:, 1:]
    elif file == "test_dualpot_lin.txt":
        wavefuncsmat = np.loadtxt\
            ("__unittestfiles__/test_dualpot_lin_wf.dat")[:, 1:]
    elif file == "test_dualpot_cspline.txt":
        wavefuncsmat = np.loadtxt\
            ("__unittestfiles__/test_dualpot_cspline_wf.dat")[:, 1:]
    elif file == "test_asympot.txt":
        wavefuncsmat = np.loadtxt\
            ("__unittestfiles__/test_asympot_wf.dat")[:, 1:]
    else:
        wavefuncsmat = np.ones((len(x_vals), eigmax - eigmin + 1))
    # checking if the norm of the wavefunctions is equal to testdata
    sol.run(filedir, "__output__")
    testwavefuncs = np.loadtxt("__output__/wavefuncs.dat")[:, 1:]
    assert np.all(np.abs(np.abs(wavefuncsmat)**2 - np.abs(testwavefuncs)**2)
                  < ERROR)


def _creat_testdata(file):
    """creates energies.dat and wavefuncs.dat for the files in the testfiles
       folder which are used in unit testing.

       Args:
           file    filename from the testfiles folder
    """
    # Creating testdata and visualizing
    filedir = "__testfiles__/" + file
    sol.run(filedir, "__output__")
    vis.run("__output__")
    plt.pause(5)
    # Saving testdata if approved
    check = input("Should this data be used as testdata [y/n]? ")
    if check == "y":
        newdirwf = "__unittestfiles__/" + file.split(".")[0] + "_wf.dat"
        newdirenergy = "__unittestfiles__/" + file.split(".")[0]\
            + "_energy.dat"
        testdatawf = np.loadtxt("__output__/wavefuncs.dat")
        np.savetxt(newdirwf, testdatawf)
        testdataenergy = np.loadtxt("__output__/energies.dat")
        np.savetxt(newdirenergy, testdataenergy)
        plt.close('all')
    if check == "n":
        plt.close('all')
