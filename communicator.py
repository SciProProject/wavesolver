#!/usr/bin/env python3
"""Retrives data for the solvers module from the files given by the user."""
import solvers
import visualization


def main():
    """Solves the Schrödingerequation from a file containing data and outputs
       files with the eigenenergies, wavefunctions, interpolated potentials
       and expected values also visualizing the solution.
    """
    filedir = input('Please input the directory of the file with the data. ')

    path = input("""Please inpute the directory where you
want to save the files.""")
    solvers.run(filedir, path)

    visualization.run(path)


if __name__ == '__main__':
    main()
