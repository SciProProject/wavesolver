Wavefuncsolver modul
Authors: Mike Pfannenstiel, Fabian Piwowarczyk

Functionality:
1. The user calls the communicater.py in the shell to start the programm.
2. The programm asks for the directory of the input file.
    - The input file (example.inp) should be formated like the example below.
3. The programm asks for the directory where the output files should be saved.
The output files are:
    - energies.dat:
        Contains the eigenvalues for the given range.
    - expvalues.dat
        Contains the expected values (first column)
        and standard diviation (second column).
    - potential.dat
        Contains the xy declarations for the interpolated potential
    - wavefuncs.dat
        Contains the x-values of the wavefunctions (first column)
        and every column after the first contains the y-values for the
        corresponding eigenvalue.
4. The programm starts visualization.
5. Stonks/Profit

example.inp:
##############################################################################
2.0             # mass
-2.0 2.0 1999   # xMin xMax nPoint
1 5             # first and last eigenvalue to print
linear          # interpolation type
2               # nr. of interpolation points and xy declarations
-2.0 0.0
2.0 0.0
##############################################################################