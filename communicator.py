"""Retrives data for the solvers module from the files given by the user."""


def input_reader(filedir):
    """Collects data from the file the user iputs whitch is used to solve
       the wavefunction in the solvers module.

       Args:
           filedir: directory of the file

       Returns:
           List of data used in solving the schroedinger equation.
    """
    with open(filedir, "r") as fp:
        initlist = fp.readlines()
        datalist = []
        for tup in enumerate(initlist):
            part1 = tup.partition('#')
            part2 = part1[0].partition('\n')
            datalist.append(part2[0])
    return datalist
