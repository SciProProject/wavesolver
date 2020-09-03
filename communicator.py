#!/usr/bin/env python3
"""Retrives data for the solvers module from the files given by the user."""
import solvers
import visualization

def main():
    filedir = input('Please input the directory of the file with the data. ')

    solvers.run(filedir)

    visualization.run()

if __name__ == '__main__':
    main()