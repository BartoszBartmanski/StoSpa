#!/usr/bin/python3

"""Plots error against kappa.

Usage:
    plot_error_against_kappa.py -h | --help
    plot_error_against_kappa.py <file_name>... [options]

Options:
    -h --help                       Show this screen
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def plot(arguments):
    labels = ["FEM", "FDM", "FVM", "FET"]

    for name in arguments["<file_name>"]:
        plt.figure(name)
        d = np.loadtxt(name)
        assert d.shape[0] == 5
        for i in range(4):
            plt.plot(d[0], d[i+1], label=labels[i])

        plt.legend()
        plt.ylabel(r"Error")
        plt.xlabel(r"$\kappa$")
        plt.tight_layout()


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    plot(args)

    plt.show()
