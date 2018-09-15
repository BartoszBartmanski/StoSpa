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


def plot(d, name):
    assert d.shape[0] == 5

    labels = ["FEM", "FDM", "FVM", "FET"]

    plt.figure(name)

    for i in range(4):
        plt.plot(d[0], d[i+1], label=labels[i])

    plt.legend()
    plt.ylabel(r"Error")
    plt.xlabel(r"$\kappa$")


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    for name in arguments["<file_name>"]:
        data = np.loadtxt(name)
        plot(data, name)

    plt.show()
