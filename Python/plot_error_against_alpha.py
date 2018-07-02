#!/usr/bin/python3

"""Plots error against alpha.

Usage:
    plot_error_against_alpha.py -h | --help
    plot_error_against_alpha.py <file_name>... [options]

Options:
    -h --help                       Show this screen
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def plot(arguments):

    for name in arguments["<file_name>"]:
        d = np.loadtxt(name)
        assert d.shape[0] == 2
        plt.plot(d[0], d[1])

    plt.ylabel(r"Error")
    plt.xlabel(r"$\alpha$")
    plt.tight_layout()


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    plot(args)

    plt.show()
