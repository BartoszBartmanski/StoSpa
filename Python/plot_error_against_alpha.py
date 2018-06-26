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

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    for name in args["<file_name>"]:
        d = np.loadtxt(name)
        assert d.shape[0] == 2
        plt.plot(d[0], d[1])

    plt.ylabel(r"Error")
    plt.xlabel(r"$\alpha$")
    plt.tight_layout()
    plt.show()
