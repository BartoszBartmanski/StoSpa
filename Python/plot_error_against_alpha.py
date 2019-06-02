#!/usr/bin/env python3

"""Plots error against alpha.

Usage:
    plot_error_against_alpha.py -h | --help
    plot_error_against_alpha.py <filename>... [options]

Options:
    -h --help                       Show this screen
"""

import matplotlib.pyplot as plt
import numpy as np

from docopt import docopt

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    for name in args["<filename>"]:
        data = np.loadtxt(name)
        plt.plot(data[0], data[1])

    plt.ylabel(r"Error")
    plt.xlabel(r"$\alpha$")
    plt.show()
