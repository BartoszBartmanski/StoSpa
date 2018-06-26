#!/usr/bin/python3

"""Plots error against beta.

Usage:
    plot_error_against_beta.py -h | --help
    plot_error_against_beta.py <file_name>... [options]

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
        # Get the size of square grid of errors

        num_points = None

        f = open(name, "r")
        for line in f:
            if "#" not in line:
                break
            if "num_points" in line and "command" not in line.lower():
                num_points = int(line.split(" ")[-1])
                break
        f.close()

        d = np.loadtxt(name)
        if num_points is None:
            num_points = np.sqrt(d.shape[1])
        d = d.reshape([3, num_points, num_points]).transpose([0, 2, 1])
        index = np.unravel_index(np.argmin(d[2]), d[2].shape)

        plt.figure(name)
        cs = plt.contour(d[0], d[1], d[2], np.linspace(np.amin(d[2]), np.amax(d[2]), 20))
        plt.clabel(cs, inline=1, fontsize=10)
        plt.ylabel(r"$\beta_y$")
        plt.xlabel(r"$\beta_x$")
        plt.colorbar(cs)
        plt.gca().set_aspect('equal')
        plt.scatter(d[0][index], d[1][index])

    plt.show()
