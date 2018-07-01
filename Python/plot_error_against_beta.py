#!/usr/bin/python3

"""Plots error against beta.

Usage:
    plot_error_against_beta.py -h | --help
    plot_error_against_beta.py <file_name>... [options]

Options:
    -h --help                       Show this screen
    --get_error=<val>		        Get value of error for given values of beta.
    --get_min                       Returns the location of the minimum point.
    --1d                            Plots the error along beta_x = beta_y
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)
    if args["--get_error"]:
        args["--get_error"] = [float(x) for x in args["--get_error"].split(",")]
        if len(args["--get_error"]) == 1:
            args["--get_error"].append(args["--get_error"][0])

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
        min_index = np.unravel_index(np.argmin(d[2]), d[2].shape)

        plt.figure(name)

        if args["--get_min"]:
            s = r"Cooridinates of the minimum for {} = ({}, {})"
            print(s.format(name, d[0][min_index], d[1][min_index]))

        if args["--get_error"]:
            s = r"Value of error at ({}, {}) = {}"
            idx_0 = (np.abs(d[0] - args["--get_error"][0])).argmin()
            idx_1 = (np.abs(d[1] - args["--get_error"][1])).argmin()
            idx = np.unravel_index([idx_0, idx_1], [num_points, num_points])
            idx = (idx[0][1], idx[1][0])
            error = d[2][idx]
            print(s.format(d[0][idx], d[1][idx], error))

        if args["--1d"]:
            plt.plot(d[1][d[0] == d[1]], d[2][d[0] == d[1]], color="black")
            plt.xlabel(r"$\beta$")
            plt.ylabel("e")
        else:
            cs = plt.contour(d[0], d[1], d[2], np.linspace(np.amin(d[2]), np.amax(d[2]), 20))
            plt.clabel(cs, inline=1, fontsize=10)
            plt.ylabel(r"$\beta_y$")
            plt.xlabel(r"$\beta_x$")
            plt.gca().set_aspect('equal')

            plt.scatter(d[0][min_index], d[1][min_index], label="Min.")
            plt.legend()

    plt.show()
