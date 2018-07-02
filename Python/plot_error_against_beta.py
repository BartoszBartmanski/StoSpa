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


def get_data(file_name):
    num_points = None

    f = open(file_name, "r")
    for line in f:
        if "#" not in line:
            break
        if "num_points" in line and "command" not in line.lower():
            num_points = int(line.split(" ")[-1])
            break
    f.close()

    data = np.loadtxt(file_name)
    if num_points is None:
        num_points = np.sqrt(data.shape[1])
    data = data.reshape([3, num_points, num_points]).transpose([0, 2, 1])

    return data


def get_min(data):
    return np.unravel_index(np.argmin(data), data.shape)


def get_index(mesh_x, mesh_y, value):
    """Returns the index of point in the mesh closest to the parameter value"""

    idx_0 = (np.abs(mesh_x - value[0])).argmin()
    idx_1 = (np.abs(mesh_y - value[1])).argmin()
    idx = np.unravel_index([idx_0, idx_1], mesh_x.shape)
    idx = (idx[0][1], idx[1][0])

    return idx


def plot(arguments):
    for name in arguments["<file_name>"]:

        d = get_data(name)

        plt.figure(name)

        if arguments["--get_min"]:
            s = r"Coordinates of the minimum for {} = ({}, {})"
            min_index = get_min(d[2])
            print(s.format(name, d[0][min_index], d[1][min_index]))

        if arguments["--get_error"]:
            s = r"Value of error at ({}, {}) = {}"
            index = get_index(d[0], d[1], arguments["--ger_error"])
            print(s.format(d[0][index], d[1][index], d[2][index]))

        if arguments["--1d"]:
            plt.plot(d[1][d[0] == d[1]], d[2][d[0] == d[1]], color="black")
            plt.xlabel(r"$\beta$")
            plt.ylabel("e")
        else:
            cs = plt.contour(d[0], d[1], d[2], np.linspace(np.amin(d[2]), np.amax(d[2]), 20))
            plt.clabel(cs, inline=1, fontsize=10)
            plt.ylabel(r"$\beta_y$")
            plt.xlabel(r"$\beta_x$")
            plt.gca().set_aspect('equal')

            min_index = get_min(d[2])
            plt.scatter(d[0][min_index], d[1][min_index], label="Min.")
            plt.legend()


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)
    if args["--get_error"]:
        args["--get_error"] = [float(x) for x in args["--get_error"].split(",")]
        if len(args["--get_error"]) == 1:
            args["--get_error"].append(args["--get_error"][0])

    plot(args)

    plt.show()
