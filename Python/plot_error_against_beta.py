#!/usr/bin/python3

"""Plots error against beta.

Usage:
    plot_error_against_beta.py -h | --help
    plot_error_against_beta.py <filename>... [options]

Options:
    -h --help                       Show this screen
    --get_error=<val>		        Get value of error for given values of beta.
    --get_min                       Returns the location of the minimum point.
    --1d                            Plots the error along beta_x = beta_y
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def get_data(filename):
    num_points = None

    f = open(filename, "r")
    for line in f:
        if "#" not in line:
            break
        if "num_points" in line and "command" not in line.lower():
            num_points = int(line.split(" ")[-1])
            break
    f.close()

    d = np.loadtxt(filename)
    if num_points is None:
        num_points = np.sqrt(d.shape[1])
    d = d.reshape([3, num_points, num_points]).transpose([0, 2, 1])

    return d


def get_min(d):
    return np.unravel_index(np.argmin(d), d.shape)


def get_index(mesh_x, mesh_y, value):
    """Returns the index of point in the mesh closest to the parameter value"""

    idx_0 = (np.abs(mesh_x - value[0])).argmin()
    idx_1 = (np.abs(mesh_y - value[1])).argmin()
    idx = np.unravel_index([idx_0, idx_1], mesh_x.shape)
    idx = (idx[0][1], idx[1][0])

    return idx


def plot(d, name, reduce=True):

    plt.figure(name)

    if reduce:
        plt.plot(d[1][d[0] == d[1]], d[2][d[0] == d[1]], color="black")
        plt.xlabel(r"$\beta$")
        plt.ylabel("e")
    else:
        num_cont = np.linspace(np.amin(d[2]), np.amax(d[2]), 20)
        cs = plt.contour(d[0], d[1], d[2], num_cont)
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
        e = args["--get_error"].split(",")
        args["--get_error"] = [float(x) for x in e]
        if len(args["--get_error"]) == 1:
            args["--get_error"].append(args["--get_error"][0])

    for name in args["<filename>"]:
        data = get_data(name)

        if args["--get_min"]:
            min_index = get_min(data[2])
            s = r"Coordinates of the minimum for {} = ({}, {})"
            print(s.format(name, data[0][min_index], data[1][min_index]))

        if args["--get_error"]:
            index = get_index(data[0], data[1], arguments["--ger_error"])
            s = r"Value of error at ({}, {}) = {}"
            print(s.format(data[0][index], data[1][index], data[2][index]))

        plot(data, name, args["--1d"])

    plt.show()
