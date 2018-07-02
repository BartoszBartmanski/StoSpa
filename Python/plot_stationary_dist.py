#!/usr/bin/python3

"""Plots the stationary distribution.

Usage:
    plot_stationary_dist.py -h | --help
    plot_stationary_dist.py <file_name>... [options]

Options:
    -h --help                           Show this screen
    --x_lim=<x_lim>                     Change the x-axis limits.
    --separate                          Whether to use a separate figure for each file.
"""

from docopt import docopt
import scipy.misc as misc
import numpy as np
import matplotlib.pyplot as plt


def get_mean(file_name):
    mean = None

    f = open(file_name, "r")
    for line in f:
        if "#" not in line:
            break
        if "means" in line:
            mean = float(line.split(" ")[-1])
            break
    f.close()

    return mean


def get_range(data_len, x_lim):
    x_lim = x_lim + [0, data_len][len(x_lim):]
    tick_step = max((x_lim[1] - x_lim[0]) // 10, 1)
    x_lim += [tick_step]
    return x_lim


def get_data(file_name):
    y = np.atleast_2d(np.loadtxt(file_name))
    x = np.arange(y.shape[1])

    return x, y


def get_analytic(mean, x_lim):
    x = np.arange(x_lim[0], x_lim[1])
    y = 1.0 / (misc.factorial(x)) * mean ** x * np.exp(-mean)

    return x, y


def plot(arguments):

    analytic = True
    for name in arguments["<file_name>"]:
        if arguments["--separate"]:
            plt.figure(name)
            analytic = True

        d = get_data(name)
        limits = get_range(len(d[0]), arguments["--x_lim"])

        if analytic:
            m = get_mean(name)
            d_a = get_analytic(m, limits)
            plt.plot(d_a[0], d_a[1], label="Analytical", color="black")
            analytic = False

        plt.bar(d[0], d[1][0] / np.sum(d[1][0]), label="Stochastic")
        for i in range(1, len(d[1])):
            plt.plot(d[0], d[1][i]/np.sum(d[1][i]))

        plt.xlim(limits[0]-0.5, limits[1]+0.5)
        plt.gca().set_xticks(np.arange(limits[0], limits[1]+1, limits[2]))
        plt.title("Stationary dist.")
        plt.xlabel("# of molecules")
        plt.ylabel("Probability")
        plt.legend()


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)
    if args["--x_lim"]:
        args["--x_lim"] = [int(x) for x in args["--x_lim"].split(",")]
    else:
        args["--x_lim"] = []

    plot(args)

    plt.show()
