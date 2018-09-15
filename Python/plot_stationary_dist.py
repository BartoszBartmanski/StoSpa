#!/usr/bin/python3

"""Plots the stationary distribution.

Usage:
    plot_stationary_dist.py -h | --help
    plot_stationary_dist.py <file_name>... [options]

Options:
    -h --help                 Show this screen
    --x_lim=<x_lim>           Change the x-axis limits.
    --sep                     Whether to use a separate figure for each file.
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

    if not mean:
        data = np.loadtxt(file_name)
        mean = np.argmax(data)

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


def plot(d, name, lim, sep_figs=False):
    if sep_figs:
        plt.figure(name)
        plt.bar(d[0], d[1][0] / np.sum(d[1][0]), label="Stochastic")
    else:
        plt.bar(d[0], d[1][0] / np.sum(d[1][0]), label=name)

    for i in range(1, len(d[1])):
        plt.plot(d[0], d[1][i]/np.sum(d[1][i]))

    plt.xlim(lim[0]-0.5, lim[1]+0.5)
    plt.gca().set_xticks(np.arange(lim[0], lim[1]+1, lim[2]))
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

    redraw = True
    for name in args["<file_name>"]:
        data = get_data(name)
        limits = get_range(len(data[0]), args["--x_lim"])
        plot(data, name, limits, args["--sep"])

        if redraw:
            d_a = get_analytic(get_mean(name), limits)
            plt.plot(d_a[0], d_a[1], label="Analytical", color="black")
            if args["--sep"]:
                redraw = True
            else:
                redraw = False

    plt.show()
