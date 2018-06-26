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

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    colors = ["r", "y", "c", "m", "g"]

    analytic = True
    for name in args["<file_name>"]:
        if args["--separate"]:
            analytic = True
            plt.figure(name)

        y = np.atleast_2d(np.loadtxt(name))
        x = np.arange(y.shape[1])

        # Get the simulation mean number of molecules for species A from the first file
        mean_a = None
        f = open(name, "r")
        for line in f:
            if "#" not in line:
                break
            if "means" in line:
                mean_a = float(line.split(" ")[-1])
                break
        f.close()

        plt.bar(x, y[0]/np.sum(y[0]), label="Stochastic")
        limits = [0, len(x)]
        if args["--x_lim"]:
            args["--x_lim"] = [int(x) for x in args["--x_lim"].split(",")]
            limits = args["--x_lim"] + limits[len(args["--x_lim"]):]
        tick_step = max((limits[1] - limits[0])//10, 1)

        if mean_a and analytic:
            x_values = np.arange(limits[0], limits[1])
            y_values = 1.0 / (misc.factorial(x_values)) * mean_a ** x_values * np.exp(-mean_a)
            plt.plot(x_values, y_values, label="Analytical", color="black")
        analytic = False

        for i in range(1, len(y)):
            plt.plot(x, y[i]/np.sum(y[i]), color=colors[i % len(colors)])

        plt.xlim(limits[0]-0.5, limits[1]+0.5)
        plt.gca().set_xticks(np.arange(limits[0], limits[1]+1, tick_step))
        plt.title("Stationary dist.")
        plt.xlabel("# of molecules")
        plt.ylabel("Probability")
        plt.legend()
    plt.show()
