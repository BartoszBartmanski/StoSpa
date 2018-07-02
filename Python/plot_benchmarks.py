#!/usr/bin/python3

"""Plots the benchmarking data.

Usage:
    plot_benchmarks.py -h | --help
    plot_benchmarks.py <file_name>... [options]

Options:
    -h --help                       Show this screen
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def plot(arguments):

    for name in arguments["<file_name>"]:
        f = open(name)
        lines = f.readlines()
        times = []

        for line in lines:
            times.append(float(line.split()[0]))
        times = np.array(times)

        plt.plot(np.arange(len(times)), times)

    plt.xlabel("Run")
    plt.ylabel("Time [s]")


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    plot(args)

    plt.show()
