#!/usr/bin/env python3

"""Plots the benchmarking data.

Usage:
    plot_benchmarks.py -h | --help
    plot_benchmarks.py <filename>... [options]

Options:
    -h --help                       Show this screen
"""

import matplotlib.pyplot as plt
import numpy as np

from docopt import docopt


def get_data(name):

    f = open(name)
    lines = f.readlines()
    times = []

    for line in lines:
        times.append(float(line.split()[0]))
    times = np.array(times)

    return times


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    for name in arg["<filename>"]:
        data = get_data(name)

        plt.plot(np.arange(len(data)), data)

    plt.xlabel("Run")
    plt.ylabel("Time [s]")
    plt.show()
