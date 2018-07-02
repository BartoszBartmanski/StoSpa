#!/usr/bin/python3

"""Plots the power spectra.

Usage:
    plot_power_spectra.py -h | --help
    plot_power_spectra.py <file_name>... [options]

Options:
    -h --help                       Show this screen
    --time=<time>                   Index of the time-step [default: -1]
    --species=<species>             Index of the species [default: 0]
"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def get_data(file_name):
    original_data_shape = None

    f = open(file_name, "r")
    for line in f:
        if "#" not in line:
            break
        if "original_data_shape" in line:
            original_data_shape = line.split(" ")[-1].split(",")
            original_data_shape = [int(x) for x in original_data_shape]
    f.close()

    data = np.loadtxt(file_name).reshape(original_data_shape)

    return data


def plot(arguments):

    for name in arguments["<file_name>"]:
        # Get the size of square grid of errors

        d = get_data(name)
        num_dims = len(d[0, 0].shape)
        d = d[arguments["--species"], arguments["--time"]]

        plt.figure(name)

        if num_dims == 1:
            plt.bar(np.arange(len(d)), d)
            plt.xlabel(r"$m$")
            plt.ylabel(r"$P_s$")
        else:
            im = plt.imshow(d, origin="lower", aspect="auto", cmap="Greys")
            plt.ylabel(r"$m_y$")
            plt.xlabel(r"$m_x$")
            plt.colorbar(im)


if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)
    args["--time"] = int(args["--time"])
    args["--species"] = int(args["--species"])

    plot(args)

    plt.show()
