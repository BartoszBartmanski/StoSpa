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

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)
    args["--time"] = int(args["--time"])
    args["--species"] = int(args["--species"])

    for name in args["<file_name>"]:
        # Get the size of square grid of errors

        original_data_shape = None
        num_dims = None

        f = open(name, "r")
        for line in f:
            if "#" not in line:
                break
            if "original_data_shape" in line:
                original_data_shape = line.split(" ")[-1].split(",")
                original_data_shape = [int(x) for x in original_data_shape]
            if "num_dims" in line and "command" not in line.lower():
                num_dims = int(line.split(" ")[-1])
        f.close()
        d = np.loadtxt(name).reshape(original_data_shape)

        plt.figure(name)

        data = d[args["--species"], args["--time"]]
        if num_dims == 1:
            plt.bar(np.arange(len(data)), data)
            plt.xlabel(r"$m$")
            plt.ylabel(r"$P_s$")
        else:
            im = plt.imshow(data, origin="lower", aspect="auto")
            plt.ylabel(r"$m_y$")
            plt.xlabel(r"$m_x$")
            plt.colorbar(im)

    plt.show()
