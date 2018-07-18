#!/usr/bin/python3

"""Visualising the data from stochastic simulations.

Usage:
    plot_peaks.py -h | --help
    plot_peaks.py <file_name> find_peaks [options]
    plot_peaks.py <file_name> triangulation [options]
    plot_peaks.py <file_name> turing_stats ( height | length | area ) [options]

Options:
    -h --help                       Show this screen.
    
"""

from docopt import docopt
from data import Data
from matplotlib import pyplot as plt
import pickle

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    # Open the data
    if ".p" in args["<file_name>"]:
        with open(args["<file_name>"], "rb") as handle:
            data = pickle.load(handle)
    else:
        data = Data(args["<file_name>"])

    if data.num_dims == 1:
        plot = Plot1d(data.domain[0], data.stochastic[0, -1], data.file_name, args["--save"], q["rep"])
    else:
        plot = Plot2d(data.stochastic[0, -1], data.file_name, args["--save"])

    if args["find_peaks"]:
        peak_indices = data.find_peaks()
        if data.num_dims == 1:
            plot = Plot1d(data.domain[0], data.stochastic[0, -1], data.file_name, args["--save"], q["rep"])
            plot.set_x_lim(data.domain_bounds)
            plot.add_plot(data.domain[0][peak_indices[0]], data.stochastic[0, -1][peak_indices[0]],
                          ls="x", color=q["tertiary_color"], rep="scatter")
            plot.set_ax_labels(data.labels[0], data.data_type)
        else:
            plot = Plot2d(data.stochastic[0, -1], data.file_name, args["--save"])
            plot.add_plot(data.domain[0][peak_indices[0][1]], data.domain[1][peak_indices[0][0]],
                          ls="x", color=q["tertiary_color"], rep="scatter")
            plot.set_ax_labels(data.labels[0], data.labels[1])
            plot.set_extent(np.array([data.domain_bounds[0], data.domain_bounds[1],
                                      data.domain_bounds[2], data.domain_bounds[3]]))

    elif args["triangulation"]:
        peak_indices = data.find_peaks()
        if data.num_dims == 1:
            plot = Plot1d(data.domain[0], data.stochastic[0, -1], data.file_name, args["--save"], q["rep"])
            plot.set_x_lim(data.domain_bounds)
            plot.add_plot(data.domain[0][peak_indices[0]], data.stochastic[0, -1][peak_indices[0]],
                          ls="x", color=q["tertiary_color"], rep="scatter")
            plot.set_ax_labels(data.labels[0], data.data_type)
        else:
            points = np.array([data.domain[0][peak_indices[0][1]], data.domain[1][peak_indices[0][0]]])

            plot = Plot2d(data.stochastic[0, -1], data.file_name, args["--save"])
            plot.add_plot(points[0], points[1], rep="triangulation", color=q["tertiary_color"])
            plot.set_ax_labels(data.labels[0], data.labels[1])
            plot.set_extent(np.array([data.domain_bounds[0], data.domain_bounds[1],
                                      data.domain_bounds[2], data.domain_bounds[3]]))

    elif args["turing_stats"]:
        stats = data.get_spot_properties()
        if args["height"]:
            plot = Plot1d(stats[0], file_name=data.file_name, save_dir=args["--save"], rep="hist")
            plot.set_title("Peak heights")
            plot.append_filename("height")
        elif args["length"]:
            plot = Plot1d(stats[1], file_name=data.file_name, save_dir=args["--save"], rep="hist")
            plot.set_title("Average peak-to-peak distance")
            plot.append_filename("length")
        elif args["area"]:
            plot = Plot1d(stats[2], file_name=data.file_name, save_dir=args["--save"], rep="hist")
            plot.set_title("Average peak area")
            plot.append_filename("area")
        plot.set_y_label("Frequency density")
