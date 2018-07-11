#!/usr/bin/python3

"""Visualising the data from stochastic simulations.

Usage:
    sto_plot.py -h | --help
    sto_plot.py <sim_name> [options]
    sto_plot.py <sim_name> (-m | --movie) [options]
    sto_plot.py <sim_name> (-p | --plot) [options]
    sto_plot.py <sim_name> --1d [options]
    sto_plot.py <sim_name> --2d [options]
    sto_plot.py <sim_name> --3d [options]
    sto_plot.py <sim_name> --error [options]
    sto_plot.py <sim_name> --get_sims [options]
    sto_plot.py <sim_name> --find_peaks [options]
    sto_plot.py <sim_name> --triangulation [options]
    sto_plot.py <sim_name> --turing_stats ( mean_height | mean_length | mean_area ) [options]
    sto_plot.py <sim_name> --stationary_dist [sum] [options]
    sto_plot.py <sim_name> -b | --benchmark [options]
    sto_plot.py <sim_name> -c | --check [options]
    sto_plot.py <sim_name> --pickle [<file_name>] [options]
    sto_plot.py <sim_name> [options]
    sto_plot.py --set <name> <value>
    sto_plot.py --info

Options:
    -h --help                       Show this screen
    -d=<dir>                        Defines the current working directory. [default: ./]
    <sim_name>                      Name of the simulation that will be plotted.
    -m --movie                      The simulation will be made into an animation.
    -p --plot                       The simulation will be plotted.
    -s --power_spectrum             Plot a power_spectrum for the data instead of actual data.
    -w --wave_mode_freq             Plot of the frequencies of each wave mode across different simulations.
    --1d                            A 1d plot of the power spectrum for a specified wave-mode.
    --2d                            A 2d plot of the simulation.
    --3d                            A 3d plot of the data.
    -b --benchmark                  Plot the results of the benchmarking tests.
    -c --check                      Checks that all the simulations are consistent.
    --error                         Calculate the error for each time step compared to the analytic solution.
    --get_sims                      Returns the files for which the dominant wavemode is specified by --mode.
    --find_peaks                    Finds the locations of any maximum in the data.
    --stationary_dist               Plots the stationary distribution.
    --sim_info                      Display the information about the simulation.
    --add_analytic=<path>           Open a file that contains the analytic solution.
    --save_analytic=<path>          Save the analytic solution to a file.
    --save_stochastic=<path>        Save the analytic solution to a file.
    --show_analytic                 Show the analytic (if the file for analytic solutions exists).
    --show_difference               Show the difference between the analytic and stochastic solutions.
    --interval=<fps>                Interval between frames in a movie produced based on the data. [default: 24]
    --frame=<frame_index>           Index of the frame to be picked from the data. [default: -1]
    --species=<species>             The index of the species to be selected.
    --bounds=<bounds>               Indices of the spatial points by which to bound the output data.
    --time=<time>                   Indices of time points by which to bound the data. For 2d separate by a comma.
    --mode=<mode>                   Select specific wavemode. For 2d separate by a comma.
    --remove_zero                   Removes the zeroth wave mode from power spectrum
    --save=<save_dir>               Save the plot of the simulation.
    --latex=<latex_font_size>       Changes the figure plot size to make sure that figures look LaTeX-native.
    --append=<name>                 Appends file name.
    --set_title=<title>             Change the title of the plot from the default.
    --set_x_label=<xlabel>          Change the label of the x-axis from the default.
    --set_y_label=<ylabel>          Change the label of the y-axis from the default.
    --add_h_line=<value>            Adds a horizontal line at the specified value.
    --add_v_line=<value>            Adds a vertical line at the specified value.
    --set_x_lim=<x_lim>             Sets the values of the x-limits.
    --set_y_lim=<y_lim>             Sets the values of the y-limits.

"""

import pickle
import os
from parameters import Parameters
from data import get_data
from plot import *
from docopt import docopt
from scipy import misc

if __name__ == '__main__':

    # Get the command line arguments
    args = docopt(__doc__)

    needed = ["primary_colormap",
              "secondary_colormap",
              "tertiary_colormap",
              "primary_color",
              "secondary_color",
              "tertiary_color",
              "rep",
              "fig_size",
              "aspect",
              "font_size",
              "dpi",
              "file_type"]

    # Get the saved settings
    q = Parameters(os.path.expanduser("~") + "/.sto_plot.json", needed)
    plot = None

    # Display the saved settings
    if args["--info"]:
        q.display()
    # Change one of the saved settings
    elif args["--set"]:
        q.change(args["<name>"], args["<value>"])
    # Otherwise, work with simulations data
    else:
        # Change the directory
        os.chdir(args["-d"])

        # Open the data
        if ".p" in args["<sim_name>"]:
            with open(args["<sim_name>"], "rb") as handle:
                data = pickle.load(handle)
        else:
            data = Data(args["<sim_name>"])

        if args["--pickle"]:
            if args["<file_name>"] is None:
                file_name = data.file_name.replace(".dat", ".p")
            else:
                file_name = args["<file_name>"]
                if ".p" not in file_name:
                    file_name += ".p"
            with open(file_name, "wb") as handle:
                pickle.dump(data, handle)
                print("File pickled as", file_name)

        # If any of these options is used then an analytic solution will be necessary
        bool_analytic = (args["--add_analytic"] or
                         args["--show_analytic"] or
                         args["--show_difference"] or
                         args["--error"] or
                         args["--save_analytic"])
        # If an analytic solution is necessary, add the analytic solution, either from function or a file
        if bool_analytic:
            data.add_analytic(get_data(args["--add_analytic"]))

        # If needs be, display data information
        if args["--sim_info"]:
            print(data)

        # If necessary transform the data
        if args["--power_spectrum"]:
            data.dct(args["--remove_zero"], bool_analytic)
            app_label = r"$P_s$"
        elif args["--wave_mode_freq"]:
            data.count_modes()
        else:
            data.average()

        # Correct the values by which to slice the data
        if args["--bounds"] is None:
            args["--bounds"] = []
        else:
            args["--bounds"] = [int(x) for x in args["--bounds"].split(",")]
        if args["--time"] is None:
            args["--time"] = []
        else:
            args["--time"] = [int(x) for x in args["--time"].split(",")]
        if args["--mode"] is None:
            args["--mode"] = []
        else:
            args["--mode"] = ([int(x) for x in args["--mode"].split(",")] + [0, 0])[:data.num_dims]
        mode_str = "_".join(map(str, args["--mode"]))

        try:
            args["--species"] = int(args["--species"])
        except TypeError:
            args["--species"] = None

        # Correct the frame argument, if plot rather than a movie
        if args["--plot"]:
            args["--time"] = [int(args["--frame"])+data.num_time_steps, int(args["--frame"])+data.num_time_steps+1]
            args["--time"] = [x % (data.num_time_steps+1) for x in args["--time"]]

        # Slice the data appropriately
        data.slice(args["--time"], args["--bounds"], args["--species"])

        if args["--save_analytic"] or args["--save_stochastic"]:
            if args["--save_analytic"]:
                data.save_analytic(data.file_name, args["--save_analytic"], args["--append"])
            if args["--save_stochastic"]:
                data.save_stochastic(data.file_name, args["--save_stochastic"], args["--append"])

        # Prints and/or saves to a file the simulations with the specified wavemode
        elif args["--get_sims"]:
            mode_sims = data.get_sim_names(args["--mode"])
            print(mode_sims)
            if args["--save"]:
                if os.path.isdir(args["--save"]):
                    args["--save"] += data.file_name + "_list.dat"
                with open(args["--save"], "w") as handle_sims:
                    handle_sims.write("# List of simulation names with a dominant wavemode {}\n".format(mode_str))
                    handle_sims.write(mode_sims)

        elif args["--find_peaks"]:
            peak_indices = data.find_peaks()
            data.reset()
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

        elif args["--triangulation"]:
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

        elif args["--turing_stats"]:
            stats = data.get_spot_properties()
            if args["mean_height"]:
                plot = Plot1d(stats[0], file_name=data.file_name, save_dir=args["--save"], rep="hist")
                plot.set_title("Peak heights")
                plot.append_filename("mean_height")
            elif args["mean_length"]:
                plot = Plot1d(stats[1], file_name=data.file_name, save_dir=args["--save"], rep="hist")
                plot.set_title("Average peak-to-peak distance")
                plot.append_filename("mean_length")
            elif args["mean_area"]:
                plot = Plot1d(stats[2], file_name=data.file_name, save_dir=args["--save"], rep="hist")
                plot.set_title("Average peak area")
                plot.append_filename("mean_area")
            plot.set_y_label("Frequency density")

        elif args["--stationary_dist"]:
            values, counts = data.get_stationary_dist(args["sum"])
            plot = Plot1d(values[0], counts[0], data.file_name, args["--save"], "bar")
            plot.set_x_label("number of molecules")
            plot.set_y_label("stationary dist.")
            plot.append_filename("stationary_dist")
            if args["--show_analytic"]:
                b_0 = np.sum(data.stochastic[1, 0])
                domain_size = data.domain_bounds[1] - data.domain_bounds[0]
                factor = data.prod[0] * domain_size**2 / (data.other_rates["TwoSpeciesDecay"] * b_0)
                analytical = 1/(misc.factorial(values[0])) * factor ** values[0] * np.exp(-factor)
                plot.add_plot(values[0], analytical, color=q["secondary_color"])

        # Plots the error, when comparing what's in the stochastic and analytic data containers
        elif args["--error"]:
            plot = Plot1d(data.t_values,
                          data.get_error(),
                          data.file_name,
                          args["--save"],
                          "line")

            plot.set_ax_labels(r"$t$", r"$e$")
            plot.append_filename("error")

        # Plots a 1d data, whatever it may be
        elif args["--1d"]:
            if data.num_dims == 1:
                wavemode_data = data.stochastic[0, :, args["--mode"][0]]
            else:
                wavemode_data = data.stochastic[0, :, args["--mode"][1], args["--mode"][0]]

            plot = Plot1d(data.t_values,
                          wavemode_data,
                          data.file_name,
                          args["--save"],
                          "line")

            plot.set_ax_labels(r"$t$", data.data_type)
            plot.append_filename("index_{}".format(mode_str))

        # Plots a 2d data, whatever it may be
        elif args["--2d"]:
            if data.num_dims == 1:
                plot = Plot2d(data.stochastic[0],
                              data.file_name,
                              args["--save"])

                plot.set_ax_labels(data.labels[0], data.labels[1])
                plot.set_extent(np.array([data.domain_bounds[0], data.domain_bounds[1],
                                          data.t_values[0], data.t_values[-1]]))
            else:
                plot = Plot2d(data.stochastic[0, -1],
                              data.file_name,
                              args["--save"])

                plot.set_ax_labels(data.labels[0], data.labels[1])
                plot.set_extent(np.array([data.domain_bounds[0], data.domain_bounds[1],
                                          data.domain_bounds[2], data.domain_bounds[3]]))
            plot.append_filename("2d")

        # Plots a 3d data, whatever it may be
        elif args["--3d"]:
            if args["--power_spectrum"] and data.num_dims == 1:
                plot = Plot3dLines(data.t_values,
                                   data.freq[0],
                                   data.stochastic[0].transpose(),
                                   data.file_name,
                                   args["--save"])
                plot.append_filename("power_spectrum_3d")
                plot.set_ax_labels(data.labels[1], data.labels[0], data.data_type)
            else:
                if data.num_dims == 1:
                    plot = Plot3d(data.t_values,
                                  data.x[0],
                                  data.stochastic[0].transpose(),
                                  data.file_name,
                                  args["--save"])
                    plot.set_ax_labels(data.labels[1], data.labels[0], data.data_type)

                else:
                    plot = Plot3d(data.x[0], data.x[1],
                                  data.stochastic[0, -1],
                                  data.file_name,
                                  args["--save"])
                    plot.set_ax_labels(data.labels[0], data.labels[1], data.data_type)
                plot.append_filename("3d")

        elif args["--plot"] or args["--movie"]:
            if data.num_dims == 1:
                plot = PlotSim1d(data,
                                 data.file_name,
                                 args["--save"],
                                 args["--show_analytic"],
                                 args["--show_difference"])
            else:
                plot = PlotSim2d(data,
                                 data.file_name,
                                 args["--save"],
                                 args["--show_analytic"],
                                 args["--show_difference"])

            if args["--power_spectrum"]:
                if data.num_dims == 1:
                    plot.set_ax_labels(r"$m$", r"$P_s$")
                    plot.set_domain(data.freq[0])
                elif data.num_dims == 2:
                    plot.set_ax_labels(r"$m_x$", r"$m_y$")
                    plot.set_domain(data.freq)

                plot.append_filename("power_spectrum")

            elif args["--wave_mode_freq"]:
                if data.num_dims == 1:
                    plot.set_ax_labels(r"$m$", r"Frequency")
                    plot.set_domain(data.freq[0])
                elif data.num_dims == 2:
                    plot.set_ax_labels(r"$m_x$", r"$m_y$")
                    plot.set_domain(data.freq)

                plot.append_filename("wavemode_freq")

            if hasattr(plot, "set_stochastic_colour"):
                plot.set_stochastic_colour(q["primary_color"])
                plot.set_analytic_colour(q["secondary_color"])
                plot.set_representation(q["rep"])

            if hasattr(plot, "set_stochastic_colourmap"):
                plot.set_stochastic_colourmap(q["primary_colormap"])
                plot.set_analytic_colourmap(q["secondary_colormap"])
                plot.set_difference_colourmap(q["tertiary_colormap"])

        if hasattr(plot, "alt_ticks") and not isinstance(plot, (Plot1d, Plot3d, Plot3dLines)):
            if args["--power_spectrum"] or args["--wave_mode_freq"]:
                if data.num_dims == 1:
                    plot.alt_ticks(0, start=data.freq[0][0], end=data.freq[0][-1])
                else:  # data.num_dims == 2
                    plot.alt_ticks(start=[data.freq[0][0], data.freq[1][0]],
                                   end=[data.freq[0][-1], data.freq[1][-1]])

        if hasattr(plot, "set_colour"):
            plot.set_colour(q["primary_color"])

        if hasattr(plot, "set_colormap"):
            plot.set_colormap(q["primary_colormap"])

        if hasattr(plot, "colorbar") and not isinstance(plot, (PlotSim1d, Plot1d)):
            plot.colorbar()

        if hasattr(plot, "set_interval"):
            plot.set_interval(int(args["--interval"]))

        # From this point on are the settings that apply for all the python scripts
        if hasattr(plot, "set_aspect"):
            plot.set_aspect(q["aspect"])

        if hasattr(plot, "set_fig_size"):
            plot.set_fig_size(float(q["fig_size"]))

        if hasattr(plot, "set_x_lim") and args["--set_x_lim"]:
            args["--set_x_lim"] = [float(x) for x in args["--set_x_lim"]]
            bounds = [data.domain_bounds[0], data.domain_bounds[1]]
            args["--set_x_lim"] += bounds[len(args["--set_x_lim"]):2]
            plot.set_x_lim(args["--set_x_lim"])

        if hasattr(plot, "set_y_lim") and args["--set_y_lim"]:
            args["--set_y_lim"] = [float(x) for x in args["--set_y_lim"]]
            bounds = [0.9 * np.min(data.min_value), 1.1 * np.max(data.max_value[0])]
            args["--set_y_lim"] += bounds[len(args["--set_y_lim"]):2]
            plot.set_y_lim(args["--set_y_lim"])

        if hasattr(plot, "append_filename") and args["--append"]:
            plot.append_filename(args["--append"])

        if hasattr(plot, "set_title") and args["--set_title"]:
            plot.set_title(r"{}".format(args["--set_title"]))

        if hasattr(plot, "set_x_label") and args["--set_x_label"]:
            plot.set_x_label(r"{}".format(args["--set_x_label"]))

        if hasattr(plot, "set_y_label") and args["--set_y_label"]:
            plot.set_y_label(r"{}".format(args["--set_y_label"]))

        if hasattr(plot, "add_h_line") and args["--add_h_line"]:
            plot.add_h_line(float(args["--add_h_line"]))

        if hasattr(plot, "add_v_line") and args["--add_v_line"]:
            plot.add_v_line(float(args["--add_v_line"]))

        if hasattr(plot, "update_all_text"):
            plot.update_all_text(float(q["font_size"]))

        if hasattr(plot, "show"):
            plot.show(float(q["dpi"]), q["file_type"])
