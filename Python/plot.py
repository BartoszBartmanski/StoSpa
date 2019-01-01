#!/usr/bin/python3

import numpy as np
from matplotlib import animation
from matplotlib import figure
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import tri
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from os import path

from data import Data

ints = (int, np.int, np.int8, np.int16, np.int32, np.int64)
floats = (float, np.float, np.float16, np.float32, np.float64, np.float128)


def get_app_string(val, val_name="t"):
    """Calculates the number of digits necessary (used in time display in figures)"""

    number_as_string = str(float(val))
    str_len = len(number_as_string[number_as_string.index(".") + 1:])
    string = val_name + " = {:." + str(str_len) + "f}"

    return string


def save_or_show(obj, file_name, save_dir, file_type="pdf", dpi=None, data_array=None, file_header=""):
    assert isinstance(file_name, str) or file_name is None
    assert isinstance(save_dir, str) or save_dir is None
    assert isinstance(data_array, (np.ndarray, list, tuple)) or data_array is None
    assert isinstance(file_header, str)

    file_name = file_name.replace(" ", "_")

    if save_dir is not None:
        print("Plot saved as", file_name)
        full_path = save_dir + "/" + file_name
        full_path = full_path.replace("//", "/")

        if path.isdir(save_dir):
            # Save the plot
            if isinstance(obj, animation.FuncAnimation):
                obj.save(full_path + ".avi")
                full_path += "_movie"
            elif isinstance(obj, figure.Figure):
                obj.savefig(full_path + "." + file_type, dpi=dpi)
                full_path += "_plot"
            else:
                print("\nERROR: Unrecognised class submitted for the obj parameter in save_or_show function!")
                quit()

            if data_array is not None:
                if isinstance(data_array, np.ndarray):
                    with open(full_path + ".dat", "wb") as handle:
                        np.savetxt(handle, data_array, header=file_header)
                elif isinstance(data_array, (list, tuple)):
                    for counter in range(len(data_array)):
                        assert isinstance(data_array[counter], np.ndarray)
                        with open(full_path + "_{}.dat".format(counter), "wb") as handle:
                            np.savetxt(handle, data_array[counter], header=file_header)

        else:
            print("\nERROR: Directory specified ({}) is not valid!".format(save_dir))
            quit()
    else:
        if isinstance(obj, figure.Figure):
            obj.canvas.set_window_title(file_name)
        plt.show()


class Axes2d(object):
    def __init__(self, num_x, num_y, array_shape=None, filename="Figure", save_dir=None, string=None):
        assert isinstance(num_x, int)
        assert isinstance(num_y, int)
        assert (isinstance(array_shape, (tuple, list, np.ndarray)) and len(array_shape) == 2) or array_shape is None
        assert isinstance(filename, str)
        assert isinstance(save_dir, str) or save_dir is None

        self.file_name = filename
        self.save_dir = save_dir
        self.num_axes_x = num_x  # Number of axes horizontally in the figure environment
        self.num_axes_y = num_y  # Number of axes vertically in the figure environment
        self.array_shape = array_shape
        self.extent = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.titles = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.x_labels = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.y_labels = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.labels = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.axes = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.size = 6
        self.font_size = np.array([[8] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.tick_size = np.array([[8] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.cbar_font_size = np.array([[8] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.cbar_tick_size = np.array([[8] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.legend_font_size = np.array([[8] * self.num_axes_x for _ in range(self.num_axes_y)])
        self.latex_ratio = 1.0
        self.representation = "line"  # Optional attribute. Useful if 1d data is plotted.
        self.desc_string = string
        self.data = []  # Container for arrays that are not of the same size, so other plotting classes cannot be used.

        ratios = self.size * np.ones(self.num_axes_y)
        self.desc = string is not None
        if self.desc:
            ratios = np.append(ratios, [0.1])

        self.fig = plt.figure(figsize=(self.size * self.num_axes_x, (self.size + 0.1*int(self.desc)) * self.num_axes_y))
        gs = gridspec.GridSpec(self.num_axes_y + int(self.desc), self.num_axes_x, height_ratios=ratios)

        for i in range(self.num_axes_y):
            for j in range(self.num_axes_x):
                self.axes[i, j] = plt.subplot(gs[i, j])

        self.text = None
        if self.desc:
            self.ax_param = plt.subplot(gs[-1, :])
            self.ax_param.axis('off')
            self.text = self.ax_param.text(0, 0, self.desc_string, size=15)

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_save_dir(self, path_to_dir):
        assert isinstance(path_to_dir, str)
        self.save_dir = path_to_dir

    def __set_fig_size__(self, width=1.0, height=1.0):
        dims = np.array([self.size * width * self.latex_ratio * self.num_axes_x,
                         self.size * height * self.latex_ratio * self.num_axes_y + int(self.desc)])
        dims /= width

        self.fig.set_size_inches(dims[0], dims[1], forward=True)

    def set_fig_size(self, size, latex_font_size=None):
        assert isinstance(size, (float, int))
        self.size = size
        # For given latex font size, there is an associated maximum available width
        latex_font_sizes = {10: 345, 11: 360, 12: 390}
        inches_per_pt = 1.0 / 72.27

        if latex_font_size:
            self.size *= 0.1  # To turn this into a fraction of the available width
            self.latex_ratio = latex_font_sizes.get(latex_font_size) * inches_per_pt / self.num_axes_x

        self.__set_fig_size__()

    def set_font_size(self, size, ax_index=None):
        assert isinstance(size, (float, int))
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.font_size[tuple(ax_index)] = size
            if self.titles[tuple(ax_index)] is not None:
                self.set_title(self.titles[tuple(ax_index)], self.font_size[tuple(ax_index)], tuple(ax_index))
            if self.x_labels[tuple(ax_index)] is not None:
                self.set_x_label(self.x_labels[tuple(ax_index)], self.font_size[tuple(ax_index)], tuple(ax_index))
            if self.y_labels[tuple(ax_index)] is not None:
                self.set_y_label(self.y_labels[tuple(ax_index)], self.font_size[tuple(ax_index)], tuple(ax_index))
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.font_size[i, j] = size
                    if self.titles[i, j] is not None:
                        self.set_title(self.titles[i, j], self.font_size[i, j], (i, j))
                    if self.x_labels[i, j] is not None:
                        self.set_x_label(self.x_labels[i, j], self.font_size[i, j], (i, j))
                    if self.y_labels[i, j] is not None:
                        self.set_y_label(self.y_labels[i, j], self.font_size[i, j], (i, j))

    def set_tick_size(self, size, ax_index=None):
        assert isinstance(size, (int, float))
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.tick_size[tuple(ax_index)] = size
            self.axes[tuple(ax_index)].tick_params(labelsize=self.tick_size[tuple(ax_index)])
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.tick_size[i, j] = size
                    self.axes[i, j].tick_params(labelsize=self.tick_size[i, j])

    def set_colorbar_font_size(self, size, ax_index=None):
        assert isinstance(size, (int, float))
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.cbar_font_size[tuple(ax_index)] = size
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.cbar_font_size[i, j] = size

    def set_colorbar_tick_size(self, size, ax_index=None):
        assert isinstance(size, (int, float))
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.cbar_tick_size[tuple(ax_index)] = size
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.cbar_tick_size[i, j] = size

    def set_legend_font_size(self, size, ax_index=None):
        assert isinstance(size, (int, float))
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.legend_font_size[tuple(ax_index)] = size
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.legend_font_size[i, j] = size

    def update_all_text(self, size):
        assert isinstance(size, (ints, floats))
        self.set_font_size(size)
        self.set_tick_size(size)
        self.set_colorbar_font_size(size)
        self.set_colorbar_tick_size(size)
        self.set_legend_font_size(size)
        if self.desc:
            self.change_string(self.desc_string, size)

    def append_filename(self, string):
        assert isinstance(string, str)
        self.file_name += ("_"+string)

    def change_string(self, string, font_size=None):
        assert isinstance(string, str)
        assert self.desc is True, "Description axis has not been initialised!"
        self.text.set_text(string)
        self.text.set_fontsize(font_size)

    def legend(self):
        for i in range(self.num_axes_y):
            for j in range(self.num_axes_x):
                if self.labels[i, j] is not None:
                    self.axes[i, j].legend(prop={"size": self.legend_font_size[i, j]})

    def set_extent(self, values, ax_index=None):
        assert isinstance(values, (tuple, list, np.ndarray))
        values = np.array(values)
        assert values.ndim == 1 and len(values) == 4

        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            ax_index = tuple(ax_index)
            self.extent[ax_index] = values
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.extent[i, j] = values

    def alt_ticks(self, axis=None, tick_freq=None, start=None, end=None, step=None, ax_index=None, ar_shape=None):
        assert isinstance(tick_freq, (ints, tuple, list, np.ndarray)) or tick_freq is None
        assert isinstance(start, (ints, floats, tuple, list, np.ndarray)) or start is None
        assert isinstance(end, (ints, floats, tuple, list, np.ndarray)) or end is None
        assert isinstance(step, (ints, floats, tuple, list, np.ndarray)) or step is None
        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        assert isinstance(ar_shape, (tuple, list, np.ndarray)) or ar_shape is None

        if ar_shape is not None:
            arr_size = (ar_shape[0], ar_shape[1])
        else:
            assert self.array_shape is not None
            arr_size = (self.array_shape[1], self.array_shape[0])

        if start is None:
            start = np.array([0, 0], dtype=int)
        elif isinstance(start, (ints, floats)):
            start = np.array([start, start], dtype=int)
        else:
            start = np.array(start, dtype=int)
        if end is None:
            end = np.array([arr_size[0], arr_size[1]], dtype=int)
        elif isinstance(end, (ints, floats)):
            end = np.array([end, end], dtype=int)
        else:
            end = np.array(end, dtype=int)
        if step is None:
            step = np.array([1, 1], dtype=int)
        elif isinstance(step, (ints, floats)):
            step = np.array([step, step], dtype=int)
        else:
            step = np.array(step, dtype=int)

        tick_labels = [np.arange(start[0], end[0]+1, step[0]), np.arange(start[1], end[1]+1, step[1])]
        if axis == 0:
            ticks = tick_labels  # for 1d
        else:
            ticks = [np.arange(len(tick_labels[0])), np.arange(len(tick_labels[1]))]  # for 2d

        if isinstance(tick_freq, int):
            tick_freq = (tick_freq, tick_freq)
        elif tick_freq is None:
            tick_freq = [0, 0]
            tick_freq[0] = max(arr_size[0] // 10, 1)
            tick_freq[1] = max(arr_size[1] // 10, 1)

        assert isinstance(ax_index, (tuple, list, np.ndarray)) or ax_index is None
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            ax_index = tuple(ax_index)
            if axis == 0:
                self.set_x_lim([ticks[0][0]-0.5*step[0], ticks[0][-1]+0.5*step[0]], ax_index)
                self.axes[ax_index].set_xticks(ticks[0][::tick_freq[0]])
                self.axes[ax_index].set_xticklabels(tick_labels[0][::tick_freq[0]])
            elif axis == 1:
                self.set_y_lim([ticks[1][0]-0.5*step[1], ticks[1][-1]+0.5*step[1]], ax_index)
                self.axes[ax_index].set_yticks(ticks[1][::tick_freq[1]])
                self.axes[ax_index].set_yticklabels(tick_labels[1][::tick_freq[1]])
            else:
                self.set_x_lim([ticks[0][0]-0.5*step[0], ticks[0][-1]+0.5*step[0]], ax_index)
                self.axes[ax_index].set_xticks(ticks[0][::tick_freq[0]])
                self.axes[ax_index].set_xticklabels(tick_labels[0][::tick_freq[0]])

                self.set_y_lim([ticks[1][0]-0.5*step[1], ticks[1][-1]+0.5*step[1]], ax_index)
                self.axes[ax_index].set_yticks(ticks[1][::tick_freq[1]])
                self.axes[ax_index].set_yticklabels(tick_labels[1][::tick_freq[1]])

        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if axis == 0:
                        self.set_x_lim([ticks[0][0]-0.5*step[0], ticks[0][-1]+0.5*step[0]], (i, j))
                        self.axes[i, j].set_xticks(ticks[0][::tick_freq[0]])
                        self.axes[i, j].set_xticklabels(tick_labels[0][::tick_freq[0]])
                    elif axis == 1:
                        self.set_y_lim([ticks[1][0]-0.5*step[1], ticks[1][-1]+0.5*step[1]], (i, j))
                        self.axes[i, j].set_yticks(ticks[1][::tick_freq[1]])
                        self.axes[i, j].set_yticklabels(tick_labels[1][::tick_freq[1]])
                    else:
                        self.set_x_lim([ticks[0][0]-0.5*step[0], ticks[0][-1]+0.5*step[0]], (i, j))
                        self.axes[i, j].set_xticks(ticks[0][::tick_freq[0]])
                        self.axes[i, j].set_xticklabels(tick_labels[0][::tick_freq[0]])

                        self.set_y_lim([ticks[1][0]-0.5*step[1], ticks[1][-1]+0.5*step[1]], (i, j))
                        self.axes[i, j].set_yticks(ticks[1][::tick_freq[1]])
                        self.axes[i, j].set_yticklabels(tick_labels[1][::tick_freq[1]])

        self.extent = np.array([[None] * self.num_axes_x for _ in range(self.num_axes_y)])

    def set_labels(self, label, ax_index=None):
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.labels[tuple(ax_index)] = label
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.labels[i, j] = label

    def colorbar(self, label="", ax_index=None):
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.labels[tuple(ax_index)] = label
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.labels[i, j] = label

    def add_h_line(self, value, ls="--", colour="black", ax_index=None):
        assert isinstance(value, (int, float))
        assert isinstance(colour, str)
        assert isinstance(ls, str)

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].axhline(value, linestyle=ls, color=colour)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].axhline(value, linestyle=ls, color=colour)

    def add_v_line(self, value, ls="--", colour="black", ax_index=None):
        assert isinstance(value, (int, float))
        assert isinstance(colour, str)
        assert isinstance(ls, str)

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].axvline(value, linestyle=ls, color=colour)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].axvline(value, linestyle=ls, color=colour)

    def set_x_lim(self, values, ax_index=None):
        assert isinstance(values, (tuple, list, np.ndarray))
        values = np.array(values)
        assert len(values) == 2 and values.ndim == 1

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].set_xlim(values)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].set_xlim(values)

    def set_y_lim(self, values, ax_index=None):
        assert isinstance(values, (tuple, list, np.ndarray))
        values = np.array(values)
        assert len(values) == 2 and values.ndim == 1

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].set_ylim(values)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].set_ylim(values)

    def set_x_label(self, label, font_size=None, ax_index=None):
        assert isinstance(label, str)
        assert isinstance(font_size, (floats, ints)) or font_size is None

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            if font_size is None:
                font_size = self.font_size[tuple(ax_index)]
            self.x_labels[tuple(ax_index)] = label
            self.axes[tuple(ax_index)].set_xlabel(label, fontsize=font_size)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if font_size is None:
                        font_size = self.font_size[i, j]
                    self.x_labels[i, j] = label
                    self.axes[i, j].set_xlabel(label, fontsize=font_size)

    def set_y_label(self, label, font_size=None, ax_index=None):
        assert isinstance(label, str)
        assert isinstance(font_size, (floats, ints)) or font_size is None

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            if font_size is None:
                font_size = self.font_size[tuple(ax_index)]
            self.y_labels[tuple(ax_index)] = label
            self.axes[tuple(ax_index)].set_ylabel(label, fontsize=font_size)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if font_size is None:
                        font_size = self.font_size[i, j]
                    self.y_labels[i, j] = label
                    self.axes[i, j].set_ylabel(label, fontsize=font_size)

    def set_ax_labels(self, x_label="", y_label="", font_size=None, ax_index=None):
        assert isinstance(x_label, str)
        assert isinstance(y_label, str)
        assert isinstance(font_size, (floats, ints)) or font_size is None

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            if font_size is None:
                font_size = self.font_size[tuple(ax_index)]
            self.y_labels[tuple(ax_index)] = y_label
            self.axes[tuple(ax_index)].set_ylabel(y_label, fontsize=font_size)
            self.x_labels[tuple(ax_index)] = x_label
            self.axes[tuple(ax_index)].set_xlabel(x_label, fontsize=font_size)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if font_size is None:
                        font_size = self.font_size[i, j]
                    self.y_labels[i, j] = y_label
                    self.axes[i, j].set_ylabel(y_label, fontsize=font_size)
                    self.x_labels[i, j] = x_label
                    self.axes[i, j].set_xlabel(x_label, fontsize=font_size)

    def set_title(self, title, font_size=None, ax_index=None):
        assert isinstance(title, str), "title = {}".format(title)
        assert isinstance(font_size, (floats, ints)) or font_size is None, "font_size = {}".format(font_size)

        if isinstance(ax_index, (tuple, list, np.ndarray)):
            if font_size is None:
                font_size = self.font_size[tuple(ax_index)]
            self.titles[tuple(ax_index)] = title
            self.axes[tuple(ax_index)].set_title(title, fontsize=font_size)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if font_size is None:
                        font_size = self.font_size[i, j]
                    self.titles[i, j] = title
                    self.axes[i, j].set_title(title, fontsize=font_size)

    def use_x_log_scale(self, ax_index=None):
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].set_xscale("log")
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].set_xscale("log")

    def use_y_log_scale(self, ax_index=None):
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].set_yscale("log")
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].set_yscale("log")

    def set_representation(self, rep):
        assert isinstance(rep, str)
        error = "Parameter representation can only be:\n- line\n- step\n- bar\n"
        assert rep == "line" or rep == "bar" or rep == "step", error
        self.representation = rep

    def grid(self, val=True, ax_index=None):
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.axes[tuple(ax_index)].grid(val)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    self.axes[i, j].grid(val)

    def add_plot(self, x, y=None, ls="-", label=None, color=None, rep="line", ax_index=None):
        if y is None:
            y = x
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            if rep == "line":
                self.axes[tuple(ax_index)].plot(x, y, ls=ls, color=color, label=label)
            elif rep == "bar":
                dx = x[1] - x[0]
                self.axes[tuple(ax_index)].bar(x, y, width=dx, edgecolor="black", color=color, label=label)
            elif rep == "step":
                self.axes[tuple(ax_index)].step(x, y, color=color, label=label)
            elif rep == "scatter":
                self.axes[tuple(ax_index)].scatter(x, y, marker=ls, c=color, label=label)
            elif rep == "triangulation":
                self.axes[tuple(ax_index)].triplot(tri.Triangulation(x, y), color=color)
            elif rep == "hist":
                self.axes[tuple(ax_index)].hist(x, color=color, label=label)
        else:
            for i in range(self.num_axes_y):
                for j in range(self.num_axes_x):
                    if rep == "line":
                        self.axes[i, j].plot(x, y, ls=ls, color=color, label=label)
                    elif rep == "bar":
                        dx = x[1] - x[0]
                        self.axes[i, j].bar(x, y, width=dx, edgecolor="black", color=color, label=label)
                    elif rep == "step":
                        self.axes[i, j].step(x, y, color=color, label=label)
                    elif rep == "scatter":
                        self.axes[i, j].scatter(x, y, marker=ls, c=color, label=label)
                    elif rep == "triangulation":
                        self.axes[i, j].triplot(tri.Triangulation(x, y), color=color)
                    elif rep == "hist":
                        self.axes[i, j].hist(x, color=color, label=label)

        self.data.append(np.array([x, y], dtype=float))

    def show(self, dpi=100, file_type="pdf"):
        self.legend()
        plt.tight_layout()
        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi, self.data)


class Axes3d(object):

    def __init__(self, array_shape=None, filename="Figure", save_dir=None):

        assert isinstance(array_shape, (tuple, list, np.ndarray)) and len(array_shape) == 2 or array_shape is None
        assert isinstance(filename, str)
        assert isinstance(save_dir, str) or save_dir is None

        self.size = 7
        self.font_size = 8
        self.tick_size = 8
        self.cbar_font_size = 8
        self.cbar_tick_size = 8
        self.legend_font_size = 8
        self.latex_ratio = 1.0
        self.title = None
        self.x_label = None
        self.y_label = None
        self.z_label = None
        self.fig = plt.figure(figsize=(self.size, self.size))
        self.ax = Axes3D(self.fig)

        self.view_elev = 20
        self.view_azim = -130
        self.view_dist = 12
        self.ax.view_init(self.view_elev, self.view_azim)
        self.ax.dist = self.view_dist

        self.array_shape = array_shape
        self.file_name = filename
        self.save_dir = save_dir

        self.lw = 0

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_save_dir(self, path_to_dir):
        assert isinstance(path_to_dir, str)
        self.save_dir = path_to_dir

    def __set_fig_size__(self, width=1.0, height=1.0):
        dims = np.array([self.size * width * self.latex_ratio, self.size * height * self.latex_ratio])
        dims /= width

        self.fig.set_size_inches(dims[0], dims[1], forward=True)

    def set_fig_size(self, size, latex_font_size=None):
        assert isinstance(size, (float, int))
        self.size = size

        latex_font_sizes = {10: 345, 11: 360, 12: 390}
        inches_per_pt = 1.0 / 72.27

        if latex_font_size:
            self.size *= 0.1  # To turn this into a fraction of the available width
            self.latex_ratio = latex_font_sizes.get(latex_font_size) * inches_per_pt

        self.__set_fig_size__()

    def set_font_size(self, size):
        assert isinstance(size, (float, int))
        self.font_size = size
        if self.title is not None:
            self.set_title(self.title, self.font_size)
        if self.x_label is not None:
            self.set_x_label(self.x_label, self.font_size)
        if self.y_label is not None:
            self.set_y_label(self.y_label, self.font_size)
        if self.z_label is not None:
            self.set_z_label(self.z_label, self.font_size)

    def set_tick_size(self, size):
        assert isinstance(size, (int, float))
        self.tick_size = size
        self.ax.tick_params(labelsize=self.tick_size)

    def set_colorbar_font_size(self, size):
        assert isinstance(size, (int, float))
        self.cbar_font_size = size

    def set_colorbar_tick_size(self, size):
        assert isinstance(size, (int, float))
        self.cbar_tick_size = size

    def set_legend_font_size(self, size):
        assert isinstance(size, (int, float))
        self.legend_font_size = size

    def update_all_text(self, size):
        self.set_font_size(size)
        self.set_tick_size(size)
        self.set_colorbar_font_size(size)
        self.set_colorbar_tick_size(size)
        self.set_legend_font_size(size)

    def legend(self):
        self.ax.legend(prop={"size": self.legend_font_size})

    def append_filename(self, string):
        assert isinstance(string, str)
        self.file_name += ("_"+string)

    def set_x_lim(self, values):
        assert len(values) == 2
        assert isinstance(values[0], (int, float)) or values[0] is None
        assert isinstance(values[1], (int, float)) or values[1] is None
        self.ax.set_xlim(values)

    def set_y_lim(self, values):
        assert len(values) == 2
        assert isinstance(values[0], (int, float)) or values[0] is None
        assert isinstance(values[1], (int, float)) or values[1] is None
        self.ax.set_ylim(values)

    def set_z_lim(self, values):
        assert len(values) == 2
        assert isinstance(values[0], (int, float)) or values[0] is None
        assert isinstance(values[1], (int, float)) or values[1] is None
        self.ax.set_zlim(values)

    def use_x_log_scale(self):
        self.ax.set_xscale("log")

    def use_y_log_scale(self):
        self.ax.set_yscale("log")

    def use_z_log_scale(self):
        self.ax.set_zscale("log")

    def set_x_label(self, label, font_size=None):
        assert isinstance(label, str)
        self.x_label = label
        self.ax.set_xlabel(label, fontsize=font_size)

    def set_y_label(self, label, font_size=None):
        assert isinstance(label, str)
        self.y_label = label
        self.ax.set_ylabel(label, fontsize=font_size)

    def set_z_label(self, label, font_size=None):
        assert isinstance(label, str)
        self.z_label = label
        self.ax.set_zlabel(label, fontsize=font_size)

    def set_ax_labels(self, x_label="", y_label="", z_label="", font_size=None):
        assert isinstance(x_label, str)
        assert isinstance(y_label, str)
        assert isinstance(font_size, (floats, ints)) or font_size is None

        if font_size is None:
            font_size = self.font_size
        self.x_label = x_label
        self.ax.set_xlabel(x_label, fontsize=font_size)
        self.y_label = y_label
        self.ax.set_ylabel(y_label, fontsize=font_size)
        self.z_label = z_label
        self.ax.set_zlabel(z_label, fontsize=font_size)

    def set_title(self, title, font_size=None):
        assert isinstance(title, str)
        self.title = title
        self.ax.set_title(title, fontsize=font_size)

    def set_view(self, elev, azim, dist):
        assert isinstance(elev, (ints, floats))
        assert isinstance(azim, (ints, floats))
        assert isinstance(dist, (ints, floats))
        self.view_elev = elev
        self.view_azim = azim
        self.view_dist = dist
        self.ax.view_init(elev, azim)
        self.ax.dist = dist

    def set_line_width(self, width):
        assert isinstance(width, (float, int))
        self.lw = width

    def alter_ticks(self):
        # Not sure what this does
        self.ax.zaxis.set_major_locator(LinearLocator(10))
        self.ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    def alt_ticks(self, axis=None, tick_freq=None, start=None, end=None, step=None, ar_shape=None):
        assert isinstance(tick_freq, (ints, tuple, list, np.ndarray)) or tick_freq is None
        assert isinstance(start, (ints, floats, tuple, list, np.ndarray)) or start is None
        assert isinstance(end, (ints, floats, tuple, list, np.ndarray)) or end is None
        assert isinstance(step, (ints, floats, tuple, list, np.ndarray)) or step is None
        assert isinstance(ar_shape, (tuple, list, np.ndarray)) or ar_shape is None

        if ar_shape is not None:
            arr_size = (ar_shape[0], ar_shape[1])
        else:
            assert self.array_shape is not None
            arr_size = (self.array_shape[1], self.array_shape[0])

        if start is None:
            start = np.array([0, 0], dtype=int)
        elif isinstance(start, (int, float)):
            start = np.array([start, start], dtype=int)
        else:
            start = np.array(start, dtype=int)
        if end is None:
            end = np.array([arr_size[1], arr_size[0]], dtype=int)
        elif isinstance(end, (int, float)):
            end = np.array([end, end], dtype=int)
        else:
            end = np.array(end, dtype=int)
        if step is None:
            step = np.array([1, 1], dtype=int)
        elif isinstance(step, (int, float)):
            step = np.array([step, step], dtype=int)
        else:
            step = np.array(step, dtype=int)

        ticks = [np.arange(arr_size[0]), np.arange(arr_size[1])]
        tick_labels = [np.arange(start[0], end[0], step[0]), np.arange(start[1], end[1], step[1])]

        if isinstance(tick_freq, int):
            tick_freq = (tick_freq, tick_freq)
        elif tick_freq is None:
            tick_freq = [0, 0]
            tick_freq[0] = max(arr_size[0] // 10, 1)
            tick_freq[1] = max(arr_size[1] // 10, 1)

        if axis == 0:
            self.ax.set_xticks(ticks[0][::tick_freq[0]])
            self.ax.set_xticklabels(tick_labels[0][::tick_freq[0]])
        elif axis == 1:
            self.ax.set_yticks(ticks[1][::tick_freq[1]])
            self.ax.set_yticklabels(tick_labels[1][::tick_freq[1]])
        else:
            self.ax.set_xticks(ticks[0][::tick_freq[0]])
            self.ax.set_xticklabels(tick_labels[0][::tick_freq[0]])

            self.ax.set_yticks(ticks[1][::tick_freq[1]])
            self.ax.set_yticklabels(tick_labels[1][::tick_freq[1]])

        # plt.xticks(rotation=-90)

    def show(self, dpi=100, file_type="pdf"):
        self.legend()
        plt.tight_layout()
        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi)


class PlotSim1d(Axes2d):
    def __init__(self,
                 sims,
                 file_name="sim",
                 save_dir=None,
                 show_analytic=False,
                 show_difference=False):

        # First check the input and assign them to member attributes
        assert isinstance(sims, Data)
        self.sims = sims
        assert sims.num_dims == 1
        assert isinstance(show_analytic, bool)
        self.show_analytic = show_analytic
        assert isinstance(show_difference, bool)
        self.show_difference = show_difference

        assert isinstance(file_name, str)
        assert isinstance(save_dir, str) or save_dir is None

        super(PlotSim1d, self).__init__(1 + int(self.show_difference),
                                        self.sims.num_species,
                                        (1, self.sims.stochastic[0, 0].shape[0]),
                                        file_name,
                                        save_dir,
                                        "")

        self.set_x_lim(self.sims.domain_bounds)
        self.set_x_label("x")
        self.set_y_label(self.sims.data_type)

        for species in range(self.sims.num_species):
            self.set_title("Simulation Results - species {}".format(self.sims.indices[species]), ax_index=(species, 0))
            self.set_y_lim((0.9 * self.sims.min_value[species], 1.1 * self.sims.max_value[species]), (species, 0))
            self.set_labels(["stochastic", "analytic"], (species, 0))

            if self.show_difference:
                self.set_title("Difference", ax_index=(species, 1))
                self.set_y_lim((-0.1 * self.sims.max_value[species], 0.1 * self.sims.max_value[species]), (species, 1))
                self.set_labels("difference", (species, 1))

        # Initialise the rest of the attributes
        self.interval = 10
        self.representation = "bar"
        self.stochastic_colour = "blue"
        self.analytic_colour = "green"
        self.x = self.sims.domain[0]

    def set_domain(self, array):
        assert isinstance(array, (list, tuple, np.ndarray))
        self.x = np.array(array)

    def set_interval(self, value):
        assert isinstance(value, int)
        self.interval = value

    def set_stochastic_colour(self, colour):
        assert isinstance(colour, str)
        self.stochastic_colour = colour

    def set_analytic_colour(self, colour):
        assert isinstance(colour, str)
        self.analytic_colour = colour

    def show(self, dpi=None, file_type="pdf"):

        stochastic_lines = []
        analytic_lines = []
        difference_lines = []

        # Calculate the spacing in the x-axis data
        if len(self.x) == 1:
            dx = self.x[0]
        else:
            dx = self.x[1] - self.x[0]
        string = get_app_string(self.sims.time_step)

        for species in range(self.sims.num_species):

            # Set up the line, bar or step objects that will have data assigned to them at each time step
            if self.representation == "line":
                stochastic_lines.append(
                    self.axes[species, 0].plot(self.x, self.sims.stochastic[species, 0],
                                               lw=2.0, color=self.stochastic_colour,
                                               label=self.labels[species, 0][0])[0])

            elif self.representation == "step":
                stochastic_lines.append(
                    self.axes[species, 0].step(self.x + 0.5 * dx, self.sims.stochastic[species, 0],
                                               lw=2.0, color=self.stochastic_colour,
                                               label=self.labels[species, 0][0])[0])

            elif self.representation == "bar":
                stochastic_lines.append(
                    self.axes[species, 0].bar(self.x, self.sims.stochastic[species, 0],
                                              width=dx, color=self.stochastic_colour,
                                              label=self.labels[species, 0][0]))

            if self.show_analytic:
                analytic_lines.append(
                    self.axes[species, 0].plot(self.x, self.sims.analytic[species, 0],
                                               lw=2.0, color=self.analytic_colour,
                                               label=self.labels[species, 0][1])[0])

            # There is an option to show the difference between the stochastic and analytic solutions
            if self.show_difference:
                diff_0 = self.sims.analytic[species, 0] - self.sims.stochastic[species, 0]
                difference_lines.append(self.axes[species, 1].plot(self.x, diff_0, lw=2.0, color="black")[0])

        # Set the layout to be tight, as it makes the figures look better.
        plt.tight_layout()
        self.legend()

        if self.save_dir is not None:
            self.sims.save_stochastic(self.file_name, self.save_dir)
            if self.show_analytic or self.show_difference:
                self.sims.save_analytic(self.file_name, self.save_dir)

        if self.sims.num_time_steps == 1:
            self.change_string(string.format(self.sims.end_time))
            save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi)

        elif self.sims.num_time_steps > 1:

            # Animation function.  This is called sequentially.
            def animate(i):
                self.change_string(string.format(self.sims.t_values[i]))

                for j in range(self.sims.num_species):
                    if self.representation == "line":
                        stochastic_lines[j].set_data(self.x, self.sims.stochastic[j, i])
                    elif self.representation == "step":
                        stochastic_lines[j].set_data(self.x + 0.5 * dx, self.sims.stochastic[j, i])
                    elif self.representation == "bar":
                        for rect, yi in zip(stochastic_lines[j], self.sims.stochastic[j, i]):
                            rect.set_height(yi)

                    if self.show_analytic and not self.show_difference:
                        analytic_lines[j].set_data(self.x, self.sims.analytic[j, i])

                    elif self.show_difference and not self.show_analytic:
                        diff = self.sims.analytic[j, i] - self.sims.stochastic[j, i]
                        difference_lines[j].set_data(self.x, diff)

                    elif self.show_analytic and self.show_difference:
                        diff = self.sims.analytic[j, i] - self.sims.stochastic[j, i]
                        analytic_lines[j].set_data(self.x, self.sims.analytic[j, i])
                        difference_lines[j].set_data(self.x, diff)

                if self.show_analytic and not self.show_difference:
                    return stochastic_lines, analytic_lines, self.text
                elif self.show_difference and not self.show_analytic:
                    return stochastic_lines, difference_lines, self.text
                elif self.show_analytic and self.show_difference:
                    return stochastic_lines, analytic_lines, difference_lines, self.text
                else:
                    return stochastic_lines, self.text

            # Call the animator. blit=True means only re-draw the parts that have changed.
            anim = animation.FuncAnimation(self.fig,
                                           animate,
                                           frames=self.sims.num_time_steps,
                                           interval=self.interval,
                                           repeat=False)

            save_or_show(anim, self.file_name, self.save_dir, file_type, dpi)

        else:
            print("\nERROR: Parameter num_time_steps is not valid! "
                  "(num_time_steps = {})".format(self.sims.num_time_steps))


class PlotSim2d(Axes2d):
    def __init__(self,
                 sims,
                 file_name="sim",
                 save_dir=None,
                 show_analytic=False,
                 show_difference=False):

        # First check the input and assign them to member attributes
        assert isinstance(sims, Data)
        self.sims = sims
        assert sims.num_dims == 2
        assert isinstance(show_analytic, bool)
        self.show_analytic = show_analytic
        assert isinstance(show_difference, bool)
        self.show_difference = show_difference

        assert isinstance(file_name, str)
        assert isinstance(save_dir, str) or save_dir is None

        super(PlotSim2d, self).__init__(1 + int(self.show_analytic) + int(self.show_difference),
                                        self.sims.num_species,
                                        self.sims.stochastic[0, 0].shape,
                                        file_name,
                                        save_dir,
                                        "")

        self.set_x_lim(self.sims.domain_bounds[0:2])
        self.set_y_lim(self.sims.domain_bounds[2:4])
        self.set_x_label("x")
        self.set_y_label("y")

        for species in range(self.sims.num_species):
            self.set_title("Simulation Results - species {}".format(self.sims.indices[species]), ax_index=(species, 0))

            if show_analytic and not show_difference:
                self.set_title("Analytic", ax_index=(species, 1))

            elif show_difference and not show_analytic:
                self.set_title("Difference", ax_index=(species, 1))

            elif show_difference and show_analytic:
                self.set_title("Analytic", ax_index=(species, 1))
                self.set_title("Difference", ax_index=(species, 2))

        # Initialise the rest of the attributes
        self.interval = 10
        self.aspect = "auto"
        self.stochastic_colormap = "Greys"
        self.analytic_colormap = "Greys"
        self.diff_colormap = "RdBu"
        self.x = self.sims.domain

        self.set_extent(self.sims.domain_bounds)
        self.set_labels(self.sims.data_type)

    def set_domain(self, array):
        assert isinstance(array, (list, tuple, np.ndarray))
        self.x = array

    def set_interval(self, value):
        assert isinstance(value, int)
        self.interval = value

    def set_aspect(self, value):
        assert isinstance(value, (str, int, float))
        if isinstance(value, str):
            assert value in ("auto", "equal")
        self.aspect = value

    def set_stochastic_colourmap(self, colormap):
        assert isinstance(colormap, str)
        self.stochastic_colormap = colormap

    def set_analytic_colourmap(self, colormap):
        assert isinstance(colormap, str)
        self.analytic_colormap = colormap

    def set_difference_colourmap(self, colormap):
        assert isinstance(colormap, str)
        self.diff_colormap = colormap

    def show(self, dpi=None, file_type="pdf"):

        stochastic_images = []
        analytic_images = []
        difference_images = []

        string = get_app_string(self.sims.time_step)

        for species in range(self.sims.num_species):
            stochastic_images.append(
                self.axes[species, 0].imshow(X=self.sims.stochastic[species, 0],
                                             interpolation='nearest', origin='lower',
                                             extent=self.extent[species, 0],
                                             vmin=0.9 * self.sims.min_value[species],
                                             vmax=1.1 * self.sims.max_value[species],
                                             aspect=self.aspect))
            stochastic_images[-1].set_cmap(self.stochastic_colormap)
            cbr_stochastic = self.fig.colorbar(stochastic_images[-1], ax=self.axes[species, 0])
            cbr_stochastic.set_label(self.labels[species, 0], fontsize=self.cbar_font_size[species, 0])
            cbr_stochastic.ax.tick_params(labelsize=self.cbar_tick_size[species, 0])

            if self.show_analytic and not self.show_difference:
                analytic_images.append(
                    self.axes[species, 1].imshow(X=self.sims.analytic[species, 0],
                                                 interpolation='nearest', origin='lower',
                                                 extent=self.extent[species, 1],
                                                 vmin=0.9 * self.sims.min_value[species],
                                                 vmax=1.1 * self.sims.max_value[species],
                                                 aspect=self.aspect))
                analytic_images[-1].set_cmap(self.analytic_colormap)
                cbr_analytic = self.fig.colorbar(analytic_images[-1], ax=self.axes[species, 1])
                cbr_analytic.set_label(self.labels[species, 1], fontsize=self.cbar_font_size[species, 1])
                cbr_analytic.ax.tick_params(labelsize=self.cbar_tick_size[species, 1])

            elif self.show_difference and not self.show_analytic:
                difference_images.append(
                    self.axes[species, 1].imshow(X=0.01 * self.sims.stochastic[species, 0],
                                                 interpolation='nearest', origin='lower',
                                                 extent=self.extent[species, 1],
                                                 vmin=-0.1 * self.sims.max_value[species],
                                                 vmax=0.1 * self.sims.max_value[species],
                                                 aspect=self.aspect))
                difference_images[-1].set_cmap(self.diff_colormap)
                cbr_diff = self.fig.colorbar(difference_images[-1], ax=self.axes[species, 1])
                cbr_diff.set_label(self.labels[species, 1], fontsize=self.cbar_font_size[species, 1])
                cbr_diff.ax.tick_params(labelsize=self.cbar_tick_size[species, 1])

            elif self.show_difference and self.show_analytic:
                analytic_images.append(
                    self.axes[species, 1].imshow(X=self.sims.analytic[species, 0],
                                                 interpolation='nearest', origin='lower',
                                                 extent=self.extent[species, 1],
                                                 vmin=0.9 * self.sims.min_value[species],
                                                 vmax=1.1 * self.sims.max_value[species],
                                                 aspect=self.aspect))
                analytic_images[-1].set_cmap(self.analytic_colormap)
                cbr_analytic = self.fig.colorbar(analytic_images[-1], ax=self.axes[species, 1])
                cbr_analytic.set_label(self.labels[species, 1], fontsize=self.cbar_font_size[species, 1])
                cbr_analytic.ax.tick_params(labelsize=self.cbar_tick_size[species, 1])

                diff_0 = self.sims.analytic[species, 0] - self.sims.stochastic[species, 0]
                difference_images.append(
                    self.axes[species, 2].imshow(X=diff_0,
                                                 interpolation='nearest', origin='lower',
                                                 extent=self.extent[species, 2],
                                                 vmin=-0.1 * self.sims.max_value[species],
                                                 vmax=0.1 * self.sims.max_value[species],
                                                 aspect=self.aspect))
                difference_images[-1].set_cmap(self.diff_colormap)
                cbr_diff = self.fig.colorbar(difference_images[-1], ax=self.axes[species, 2])
                cbr_diff.set_label(self.labels[species, 2], fontsize=self.cbar_font_size[species, 2])
                cbr_diff.ax.tick_params(labelsize=self.cbar_tick_size[species, 2])

        self.__set_fig_size__(1.25 * self.size, self.size)

        # Set the layout to be tight, as it makes the figures look better.
        plt.tight_layout()

        if self.save_dir is not None:
            self.sims.save_stochastic(self.file_name, self.save_dir)
            if self.show_analytic or self.show_difference:
                self.sims.save_analytic(self.file_name, self.save_dir)

        if self.sims.num_time_steps == 1:
            self.change_string(string.format(self.sims.end_time))
            save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi)

        elif self.sims.num_time_steps > 1:
            # Animation function.  This is called sequentially.
            def animate(i):
                self.change_string(string.format(self.sims.t_values[i]))

                for j in range(self.sims.num_species):
                    stochastic_images[j].set_array(self.sims.stochastic[j, i])

                if self.show_analytic and not self.show_difference:
                    for j in range(self.sims.num_species):
                        analytic_images[j].set_array(self.sims.analytic[j, i])

                elif self.show_difference and not self.show_analytic:
                    for j in range(self.sims.num_species):
                        diff = self.sims.analytic[j, i] - self.sims.stochastic[j, i]
                        difference_images[j].set_array(diff)

                elif self.show_analytic and self.show_difference:
                    for j in range(self.sims.num_species):
                        diff = self.sims.analytic[j, i] - self.sims.stochastic[j, i]
                        analytic_images[j].set_data(self.sims.analytic[j, i])
                        difference_images[j].set_data(diff)

                if self.show_analytic and not self.show_difference:
                    return stochastic_images, analytic_images, self.text
                elif self.show_difference and not self.show_analytic:
                    return stochastic_images, difference_images, self.text
                elif self.show_analytic and self.show_difference:
                    return stochastic_images, analytic_images, difference_images, self.text
                else:
                    return stochastic_images, self.text

            # Call the animator. blit=True means only re-draw the parts that have changed.
            anim = animation.FuncAnimation(self.fig,
                                           animate,
                                           frames=self.sims.num_time_steps,
                                           interval=self.interval,
                                           repeat=False)

            save_or_show(anim, self.file_name, self.save_dir, file_type, dpi)

        else:
            print("\nERROR: Parameter num_time_steps is not valid! "
                  "(num_time_steps = {})".format(self.sims.num_time_steps))


class Plot1d(Axes2d):

    def __init__(self, x, y=None, file_name="Figure", save_dir=None, rep="line", string=None):

        assert isinstance(rep, str) or rep is None
        assert isinstance(x, (list, tuple, np.ndarray))
        assert isinstance(y, (list, tuple, np.ndarray)) or y is None
        assert isinstance(file_name, str)
        assert isinstance(save_dir, str) or save_dir is None

        x = np.array(x)
        if y is None:
            y = x
        else:
            y = np.array(y)
            message = "Sizes:{}, {}".format(x.shape, y.shape)
            assert x.shape == np.array(y).shape, message

        assert x.ndim >= 1

        if x.ndim == 1:
            x = [[x]]
            y = [[y]]
        elif x.ndim == 2:
            x = [x]
            y = [y]

        self.x = np.array(x)
        self.y = np.array(y)
        self.num_axes = self.x.shape[0]
        self.num_arr = self.x.shape[1]
        self.rep = rep.lower()
        self.array_len = self.x[0, 0].shape[0]

        super(Plot1d, self).__init__(self.num_axes, 1,
                                     (1, self.array_len),
                                     file_name,
                                     save_dir,
                                     string)

        self.colours = np.array([[None]*self.num_arr for _ in range(self.num_axes)])
        self.ls = np.array([["-"]*self.num_arr for _ in range(self.num_axes)])
        for i in range(self.num_axes_y):
            for j in range(self.num_axes_x):
                self.labels[i, j] = [None]*self.num_arr

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_linestyle(self, ls, ax_index=None):
        if isinstance(ls, str):
            ls = [ls] * self.num_arr
        assert isinstance(ls, (list, tuple, np.ndarray))
        assert len(ls) == self.num_arr
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.ls[tuple(ax_index)] = ls
        else:
            self.ls = np.array([ls for _ in range(self.num_axes)])

    def set_colour(self, colour, ax_index=None):
        if isinstance(colour, str):
            colour = [colour] * self.num_arr
        assert isinstance(colour, (list, tuple, np.ndarray))
        assert len(colour) == self.num_arr
        if isinstance(ax_index, (tuple, list, np.ndarray)):
            self.colours[tuple(ax_index)] = colour
        else:
            self.colours = np.array([colour for _ in range(self.num_axes)])

    def show(self, dpi=None, file_type="pdf"):
        for ax in range(self.num_axes):
            assert len(self.labels[0, ax]) == self.num_arr
            for arr in range(self.num_arr):
                if self.rep == "line":
                    self.axes[0, ax].plot(self.x[ax, arr], self.y[ax, arr],
                                          ls=self.ls[ax, arr],
                                          color=self.colours[ax, arr],
                                          label=self.labels[0, ax][arr])
                elif self.rep == "bar":
                    if len(self.x[ax, arr]) == 1:
                        dx = self.x[ax, arr][0]
                    else:
                        dx = self.x[ax, arr][1] - self.x[ax, arr][0]
                    self.axes[0, ax].bar(self.x[ax, arr], self.y[ax, arr],
                                         width=dx,
                                         color=self.colours[ax, arr],
                                         label=self.labels[0, ax][arr])
                elif self.rep == "step":
                    self.axes[0, ax].step(self.x[ax, arr], self.y[ax, arr],
                                          color=self.colours[ax, arr],
                                          label=self.labels[0, ax][arr])
                elif self.rep == "hist":
                    self.axes[0, ax].hist(self.x[ax, arr],
                                          color=self.colours[ax, arr],
                                          label=self.labels[0, ax][arr])

        self.legend()

        data_plotted = np.concatenate([np.concatenate(self.x), np.concatenate(self.y)])
        if len(self.data):
            data_plotted = [data_plotted]
            for i in range(len(self.data)):
                data_plotted.append(self.data[i])
        original_data_shape = (2, self.x.shape[0], self.x.shape[1], self.x.shape[2])
        hdr = "original_data_shape = {}".format(str(original_data_shape).strip("()").replace(" ", ""))

        plt.tight_layout()

        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi, data_plotted, hdr)


class Plot2d(Axes2d):

    def __init__(self, x, file_name="Figure", save_dir=None, string=None):

        assert isinstance(x, (tuple, list, np.ndarray))

        x = np.array(x)
        assert x.ndim >= 2

        if x.ndim == 2:
            self.x = np.array([[x]])
        elif x.ndim == 3:
            self.x = np.array([x])
        else:
            self.x = np.array(x)

        super(Plot2d, self).__init__(self.x.shape[0],
                                     self.x.shape[1],
                                     self.x[0, 0].shape,
                                     file_name,
                                     save_dir,
                                     string)

        self.colormap = "Greys"
        self.interpolation = "nearest"
        self.aspect = "auto"

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_colormap(self, colormap):
        assert isinstance(colormap, str)
        self.colormap = colormap

    def set_aspect(self, value):
        assert isinstance(value, (str, int, float))
        if isinstance(value, str):
            assert value in ("auto", "equal")
        self.aspect = value

    def set_interpolation(self, name):
        assert isinstance(name, str)
        self.interpolation = name

    def show(self, dpi=None, file_type="pdf"):
        padding = 0
        for i in range(self.num_axes_y):
            for j in range(self.num_axes_x):
                im = self.axes[i, j].imshow(self.x[j, i],
                                            origin="lower",
                                            aspect=self.aspect,
                                            interpolation=self.interpolation,
                                            extent=self.extent[i, j])

                if self.labels[i, j] is not None:
                    cbr = self.fig.colorbar(im, ax=self.axes[i, j])
                    cbr.set_label(self.labels[i, j], fontsize=self.cbar_font_size[i, j])
                    cbr.ax.tick_params(labelsize=self.cbar_tick_size[i, j])

                im.set_cmap(self.colormap)

            if not all(v is None for v in self.labels[i, :]):
                padding += (0.25 * self.size)

        # Update the figure size, if colorbars were added
        self.__set_fig_size__((self.size + padding), self.size)

        data_plotted = np.concatenate(np.concatenate(self.x))
        if len(self.data):
            data_plotted = [data_plotted]
            for i in range(len(self.data)):
                data_plotted.append(self.data[i])
        original_data_shape = list(self.x.shape)
        hdr = "original_data_shape = {}".format(str(original_data_shape).strip("()").replace(" ", ""))

        plt.tight_layout()

        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi, data_plotted, hdr)


class Plot3dLines(Axes3d):

    def __init__(self, x, y, z, file_name="Figure", save_dir=None):
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
        assert isinstance(file_name, str)
        assert isinstance(save_dir, str) or save_dir is None

        self.x = np.ones(len(y)).reshape([len(y), 1]) * x
        self.y = y.reshape([len(y), 1])*np.ones(len(x))
        self.z = z
        assert self.x.shape == self.y.shape

        super(Plot3dLines, self).__init__(self.x.shape, file_name, save_dir)

        self.colours = ["black" for _ in range(len(self.x))]

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_colour(self, colour):
        assert isinstance(colour, str) or isinstance(colour, (list, tuple))
        if isinstance(colour, str):
            self.colours = [colour for _ in range(len(self.x))]

        if isinstance(colour, (list, tuple)):
            assert len(colour) == len(self.x)
            self.colours = colour

    def show(self, dpi=None, file_type="pdf"):
        for i in range(len(self.x)):
            self.ax.plot(self.x[i], self.y[i], self.z[i], color=self.colours[i])

        data_plotted = np.concatenate([self.x, self.y, self.z])
        original_data_shape = (3, self.x.shape[0], self.x.shape[1])
        hdr = "original_data_shape = {}".format(str(original_data_shape).strip("()").replace(" ", ""))

        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi, data_plotted, hdr)


class Plot3d(Axes3d):

    def __init__(self, x, y, z, file_name="Figure", save_dir=None):
        assert isinstance(x, np.ndarray)
        assert x.ndim == 1
        assert isinstance(y, np.ndarray)
        assert y.ndim == 1
        assert isinstance(z, np.ndarray)
        assert isinstance(file_name, str)
        assert isinstance(save_dir, str) or save_dir is None

        self.x, self.y = np.meshgrid(x, y)
        self.z = z

        message = "Sizes don't match: {}, {}, {}".format(self.x.shape, self.y.shape, self.z.shape)
        assert self.x.shape == self.y.shape == self.z.shape, message

        super(Plot3d, self).__init__(self.x.shape, file_name, save_dir)

        self.colourmap = "Greys"
        self.cbar_label = None
        self.cbar_aspect = None
        self.cbar_shrink = None

    def __del__(self):
        return "{} deleted".format(self.__class__.__name__)

    def set_colourmap(self, colourmap):
        assert isinstance(colourmap, str)
        self.colourmap = colourmap

    def colorbar(self, label="", shrink=0.5, aspect=10):
        assert isinstance(label, str)
        assert isinstance(shrink, (float, int))
        assert isinstance(aspect, (float, int))

        self.cbar_label = label
        self.cbar_shrink = shrink
        self.cbar_aspect = aspect

    def show(self, dpi=None, file_type="pdf"):
        surf = self.ax.plot_surface(self.x, self.y, self.z, cmap=self.colourmap, linewidth=self.lw, antialiased=False)

        if self.cbar_label is not None:
            cbr = self.fig.colorbar(surf, shrink=self.cbar_shrink, aspect=self.cbar_aspect, pad=-0.05)
            cbr.set_label(self.cbar_label, fontsize=self.cbar_font_size)
            cbr.ax.tick_params(labelsize=self.cbar_tick_size)
            # Update the figure size, if colorbars were added
            self.__set_fig_size__(1.1*self.size, self.size)

        data_plotted = np.concatenate([self.x, self.y, self.z])
        original_data_shape = (3, self.x.shape[0], self.x.shape[1])
        hdr = "original_data_shape = {}".format(str(original_data_shape).strip("()").replace(" ", ""))

        save_or_show(self.fig, self.file_name, self.save_dir, file_type, dpi, data_plotted, hdr)
