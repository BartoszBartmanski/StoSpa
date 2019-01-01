#!/usr/bin/python3

import numpy as np
import os
import sys
from matplotlib import tri
from scipy import fftpack as fft
from scipy import ndimage


def analytic_solution(num_dims,
                      t_val,
                      x_val=None,
                      domain_bounds=(0.0, 1.0),
                      x_0=(0.5, 0.5),
                      d=1.0,
                      k_decay=0.0,
                      k_influx=0.0,
                      trunc_order=100,
                      num_points=None):
    """This function returns the analytic solution to the heat equation with decay i.e. du/dt = nabla^2 u + k_1 - k_2 u
    k_1 is the production rate, k_2 is the decay rate
    Returns x-axis values, followed by an array of the solutions at different time points"""
    if isinstance(t_val, (int, float)):
        t_val = np.array([t_val])
    if isinstance(num_points, (int, float)):
        num_points = [num_points, num_points]
    if isinstance(x_0, (int, float)):
        x_0 = np.array([x_0, x_0])
    if len(domain_bounds) < 4:
        domain_bounds = (domain_bounds[0], domain_bounds[1], domain_bounds[0], domain_bounds[1])

    assert isinstance(t_val, (list, tuple, np.ndarray))
    assert isinstance(x_val, (tuple, list, np.ndarray)) or x_val is None
    assert isinstance(domain_bounds, (list, tuple, np.ndarray))
    assert isinstance(x_0, (tuple, list, np.ndarray))
    assert isinstance(d, (int, float))
    assert isinstance(k_decay, (int, float))
    assert isinstance(k_influx, (int, float))
    assert isinstance(trunc_order, int)

    length = float(domain_bounds[1] - domain_bounds[0])
    t = np.array(t_val)
    if x_val is None:
        assert num_points is not None
        x_val = [np.linspace(domain_bounds[0], domain_bounds[1], num_points[0]),
                 np.linspace(domain_bounds[0], domain_bounds[1], num_points[1])]

    if num_dims == 1:
        if isinstance(x_val[0], (tuple, list, np.ndarray)):
            x = np.array(x_val[0])
            y = np.array(x_val[0])
        else:
            x = np.array(x_val)
            y = np.array(x_val)

        assert t.ndim == 1
        t = t.reshape([t.shape[0], 1])

        u = 1.0 / length

        for n in range(1, trunc_order):
            u += (2/length)*np.cos((n*np.pi/length)*x_0[0])*np.cos((n*np.pi/length)*x)*np.exp(-d*(n*np.pi/length)**2*t)

    else:
        assert isinstance(x_val[0], (tuple, list, np.ndarray))
        assert isinstance(x_val[1], (tuple, list, np.ndarray))
        x = np.array(x_val[0])
        y = np.array(x_val[1])

        xx, yy = np.meshgrid(x, y)

        assert t.ndim == 1
        t = t.reshape([t.shape[0], 1, 1])

        u = 1.0 / length ** 2

        for k in range(1, trunc_order):
            u += (2.0 / length ** 2) * np.cos(k * np.pi * x_0[1] / length) * np.cos(k * np.pi * yy / length) * np.exp(
                -d * t * (k * np.pi / length) ** 2)
        for j in range(1, trunc_order):
            u += (2.0 / length ** 2) * np.cos(j * np.pi * x_0[0] / length) * np.cos(j * np.pi * xx / length) * np.exp(
                -d * t * (j * np.pi / length) ** 2)
        for j in range(1, trunc_order):
            for k in range(1, trunc_order):
                u += (4.0 / length ** 2) * np.cos(j * np.pi * x_0[0] / length) * np.cos(k * np.pi * x_0[1] / length) * \
                     np.cos(j * np.pi * xx / length) * np.cos(k * np.pi * yy / length) * \
                     np.exp(-d * t * ((j * np.pi / length) ** 2 + (k * np.pi / length) ** 2))

    if k_decay > 0.0 and k_influx == 0.0:
        u *= np.exp(- k_decay * t)
    elif k_decay == 0.0 and k_influx > 0.0:
        u += k_influx * t
    elif k_decay > 0.0 and k_influx > 0.0:
        u += k_influx * (1.0 - np.exp(-k_decay * t)) / k_decay

    if num_dims == 1:
        return u, x
    else:
        return u, x, y


def power_spec(a, factor=1.0, dim=1):
    if dim == 1:
        a = np.sum(np.absolute(factor * fft.dct(a, type=1, axis=-1)) ** 2, axis=0)

    elif dim == 2:  # for 2d
        a = np.sum(np.absolute(factor * fft.dct(fft.dct(a, type=1, axis=-2), type=1, axis=-1)) ** 2, axis=0)

    return a


def get_data(file_name):
    if file_name is None:
        return None

    desc = ""
    data_array = np.array([])
    data_shape = []
    try:
        data_array = np.loadtxt(file_name)
        f = open(file_name, "r")
    except IOError:
        print("Could not find the file {}".format(file_name))
    else:

        for line in f:
            if "#" not in line:
                break
            elif "original_data_shape" in line:
                data_shape = [int(x) for x in line.split()[-1].split(",")]
            else:
                desc += line
        f.close()
        print(desc)

    if len(data_shape) > 0:
        data_array = data_array.reshape(data_shape)
    else:
        print("Could not find the original dimensions of this data!\n")

    return data_array


class Sim(object):
    """Class to open the parameters file and the data files associated with a specific simulation."""

    def __init__(self, file_name):

        assert isinstance(file_name, str)

        try:
            f = open(file_name, "r")
        except IOError:
            print("Could not find the file: {}\n".format(file_name))
            quit()
        else:
            # Assign all the parameters
            self.file_name = file_name
            self.additional_reactions = {}
            self.info = ""

            # Declare attribues some attributes
            self.h = None
            self.kappa = None
            self.alpha = None
            self.beta = None
            self.data_type = "molecules"

            sys.stdout.write("\r{}".format("Opening file: " + self.file_name))
            sys.stdout.flush()

            for line in f:

                if "#" not in line:
                    break

                elif "command" in line.lower():
                    pass

                elif "num_runs" in line:
                    self.num_runs = int(line.split()[-1])

                elif "num_species" in line:
                    self.num_species = int(line.split()[-1])

                elif "num_method" in line:
                    self.num_method = line.split()[-1].lower()

                elif "num_dims" in line:
                    self.num_dims = int(line.split()[-1])

                elif "num_voxels" in line:
                    self.num_voxels = int(line.split()[-1])

                elif "domain_bounds" in line:
                    self.domain_bounds = np.array(line.split()[-1].split(",")).astype(float)

                elif " h " in line:
                    self.h = float(line.split()[-1])

                elif "kappa" in line:
                    self.kappa = float(line.split()[-1])

                elif "time_step" in line:
                    self.time_step = float(line.split()[-1])

                elif "end_time" in line:
                    self.end_time = float(line.split()[-1])

                elif "data_type" in line:
                    self.data_type = line.split()[-1]

                elif "bc" in line:
                    self.bc = line.split()[-1]

                elif "alpha" in line:
                    self.alpha = float(line.split()[-1])

                elif "beta" in line:
                    self.beta = np.array(line.split()[-1].split(",")).astype(float)

                elif "diffusion_constants" in line:
                    self.diffusion_constants = np.array(line.split()[-1].split(",")).astype(float)

                elif "decay_constants" in line:
                    self.decay_constants = np.array(line.split()[-1].split(",")).astype(float)

                elif "prod_constants" in line:
                    self.prod_constants = np.array(line.split()[-1].split(",")).astype(float)

                elif "\t" in line:
                    self.additional_reactions[line.split()[1]] = float(line.split()[-1])

                self.info += line.replace("# ", "").replace("#\t ", "\t")

        if not hasattr(self, "diffusion_constants"):
            self.diffusion_constants = np.zeros(self.num_species)
        if not hasattr(self, "decay_constants"):
            self.decay_constants = np.zeros(self.num_species)
        if not hasattr(self, "prod_constants"):
            self.prod_constants = np.zeros(self.num_species)

        if self.kappa is not None:
            assert self.num_voxels * self.kappa == np.floor(self.num_voxels * self.kappa)  # check kappa is correct

        # Get the stochastic simulation data
        stochastic = np.loadtxt(self.file_name)
        self.num_time_steps = int(stochastic.shape[0]/self.num_species)
        if self.num_dims == 1:
            self.num_voxels = [self.num_voxels, self.num_voxels]
        else:  # self.num_dims == 2
            self.num_voxels = [self.num_voxels, int(self.num_voxels * self.kappa)]
            stochastic = stochastic.reshape([stochastic.shape[0], self.num_voxels[1], self.num_voxels[0]])

        self.stochastic = []
        for species in range(self.num_species):
            self.stochastic.append(stochastic[species::self.num_species])

        # At this point the data is organised as follows
        # 0th index - species
        # 1st index - time point
        # rest - spatial indices
        new_shape = [self.num_species, self.num_time_steps] + [self.num_voxels[1], self.num_voxels[0]][:self.num_dims]
        self.stochastic = np.array(self.stochastic).reshape(new_shape)

        if self.h is None:
            self.h = (self.domain_bounds[1] - self.domain_bounds[0]) / self.num_voxels[-1]

    def __eq__(self, other):
        ignored = ("info", "stochastic", "file_name")
        for k1 in self.__dict__:
            if k1 not in ignored and \
                    (k1 not in other.__dict__ or not np.array_equal(other.__dict__[k1], self.__dict__[k1])):
                return False
        for k2 in other.__dict__:
            if k2 not in ignored and k2 not in self.__dict__:
                return False
        return True

    def diff(self, other):
        difference = []
        ignored = ("info", "stochastic", "file_name")
        for k1 in self.__dict__:
            if k1 not in ignored and \
                    (k1 not in other.__dict__ or not np.array_equal(other.__dict__[k1], self.__dict__[k1])):
                difference.append(k1)
        for k2 in other.__dict__:
            if k2 not in ignored and k2 not in self.__dict__:
                difference.append(k2)
        return difference

    def __str__(self):
        return self.info

    def species(self, index):
        assert isinstance(index, int)
        return self.stochastic[index % self.num_species]


class Data(object):
    """Class to store the data from multiple simulations
        0th index - simulation index
        1st index - species index
        2nd index - time index
        3rd/4th index - spatial index
    """

    def __init__(self, file_name):
        """Constructor method."""

        self.__stochastic__ = []
        self.sim_names = []

        # Check whether multiple files will be read (this is indicated by hashtag symbol at the end of file name)
        if "#.dat" in file_name:
            self.file_name = file_name.split("#.dat")[0]
            # Check how many files there are of the similar name
            files = os.listdir("./")
            for i in range(1000):
                specific_filename = self.file_name + str(i).zfill(3) + ".dat"
                if specific_filename in files and specific_filename not in self.sim_names:
                    self.sim_names.append(specific_filename)
        else:
            self.file_name = file_name.split(".dat")[0]
            self.sim_names.append(file_name)

        self.num_sims = len(self.sim_names)
        if self.num_sims == 0:
            print("No simulations with such a name ({})\n".format(file_name))
            quit()

        # First get one simulation and then compare all the parameters to this
        sim_0 = Sim(self.sim_names[0])
        for name in self.sim_names:
            next_sim = Sim(name)
            if next_sim == sim_0:
                self.__stochastic__.append(next_sim.stochastic)
            else:
                message = "\n{} != {} => The following attributes don't agree: {}"
                differences = str(sim_0.diff(next_sim)).strip("[]").replace("'", "")
                print(message.format(next_sim.file_name, sim_0.file_name, differences))
        sys.stdout.write("\n")

        # Keep only one set of parameters
        self.num_species = sim_0.num_species
        self.num_method = sim_0.num_method
        self.num_dims = sim_0.num_dims
        self.num_voxels = sim_0.num_voxels
        self.h = sim_0.h
        self.kappa = sim_0.kappa
        self.domain_bounds = sim_0.domain_bounds
        self.alpha = sim_0.alpha
        self.beta = sim_0.beta
        self.data_type = sim_0.data_type
        self.time_step = sim_0.time_step
        self.num_time_steps = sim_0.num_time_steps
        self.end_time = sim_0.end_time
        self.bc = sim_0.bc
        self.info = sim_0.info

        self.diffusion_constants = sim_0.diffusion_constants
        self.decay_constants = sim_0.decay_constants
        self.prod_constants = sim_0.prod_constants
        self.additional_reactions = sim_0.additional_reactions

        # 0th index - simulation index
        # 1st index - species index
        # 2nd index - time index
        # 3rd/4th index - spatial index
        self.__stochastic__ = np.array(self.__stochastic__)
        # 0th index - species index
        # 1st index - time index
        # 2nd/3rd index - spatial index
        self.stochastic = self.__stochastic__[0]  # first sim stochastic data
        self.analytic = np.zeros(self.stochastic.shape)

        self.t_values = np.arange(0, self.end_time, self.time_step)
        self.total_num_molecules = np.sum(self.stochastic[:, 0], axis=(1, 2)[:self.num_dims])
        self.min_value = np.zeros(self.num_species)
        self.max_value = np.ones(self.num_species)
        self.indices = np.arange(self.num_species)

        self.__update_max_min__()

        if self.num_dims == 1:
            # For 1d data will need the x-axis data
            self.domain = [np.linspace(self.domain_bounds[0] + 0.5 * self.h, self.domain_bounds[1] - 0.5 * self.h,
                                       self.num_voxels[0])]

            self.freq = [np.arange(self.num_voxels[0])]

            self.domain_bounds = [self.domain_bounds[0], self.domain_bounds[1]]

            self.voxel_size = self.h

            self.x_0 = [self.domain[0][np.argmax(self.stochastic[0, 0])]]

            self.corners = [self.domain[0][0], self.domain[0][-1]]

            self.labels = [r"$x$", r"$t$"]

        else:
            # For 2d just provide two domain - along x-axis and along y-axis
            domain_x = np.linspace(self.domain_bounds[0] + 0.5 * self.kappa * self.h,
                                   self.domain_bounds[1] - 0.5 * self.kappa * self.h,
                                   self.num_voxels[0])
            domain_y = np.linspace(self.domain_bounds[0] + 0.5 * self.h,
                                   self.domain_bounds[1] - 0.5 * self.h,
                                   self.num_voxels[1])

            self.domain = [domain_x, domain_y]

            self.freq = [np.arange(self.num_voxels[0]), np.arange(self.num_voxels[1])]

            self.domain_bounds = [self.domain_bounds[0], self.domain_bounds[1],
                                  self.domain_bounds[0], self.domain_bounds[1]]

            self.voxel_size = self.kappa * self.h ** 2

            idx = tuple(np.unravel_index(np.argmax(self.stochastic[0, 0]), self.num_voxels))
            self.x_0 = [self.domain[0][idx[0]], self.domain[0][idx[0]]]

            self.corners = [self.domain[0][0], self.domain[0][-1], self.domain[1][0], self.domain[1][-1]]

            self.labels = [r"$x$", r"$y$"]

        self.x = self.domain  # Container for the domain or frequency

    def __repr__(self):
        return self.info

    def __str__(self):
        return self.info

    def __update_max_min__(self):
        self.min_value = np.amin(self.stochastic, axis=(1, 2, 3)[:1+self.num_dims])
        self.max_value = np.amax(self.stochastic, axis=(1, 2, 3)[:1+self.num_dims])

    def __save__(self, att_str, file_name=None, save_dir="./", append=None):
        assert att_str == "stochastic" or att_str == "analytic"
        assert isinstance(file_name, str) or file_name is None
        assert isinstance(append, str) or append is None
        assert os.path.isdir(save_dir), "Not a valid path {}".format(save_dir)
        if file_name is None:
            file_name = self.file_name

        data = getattr(self, att_str)
        save_data = np.reshape(data, np.prod(data.shape))

        header = "Data output from data.py module.\n"
        header += "original_data_shape = {}\n".format(str(data.shape).strip("()").replace(" ", ""))

        file_name += "_" + att_str
        if append:
            file_name += "_" + append
        path_to_file = save_dir + "/" + file_name
        path_to_file = path_to_file.replace("//", "/")

        print("Data saved as", file_name)
        with open(path_to_file + ".dat", "wb") as handle:
            np.savetxt(handle, save_data, header=str(header + self.info))

    def save_stochastic(self, file_name=None, save_dir="./", append=None):
        self.__save__("stochastic", file_name, save_dir, append)

    def save_analytic(self, file_name=None, save_dir="./", append=None):
        self.__save__("analytic", file_name, save_dir, append)

    def add_analytic(self, analytic_array=None):
        if analytic_array is not None:
            try:
                analytic_array = np.array(analytic_array).reshape(self.stochastic.shape)
            except ValueError:
                print("The shape of the analytic array does not match the shape of the stochastic array! "
                      "{} != {}".format(self.stochastic.shape, analytic_array.shape))

        else:
            analytic_array = []
            for idx in range(self.num_species):
                analytic_array.append(analytic_solution(self.num_dims,
                                                        self.t_values,
                                                        self.domain,
                                                        self.domain_bounds,
                                                        self.x_0,
                                                        self.diffusion_constants[idx],
                                                        self.decay_constants[idx],
                                                        self.prod_constants[idx] / self.total_num_molecules[idx])[0])
                if self.data_type == "molecules":
                    analytic_array[-1] *= (self.total_num_molecules[idx] * self.voxel_size)
            analytic_array = np.array(analytic_array)

        self.analytic = analytic_array

    def reset(self):
        """Resets the stochastic container to the original state."""
        self.stochastic = self.__stochastic__[0]
        self.__update_max_min__()

    def slice(self, t, x, species=None):
        """Takes an appropriate slice of the data both through the time series and the spatial data."""
        # Make sure that the list of bounds is of the appropriate length

        if not isinstance(t, list):
            t = []
        t += [0, self.num_time_steps + 1][len(t):2]

        if not isinstance(x, list):
            x = []
        x += [0, self.num_voxels[0], 0, self.num_voxels[1]][len(x):4]

        if isinstance(species, int):
            species = [species % (self.num_species + 1), species+1 % (self.num_species + 1)]
            self.num_species = 1
        else:
            species = [0, self.num_species + 1]

        for elem in t:
            assert isinstance(elem, int)
        for elem in x:
            assert isinstance(elem, int)

        # Correct indices of the bounds
        self.t_values = self.t_values[t[0]:t[1]]
        self.end_time = self.t_values[-1]
        self.num_time_steps = len(self.t_values)
        self.indices = self.indices[species[0]:species[1]]

        if self.num_dims == 1:
            # Take an appropriate slice
            self.domain[0] = self.domain[0][x[0]:x[1]]

            self.domain_bounds = [self.domain[0][0] - 0.5 * self.h, self.domain[0][-1] + 0.5 * self.h]
            self.freq[0] = self.freq[0][x[0]:x[1]]

            self.corners = [self.domain[0][0], self.domain[0][-1]]

            self.stochastic = self.stochastic[species[0]:species[1], t[0]:t[1], x[0]:x[1]]
            self.analytic = self.analytic[species[0]:species[1], t[0]:t[1], x[0]:x[1]]

        else:
            # Take an appropriate slice of the domain
            self.domain[0] = self.domain[0][x[0]:x[1]]
            self.domain[1] = self.domain[1][x[2]:x[3]]

            # Update the domain bounds
            self.domain_bounds = [self.domain[0][0] - 0.5 * self.kappa * self.h,
                                  self.domain[0][-1] + 0.5 * self.kappa * self.h,
                                  self.domain[1][0] - 0.5 * self.h,
                                  self.domain[1][-1] + 0.5 * self.h]

            # Take an appropriate slice of the frequency spectrum
            self.freq[0] = self.freq[0][x[0]:x[1]]
            self.freq[1] = self.freq[1][x[2]:x[3]]

            self.corners = [self.domain[0][0], self.domain[0][-1], self.domain[1][0], self.domain[1][-1]]

            self.stochastic = self.stochastic[species[0]:species[1], t[0]:t[1], x[2]:x[3], x[0]:x[1]]
            self.analytic = self.analytic[species[0]:species[1], t[0]:t[1], x[2]:x[3], x[0]:x[1]]

        self.__update_max_min__()

    def average(self):
        """Averages over all the simulations."""
        self.stochastic = np.sum(self.__stochastic__, axis=0) / self.num_sims
        self.__update_max_min__()

    def fft(self, remove_zero=False, analytic=False):
        """Performs the discrete Fourier transform."""

        if self.num_dims == 1:  # for 1d
            self.freq = np.fft.rfftfreq(self.num_voxels[0], 0.5 / self.num_voxels[0])
            self.stochastic = np.sum(np.absolute(np.fft.rfft(self.__stochastic__, axis=3, norm='ortho')) ** 2, axis=0)
            if remove_zero:
                self.stochastic[:, :, np.where(self.freq[0] == 0)] = 0

            if analytic:
                self.analytic = np.absolute(np.fft.rfft(self.analytic, axis=2, norm='ortho')) ** 2
                if remove_zero:
                    self.analytic[:, :, np.where(self.freq[0] == 0)] = 0

            self.labels = [r"$m$", r"$t$"]

        else:  # for 2d
            self.freq[0] = np.fft.rfftfreq(self.num_voxels[0], 0.5 / self.num_voxels[0])
            self.freq[1] = np.fft.rfftfreq(self.num_voxels[1], 0.5 / self.num_voxels[1])
            self.stochastic = np.sum(
                np.absolute(np.fft.rfft2(self.__stochastic__, axes=(3, 4), norm='ortho')) ** 2, axis=0)
            if remove_zero:
                self.stochastic[:, :, np.where(self.freq[1] == 0), np.where(self.freq[0] == 0)] = 0

            if analytic:
                self.analytic = np.absolute(np.fft.rfft2(self.analytic, axes=(2, 3), norm='ortho')) ** 2
                if remove_zero:
                    self.analytic[:, :, np.where(self.freq[1] == 0), np.where(self.freq[0] == 0)] = 0

            self.labels = [r"$m_x$", r"$m_y$"]

        self.__update_max_min__()

        self.data_type = r"$P_s$"
        self.x = self.freq

    def dct(self, remove_zero=False, analytic=False):
        """Performs discrete Fourier cosine transform."""

        self.stochastic = power_spec(self.__stochastic__, self.voxel_size, dim=self.num_dims)

        if analytic:
            self.analytic = power_spec(self.analytic, self.voxel_size, dim=self.num_dims)

        if self.num_dims == 1:  # for 1d

            if remove_zero:
                loc = np.where(self.freq[0] == 0)
                self.stochastic[:, :, loc] = 0
                self.analytic[:, :, loc] = 0

            self.labels = [r"$m$", r"$t$"]

        else:  # for 2d

            if remove_zero:
                loc = [np.where(self.freq[1] == 0), np.where(self.freq[0] == 0)]
                self.stochastic[:, :, loc[0], loc[1]] = 0
                self.analytic[:, :, loc[0], loc[1]] = 0

            self.labels = [r"$m_x$", r"$m_y$"]

        self.__update_max_min__()

        self.data_type = r"$P_s$"
        self.x = self.freq

    def find_dominant_modes(self):
        """Determines which wavemode is the dominant one by finding the maximum of the power spectrum."""

        if self.num_dims == 1:  # for 1d
            ps_s = np.absolute(fft.dct(self.__stochastic__, type=1, axis=3))**2
            ps_s[:, :, :, np.where(self.freq[0] == 0)] = 0

        else:  # for 2d
            ps_s = np.absolute(fft.dct(fft.dct(self.__stochastic__, type=1, axis=3), type=1, axis=4))**2
            ps_s[:, :, :, np.where(self.freq[1] == 0), np.where(self.freq[0] == 0)] = 0
            ps_s = ps_s.reshape(ps_s.shape[0], ps_s.shape[1], ps_s.shape[2], ps_s.shape[3] * ps_s.shape[4])

        wave_mode_freq = np.zeros(ps_s.shape)
        i, j, k = ps_s.shape[:3]
        id_0, id_1, id_2 = np.ogrid[:i, :j, :k]

        wave_mode_freq[id_0, id_1, id_2, np.argmax(ps_s, axis=3)] = 1
        wave_mode_freq = wave_mode_freq.reshape(self.__stochastic__.shape)

        return wave_mode_freq

    def count_modes(self):
        """Counts the number of simulations in which each wave mode is the dominant one
        (ignoring the zeroth mode)"""

        wave_mode_freq = self.find_dominant_modes()

        self.stochastic = np.sum(wave_mode_freq, axis=0)

        self.__update_max_min__()

        if self.num_dims == 1:
            self.labels = [r"$m$", r"$t$"]
        else:  # self.num_dimensions == 2
            self.labels = [r"$m_x$", r"$m_y$"]
        self.data_type = "Frequency"
        self.x = self.freq

    def get_sim_names(self, wavemode, species=0, time=-1):
        """Returns the simulation names with the specified dominant wavemode."""

        assert isinstance(species, int)
        if self.num_dims == 1:
            if isinstance(wavemode, (tuple, list, np.ndarray)):
                wavemode = wavemode[0]
            assert isinstance(wavemode, int)
            wave_mode_freq = self.find_dominant_modes()[:, species, time, wavemode]
        else:  # self.num_dimensions == 2
            if isinstance(wavemode, int):
                wavemode = (wavemode, wavemode)
            assert isinstance(wavemode, (tuple, list, np.ndarray))
            wave_mode_freq = self.find_dominant_modes()[:, species, time, wavemode[1], wavemode[0]]

        specific_sims = "\n".join(np.array(self.sim_names)[wave_mode_freq.astype(bool)])

        return specific_sims

    def find_peaks(self, size=None, cutoff=0.5, species=0, time=-1):
        """Returns the indices and the values of the local maxima of the data."""

        assert isinstance(species, int)
        assert isinstance(time, int)
        if isinstance(size, int):
            search = [(size, size)[:self.num_dims]] * self.num_sims
        elif isinstance(size, (tuple, list, np.ndarray)):
            assert isinstance(size, (list, tuple, np.ndarray)) and len(size) == 2
            search = [size[:self.num_dims]] * self.num_sims
        else:
            # Get the values of the dominant wavemodes to determine the search neighbourhood for each image
            wave_mode_freq = self.find_dominant_modes()[:, species, time]
            s = wave_mode_freq.shape
            if self.num_dims == 2:
                wave_mode_freq = wave_mode_freq.reshape(s[0], -1)

            # Get the dominant modes for each simulation
            m = np.sqrt(np.sum(np.array(np.unravel_index(np.argmax(wave_mode_freq, axis=-1), s[1:]))**2, axis=0))

            search_kappa = 0.5
            if self.num_dims == 1:
                search = search_kappa * (2 * self.domain_bounds[1]) / (m * self.h)
            else:  # self.num_dims == 2
                search = np.ones([self.num_sims, 2])
                search[:, 1] = search_kappa * (2 * self.domain_bounds[1]) / (m * self.kappa * self.h)
                search[:, 0] = search_kappa * (2 * self.domain_bounds[3]) / (m * self.h)

        indices_of_peaks = []
        for i in range(self.num_sims):
            filtered_data = ndimage.maximum_filter(self.__stochastic__[i, species, time], search[i], mode="constant")
            filtered_data[filtered_data < cutoff * np.max(filtered_data)] = 0
            indices_of_peaks.append(np.array(np.where(filtered_data == self.__stochastic__[i, species, time])))

        return indices_of_peaks

    def get_spot_properties(self, size=None, cutoff=0.5, species=0, time=-1):
        """Returns the mean spot area. For Turing patterns."""
        centers_of_spots = self.find_peaks(size, cutoff, species, time)
        indices_of_spots = self.find_peaks(1, cutoff, species, time)
        area_of_spots = []
        values_of_peaks = []
        dist_between_peaks = []
        if self.num_dims == 1:
            for i in range(self.num_sims):
                area_of_spots.append(len(indices_of_spots[i][0]) * self.voxel_size / len(centers_of_spots[i][0]))
                values_of_peaks += list(self.__stochastic__[i, species, time][centers_of_spots[i]][0])
                loc = self.domain[0][centers_of_spots[i]]
                dist_between_peaks.append(np.sum(np.diff(loc)) / len(np.diff(loc)))
        else:  # self.num_dims == 2
            for i in range(self.num_sims):
                area_of_spots.append(len(indices_of_spots[i][0]) * self.voxel_size / len(centers_of_spots[i][0]))
                sto = self.__stochastic__[i, species, time]
                values_of_peaks += list(sto.flatten()[np.ravel_multi_index(centers_of_spots[i], sto.shape)])
                loc = np.array([self.domain[0][centers_of_spots[i][1]], self.domain[1][centers_of_spots[i][0]]])
                triag = tri.Triangulation(loc[0], loc[1])
                mean = 0
                for edge in triag.edges:
                    mean += np.sqrt(np.sum((loc.transpose()[edge[0]] - loc.transpose()[edge[1]])**2))
                dist_between_peaks.append(mean / len(triag.edges))
        return values_of_peaks, dist_between_peaks, area_of_spots

    def get_error(self, species=0):

        axes = (1, 2)[:self.num_dims]

        error = self.stochastic[species] - self.analytic[species]

        if self.data_type == "molecules":
            error /= (self.total_num_molecules[species] * self.voxel_size)

        error = np.sqrt(self.voxel_size * np.sum(error ** 2, axes))

        return error

    def get_stationary_dist(self, take_sum=False):

        all_values = []
        all_counts = []

        for j in range(self.num_species):
            if take_sum:
                values, counts = np.unique(np.sum(self.__stochastic__[:, j], axis=(-1, -2)[:self.num_dims]),
                                           return_counts=True)
            else:
                values, counts = np.unique(self.__stochastic__[:, j], return_counts=True)
            counts = np.array(counts) / (self.num_time_steps * self.num_sims)
            all_values.append(values)
            all_counts.append(counts)

        return all_values, all_counts
