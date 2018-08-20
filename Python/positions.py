#!/usr/bin/python3

import numpy as np


def flip_array(vec, new_shape):
    if isinstance(new_shape, int):
        new_shape = [1, new_shape]

    # Turn vec into numpy array and get its shape
    vec = np.array(vec)
    shape = vec.shape

    # Get a temporary shape to flip the array
    assert np.prod(shape) % np.prod(new_shape) == 0, "Can't reorder the array!"
    temp_shape = new_shape + [int(np.prod(shape)/np.prod(new_shape))]

    # Cut the list
    vec = np.array(vec[:np.prod(new_shape)])
    # Flip the order
    vec = np.flip(vec.reshape(temp_shape), axis=0).reshape(shape)

    return vec


class Positions(object):
    def __init__(self, num_axes, fig_params):

        # Check that fig_params contains all the necessary info
        needed = ["fig_width", "fig_height", "free_width", "free_height",
                  "left", "bottom", "width", "height"]
        assert isinstance(fig_params, dict)
        for x in needed:
            assert x in fig_params, "{} is missing!".format(x)

        # Make sure that num_axes is a tuple
        if isinstance(num_axes, int):
            self.__num_axes__ = [1, num_axes]
        else:
            self.__num_axes__ = num_axes
        assert isinstance(self.__num_axes__, (tuple, list))
        assert len(self.__num_axes__) == 2

        self.__h__ = fig_params["height"] / fig_params["fig_height"]
        self.__w__ = fig_params["width"] / fig_params["fig_width"]
        self.__free_w__ = fig_params["free_width"] / fig_params["fig_width"]
        self.__free_h__ = fig_params["free_height"] / fig_params["fig_height"]

        self.__positions__ = []
        pos0 = np.array([fig_params["left"]/fig_params["fig_width"],
                         fig_params["bottom"]/fig_params["fig_height"],
                         fig_params["width"]/fig_params["fig_width"],
                         fig_params["height"]/fig_params["fig_height"]])

        w_space = fig_params["width"] + fig_params["free_width"]
        h_space = fig_params["height"] + fig_params["free_height"]
        add_x = np.array([w_space / fig_params["fig_width"], 0, 0, 0])
        add_y = np.array([0, h_space / fig_params["fig_height"], 0, 0])
        for i in range(num_axes[0]):
            for j in range(num_axes[1]):
                self.__positions__.append(pos0 + j*add_x + i*add_y)

    def get_positions(self, from_top=False):
        if from_top:
            return flip_array(self.__positions__, self.__num_axes__)
        else:
            return self.__positions__

    def get_width(self):
        return np.array([self.__w__, 0, 0, 0])

    def get_height(self):
        return np.array([0, self.__h__, 0, 0])

    def get_free_width(self):
        return np.array([self.__free_w__, 0, 0, 0])

    def get_free_height(self):
        return np.array([0, self.__free_h__, 0, 0])
