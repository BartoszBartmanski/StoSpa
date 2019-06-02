#!/usr/bin/env python3

import numpy as np


class Positions(object):
    def __init__(self, num_axes, fig_params, from_top=False):

        # Check that fig_params contains all the necessary info
        needed = ["fig_width", "fig_height", "free_width", "free_height",
                  "left", "bottom", "width", "height"]
        assert isinstance(fig_params, dict)
        for x in needed:
            assert x in fig_params, "{} is missing!".format(x)

        # Make sure that num_axes is a tuple
        # self.__num_axes__ = [num_rows, num_cols]
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
        self.__fig_w__ = fig_params["fig_width"]  # in inches
        self.__fig_h__ = fig_params["fig_height"]  # in inches

        self.__positions__ = []
        self.__rows__ = []
        self.__cols__ = []
        pos0 = np.array([fig_params["left"]/fig_params["fig_width"],
                         fig_params["bottom"]/fig_params["fig_height"],
                         fig_params["width"]/fig_params["fig_width"],
                         fig_params["height"]/fig_params["fig_height"]])

        w_space = fig_params["width"] + fig_params["free_width"]
        h_space = fig_params["height"] + fig_params["free_height"]
        add_x = np.array([w_space / fig_params["fig_width"], 0, 0, 0])
        add_y = np.array([0, h_space / fig_params["fig_height"], 0, 0])

        for i in range(self.__num_axes__[0]):
            for j in range(self.__num_axes__[1]):
                self.__positions__.append(pos0 + j*add_x + i*add_y)
                self.__rows__.append(i)
                self.__cols__.append(j)

        if from_top:
            self.__positions__ = self.flip_array(self.__positions__,
                                                 self.__num_axes__)
        else:
            self.__rows__ = self.flip_array(self.__rows__, self.__num_axes__)
            self.__cols__ = self.flip_array(self.__cols__, self.__num_axes__)

    @staticmethod
    def flip_array(vec, new_shape):
        if isinstance(new_shape, int):
            new_shape = [1, new_shape]
        else:
            new_shape = list(new_shape)

        # Turn vec into numpy array and get its shape
        vec = np.array(vec)
        shape = vec.shape

        # Get a temporary shape to flip the array
        m = "Can't reorder the array!"
        assert np.prod(shape) % np.prod(new_shape) == 0, m
        temp_shape = new_shape + [int(np.prod(shape)/np.prod(new_shape))]

        # Cut the list
        vec = np.array(vec[:np.prod(new_shape)])
        # Flip the order
        vec = np.flip(vec.reshape(temp_shape), axis=0).reshape(shape)

        return vec

    def get_position(self, idx=0):
        return self.__positions__[idx]

    def get_width(self):
        return np.array([self.__w__, 0, 0, 0])

    def get_height(self):
        return np.array([0, self.__h__, 0, 0])

    def get_free_width(self):
        return np.array([self.__free_w__, 0, 0, 0])

    def get_free_height(self):
        return np.array([0, self.__free_h__, 0, 0])

    def get_fig_size(self):
        return [self.__fig_w__, self.__fig_h__]

    def get_row(self, idx):
        return self.__rows__[idx]

    def get_col(self, idx):
        return self.__cols__[idx]

    def get_num_axes(self):
        return np.prod(self.__num_axes__)

    def shift_pos_x(self, value, cols):
        """Value in inches"""
        if isinstance(cols, int):
            cols = [cols]
        assert isinstance(cols, (tuple, list))
        shift = np.array([value/self.__fig_w__, 0, 0, 0])

        for i in range(self.get_num_axes()):
            if self.__cols__[i] in cols:
                self.__positions__[i] += shift

    def shift_pos_y(self, value, rows):
        """Value in inches"""
        if isinstance(rows, int):
            rows = [rows]
        assert isinstance(rows, (tuple, list))
        shift = np.array([0, value/self.__fig_h__, 0, 0])

        for i in range(self.get_num_axes()):
            if self.__rows__[i] in rows:
                self.__positions__[i] += shift

    def get_adj_position(self, idx, shift_x=0.0, shift_y=0.0,
                         width=None, height=None):
        new_position = self.__positions__[idx].copy()
        if width is None:
            width = self.__w__ * self.__fig_w__
        if height is None:
            height = self.__h__ * self.__fig_h__

        width /= self.__fig_w__
        height /= self.__fig_h__
        shift_x /= self.__fig_w__
        shift_y /= self.__fig_h__

        new_position[2] = width
        new_position[3] = height

        if shift_x > 0.0:
            new_position[0] += (self.__w__ + abs(shift_x))
        elif shift_x < 0.0:
            new_position[0] -= (abs(shift_x) + width)

        if shift_y > 0.0:
            new_position[1] += (self.__h__ + abs(shift_y))
        elif shift_y < 0.0:
            new_position[1] -= (abs(shift_y) + height)

        return new_position
