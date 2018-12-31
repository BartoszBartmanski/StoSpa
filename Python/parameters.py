#!/usr/bin/python3

"""Visualising the data from stochastic simulations.

Usage:
    parameters.py -h | --help
    parameters.py <filename>
    parameters.py <filename> change [--] <key> <value> [<type>]
    parameters.py <filename> add [--] <key> <value> [<type>]
    parameters.py <filename> delete [--] <key>

Options:
    -h --help                       Show this screen
    change                          Changes the specified entry.
    add                             Adds an entry to the file.
    delete                          Deletes an entry in the file.
"""

import json
from docopt import docopt


def change_item(value, t):

    if t == "int" or t == "i":
        value = int(value)
    elif t == "float" or t == "f":
        value = float(value)
    else:
        value = str(value)

    return value


def change_type(value, t):
    assert isinstance(value, str), type(value)
    value = value.split(",")

    for i in range(len(value)):
        value[i] = change_item(value[i], t)

    if len(value) == 1:
        value = value[0]

    return value


class Parameters(object):
    """Class to store, change and get parameters in json'd dictionaries."""

    def __init__(self, path_to_file, keys=None, params=None):

        assert isinstance(path_to_file, str), type(path_to_file)
        if keys is None:
            keys = []
        if params is None:
            params = {}
        assert isinstance(params, dict)
        assert isinstance(keys, (list, tuple))

        self.__path_to_file__ = path_to_file
        self.__keys__ = keys
        self.__params__ = params

        try:
            f = open(self.__path_to_file__, "r")
        except IOError:
            # Ask for the parameters
            self.__check_keys__()

            # Save the parameters
            self.__save__()
        else:
            # Open the dictionary
            self.__params__ = json.load(f)
            f.close()

            # Check that all the necessary keys are present
            if len(self.__keys__) == 0:
                self.__keys__ = list(self.__params__.keys())
            else:
                self.__check_keys__()

            # Save the parameters
            self.__save__()

    def __getitem__(self, item):
        if item not in self.__params__.keys():
            self.__get_value__(item)
            self.__save__()
        return self.__params__[item]

    def __setitem__(self, item, value):
        self.__params__[item] = value

    def __get_value__(self, key):
        val = input("{} [<input_type>] = ".format(key))
        if " " in val:
            val, t = val.split(" ")[:2]
        else:
            t = "str"
        self.__params__[key] = change_type(val, t)

    def __check_keys__(self):
        for key in self.__keys__:
            if key not in self.__params__.keys():
                self.__get_value__(key)

    def __save__(self):
        with open(self.__path_to_file__, "w") as outfile:
            json.dump(self.__params__, outfile, indent=4, sort_keys=True)

    def get_dict(self):
        return self.__params__

    def update(self, another_obj):
        assert isinstance(another_obj, Parameters)
        self.__params__.update(another_obj.get_dict())

    def display(self):
        self.__keys__.sort()
        for parameter in self.__keys__:
            print("  {} - {}".format(parameter, self.__params__[parameter]))

    def change(self, name, val):
        m = "Parameter {} not found!".format(name)
        assert name in self.__params__.keys(), m
        self.__params__[name] = val

        # Save the parameters
        self.__save__()

    def add(self, name, val):
        m = "Parameter {} already in the dictionary!".format(name)
        assert name not in self.__params__.keys(), m
        self.__params__[name] = val

        # Save the parameters
        self.__save__()

    def delete(self, name):
        if name not in self.__params__.keys():
            print("No such entry!")
        else:
            del self.__params__[name]

        # Save the parameters
        self.__save__()


if __name__ == '__main__':
    args = docopt(__doc__)

    parameters = Parameters(args["<filename>"])
    if args["change"]:
        args["<value>"] = change_type(args["<value>"], args["<type>"])
        parameters.change(args["<key>"], args["<value>"])
    elif args["add"]:
        args["<value>"] = change_type(args["<value>"], args["<type>"])
        parameters.add(args["<key>"], args["<value>"])
    elif args["delete"]:
        parameters.delete(args["<key>"])
    else:
        parameters.display()
