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
    add                             Adds an entry to the class stored in the file.
    delete                          Deletes an entry in the file.
"""

import json
from os import path
from docopt import docopt


def _change_type(value, t):

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
        value[i] = _change_type(value[i], t)

    if len(value) == 1:
        value = value[0]

    return value


def get_dict_values(a_dict, keys):
    for key in keys:
        if key not in a_dict:
            val = input("{} [<input_type>] = ".format(key))
            if " " in val:
                val, t = val.split(" ")[:2]
            else:
                t = "str"
            a_dict[key] = change_type(val, t)
    return a_dict


class Parameters(object):
    """Class to store, change and get parameters in pickled dictionaries."""

    def __init__(self, path_to_file, keys=None, params=None):

        assert isinstance(path_to_file, str), type(path_to_file)
        assert path.exists(path_to_file), "{} file does not exist!".format(path_to_file)
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
            self.__params__ = get_dict_values(self.__params__, self.__keys__)

            # Save the parameters
            with open(self.__path_to_file__, "w") as outfile:
                json.dump(self.__params__, outfile, indent=4, sort_keys=True)
        else:
            # Unpickle the dictionary
            self.__params__ = json.load(f)
            f.close()

            if len(self.__keys__) == 0:
                self.__keys__ = list(self.__params__.keys())

            # Check that all the necessary parameters are present in this dictionary
            self.__params__ = get_dict_values(self.__params__, self.__keys__)

            # Save the parameters
            with open(self.__path_to_file__, "w") as outfile:
                json.dump(self.__params__, outfile, indent=4, sort_keys=True)

    def __getitem__(self, item):
        return self.__params__[item]

    def __setitem__(self, item, value):
        self.__params__[item] = value

    def get_dict(self):
        return self.__params__

    def display(self):
        self.__keys__.sort()
        for parameter in self.__keys__:
            print("  {} - {}".format(parameter, self.__params__[parameter]))

    def change(self, name, val):
        assert name in self.__params__.keys(), "Parameter {} not found!".format(name)
        self.__params__[name] = val

        # Save the parameters
        with open(self.__path_to_file__, "w") as outfile:
            json.dump(self.__params__, outfile, indent=4, sort_keys=True)

    def add(self, name, val):
        assert name not in self.__params__.keys(), "Parameter {} already in the dictionary!".format(name)
        self.__params__[name] = val

        # Save the parameters
        with open(self.__path_to_file__, "w") as outfile:
            json.dump(self.__params__, outfile, indent=4, sort_keys=True)

    def delete(self, name):
        if name not in self.__params__.keys():
            print("No such entry!")
        else:
            del self.__params__[name]

        # Save the parameters
        with open(self.__path_to_file__, "w") as outfile:
            json.dump(self.__params__, outfile, indent=4, sort_keys=True)


if __name__ == '__main__':
    args = docopt(__doc__)

    parameters = Parameters(args["<filename>"])
    if args["change"]:
        parameters.change(args["<key>"], change_type(args["<value>"], args["<type>"]))
    elif args["add"]:
        parameters.add(args["<key>"], change_type(args["<value>"], args["<type>"]))
    elif args["delete"]:
        parameters.delete(args["<key>"])
    else:
        parameters.display()
