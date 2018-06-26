#!/usr/bin/python3


"""Visualising the data from stochastic simulations.

Usage:
    parameters.py -h | --help
    parameters.py display <filename>...
    parameters.py <filename> change <key> <value>
    parameters.py <filename> add <key> <value>
    parameters.py <filename> delete <key>

Options:
    -h --help                       Show this screen
    change                          Changes the specified entry.
    add                             Adds an entry to the class stored in the file.
    delete                          Deletes an entry in the file.
"""

import pickle
from os import path
from docopt import docopt


class Parameters(object):
    """Class to store, change and get parameters in pickled dictionaries."""

    def __init__(self, file_name, save_dir="./", keys=None, params=None):
        assert isinstance(file_name, str)
        assert isinstance(save_dir, str)
        assert path.isdir(save_dir)
        if keys is None:
            keys = []
        if params is None:
            params = {}
        assert isinstance(params, dict)
        assert isinstance(keys, (list, tuple))

        self.file_name = file_name
        self.save_dir = save_dir
        self.keys = keys
        self.params = params

        try:
            f = open(self.save_dir + "/" + self.file_name, "rb")
        except IOError:
            # Ask for the parameters
            for key in keys:
                if key not in self.params:
                    val = input("{} = ".format(key))
                    self.params[key] = val

            # Save the parameters
            pickle.dump(self.params, open(self.save_dir + "/" + self.file_name, "wb"))
        else:
            # Unpickle the dictionary
            self.params = pickle.load(f)
            f.close()

            if len(self.keys) == 0:
                self.keys = self.params.keys()

            # Check that all the necessary parameters are present in this dictionary
            for key in self.keys:
                if key not in self.params.keys():
                    self.params[key] = input("Value for {} not found!\n{} = ".format(key, key))

            pickle.dump(self.params, open(self.save_dir + "/" + self.file_name, "wb"))

    def display(self):
        for parameter in self.keys:
            print("{} - {}".format(parameter, self.params[parameter]))

    def get_parameters(self):
        return self.params

    def change(self, name, val):
        assert name in self.params.keys(), "Parameter {} not found!".format(name)
        self.params[name] = str(val)
        # Save the parameters
        pickle.dump(self.params, open(self.save_dir + "/" + self.file_name, "wb"))

    def add(self, name, val):
        assert name not in self.params.keys(), "Parameter {} already in the dictionary!".format(name)
        self.params[name] = str(val)
        # Save the parameters
        pickle.dump(self.params, open(self.save_dir + "/" + self.file_name, "wb"))

    def delete(self, name):
        if name not in self.params.keys():
            print("No such entry!")
        else:
            del self.params[name]
        # Save the parameters
        pickle.dump(self.params, open(self.save_dir + "/" + self.file_name, "wb"))


if __name__ == '__main__':
    args = docopt(__doc__)

    parameters = Parameters(args["<filename>"][0], "./")
    if args["change"]:
        parameters.change(args["<key>"], args["<val>"])
    elif args["add"]:
        parameters.add(args["<key>"], args["<val>"])
        parameters = Parameters(args["<filename>"], "./")
    elif args["delete"]:
        parameters.delete(args["<key>"])
    elif args["display"]:
        green = "\033[1;32m"
        nc = "\033[0m"
        for arg in args["<filename>"]:
            parameters = Parameters(arg, "./")
            print(green + arg + nc)
            parameters.display()
