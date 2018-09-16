#!/usr/bin/python3

"""
Script to plot any miscellaneous data.

Usage:
    jump_rates.py -h | --help
    jump_rates.py lambda [<index>] [options]
    jump_rates.py theta [<index>] [options]

Options:
    -h --help               Show this screen.
    --length=<val>          Vertical voxel length. [default: 0.025]
    --kappa=<val>           Aspect ratio of the voxels. [default: 1.0]
    --diff=<val>            Macro-scale diffusion constant. [default: 1.0]

"""

from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt


def get_parameter(param):
    if isinstance(param, (int, float)):
        param = float(param)
    elif isinstance(param, (tuple, list, np.ndarray)):
        param = np.array(param)
    return param


class JumpRate(object):
    def __init__(self, diff, h, kappa):
        """
        diff - diffusion coefficient
        h - compartment size (in the y-direction in 2d)
        kappa - compartment aspect ratio (horizontal/vertical)
        """
        assert isinstance(diff, (float, int))
        assert isinstance(h, (float, int))
        assert isinstance(kappa, (float, int))

        self.__diff__ = float(diff)
        self.__h__ = float(h)
        self.__kappa__ = float(kappa)

    def get_diff(self):
        return self.__diff__

    def get_h(self):
        return self.__h__

    def get_kappa(self):
        return self.__kappa__

    def get_lambda(self, index, parameter, num_dims=2):
        assert isinstance(index, int)
        if index == 1:
            return self.__lambda_1__(parameter)
        elif index == 2:
            return self.__lambda_2__(parameter)
        elif index == 3:
            return self.__lambda_3__(parameter)
        else:
            l_1 = self.__lambda_1__(parameter)
            l_2 = self.__lambda_2__(parameter)
            l_3 = self.__lambda_3__(parameter)
            if num_dims == 2:
                return (2 * l_1 + 4 * l_2 + 2 * l_3)
            else:
                return l_1 + l_2

    def get_theta(self, index, parameter, num_dims=2):
        assert isinstance(index, int)
        l_1 = self.__lambda_1__(parameter)
        l_2 = self.__lambda_2__(parameter)
        l_3 = self.__lambda_3__(parameter)
        if num_dims == 2:
            l_0 = 2 * l_1 + 4 * l_2 + 2 * l_3
        else:
            l_0 = l_1 + l_2
        if index == 1:
            return l_1 / l_0
        elif index == 2:
            return l_2 / l_0
        elif index == 3:
            return l_3 / l_0
        else:
            return l_0 / l_0

    def __lambda_1__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * np.ones(len(parameter))

    def __lambda_2__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * np.ones(len(parameter))

    def __lambda_3__(self, parameter):
        parameter = get_parameter(parameter)
        return 0 * np.ones(len(parameter))


class FDM(JumpRate):

    def __lambda_1__(self, alpha):
        alpha = get_parameter(alpha)
        val = self.__diff__ * (1.0 - self.__kappa__ * alpha)
        val /= (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val

    def __lambda_2__(self, alpha):
        alpha = get_parameter(alpha)
        val = self.__diff__ * alpha
        val /= (2.0 * self.__kappa__ * self.__h__ ** 2)
        return val

    def __lambda_3__(self, alpha):
        alpha = get_parameter(alpha)
        val = self.__diff__ * (self.__kappa__**2 - self.__kappa__ * alpha)
        val /= (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val


class FEM(JumpRate):

    def __lambda_1__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ * (2.0 - self.__kappa__ ** 2)
        val /= (3.0 * self.__kappa__ ** 2 * self.__h__ ** 2)
        return val * np.ones(len(parameter))

    def __lambda_2__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ * (self.__kappa__ ** 2 + 1.0)
        val /= (6.0 * self.__kappa__ ** 2 * self.__h__ ** 2)
        return val * np.ones(len(parameter))

    def __lambda_3__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ * (2.0 * self.__kappa__ ** 2 - 1.0)
        val /= (3.0 * self.__kappa__ ** 2 * self.__h__ ** 2)
        return val * np.ones(len(parameter))


class FVM(JumpRate):

    def __lambda_1__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ / (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val * np.ones(len(parameter))

    def __lambda_2__(self, parameter):
        parameter = get_parameter(parameter)
        return 0 * np.ones(len(parameter))

    def __lambda_3__(self, parameter):
        parameter = get_parameter(parameter)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * np.ones(len(parameter))


class FET(JumpRate):

    def __init__(self, diff, h, kappa, trunc_order=100):
        super().__init__(diff, h, kappa)
        assert isinstance(trunc_order, (float, int))
        self.__trunc_order__ = int(trunc_order)

    def __lambda_1__(self, beta_y):
        beta_y = get_parameter(beta_y)
        theta_1 = 0
        lambda_0 = 0
        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                lambda_0 += (-1)**(j + k) / ((2 * k - 1) * (2 * j - 1) * temp)

                val = 8.0 * (-1)**(k + 1) * (2.0 * k - 1.0)
                val /= (np.pi**2 * temp * (2.0 * j - 1.0))
                theta_1 += np.sin(0.5 * (2 * j - 1) * np.pi * beta_y) * val

        kappa2h2 = self.__kappa__ ** 2 * self.__h__ ** 2
        lambda_0 = np.pi**4 * self.__diff__ / (64 * kappa2h2 * lambda_0)

        return theta_1 * lambda_0

    def __lambda_2__(self, beta_x, beta_y=None):
        if beta_y is None:
            beta_y = beta_x
        beta_x = get_parameter(beta_x)
        beta_y = get_parameter(beta_y)
        theta_2 = 0.0
        lambda_0 = 0.0
        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                lambda_0 += (-1)**(j+k) / ((2 * k - 1) * (2 * j - 1) * temp)

                val1 = (2*k-1)*(np.sin(np.pi*(beta_y*j-0.5*beta_y+j))+1.0)
                val1 /= (2 * j - 1)
                val2 = (2*j-1)*(np.sin(np.pi*(beta_x*k-0.5*beta_x+k))+1.0)
                val2 /= (self.__kappa__**2 * (2 * k - 1))
                theta_2 += 4.0*(-1)**(j+k) * (val1 + val2) / (np.pi**2 * temp)

        kappa2h2 = self.__kappa__ ** 2 * self.__h__ ** 2
        lambda_0 = np.pi**4 * self.__diff__ / (64 * kappa2h2 * lambda_0)

        return theta_2 * lambda_0

    def __lambda_3__(self, beta_x):
        beta_x = get_parameter(beta_x)
        theta_3 = 0
        lambda_0 = 0
        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                lambda_0 += (-1)**(j+k) / ((2 * k - 1) * (2 * j - 1) * temp)

                val = 8.0 * (-1) ** (j+1) * (2 * j - 1) * self.__kappa__**2
                val /= (np.pi**2 * temp * (2 * k - 1))
                theta_3 += np.sin(0.5 * (2 * k - 1) * np.pi * beta_x) * val

        kappa2h2 = self.__kappa__ ** 2 * self.__h__ ** 2
        lambda_0 = np.pi**4 * self.__diff__ / (64 * kappa2h2 * lambda_0)

        return theta_3 * lambda_0


if __name__ == '__main__':
    args = docopt(__doc__)

    args["--length"] = float(args["--length"])
    args["--kappa"] = float(args["--kappa"])
    args["--diff"] = float(args["--diff"])
    if args["<index>"] is None:
        args["<index>"] = 0
    else:
        args["<index>"] = int(args["<index>"])

    fdm = FDM(args["--diff"], args["--length"], args["--kappa"])
    fem = FEM(args["--diff"], args["--length"], args["--kappa"])
    fvm = FVM(args["--diff"], args["--length"], args["--kappa"])
    fet = FET(args["--diff"], args["--length"], args["--kappa"])

    x = np.linspace(0, 1, 50)

    if args["lambda"]:

        plt.plot(x, fdm.get_lambda(args["<index>"], x), label="FDM")
        plt.plot(x, fem.get_lambda(args["<index>"], x), label="FEM")
        plt.plot(x, fvm.get_lambda(args["<index>"], x), label="FVM")
        plt.plot(x, fet.get_lambda(args["<index>"], x), label="FET")
        plt.ylabel(r"$\lambda_{}$".format(args["<index>"]))

    elif args["theta"]:
        plt.plot(x, fdm.get_theta(args["<index>"], x), label="FDM")
        plt.plot(x, fem.get_theta(args["<index>"], x), label="FEM")
        plt.plot(x, fvm.get_theta(args["<index>"], x), label="FVM")
        plt.plot(x, fet.get_theta(args["<index>"], x), label="FET")
        plt.ylabel(r"$\theta_{}$".format(args["<index>"]))

    plt.xlabel(r"$\alpha$, $\beta$")
    plt.legend()
    plt.show()
