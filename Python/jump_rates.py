#!/usr/bin/env python3

"""
Plots jump rates against parameters.

Usage:
    jump_rates.py -h | --help
    jump_rates.py lambda [<index>] [options]
    jump_rates.py theta [<index>] [options]
    jump_rates.py check [options]

Options:
    -h --help               Show this screen.
    --length=<val>          Vertical voxel length. [default: 0.025]
    --kappa=<val>           Aspect ratio of the voxels. [default: 1.0]
    --diff=<val>            Macro-scale diffusion constant. [default: 1.0]
    --save=<path>           Save the figure in <path>.

"""

import numpy as np
from matplotlib import pyplot as plt

from docopt import docopt


class JumpRate(object):
    def __init__(self, diff, h):
        """
        diff - diffusion coefficient
        h - compartment size (in the y-direction in 2d)
        """
        assert isinstance(diff, (float, int))
        assert isinstance(h, (float, int))

        self.__diff__ = float(diff)
        self.__h__ = float(h)
        self.__num_dims__ = 1

    @staticmethod
    def __get_parameter__(param, const=False):
        if param is None:
            param = 1
        if isinstance(param, (int, float)):
            param = float(param)
            if const:
                param = 1
        elif isinstance(param, (tuple, list, np.ndarray)):
            param = np.array(param)
            if const:
                param = np.ones(len(param))
        return param

    def get_diff(self):
        return self.__diff__

    def get_h(self):
        return self.__h__

    def get_lambda(self, index, parameter=None):
        assert isinstance(index, int)
        if index == 1:
            return self.__lambda_1__(parameter)
        elif index == 2:
            return self.__lambda_2__(parameter)
        else:
            return self.__lambda_0__(parameter)

    def get_theta(self, index, parameter=None):
        assert isinstance(index, int)
        l_1 = self.__lambda_1__(parameter)
        l_2 = self.__lambda_2__(parameter)
        l_0 = self.__lambda_0__(parameter)

        if index == 1:
            return l_1 / l_0
        elif index == 2:
            return l_2 / l_0
        else:
            return (l_1 + l_2) / l_0

    def __lambda_1__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * parameter

    def __lambda_2__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * parameter

    def __lambda_0__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = 2.0 * self.__diff__ / (self.__h__ ** 2)
        return val * parameter


class JumpRate2D(JumpRate):
    def __init__(self, diff, h, kappa):
        """
        diff - diffusion coefficient
        h - compartment size (in the y-direction in 2d)
        kappa - compartment aspect ratio (horizontal/vertical)
        """
        super().__init__(diff, h)
        self.__kappa__ = self.__get_parameter__(kappa)
        self.__num_dims__ = 2

    def get_kappa(self):
        return self.__kappa__

    def get_lambda(self, index, parameter=None):
        assert isinstance(index, int)
        if index == 1:
            return self.__lambda_1__(parameter)
        elif index == 2:
            return self.__lambda_2__(parameter)
        elif index == 3:
            return self.__lambda_3__(parameter)
        else:
            return self.__lambda_0__(parameter)

    def get_theta(self, index, parameter=None):
        assert isinstance(index, int)
        l_1 = self.__lambda_1__(parameter)
        l_2 = self.__lambda_2__(parameter)
        l_3 = self.__lambda_3__(parameter)
        l_0 = self.__lambda_0__(parameter)

        if index == 1:
            return l_1 / l_0
        elif index == 2:
            return l_2 / l_0
        elif index == 3:
            return l_3 / l_0
        else:
            return (2*l_1 + 4*l_2 + 2*l_3) / l_0

    def __lambda_1__(self, parameter):
        raise NotImplementedError

    def __lambda_2__(self, parameter):
        raise NotImplementedError

    def __lambda_3__(self, parameter):
        raise NotImplementedError

    def __lambda_0__(self, parameter):
        l_1 = self.__lambda_1__(parameter)
        l_2 = self.__lambda_2__(parameter)
        l_3 = self.__lambda_3__(parameter)
        return 2 * l_1 + 4 * l_2 + 2 * l_3

    def check(self, parameter):
        l1 = self.__lambda_1__(parameter)
        l2 = self.__lambda_2__(parameter)
        l3 = self.__lambda_3__(parameter)
        l0 = self.__lambda_0__(parameter)
        return np.absolute(((2 * l1 + 4*l2 + 2*l3) - l0)/l0)


class FDM(JumpRate2D):

    def __lambda_1__(self, alpha):
        alpha = self.__get_parameter__(alpha)
        val = self.__diff__ * (1.0 - self.__kappa__ * alpha)
        val /= (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val

    def __lambda_2__(self, alpha):
        alpha = self.__get_parameter__(alpha)
        val = self.__diff__ * alpha
        val /= (2.0 * self.__kappa__ * self.__h__ ** 2)
        return val

    def __lambda_3__(self, alpha):
        alpha = self.__get_parameter__(alpha)
        val = self.__diff__ * (self.__kappa__**2 - self.__kappa__ * alpha)
        val /= (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val


class FEM(JumpRate2D):
    def __init__(self, diff, h, kappa, alt_lambda_0=False):
        super().__init__(diff, h, kappa)
        self.__a__ = 2 * 2.0
        self.__b__ = 2 * 1.0
        self.__c__ = 1.0
        self.__d__ = 2 * 4.0
        self.__denominator__ = 6.0 * self.__kappa__ ** 2 * self.__h__ ** 2
        self.__alt_lambda_0__ = alt_lambda_0  # Use sum of lambdas for lambda_0

    def __lambda_1__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ * (self.__a__ - self.__b__ * self.__kappa__ ** 2)
        val /= self.__denominator__
        return val * parameter

    def __lambda_2__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ * self.__c__ * (self.__kappa__ ** 2 + 1.0)
        val /= self.__denominator__
        return val * parameter

    def __lambda_3__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ * (self.__a__ * self.__kappa__ ** 2 - self.__b__)
        val /= self.__denominator__
        return val * parameter

    def __lambda_0__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        if self.__alt_lambda_0__:
            val = 2 * self.__lambda_1__(parameter)
            val += 4 * self.__lambda_2__(parameter)
            val += 2 * self.__lambda_3__(parameter)
        else:
            val = self.__diff__ * self.__d__ * (self.__kappa__ ** 2 + 1)
            val /= self.__denominator__
        return val * parameter


class FVM(JumpRate2D):

    def __lambda_1__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ / (self.__kappa__ ** 2 * self.__h__ ** 2)
        return val * parameter

    def __lambda_2__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        return 0 * parameter

    def __lambda_3__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        val = self.__diff__ / (self.__h__ ** 2)
        return val * parameter


class FET(JumpRate2D):

    def __init__(self, diff, h, kappa, trunc_order=100):
        super().__init__(diff, h, kappa)
        assert isinstance(trunc_order, (float, int))
        self.__trunc_order__ = int(trunc_order)

    def __lambda_1__(self, beta_y):
        beta_y = self.__get_parameter__(beta_y)
        theta_1 = 0

        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                val = 8.0 * (-1)**(k + 1) * (2.0 * k - 1.0)
                val /= (np.pi**2 * temp * (2.0 * j - 1.0))
                theta_1 += np.sin((j - 0.5) * np.pi * beta_y) * val

        return theta_1 * self.__lambda_0__(beta_y)

    def __lambda_2__(self, beta_x, beta_y=None):
        if beta_y is None:
            beta_y = beta_x
        beta_x = self.__get_parameter__(beta_x)
        beta_y = self.__get_parameter__(beta_y)
        assert len(beta_x) == len(beta_y)

        theta_2 = 0.0
        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                val1 = (2*j-1)*(np.sin(np.pi*(beta_x*k-0.5*beta_x+k))+1.0)
                val1 *= self.__kappa__**2
                val1 /= (2 * k - 1)
                val2 = (2*k-1)*(np.sin(np.pi*(beta_y*j-0.5*beta_y+j))+1.0)
                val2 /= (2 * j - 1)
                theta_2 += 4.0*(-1)**(j+k) * (val1 + val2) / (np.pi**2 * temp)

        return theta_2 * self.__lambda_0__(beta_x)

    def __lambda_3__(self, beta_x):
        beta_x = self.__get_parameter__(beta_x)

        theta_3 = 0
        for k in range(1, self.__trunc_order__+1):
            for j in range(1, self.__trunc_order__+1):
                temp = (self.__kappa__**2 * (2*j-1)**2 + (2*k-1)**2)

                val = 8.0 * (-1)**(j+1) * (2 * j - 1) * self.__kappa__**2
                val /= (np.pi**2 * temp * (2 * k - 1))
                theta_3 += np.sin((k - 0.5) * np.pi * beta_x) * val

        return theta_3 * self.__lambda_0__(beta_x)

    def __lambda_0__(self, parameter):
        parameter = self.__get_parameter__(parameter, const=True)
        k = np.arange(1, self.__trunc_order__ + 1)
        kk, jj = np.meshgrid(k, k)

        top = 64 * (-1)**(jj + kk) * self.__kappa__**2 * self.__h__**2
        bottom_1 = self.__diff__ * np.pi**4 * (2*kk-1)*(2*jj-1)
        bottom_2 = (self.__kappa__**2 * (2*jj-1)**2 + (2*kk-1)**2)
        return 1.0/np.sum(top/(bottom_1*bottom_2)) * parameter


if __name__ == '__main__':
    args = docopt(__doc__)

    args["--length"] = float(args["--length"])
    args["--kappa"] = float(args["--kappa"])
    args["--diff"] = float(args["--diff"])
    if args["<index>"] in ["1", "2", "3"]:
        args["<index>"] = int(args["<index>"])
    else:
        args["<index>"] = 0

    fdm = FDM(args["--diff"], args["--length"], args["--kappa"])
    fem = FEM(args["--diff"], args["--length"], args["--kappa"])
    fvm = FVM(args["--diff"], args["--length"], args["--kappa"])
    fet = FET(args["--diff"], args["--length"], args["--kappa"])

    x = np.linspace(0, 1, 50)
    fig = plt.figure()

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

    elif args["check"]:
        plt.plot(x, fdm.check(x), label="FDM")
        plt.plot(x, fem.check(x), label="FEM")
        plt.plot(x, fvm.check(x), label="FVM")
        plt.plot(x, fet.check(x), label="FET")
        plt.ylabel(r"$\left|\sum_i \lambda_i - \lambda_0\right|/\lambda_0$")

    plt.xlabel(r"$\alpha$, $\beta$")
    plt.legend()
    plt.tight_layout()

    if args["--save"]:
        plt.savefig(args["--save"])
    else:
        plt.show()
