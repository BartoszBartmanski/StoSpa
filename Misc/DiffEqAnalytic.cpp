//
// Created by bartosz on 24/02/18.
//

#include "DiffEqAnalytic.hpp"

DiffEqAnalytic::DiffEqAnalytic(unsigned num_dims,
                               double time,
                               vector<double> x_0,
                               vector<double> domain,
                               vector<unsigned> num_points,
                               unsigned total_molecules,
                               double diff,
                               double decay,
                               double prod,
                               unsigned trunc_order,
                               bool voxels)
{
    // Check the input is sensible
    assert(num_dims == 1 or num_dims == 2);
    assert(time > 0);
    assert(total_molecules > 0);
    assert(diff > 0);
    assert(decay >= 0);
    assert(prod >= 0);
    assert(trunc_order > 0);
    for (unsigned i=0; i<num_dims; i++)
    {
        assert(domain[0] < x_0[i] and x_0[i] < domain[1]);
    }

    // This function returns the analytic solution to the heat equation with decay i.e. du/dt = nabla^2 u + k_prod - k_decay u
    // k is the decay rate

    prod /= total_molecules;
    vector<double> u;
    double l = domain[1] - domain[0];
    vector<double> spacing = {0, 0};
    if (voxels)
    {
        spacing = {l / num_points[0], l / num_points[1]};
    }

    // To save computation time - probably
    double pi_l = M_PI / l;

    if (num_dims == 1)
    {
        double factor = 2.0 / l;
        vector<double> x = linspace(domain[0]+0.5*spacing[0], domain[1]-0.5*spacing[0], num_points[0]);
        u = vector<double>(num_points[0], 1.0 / l);

        for (unsigned i = 0; i < num_points[0]; i++)
        {
            for (unsigned n = 1; n < trunc_order + 1; n++)
            {
                double val = factor * cos(n * pi_l * x_0[0]) * cos(n * pi_l * x[i]) * exp(-diff * pow(n * pi_l, 2) * time);
                u[i] += val;
            }
        }
    }
    else if (num_dims == 2)
    {
        double factor_0 = 1.0 / pow(l, 2);
        double factor_1 = 2.0 * factor_0;
        double factor_2 = 2.0 * factor_1;
        vector<double> x = linspace(domain[0]+0.5*spacing[0], domain[1]-0.5*spacing[0], num_points[0]);
        vector<double> y = linspace(domain[0]+0.5*spacing[1], domain[1]-0.5*spacing[1], num_points[1]);
        u = vector<double>(num_points[0]*num_points[1], factor_0);

        for (unsigned j=0; j < num_points[1]; j++)
        {
            for (unsigned i = 0; i < num_points[0]; i++)
            {
                for (unsigned m = 1; m < trunc_order + 1; m++)
                {
                    u[j*num_points[0] + i] += (factor_1 * cos(m * pi_l * x_0[1]) * cos(m * pi_l * y[j]) * exp(-diff * time * pow(m * pi_l, 2)));
                }

                for (unsigned n = 1; n < trunc_order + 1; n++)
                {
                    u[j*num_points[0] + i] += (factor_1 * cos(n * pi_l * x_0[0]) * cos(n * pi_l * x[i]) * exp(-diff * time * pow(n * pi_l, 2)));
                }

                for (unsigned m = 1; m < trunc_order + 1; m++)
                {
                    for (unsigned n = 1; n < trunc_order + 1; n++)
                    {
                        u[j * num_points[0] + i] += (factor_2 * cos(n * pi_l * x_0[0]) * cos(n * pi_l * x[i]) * cos(m * pi_l * x_0[1]) * cos(m * pi_l * y[j]) * exp(-diff * time * (pow(m * pi_l, 2) + pow(n * pi_l, 2))));
                    }
                }
            }
        }

    }

    // Additional terms and factors for when there is decay or increase in the number of molecules
    if (decay > 0.0 and prod == 0.0)
    {
        u = exp(-decay * time) * u;
    }
    else if (decay == 0.0 and prod > 0.0)
    {
        u = u + prod * time;
    }
    else if (decay > 0.0 and prod > 0.0)
    {
        u = exp(-decay * time) * u;
        u = u + prod * (1.0 - exp(-decay * time)) / decay;
    }

    mSolution = u;

}

vector<double> DiffEqAnalytic::GetAnalytic()
{
    return mSolution;
}
