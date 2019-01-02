//
// Created by bartosz on 31/12/18.
//

#ifndef STOSPA_SIMFUNCTIONS_HPP
#define STOSPA_SIMFUNCTIONS_HPP

#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"

/**
 * Brings together both constructors into a single line function.
 * (Templates are not used as they do not allow to specify dimension at runtime.)
 * @param params - parameters object that contains simulation parameters
 * @return - a pointer to the AbstractSimulation object
 */
unique_ptr<AbstractSimulation> simulator(Parameters params)
{
    unique_ptr<AbstractSimulation> sim;
    if (params.GetNumDims() == 1)
    {
        sim = make_unique<Simulation_1d>(params);
    }
    else
    {
        sim = make_unique<Simulation_2d>(params);
    }
    return move(sim);
}

/**
 * Calculates the midpoint of the middle voxel (problems arise when even number of voxels)
 * @param domain_bounds - bounds of the domain
 * @param num_voxels - number of voxels in the x-direction
 * @param kappa - voxel aspect ratio
 * @return vector<double>
 */
vector<double> get_midpoint(const vector<double>& domain_bounds, const unsigned& num_voxels, const double& kappa)
{
    vector<unsigned> voxels = {num_voxels, unsigned(kappa*num_voxels)};
    vector<double> midpoint(2);

    for (unsigned i=0; i < 2; i++)
    {
        double size = (domain_bounds[1] - domain_bounds[0]) / voxels[i];
        midpoint[i] = (voxels[i]/2 + 0.5 ) * size;
    }
    return midpoint;
}

/**
 * Returns pointer to JumpRate object that can give jump rates in specified directions.
 * @return unique_ptr<JumpRate>
 */

unique_ptr<JumpRate> get_jump_rate(unsigned dim, const string& num_method, double kappa, double h, double alpha, vector<double> beta)
{
    unique_ptr<JumpRate> jump_rate;

    if (dim == 1)
    {
        jump_rate = make_unique<JumpRate1d>(h);
    }
    else
    {
        if (num_method == "fdm")
        {
            jump_rate = make_unique<FDM>(kappa, h, alpha);
        }
        else if (num_method == "fem")
        {
            jump_rate = make_unique<FEM>(kappa, h);
        }
        else if (num_method == "fvm")
        {
            jump_rate = make_unique<FVM>(kappa, h);
        }
        else if (num_method == "fet")
        {
            jump_rate = make_unique<FET>(kappa, h, beta, 1000);
        }
        else if (num_method == "fetU")
        {
            jump_rate = make_unique<FETUniform>(kappa, h, beta, 1000);
        }
        else
        {
            throw runtime_error("Unknown input for numerical method from which to derive the jump coefficients!");
        }
    }

    return move(jump_rate);
}
unique_ptr<JumpRate> get_jump_rate(Parameters params)
{
    return move(get_jump_rate(params.GetNumDims(), params.GetNumMethod(), params.GetKappa(), params.GetH(), params.GetAlpha(), params.GetBeta()));
}

#endif //STOSPA_SIMFUNCTIONS_HPP
