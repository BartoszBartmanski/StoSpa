//
// Created by bartosz on 31/12/18.
//

#ifndef STOSPA_SIMFUNCTIONS_HPP
#define STOSPA_SIMFUNCTIONS_HPP

#include "Simulation1d.hpp"
#include "Simulation2d.hpp"

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
        sim = make_unique<Simulation1d>(params);
    }
    else
    {
        sim = make_unique<Simulation2d>(params);
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

unique_ptr<JumpRate> get_jump_rates(vector<double> voxel_dims, const string &num_method="fdm", double alpha={}, vector<double> beta={})
{
    assert(!voxel_dims.empty());
    unique_ptr<JumpRate> jump_rate;

    if (voxel_dims.size() == 1)
    {
        jump_rate = make_unique<JumpRate1d>(voxel_dims);
    }
    else
    {
        if (num_method == "fdm")
        {
            jump_rate = make_unique<FDM>(voxel_dims, alpha);
        }
        else if (num_method == "fem")
        {
            jump_rate = make_unique<FEM>(voxel_dims);
        }
        else if (num_method == "fvm")
        {
            jump_rate = make_unique<FVM>(voxel_dims);
        }
        else if (num_method == "fet")
        {
            jump_rate = make_unique<FET>(voxel_dims, beta, 1000);
        }
        else if (num_method == "fetU")
        {
            jump_rate = make_unique<FETUniform>(voxel_dims, beta, 1000);
        }
        else
        {
            throw runtime_error("Unknown input for numerical method from which to derive the jump coefficients!");
        }
    }

    return move(jump_rate);
}
unique_ptr<JumpRate> get_jump_rates(Parameters &params)
{
    double domain_size = params.GetDomainBounds()[1] - params.GetDomainBounds()[0];
    vector<double> voxel_dims;
    double h_x = domain_size/params.GetNumVoxels();

    if (params.GetNumDims() == 1)
    {
        voxel_dims = {h_x};
    }
    else
    {
        voxel_dims = {h_x, h_x / params.GetKappa()};
    }

    return move(get_jump_rates(voxel_dims, params.GetNumMethod(), params.GetAlpha(), params.GetBeta()));
}

#endif //STOSPA_SIMFUNCTIONS_HPP
