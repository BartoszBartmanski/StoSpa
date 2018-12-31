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
 * Gets the simulation to the specified time point and calculates the error.
 * @param sim - reference to a simulation object
 * @param analytic - vector of analytic solution
 * @param end_time - end time of the the simulation
 * @return double
 */
double get_error(AbstractSimulation& sim, const vector<double>& analytic, double end_time)
{
    sim.Advance(end_time);
    return sim.GetError(analytic);
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

#endif //STOSPA_SIMFUNCTIONS_HPP
