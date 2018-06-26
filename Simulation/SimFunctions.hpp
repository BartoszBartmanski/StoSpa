//
// Created by bartmanski on 27/03/18.
//

#ifndef STOSPA_SIMFUNCTIONS_HPP
#define STOSPA_SIMFUNCTIONS_HPP

#include "AbstractSimulation.hpp"

/**
 * Gets the simulation to the specified time point and calculates the error.
 * @param sim - reference to a simulation object
 * @param analytic - vector of analytic solution
 * @param end_time - end time of the the simulation
 * @return double
 */
double get_error(AbstractSimulation& sim, const vector<double>& analytic, double end_time);

/**
 * Calculates the size of the voxel
 * @param domain_bounds - bounds of the domain
 * @param num_voxels - number of voxels in the x-direction
 * @param kappa - voxel aspect ratio
 * @return double
 */
double get_voxel_size(const vector<double>& domain_bounds, const unsigned& num_voxels, const double& kappa);

/**
 * Calculates the midpoint of the middle voxel (problems arise when even number of voxels)
 * @param domain_bounds - bounds of the domain
 * @param num_voxels - number of voxels in the x-direction
 * @param kappa - voxel aspect ratio
 * @return vector<double>
 */
vector<double> get_midpoint(const vector<double>& domain_bounds, const unsigned& num_voxels, const double& kappa);

#endif //STOSPA_SIMFUNCTIONS_HPP
