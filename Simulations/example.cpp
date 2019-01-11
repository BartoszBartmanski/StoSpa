//
// Created by bartosz on 02/01/19.
//

#include "Simulator.hpp"

int main(int argc, const char** argv)
{
    // Declare the simulation name
    string sim_name = "example.dat";

    // Define all the variables
    unsigned num_runs = 1;
    unsigned num_species = 1;
    unsigned num_voxels = 10;
    vector<double> domain_bounds = {0, 1.0};
    string bc = "reflective";
    double diff = 0.1;
    double end_time = 5.0;
    double time_step = 0.01;

    // Declare a pointer for the simulation object
    Simulation1d sim(num_runs, num_species, num_voxels, domain_bounds, bc);

    // Setup the number of molecules
    sim.SetVoxels({0}, 1000, 0);

    // Setup the reaction rates (diffusion, decay and production)
    sim.SetDiffusionRate(get_jump_rates(sim.GetVoxelDims()), diff, 0);

    // Check whether this file name already exists and alter it appropriately
    sim.Run(sim_name, end_time, time_step);

    return 0;
}

