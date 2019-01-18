//
// Created by bartosz on 11/12/18.
//

#include <iostream>
#include "docopt.h"
#include "Utilities.hpp"
#include "Parameters.hpp"
#include "Simulator.hpp"
#include "GrowthRates.hpp"


static const char USAGE[] =
        R"(Running stochastic simulations.
If couple of inputs are necessary for one argument, separate them by a comma.

    Usage:
      growing_domain [options]

    Options:
      -h --help                                     Show this screen.
      -d --num_dims=<num_dims>                      The number of dimensions of the system. [default: 1]
      --end_time=<end_time>                         Time until to run the simulation [default: 5.0].
      --time_step=<time_step>                       Length of the time step [default: 0.01].
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_runs=<num_runs>                         Number of realisations [default: 1].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 21]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,20.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               The value of alpha [default: 0.0].
      --beta=<beta>                                 The value of beta [default: 0.0,0.0].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --diff=<diff>                                 Diffusion coefficient [default: 1.0].
      --initial_num=<num>                           Initial number of molecules [default: 1000].
      --growth_rate=<rate>                          Growth rate of the domain. [default: 1.0]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters p(args);
    p.SetComments("Data for the simulation.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumSpecies(1);
    p.Add("data_type", "molecules");
    p.Add("growth_rate", args["--growth_rate"].asString());
    double growth = stod(args["--growth_rate"].asString());

    // Declare the simulation name
    string sim_name = "growing_domain";
    if (args["--append"]) { sim_name = "_" + args["--append"].asString(); }

    // Declare a pointer for the simulation object
    Simulation1d sim(p);

    sim.SetVoxels({0}, p.GetInitialNum()[0], 0);

    sim.SetDiffusionRate(get_jump_rates(p), p.GetDiff()[0], 0);

    sim.SetGrowth(make_unique<Exponential>(growth));

    string path_to_file = update_path(p.GetSaveDir(), sim_name, p.GetStartIndex());  // Get appropriate filename
    p.Save(path_to_file);  // Save parameters
    sim.Run(path_to_file, p.GetEndTime(), p.GetTimeStep());

    return 0;
}