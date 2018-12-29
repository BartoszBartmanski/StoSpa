//
// Created by bartosz on 11/12/18.
//

#include <iostream>
#include "docopt.h"
#include "Utilities.hpp"
#include "Parameters.hpp"
#include "Simulation_1d.hpp"


static const char USAGE[] =
        R"(Running stochastic simulations.
If couple of inputs are necessary for one argument, separate them by a comma.

    Usage:
      growing_domain [options]
      growing_domain -h | --help

    Options:
      -h --help                                     Show this screen.
      -d --num_dims=<num_dims>                      The number of dimensions of the system. [default: 1]
      --end_time=<end_time>                         Time until to run the simulation [default: 5.0].
      --time_step=<time_step>                       Length of the time step [default: 0.01].
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_runs=<num_runs>                         Number of realisations [default: 1].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 21]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,20.0]
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               The value of alpha [default: 0.0].
      --beta=<beta>                                 The value of beta [default: 0.0,0.0].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --diff=<diff>                                 Diffusion coefficient [default: 1.0].
      --initial_num=<num>                           Initial number of molecules [default: 1000].
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters p(args);
    p.SetComments("Data for the simulation.");
    p.SetCommand(arr_to_str(argc, argv));

    // Declare the simulation name
    string sim_name = "growing_domain";
    if (args["--append"]) { sim_name = "_" + args["--append"].asString(); }

    p.SetBC("Exponential");

    // Calculate the number of steps that will need to be taken to reach the desired point in time
    auto num_steps = unsigned(p.GetEndTime()/p.GetTimeStep());

    // The data saved will be number of molecules, so save this as well
    p.Add("data_type", "molecules");

    // Declare a pointer for the simulation object
    Simulation_1d sim(p.GetNumRuns(), p.GetNumSpecies(), p.GetNumMethod(), p.GetNumVoxels(), p.GetDomainBounds(), p.GetBC());
    sim.UseExtrande();

    sim.SetInitialNumMolecules({0}, p.GetInitialNum()[0], 0);

    sim.SetDiffusionRate(p.GetDiff()[0], 0);

    string path_to_file = update_path(p.GetSaveDir(), sim_name, p.GetStartIndex());
    p.Save(path_to_file);

    unique_ptr<ofstream> output = make_unique<ofstream>(path_to_file, ios::app);

    // Initialise progress object
    Progress prog(num_steps);

    // Run the SSA
    for (unsigned i=0; i < num_steps; i++)
    {
        // Move to the next time step
        sim.Advance(p.GetTimeStep(), i);

        // Save the stochastic simulation
        save_vector(sim.GetAverageNumMolecules(), output);

        // Print the progress of the simulation
        prog.Show();
    }

    cout << "Data saved in " << path_to_file << endl;


}