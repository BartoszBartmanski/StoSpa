//
// Created by bartosz on 28/12/18.
//

#include <iostream>
#include "docopt.h"
#include "Parameters.hpp"
#include "Simulator.hpp"
#include "Decay.hpp"
#include "Production.hpp"

static const char USAGE[] =
        R"(Running stochastic simulations.
If couple of inputs are necessary for one argument, separate them by a comma.

    Usage:
      simple_diffusion [options]

    Options:
      -h --help                                     Show this screen.
      -d --num_dims=<num_dims>                      The number of dimensions of the system. [default: 1]
      --end_time=<end_time>                         Time until to run the simulation [default: 1.0].
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
      --decay=<decay>                               Rate of decay [default: 0.0].
      --prod=<prod>                                 Rate of production [default: 0.0].
      --initial_num=<num>                           Initial number of molecules [default: 1000].
      --initial_pos=<pos>                           Initial position (in indices) of the molecules.

)";


int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);
    Parameters p(args);
    p.SetComments("Data for the simulation.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumSpecies(1);
    p.Add("data_type", "molecules");

    // Declare the simulation name
    string sim_name = "sim_" + to_string(p.GetNumDims()) + "d";
    if (args["--append"]) { sim_name = "_" + args["--append"].asString(); }

    // Get the initial position of the molecules
    vector<unsigned> initial_pos;
    if (args["--initial_pos"])
    {
        initial_pos = split<unsigned>(args["--initial_pos"].asString());
        p.Add("initial_pos", args["--initial_pos"].asString());
    }

    // Declare a pointer for the simulation object
    auto sim = simulator(p);

    // Setup the number of molecules
    if (initial_pos.empty()) { initial_pos = floor_div(sim->GetNumVoxels(), 2); }
    sim->SetInitialNumMolecules(initial_pos, p.GetInitialNum()[0], 0);

    // Setup the reaction rates (diffusion, decay and production)
    sim->SetDiffusionRate(p.GetDiff()[0], 0);
    sim->AddReaction(make_unique<Decay>(p.GetDecay()[0], 0));
    sim->AddReaction(make_unique<Production>(p.GetProd()[0], 0));

    // Check whether this file name already exists and alter it appropriately
    string path_to_file = update_path(p.GetSaveDir(), sim_name, p.GetStartIndex());
    p.Save(path_to_file);
    sim->Run(path_to_file, p.GetEndTime(), p.GetTimeStep());

    return 0;
}
