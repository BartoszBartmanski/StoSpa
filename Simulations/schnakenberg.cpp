//
// Created by bartosz on 30/12/18.
//

#include <iostream>
#include "docopt.h"
#include "Parameters.hpp"
#include "Simulator.hpp"
#include "Decay.hpp"
#include "Production.hpp"
#include "Schnakenberg.hpp"

static const char USAGE[] =
        R"(Running stochastic simulations.
If couple of inputs are necessary for one argument, separate them by a comma.

    Usage:
      schnakenberg [options]

    Options:
      -h --help                                     Show this screen.
      -d --num_dims=<num_dims>                      The number of dimensions of the system. [default: 1]
      --end_time=<end_time>                         Time until to run the simulation [default: 1800.0].
      --time_step=<time_step>                       Length of the time step [default: 10.0].
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_runs=<num_runs>                         Number of realisations [default: 1].p.GetTimeStep()
      --num_voxels=<num_voxels>                     Number of voxels. [default: 40]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,1.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               The value of alpha [default: 0.0].
      --beta=<beta>                                 The value of beta [default: 0.0,0.0].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --diff=<diff>                                 Diffusion coefficients [default: 0.00001,0.001].
      --decay=<decay>                               Rate of decay [default: 0.02,0.0].
      --prod=<prod>                                 Rate of production [default: 1.0,3.0].
      --rate=<rate>                                 Rate of schnakenberg reaction [default: 0.000001]
      --initial_num=<num>                           Initial number of molecules [default: 200,75].
      --growth_rate=<rate>                          Growth rate of the domain. [default: 0.0]
)";


int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "StoSpa 0.1");
    Parameters p(args);
    p.SetComments("Data for the simulation.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumSpecies(2);
    p.Add("data_type", "molecules");
    double rate = stod(args["--rate"].asString());

    // Declare the simulation name
    string sim_name = "schnakenberg_" + to_string(p.GetNumDims()) + "d";
    if (args["--append"]) { sim_name = "_" + args["--append"].asString(); }

    // Initialise the simulation object
    auto sim = simulator(p);

    // Setup the number of molecules
    sim->SetVoxels(p.GetInitialNum()[0] * ones(sim->GetNumVoxels()), 0);
    sim->SetVoxels(p.GetInitialNum()[1] * ones(sim->GetNumVoxels()), 1);

    // Setup the reaction rates
    sim->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[0], 0);
    sim->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[1], 1);
    sim->AddReaction(make_unique<Decay>(p.GetDecay()[0], 0));

    p.SetProd(p.GetProd()/sim->GetVoxelSize());
    sim->AddReaction(make_unique<Production>(p.GetProd()[0], 0));
    sim->AddReaction(make_unique<Production>(p.GetProd()[1], 1));

    unique_ptr<AbstractReaction> schnakenberg = make_unique<Schnakenberg>(rate * pow(sim->GetVoxelSize(), 2));
    p.AddReaction(schnakenberg);
    sim->AddReaction(move(schnakenberg));

    double growth = stod(args["--growth_rate"].asString());
    if (growth > 0)
    {
        p.Add("growth_rate", args["--growth_rate"].asString());
        sim->SetGrowth(make_unique<Exponential>(p.GetNumDims(), growth));
    }

    // Check whether this file name already exists and alter it appropriately
    string path_to_file = update_path(p.GetSaveDir(), sim_name, p.GetStartIndex());
    p.Save(path_to_file);
    sim->Run(path_to_file, p.GetEndTime(), p.GetTimeStep());

    return 0;
}
