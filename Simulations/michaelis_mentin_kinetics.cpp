#include <iostream>
#include "docopt.h"
#include "Parameters.hpp"
#include "Simulator.hpp"
#include "Production.hpp"
#include "MichaelisMentin.hpp"
#include "MichaelisMentinReduced.hpp"

static const char USAGE[] =
        R"(Running Michealis-Mentin kinetics E + S [k_2]<->[k_1] C ->[k_3] E + P and
    the reduced kinetics S ->[k_3 E_t S / (S + k_m)] P
    Indexes of the species:
        Species E - index 0
        Species S - index 1
        Species C - index 2
        Species P - index 3

    Usage:
      MichaelisMentinKinetics [options]

    Options:
      -h --help                                     Show this screen.
      --info                                        Show all the parameters.
      --end_time=<end_time>                         Time until to run the simulation [default: 500].
      --time_step=<time_step>                       Length of the time step [default: 0.5].
      -d --num_dims=<num_dims>                      Number of dimensions of the domain [default: 1]
      --num_runs=<num_runs>                         Number of runs [default: 1].
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 1]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,100.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               Value of alpha [default: 0.0].
      --beta=<beta>                                 Value of beta [default: 0.0,0.0].
      --diff=<diff>                                 The diffusion coefficients [default: 1.0,1.0,1.0,1.0].
      --k_0=<k_0>                                   0 -> S reaction rate [default: 0.25].
      --k_1=<k_1>                                   E + S -> C reaction rate [default: 1.0].
      --k_2=<k_2>                                   C -> E + S reaction rate [default: 1.0].
      --k_3=<k_3>                                   C -> E + P reaction rate [default: 5.0].
      --initial_num=<val>                           Initial number of molecules [default: 10,0,0,0].
      --mean_and_var                                Output just the mean and the variance for the data.
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);
    Parameters p(args);
    p.SetComments("Data for a Michaelis-Mentin kinetics simulation.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumThreads(1);

    // Set up the reactions
    double k_0 = stod(args["--k_0"].asString());
    p.Add("k_0", to_string(k_0));

    double k_1 = stod(args["--k_1"].asString());
    p.Add("k_1", to_string(k_1));

    double k_2 = stod(args["--k_2"].asString());
    p.Add("k_2", to_string(k_2));

    double k_3 = stod(args["--k_3"].asString());
    p.Add("k_3", to_string(k_3));

    unsigned e_t = p.GetInitialNum()[0];
    p.Add("e_t", to_string(e_t));

    // Get the file name
    string file_name_full = "michaelis_mentin_full";
    string file_name_reduced = "michaelis_mentin_reduced";
    if (args["--append"])
    {
        file_name_full += ("_" + args["--append"].asString());
        file_name_reduced += ("_" + args["--append"].asString());
    }

    /*
     * First we setup the full system.
     */

    // Point the pointer to the object
    p.SetNumSpecies(4);
    auto sim_full = simulator(p);

    for (unsigned i=0; i < p.GetNumSpecies(); i++)
    {
        // Setup the number of molecules
        sim_full->SetVoxels(p.GetInitialNum()[i] * ones(sim_full->GetNumVoxels()), i);
        // Setup the reaction rates
        sim_full->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[i], i);
    }

    sim_full->AddReaction(make_unique<MichaelisMentin_I>(k_1));
    sim_full->AddReaction(make_unique<MichaelisMentin_II>(k_2));
    sim_full->AddReaction(make_unique<MichaelisMentin_III>(k_3));
    sim_full->AddReaction(make_unique<Production>(k_0, 1));
    p.SetProd({0, k_0, 0, 0});

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file_full = update_path(p.GetSaveDir(), file_name_full, p.GetStartIndex());
    p.Save(path_to_file_full);
    sim_full->Run(path_to_file_full, p.GetEndTime(), p.GetTimeStep());

    /*
     * Now we setup the reduced system.
     */

    // Point the pointer to the object
    p.SetNumSpecies(2);
    auto sim_reduced = simulator(p);

    // Setup the number of molecules
    sim_reduced->SetVoxels(p.GetInitialNum()[1] * ones(sim_reduced->GetNumVoxels()), 0);
    sim_reduced->SetVoxels(p.GetInitialNum()[3] * ones(sim_reduced->GetNumVoxels()), 1);
    // Setup the reaction rates
    sim_reduced->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[1], 0);
    sim_reduced->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[3], 1);

    sim_reduced->AddReaction(make_unique<Production>(k_0, 0));

    double k_m = (k_2 + k_3) / k_1;
    sim_reduced->AddReaction(make_unique<MichaelisMentinReduced>(k_3, e_t, k_m));

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file_reduced = update_path(p.GetSaveDir(), file_name_reduced, p.GetStartIndex());
    p.Save(path_to_file_reduced);        // Setup the reaction rates
    sim_reduced->Run(path_to_file_reduced, p.GetEndTime(), p.GetTimeStep());

    return 0;
}
