#include <iostream>
#include <Production.hpp>
#include "docopt.h"
#include "Parameters.hpp"
#include "Simulator.hpp"
#include "TwoSpeciesDecay.hpp"

static const char USAGE[] =
        R"(Produces the stationary distribution for the following reaction kinetics: A + B -> B and 0 -> A.

    Usage:
      TwoSpeciesDecayDist [options]

    Options:
      -h --help                                     Show this screen.
      --end_time=<end_time>                         Time until to run the simulation [default: 500000].
      --time_step=<time_step>                       Length of the time step [default: 0.5].
      -d --num_dims=<num_dims>                      Number of dimensions of the domain [default: 1]
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 50]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,1.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               Value of alpha [default: 0.0].
      --beta=<beta>                                 Value of beta [default: 0.0,0.0].
      --diff=<diff>                                 The diffusion coefficient for species A [default: 1.0,1.0].
      --k_1=<k_1>                                   The rate constant for the two species decay reaction. [default: 0.2]
      --k_2=<k_2>                                   The production rate of species A. [default: 1.0]
      --initial_num=<val>                           Initial number of molecules. [default: 5,1]
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);
    Parameters p(args);
    p.SetComments("Data for a stationary distribution for two-species-decay reaction given the following parameters.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumSpecies(2);

    // Get the file name
    string file_name = "two_species_decay_dist";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    double k_2 = stod(args["--k_2"].asString());                        // production of species 0
    double k_1 = stod(args["--k_1"].asString());                        // decay of species 0

    vector<unsigned> stationary_dist(100 * p.GetInitialNum()[0], 0);

    // Point the pointer to the object
    auto sim = simulator(p);

    // Setup the number of molecules
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<unsigned> rand_x(0, p.GetNumVoxels() - 1);
    uniform_int_distribution<unsigned> rand_y(0, unsigned(p.GetNumVoxels()/p.GetKappa()) - 1);
    for (unsigned i=0; i < p.GetInitialNum()[0]; i++)
    {
        sim->SetVoxels({rand_x(gen), rand_y(gen)}, 1, 0);
    }
    sim->SetVoxels({0, 0}, p.GetInitialNum()[1], 1);

    // Setup the reaction rates
    sim->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[0], 0);
    sim->SetDiffusionRate(get_jump_rates(p), p.GetDiff()[1], 1);

    p.SetProd({k_2, 0});
    sim->AddReaction(make_unique<Production>(k_2, 0));

    unique_ptr<AbstractReaction> two_species_decay = make_unique<TwoSpeciesDecay>(k_1);
    p.AddReaction(two_species_decay);
    sim->AddReaction(move(two_species_decay));

    // Initialise progress object
    auto num_steps = unsigned(p.GetEndTime() / p.GetTimeStep());
    Progress prog(num_steps-1);

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(p.GetSaveDir(), file_name, p.GetStartIndex());
    string path_to_raw_data = path_to_file.substr(0, path_to_file.size()-4)+ "_raw.dat";

    // Save the information about this plot in the same file as the data
    p.Save(path_to_raw_data);
    p.Add("row 1", "Number of times that the system has been in the state given by the index");
    p.Save(path_to_file);

    for (unsigned i = 1; i < num_steps; i++)
    {
        sim->Advance(i*p.GetTimeStep());
        unsigned total = sim->GetTotalMolecules();
        save(total, path_to_raw_data);
        if (total < stationary_dist.size())
        {
            stationary_dist[total] += 1;
        }

        // Show progress of the simulation
        prog.Show();
    }

    // Save the data
    save(stationary_dist, path_to_file);

    // Add the simulation name to the log file
    prog.End(path_to_file);

    return 0;
}
