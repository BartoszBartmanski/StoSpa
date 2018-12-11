#include <iostream>
#include <algorithm>
#include <Production.hpp>
#include "docopt.h"
#include "SimFunctions.hpp"
#include "Parameters.hpp"
#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"
#include "TwoSpeciesDecay.hpp"

static const char USAGE[] =
        R"(Produces the stationary distribution for the following reaction kinetics: A + B -> B and 0 -> A.

    Usage:
      TwoSpeciesDecayDist [options]
      TwoSpeciesDecayDist -h | --help

    Options:
      -h --help                                     Show this screen.
      --info                                        Show all the parameters.
      --end_time=<end_time>                         Time until to run the simulation [default: 500000].
      --time_step=<time_step>                       Length of the time step [default: 0.5].
      -d --num_dims=<num_dims>                      Number of dimensions of the domain [default: 1]
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 1]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,1.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               Value of alpha [default: 0.0].
      --beta=<beta>                                 Value of beta [default: 0.0,0.0].
      --D_a=<D_a>                                   The diffusion coefficient for species A [default: 1.0].
      --D_b=<D_b>                                   The diffusion coefficient for species B [default: 1.0].
      --k_1=<k_1>                                   The rate constant for the two species decay reaction. [default: 0.2]
      --k_2=<k_2>                                   The production rate of species A. [default: 1.0]
      --initial_num=<val>                           Initial number of molecules.
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for a stationary distribution for two-species-decay reaction given the following parameters.");
    params.SetCommand(arr_to_str(argc, argv));
    params.SetNumThreads(1);

    // Get the save directory
    string dir_name = SAVE_DIR;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }

    // Get the file name
    string file_name = "two_species_decay_dist";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Get the start index of the file extension
    auto start_index = unsigned(stoi(args["--start_index"].asString()));

    unsigned num_runs = 1;
    params.SetNumRuns(num_runs);

    auto num_dims = unsigned(stoi(args["--num_dims"].asString()));
    params.SetNumDims(num_dims);

    unsigned num_species = 2;
    params.SetNumSpecies(num_species);

    // Get the number of voxels
    auto num_voxels = unsigned(stoi(args["--num_voxels"].asString()));
    params.SetNumVoxels(num_voxels);

    // Get the domain bounds
    vector<double> domain_bounds = split_stod(args["--domain_bounds"].asString());
    params.SetDomainBounds(domain_bounds);

    // Get the numerical method that will be used to derive the jump coefficients
    string num_method = args["--num_method"].asString();
    params.SetNumMethod(num_method);

    // Get the boundary condition for the domain
    string bc = args["--bc"].asString();
    params.SetBC(bc);

    // Get the end time of the simulation
    double end_time = stod(args["--end_time"].asString());
    params.SetEndTime(end_time);

    // Get the time step for the simulation
    double time_step = stod(args["--time_step"].asString());
    params.SetTimeStep(time_step);

    // Get the voxel aspect ratio
    double kappa = unsigned(stod(args["--kappa"].asString())*num_voxels)/double(num_voxels);
    params.SetKappa(kappa);

    // Get the value of alpha
    double alpha = stod(args["--alpha"].asString());
    params.SetAlpha(alpha);

    // Get the value of beta
    vector<double> beta = split_stod(args["--beta"].asString());
    if (beta.size() == 1) { beta.push_back(beta[0]); }
    params.SetBeta(beta);

    // Get the diffusion coefficients and reaction rates
    double D_a = stod(args["--D_a"].asString());                        // diffusion of species 0
    double D_b = stod(args["--D_b"].asString());                        // diffusion of species 1
    params.SetDiff({D_a, D_b});

    double k_2 = stod(args["--k_2"].asString());                        // production of species 0
    auto prod = make_shared<Production>(k_2, 0);
    params.SetProd({prod->GetRateConstant(), 0});

    double k_1 = stod(args["--k_1"].asString());                        // decay of species 0
    auto two_species_decay = make_shared<TwoSpeciesDecay>(k_1);
    params.AddAdditionalReactions(two_species_decay->GetReactionName(), two_species_decay->GetRateConstant());

    auto mean_a = unsigned(k_2 * pow(domain_bounds[1] - domain_bounds[0], 2) / k_1);
    params.Add("means", to_string(mean_a));

    // Get the initial number of molecules
    vector<unsigned> initial_num;
    if (args["--initial_num"]) { initial_num = split_stou(args["--initial_num"].asString()); }
    else { initial_num = {mean_a/num_voxels, 1}; }
    params.SetInitialNumMolecules(initial_num);

    // Calculate the number of steps that will need to be taken
    auto num_steps = unsigned(end_time / time_step);

    // Create a vector to store the values
    vector<unsigned> stationary_dist(5 * mean_a, 0);

    // Create a pointer to an AbstractSimulation
    unique_ptr<AbstractSimulation> sim;

    // Point the pointer to the object
    if (num_dims == 1)
    {
        sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
    }
    else
    {
        sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
    }

    // Setup the number of molecules
    sim->SetInitialState(initial_num[0] * ones(sim->GetNumVoxels()), 0);
    sim->SetInitialNumMolecules({0, 0}, initial_num[1], 1);

    // Setup the reaction rates
    sim->SetDiffusionRate(D_a, 0);
    sim->SetDiffusionRate(D_b, 1);
    sim->AddReaction(prod);
    sim->AddReaction(two_species_decay);

    // Initialise progress object
    Progress p(num_steps);

    for (unsigned i = 0; i < num_steps; i++)
    {
        sim->Advance(time_step, i);
        unsigned total = sim->GetTotalMolecules();
        if (total < stationary_dist.size())
        {
            stationary_dist[total] += 1;
        }

        // Show progress of the simulation
        p.Show();
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(dir_name, file_name, start_index);

    // Save the information about this plot in the same file as the data
    params.Add("row 1", "Number of times that the system has been in the state given by the index");
    params.Save(path_to_file);

    // Save the data
    save_vector(stationary_dist, path_to_file);

    // Add the simulation name to the log file
    cout << "Data saved in " << path_to_file << endl;

    return 0;
}
