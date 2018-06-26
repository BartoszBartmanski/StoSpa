#include <iostream>
#include "docopt.h"
#include "SimFunctions.hpp"
#include "Parameters.hpp"
#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"
#include "Decay.hpp"
#include "Production.hpp"
#include "Schnakenberg.hpp"
#include "Dimerisation.hpp"
#include "TwoSpeciesDecay.hpp"

static const char USAGE[] =
        R"(Running stochastic simulations.
If couple of inputs are necessary for one argument, separate them by a comma.

    Usage:
      StoSpa sim [--initial_num=<num>] [--initial_pos <pos>] [options]
      StoSpa schnakenberg [options]
      StoSpa two_species_decay [--initial_num=<val>] [--k_1=<k>] [--k_2=<k>] [options]
      StoSpa dimerisation [options]
      StoSpa -h | --help
      StoSpa --version

    Options:
      -h --help                                     Show this screen.
      --version                                     Show version.
      --info                                        Show all the parameters.
      -d --num_dims=<num_dims>                      The number of dimensions of the system. [default: 1]
      --end_time=<end_time>                         Time until to run the simulation [default: 1.0].
      --time_step=<time_step>                       Length of the time step [default: 0.01].
      --num_method=<num_method>                     The method of derivation of jump coefficients [default: fem].
      --num_runs=<num_runs>                         Number of realisations [default: 1].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 40]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,1.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.0].
      --alpha=<alpha>                               The value of alpha [default: 0.0].
      --beta=<beta>                                 The value of beta [default: 0.0,0.0].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --average                                     Whether to average the simulations straight away.
      --diff=<diff>                                 Diffusion coefficient [default: 1.0,1.0].
      --decay=<decay>                               Rate of decay [default: 0.0,0.0].
      --prod=<prod>                                 Rate of production [default: 0.0,0.0].
      --initial_num=<num>                           Initial number of molecules [default: 1000,1].
      --initial_pos                                 Initial position (in indices) of the molecules.
      --k_1=<k>                                     Rate of two-species-decay reaction [default: 0.2].
      --k_2=<k>                                     Rate of production of species A [default: 1.0].
)";

void Run(const unique_ptr<AbstractSimulation>& sim, const string& path_to_file, unsigned num_steps, double time_step)
{
    unsigned num_species = sim->GetNumSpecies();

    // Create a file handle
    unique_ptr<ofstream> output = make_unique<ofstream>(path_to_file, ios::app);

    // Initialise progress object
    Progress p(num_steps+1);

    // Run the SSA
    for (unsigned i=0; i <= num_steps; i++)
    {
        // Move to the next time step
        sim->Advance(time_step, i);

        for (unsigned species=0; species < num_species; species++)
        {
            // Save the stochastic simulation
            save_vector(sim->GetAverageNumMolecules(species), output);
        }

        // Print the progress of the simulation
        p.Show(i);
    }

    cout << "Data saved in " << path_to_file << endl;
}

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true, "StoSpa 0.1");

    // Show the raw command line input if prompted
    if (args["--info"].asBool())
    {
        for (auto const &arg : args)
        {
            std::cout << arg.first << ": " << arg.second << std::endl;
        }
    }

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for the simulation.");
    params.SetCommand(arr_to_str(argc, argv));
    params.SetNumThreads(1);

    // Get the save directory
    string dir_name;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }
    else { dir_name = get_dir(); }

    // Declare the simulation name
    string sim_name;
    if (args["--append"]) { sim_name = "_" + args["--append"].asString(); }

    // Get the start index to be appended at the end of simulation file name
    auto start_index = unsigned(stoi(args["--start_index"].asString()));

    // Get number of dimensions of the domain
    unsigned num_dims = unsigned(stoi(args["--num_dims"].asString()));
    params.SetNumDims(num_dims);

    // Get the number of voxels in the x-direction
    unsigned num_voxels = unsigned(stoi(args["--num_voxels"].asString()));
    params.SetNumVoxels(num_voxels);

    // Get the domain bounds
    vector<double> domain_bounds = split_stod(args["--domain_bounds"].asString());
    params.SetDomainBounds(domain_bounds);

    // Get initial number of molecules
    vector<unsigned> initial_num = split_stou(args["--initial_num"].asString());
    params.SetInitialNumMolecules(initial_num);

    // Get the initial position of the molecules
    vector<unsigned> initial_pos;
    if (args["--initial_pos"].asBool())
    {
        initial_pos = split_stou(args["<pos>"].asString());
        params.Add("initial_pos", args["<pos>"].asString());
    }

    // Get the numerical method of derivation for the jump coefficients
    string num_method = args["--num_method"].asString();
    params.SetNumMethod(num_method);

    // Get the boundary condition
    string bc = args["--bc"].asString();
    params.SetBC(bc);

    // Get the end time of the simulation
    double end_time = stod(args["--end_time"].asString());
    params.SetEndTime(end_time);

    // Get the time step
    double time_step = stod(args["--time_step"].asString());
    params.SetTimeStep(time_step);

    // Calculate the number of steps that will need to be taken to reach the desired point in time
    auto num_steps = unsigned(end_time/time_step);

    // Get the diffusion coefficients
    vector<double> diff = split_stod(args["--diff"].asString());
    params.SetDiff(diff);

    // Get the decay constants
    vector<double> decay = split_stod(args["--decay"].asString());
    params.SetDecay(decay);

    // Get the production constants
    vector<double> prod = split_stod(args["--prod"].asString());
    params.SetProd(prod);

    // Get the voxel aspect ratio - kappa
    double kappa = stod(args["--kappa"].asString());
    if (num_dims == 1) { kappa = 0; }
    params.SetKappa(kappa);

    // Get the FDM derivation parameter - alpha
    double alpha = stod(args["--alpha"].asString());
    params.SetAlpha(alpha);

    // Get the FET derivation parameter - beta
    vector<double> beta = split_stod(args["--beta"].asString());
    if (beta.size() == 1) { beta.push_back(beta[0]); }
    params.SetBeta(beta);

    // Declare the num_species variable, which will be set later
    unsigned num_species;

    // Calculate the voxel size
    double voxel_size = get_voxel_size(domain_bounds, num_voxels, kappa);

    // The data saved will be number of molecules, so save this as well
    params.Add("data_type", "molecules");

    // Initialise the both number of runs variables
    unsigned num_runs = 1;
    unsigned num_separate_runs = 1;

    // Determine whether the simulations will be averaged or save separately
    if (args["--average"].asBool()) { num_runs = unsigned(stoi(args["--num_runs"].asString())); }
    else { num_separate_runs = unsigned(stoi(args["--num_runs"].asString())); }
    params.SetNumRuns(num_runs);

    // Declare a pointer for the simulation object
    unique_ptr<AbstractSimulation> sim;

    if (args["sim"].asBool())
    {
        sim_name.insert(0, "sim_" + to_string(num_dims) + "d");

        num_species = 1;
        params.SetNumSpecies(num_species);

        auto decay_reaction = make_shared<Decay>(decay[0], 0);
        params.SetDecay(decay[0]);

        auto prod_reaction = make_shared<Production>(prod[0], 0);
        params.SetProd(prod[0]);
        
        for (unsigned i = 0; i < num_separate_runs; i++)
        {
            if (num_dims == 1)
            {
                sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
            }
            else
            {
                sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
            }

            // Check that the voxel size has been correctly calculated
            assert(sim->GetVoxelSize() == voxel_size);

            // Setup the number of molecules
            if (initial_pos.empty()) { initial_pos = floor_div(sim->GetNumVoxels(), 2); }
            sim->SetInitialNumMolecules(initial_pos, initial_num[0], 0);

            // Setup the reaction rates (diffusion, decay and production)
            sim->SetDiffusionRate(diff[0], 0);
            sim->AddReaction(decay_reaction);
            sim->AddReaction(prod_reaction);

            // Check whether this file name already exists and alter it appropriately
            string path_to_file = update_path(dir_name, sim_name, start_index);
            params.Save(path_to_file);

            // Run the simulation
            Run(sim, path_to_file, num_steps, time_step);
        }
    }
    else if (args["schnakenberg"].asBool())
    {
        sim_name.insert(0, "schnakenberg_" + to_string(num_dims) + "d");

        num_species = 2;
        params.SetNumSpecies(num_species);

        double D_u = 0.00001;                                                   // diffusion of species 0
        double D_v = 0.001;                                                     // diffusion of species 1
        params.SetDiff({D_u, D_v});        

        // Create the additional reaction object
        auto schnakenberg = make_shared<Schnakenberg>(0.000001 * pow(voxel_size, 2));
        params.AddAdditionalReactions(schnakenberg->GetReactionName(), schnakenberg->GetRateConstant());

        auto decay_reaction = make_shared<Decay>(0.02, 0);
        params.SetDecay({decay_reaction->GetRateConstant(), 0});
        
        auto prod_species_0 = make_shared<Production>(1.0 / voxel_size, 0);
        auto prod_species_1 = make_shared<Production>(3.0 / voxel_size, 1);
        params.SetProd({prod_species_0->GetRateConstant(), prod_species_1->GetRateConstant()});

        // Initialise the steady state values
        initial_num = {200, 75};
        params.SetInitialNumMolecules(initial_num);

        for (unsigned i = 0; i < num_separate_runs; i++)
        {
            if (num_dims == 1)
            {
                sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
            }
            else
            {
                sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
            }

            // Check that the voxel size has been correctly calculated
            assert(sim->GetVoxelSize() == voxel_size);

            // Setup the number of molecules
            sim->SetInitialState(initial_num[0] * ones(sim->GetNumVoxels()), 0);
            sim->SetInitialState(initial_num[1] * ones(sim->GetNumVoxels()), 1);

            // Setup the reaction rates
            sim->SetDiffusionRate(D_u, 0);
            sim->SetDiffusionRate(D_v, 1);
            sim->AddReaction(decay_reaction);
            sim->AddReaction(prod_species_0);
            sim->AddReaction(prod_species_1);
            sim->AddReaction(schnakenberg);

            // Check whether this file name already exists and alter it appropriately
            string path_to_file = update_path(dir_name, sim_name, start_index);
            params.Save(path_to_file);

            // Run the simulation
            Run(sim, path_to_file, num_steps, time_step);
        }
    }
    else if (args["two_species_decay"].asBool())
    {
        sim_name.insert(0, "two_species_decay_" + to_string(num_dims) + "d");

        num_species = 2;
        params.SetNumSpecies(num_species);

        // Get the production rate of species 0
        double k_2 = stod(args["--k_2"].asString());                        // production of species 0
        auto prod_reaction = make_shared<Production>(k_2, 0);
        params.SetProd({k_2, 0});

        // Get the rate of the two species decay reaction
        double k_1 = stod(args["--k_1"].asString());                        // decay of species 0
        auto two_species_decay = make_shared<TwoSpeciesDecay>(k_1);
        params.AddAdditionalReactions(two_species_decay->GetReactionName(), two_species_decay->GetRateConstant());

        for (unsigned i = 0; i < num_separate_runs; i++)
        {
            if (num_dims == 1)
            {
                sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
            }
            else
            {
                sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
            }

            // Setup the reaction rates (diffusion, decay and production)
            sim->SetDiffusionRate(diff[0], 0);
            sim->SetDiffusionRate(diff[1], 1);
            sim->AddReaction(prod_reaction);
            sim->AddReaction(two_species_decay);

            // Setup the number of molecules
            sim->SetInitialState(initial_num[0] * ones(sim->GetNumVoxels()), 0);
            sim->SetInitialNumMolecules({0, 0}, 1, 1);

            // Check whether this file name already exists and alter it appropriately
            string path_to_file = update_path(dir_name, sim_name, start_index);
            params.Save(path_to_file);

            // Run the simulation
            Run(sim, path_to_file, num_steps, time_step);
        }
    }
    else if (args["dimerisation"].asBool())
    {
        sim_name.insert(0, "dimerisation_" + to_string(num_dims) + "d");

        num_species = 1;
        params.SetNumSpecies(num_species);

        // Rate constants
        double k_1 = 0.1;  // rate of dimerisation A + A -> 0
        auto dimerisation = make_shared<Dimerisation>(k_1, 0);
        params.AddAdditionalReactions(dimerisation->GetReactionName(), dimerisation->GetRateConstant());

        auto prod_reaction = make_shared<Production>(prod[0], 0);
        params.SetProd(prod[0]);

        for (unsigned i = 0; i < num_separate_runs; i++)
        {
            if (num_dims == 1)
            {
                sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
            }
            else
            {
                sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
            }

            // Check that the voxel size has been correctly calculated
            assert(sim->GetVoxelSize() == voxel_size);

            // Setup the number of molecules
            sim->SetInitialState(ones(sim->GetNumVoxels()), 0);

            // Setup the reaction rates
            sim->SetDiffusionRate(diff[0], 0);
            sim->AddReaction(prod_reaction);
            sim->AddReaction(dimerisation);

            // Check whether this file name already exists and alter it appropriately
            string path_to_file = update_path(dir_name, sim_name, start_index);
            params.Save(path_to_file);

            // Run the simulation
            Run(sim, path_to_file, num_steps, time_step);
        }
    }

    return 0;
}
