#include <iostream>
#include <algorithm>
#include "docopt.h"
#include "SimFunctions.hpp"
#include "Parameters.hpp"
#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"
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
      MichaelisMentinKinetics -h | --help

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
      --Diff=<Diff>                                 The diffusion coefficients [default: 1.0,1.0,1.0,1.0].
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

void Run(const unique_ptr<AbstractSimulation>& sim, const string& path_to_file, unsigned num_steps, double time_step)
{
    unsigned num_species = sim->GetNumSpecies();

    // Initialise progress object
    Progress p(num_steps);

    for (unsigned i = 0; i < num_steps; i++)
    {
        sim->Advance(time_step, i);

        for (unsigned species=0; species < num_species; species++)
        {
            // Save the stochastic simulation
            save_vector(sim->GetAverageNumMolecules(species), path_to_file);
        }

        // Show progress of the simulation
        p.Show();
    }

    // Add the simulation name to the log file
    cout << "Data saved in " << path_to_file << endl;
}

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for a Michaelis-Mentin kinetics simulation.");
    params.SetCommand(arr_to_str(argc, argv));
    params.SetNumThreads(1);

    // Get the save directory
    string dir_name = SAVE_DIR;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }

    // Get the start index of the file extension
    auto start_index = unsigned(stoi(args["--start_index"].asString()));

    unsigned num_runs = unsigned(stoi(args["--num_runs"].asString()));
    params.SetNumRuns(num_runs);

    auto num_dims = unsigned(stoi(args["--num_dims"].asString()));
    params.SetNumDims(num_dims);

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
    if (num_dims == 1) { kappa = 0.0; }
    params.SetKappa(kappa);

    // Get the value of alpha
    double alpha = stod(args["--alpha"].asString());
    params.SetAlpha(alpha);

    // Get the value of beta
    vector<double> beta = split_stod(args["--beta"].asString());
    if (beta.size() == 1) { beta.push_back(beta[0]); }
    params.SetBeta(beta);

    // Get the initial number of molecules
    vector<unsigned> initial_num = split_stou(args["--initial_num"].asString());

    // Get the diffusion coefficients and reaction rates
    vector<double> Diff = split_stod(args["--Diff"].asString());
    params.SetDiff(Diff);

    // Set up the reactions
    double k_0 = stod(args["--k_0"].asString());
    params.Add("k_0", to_string(k_0));

    double k_1 = stod(args["--k_1"].asString());
    params.Add("k_1", to_string(k_1));

    double k_2 = stod(args["--k_2"].asString());
    params.Add("k_2", to_string(k_2));

    double k_3 = stod(args["--k_3"].asString());
    params.Add("k_3", to_string(k_3));

    unsigned e_t = initial_num[0];
    params.Add("e_t", to_string(e_t));

    // Calculate the number of steps that will need to be taken
    auto num_steps = unsigned(end_time / time_step);

    // Create a pointer to an AbstractSimulation
    unique_ptr<AbstractSimulation> sim;
    string file_name;
    string path_to_file;
    unsigned num_species;

    /*
     * First we simulate the full system.
     */

    // Get the file name
    file_name = "michaelis_mentin_full";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    num_species = 4;
    params.SetNumSpecies(num_species);

    // Point the pointer to the object
    if (num_dims == 1)
    {
        sim = make_unique<Simulation_1d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
    }
    else
    {
        sim = make_unique<Simulation_2d>(num_runs, num_species, num_method, num_voxels, domain_bounds, bc, kappa, alpha, beta[0], beta[1]);
    }

    for (unsigned i=0; i < num_species; i++)
    {
        // Setup the number of molecules
        sim->SetInitialState(initial_num[i] * ones(sim->GetNumVoxels()), i);
        // Setup the reaction rates
        sim->SetDiffusionRate(Diff[i], i);
    }

    auto e_s_to_c = make_shared<MichaelisMentin_I>(k_1);
    sim->AddReaction(e_s_to_c);

    auto c_to_e_s = make_shared<MichaelisMentin_II>(k_2);
    sim->AddReaction(c_to_e_s);

    auto c_to_e_p = make_shared<MichaelisMentin_III>(k_3);
    sim->AddReaction(c_to_e_p);

    auto prod_full = make_shared<Production>(k_0, 1);
    sim->AddReaction(prod_full);
    params.SetProd({0, k_0, 0, 0});

    // Check that the directory exists and that the no file is being over-written.
    path_to_file = update_path(dir_name, file_name, start_index);

    // Save the information about this plot in the same file as the data
    params.Save(path_to_file);

    Run(sim, path_to_file, num_steps, time_step);

    /*
     * Now we simulate the reduced system.
     */

    // Get the file name
    file_name = "michaelis_mentin_reduced";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Update the number of species
    num_species = 2;
    params.SetNumSpecies(num_species);

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
    sim->SetInitialState(initial_num[1] * ones(sim->GetNumVoxels()), 0);
    sim->SetInitialState(initial_num[3] * ones(sim->GetNumVoxels()), 1);
    // Setup the reaction rates
    sim->SetDiffusionRate(Diff[1], 0);
    sim->SetDiffusionRate(Diff[3], 1);

    auto prod_reduced = make_shared<Production>(k_0, 0);
    params.SetProd({k_0, 0});
    sim->AddReaction(prod_reduced);

    double k_m = (k_2 + k_3) / k_1;
    auto s_to_p = make_shared<MichaelisMentinReduced>(k_3, e_t, k_m);
    sim->AddReaction(s_to_p);

    // Check that the directory exists and that the no file is being over-written.
    path_to_file = update_path(dir_name, file_name, start_index);

    // Save the information about this plot in the same file as the data
    params.Save(path_to_file);

    Run(sim, path_to_file, num_steps, time_step);
}
