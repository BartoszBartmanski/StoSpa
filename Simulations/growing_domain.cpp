//
// Created by bartosz on 11/12/18.
//

#include <iostream>
#include "docopt.h"
#include "SimFunctions.hpp"
#include "Parameters.hpp"
#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"

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
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for the simulation.");
    params.SetCommand(arr_to_str(argc, argv));

    // Get the save directory
    string dir_name = SAVE_DIR;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }

    // Declare the simulation name
    string sim_name = "growing_domain";
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

    // Get the numerical method of derivation for the jump coefficients
    string num_method = args["--num_method"].asString();
    params.SetNumMethod(num_method);

    // Get the boundary condition
    string bc = "Exponential";
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
    double diff = stod(args["--diff"].asString());
    params.SetDiff(diff);

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
    unsigned num_species = 1;
    params.SetNumSpecies(num_species);

    // Calculate the voxel size
    double voxel_size = get_voxel_size(domain_bounds, num_voxels, kappa);

    // The data saved will be number of molecules, so save this as well
    params.Add("data_type", "molecules");

    // Initialise the both number of runs variables
    unsigned num_runs = unsigned(stoi(args["--num_runs"].asString()));
    params.SetNumRuns(num_runs);

    // Declare a pointer for the simulation object
    Simulation_1d sim = Simulation_1d(num_runs, num_species, num_method, num_voxels, domain_bounds, bc);
    sim.UseExtrande();

    vector<unsigned> initial_pos = floor_div(sim.GetNumVoxels(), 2);
    sim.SetInitialNumMolecules(initial_pos, initial_num[0], 0);

    sim.SetDiffusionRate(diff, 0);

    string path_to_file = update_path(dir_name, sim_name, start_index);
    params.Save(path_to_file);

    unique_ptr<ofstream> output = make_unique<ofstream>(path_to_file, ios::app);

    // Initialise progress object
    Progress p(num_steps);

    // Run the SSA
    for (unsigned i=0; i < num_steps; i++)
    {
        // Move to the next time step
        sim.Advance(time_step, i);

        for (unsigned species=0; species < num_species; species++)
        {
            // Save the stochastic simulation
            save_vector(sim.GetAverageNumMolecules(species), output);
        }

        // Print the progress of the simulation
        p.Show();
    }

    cout << "Data saved in " << path_to_file << endl;


}