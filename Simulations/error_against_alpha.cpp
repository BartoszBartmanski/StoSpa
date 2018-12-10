#include <iostream>
#include <future>
#include "docopt.h"
#include "Utilities.hpp"
#include "Parameters.hpp"
#include "SimFunctions.hpp"
#include "DiffEqAnalytic.hpp"
#include "Simulation_2d.hpp"
#include "Decay.hpp"
#include "Production.hpp"

static const char USAGE[] =
        R"(Produces a data for the plot of error against value of alpha.

    Usage:
      error_against_alpha [options]
      error_against_alpha -h | --help

    Options:
      -h --help                                     Show this screen.
      --info                                        Show all the parameters.
      --end_time=<end_time>                         Time until to run the simulation [default: 5.0].
      --num_runs=<num_runs>                         Number of realisations [default: 1].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 21]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,20.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --kappa=<kappa>                               The voxel aspect ratio [default: 1.4].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --diff=<diff>                                 DiffusionReflective coefficient [default: 1.0].
      --decay=<decay>                               Rate of decay [default: 0.0].
      --prod=<prod>                                 Rate of production [default: 0.0].
      --initial_num=<val>                           Initial number of molecules [default: 5000000].
      --num_points=<points>                         Number of points of alpha for which to run the simulation. [default: 50]
      --trunc_order=<trunc>                         The truncation order of analytic solution. [default: 1000]
      -j --num_threads=<num_threads>                Number of threads to be used. [default: 1]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for a plot of values of error against values of alpha in the 2d simulations with alpha being varied.");
    params.SetCommand(arr_to_str(argc, argv));
    params.SetNumDims(2);

    // Get the directory name
    string dir_name;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }
    else { dir_name = get_dir(); }

    // Get the start index of the file extension
    auto start_index = unsigned(stoi(args["--start_index"].asString()));

    // Name the file
    string file_name = "error_against_alpha";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Get the number of runs
    auto num_runs = unsigned(stoi(args["--num_runs"].asString()));
    params.SetNumRuns(num_runs);

    unsigned num_species = 1;
    params.SetNumSpecies(num_species);

    string method = "fdm";
    params.SetNumMethod(method);

    // Get the number of voxels
    auto num_voxels = unsigned(stoi(args["--num_voxels"].asString()));
    params.SetNumVoxels(num_voxels);

    // Get the domain bounds
    vector<double> domain_bounds = split_stod(args["--domain_bounds"].asString());
    params.SetDomainBounds(domain_bounds);

    // Get the initial number of molecules
    auto initial_num = unsigned(stoi(args["--initial_num"].asString()));
    params.SetInitialNumMolecules(initial_num);

    // Get the boundary condition for the domain
    string bc = args["--bc"].asString();
    params.SetBC(bc);

    // Get the end time of the simulation
    double end_time = stod(args["--end_time"].asString());
    params.SetEndTime(end_time);

    // Get the diffusion coefficient
    double diff = stod(args["--diff"].asString());
    params.SetDiff(diff);

    // Get the decay rate
    double decay = stod(args["--decay"].asString());
    params.SetDecay(decay);

    // Get the production rate
    double prod = stod(args["--prod"].asString());
    params.SetProd(prod);

    // Get the voxel aspect ratio
    double kappa = unsigned(stod(args["--kappa"].asString())*num_voxels)/double(num_voxels);
    params.SetKappa(kappa);

    // Get the number of points of alpha for which to run the simulation
    unsigned num_points = unsigned(stoi(args["--num_points"].asString()));
    params.SetNumPoints(num_points);

    // Get the truncation order for the analytic solution
    unsigned trunc_order = unsigned(stoi(args["--trunc_order"].asString()));
    params.SetTruncOrder(trunc_order);

    // Get the number of cores
    unsigned num_threads = unsigned(stoi(args["--num_threads"].asString()));
    num_threads = unsigned(gcd(num_points, num_threads));
    params.SetNumThreads(num_threads);

    // Create a vector of alpha values
    vector<double> alpha = linspace(0.0, 1.0/kappa, num_points);

    // Create a vector to hold the results
    vector<long double> error(num_points, 0);

    // Calculate the initial point of the molecules
    vector<double> mid_point = get_midpoint(domain_bounds, num_voxels, kappa);

    // Calculate the analytical solution
    DiffEqAnalytic analytic = DiffEqAnalytic(2,
                                             end_time,
                                             mid_point,
                                             domain_bounds,
                                             {num_voxels, unsigned(kappa*num_voxels)},
                                             initial_num,
                                             diff,
                                             decay,
                                             prod,
                                             trunc_order);
    vector<double> sol = analytic.GetAnalytic();

    unsigned num_sets = num_points / num_threads;

    // Create vector of futures
    vector<future<double>> futures(num_threads);
    vector<Simulation_2d> sims(num_threads);

    // Initialise progress object
    Progress p(num_points);

    // Collect the data
    for (unsigned set=0; set < num_sets; set++)
    {
        for (unsigned i=0; i < num_threads; i++)
        {
            unsigned k = set * num_threads + i;
            sims[i] = Simulation_2d(num_runs, num_species, method, num_voxels, domain_bounds, bc, kappa, alpha[k], 0.0, 0.0);
            sims[i].SetDiffusionRate(diff, 0);
            sims[i].AddReaction(make_shared<Decay>(decay, 0));
            sims[i].AddReaction(make_shared<Production>(prod, 0));
            sims[i].SetInitialNumMolecules(floor_div(sims[i].GetNumVoxels(), 2), initial_num, 0);
        }

        for (unsigned i=0; i < num_threads; i++)
        {
            futures[i] = async(launch::async, get_error, ref(sims[i]), sol, end_time);
        }

        for (unsigned i=0; i < num_threads; i++)
        {
            unsigned k = set * num_threads + i;
            error[k] = futures[i].get();

            // Show progress of the simulation
            p.Show();
        }
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(dir_name, file_name, start_index);

    // Save the information about this plot in the same file as the data
    params.Add("row 1", "values of alpha");
    params.Add("row 2", "values of error");
    params.Save(path_to_file);

    // Save the data
    save_vector(alpha, path_to_file);
    save_vector(error, path_to_file);

    // Add the simulation name to the log file
    cout << "Data saved in " << path_to_file << endl;

    return 0;
}
