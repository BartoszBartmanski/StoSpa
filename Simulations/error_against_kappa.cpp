#include <iostream>
#include <future>
#include "docopt.h"
#include "SimFunctions.hpp"
#include "Parameters.hpp"
#include "DiffEqAnalytic.hpp"
#include "Simulation_2d.hpp"
#include "Decay.hpp"
#include "Production.hpp"

static const char USAGE[] =
        R"(Produces a data for the plot of error against value of alpha.

    Usage:
      error_against_kappa [options]
      error_against_kappa -h | --help

    Options:
      -h --help                                     Show this screen.
      --info                                        Show all the parameters.
      --end_time=<end_time>                         Time until to run the simulation [default: 5.0].
      --num_runs=<num_runs>                         Number of realisations [default: 1].
      --num_voxels=<num_voxels>                     Number of voxels. [default: 21]
      --domain_bounds=<domain_bounds>               The bounds of the domain. [default: 0.0,20.0]
      --bc=<bc>                                     The boundary condition [default: reflective].
      --alpha=<alpha>                               Value of alpha [default: 0.0].
      --beta=<beta>                                 Value of beta [default: 0.0,0.0].
      --dir_name=<dir_name>                         The name of save dir.
      --append=<append>                             Append filename.
      --start_index=<index>                         The start index of the simulations. [default: 0]
      --diff=<diff>                                 DiffusionReflective coefficient [default: 1.0].
      --decay=<decay>                               Rate of decay [default: 0.0].
      --prod=<prod>                                 Rate of production [default: 0.0].
      --initial_num=<val>                           Initial number of molecules [default: 5000000].
      --num_points=<points>                         Number of points of kappa for which to run the simulation. [default: 10]
      --trunc_order=<trunc>                         The truncation order of analytic solution. [default: 1000]
)";

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

    // Create a parameters object - used to order the parameters order in the comments of the data files
    Parameters params;
    params.SetComments("Data for a plot of values of error against values of kappa in the 2d simulations with kappa being varied.");
    params.SetCommand(arr_to_str(argc, argv));
    params.SetNumDims(2);
    params.SetNumThreads(4);

    // Get the save directory
    string dir_name = SAVE_DIR;
    if (args["--dir_name"]) { dir_name = args["--dir_name"].asString(); }

    // Name the file
    string file_name = "error_against_kappa";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Get the start index of the file extension
    auto start_index = unsigned(stoi(args["--start_index"].asString()));

    // Get the number of runs
    auto num_runs = unsigned(stoi(args["--num_runs"].asString()));
    params.SetNumRuns(num_runs);

    unsigned num_species = 1;
    params.SetNumSpecies(num_species);

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

    // Get the value of alpha
    double alpha = stod(args["--alpha"].asString());
    params.SetAlpha(alpha);

    // Get the value of beta
    vector<double> beta = split_stod(args["--beta"].asString());
    if (beta.size() == 1) { beta.push_back(beta[0]); }
    params.SetBeta(beta);

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

    // Get the truncation order for the analytic solution
    unsigned trunc_order = unsigned(stoi(args["--trunc_order"].asString()));
    params.SetTruncOrder(trunc_order);

    // Get the number of points of kappa for which to run the simulation
    unsigned num_points = unsigned(stoi(args["--num_points"].asString()));

    // Get rid of the duplicate kappa values in the vector that holds the values
    vector<double> tmp_kappa = linspace(1.0, 1.4, num_points);
    vector<double> kappa(1, 1.0);

    for (unsigned i=1; i < tmp_kappa.size(); i++)
    {
        double next_kappa =  unsigned(tmp_kappa[i]*num_voxels)/double(num_voxels);
        if (next_kappa != kappa.back())
        {
            kappa.push_back(next_kappa);
        }
    }

    // Correct the num_points variable
    num_points = unsigned(kappa.size());
    params.SetNumPoints(num_points);

    vector<string> methods = {"fem", "fdm", "fvm", "fet"};
    map<string, vector<double>> error;
    map<string, future<double>> futures;
    map<string, Simulation_2d> sims;

    // Initialise progress object
    Progress p(num_points);

    // Collect the data
    for (unsigned i=0; i < num_points; i++)
    {

        vector<double> mid_point = get_midpoint(domain_bounds, num_voxels, kappa[i]);

        DiffEqAnalytic analytic = DiffEqAnalytic(2,
                                                 end_time,
                                                 mid_point,
                                                 domain_bounds,
                                                 {num_voxels, unsigned(kappa[i]*num_voxels)},
                                                 initial_num,
                                                 diff,
                                                 decay,
                                                 prod,
                                                 trunc_order);
        vector<double> sol = analytic.GetAnalytic();

        for (const string& method : methods)
        {
            sims[method] = Simulation_2d(num_runs, num_species, method, num_voxels, domain_bounds, bc, kappa[i], alpha, beta[0], beta[1]);
            sims[method].SetDiffusionRate(diff, 0);
            sims[method].AddReaction(make_shared<Decay>(decay, 0));
            sims[method].AddReaction(make_shared<Production>(prod, 0));
            sims[method].SetInitialNumMolecules(floor_div(sims[method].GetNumVoxels(), 2), initial_num, 0);
        }

        for (const string& method : methods)
        {
            futures[method] = async(launch::async, get_error, ref(sims[method]), sol, end_time);
        }

        for (const string& method : methods)
        {
            error[method].push_back(futures[method].get());
        }

        // Show progress of the simulation
        p.Show();
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(dir_name, file_name, start_index);

    // Save the information about this plot in the same file as the data
    params.Add("row 1", "values of kappa");
    params.Add("row 2", "values of error using fem");
    params.Add("row 3", "values of error using fdm");
    params.Add("row 4", "values of error using fvm");
    params.Add("row 5", "values of error using fet");
    params.Save(path_to_file);

    // Save the data
    save_vector(kappa, path_to_file);
    for (const string& method : methods)
    {
        save_vector(error[method], path_to_file);
    }

    // Add the simulation name to the log file
    cout << "Data saved in " << path_to_file << endl;

    return 0;
}

