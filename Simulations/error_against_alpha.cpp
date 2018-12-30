#include <iostream>
#include <future>
#include "docopt.h"
#include "Utilities.hpp"
#include "Parameters.hpp"
#include "DiffEqAnalytic.hpp"
#include "Simulation_2d.hpp"
#include "Decay.hpp"
#include "Production.hpp"

static const char USAGE[] =
        R"(Produces a data for the plot of error against value of alpha.

    Usage:
      error_against_alpha [options]

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
    Parameters p(args);
    p.SetComments("Data for a plot of values of error against values of alpha in the 2d simulations with alpha being varied.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumDims(2);
    p.SetNumMethod("fdm");

    // Name the file
    string file_name = "error_against_alpha";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Create a vector of alpha values
    vector<double> alpha = linspace(0.0, 1.0/p.GetKappa(), p.GetNumPoints());

    // Create a vector to hold the results
    vector<long double> error(p.GetNumPoints(), 0);

    // Calculate the initial point of the molecules
    vector<double> mid_point = get_midpoint(p.GetDomainBounds(), p.GetNumVoxels(), p.GetKappa());

    // Calculate the analytical solution
    DiffEqAnalytic analytic(2,
                            p.GetEndTime(),
                            mid_point,
                            p.GetDomainBounds(),
                            {p.GetNumVoxels(), unsigned(p.GetKappa()*p.GetNumVoxels())},
                            p.GetInitialNum()[0],
                            p.GetDiff()[0],
                            p.GetDecay()[0],
                            p.GetProd()[0],
                            p.GetTruncOrder());
    vector<double> sol = analytic.GetAnalytic();

    unsigned num_points = p.GetNumPoints();
    unsigned num_threads = p.GetNumThreads();
    unsigned num_sets = num_points / num_threads;

    // Create vector of futures
    vector<future<double>> futures(num_threads);
    vector<Simulation_2d> sims;

    // Initialise progress object
    Progress prog(p.GetNumPoints());

    // Collect the data
    for (unsigned set=0; set < num_sets; set++)
    {
        sims.clear();
        for (unsigned i=0; i < num_threads; i++)
        {
            unsigned k = set * num_threads + i;
            sims.emplace_back(Simulation_2d(p));
            sims.back().SetAlpha(alpha[k]);
            sims.back().SetDiffusionRate(p.GetDiff()[0], 0);
            sims.back().AddReaction(make_unique<Decay>(p.GetDecay()[0], 0));
            sims.back().AddReaction(make_unique<Production>(p.GetProd()[0], 0));
            sims.back().SetInitialNumMolecules(floor_div(sims[i].GetNumVoxels(), 2), p.GetInitialNum()[0], 0);
        }

        for (unsigned i=0; i < num_threads; i++)
        {
            futures[i] = async(launch::async, get_error, ref(sims[i]), sol, p.GetEndTime());
        }

        for (unsigned i=0; i < num_threads; i++)
        {
            unsigned k = set * num_threads + i;
            error[k] = futures[i].get();

            // Show progress of the simulation
            prog.Show();
        }
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(p.GetSaveDir(), file_name, p.GetStartIndex());

    // Save the information about this plot in the same file as the data
    p.Add("row 1", "values of alpha");
    p.Add("row 2", "values of error");
    p.Save(path_to_file);

    // Save the data
    save_vector(alpha, path_to_file);
    save_vector(error, path_to_file);

    // Add the simulation name to the log file
    prog.End(path_to_file);

    return 0;
}
