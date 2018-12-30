#include <iostream>
#include <future>
#include "docopt.h"
#include "Parameters.hpp"
#include "DiffEqAnalytic.hpp"
#include "Simulation_2d.hpp"
#include "Decay.hpp"
#include "Production.hpp"

static const char USAGE[] =
        R"(Produces a data for the plot of error against value of alpha.

    Usage:
      error_against_kappa [options]

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
    Parameters p(args);
    p.SetComments("Data for a plot of values of error against values of kappa in the 2d simulations with kappa being varied.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumDims(2);
    p.SetNumThreads(4);

    // Name the file
    string file_name = "error_against_kappa";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Get the number of points of kappa for which to run the simulation
    unsigned num_points = p.GetNumPoints();

    // Get rid of the duplicate kappa values in the vector that holds the values
    vector<double> tmp_kappa = linspace(1.0, 1.4, num_points);
    vector<double> kappa(1, 1.0);

    for (unsigned i=1; i < tmp_kappa.size(); i++)
    {
        double next_kappa =  unsigned(tmp_kappa[i]*p.GetNumVoxels())/double(p.GetNumVoxels());
        if (next_kappa != kappa.back())
        {
            kappa.push_back(next_kappa);
        }
    }

    // Correct the num_points variable
    num_points = unsigned(kappa.size());
    p.SetNumPoints(num_points);

    vector<string> methods = {"fem", "fdm", "fvm", "fet"};
    vector<vector<double>> error(methods.size());
    vector<future<double>> futures(methods.size());
    vector<Simulation_2d> sims;

    // Initialise progress object
    Progress prog(num_points);

    // Collect the data
    for (unsigned i=0; i < num_points; i++)
    {

        vector<double> mid_point = get_midpoint(p.GetDomainBounds(), p.GetNumVoxels(), kappa[i]);

        DiffEqAnalytic analytic(2,
                                p.GetEndTime(),
                                mid_point,
                                p.GetDomainBounds(),
                                {p.GetNumVoxels(), unsigned(kappa[i]*p.GetNumVoxels())},
                                p.GetInitialNum()[0],
                                p.GetDiff()[0],
                                p.GetDecay()[0],
                                p.GetProd()[0],
                                p.GetTruncOrder());
        vector<double> sol = analytic.GetAnalytic();

        sims.clear();
        for (const string& method : methods)
        {
            sims.emplace_back(Simulation_2d(p.GetNumRuns(), p.GetNumSpecies(), method, p.GetNumVoxels(), p.GetDomainBounds(), p.GetBC(), kappa[i]));
            sims.back().SetAlpha(p.GetAlpha());
            sims.back().SetBeta(p.GetBeta());
            sims.back().SetDiffusionRate(p.GetDiff()[0], 0);
            sims.back().AddReaction(make_unique<Decay>(p.GetDecay()[0], 0));
            sims.back().AddReaction(make_unique<Production>(p.GetProd()[0], 0));
            sims.back().SetInitialNumMolecules(floor_div(sims.back().GetNumVoxels(), 2), p.GetInitialNum()[0], 0);
        }

        for (unsigned j=0; j< methods.size(); j++)
        {
            futures[j] = async(launch::async, get_error, ref(sims[j]), sol, p.GetEndTime());
        }

        for (unsigned j=0; j< methods.size(); j++)
        {
            error[j].push_back(futures[j].get());
        }

        // Show progress of the simulation
        prog.Show();
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(p.GetSaveDir(), file_name, p.GetStartIndex());

    // Save the information about this plot in the same file as the data
    p.Add("row 1", "values of kappa");
    for (unsigned j=0; j< methods.size(); j++)
    {
        p.Add("row " + to_string(j+2), "values of error using " + methods[j]);
    }
    p.Save(path_to_file);

    // Save the data
    save_vector(kappa, path_to_file);
    for (unsigned j=0; j< methods.size(); j++)
    {
        save_vector(error[j], path_to_file);
    }

    // Add the simulation name to the log file
    prog.End(path_to_file);

    return 0;
}

