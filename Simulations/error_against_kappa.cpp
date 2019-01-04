#include <iostream>
#include "docopt.h"
#include "Parameters.hpp"
#include "DiffEqAnalytic.hpp"
#include "Simulation2d.hpp"
#include "Simulator.hpp"
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
      -j --num_threads=<num_threads>                Number of threads to be used. [default: 1]
)";

vector<double> get_kappa(unsigned num_points, unsigned num_voxels)
{
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

    return kappa;
}

int main(int argc, const char** argv)
{
    // Get command line input
    map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);
    Parameters p(args);
    p.SetComments("Data for a plot of values of error against values of kappa in the 2d simulations with kappa being varied.");
    p.SetCommand(arr_to_str(argc, argv));
    p.SetNumDims(2);

    // Name the file
    string file_name = "error_against_kappa";
    if (args["--append"]) { file_name += ("_" + args["--append"].asString()); }

    // Get vector of values of kappa and correct the number of points variable
    vector<double> kappa = get_kappa(p.GetNumPoints(), p.GetNumVoxels());
    auto num_points = unsigned(kappa.size());
    p.SetNumPoints(num_points);

    // Initialise the containers
    vector<string> methods = {"fem", "fdm", "fvm", "fet"};
    map<string, vector<double>> error;
    for (const auto& method : methods)
    {
        error[method] = vector<double>(num_points, 0);
    }

    // Initialise progress object
    Progress prog(num_points);

    // Collect the data
    #pragma omp parallel for num_threads(p.GetNumThreads())
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

        for (const auto& method : methods)
        {
            Simulation2d sim(p.GetNumRuns(), p.GetNumSpecies(), p.GetNumVoxels(), p.GetDomainBounds(), p.GetBC(),
                              kappa[i]);
            unique_ptr<JumpRate> jump_rate = get_jump_rates(sim.GetVoxelDims(), method, p.GetAlpha(), p.GetBeta());
            sim.SetDiffusionRate(move(jump_rate), p.GetDiff()[0], 0);
            sim.AddReaction(make_unique<Decay>(p.GetDecay()[0], 0));
            sim.AddReaction(make_unique<Production>(p.GetProd()[0], 0));
            sim.SetInitialNumMolecules(floor_div(sim.GetNumVoxels(), 2), p.GetInitialNum()[0], 0);
            sim.Advance(p.GetEndTime());
            error[method][i] = sim.GetError(sol);
        }

        // Show progress of the simulation
        prog.Show();
    }

    // Check that the directory exists and that the no file is being over-written.
    string path_to_file = update_path(p.GetSaveDir(), file_name, p.GetStartIndex());

    // Save the information about this plot in the same file as the data
    p.Add("row 1", "values of kappa");
    p.Add("row 2", "values of error using fem");
    p.Add("row 3", "values of error using fdm");
    p.Add("row 4", "values of error using fvm");
    p.Add("row 5", "values of error using fet");
    p.Save(path_to_file);

    // Save the data
    save_vector(kappa, path_to_file);
    save_vector(error["fem"], path_to_file);
    save_vector(error["fdm"], path_to_file);
    save_vector(error["fvm"], path_to_file);
    save_vector(error["fet"], path_to_file);

    // Add the simulation name to the log file
    prog.End(path_to_file);

    return 0;
}

