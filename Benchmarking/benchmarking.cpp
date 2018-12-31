
#include <chrono>
#include "Version.hpp"
#include "Simulation_1d.hpp"
#include "Simulation_2d.hpp"

void benchmark(unsigned num_dims, unsigned num_runs, const string& file_name)
{
    // Create a pointer to an AbstractSimulation
    unique_ptr<AbstractSimulation> sim;

    unsigned num_species = 1;
    string method = "fdm";
    unsigned num_voxels = 21;
    vector<double> bounds = {0.0, 20.0};
    string bc = "reflective";
    double kappa = 1.0;

    // Point the pointer to the object
    if (num_dims == 1)
    {
        sim = make_unique<Simulation_1d>(num_runs, num_species, method, num_voxels, bounds, bc);
    }
    else
    {
        sim = make_unique<Simulation_2d>(num_runs, num_species, method, num_voxels, bounds, bc, kappa);
    }

    sim->SetDiffusionRate(1.0, 0);
    sim->SetInitialNumMolecules(floor_div(sim->GetNumVoxels(), 2), 1000, 0);

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    sim->Advance(5.0);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

    // Calculate the elapsed time
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    time_t end_time = chrono::system_clock::to_time_t(t2);

    // Print the time it took to finish the simulation
    string runs = " run";
    if (num_runs > 1) { runs += "s"; }
    cout << "It took " << time_span.count() << " to finish running " << num_runs << runs << " of " << num_dims << "d simulation" << endl;

    // Save the results
    fstream outfile;
    outfile.open(string(PATH_TO_SOURCE) + "/Benchmarking/" + file_name, ios_base::app);
    outfile << time_span.count() << " on " << ctime(&end_time);
    outfile.close();
}

int main()
{
    benchmark(1, 1, "benchmark_1d.dat");
    benchmark(1, 50, "benchmark_many_runs_1d.dat");
    benchmark(2, 1, "benchmark_2d.dat");
    benchmark(2, 50, "benchmark_many_runs_2d.dat");

    return 0;
}
