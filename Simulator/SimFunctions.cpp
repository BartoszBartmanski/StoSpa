#include "SimFunctions.hpp"

double get_error(AbstractSimulation& sim, const vector<double>& analytic, double end_time)
{
    sim.Advance(end_time);
    return sim.GetError(analytic);
}

double get_voxel_size(const vector<double>& domain_bounds, const unsigned& num_voxels, const double& kappa)
{
    double voxel_size;

    if (!kappa) // i.e. in one dimension
    {
        voxel_size = (domain_bounds[1] - domain_bounds[0]) / double(num_voxels);
    }
    else  // in two dimensions
    {
        double h = (domain_bounds[1] - domain_bounds[0]) / floor(kappa * num_voxels);
        voxel_size = pow(h, 2) * kappa;
    }

    return voxel_size;
}

vector<double> get_midpoint(const vector<double>& domain_bounds, const unsigned& num_voxels, const double& kappa)
{
    vector<unsigned> voxels = {num_voxels, unsigned(kappa*num_voxels)};
    vector<double> midpoint(2);

    for (unsigned i=0; i < 2; i++)
    {
        double size = (domain_bounds[1] - domain_bounds[0]) / voxels[i];
        midpoint[i] = (voxels[i]/2 + 0.5 ) * size;
    }
    return midpoint;
}
