//
// Created by bartosz on 24/02/18.
//

#ifndef STOSPA_DIFFEQANALYTIC_HPP
#define STOSPA_DIFFEQANALYTIC_HPP

#include <vector>
#include "VectorFunctions.hpp"

using namespace std;

class DiffEqAnalytic
{
protected:
    vector<double> mSolution;

public:
    /**
     * Default constructor
     */
    DiffEqAnalytic(unsigned num_dims,
                   double time,
                   vector<double> x_0,
                   vector<double> domain,
                   vector<unsigned> num_points,
                   unsigned total_molecules,
                   double diff,
                   double decay,
                   double prod,
                   unsigned trunc_order=1000,
                   bool voxels=true);

    ~DiffEqAnalytic()= default;

    vector<double> GetAnalytic();
};


#endif //STOSPA_DIFFEQANALYTIC_HPP
