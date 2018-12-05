//
// Created by bartosz on 19/04/17.
//

#ifndef STOSPA_GRIDS_HPP
#define STOSPA_GRIDS_HPP

#include <vector>
#include <limits>
#include <cassert>

using namespace std;

class Grid
{
public:
    /** Infinity. */
    double inf = numeric_limits<double>::infinity();

    double voxelSize=1.0;

    vector<unsigned> numVoxels;

    /** Matrix of integers, representing the number of molecules within each voxel.
     * There is a vector of integers for species. */
    vector<vector<unsigned>> voxels;

    /** Vector of times until the next reaction for each voxel. */
    vector<double> time_increments;

    /** Vector of bounds for the total propensities. */
    vector<double> a_0;

    /** Constructor. */
    Grid(unsigned num_species, double voxel_size, unsigned num_voxels_x, unsigned int num_voxels_y=1);

    /** Default constructor. */
    Grid()=default;

};


#endif //STOSPA_GRIDS_HPP
