//
// Created by bartosz on 19/04/17.
//

#ifndef STOSPA_GRIDS_HPP
#define STOSPA_GRIDS_HPP

#include <vector>
#include <limits>
#include <cassert>
#include <boost/heap/fibonacci_heap.hpp>

using namespace std;

class Grid
{
public:

    double scale=1.0;

    double voxelSize=1.0;

    double time=0.0;

    vector<unsigned> numVoxels;

    /** Matrix of integers, representing the number of molecules within each voxel.
     * There is a vector of integers for species. */
    vector<vector<unsigned>> voxels;

    /** Vector of times until the next reaction for each voxel. */
    boost::heap::fibonacci_heap<pair<double, unsigned>> next_reaction_time;

    vector<boost::heap::fibonacci_heap<pair<double, unsigned>>::handle_type> handles;

    /** Vector of bounds for the total propensities. */
    vector<double> a_0;

    /** Constructor. */
    Grid(unsigned num_species, double voxel_size, unsigned num_voxels_x, unsigned int num_voxels_y=1);

    /** Default constructor. */
    Grid()=default;

};


#endif //STOSPA_GRIDS_HPP
