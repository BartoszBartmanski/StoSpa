//
// Created by bartosz on 19/04/17.
//

#include "Grid.hpp"

Grid::Grid(unsigned num_species, double voxel_size, unsigned num_voxels_x, unsigned int num_voxels_y)
{
    // Check for sensible input
    assert(num_species > 0);
    assert(voxel_size > 0);
    assert(num_voxels_x > 0);
    assert(num_voxels_y > 0);

    numVoxels = {num_voxels_x, num_voxels_y};
    voxelSize = voxel_size;
    voxels.resize(num_species);
    time = 0.0;

    for (unsigned i=0; i<num_species; i++)
    {
        voxels[i] = vector<unsigned>(num_voxels_x*num_voxels_y, 0);
    }
    a_0 = vector<double>(num_voxels_x*num_voxels_y, 0);

    handles.resize(num_voxels_x * num_voxels_y);
}
