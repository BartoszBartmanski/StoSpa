//
// Created by bartosz on 23/06/17.
//

#include "Schnakenberg.hpp"

Schnakenberg::Schnakenberg(double reaction_rate)
{
    // Schnakenberg kinetics: 2A + B -> 3A
    // Species A - index 0
    // Species B - index 1

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "schnakenberg";
}

void Schnakenberg::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void Schnakenberg::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 2);
}

double Schnakenberg::GetPropensity(Grid& grid, const int& voxel_index)
{
    return mRateConstant * grid.voxels[0][voxel_index] * (grid.voxels[0][voxel_index] - 1) * grid.voxels[1][voxel_index] / pow(grid.voxelSize, 2);
}

int Schnakenberg::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] += 1;
    grid.voxels[1][voxel_index] -= 1;
    return voxel_index;
}