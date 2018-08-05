//
// Created by bartosz on 05/07/17.
//

#include "GrayScott.hpp"

GrayScott_I::GrayScott_I(double reaction_rate)
{
    // Gray-scott reaction I kinetics: A + 2B -> 3B
    // Species A - index 0
    // Species B - index 1

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "GrayScott_I";
}

void GrayScott_I::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void GrayScott_I::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 3);
}

double GrayScott_I::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] * (grid.voxels[1][voxel_index] - 1) / pow(grid.voxelSize, 2);

    return propensity;
}

int GrayScott_I::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] -= 1;
    grid.voxels[1][voxel_index] += 1;
    return voxel_index;
}

GrayScott_II::GrayScott_II(double reaction_rate)
{
    // Gray-scott reaction II kinetics: B -> C
    // Species B - index 1
    // Species C - index 2

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "GrayScott_II";
}

void GrayScott_II::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void GrayScott_II::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 3);
}

double GrayScott_II::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[1][voxel_index];

    return propensity;
}

int GrayScott_II::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[1][voxel_index] -= 1;
    grid.voxels[2][voxel_index] += 1;
    return voxel_index;
}