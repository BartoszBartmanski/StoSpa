//
// Created by bartosz on 18/07/18.
//

#include "MichaelisMentin.hpp"

MichaelisMentin_I::MichaelisMentin_I(double reaction_rate)
{
    // Michaelis-Mentin kinetics: E + S -> C
    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "MichaelisMentin_I";
}

void MichaelisMentin_I::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void MichaelisMentin_I::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 4);
}

double MichaelisMentin_I::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] / grid.voxelSize;

    return propensity;
}

int MichaelisMentin_I::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] -= 1;
    grid.voxels[1][voxel_index] -= 1;
    grid.voxels[2][voxel_index] += 1;
    return voxel_index;
}

MichaelisMentin_II::MichaelisMentin_II(double reaction_rate)
{
    // Michaelis-Mentin kinetics: E + S <- C

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "MichaelisMentin_II";
}

void MichaelisMentin_II::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void MichaelisMentin_II::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 4);
}

double MichaelisMentin_II::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[2][voxel_index];

    return propensity;
}

int MichaelisMentin_II::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] += 1;
    grid.voxels[1][voxel_index] += 1;
    grid.voxels[2][voxel_index] -= 1;
    return voxel_index;
}

MichaelisMentin_III::MichaelisMentin_III(double reaction_rate)
{
    // Michaelis-Mentin kinetics: C -> E + P

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "MichaelisMentin_III";
}

void MichaelisMentin_III::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void MichaelisMentin_III::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 4);
}

double MichaelisMentin_III::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[2][voxel_index];

    return propensity;
}

int MichaelisMentin_III::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[2][voxel_index] -= 1;
    grid.voxels[0][voxel_index] += 1;
    grid.voxels[3][voxel_index] += 1;
    return voxel_index;
}