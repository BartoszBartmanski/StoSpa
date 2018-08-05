//
// Created by bartosz on 18/07/18.
//

#include "MichaelisMentinReduced.hpp"

MichaelisMentinReduced::MichaelisMentinReduced(double reaction_rate, double e_t, double k_m)
{
    // Michaelis-Mentin kinetics: S -> P
    assert(reaction_rate >= 0);
    assert(e_t >= 0);
    assert(k_m > 0);

    mRateConstant = reaction_rate;
    mTotalEnzyme = e_t;
    mKm = k_m;
    mReactionName = "MichaelisMentinReduced";
}

void MichaelisMentinReduced::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void MichaelisMentinReduced::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 2);
}

double MichaelisMentinReduced::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * mTotalEnzyme * grid.voxels[0][voxel_index] / (grid.voxels[0][voxel_index] + grid.voxelSize * mKm);

    return propensity;
}

int MichaelisMentinReduced::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] -= 1;
    grid.voxels[1][voxel_index] += 1;
    return voxel_index;
}