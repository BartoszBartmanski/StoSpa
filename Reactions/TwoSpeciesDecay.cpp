//
// Created by bartmanski on 22/11/17.
//

#include "TwoSpeciesDecay.hpp"

TwoSpeciesDecay::TwoSpeciesDecay(double reaction_rate)
{
    // Schnakenberg kinetics: A + B -> B
    // Species A - index 0
    // Species B - index 1

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mReactionName = "TwoSpeciesDecay";
}

void TwoSpeciesDecay::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void TwoSpeciesDecay::CheckNumSpecies(unsigned num_species)
{
    assert(num_species == 2);
}

double TwoSpeciesDecay::GetPropensity(Grid& grid, const int& voxel_index)
{
    return mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] / grid.voxelSize;
}

int TwoSpeciesDecay::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[0][voxel_index] -= 1;
    return voxel_index;
}