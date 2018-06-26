//
// Created by bartosz on 20/04/18.
//

#include "Decay.hpp"

Decay::Decay(double reaction_rate, unsigned species)
{
    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mSpeciesIndex = species;

    mReactionName = "Decay";
}

void Decay::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void Decay::CheckNumSpecies(unsigned num_species)
{
    assert(num_species > 0);
}

double Decay::GetPropensity(Grid& grid, const int& voxel_index)
{
    return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
}

int Decay::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[mSpeciesIndex][voxel_index] -= 1;
    return voxel_index;
}