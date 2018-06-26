//
// Created by bartosz on 20/04/18.
//

#include "Production.hpp"

Production::Production(double reaction_rate, unsigned species)
{
    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mSpeciesIndex = species;
    mReactionName = "Production";
}

void Production::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void Production::CheckNumSpecies(unsigned num_species)
{
    assert(num_species > 0);
}

double Production::GetPropensity(Grid& grid, const int& voxel_index)
{
    return mRateConstant * grid.voxelSize;
}

int Production::UpdateGrid(Grid& grid, const int& voxel_index)
{
    grid.voxels[mSpeciesIndex][voxel_index] += 1;
    return voxel_index;
}