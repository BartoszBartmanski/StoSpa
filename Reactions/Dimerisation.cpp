//
// Created by bartosz on 02/07/17.
//

#include "Dimerisation.hpp"

Dimerisation::Dimerisation(double reaction_rate, unsigned int species)
{
    // Dimerisation kinetics 2A -> 0

    assert(reaction_rate >= 0);

    mRateConstant = reaction_rate;
    mSpeciesIndex = species;
    mReactionName = "Dimerisation";
}

void Dimerisation::SetRateConstant(double rate_constant)
{
    mRateConstant = rate_constant;
}

void Dimerisation::CheckNumSpecies(unsigned num_species)
{
    assert(num_species > 0);
}

double Dimerisation::GetPropensity(Grid& grid, const int& voxel_index)
{
    double propensity = mRateConstant * grid.voxels[mSpeciesIndex][voxel_index] * (grid.voxels[mSpeciesIndex][voxel_index] - 1) / grid.voxelSize;

    return propensity;
}

int Dimerisation::UpdateGrid(Grid& grid, const int& voxel_index)
{
    if (grid.voxels[mSpeciesIndex][voxel_index] >= 2)
    {
        grid.voxels[mSpeciesIndex][voxel_index] -= 2;
    }
    return voxel_index;
}
