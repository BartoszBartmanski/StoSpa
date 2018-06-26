//
// Created by bartosz on 25/04/18.
//

#include "DiffusionPeriodic.hpp"

DiffusionPeriodic::DiffusionPeriodic(double reaction_rate, unsigned species, vector<int> direction)
{
    assert(reaction_rate >= 0);
    mReactionName = "DiffusionPeriodic";

    mRateConstant = reaction_rate;
    mSpeciesIndex = species;

    if (direction.size() == 1)
    {
        direction.push_back(0);
    }
    mDirection = direction;
    mUnflattenedIndex.reserve(2);
}

void DiffusionPeriodic::SetRateConstant(double rate_constant)
{
    assert(rate_constant > 0);
    mRateConstant = rate_constant;
}

void DiffusionPeriodic::CheckNumSpecies(unsigned num_species)
{
    assert(num_species > 0);
}

double DiffusionPeriodic::GetPropensity(Grid& grid, const int& voxel_index)
{
    return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
}

int DiffusionPeriodic::UpdateGrid(Grid& grid, const int& voxel_index)
{
    mUnflattenedIndex[0] = mod((voxel_index % grid.numVoxels[0]) + mDirection[0], grid.numVoxels[0]);
    mUnflattenedIndex[1] = mod((voxel_index / grid.numVoxels[0]) + mDirection[1], grid.numVoxels[1]);

    int jump_index = mUnflattenedIndex[1] * grid.numVoxels[0] + mUnflattenedIndex[0];
    grid.voxels[mSpeciesIndex][voxel_index] -= 1;
    grid.voxels[mSpeciesIndex][jump_index] += 1;

    return jump_index;
}