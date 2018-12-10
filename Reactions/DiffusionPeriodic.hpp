//
// Created by bartosz on 25/04/18.
//

#ifndef STOSPA_DIFFUSIONPERIODIC_HPP
#define STOSPA_DIFFUSIONPERIODIC_HPP

#include "AbstractReaction.hpp"
#include <Utilities.hpp>

class DiffusionPeriodic : public AbstractReaction
{
    unsigned mSpeciesIndex;

    vector<int> mUnflattenedIndex;

    vector<int> mDirection;

public:
    DiffusionPeriodic(double reaction_rate, unsigned species, vector<int> direction)
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

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species > 0);
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
    }


    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        mUnflattenedIndex[0] = mod((voxel_index % grid.numVoxels[0]) + mDirection[0], grid.numVoxels[0]);
        mUnflattenedIndex[1] = mod((voxel_index / grid.numVoxels[0]) + mDirection[1], grid.numVoxels[1]);

        int jump_index = mUnflattenedIndex[1] * grid.numVoxels[0] + mUnflattenedIndex[0];
        grid.voxels[mSpeciesIndex][voxel_index] -= 1;
        grid.voxels[mSpeciesIndex][jump_index] += 1;

        return jump_index;
    }

};


#endif //STOSPA_DIFFUSIONPERIODIC_HPP
