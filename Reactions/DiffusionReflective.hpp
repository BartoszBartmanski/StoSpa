//
// Created by bartosz on 21/04/18.
//

#ifndef STOSPA_DIFFUSIONREFLECTIVE_HPP
#define STOSPA_DIFFUSIONREFLECTIVE_HPP

#include <AbstractReaction.hpp>

class DiffusionReflective : public AbstractReaction
{
    unsigned mSpeciesIndex;

    vector<int> mUnflattenedIndex;

    vector<int> mDirection;

public:
    DiffusionReflective(double reaction_rate, unsigned species, vector<int> direction)
    {
        assert(reaction_rate >= 0);
        mReactionName = "DiffusionReflective";

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

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        int jump_index = voxel_index;
        mUnflattenedIndex[0] = (voxel_index % grid.numVoxels[0]) + mDirection[0];
        mUnflattenedIndex[1] = (voxel_index / grid.numVoxels[0]) + mDirection[1];

        unsigned check = 0;
        for (unsigned i=0; i < 2; i++)
        {
            check += (mUnflattenedIndex[i] < 0 or mUnflattenedIndex[i] >= int(grid.numVoxels[i]));
        }

        if (check == 0)
        {
            jump_index = mUnflattenedIndex[1] * grid.numVoxels[0] + mUnflattenedIndex[0];
            grid.voxels[mSpeciesIndex][voxel_index] -= 1;
            grid.voxels[mSpeciesIndex][jump_index] += 1;
        }
        return jump_index;
    }

};


#endif //STOSPA_DIFFUSIONREFLECTIVE_HPP
