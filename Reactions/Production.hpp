//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_PRODUCTION_HPP
#define STOSPA_PRODUCTION_HPP

#include <AbstractReaction.hpp>

class Production : public AbstractReaction
{
private:
    unsigned mSpeciesIndex;

public:
    Production(double reaction_rate, unsigned species)
    {
        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mSpeciesIndex = species;
        mReactionName = "Production";
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        (void)voxel_index;
        return mRateConstant * grid.voxelSize;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        (void)voxel_index;
        return mRateConstant * grid.voxelSize;
    }


    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[mSpeciesIndex][voxel_index] += 1;
        return voxel_index;
    }
};


#endif //STOSPA_PRODUCTION_HPP
