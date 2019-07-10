//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_DECAY_HPP
#define STOSPA_DECAY_HPP

#include <AbstractReaction.hpp>

class Decay : public AbstractReaction
{
    unsigned mSpeciesIndex;
public:
    explicit Decay(double reaction_rate, unsigned species)
    {
        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mSpeciesIndex = species;

        mReactionName = "Decay";
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
	unsigned num = grid.voxels[mSpeciesIndex][voxel_index];
	if (num > 0) { num -= 1; }
	return mRateConstant * num;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[mSpeciesIndex][voxel_index] -= 1;
        return voxel_index;
    }
};


#endif //STOSPA_DECAY_HPP
