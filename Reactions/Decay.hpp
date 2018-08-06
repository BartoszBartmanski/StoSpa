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
        return mRateConstant * grid.voxels[mSpeciesIndex][voxel_index];
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[mSpeciesIndex][voxel_index] -= 1;
        return voxel_index;
    }
};


#endif //STOSPA_DECAY_HPP
