//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_NONE_HPP
#define STOSPA_NONE_HPP

#include <AbstractReaction.hpp>

class None : public AbstractReaction
{
    unsigned mSpeciesIndex;
public:
    explicit None(double reaction_rate, unsigned species)
    {
        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mSpeciesIndex = species;

        mReactionName = "None";
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
        return mRateConstant;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        return voxel_index;
    }
};


#endif //STOSPA_DECAY_HPP
