//
// Created by bartosz on 02/07/17.
//

#ifndef STOSPA_DIMERISATION_HPP
#define STOSPA_DIMERISATION_HPP

#include "AbstractReaction.hpp"

class Dimerisation : public AbstractReaction
{
protected:
    /** Species for which this reaction will take place. */
    unsigned mSpeciesIndex;

public:

    Dimerisation(double reaction_rate, unsigned int species_index)
    {
        // Dimerisation kinetics 2A -> 0

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mSpeciesIndex = species_index;
        mReactionName = "Dimerisation";
    }

    void SetRateConstant(double rate_constant) override
    {
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species > 0);
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[mSpeciesIndex][voxel_index] * (grid.voxels[mSpeciesIndex][voxel_index] - 1) / grid.voxelSize;

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        if (grid.voxels[mSpeciesIndex][voxel_index] >= 2)
        {
            grid.voxels[mSpeciesIndex][voxel_index] -= 2;
        }
        return voxel_index;
    }

};


#endif //STOSPA_DIMERISATION_HPP
