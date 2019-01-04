//
// Created by bartmanski on 22/11/17.
//

#ifndef STOSPA_TWOSPECIESDECAY_HPP
#define STOSPA_TWOSPECIESDECAY_HPP

#include "AbstractReaction.hpp"

class TwoSpeciesDecay : public AbstractReaction
{
public:
    explicit TwoSpeciesDecay(double reaction_rate)
    {
        // Schnakenberg kinetics: A + B -> B
        // Species A - index 0
        // Species B - index 1

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "TwoSpeciesDecay";
    }

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 2);
        (void)num_species;
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] / grid.voxelSize;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        unsigned future_0 = grid.voxels[0][voxel_index] - 1;
        unsigned future_1 = grid.voxels[1][voxel_index];

        return mRateConstant * future_0 * future_1 / grid.voxelSize;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] -= 1;
        return voxel_index;
    }

};


#endif //STOSPA_TWOSPECIESDECAY_HPP
