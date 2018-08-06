//
// Created by bartosz on 23/06/17.
//

#ifndef STOSPA_SCHNAKENBERG_HPP
#define STOSPA_SCHNAKENBERG_HPP


#include "AbstractReaction.hpp"

class Schnakenberg : public AbstractReaction
{
public:
    explicit Schnakenberg(double reaction_rate)
    {
        // Schnakenberg kinetics: 2A + B -> 3A
        // Species A - index 0
        // Species B - index 1

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "schnakenberg";
    }

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 2);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        return mRateConstant * grid.voxels[0][voxel_index] * (grid.voxels[0][voxel_index] - 1) * grid.voxels[1][voxel_index] / pow(grid.voxelSize, 2);
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] += 1;
        grid.voxels[1][voxel_index] -= 1;
        return voxel_index;
    }

};


#endif //STOSPA_SCHNAKENBERG_HPP
