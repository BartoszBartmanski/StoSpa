//
// Created by bartosz on 05/07/17.
//

#ifndef STOSPA_GRAYSCOTT_I_HPP
#define STOSPA_GRAYSCOTT_I_HPP

#include "AbstractReaction.hpp"

class GrayScott_I : public AbstractReaction
{
public:
    explicit GrayScott_I(double reaction_rate)
    {
        // Gray-scott reaction I kinetics: A + 2B -> 3B
        // Species A - index 0
        // Species B - index 1

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "GrayScott_I";
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 3);
        (void)num_species;
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] * (grid.voxels[1][voxel_index] - 1) / pow(grid.voxelSize, 2);

        return propensity;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        unsigned future_0 = grid.voxels[0][voxel_index] - 1;
        unsigned future_1 = grid.voxels[1][voxel_index] + 1;

        double propensity = mRateConstant * future_0 * future_1 * (future_1 - 1) / pow(grid.voxelSize, 2);

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] -= 1;
        grid.voxels[1][voxel_index] += 1;
        return voxel_index;
    }
};

class GrayScott_II : public AbstractReaction
{
public:
    explicit GrayScott_II(double reaction_rate)
    {
        // Gray-scott reaction II kinetics: B -> C
        // Species B - index 1
        // Species C - index 2

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "GrayScott_II";
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 3);
        (void)num_species;
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[1][voxel_index];

        return propensity;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        unsigned future = grid.voxels[1][voxel_index] - 1;
        double propensity = mRateConstant * future;

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[1][voxel_index] -= 1;
        grid.voxels[2][voxel_index] += 1;
        return voxel_index;
    }

};


#endif //STOSPA_GRAYSCOTT_I_HPP
