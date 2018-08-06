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

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 3);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] * (grid.voxels[1][voxel_index] - 1) / pow(grid.voxelSize, 2);

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

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 3);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[1][voxel_index];

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
