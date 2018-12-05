//
// Created by bartosz on 18/07/18.
//

#ifndef STOSPA_MICHAELISMENTIN_HPP
#define STOSPA_MICHAELISMENTIN_HPP

#include "AbstractReaction.hpp"

// E + S [k_2]<->[k_1] C ->[k_3] E + P
// Species E - index 0
// Species S - index 1
// Species C - index 2
// Species P - index 3

class MichaelisMentin_I : public AbstractReaction
{
public:
    explicit MichaelisMentin_I(double reaction_rate)
    {
        // Michaelis-Mentin kinetics: E + S -> C
        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "MichaelisMentin_I";
    }

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 4);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[0][voxel_index] * grid.voxels[1][voxel_index] / grid.voxelSize;

        return propensity;
    }

    double GetFuturePropensity(Grid& grid, const int& voxel_index) override
    {
        unsigned future = grid.voxels[2][voxel_index] + 1;
        double propensity = mRateConstant * future;

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] -= 1;
        grid.voxels[1][voxel_index] -= 1;
        grid.voxels[2][voxel_index] += 1;
        return voxel_index;
    }
};

class MichaelisMentin_II : public AbstractReaction
{
public:
    explicit MichaelisMentin_II(double reaction_rate)
    {
        // Michaelis-Mentin kinetics: E + S <- C

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "MichaelisMentin_II";
    }

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 4);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[2][voxel_index];

        return propensity;
    }

    double GetFuturePropensity(Grid& grid, const int& voxel_index) override
    {
        unsigned future = grid.voxels[2][voxel_index] - 1;
        double propensity = mRateConstant * future;

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] += 1;
        grid.voxels[1][voxel_index] += 1;
        grid.voxels[2][voxel_index] -= 1;
        return voxel_index;
    }

};

class MichaelisMentin_III : public AbstractReaction
{
public:
    explicit MichaelisMentin_III(double reaction_rate)
    {
        // Michaelis-Mentin kinetics: C -> E + P

        assert(reaction_rate >= 0);

        mRateConstant = reaction_rate;
        mReactionName = "MichaelisMentin_III";
    }

    void SetRateConstant(double rate_constant) override
    {
        assert(rate_constant > 0);
        mRateConstant = rate_constant;
    }

    void CheckNumSpecies(unsigned num_species) override
    {
        assert(num_species == 4);
    }

    double GetPropensity(Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * grid.voxels[2][voxel_index];

        return propensity;
    }

    double GetFuturePropensity(Grid& grid, const int& voxel_index) override
    {
        unsigned future = grid.voxels[2][voxel_index] - 1;
        double propensity = mRateConstant * future;

        return propensity;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[2][voxel_index] -= 1;
        grid.voxels[0][voxel_index] += 1;
        grid.voxels[3][voxel_index] += 1;
        return voxel_index;
    }

};

#endif //STOSPA_MICHAELISMENTIN_HPP
