//
// Created by bartosz on 18/07/18.
//

#ifndef STOSPA_MICHAELISMENTINREDUCED_HPP
#define STOSPA_MICHAELISMENTINREDUCED_HPP

#include "AbstractReaction.hpp"

// S ->[k_3 E_t/(S+K_m)] P
// Species S - index 0
// Species P - index 1

class MichaelisMentinReduced : public AbstractReaction
{
protected:
    double mTotalEnzyme;

    double mKm;

public:
    explicit MichaelisMentinReduced(double reaction_rate, double e_t, double k_m)
    {
        // Michaelis-Mentin kinetics: S -> P
        assert(reaction_rate >= 0);
        assert(e_t >= 0);
        assert(k_m > 0);

        mRateConstant = reaction_rate;
        mTotalEnzyme = e_t;
        mKm = k_m;
        mReactionName = "MichaelisMentinReduced";
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

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        double propensity = mRateConstant * mTotalEnzyme * grid.voxels[0][voxel_index] / (grid.voxels[0][voxel_index] + grid.voxelSize * mKm);

        return propensity;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        unsigned future_0 = grid.voxels[0][voxel_index] - 1;

        double propensity = mRateConstant * mTotalEnzyme * future_0 / (future_0 + grid.voxelSize * mKm);

        return propensity;
    }


    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        grid.voxels[0][voxel_index] -= 1;
        grid.voxels[1][voxel_index] += 1;
        return voxel_index;
    }
};


#endif //STOSPA_MICHAELISMENTINREDUCED_HPP
