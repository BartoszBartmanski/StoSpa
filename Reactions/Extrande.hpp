//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_NONE_HPP
#define STOSPA_NONE_HPP

#include <AbstractReaction.hpp>

class Extrande : public AbstractReaction
{
public:
    explicit Extrande()
    {
        // Rate constant here being actual total propensity in a voxel
        mRateConstant = 0;

        mReactionName = "Extrande";
    }

    double GetPropensity(const Grid& grid, const int& voxel_index) override
    {
        (void)grid;
        (void)voxel_index;
        return mRateConstant;
    }

    double GetFuturePropensity(const Grid& grid, const int& voxel_index) override
    {
        (void)grid;
        (void)voxel_index;
        return 0;
    }

    int UpdateGrid(Grid& grid, const int& voxel_index) override
    {
        (void)grid;
        return voxel_index;
    }
};


#endif //STOSPA_DECAY_HPP
