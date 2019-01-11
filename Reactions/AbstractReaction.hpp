//
// Created by bartosz on 23/06/17.
//

#ifndef STOSPA_ABSTRACTREACTION_HPP
#define STOSPA_ABSTRACTREACTION_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "Grid.hpp"

class AbstractReaction
{
protected:

    /** The name of the reaction. */
    string mReactionName="";

    /** The input rate constant. */
    double mRateConstant=0.0;

public:

    AbstractReaction()= default;

    /** Default destructor. */
    virtual ~AbstractReaction()= default;

    /** Set the rate constant for this reaction. */
    void SetRateConstant(double rate_constant)
    {
        assert(rate_constant >= 0);
        mRateConstant = rate_constant;
    }

    /** Get the rate constant. */
    double GetRateConstant()
    {
        return mRateConstant;
    }

    /** Set the reaction name. */
    void SetReactionName(string reaction_name)
    {
        mReactionName = move(reaction_name);
    }

    /** Get the reaction name. */
    string GetReactionName()
    {
        return mReactionName;
    }

    /** A check that the simulation has the appropriate number of species. */
    virtual void CheckNumSpecies(unsigned num_species)
    {
        assert(num_species > 0);
        (void)num_species;
    }

    /** Get the propensity (this will be used within a simulation class). */
    virtual double GetPropensity(const Grid& grid, const int& voxel_index)=0;

    virtual double GetFuturePropensity(const Grid& grid, const int& voxel_index)=0;

    /** Update the grids that hold the molecules (this will be used within a simulation class). */
    virtual int UpdateGrid(Grid& grid, const int& voxel_index)=0;

};


#endif //STOSPA_ABSTRACTREACTION_HPP
