//
// Created by bartosz on 25/04/18.
//

#ifndef STOSPA_DIFFUSIONPERIODIC_HPP
#define STOSPA_DIFFUSIONPERIODIC_HPP

#include "AbstractReaction.hpp"
#include <Utilities.hpp>

class DiffusionPeriodic : public AbstractReaction
{
    unsigned mSpeciesIndex;

    vector<int> mUnflattenedIndex;

    vector<int> mDirection;

public:
    DiffusionPeriodic(double reaction_rate, unsigned species, vector<int> direction);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_DIFFUSIONPERIODIC_HPP
