//
// Created by bartosz on 21/04/18.
//

#ifndef STOSPA_DIFFUSIONREFLECTIVE_HPP
#define STOSPA_DIFFUSIONREFLECTIVE_HPP

#include <AbstractReaction.hpp>

class DiffusionReflective : public AbstractReaction
{
    unsigned mSpeciesIndex;

    vector<int> mUnflattenedIndex;

    vector<int> mDirection;

public:
    DiffusionReflective(double reaction_rate, unsigned species, vector<int> direction);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_DIFFUSIONREFLECTIVE_HPP
