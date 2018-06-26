//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_PRODUCTION_HPP
#define STOSPA_PRODUCTION_HPP

#include <AbstractReaction.hpp>

class Production : public AbstractReaction
{
private:
    unsigned mSpeciesIndex;

public:
    Production(double reaction_rate, unsigned species);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;
};


#endif //STOSPA_PRODUCTION_HPP
