//
// Created by bartosz on 20/04/18.
//

#ifndef STOSPA_DECAY_HPP
#define STOSPA_DECAY_HPP

#include <AbstractReaction.hpp>

class Decay : public AbstractReaction
{
    unsigned mSpeciesIndex;
public:
    explicit Decay(double reaction_rate, unsigned species);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;
};


#endif //STOSPA_DECAY_HPP
