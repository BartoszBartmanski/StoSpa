//
// Created by bartmanski on 22/11/17.
//

#ifndef STOSPA_TWOSPECIESDECAY_HPP
#define STOSPA_TWOSPECIESDECAY_HPP

#include "AbstractReaction.hpp"

class TwoSpeciesDecay : public AbstractReaction
{
public:
    explicit TwoSpeciesDecay(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_TWOSPECIESDECAY_HPP
