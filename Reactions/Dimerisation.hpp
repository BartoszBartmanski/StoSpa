//
// Created by bartosz on 02/07/17.
//

#ifndef STOSPA_DIMERISATION_HPP
#define STOSPA_DIMERISATION_HPP

#include "AbstractReaction.hpp"

class Dimerisation : public AbstractReaction
{
protected:
    /** Species for which this reaction will take place. */
    unsigned mSpeciesIndex;

public:

    Dimerisation(double reaction_rate, unsigned int species_index);

    void SetRateConstant(double rate_constant) override ;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_DIMERISATION_HPP
