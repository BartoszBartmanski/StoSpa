//
// Created by bartosz on 23/06/17.
//

#ifndef STOSPA_SCHNAKENBERG_HPP
#define STOSPA_SCHNAKENBERG_HPP


#include "AbstractReaction.hpp"

class Schnakenberg : public AbstractReaction
{
public:
    explicit Schnakenberg(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_SCHNAKENBERG_HPP
