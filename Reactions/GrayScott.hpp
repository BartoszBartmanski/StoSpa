//
// Created by bartosz on 05/07/17.
//

#ifndef STOSPA_GRAYSCOTT_I_HPP
#define STOSPA_GRAYSCOTT_I_HPP

#include "AbstractReaction.hpp"

class GrayScott_I : public AbstractReaction
{
public:
    explicit GrayScott_I(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;
};

class GrayScott_II : public AbstractReaction
{
public:
    explicit GrayScott_II(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};


#endif //STOSPA_GRAYSCOTT_I_HPP
