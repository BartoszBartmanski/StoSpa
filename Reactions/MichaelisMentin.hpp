//
// Created by bartosz on 18/07/18.
//

#ifndef STOSPA_MICHAELISMENTIN_HPP
#define STOSPA_MICHAELISMENTIN_HPP

#include "AbstractReaction.hpp"

// E + S [k_2]<->[k_1] C ->[k_3] E + P
// Species E - index 0
// Species S - index 1
// Species C - index 2
// Species P - index 3

class MichaelisMentin_I : public AbstractReaction
{
public:
    explicit MichaelisMentin_I(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;
};

class MichaelisMentin_II : public AbstractReaction
{
public:
    explicit MichaelisMentin_II(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};

class MichaelisMentin_III : public AbstractReaction
{
public:
    explicit MichaelisMentin_III(double reaction_rate);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;

};

#endif //STOSPA_MICHAELISMENTIN_HPP
