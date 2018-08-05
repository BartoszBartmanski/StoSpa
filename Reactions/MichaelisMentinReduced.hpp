//
// Created by bartosz on 18/07/18.
//

#ifndef STOSPA_MICHAELISMENTINREDUCED_HPP
#define STOSPA_MICHAELISMENTINREDUCED_HPP

#include "AbstractReaction.hpp"

// S ->[k_3 E_t/(S+K_m)] P
// Species S - index 0
// Species P - index 1

class MichaelisMentinReduced : public AbstractReaction
{
protected:
    double mTotalEnzyme;

    double mKm;

public:
    explicit MichaelisMentinReduced(double reaction_rate, double e_t, double k_m);

    void SetRateConstant(double rate_constant) override;

    void CheckNumSpecies(unsigned num_species) override;

    double GetPropensity(Grid& grid, const int& voxel_index) override;

    int UpdateGrid(Grid& grid, const int& voxel_index) override;
};


#endif //STOSPA_MICHAELISMENTINREDUCED_HPP
