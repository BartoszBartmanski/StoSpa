//
// Created by bartosz on 26/04/18.
//

#ifndef STOSPA_VOXEL_HPP
#define STOSPA_VOXEL_HPP

#include <vector>
#include <limits>
#include <cassert>

using namespace std;

class Voxel
{
protected:
    double inf = numeric_limits<double>::infinity();

    vector<unsigned> mNumMolecules;

    double mNextReactionTime;

    double mVoxelSize;

    vector<double> mPosition;

    bool mOnBoundary;

public:
    Voxel(unsigned num_species, double voxel_size, vector<double> position, bool on_boundary);

    void UpdateNextReactionTime(double time);

    double GetNextReactionTime();

    double GetVoxelSize();

    vector<double> GetPosition();

    bool OnBoundary();

};


#endif //STOSPA_VOXEL_HPP
