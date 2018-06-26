//
// Created by bartosz on 26/04/18.
//

#include "Voxel.hpp"

Voxel::Voxel(unsigned num_species, double voxel_size, vector<double> position, bool on_boundary)
{
    mNextReactionTime = inf;
    mNumMolecules.resize(num_species);
    mVoxelSize = voxel_size;
    mPosition = move(position);
    mOnBoundary = on_boundary;
}

void Voxel::UpdateNextReactionTime(double time)
{
    assert(time >= 0);
    mNextReactionTime = time;
}

double Voxel::GetNextReactionTime()
{
    return mNextReactionTime;
}

double Voxel::GetVoxelSize()
{
    return mVoxelSize;
}

vector<double> Voxel::GetPosition()
{
    return mPosition;
}

bool Voxel::OnBoundary()
{
    return mOnBoundary;
}