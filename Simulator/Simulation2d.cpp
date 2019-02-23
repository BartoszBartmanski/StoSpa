#include "Simulation2d.hpp"

Simulation2d::Simulation2d(unsigned num_runs, unsigned num_species, unsigned num_voxels, vector<double> domain_bounds,
                             string boundary_condition, double ratio)
{
    // First check the input parameters
    assert(num_runs > 0);
    assert(num_species > 0);
    assert(num_voxels > 0);
    assert(domain_bounds.size() == 2);
    assert(ratio > 0.0);

    // Input parameters
    mDim = 2;
    mNumRuns = num_runs;
    mNumSpecies = num_species;
    mNumVoxels = {num_voxels, unsigned(ratio * num_voxels)};
    mTotalNumVoxels = mNumVoxels[0] * mNumVoxels[1];
    mDomainBounds = move(domain_bounds);
    mBC = move(boundary_condition);

    // Defines the directions in which molecules can jump
    mJumpDirections = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};

    // Correction to mRatio
    mRatio = mNumVoxels[1]/double(mNumVoxels[0]);  // correction for some values of aspect ratio not being possible

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    double domain_size = mDomainBounds[1] - mDomainBounds[0];
    mVoxelDims = {domain_size/mNumVoxels[0], domain_size/mNumVoxels[1]};
    mInitialVoxelSize = mVoxelDims[0] * mVoxelDims[1];

    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mInitialVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }

}

Simulation2d::Simulation2d(Parameters params)
{
    // First check the input parameters
    assert(params.GetNumRuns() > 0);
    assert(params.GetNumSpecies() > 0);
    assert(params.GetNumVoxels() > 0);
    assert(params.GetDomainBounds().size() == 2);
    assert(params.GetKappa() > 0.0);
    assert(params.GetAlpha() >= 0.0);
    assert(params.GetBeta().size() == 2);
    assert(params.GetBeta()[0] >= 0.0 and params.GetBeta()[0] <= 1.0);
    assert(params.GetBeta()[1] >= 0.0 and params.GetBeta()[1] <= 1.0);

    // Input parameters
    mDim = 2;
    mNumRuns = params.GetNumRuns();
    mNumSpecies = params.GetNumSpecies();
    mNumVoxels = {params.GetNumVoxels(), unsigned(params.GetKappa() * params.GetNumVoxels())};
    mTotalNumVoxels = mNumVoxels[0] * mNumVoxels[1];
    mDomainBounds = params.GetDomainBounds();
    mBC = params.GetBC();

    // Defines the directions in which molecules can jump
    mJumpDirections = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};

    // Additional input parameters (due to working in 2d)
    mRatio = mNumVoxels[1]/double(mNumVoxels[0]);  // correction for some values of aspect ratio not being possible

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    double domain_size = mDomainBounds[1] - mDomainBounds[0];
    mVoxelDims = {domain_size/mNumVoxels[0], domain_size/mNumVoxels[1]};
    mInitialVoxelSize = mVoxelDims[0] * mVoxelDims[1];

    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mInitialVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }

}

double Simulation2d::GetVoxelRatio()
{
    return mRatio;
}

void Simulation2d::SetVoxels(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species)
{
    // Check that the input is sensible
    assert(species < mNumSpecies);
    assert(voxel_index.size() == 2);  // We expect location of a single voxel in the input
    assert(voxel_index[0] < mNumVoxels[0]);
    assert(voxel_index[1] < mNumVoxels[1]);

    // Populate the vector mVoxels at the specified position
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run].voxels[species][voxel_index[1]*mNumVoxels[0] + voxel_index[0]] = num_molecules;
    }
    mTotalNumMolecules[species] = vec_sum(mGrids[0].voxels[species]);
}

void Simulation2d::SetVoxels(vector<vector<unsigned> > initial_state, unsigned species)
{
    // Check that the input is sensible
    assert(species < mNumSpecies);
    assert(initial_state[0].size() == mNumVoxels[0]);
    assert(initial_state.size() == mNumVoxels[1]);

    // Now change the array of number of molecules to a single vector, as it is stored in the mGrids
    vector<unsigned> vec;
    for (unsigned vertical_index = 0; vertical_index < mNumVoxels[1]; vertical_index++)
    {
        vec.insert(vec.end(), initial_state[vertical_index].begin(), initial_state[vertical_index].end());
    }

    for (unsigned i=0; i < mNumVoxels[0]*mNumVoxels[1]; i++)
    {
        for (unsigned run=0; run < mNumRuns; run++)
        {
            mGrids[run].voxels[species][i] = vec[i];
        }
    }

    mTotalNumMolecules[species] = vec_sum(mGrids[0].voxels[species]);
}
