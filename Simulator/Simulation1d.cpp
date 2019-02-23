
#include "Simulation1d.hpp"

Simulation1d::Simulation1d(unsigned num_runs, unsigned num_species, unsigned num_voxels, vector<double> domain_bounds,
                             string boundary_condition)
{
    // First check the input parameters
    assert(num_runs > 0);
    assert(num_species > 0);
    assert(num_voxels > 0);
    assert(domain_bounds.size() == 2);

    // Input parameters
    mDim = 1;
    mNumRuns = num_runs;
    mNumSpecies = num_species;
    mNumVoxels = {num_voxels, 1};
    mTotalNumVoxels = num_voxels;
    mDomainBounds = move(domain_bounds);
    mBC = move(boundary_condition);

    // Defines the directions in which molecules can jump
    mJumpDirections = {{-1,0}, {1,0}};

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    mInitialVoxelSize = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[0]);
    mVoxelDims = {mInitialVoxelSize};

    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mInitialVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }
}

Simulation1d::Simulation1d(Parameters params)
{
    // First check the input parameters
    assert(params.GetNumRuns() > 0);
    assert(params.GetNumSpecies() > 0);
    assert(params.GetNumVoxels() > 0);
    assert(params.GetDomainBounds().size() == 2);

    // Input parameters
    mDim = 1;
    mNumRuns = params.GetNumRuns();
    mNumSpecies = params.GetNumSpecies();
    mNumVoxels = {params.GetNumVoxels(), 1};
    mTotalNumVoxels = mNumVoxels[0];
    mDomainBounds = params.GetDomainBounds();
    mBC = params.GetBC();

    // Defines the directions in which molecules can jump
    mJumpDirections = {{-1,0}, {1,0}};

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    mInitialVoxelSize = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[0]);
    mVoxelDims = {mInitialVoxelSize};

    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mInitialVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }
}

void Simulation1d::SetVoxels(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species)
{
    // Check that the input is sensible
    assert(species < mNumSpecies);
    assert(voxel_index[0] < mNumVoxels[0]);
    assert(0 < num_molecules);

    // Populate the vector mVoxels at the specified position
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run].voxels[species][voxel_index[0]] = num_molecules;
    }
    mTotalNumMolecules[species] = vec_sum(mGrids[0].voxels[species]);
}

void Simulation1d::SetVoxels(vector<vector<unsigned>> initial_state, unsigned species)
{
    // Check that the input is sensible
    assert(species < mNumSpecies);
    assert(initial_state[0].size() == mNumVoxels[0]);

    for (unsigned i=0; i < initial_state[0].size(); i++)
    {
        for (unsigned run=0; run < mNumRuns; run++)
        {
            mGrids[run].voxels[species][i] = initial_state[0][i];
        }
    }

    mTotalNumMolecules[species] = vec_sum(mGrids[0].voxels[species]);
}
