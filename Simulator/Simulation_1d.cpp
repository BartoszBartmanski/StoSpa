
#include "Simulation_1d.hpp"

Simulation_1d::Simulation_1d(unsigned num_runs,
                             unsigned num_species,
                             string num_method,
                             unsigned num_voxels,
                             vector<double> domain_bounds,
                             string boundary_condition)
{
    // First check the input parameters
    assert(num_runs > 0);
    assert(num_species > 0);
    assert(num_voxels > 0);
    assert(domain_bounds.size() == 2);

    // Input parameters
    mNumRuns = num_runs;
    mNumSpecies = num_species;
    mNumMethod = move(num_method);
    mNumVoxels = {num_voxels, 1};
    mTotalNumVoxels = num_voxels;
    mDomainBounds = domain_bounds;
    mBC = move(boundary_condition);

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    m_h = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[0]);
    mVoxelSize = m_h;
    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }
}

Simulation_1d::Simulation_1d(Parameters params)
{
    // First check the input parameters
    assert(params.GetNumRuns() > 0);
    assert(params.GetNumSpecies() > 0);
    assert(params.GetNumVoxels() > 0);
    assert(params.GetDomainBounds().size() == 2);

    // Input parameters
    mNumRuns = params.GetNumRuns();
    mNumSpecies = params.GetNumSpecies();
    mNumMethod = params.GetNumMethod();
    mNumVoxels = {params.GetNumVoxels(), 1};
    mTotalNumVoxels = mNumVoxels[0];
    mDomainBounds = params.GetDomainBounds();
    mBC = params.GetBC();

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    m_h = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[0]);
    mVoxelSize = m_h;
    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }
}

void Simulation_1d::SetDiffusionRate(double diff, unsigned int species)
{
    // Check for sensible input
    assert(diff >= 0.0);
    assert(species < mNumSpecies);

    mDiffusionCoefficients[species] = diff;
    double lambda = diff / pow(m_h, 2);

    vector<vector<int>> directions = {{-1}, {1}};
    if (mBC == "reflective")
    {
        mReactions.emplace_back(make_unique<DiffusionReflective>(lambda, species, directions[0]));
        mReactions.emplace_back(make_unique<DiffusionReflective>(lambda, species, directions[1]));
    }
    else if (mBC == "periodic")
    {
        mReactions.emplace_back(make_unique<DiffusionPeriodic>(lambda, species, directions[0]));
        mReactions.emplace_back(make_unique<DiffusionPeriodic>(lambda, species, directions[1]));
    }
    else if (mBC == "Exponential")
    {
        mReactions.emplace_back(make_unique<DiffusionReflectiveExp>(lambda, species, directions[0], 1.0));
        mReactions.emplace_back(make_unique<DiffusionReflectiveExp>(lambda, species, directions[1], 1.0));
    }
    else
    {
        throw runtime_error("Boundary condition can only be one of the following: reflective, periodic");
    }
}

void Simulation_1d::SetInitialNumMolecules(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species)
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

void Simulation_1d::SetInitialState(vector<vector<unsigned>> initial_state, unsigned species)
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

void Simulation_1d::SetInitialState(vector<vector<int>> initial_state, unsigned species)
{
    // Check that the input is sensible
    assert(species < mNumSpecies);
    assert(initial_state[0].size() == mNumVoxels[0]);

    for (unsigned i=0; i < initial_state[0].size(); i++)
    {
        for (unsigned run=0; run < mNumRuns; run++)
        {
            mGrids[run].voxels[species][i] = unsigned(initial_state[0][i]);
        }
    }

    mTotalNumMolecules[species] = vec_sum(mGrids[0].voxels[species]);
}