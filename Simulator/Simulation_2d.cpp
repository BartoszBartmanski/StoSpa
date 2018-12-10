#include "Simulation_2d.hpp"

Simulation_2d::Simulation_2d(unsigned num_runs, unsigned num_species, string num_method, unsigned num_voxels,
                             vector<double> domain_bounds, string boundary_condition, double ratio, double alpha, double beta_x,
                             double beta_y)
{
    // First check the input parameters
    assert(num_runs > 0);
    assert(num_species > 0);
    assert(num_voxels > 0);
    assert(domain_bounds.size() == 2);
    assert(ratio > 0.0);
    assert(alpha >= 0.0);
    assert(beta_x >= 0.0 and beta_x <= 1.0);
    assert(beta_y >= 0.0 and beta_y <= 1.0);

    // Input parameters
    mNumRuns = num_runs;
    mNumSpecies = num_species;
    mNumVoxels = {num_voxels, unsigned(ratio * num_voxels)};
    mTotalNumVoxels = mNumVoxels[0] * mNumVoxels[1];

    mDomainBounds = domain_bounds;
    mBC = move(boundary_condition);

    // Additional input parameters (due to working in 2d)
    mRatio = mNumVoxels[1]/double(mNumVoxels[0]);  // correction for some values of aspect ratio not being possible
    mAlpha = alpha;
    mBetaX = beta_x;
    mBetaY = beta_y;

    // Simulation attributes that will change with each time step
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    m_h = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[1]);
    mVoxelSize = pow(m_h, 2) * mRatio;
    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);
    mDiffusionCoefficients.resize(mNumSpecies);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mVoxelSize, mNumVoxels[0], mNumVoxels[1]);
    }

    mNumMethod = move(num_method);
    if (mNumMethod == "fdm")
    {
        mJumpRates = make_shared<FDM>(mRatio, m_h, mAlpha);
    }
    else if (mNumMethod == "fem")
    {
        mJumpRates = make_shared<FEM>(mRatio, m_h);
    }
    else if (mNumMethod == "fvm")
    {
        mJumpRates = make_shared<FVM>(mRatio, m_h);
    }
    else if (mNumMethod == "fet")
    {
        mJumpRates = make_shared<FET>(mRatio, m_h, mBetaX, mBetaY, 1000);
    }
    else if (mNumMethod == "fetU")
    {
        mJumpRates = make_shared<FETUniform>(mRatio, m_h, mBetaX, mBetaY, 1000);
    }
    else
    {
        throw runtime_error("Unknown input for numerical method from which to derive the jump coefficients!");
    }

}

double Simulation_2d::GetVoxelRatio()
{
    return mRatio;
}

double Simulation_2d::GetAlpha()
{
    return mAlpha;
}

double Simulation_2d::GetBetaX()
{
    return mBetaX;
}

double Simulation_2d::GetBetaY()
{
    return mBetaY;
}

void Simulation_2d::SetDiffusionRate(double diff, unsigned int species)
{
    // Check for sensible input
    assert(diff >= 0.0);
    assert(species < mNumSpecies);

    mDiffusionCoefficients[species] = diff;

    vector<vector<int>> directions = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
    if (mBC == "reflective")
    {
        for (auto direction : directions)
        {
            shared_ptr<DiffusionReflective> diff_class;
            if (direction[1] == 0)  // Horizontal jumps
            {
                diff_class = make_shared<DiffusionReflective>(diff * mJumpRates->GetLambda1(), species, direction);
            }
            else if (direction[0] == 0)  // Vertical jumps
            {
                diff_class = make_shared<DiffusionReflective>(diff * mJumpRates->GetLambda3(), species, direction);
            }
            else  // Diagonal jumps
            {
                diff_class = make_shared<DiffusionReflective>(diff * mJumpRates->GetLambda2(), species, direction);
            }
            this->AddReaction(diff_class);
        }
    }
    else if (mBC == "periodic")
    {
        for (auto direction : directions)
        {
            shared_ptr<DiffusionPeriodic> diff_class;
            if (direction[1] == 0)  // Horizontal jumps
            {
                diff_class = make_shared<DiffusionPeriodic>(diff * mJumpRates->GetLambda1(), species, direction);
            }
            else if (direction[0] == 0)  // Vertical jumps
            {
                diff_class = make_shared<DiffusionPeriodic>(diff * mJumpRates->GetLambda3(), species, direction);
            }
            else  // Diagonal jumps
            {
                diff_class = make_shared<DiffusionPeriodic>(diff * mJumpRates->GetLambda2(), species, direction);
            }
            this->AddReaction(diff_class);
        }
    }
    else
    {
        throw runtime_error("Boundary condition can only be one of the following: reflective, periodic");
    }

}

void Simulation_2d::SetInitialNumMolecules(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species)
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

void Simulation_2d::SetInitialState(vector<vector<unsigned> > initial_state, unsigned species)
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

void Simulation_2d::SetInitialState(vector<vector<int> > initial_state, unsigned species)
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
