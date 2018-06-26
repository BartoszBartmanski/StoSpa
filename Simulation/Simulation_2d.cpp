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

    if (num_method != "fdm" and num_method != "fvm" and num_method != "fem" and num_method != "fet")
    {
        throw runtime_error("Parameter num_method can only be one of the following:\n- fdm\n- fvm\n- fem\n- fet");
    }

    if (boundary_condition != "reflective" and boundary_condition != "periodic")
    {
        throw runtime_error("Parameter boundary_condition can only be reflective or periodic");
    }

    // Input parameters
    mNumRuns = num_runs;
    mNumSpecies = num_species;
    mNumMethod = num_method;
    mNumVoxels = {num_voxels, unsigned(ratio * num_voxels)};
    mTotalNumVoxels = mNumVoxels[0] * mNumVoxels[1];

    mDomainBounds = domain_bounds;
    mBC = boundary_condition;

    // Additional input parameters (due to working in 2d)
    mRatio = mNumVoxels[1]/double(mNumVoxels[0]);  // correction for some values of aspect ratio not being possible
    mAlpha = alpha;
    mBetaX = beta_x;
    mBetaY = beta_y;

    // Simulation attributes that will change with each time step
    mCurrentTime = vector<double>(mNumRuns, 0);
    mNumJumps = 0;
    mTime = 0.0;

    // Simulation attributes that will remain constant throughout the simulation
    m_h = (mDomainBounds[1] - mDomainBounds[0]) / double(mNumVoxels[1]);
    mVoxelSize = pow(m_h, 2) * mRatio;
    mTotalNumMolecules = vector<unsigned>(mNumSpecies, 0);

    mGrids.resize(mNumRuns);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        mGrids[run] = Grid(mNumSpecies, mVoxelSize, mNumVoxels[0], mNumVoxels[1]);
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

double Simulation_2d::GetLambdaDiffusion(unsigned int species, unsigned truncation_order)
{
    double lambda_diffusion = 0;
    // First we calculate the expected time to leave the voxel, then we invert this to calculate the propensity
    for (unsigned j = 1; j < truncation_order + 1; j++)
    {
        for (unsigned k = 1; k < truncation_order + 1; k++)
        {
            lambda_diffusion += (pow(-1, j + k)/((2.0*j-1.0) * (2.0*k-1.0) * (pow(mRatio*(2.0*j-1.0), 2) + pow((2.0*k-1.0), 2))));
        }
    }

    lambda_diffusion *= (64.0 * pow(mRatio * m_h, 2) / (pow(M_PI, 4) * mDiffusionCoefficients[species]));

    lambda_diffusion = 1.0 / lambda_diffusion;

    return lambda_diffusion;
}

double Simulation_2d::GetTheta1(unsigned truncation_order)
{
    double theta1 = 0.0;

    for (unsigned j = 1; j < truncation_order + 1; j++)
    {
        for (unsigned k = 1; k < truncation_order + 1; k++)
        {
            double factor = 4.0 * pow(mRatio * m_h, 2) / (pow(M_PI, 2) * (pow(mRatio * (2.0*j-1.0), 2) + pow((2.0*k-1.0), 2)));
            theta1 += ((2.0 * pow(-1, k+1) * (2.0*k-1.0) * sin(0.5*(2.0*j-1.0)*M_PI*mBetaY) * factor) / (pow(mRatio * m_h, 2) * (2.0*j-1.0)));
        }
    }

    return theta1;
}

double Simulation_2d::GetTheta2(unsigned truncation_order)
{
    double theta2 = 0.0;

    for (unsigned j = 1; j < truncation_order + 1; j++)
    {
        for (unsigned k = 1; k < truncation_order + 1; k++)
        {
            double factor = 4.0 * pow(mRatio * m_h, 2) / (pow(M_PI, 2) * (pow(mRatio * (2.0*j-1.0), 2) + pow((2.0*k-1.0), 2)));
            theta2 += ((pow(-1, j+k) * (2.0*k-1.0) * (sin(M_PI*(mBetaY * j - 0.5 * mBetaY + j)) + 1.0) * factor) / (pow(m_h * mRatio, 2) * (2.0*j-1.0)));
            theta2 += ((pow(-1, j+k) * (2.0*j-1.0) * (sin(M_PI*(mBetaX * k - 0.5 * mBetaX + k)) + 1.0) * factor) / (pow(m_h, 2) * (2.0*k-1.0)));
        }
    }

    return theta2;
}

double Simulation_2d::GetTheta3(unsigned truncation_order)
{
    double theta3 = 0.0;

    for (unsigned j = 1; j < truncation_order + 1; j++)
    {
        for (unsigned k = 1; k < truncation_order + 1; k++)
        {
            double factor = 4.0 * pow(mRatio * m_h, 2) / (pow(M_PI, 2) * (pow(mRatio * (2.0*j-1.0), 2) + pow((2.0*k-1.0), 2)));
            theta3 += ((2.0 * pow(-1, j+1) * (2.0*j-1.0) * sin(0.5*(2.0*k-1.0)*M_PI*mBetaX) * factor) / (pow(m_h, 2) * (2.0*k-1.0)));
        }
    }

    return theta3;
}

void Simulation_2d::SetDiffusionRate(double diffusion_coefficient, unsigned int species)
{
    // Check for sensible input
    assert(diffusion_coefficient >= 0.0);
    assert(species < mNumSpecies);

    mDiffusionCoefficients[species] = diffusion_coefficient;

    vector<double> lambda(3, 0);
    if (mNumMethod == "fdm")
    {
        // Check that the parameters satisfy the discrete maximum principle (which means that the probabilities are non-negative)
        double upper_bound = min(mRatio, 1.0/mRatio);
        assert(mAlpha >= 0.0 and mAlpha <= upper_bound and upper_bound <= 1.0);

        lambda[0] = diffusion_coefficient * (1.0 - mRatio * mAlpha) / pow(mRatio * m_h, 2);
        lambda[1] = diffusion_coefficient * (0.5 * mRatio *mAlpha) / pow(mRatio * m_h, 2);
        lambda[2] = diffusion_coefficient * (pow(mRatio, 2) - mRatio * mAlpha) / pow(mRatio * m_h, 2);

    }
    else if (mNumMethod == "fem")
    {
        // Check that the parameters satisfy the discrete maximum principle
        assert(mRatio >= 1.0/sqrt(2) and mRatio <= sqrt(2));

        // This method is equivalent to fdm when the free parameter is as follows
        mAlpha = (pow(mRatio, 2) + 1.0) / (3.0 * mRatio);

        lambda[0] = diffusion_coefficient * (2.0 - pow(mRatio, 2)) / (3.0 * pow(mRatio * m_h, 2));
        lambda[1] = diffusion_coefficient * (pow(mRatio, 2) + 1.0) / (6.0 * pow(mRatio * m_h, 2));
        lambda[2] = diffusion_coefficient * (2.0 * pow(mRatio, 2) - 1.0) / (3.0 * pow(mRatio * m_h, 2));
    }
    else if (mNumMethod == "fvm")
    {
        lambda[0] = diffusion_coefficient / pow(mRatio * m_h, 2);
        lambda[1] = 0.0;
        lambda[2] = diffusion_coefficient / pow(m_h, 2);
    }
    else if (mNumMethod == "fet")
    {
        double lambda_diffusion = GetLambdaDiffusion(species, 1000);
        double theta_1 = GetTheta1(1000);
        lambda[0] = lambda_diffusion * theta_1;
        double theta_3 = GetTheta3(1000);
        lambda[2] = lambda_diffusion * theta_3;
        double theta_2 = 0.25 * (1.0 - 2.0*theta_1 - 2.0 * theta_3);
        lambda[1] = lambda_diffusion * theta_2;
    }
    else
    {
        throw runtime_error("Unknown input for numerical method from which to derive the jump coefficients!");
    }

    vector<vector<int>> directions = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
    if (mBC == "reflective")
    {
        for (auto direction : directions)
        {
            shared_ptr<DiffusionReflective> diff;
            if (direction[1] == 0)  // Horizontal jumps
            {
                diff = make_shared<DiffusionReflective>(lambda[0], species, direction);
            }
            else if (direction[0] == 0)  // Vertical jumps
            {
                diff = make_shared<DiffusionReflective>(lambda[2], species, direction);
            }
            else  // Diagonal jumps
            {
                diff = make_shared<DiffusionReflective>(lambda[1], species, direction);
            }
            this->AddReaction(diff);
        }
    }
    else
    {
        for (auto direction : directions)
        {
            shared_ptr<DiffusionPeriodic> diff;
            if (direction[1] == 0)  // Horizontal jumps
            {
                diff = make_shared<DiffusionPeriodic>(lambda[0], species, direction);
            }
            else if (direction[0] == 0)  // Vertical jumps
            {
                diff = make_shared<DiffusionPeriodic>(lambda[2], species, direction);
            }
            else  // Diagonal jumps
            {
                diff = make_shared<DiffusionPeriodic>(lambda[1], species, direction);
            }
            this->AddReaction(diff);
        }
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
