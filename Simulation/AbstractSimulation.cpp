//
// Created by bartosz on 24/06/17.
//

#include "AbstractSimulation.hpp"

AbstractSimulation::AbstractSimulation() : inf(numeric_limits<double>::infinity()),
                                           mTimesSet(false),
                                           mNumRuns(1),
                                           mNumSpecies(1),
                                           mNumMethod("fdm"),
                                           mTotalNumVoxels(1),
                                           mTime(0),
                                           mNumJumps(0),
                                           m_h(0),
                                           mVoxelSize(0)
{
    // TODO: add a seed to rd
    random_device rd;
    mSeed = rd();
    mGen = mt19937(mSeed);
    mUniform = uniform_real_distribution<double>(0.0, 1.0);

    mDiffusionCoefficients.resize(mNumSpecies);

//    shared_ptr<None> none = make_shared<None>(1.0, 0);
//    mReactions.push_back(none);
}

void AbstractSimulation::AddReaction(shared_ptr<AbstractReaction> reaction)
{
    reaction->CheckNumSpecies(mNumSpecies);
    if (reaction->GetRateConstant() > 0)
    {
        mReactions.push_back(reaction);
    }
}

vector<shared_ptr<AbstractReaction>> AbstractSimulation::GetReactions()
{
    return mReactions;
}

inline double AbstractSimulation::GetTotalPropensity(const unsigned& run, const int& voxel_index)
{
    double propensity = 0;

    for (const shared_ptr<AbstractReaction>& reaction : mReactions)
    {
        propensity += reaction->GetPropensity(mGrids[run], voxel_index);
    }

    return propensity;
}

unsigned AbstractSimulation::NextReaction(const unsigned& run, const int& voxel_index)
{
    double r_a_0 = 0;
    for (unsigned i=0; i < mReactions.size(); i++)
    {
        double propensity = mReactions[i]->GetPropensity(mGrids[run], voxel_index);
        mPropensities[i] = propensity;
        r_a_0 += propensity;
    }
    r_a_0 *= mUniform(mGen);

    unsigned reaction = 0;
    double lower_bound = 0;
    for (const double& propensity : mPropensities)
    {
        if (r_a_0 > lower_bound and r_a_0 < lower_bound + propensity)
        {
            break;
        }
        else
        {
            reaction += 1;
            lower_bound += propensity;
        }
    }

    return reaction;
}

double AbstractSimulation::Exponential(double propensity)
{
    return (-1.0/propensity) * log(mUniform(mGen));
}

void AbstractSimulation::SetupTimeIncrements()
{
    mTimesSet = true;
    mPropensities.resize(mReactions.size());

    pair<double, unsigned> a_pair;
    for (unsigned run=0; run < mNumRuns; run++)
    {
        for (unsigned i = 0; i < mTotalNumVoxels; i++)
        {
            double t_0 = 1.0 / Exponential(GetTotalPropensity(run, i));
            a_pair = make_pair(t_0, i);
            mGrids[run].handles[i] = mGrids[run].time_increments.push(a_pair);
        }
    }
}

void AbstractSimulation::SetSeed(unsigned number)
{
    mSeed = number;
    mGen = mt19937(number);
}

unsigned AbstractSimulation::GetSeed()
{
    return mSeed;
}

void AbstractSimulation::SSA_loop(const unsigned& run)
{
    // Find the smallest time until the next reaction
    double inv_time = mGrids[run].time_increments.top().first;
    mCurrentTime[run] = 1.0 / inv_time;
    int voxel_index = mGrids[run].time_increments.top().second;

    if (mCurrentTime[run] < inf)
    {
        // Determine which reaction happens next and update the molecule numbers accordingly
        unsigned reaction = NextReaction(run, voxel_index);
        int jump_index = mReactions[reaction]->UpdateGrid(mGrids[run], voxel_index);

        // Update the times until the next reaction
        inv_time = 1.0 / (mCurrentTime[run] + Exponential(GetTotalPropensity(run, voxel_index)));
        pair<double, unsigned> new_pair = make_pair(inv_time, voxel_index);
        *mGrids[run].handles[voxel_index] = new_pair;
        mGrids[run].time_increments.update(mGrids[run].handles[voxel_index]);

        if (jump_index != voxel_index)
        {
            inv_time = 1.0/(mCurrentTime[run] + Exponential(GetTotalPropensity(run, jump_index)));
            new_pair = make_pair(inv_time, jump_index);
            *mGrids[run].handles[jump_index] = new_pair;
            mGrids[run].time_increments.update(mGrids[run].handles[jump_index]);
        }

        // Update the number of jumps variable
        mNumJumps += 1;
    }

}

void AbstractSimulation::Advance(const double& time_step, const unsigned& iterator)
{
    if (!mTimesSet)
    {
        SetupTimeIncrements();
    }

    for (unsigned run=0; run < mNumRuns; run++)
    {
        while (mCurrentTime[run] < iterator * time_step)
        {
            SSA_loop(run);
        }
    }
    mTime = iterator * time_step;
}

vector<unsigned> AbstractSimulation::GetVoxels(unsigned int species, unsigned int run)
{
    return mGrids[run].voxels[species];
}

vector<double> AbstractSimulation::GetConcentration(unsigned int species)
{
    vector<double> concentration(mTotalNumVoxels);

    for (unsigned run=0; run < mNumRuns; run++)
    {
        concentration = concentration + mGrids[run].voxels[species];
    }
    concentration = concentration / (mNumRuns * mTotalNumMolecules[species] * mVoxelSize);

    return concentration;
}

vector<double> AbstractSimulation::GetAverageNumMolecules(unsigned int species)
{
    vector<double> molecules(mTotalNumVoxels);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        molecules = molecules + mGrids[run].voxels[species];
    }

    molecules = molecules / mNumRuns;

    return molecules;
}

vector<unsigned> AbstractSimulation::GetNumMolecules(unsigned species, unsigned run)
{
    assert(species < mNumSpecies);
    return mGrids[run].voxels[species];
}

double AbstractSimulation::GetCurrentTime()
{
    return mTime;
}

double AbstractSimulation::GetSpacing()
{
    return m_h;
}

double AbstractSimulation::GetVoxelSize()
{
    return mVoxelSize;
}

unsigned AbstractSimulation::GetInitialTotalMolecules(unsigned int species)
{
    return mTotalNumMolecules[species];
}

unsigned AbstractSimulation::GetTotalMolecules(unsigned int species, unsigned int run)
{
    unsigned total = 0;
    for (const unsigned& voxel : mGrids[run].voxels[species])
    {
        total += voxel;
    }
    return total;
}

unsigned AbstractSimulation::GetNumJumps()
{
    return mNumJumps;
}

vector<unsigned> AbstractSimulation::GetNumVoxels()
{
    return mNumVoxels;
}

unsigned AbstractSimulation::GetNumRuns()
{
    return mNumRuns;
}

unsigned AbstractSimulation::GetNumSpecies()
{
    return mNumSpecies;
}

string AbstractSimulation::GetNumMethod()
{
    return mNumMethod;
}

vector<double> AbstractSimulation::GetDomainBounds()
{
    return mDomainBounds;
}

string AbstractSimulation::GetBoundaryCondition()
{
    return mBC;
}

double AbstractSimulation::GetDiffusionCoefficient(unsigned int species)
{
    assert(species < mNumSpecies);
    return mDiffusionCoefficients[species];
}

double AbstractSimulation::GetError(const vector<double>& analytic, unsigned species)
{
    vector<double> u_stochastic = GetConcentration(species);

    double error = 0.0;

    for (unsigned i=0; i < u_stochastic.size(); i++)
    {
        error += (u_stochastic[i] - analytic[i]) * (u_stochastic[i] - analytic[i]);
    }

    error = sqrt(mVoxelSize * error);

    return error;
}

double AbstractSimulation::GetRelativeError(const vector<double>& analytic, unsigned species)
{

    double error = 0.0;
    double mod_u = 0.0;

    vector<double> u_stochastic = GetConcentration(species);

    for (unsigned i=0; i < u_stochastic.size(); i++)
    {
        mod_u += analytic[i] * analytic[i];
        error += (u_stochastic[i] - analytic[i]) * (u_stochastic[i] - analytic[i]);
    }

    mod_u = sqrt(mVoxelSize * mod_u);
    error = sqrt(mVoxelSize * error) / mod_u;

    return error;
}
