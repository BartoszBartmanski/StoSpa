//
// Created by bartosz on 24/06/17.
//

#include "AbstractSimulation.hpp"

AbstractSimulation::AbstractSimulation() : inf(numeric_limits<double>::infinity()),
                                           mTimesSet(false),
                                           mExtrande(false),
                                           mExtrandeIndex(0),
                                           mNumRuns(1),
                                           mNumSpecies(1),
                                           mNumMethod("fdm"),
                                           mTotalNumVoxels(1),
                                           mTime(0),
                                           m_h(0),
                                           mVoxelSize(0)
{
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    random_device rd;
    mSeed = rd();
    mGen = mt19937(mSeed);
    mUniform = uniform_real_distribution<double>(0.0, 1.0);

    mDiffusionCoefficients.resize(mNumSpecies);

}

void AbstractSimulation::AddReaction(unique_ptr<AbstractReaction>&& reaction)
{
    reaction->CheckNumSpecies(mNumSpecies);
    if (reaction->GetRateConstant() > 0)
    {
        mReactions.emplace_back(move(reaction));
    }
}

unsigned AbstractSimulation::NextReaction(const unsigned& run, const int& voxel_index)
{
    double r_a_0 = mUniform(mGen) * mGrids[run].a_0[voxel_index];

    unsigned reaction_idx = 0;
    double lower_bound = 0;
    for (const unique_ptr<AbstractReaction>& reaction : mReactions)
    {
        double propensity = reaction->GetPropensity(mGrids[run], voxel_index);
        if (r_a_0 > lower_bound and r_a_0 < lower_bound + propensity)
        {
            break;
        }
        else
        {
            reaction_idx += 1;
            lower_bound += propensity;
        }
    }

    if (mExtrande and reaction_idx == mReactions.size())
    {
        reaction_idx = mExtrandeIndex;
    }

    if (reaction_idx >= mReactions.size())
    {
        throw runtime_error("NextReaction function returns index outside of possible range");
    }

    return reaction_idx;
}

double AbstractSimulation::Exponential(const double& propensity)
{
    return (-1.0/propensity) * log(mUniform(mGen));
}

void AbstractSimulation::UpdateTotalPropensity(const unsigned& run, const int& voxel_index)
{
    if (mExtrande)
    {
        double max_next_prop = 0;
        for (const auto& reaction : mReactions)
        {
            double future = reaction->GetFuturePropensity(mGrids[run], voxel_index);
            if (future > max_next_prop)
            {
                max_next_prop = future;
            }
        }
        mReactions[mExtrandeIndex]->SetRateConstant(max_next_prop);
    }

    double total = 0;
    for (const unique_ptr<AbstractReaction>& reaction : mReactions)
    {
        total += reaction->GetPropensity(mGrids[run], voxel_index);
    }

   mGrids[run].a_0[voxel_index] = total;
}

void AbstractSimulation::UpdateTime(const unsigned& run, const int& voxel_index)
{
    UpdateTotalPropensity(run, voxel_index);
    double inv_time = 1.0 / (mGrids[run].time + Exponential(mGrids[run].a_0[voxel_index]));
    pair<double, unsigned> new_pair = make_pair(inv_time, voxel_index);
    *mGrids[run].handles[voxel_index] = new_pair;
    mGrids[run].next_reaction_time.update(mGrids[run].handles[voxel_index]);
}

void AbstractSimulation::SetupTimeIncrements()
{
    if (!mTimesSet)
    {
        pair<double, unsigned> a_pair;
        for (unsigned run = 0; run < mNumRuns; run++)
        {
            for (unsigned i = 0; i < mTotalNumVoxels; i++)
            {
                UpdateTotalPropensity(run, i);
                double t_0 = 1.0 / Exponential(mGrids[run].a_0[i]);
                a_pair = make_pair(t_0, i);
                mGrids[run].handles[i] = mGrids[run].next_reaction_time.push(a_pair);
            }
        }
        mTimesSet = true;
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

void AbstractSimulation::UseExtrande()
{
    mExtrande = true;
    // Add the (none -> none) reaction
    unique_ptr<Extrande> none = unique_ptr<Extrande>();
    mReactions.emplace_back(move(none));
    mExtrandeIndex = unsigned(mReactions.size()) - 1;
}

void AbstractSimulation::SSA_loop(const unsigned& run)
{
    // Find the smallest time until the next reaction
    double inv_time = mGrids[run].next_reaction_time.top().first;
    mGrids[run].time = 1.0 / inv_time;
    unsigned voxel_index = mGrids[run].next_reaction_time.top().second;

    if (mGrids[run].time < inf)
    {
        // Determine which reaction happens next and update the molecule numbers accordingly
        unsigned reaction = NextReaction(run, voxel_index);
        auto jump_index = unsigned(mReactions[reaction]->UpdateGrid(mGrids[run], voxel_index));

        // Update the times until the next reaction
        UpdateTime(run, voxel_index);
        if (jump_index != voxel_index)
        {
            UpdateTime(run, jump_index);
        }

        // Update the number of jumps variable
        mNumJumps[run] += 1;
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
        while (mGrids[run].time < iterator * time_step)
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

unsigned AbstractSimulation::GetNumJumps(unsigned run)
{
    return mNumJumps[run];
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
