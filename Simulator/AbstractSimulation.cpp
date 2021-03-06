//
// Created by bartosz on 24/06/17.
//

#include "AbstractSimulation.hpp"

AbstractSimulation::AbstractSimulation() : inf(numeric_limits<double>::infinity()),
                                           mTimesSet(false),
                                           mExtrande(false),
                                           mExtrandeIndex(0),
                                           mGrowth(false),
                                           mNumRuns(1),
                                           mNumSpecies(1),
                                           mTotalNumVoxels(1),
                                           mTime(0),
                                           mInitialVoxelSize(0)
{
    mNumJumps = vector<unsigned>(mNumRuns, 0);
    random_device rd;
    mSeed = rd();
    mGen = mt19937(mSeed);
    mUniform = uniform_real_distribution<double>(0.0, 1.0);

    mDiffusionCoefficients.resize(mNumSpecies);

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

void AbstractSimulation::UpdateVoxelSize(const unsigned& run)
{
    double growth = mGrowthRate->GetFunction(mGrids[run].time);

    if (mDim == 1)
    {
        mGrids[run].voxelSize = mInitialVoxelSize * growth;
        mGrids[run].diffScaleFactor = growth * growth;
    }
    else
    {
        growth *= growth;
        mGrids[run].voxelSize = mInitialVoxelSize * growth;
        mGrids[run].diffScaleFactor = growth;
    }
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

void AbstractSimulation::SSA_loop(const unsigned& run)
{
    // Find the smallest time until the next reaction
    double inv_time = mGrids[run].next_reaction_time.top().first;
    mGrids[run].time = 1.0 / inv_time;
    unsigned voxel_index = mGrids[run].next_reaction_time.top().second;

    if (mGrowth)
    {
        this->UpdateVoxelSize(run);
    }

    if (mGrids[run].time < inf)
    {
        // Determine which reaction happens next and update the molecule numbers accordingly
        unsigned reaction = this->NextReaction(run, voxel_index);
        auto jump_index = unsigned(mReactions[reaction]->UpdateGrid(mGrids[run], voxel_index));

        // Update the times until the next reaction
        this->UpdateTime(run, voxel_index);
        if (jump_index != voxel_index)
        {
            this->UpdateTime(run, jump_index);
        }

        // Update the number of jumps variable
        mNumJumps[run] += 1;
    }
}


void AbstractSimulation::Advance(const double &time_point)
{
    if (!mTimesSet)
    {
        SetupTimeIncrements();
    }

    for (unsigned run=0; run < mNumRuns; run++)
    {
        while (mGrids[run].time < time_point)
        {
            SSA_loop(run);
        }
    }
    mTime = time_point;
}

void AbstractSimulation::SetDiffusionRate(unique_ptr<JumpRate>&& method, double diff, unsigned species)
{
    // Check for sensible input
    assert(diff >= 0.0);
    assert(species < mNumSpecies);

    mDiffusionCoefficients[species] = diff;

    for (auto direction : mJumpDirections)
    {
        if (mBC == "reflective")
        {
            mReactions.emplace_back(make_unique<DiffusionReflective>(diff * method->GetLambda(direction), species, direction));
        }
        else if (mBC == "periodic")
        {
            mReactions.emplace_back(make_unique<DiffusionPeriodic>(diff * method->GetLambda(direction), species, direction));
        }
        else
        {
            throw runtime_error("Boundary condition can only be one of the following: reflective, periodic");
        }
    }
}

void AbstractSimulation::AddReaction(unique_ptr<AbstractReaction>&& reaction)
{
    reaction->CheckNumSpecies(mNumSpecies);
    if (reaction->GetRateConstant() > 0)
    {
        mReactions.emplace_back(move(reaction));
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
    if (!mExtrande)
    {
        // Add the (none -> none) reaction
        mReactions.emplace_back(make_unique<Extrande>());
        mExtrandeIndex = unsigned(mReactions.size()) - 1;

        mExtrande = true;
    }
}

void AbstractSimulation::SetGrowth(unique_ptr<GrowthRate>&& growth)
{
    mGrowth = true;
    mGrowthRate = move(growth);
    this->UseExtrande();
}

vector<unsigned> AbstractSimulation::GetVoxels(unsigned species, unsigned int run)
{
    return mGrids[run].voxels[species];
}

vector<double> AbstractSimulation::GetConcentration(unsigned species)
{
    vector<double> concentration(mTotalNumVoxels);

    for (unsigned run=0; run < mNumRuns; run++)
    {
        concentration = concentration + mGrids[run].voxels[species];
    }
    concentration = concentration / (mNumRuns * mTotalNumMolecules[species] * mInitialVoxelSize);

    return concentration;
}

vector<double> AbstractSimulation::GetAverageNumMolecules(unsigned species)
{
    vector<double> molecules(mTotalNumVoxels);
    for (unsigned run=0; run < mNumRuns; run++)
    {
        molecules = molecules + mGrids[run].voxels[species];
    }

    molecules = molecules / mNumRuns;

    return molecules;
}

double AbstractSimulation::GetCurrentTime()
{
    return mTime;
}

double AbstractSimulation::GetVoxelSize()
{
    double value = mInitialVoxelSize;
    if (mGrowth)
    {
        value = mInitialVoxelSize * mGrowthRate->GetFunction(mTime);
    }

    return value;
}

vector<double> AbstractSimulation::GetVoxelDims()
{
    auto vec = mVoxelDims;
    if (mGrowth)
    {
        auto vec = mGrowthRate->GetFunction(mTime) * mVoxelDims;
    }
    return vec;
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

    error = sqrt(mInitialVoxelSize * error);

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

    mod_u = sqrt(mInitialVoxelSize * mod_u);
    error = sqrt(mInitialVoxelSize * error) / mod_u;

    return error;
}

void AbstractSimulation::Run(const string &output, const double &endtime, const double &timestep)
{
    // Initialise progress object
    auto num_steps = unsigned(endtime/timestep);
    Progress prog(num_steps);

    unique_ptr<ofstream> p_output = make_unique<ofstream>(output, ios::app);

    // Run the SSA
    for (unsigned i=0; i < num_steps; i++)
    {
        // Move to the next time step
        this->Advance(i * timestep);

        for (unsigned species=0; species < mNumSpecies; species++)
        {
            // Save the stochastic simulation state
            save(this->GetAverageNumMolecules(species), p_output);
        }
        prog.Show();
    }

    prog.End(output);
}