//
// Created by bartmanski on 26/03/18.
//

#include "Parameters.hpp"

Parameters::Parameters()
{
    mComments = "Below are all the parameters used to generate the data contained in this file.";
    mNumThreads = 0;
    mNumRuns = 0;
    mNumDims = 0;
    mNumSpecies = 0;
    mNumVoxels = 0;
    m_h = 0;
    mKappa = 0;
    mEndTime = 0;
    mTimeStep = 0;
    mNumPoints = 0;
    mTruncOrder = 0;
}

void Parameters::Add(const string& name, const string& value)
{
    mOther[name] = value;
}

void Parameters::Save(const string& path_to_file)
{
    ostringstream oss;

    ofstream info(path_to_file);
    if (!mComments.empty()) { info << "# " << mComments << endl; }
    info << "# git_version = " << GIT_VERSION << endl;
    if (!mCommand.empty()) { info << "# command = " << mCommand << endl; }
    if (mNumThreads) { info << "# num_threads = " << mNumThreads << endl; }
    if (mNumDims)
    {
        info << "# num_dims = " << mNumDims << endl;
        if (mNumDims==1)
        {
            mAlpha.clear();
            mBeta.clear();
            mKappa=0;
        }
    }
    if (mNumRuns) { info << "# num_runs = " << mNumRuns << endl; }
    if (mNumSpecies) { info << "# num_species = " << mNumSpecies << endl; }
    if (!mNumMethod.empty()) { info << "# num_method = " << mNumMethod << endl; }
    if (mNumVoxels) { info << "# num_voxels = " << mNumVoxels << endl; }
    if (!mDomainBounds.empty()) { info << "# domain_bounds = " << mDomainBounds[0] << "," << mDomainBounds[1] << endl; }
    if (m_h) { info << "# h = " << m_h << endl; }
    if (mKappa) { info << "# kappa = " << mKappa << endl; }
    if (!mAlpha.empty()) { info << "# alpha = " << mAlpha[0] << endl; }
    if (!mBeta.empty())
    {
        if (mBeta.size() == 1) { mBeta.push_back(mBeta[0]); }
        info << "# beta = " << mBeta[0] << "," << mBeta[1] << endl;
    }
    if (!mBC.empty()) { info << "# bc = " << mBC << endl; }
    if (mEndTime) { info << "# end_time = " << mEndTime << endl; }
    if (mTimeStep) { info << "# time_step = " << mTimeStep << endl; }
    if (mNumPoints) { info << "# num_points = " << mNumPoints << endl; }
    if (!mInitialNumMolecules.empty())
    {
        oss << mInitialNumMolecules[0];
        for (unsigned i=1; i < mNumSpecies; i++)
        {
            oss << "," << mInitialNumMolecules[i];
        }
        info << "# initial_num_molecules = " << oss.str() << endl;
        oss.clear();
        oss.str("");
    }
    if (mTruncOrder) { info << "# truncation_order = " << mTruncOrder << endl; }
    if (!mDiff.empty())
    {
        oss << mDiff[0];
        for (unsigned i=1; i < mNumSpecies; i++)
        {
            oss << "," << mDiff[i];
        }
        info << "# diffusion_constants = " << oss.str() << endl;
        oss.clear();
        oss.str("");
    }
    if (!mDecay.empty())
    {
        oss << mDecay[0];
        for (unsigned i=1; i < mNumSpecies; i++)
        {
            oss << "," << mDecay[i];
        }
        info << "# decay_constants = " << oss.str() << endl;
        oss.clear();
        oss.str("");
    }
    if (!mProd.empty())
    {
        oss << mProd[0];
        for (unsigned i=1; i < mNumSpecies; i++)
        {
            oss << "," << mProd[i];
        }
        info << "# prod_constants = " << oss.str() << endl;
        oss.clear();
        oss.str("");
    }
    if (!mAdditionalReactions.empty())
    {
        info << "# additional_reactions" << endl;
        for (const auto& entry : mAdditionalReactions)
        {
            info << "#\t " << entry.first << " = " << entry.second << endl;
        }
    }
    for (const auto& entry : mOther)
    {
        info << "# " << entry.first << " = " << entry.second << endl;
    }

    info.close();
}

void Parameters::SetComments(string comment)
{
    mComments = move(comment);
}

void Parameters::SetCommand(string command)
{
    mCommand = move(command);
}

void Parameters::SetNumThreads(unsigned num_threads)
{
    mNumThreads = num_threads;
}

void Parameters::SetNumRuns(unsigned num_runs)
{
    mNumRuns = num_runs;
}

void Parameters::SetNumDims(unsigned num_dims)
{
    mNumDims = num_dims;
}

void Parameters::SetNumSpecies(unsigned num_species)
{
    mNumSpecies = num_species;
}

void Parameters::SetNumMethod(string num_method)
{
    mNumMethod = move(num_method);
}

void Parameters::SetNumVoxels(unsigned num_voxels)
{
    mNumVoxels = num_voxels;
}

void Parameters::SetDomainBounds(vector<double> domain_bounds)
{
    mDomainBounds = move(domain_bounds);
}

void Parameters::SetBC(string bc)
{
    mBC = move(bc);
}

void Parameters::SetKappa(double kappa)
{
    mKappa = kappa;
}

void Parameters::SetAlpha(double alpha)
{
    mAlpha = {alpha};
}

void Parameters::SetBeta(vector<double> beta)
{
    mBeta = move(beta);
}

void Parameters::SetEndTime(double end_time)
{
    mEndTime = end_time;
}

void Parameters::SetTimeStep(double time_step)
{
    mTimeStep = time_step;
}

void Parameters::SetNumPoints(unsigned num_points)
{
    mNumPoints = num_points;
}

void Parameters::SetDiff(vector<double> diff)
{
    mDiff = move(diff);
}

void Parameters::SetDiff(double diff)
{
    mDiff = {diff};
}

void Parameters::SetDecay(vector<double> decay)
{
    mDecay = move(decay);
}

void Parameters::SetDecay(double decay)
{
    mDecay = {decay};
}

void Parameters::SetProd(vector<double> prod)
{
    mProd = move(prod);
}

void Parameters::SetProd(double prod)
{
    mProd = {prod};
}

void Parameters::SetInitialNumMolecules(unsigned num_molecules)
{
    mInitialNumMolecules = {num_molecules};
}

void Parameters::SetInitialNumMolecules(vector<unsigned> num_molecules)
{
    mInitialNumMolecules = move(num_molecules);
}

void Parameters::SetTruncOrder(unsigned trunc_order)
{
    mTruncOrder = trunc_order;
}

void Parameters::AddAdditionalReactions(string name, double rate)
{
    mAdditionalReactions[name] = rate;
}
