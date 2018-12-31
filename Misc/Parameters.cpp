//
// Created by bartmanski on 26/03/18.
//

#include "Parameters.hpp"

Parameters::Parameters(map<string, docopt::value> cl_input)
{
    if (cl_input.count("--dir_name"))
    {
        if (cl_input["--dir_name"])
            mSaveDir = cl_input["--dir_name"].asString();
    }

    if (cl_input.count("--start_index"))
    {
        if (cl_input["--start_index"])
            mStartIndex = unsigned(stoi(cl_input["--start_index"].asString()));
    }

    if (cl_input.count("--num_dims"))
    {
        if (cl_input["--num_dims"])
            mNumDims = unsigned(stoi(cl_input["--num_dims"].asString()));
    }

    if (cl_input.count("--num_voxels"))
    {
        if (cl_input["--num_voxels"])
            mNumVoxels = unsigned(stoi(cl_input["--num_voxels"].asString()));
    }

    if (cl_input.count("--num_species"))
    {
        if (cl_input["--num_species"])
            mNumSpecies = unsigned(stoi(cl_input["--num_species"].asString()));
    }

    if (cl_input.count("--domain_bounds"))
    {
        if (cl_input["--domain_bounds"])
            mDomainBounds = split<double>(cl_input["--domain_bounds"].asString());
    }

    if (cl_input.count("--initial_num"))
    {
        if (cl_input["--initial_num"])
            mInitialNum = split<unsigned>(cl_input["--initial_num"].asString());
    }

    if (cl_input.count("--num_method"))
    {
        if (cl_input["--num_method"])
            mNumMethod = cl_input["--num_method"].asString();
    }

    if (cl_input.count("--bc"))
    {
        if (cl_input["--bc"])
            mBC = cl_input["--bc"].asString();
    }

    if (cl_input.count("--end_time"))
    {
        if (cl_input["--end_time"])
           mEndTime = stod(cl_input["--end_time"].asString());
    }

    if(cl_input.count("--time_step"))
    {
        if (cl_input["--time_step"])
            mTimeStep = stod(cl_input["--time_step"].asString());
    }

    if (cl_input.count("--diff"))
    {
        if (cl_input["--diff"])
            mDiff = split<double>(cl_input["--diff"].asString());
    }

    if (cl_input.count("--decay"))
    {
        if (cl_input["--decay"])
            mDecay = split<double>(cl_input["--decay"].asString());
    }

    if (cl_input.count("--prod"))
    {
        if (cl_input["--prod"])
            mProd = split<double>(cl_input["--prod"].asString());
    }

    if (cl_input.count("--kappa"))
    {
        if (cl_input["--kappa"])
            mKappa = unsigned(stod(cl_input["--kappa"].asString())*mNumVoxels)/double(mNumVoxels);
    }

    if (cl_input.count("--alpha"))
    {
        if (cl_input["--alpha"])
            mAlpha = stod(cl_input["--alpha"].asString());
    }

    if (cl_input.count("--beta"))
    {
        if (cl_input["--beta"])
        {
            mBeta = split<double>(cl_input["--beta"].asString());
            if (mBeta.size() == 1) { mBeta.push_back(mBeta[0]); }
        }
    }

    if (cl_input.count("--num_runs"))
    {
        if (cl_input["--num_runs"])
            mNumRuns = unsigned(stoi(cl_input["--num_runs"].asString()));
    }

    if (cl_input.count("--trunc_order"))
    {
        if (cl_input["--trunc_order"])
            mTruncOrder = unsigned(stoi(cl_input["--trunc_order"].asString()));
    }

    if (cl_input.count("--num_points"))
    {
        if (cl_input["--end_time"])
            mNumPoints = unsigned(stoi(cl_input["--num_points"].asString()));
    }

    if (cl_input.count("--num_threads"))
    {
        if (cl_input["--num_threads"])
        {
            mNumThreads = unsigned(stoi(cl_input["--num_threads"].asString()));
            mNumThreads = unsigned(gcd(mNumPoints, mNumThreads));
        }
    }
}

void Parameters::Add(string name, string value)
{
    mOther[name] = move(value);
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
            mBeta.clear();
            mKappa=-1;
        }
    }
    if (mNumRuns) { info << "# num_runs = " << mNumRuns << endl; }
    if (mNumSpecies) { info << "# num_species = " << mNumSpecies << endl; }
    if (!mNumMethod.empty()) { info << "# num_method = " << mNumMethod << endl; }
    if (mNumVoxels) { info << "# num_voxels = " << mNumVoxels << endl; }
    if (!mDomainBounds.empty()) { info << "# domain_bounds = " << mDomainBounds[0] << "," << mDomainBounds[1] << endl; }
    if (mKappa > 0) { info << "# kappa = " << mKappa << endl; }
    if (mAlpha >= 0) { info << "# alpha = " << mAlpha << endl; }
    if (!mBeta.empty())
    {
        if (mBeta.size() == 1) { mBeta.push_back(mBeta[0]); }
        info << "# beta = " << mBeta[0] << "," << mBeta[1] << endl;
    }
    if (!mBC.empty()) { info << "# bc = " << mBC << endl; }
    if (mEndTime != 0) { info << "# end_time = " << mEndTime << endl; }
    if (mTimeStep != 0) { info << "# time_step = " << mTimeStep << endl; }
    if (mNumPoints) { info << "# num_points = " << mNumPoints << endl; }
    if (!mInitialNum.empty())
    {
        oss << mInitialNum[0];
        for (unsigned i=1; i < mNumSpecies; i++)
        {
            oss << "," << mInitialNum[i];
        }
        info << "# initial_num_molecules = " << oss.str() << endl;
        oss.clear();
        oss.str("");
    }
    if (mTruncOrder > 0) { info << "# truncation_order = " << mTruncOrder << endl; }
    if (!mDiff.empty())
    {
        mDiff.resize(mNumSpecies);
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
        mDecay.resize(mNumSpecies);
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
        mProd.resize(mNumSpecies);
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

void Parameters::SetNumPoints(unsigned value)
{
    mNumPoints = value;
    mNumThreads = unsigned(gcd(mNumPoints, mNumThreads));
}

unsigned Parameters::GetNumPoints()
{
    return mNumPoints;
}

void Parameters::SetNumThreads(unsigned value)
{
    mNumThreads = value;
    mNumThreads = unsigned(gcd(mNumPoints, mNumThreads));
}

unsigned Parameters::GetNumThreads()
{
    return mNumThreads;
}

void Parameters::SetNumRuns(unsigned value)
{
    mNumRuns = value;
}

unsigned Parameters::GetNumRuns()
{
    return mNumRuns;
}

void Parameters::SetNumDims(unsigned value)
{
    mNumDims = value;
}

unsigned Parameters::GetNumDims()
{
    return mNumDims;
}

void Parameters::SetNumSpecies(unsigned value)
{
    mNumSpecies = value;
}

unsigned Parameters::GetNumSpecies()
{
    return mNumSpecies;
}

void Parameters::SetNumMethod(string value)
{
    mNumMethod = move(value);
}

string Parameters::GetNumMethod()
{
    return mNumMethod;
}

void Parameters::SetNumVoxels(unsigned value)
{
    mNumVoxels = value;
    mKappa = unsigned(mKappa*mNumVoxels)/double(mNumVoxels);
}

unsigned Parameters::GetNumVoxels()
{
    return mNumVoxels;
}

void Parameters::SetDomainBounds(vector<double> values)
{
    mDomainBounds = move(values);
}

vector<double> Parameters::GetDomainBounds()
{
    return mDomainBounds;
}

void Parameters::SetBC(string value)
{
    mBC = move(value);
}

string Parameters::GetBC()
{
    return mBC;
}

void Parameters::SetKappa(double value)
{
    mKappa = value;
    mKappa = unsigned(mKappa*mNumVoxels)/double(mNumVoxels);
}

double Parameters::GetKappa()
{
    return mKappa;
}

void Parameters::SetAlpha(double value)
{
    mAlpha = value;
}

double Parameters::GetAlpha()
{
    return mAlpha;
}

void Parameters::SetBeta(vector<double> values)
{
    mBeta = move(values);
}

vector<double> Parameters::GetBeta()
{
    return mBeta;
}

void Parameters::SetEndTime(double value)
{
    mEndTime = value;
}

double Parameters::GetEndTime()
{
    return mEndTime;
}

void Parameters::SetTimeStep(double value)
{
    mTimeStep = value;
}

double Parameters::GetTimeStep()
{
    return mTimeStep;
}

void Parameters::SetDiff(vector<double> values)
{
    mDiff = move(values);
}

vector<double> Parameters::GetDiff()
{
    return mDiff;
}

void Parameters::SetDecay(vector<double> values)
{
    mDecay = move(values);
}

vector<double> Parameters::GetDecay()
{
    return mDecay;
}

void Parameters::SetProd(vector<double> values)
{
    mProd = move(values);
}

vector<double> Parameters::GetProd()
{
    return mProd;
}

void Parameters::SetInitailNum(vector<unsigned> values)
{
    mInitialNum = move(values);
}

vector<unsigned> Parameters::GetInitialNum()
{
    return mInitialNum;
}

void Parameters::SetTruncOrder(unsigned value)
{
    mTruncOrder = value;
}

unsigned Parameters::GetTruncOrder()
{
    return mTruncOrder;
}

void Parameters::SetSaveDir(string value)
{
    mSaveDir = value;
}

string Parameters::GetSaveDir()
{
    return mSaveDir;
}

void Parameters::SetStartIndex(unsigned value)
{
    mStartIndex = value;
}

unsigned Parameters::GetStartIndex()
{
    return mStartIndex;
}

void Parameters::AddReaction(string name, double rate)
{
    mAdditionalReactions[name] = rate;
}

void Parameters::AddReaction(const unique_ptr<AbstractReaction> &p_reaction)
{
    mAdditionalReactions[p_reaction->GetReactionName()] = p_reaction->GetRateConstant();
}