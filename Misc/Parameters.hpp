//
// Created by bartmanski on 26/03/18.
//

#ifndef STOSPA_PARAMETERS_HPP
#define STOSPA_PARAMETERS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iterator>
#include "Version.hpp"
#include "Utilities.hpp"
#include "docopt.h"

using namespace std;

class Parameters
{
private:
    string mComments;

    string mCommand;

    unsigned mNumPoints=1;

    unsigned mNumThreads=1;

    unsigned mNumRuns=1;

    unsigned mNumDims=1;

    unsigned mNumSpecies=1;

    string mNumMethod;

    unsigned mNumVoxels=1;

    vector<double> mDomainBounds;

    string mBC;

    double mKappa=-1;

    double mAlpha=-1;

    vector<double> mBeta;

    double mEndTime=0;

    double mTimeStep=0;

    vector<double> mDiff;

    vector<double> mDecay;

    vector<double> mProd;

    vector<unsigned> mInitialNum;

    unsigned mTruncOrder=100;

    string mSaveDir=SAVE_DIR;

    unsigned mStartIndex=0;

    map<string, double> mAdditionalReactions;

    map<string, string> mOther;

public:
    Parameters() = default;
    explicit Parameters(map<string, docopt::value> cl_input);

    ~Parameters() = default;

    void Add(string name, string value);

    void Save(const string& path_to_file);

    void SetComments(string comment);

    void SetCommand(string command);

    void SetNumPoints(unsigned value);

    unsigned GetNumPoints();

    void SetNumThreads(unsigned value);

    unsigned GetNumThreads();

    void SetNumRuns(unsigned value);

    unsigned GetNumRuns();

    void SetNumDims(unsigned value);

    unsigned GetNumDims();

    void SetNumSpecies(unsigned value);

    unsigned GetNumSpecies();

    void SetNumMethod(string value);

    string GetNumMethod();

    void SetNumVoxels(unsigned value);

    unsigned GetNumVoxels();

    void SetDomainBounds(vector<double> values);

    vector<double> GetDomainBounds();

    void SetBC(string value);

    string GetBC();

    void SetKappa(double value);

    double GetKappa();

    void SetAlpha(double value);

    double GetAlpha();

    void SetBeta(vector<double> values);

    vector<double> GetBeta();

    void SetEndTime(double value);

    double GetEndTime();

    void SetTimeStep(double value);

    double GetTimeStep();

    void SetDiff(vector<double> values);

    vector<double> GetDiff();

    void SetDecay(vector<double> values);

    vector<double> GetDecay();

    void SetProd(vector<double> values);

    vector<double> GetProd();

    void SetInitailNum(vector<unsigned> values);

    vector<unsigned> GetInitialNum();

    void SetTruncOrder(unsigned value);

    unsigned GetTruncOrder();

    void SetSaveDir(string value);

    string GetSaveDir();

    void SetStartIndex(unsigned value);

    unsigned GetStartIndex();

    void AddAdditionalReactions(string name, double rate);

};


#endif //STOSPA_PARAMETERS_HPP
