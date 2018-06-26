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


using namespace std;

class Parameters
{
private:
    string mComments;

    string mCommand;

    unsigned mNumThreads;

    unsigned mNumRuns;

    unsigned mNumDims;

    unsigned mNumSpecies;

    string mNumMethod;

    unsigned mNumVoxels;

    vector<double> mDomainBounds;

    double m_h;

    string mBC;

    double mKappa;

    vector<double> mAlpha;

    vector<double> mBeta;

    double mEndTime;

    double mTimeStep;

    unsigned mNumPoints;

    vector<double> mDiff;

    vector<double> mDecay;

    vector<double> mProd;

    vector<unsigned> mInitialNumMolecules;

    unsigned mTruncOrder;

    map<string, double> mAdditionalReactions;

    map<string, string> mOther;


public:

    Parameters();

    ~Parameters() = default;

    void Add(const string& name, const string& value);

    void Save(const string& path_to_file);

    void SetComments(string comment);

    void SetCommand(string command);

    void SetNumThreads(unsigned num_threads);

    void SetNumRuns(unsigned num_runs);

    void SetNumDims(unsigned num_dims);

    void SetNumSpecies(unsigned num_species);

    void SetNumMethod(string num_method);

    void SetNumVoxels(unsigned num_voxels);

    void SetDomainBounds(vector<double> domain_bounds);

    void SetBC(string bc);

    void SetKappa(double kappa);

    void SetAlpha(double alpha);

    void SetBeta(vector<double> beta);

    void SetEndTime(double end_time);

    void SetTimeStep(double time_step);

    void SetNumPoints(unsigned num_points);

    void SetDiff(vector<double> diff);
    void SetDiff(double diff);

    void SetDecay(vector<double> decay);
    void SetDecay(double decay);

    void SetProd(vector<double> prod);
    void SetProd(double prod);

    void SetInitialNumMolecules(unsigned num_molecules);
    void SetInitialNumMolecules(vector<unsigned> num_molecules);

    void SetTruncOrder(unsigned trunc_order);

    void AddAdditionalReactions(string name, double rate);

};


#endif //STOSPA_PARAMETERS_HPP
