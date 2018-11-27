//
// Created by bartosz on 24/06/17.
//

#ifndef STOSPA_ABSTRACTSIMULATION_HPP
#define STOSPA_ABSTRACTSIMULATION_HPP

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <random>
#include <cassert>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <map>
#include <memory>
#include <queue>
#include "AbstractReaction.hpp"
#include "DiffusionReflective.hpp"
#include "DiffusionPeriodic.hpp"
#include "VectorFunctions.hpp"
#include "Grid.hpp"
#include "Utilities.hpp"
#include "None.hpp"

using namespace std;

class AbstractSimulation
{
protected:

    /** Helper constants */
    double inf;

    /** Switch for setting times until next reactions vector. */
    bool mTimesSet;

    /** Number of runs of this simulation */
    unsigned mNumRuns;

    /** Number of species in the simulation. */
    unsigned mNumSpecies;

    /** Numerical method from which jump coefficients are derived. */
    string mNumMethod;

    /** Number of voxels in the simulation. */
    vector<unsigned> mNumVoxels;

    unsigned mTotalNumVoxels;

    /** Vector specifying the dimensions of the simulation domain. */
    vector<double> mDomainBounds;

    /** Map that contains additional reactions. */
    vector<shared_ptr<AbstractReaction>> mReactions;

    /** Boundary condition. */
    string mBC = "reflective";

    /** Current time of the specific run of the simulation. */
    vector<double> mCurrentTime;

    /** Current time for all the runs. */
    double mTime;

    /** Number of jumps at the at time in the simulation. */
    unsigned mNumJumps;

    /** Voxel width in one dimension and voxel height in two dimensions. */
    double m_h;

    /** Voxel size (not necessarily the same as the voxel spacing) */
    double mVoxelSize;

    /** Total number of molecules of each species. Each entry corresponds to a species in the simulation. */
    vector<unsigned> mTotalNumMolecules;

    /** A vector of grids, each of which consists of an array for voxels and an array for time increments. */
    vector<Grid> mGrids;

    /** Seed used for generating a random number. */
    unsigned mSeed;

    /** For generating random numbers */
    mt19937 mGen;

    /** Uniform distribution. */
    uniform_real_distribution<double> mUniform;

    /** vector of all propensities. */
    vector<double> mPropensities;

    /** Vector of the diffusion coefficients. */
    vector<double> mDiffusionCoefficients;

    /** Calculates the total propensity for a specified voxel. */
    inline double GetTotalPropensity(const unsigned& run, const int& voxel_index);

    inline unsigned NextReaction(const unsigned& run, const int& voxel_index);

    inline double Exponential(double propensity);

public:

    /** Default constructor. */
    AbstractSimulation();


    /** Default destructor. */
    virtual ~AbstractSimulation()= default;

    void SetSeed(unsigned number);

    unsigned GetSeed();

    inline void UpdateTime(unsigned run, unsigned voxel_index);

    /**
     * SSA loop. A single molecule jump or a single reaction.
     */
    void SSA_loop(const unsigned& run);

    /**
     * Method that will invoke SSA loop
     * @param iterator
     * @param time_step
     */
    void Advance(const double& time_step, const unsigned& iterator=1);

    /** Method to occupy the grid with the time increments at the beginning of the simulation. */
    void SetupTimeIncrements();

    /**
     * Method that populates the mLambdas vector (vector of propensities).
     * @param diffusion
     * @param decay
     * @param production
     * @param species
     */
    virtual void SetDiffusionRate(double diffusion_coefficient, unsigned species)=0;

    /**
     * Method to place the specified number of molecules of the specified species at the specified voxel index
     * @param voxel_index - index of the voxel where the molecules will be placed
     * @param num_molecules - number of molecules that will be placed
     * @param species - index of the species of the molecules
     */
    virtual void SetInitialNumMolecules(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species)=0;

    /**
     * Another method of initialising the number of molecules in the voxels
     * @param initial_state - state of mVoxels initially
     * @param species - species in the give state
     */
    virtual void SetInitialState(vector<vector<unsigned> > initial_state, unsigned species)=0;
    virtual void SetInitialState(vector<vector<int> > initial_state, unsigned species)=0;

    /**
     * Method that adds different types of reactions from diffusion, decay and production.
     * @param reaction_name
     * @param rate_constant
     */
    void AddReaction(shared_ptr<AbstractReaction> reaction);

    /**
     * Returns a vector of pointers to additional reactions.
     * @return mAdditionalReactions
     */
    vector<shared_ptr<AbstractReaction>> GetReactions();

    /**
     * Returns the current state of the mVoxels
     * @param species - index of the species
     * @return mVoxels associated with the specified species
     */
    vector<unsigned> GetVoxels(unsigned species=0, unsigned run=0);

    /**
     * Returns the concentration of the specified species.
     * @param species - index of the species
     * @return concentration[species]
     */
    vector<double> GetConcentration(unsigned species=0);

    /**
     * Returns the number of molecules averaged over all the runs.
     * @param species - index of the species
     * @return average_molecules
     */
    vector<double> GetAverageNumMolecules(unsigned species = 0);

    /**
     * Returns the number of molecules of specified species and in a specified run
     * @param species - index of the species
     * @param run - index of the run
     * @return vector<unsigned> for number of molecules.
     */
    vector<unsigned> GetNumMolecules(unsigned species = 0, unsigned run = 0);

    /**
     * Returns current time.
     * @return mTime
     */
    double GetCurrentTime();

    /**
     * Returns spacing in the simulation (either the voxel size in 1d or the vertical voxel size in 2d).
     * @return mSpacing
     */
    double GetSpacing();

    /**
     * Returns voxel size
     * @return mVoxelSize
     */
    double GetVoxelSize();

    /**
     * Returns the total number of molecules in the simulation of the specified species.
     * @param species - index of the species
     * @return mTotalNumMolecules[species]
     */
    unsigned GetInitialTotalMolecules(unsigned species = 0);

    /**
     * Returns the current total number of molecules for given species index and given run.
     * @param species - index of the species
     * @param run - inde
     * @return total number of molecules of specified species in the current state of the simulation.
     */
    unsigned GetTotalMolecules(unsigned species = 0, unsigned run = 0);

    /**
     * Returns the total number of reactions that have taken place.
     * @return mNumJumps
     */
    unsigned GetNumJumps();

    /**
     * Returns the number of voxels along each direction.
     * @return mNumVoxels
     */
    vector<unsigned> GetNumVoxels();

    /**
     * Returns the number of runs within the simulation.
     * @return mNumRuns
     */
    unsigned GetNumRuns();

    /**
     * Returns the number of species.
     * @return mNumSpecies
     */
    unsigned GetNumSpecies();

    /**
     * Returns the numerical method used to derive the values of lambda associated with diffusion.
     * @return mNumMethod
     */
    string GetNumMethod();

    /**
     * Returns the domain bounds.
     * @return mDomainBounds
     */
    vector<double> GetDomainBounds();

    /**
     * Returns the boundary condition used.
     * @return mBoundaryCondtion
     */
    string GetBoundaryCondition();

    /**
     * Returns the diffusion coefficient for a given species.
     * @param species - index of the species
     * @return mDiffusionCoefficient[species]
     */
    double GetDiffusionCoefficient(unsigned species=0);

    /**
     * Returns the error for the specified species
     * @param species - index of the species
     * @return error
     */
    double GetError(const vector<double>& analytic, unsigned species=0);

    /**
     * Returns the relative error
     * @param species - index of the species
     * @return relative_error
     */
    double GetRelativeError(const vector<double>& analytic, unsigned species=0);

};

#endif //STOSPA_ABSTRACTSIMULATION_HPP
