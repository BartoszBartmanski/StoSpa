//
// Created by bartosz on 05/05/17.
//

#ifndef STOSPA_SIMULATION_2D_HPP
#define STOSPA_SIMULATION_2D_HPP

#include "AbstractSimulation.hpp"
#include "JumpRates.hpp"

using namespace std;

class Simulation_2d : public AbstractSimulation
{
protected:

    /** Ratio of horizontal spacing to the vertical spacing between voxels. */
    double mRatio=1.0;

    /** Free parameter for FDM method derivation of jump coefficients. */
    double mAlpha=0.0;

    /** Free parameter in the x-direction for FET method derivation of jump coefficients. */
    double mBetaX=0.5;

    /** Free parameter in the y-direction for FET method derivation of jump coefficients. */
    double mBetaY=0.5;

public:

    /**
     * Constructor.
     * @param num_species - number of species
     * @param num_method - numerical method from which the jump coefficients are derived
     * @param num_voxels - number of voxels present in the simulation
     * @param domain_bounds - bounds of the domain
     * @param boundary_condition - reflective or periodic boundary condition
     */
    Simulation_2d(unsigned num_runs, unsigned num_species, string num_method, unsigned num_voxels,
                  vector<double> domain_bounds, string boundary_condition, double ratio);
    explicit Simulation_2d(Parameters params);

    Simulation_2d(const Simulation_2d&) = delete; //move only type
    Simulation_2d& operator=(const Simulation_2d&) = delete; //move only type
    Simulation_2d(Simulation_2d&&) = default;
    Simulation_2d& operator=(Simulation_2d&&) = default;
    ~Simulation_2d() override = default;

    /**
     * Returns the aspect ratio of voxels (horizontal length divided by the vertical voxel length).
     * @return a double
     */
    double GetVoxelRatio();

    /**
     * Set the value of alpha (parameter in FDM derivation)
     * @param value
     */
    void SetAlpha(double value);

    /**
     * Returns the value of alpha (used in fdm derived jump rates).
     * @return a double
     */
    double GetAlpha();

    /**
     * Set the values of beta_x and beta_y (parameters in FET derivation)
     * @param values
     */
    void SetBeta(vector<double> values);

    /**
     * Returns the values of beta_x and beta_y (parameters in FET derivation)
     * @return a double
     */
    vector<double> GetBeta();

    /**
     * Returns pointer to JumpRate object that can give jump rates in specified directions.
     * @return unique_ptr<JumpRate>
     */
    unique_ptr<JumpRate> GetJumpRates();

    /**
     * Method that populates the mLambdas vector (vector of propensities).
     * @param diff
     * @param species
     */
    void SetDiffusionRate(double diff, unsigned species) override;

    /**
     * Method to place the specified number of molecules of the specified species at the specified voxel index
     * @param voxel_index - index of the voxel where the molecules will be placed
     * @param num_molecules - number of molecules that will be placed
     * @param species - index of the species of the molecules
     */
    void SetInitialNumMolecules(vector<unsigned> voxel_index, unsigned num_molecules, unsigned species) override;

    /**
     * Another method of initialising the number of molecules in the voxels
     * @param initial_state - state of mVoxels initially
     * @param species - species in the give state
     */
    void SetInitialState(vector<vector<unsigned> > initial_state, unsigned species) override;
    void SetInitialState(vector<vector<int> > initial_state, unsigned species) override;
};

#endif //STOSPA_SIMULATION_2D_HPP
