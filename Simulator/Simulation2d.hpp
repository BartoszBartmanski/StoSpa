//
// Created by bartosz on 05/05/17.
//

#ifndef STOSPA_SIMULATION_2D_HPP
#define STOSPA_SIMULATION_2D_HPP

#include "AbstractSimulation.hpp"
#include "JumpRates.hpp"

using namespace std;

class Simulation2d : public AbstractSimulation
{
protected:

    /** Ratio of horizontal spacing to the vertical spacing between voxels. */
    double mRatio=1.0;

public:

    /**
     * Constructor.
     * @param num_species - number of species
     * @param num_voxels - number of voxels present in the simulation
     * @param domain_bounds - bounds of the domain
     * @param boundary_condition - reflective or periodic boundary condition
     */
    Simulation2d(unsigned num_runs, unsigned num_species, unsigned num_voxels, vector<double> domain_bounds,
                  string boundary_condition, double ratio);
    explicit Simulation2d(Parameters params);

    Simulation2d(const Simulation2d&) = delete; //move only type
    Simulation2d& operator=(const Simulation2d&) = delete; //move only type
    Simulation2d(Simulation2d&&) = default;
    Simulation2d& operator=(Simulation2d&&) = default;
    ~Simulation2d() override = default;

    /**
     * Returns the aspect ratio of voxels (horizontal length divided by the vertical voxel length).
     * @return a double
     */
    double GetVoxelRatio();

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
