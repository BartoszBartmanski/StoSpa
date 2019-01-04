//
// Created by bartmanski on 07/04/17.
//

#ifndef STOSPA_SIMULATION_1D_HPP
#define STOSPA_SIMULATION_1D_HPP

#include "AbstractSimulation.hpp"

using namespace std;

class Simulation1d : public AbstractSimulation
{
public:

    /**
     * Default constructor
     * @param num_species - number of species
     * @param num_voxels - number of voxels present in the simulation
     * @param domain_bounds - bounds of the domain
     * @param boundary_condition - reflective or periodic boundary condition
     */
    Simulation1d(unsigned num_runs, unsigned num_species, unsigned num_voxels, vector<double> domain_bounds,
                  string boundary_condition);
    explicit Simulation1d(Parameters params);

    Simulation1d(const Simulation1d&) = delete; //move only type
    Simulation1d& operator=(const Simulation1d&) = delete; //move only type
    Simulation1d(Simulation1d&&) = default;
    Simulation1d& operator=(Simulation1d&&) = default;
    ~Simulation1d() override = default;

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


#endif //STOSPA_SIMULATION_1D_HPP
