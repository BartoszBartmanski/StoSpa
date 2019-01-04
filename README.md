# StoSpa - stochastic simulation of spatially extended systems

A C++ package for running stochastic simulations to generate sample paths 
for reaction-diffusion master equation (RDME). Reactions happen within a compartment,
and the collection of compartments makes up the domain of the system. Diffusion of 
molecules is modelled as jumps between adjacent voxels.

This package is able to simulate a system with time-dependant propensities
using the Extrande algorithm.

## Installation

First, clone this repository to your computer by running.
```
git clone https://github.com/BartoszBartmanski/StoSpa.git
```

To install we have to make a directory for the executables. Go to the directory
where StoSpa was cloned to. Then execute the following commands
```
cd StoSpa
mkdir build
cd build
```
which will create a directory for the executables.

Now we use cmake to build a makefile.
```
cmake ../
```
Next, we compile whichever simulation we want or we can compile all of the executables
```
make all
```
which will create an executable for each of the files present in Simulations directory.


## Example

Below we go through a simulation example.cpp (present in Simulations directory).

First, we add the header that enables us to run simulations.
```
#include "Simulator.hpp"
```
Next, at the start of the main function, we define the variables needed 
for the simulation.
```
string sim_name = "example.dat";
unsigned num_runs = 1;
unsigned num_species = 1;
unsigned num_voxels = 10;
vector<double> domain_bounds = {0, 1.0};
string bc = "reflective";
double diff = 0.1;
double end_time = 5.0;
double time_step = 0.01;
```
To run a simulation, we need to initialise a Simulation object (either 
Simulation1d or Simulation2d).
```
Simulation1d sim(num_runs, num_species, num_voxels, domain_bounds, bc);
```
We can set a number of molecules in the first compartment as follows
```
sim.SetInitialNumMolecules({0}, 1000, 0);
```
where we place a thousand molecules (second argument) of first species (third argument) 
into the first voxel (first argument).
Lastly, we run the simulation using the member function Run.
```
sim.Run(sim_name, end_time, time_step);
```

After we run this simulation we get the following output

![](http://users.ox.ac.uk/~shil4444/misc/sim.gif)


