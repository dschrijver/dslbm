# DSLBM
DSLBM is an MPI parallelized binary fluid code based on the Lattice Boltzmann Method. It makes use of the Color-Gradient Method for phase separation and interface tension. 

## Local installation on Ubuntu
Run the following commands to install make, gcc, openmpi, hdf5 and zlib, 
```
sudo apt install binutils openmpi-bin libopenmpi-dev zlib1g zlib1g-dev make g++
make install_hdf5
```

## Cluster installation
Load the modules corresponding to gcc, openmpi and hdf5. 

## Program settings
Parameters can be tuned in the file [params.h](params.h) before compilation. Boundary conditions and initial conditions, along with their corresponding parameters, can be tuned in [definitions.h](definitions.h). There should be some notes about the program settings:

- The number of nodes in the $x$-direction, NX, should be at least twice as big as the number of MPI processes. 
- The improved Non-Equilibrium Bounce-Back scheme in this code performs more accurately for Static Contact-Angle simulations. 
- Corners between NEBB boundaries are not implemented.  
- Custom initial conditions can be added to [definitions.h](definitions.h) and then to [initialize.c](src/initialize.c) using a newly defined flag. 
- No logic checks are performed when multiple boundary condition methods are defined for the same boundary. 
- You should not enable both the Color-Gradient Method and Shan-Chen Method in [definitions.h](definitions.h).
- The flags ```THETA_C_[BOUNDARY]``` in [definitions.h](definitions.h) are used to set the contact angles when the Color-Gradient Method is enabled, and the flags ```XI_[BOUNDARY]``` are used to tune the contact angles when the Shan-Chen Method is enabled. 

## Local compilation and running
Before compilation, copy [make.def.local](makedefs/make.def.local) to the root directory. Run the following command to compile the code,
```
make
```
The program can then be ran on $n$ processors using the following command,
```
make run_local n=n
```

## Cluster compilation and running
Before compiliation, create a make.def file corresponding to the specific cluster architecture. An example can be found in [make.def.example_cluster](makedefs/make.def.example_cluster). Copy this new file to the root directory and compile the code with the following command,
```
make
```
The ```dslbm``` executable will be generated, after which a custom jobscript can be submitted. 
