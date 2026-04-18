#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "datatypes.h"

void initialize_MPI(ParamBag *params);
void initialize_HDF5(ParamBag *params);
void initialize_fields(SimulationBag *sim);
void initialize_flags(SimulationBag *sim);
void initialize_distributions(SimulationBag *sim);
double physx(int i);
double physy(int j);
double physz(int k);

#endif