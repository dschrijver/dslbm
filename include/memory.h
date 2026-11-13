#ifndef MEMORY_H
#define MEMORY_H

#include "datatypes.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil, MatrixBag **matrices);
void allocate_stencil(SimulationBag *sim);
void allocate_distributions(SimulationBag *sim);
void allocate_fields(SimulationBag *sim);
void allocate_matrices(SimulationBag *sim);
void free_all(SimulationBag *sim);
int get_VmRSS(ParamBag *params) ;

#endif