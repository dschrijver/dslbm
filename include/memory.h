#ifndef MEMORY_H
#define MEMORY_H

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil);
void allocate_stencil(SimulationBag *sim);
void allocate_distributions(SimulationBag *sim);
void allocate_fields(SimulationBag *sim);
void free_all(SimulationBag *sim);

#endif