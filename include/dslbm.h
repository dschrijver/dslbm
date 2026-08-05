#ifndef DSLBM_H
#define DSLBM_H

#include "datatypes.h"
#include "stencil.h"

// collide.c
void collide(SimulationBag *sim);
void evaluate_color_gradients(SimulationBag *sim);
void compute_Q_corrections(SimulationBag *sim);

// communicate.c
void communicate_field(double *field, SimulationBag *sim);
void communicate_dist(double *dist, SimulationBag *sim);

// fields.c
void extract_moments(SimulationBag *sim);
void update_final_velocity(SimulationBag *sim);
void evaluate_pressure(SimulationBag *sim);

// forcing.c
void evaluate_forces(SimulationBag *sim);
void evaluate_force(int i, int j, int k, SimulationBag *sim);

// initialize.c
void initialize_MPI(ParamBag *params);
void initialize_HDF5(ParamBag *params);
void initialize_fields(SimulationBag *sim);
void initialize_flags(SimulationBag *sim);
void initialize_distributions(SimulationBag *sim);

// memory.c
void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil);
void allocate_stencil(SimulationBag *sim);
void allocate_distributions(SimulationBag *sim);
void allocate_fields(SimulationBag *sim);
void free_all(SimulationBag *sim);

// output.c
void output_data(SimulationBag *sim);

// stream.c
void stream_distributions(SimulationBag *sim);

// wetnode.c
void wetnode_macroscopic_fields(SimulationBag *sim);
void wetnode_distributions(SimulationBag *sim);

#endif