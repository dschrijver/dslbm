#ifndef FIELDS_H
#define FIELDS_H

#include "datatypes.h"

void compute_gradient(double *input, double *output_x, double *output_y, double *output_z, SimulationBag *sim);
void compute_laplacian(double *input, double *output, SimulationBag *sim);
void compute_chemical_potentials(SimulationBag *sim);
void extract_moments(SimulationBag *sim);
void compute_final_velocity(SimulationBag *sim);

#endif