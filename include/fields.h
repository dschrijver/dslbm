#ifndef FIELDS_H
#define FIELDS_H

#include "datatypes.h"

void compute_gradient(SimulationBag *sim);
void compute_laplacian(SimulationBag *sim);
void compute_chemical_potentials(SimulationBag *sim);
void extract_moments(SimulationBag *sim);
void compute_final_velocity(SimulationBag *sim);

#endif