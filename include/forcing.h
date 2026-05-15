#ifndef FORCING_H
#define FORCING_H

#include "datatypes.h"

void evaluate_forces(SimulationBag *sim);
void evaluate_force(int i, int j, int k, SimulationBag *sim);
void evaluate_surface_force(int i, int j, int k, SimulationBag *sim);
double extrapolate_wall_n(int i, int j, int k, int alpha, SimulationBag *sim);

#endif