#ifndef FORCING_H
#define FORCING_H

#include "datatypes.h"

void evaluate_forces(SimulationBag *sim);
void evaluate_force(int i, int j, int k, SimulationBag *sim);

#endif