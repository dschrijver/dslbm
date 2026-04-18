#ifndef FIELDS_H
#define FIELDS_H

#include "datatypes.h"

void extract_moments(SimulationBag *sim);
void update_final_velocity(SimulationBag *sim);
double evaluate_mass(int n, SimulationBag *sim);

#endif