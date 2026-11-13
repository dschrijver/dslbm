#ifndef ACTIVE_H
#define ACTIVE_H

void compute_molecular_field(SimulationBag *sim);
void compute_nematic_stress(SimulationBag *sim);
void compute_matrix_gradient(matrix *input, SimulationBag *sim);
void compute_matrix_divergence(matrix *input, SimulationBag *sim);

#endif