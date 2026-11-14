#ifndef ACTIVE_H
#define ACTIVE_H

void compute_molecular_field(SimulationBag *sim);
void compute_nematic_stress(SimulationBag *sim);
void compute_matrix_gradient(matrix *input, matrix *output_x, matrix *output_y, matrix *output_z, SimulationBag *sim);
void compute_matrix_laplacian(matrix *input, matrix *output, SimulationBag *sim);
void compute_matrix_divergence(matrix *input, double *output_x, double *output_y, double *output_z, SimulationBag *sim);
void evolve_Q(SimulationBag *sim);

#endif