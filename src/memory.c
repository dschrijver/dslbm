#include <stdlib.h>
#include <string.h>

#include "../include/datatypes.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil, MatrixBag **matrices)
{
    *sim = (SimulationBag *)malloc(sizeof(SimulationBag));
    *dists = (DistributionBag *)malloc(sizeof(DistributionBag));
    *fields = (FieldBag *)malloc(sizeof(FieldBag));
    *params = (ParamBag *)malloc(sizeof(ParamBag));
    *stencil = (Stencil *)malloc(sizeof(Stencil));
    *matrices = (MatrixBag *)malloc(sizeof(MatrixBag));

    (*sim)->params = (*params);
    (*sim)->dists = (*dists);
    (*sim)->fields = (*fields);
    (*sim)->stencil = (*stencil);
    (*sim)->matrices = (*matrices);
}

void allocate_stencil(SimulationBag *sim)
{
    Stencil *stencil = sim->stencil;

    int NP = stencil->NP;

    stencil->cx = (int *)malloc(NP * sizeof(int));
    stencil->cy = (int *)malloc(NP * sizeof(int));
    stencil->cz = (int *)malloc(NP * sizeof(int));
    stencil->wp = (double *)malloc(NP * sizeof(double));
    stencil->p_bounceback = (int *)malloc(NP * sizeof(int));
}

void allocate_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    long malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * NP * sizeof(double);

    dists->f1_a = (double *)malloc(malloc_size);
    dists->f2_a = (double *)malloc(malloc_size);
    dists->f1_p = (double *)malloc(malloc_size);
    dists->f2_p = (double *)malloc(malloc_size);

    // Distribution buffers
    malloc_size = (long)NY_proc * NZ_proc * NP * sizeof(double);
    dists->send_buffer_x = (double *)malloc(malloc_size);
    dists->recv_buffer_x = (double *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * NZ_proc * NP * sizeof(double);
    dists->send_buffer_y = (double *)malloc(malloc_size);
    dists->recv_buffer_y = (double *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * NP * sizeof(double);
    dists->send_buffer_z = (double *)malloc(malloc_size);
    dists->recv_buffer_z = (double *)malloc(malloc_size);

}

void allocate_fields(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    long malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * sizeof(double);

    // Densities and composition
    fields->rho = (double *)malloc(malloc_size);
    fields->phi = (double *)malloc(malloc_size);
    fields->rho_a = (double *)malloc(malloc_size);
    fields->rho_p = (double *)malloc(malloc_size);

    // Velocities
    fields->u = (double *)malloc(malloc_size);
    fields->v = (double *)malloc(malloc_size);
    fields->w = (double *)malloc(malloc_size);
    fields->u_a = (double *)malloc(malloc_size);
    fields->v_a = (double *)malloc(malloc_size);
    fields->w_a = (double *)malloc(malloc_size);
    fields->u_p = (double *)malloc(malloc_size);
    fields->v_p = (double *)malloc(malloc_size);
    fields->w_p = (double *)malloc(malloc_size);

    // Forces
    fields->Fx_a = (double *)malloc(malloc_size);
    fields->Fy_a = (double *)malloc(malloc_size);
    fields->Fz_a = (double *)malloc(malloc_size);
    fields->Fx_p = (double *)malloc(malloc_size);
    fields->Fy_p = (double *)malloc(malloc_size);
    fields->Fz_p = (double *)malloc(malloc_size);

    // Chemical potentials
    fields->mu_a = (double *)malloc(malloc_size);
    fields->mu_p = (double *)malloc(malloc_size);

    // Helper fields
    fields->gradient_input = (double *)malloc(malloc_size);
    fields->gradient_output_x = (double *)malloc(malloc_size);
    fields->gradient_output_y = (double *)malloc(malloc_size);
    fields->gradient_output_z = (double *)malloc(malloc_size);
    fields->laplacian_input = (double *)malloc(malloc_size);
    fields->laplacian_output = (double *)malloc(malloc_size);

    // Scalar buffers
    malloc_size = (long)NY_proc * NZ_proc * sizeof(double);
    fields->send_buffer_x = (double *)malloc(malloc_size);
    fields->recv_buffer_x = (double *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * NZ_proc * sizeof(double);
    fields->send_buffer_y = (double *)malloc(malloc_size);
    fields->recv_buffer_y = (double *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * sizeof(double);
    fields->send_buffer_z = (double *)malloc(malloc_size);
    fields->recv_buffer_z = (double *)malloc(malloc_size);
}

void allocate_matrices(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    long malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * sizeof(matrix);

    matrices->Q = (matrix *)malloc(malloc_size);
    matrices->H = (matrix *)malloc(malloc_size);
    matrices->sigma_nem = (matrix *)malloc(malloc_size);
    matrices->matrix_gradient_output_x = (matrix *)malloc(malloc_size);
    matrices->matrix_gradient_output_y = (matrix *)malloc(malloc_size);
    matrices->matrix_gradient_output_z = (matrix *)malloc(malloc_size);
    matrices->matrix_divergence_output_x = (double *)malloc((long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * sizeof(double));
    matrices->matrix_divergence_output_y = (double *)malloc((long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * sizeof(double));
    matrices->matrix_divergence_output_z = (double *)malloc((long)(NX_proc + 2) * (NY_proc + 2) * (NZ_proc + 2) * sizeof(double));

    // Matrix buffers
    malloc_size = (long)NY_proc * NZ_proc * sizeof(matrix);
    matrices->send_buffer_x = (matrix *)malloc(malloc_size);
    matrices->recv_buffer_x = (matrix *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * NZ_proc * sizeof(matrix);
    matrices->send_buffer_y = (matrix *)malloc(malloc_size);
    matrices->recv_buffer_y = (matrix *)malloc(malloc_size);

    malloc_size = (long)(NX_proc + 2) * (NY_proc + 2) * sizeof(matrix);
    matrices->send_buffer_z = (matrix *)malloc(malloc_size);
    matrices->recv_buffer_z = (matrix *)malloc(malloc_size);
}

void free_all(SimulationBag *sim)
{
    UNPACK_BAGS

    H5Pclose(params->fapl_id);
    H5Sclose(params->scalar_space);
    H5Sclose(params->filespace);
    H5Sclose(params->memspace);
    H5Pclose(params->dcpl_id);
    H5Pclose(params->dxpl_id);

    free(stencil->cx);
    free(stencil->cy);
    free(stencil->cz);
    free(stencil->wp);
    free(stencil->p_bounceback);

    free(dists->f1_a);
    free(dists->f1_p);
    free(dists->f2_a);
    free(dists->f2_p);

    // Distribution buffers
    free(dists->send_buffer_x);
    free(dists->recv_buffer_x);
    free(dists->send_buffer_y);
    free(dists->recv_buffer_y);
    free(dists->send_buffer_z);
    free(dists->recv_buffer_z);

    // Densities and composition
    free(fields->rho);
    free(fields->phi);
    free(fields->rho_a);
    free(fields->rho_p);

    // Velocities
    free(fields->u);
    free(fields->v);
    free(fields->w);

    // Forces
    free(fields->Fx_a);
    free(fields->Fy_a);
    free(fields->Fz_a);
    free(fields->Fx_p);
    free(fields->Fy_p);
    free(fields->Fz_p);

    // Chemical potentials
    free(fields->mu_a);
    free(fields->mu_p);

    // Helper fields
    free(fields->gradient_input);
    free(fields->gradient_output_x);
    free(fields->gradient_output_y);
    free(fields->gradient_output_z);
    free(fields->laplacian_input);
    free(fields->laplacian_output);

    // Scalar buffers
    free(fields->send_buffer_x);
    free(fields->recv_buffer_x);
    free(fields->send_buffer_y);
    free(fields->recv_buffer_y);
    free(fields->send_buffer_z);
    free(fields->recv_buffer_z);

    free(matrices->Q);
    free(matrices->H);
    free(matrices->sigma_nem);
    free(matrices->matrix_gradient_output_x);
    free(matrices->matrix_gradient_output_y);
    free(matrices->matrix_gradient_output_z);
    free(matrices->matrix_divergence_output_x);
    free(matrices->matrix_divergence_output_y);
    free(matrices->matrix_divergence_output_z);

    // Matrix buffers
    free(matrices->send_buffer_x);
    free(matrices->recv_buffer_x);
    free(matrices->send_buffer_y);
    free(matrices->recv_buffer_y);
    free(matrices->send_buffer_z);
    free(matrices->recv_buffer_z);

    free(dists);
    free(fields);
    free(stencil);
    free(matrices);
    free(params);
}

int get_VmRSS(ParamBag *params) 
{
    FILE *f = fopen("/proc/self/status", "r");
    char line[128];

    int VmRSS, total_VmRSS = 0;
    while (fgets(line, 128, f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            sscanf(line, "VmRSS: %d", &VmRSS);      
            break;
        }
    }
    fclose(f);

    MPI_Reduce(&VmRSS, &total_VmRSS, 1, MPI_INT, MPI_SUM, 0, params->comm_cart);

    return total_VmRSS;
}