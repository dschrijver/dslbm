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
    fields->Fx_nem = (double *)malloc(malloc_size);
    fields->Fy_nem = (double *)malloc(malloc_size);
    fields->Fz_nem = (double *)malloc(malloc_size);
    fields->dx_rho = (double *)malloc(malloc_size);
    fields->dy_rho = (double *)malloc(malloc_size);
    fields->dz_rho = (double *)malloc(malloc_size);

    // Chemical potentials
    fields->mu_a = (double *)malloc(malloc_size);
    fields->mu_p = (double *)malloc(malloc_size);

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
    matrices->Q_new = (matrix *)malloc(malloc_size);
    matrices->H = (matrix *)malloc(malloc_size);
    matrices->sigma_nem = (matrix *)malloc(malloc_size);
    matrices->nabla_u = (matrix *)malloc(malloc_size);
    matrices->E = (matrix *)malloc(malloc_size);
    matrices->Omega = (matrix *)malloc(malloc_size);
    matrices->dx_Q = (matrix *)malloc(malloc_size);
    matrices->dy_Q = (matrix *)malloc(malloc_size);
    matrices->dz_Q = (matrix *)malloc(malloc_size);

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