#include <stdlib.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../include/memory.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, FieldBag **fields, ParamBag **params, Stencil **stencil)
{
    *sim = (SimulationBag *)malloc(sizeof(SimulationBag));
    *dists = (DistributionBag *)malloc(sizeof(DistributionBag));
    *fields = (FieldBag *)malloc(sizeof(FieldBag));
    *params = (ParamBag *)malloc(sizeof(ParamBag));
    *stencil = (Stencil *)malloc(sizeof(Stencil));

    (*sim)->params = (*params);
    (*sim)->dists = (*dists);
    (*sim)->fields = (*fields);
    (*sim)->stencil = (*stencil);
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

    dists->f1 = (double *)malloc(malloc_size);
    dists->f2 = (double *)malloc(malloc_size);

    dists->g1 = (double *)malloc(malloc_size);
    dists->g2 = (double *)malloc(malloc_size);

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

    fields->rho = (double *)malloc(malloc_size);
    fields->phi = (double *)malloc(malloc_size);
    fields->pressure = (double *)malloc(malloc_size);
    fields->u = (double *)malloc(malloc_size);
    fields->v = (double *)malloc(malloc_size);
    fields->w = (double *)malloc(malloc_size);
    fields->Fx = (double *)malloc(malloc_size);
    fields->Fy = (double *)malloc(malloc_size);
    fields->Fz = (double *)malloc(malloc_size);

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

    free(dists->f1);
    free(dists->f2);
    free(dists->g1);
    free(dists->g2);
    free(dists->send_buffer_x);
    free(dists->recv_buffer_x);
    free(dists->send_buffer_y);
    free(dists->recv_buffer_y);
    free(dists->send_buffer_z);
    free(dists->recv_buffer_z);

    free(fields->rho);
    free(fields->phi);
    free(fields->pressure);
    free(fields->u);
    free(fields->v);
    free(fields->w);
    free(fields->Fx);
    free(fields->Fy);
    free(fields->Fz);
    free(fields->send_buffer_x);
    free(fields->recv_buffer_x);
    free(fields->send_buffer_y);
    free(fields->recv_buffer_y);
    free(fields->send_buffer_z);
    free(fields->recv_buffer_z);

    free(dists);
    free(fields);
    free(stencil);
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