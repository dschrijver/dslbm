#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/memory.h"

void allocate_bags(SimulationBag **sim, DistributionBag **dists, GlobalFieldBag **glob_fields, ComponentFieldBag **comp_fields, ParamBag **params, Stencil **stencil)
{
    *sim = (SimulationBag *)malloc(sizeof(SimulationBag));
    *dists = (DistributionBag *)malloc(sizeof(DistributionBag));
    *glob_fields = (GlobalFieldBag *)malloc(sizeof(GlobalFieldBag));
    *comp_fields = (ComponentFieldBag *)malloc(sizeof(ComponentFieldBag));
    *params = (ParamBag *)malloc(sizeof(ParamBag));
    *stencil = (Stencil *)malloc(sizeof(Stencil));

    (*sim)->params = (*params);
    (*sim)->dists = (*dists);
    (*sim)->glob_fields = (*glob_fields);
    (*sim)->comp_fields = (*comp_fields);
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
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;
    int malloc_size = NCOMP * (NX_proc + 2) * NY * NZ * NP * sizeof(double);

    dists->f1 = (double *)malloc(malloc_size);
    dists->f2 = (double *)malloc(malloc_size);
    dists->meq = (double *)malloc(NP * sizeof(double));
    dists->t_star = (double *)malloc(NP * sizeof(double));
    dists->m_star = (double *)malloc(NP * sizeof(double));
    dists->f_star = (double *)malloc(NP * sizeof(double));

    malloc_size = NCOMP * NY * NZ * NP * sizeof(double);
    dists->send_buffer = (double *)malloc(malloc_size);
    dists->recv_buffer = (double *)malloc(malloc_size);
}

void allocate_fields(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;
    int global_malloc_size = (NX_proc + 4) * NY * NZ * sizeof(double);
    int component_malloc_size = NCOMP * global_malloc_size;

    glob_fields->rho = (double *)malloc(global_malloc_size);
    glob_fields->pressure = (double *)malloc(global_malloc_size);
    glob_fields->u = (double *)malloc(global_malloc_size);
    glob_fields->v = (double *)malloc(global_malloc_size);
    glob_fields->w = (double *)malloc(global_malloc_size);
    glob_fields->rho_N = (double *)malloc(global_malloc_size);
    glob_fields->G_norm = (double *)malloc(global_malloc_size);
    glob_fields->Gx = (double *)malloc(global_malloc_size);
    glob_fields->Gy = (double *)malloc(global_malloc_size);
    glob_fields->Gz = (double *)malloc(global_malloc_size);
    glob_fields->nx = (double *)malloc(global_malloc_size);
    glob_fields->ny = (double *)malloc(global_malloc_size);
    glob_fields->nz = (double *)malloc(global_malloc_size);
    glob_fields->Qx = (double *)malloc(global_malloc_size);
    glob_fields->Qy = (double *)malloc(global_malloc_size);
    glob_fields->Qz = (double *)malloc(global_malloc_size);
    glob_fields->Fx = (double *)malloc(global_malloc_size);
    glob_fields->Fy = (double *)malloc(global_malloc_size);
    glob_fields->Fz = (double *)malloc(global_malloc_size);
    glob_fields->flag = (int *)malloc((NX_proc + 4) * (NY + 4) * (NZ + 4) * sizeof(int));

    int malloc_size = 2 * NY * NZ * sizeof(double);
    glob_fields->send_buffer = (double *)malloc(malloc_size);
    glob_fields->recv_buffer = (double *)malloc(malloc_size);

    comp_fields->rho_comp = (double *)malloc(component_malloc_size);
    comp_fields->u_comp = (double *)malloc(component_malloc_size);
    comp_fields->v_comp = (double *)malloc(component_malloc_size);
    comp_fields->w_comp = (double *)malloc(component_malloc_size);
}

void free_all(SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;
    ParamBag *params = sim->params;

    free(dists->f1);
    free(dists->f2);
    free(dists->meq);
    free(dists->t_star);
    free(dists->m_star);
    free(dists->f_star);
    free(dists->send_buffer);
    free(dists->recv_buffer);

    free(glob_fields->rho);
    free(glob_fields->pressure);
    free(glob_fields->u);
    free(glob_fields->v);
    free(glob_fields->w);
    free(glob_fields->rho_N);
    free(glob_fields->G_norm);
    free(glob_fields->Gx);
    free(glob_fields->Gy);
    free(glob_fields->Gz);
    free(glob_fields->nx);
    free(glob_fields->ny);
    free(glob_fields->nz);
    free(glob_fields->Qx);
    free(glob_fields->Qy);
    free(glob_fields->Qz);
    free(glob_fields->Fx);
    free(glob_fields->Fy);
    free(glob_fields->Fz);
    free(glob_fields->flag);
    free(glob_fields->send_buffer);
    free(glob_fields->recv_buffer);

    free(comp_fields->rho_comp);
    free(comp_fields->u_comp);
    free(comp_fields->v_comp);
    free(comp_fields->w_comp);

    free(stencil->cx);
    free(stencil->cy);
    free(stencil->cz);
    free(stencil->wp);
    free(stencil->p_bounceback);

    free(dists);
    free(glob_fields);
    free(comp_fields);
    free(stencil);
    free(params);
}