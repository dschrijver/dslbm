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

    stencil->M = (double *)malloc(NP * NP * sizeof(double));
    stencil->M_inv = (double *)malloc(NP * NP * sizeof(double));
    stencil->omega[RED] = (double *)malloc(NP * sizeof(double));
    stencil->omega[BLUE] = (double *)malloc(NP * sizeof(double));

    stencil->B = (double *)malloc(NP * sizeof(double));
    stencil->phi_eq[RED] = (double *)malloc(NP * sizeof(double));
    stencil->phi_eq[BLUE] = (double *)malloc(NP * sizeof(double));
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

    comp_fields->rho_comp = (double *)malloc(component_malloc_size);
    comp_fields->Fx = (double *)malloc(component_malloc_size);
    comp_fields->Fy = (double *)malloc(component_malloc_size);
    comp_fields->Fz = (double *)malloc(component_malloc_size);
}