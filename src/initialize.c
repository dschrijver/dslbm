#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/initialize.h"

void initialize_MPI(ParamBag *params)
{
    MPI_Comm_size(MPI_COMM_WORLD, &params->number_of_processes);
    int dims[3] = {params->number_of_processes, 1, 1};
    int xperiodic = 0;
#ifdef XPERIODIC
    xperiodic = 1;
#endif
    int periods[3] = {xperiodic, 0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &params->comm_xslices);
    MPI_Comm_rank(params->comm_xslices, &params->process_rank);
    MPI_Cart_get(params->comm_xslices, 3, dims, periods, params->process_coords);
    MPI_Cart_shift(params->comm_xslices, 0, 1, &params->process_neighbors[0], &params->process_neighbors[1]);
    params->i_start = (float)(params->process_coords[0]) / (float)params->number_of_processes * params->NX;
    params->i_end = (float)(params->process_coords[0] + 1) / (float)params->number_of_processes * params->NX;
    params->NX_proc = params->i_end - params->i_start;

    if (params->NX < 2 * params->number_of_processes)
    {
        if (params->process_rank == 0)
        {
            printf("NX (%d) must be at least twice as big as number of processes (%d)!\n", params->NX, params->number_of_processes);
            fflush(stdout);
        }

        MPI_Finalize();
        exit(1);
    }
}

void initialize_fields(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double rho_0_RED = params->rho_0_RED;
    double rho_0_BLUE = params->rho_0_BLUE;

    double zeta = stencil->zeta;
    double alpha_RED = params->alpha_RED;
    double alpha_BLUE = params->alpha_BLUE;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *rho_comp = comp_fields->rho_comp;

#ifdef INI_FLOATING_DROPLET
    double x, y, z, r2;
    double R_droplet = 15.0;
    FOR_DOMAIN
    {
        x = (double)i + 0.5 - 0.5*(double)NX;
        y = (double)j + 0.5 - 0.5*(double)NY;
        z = (double)k + 0.5 - 0.5*(double)NZ;
        r2 = x * x + y * y + z * z;

        u[INDEX_GLOB(i, j, k)] = 0.0;
        v[INDEX_GLOB(i, j, k)] = 0.0;
        w[INDEX_GLOB(i, j, k)] = 0.0;

        if (r2 < R_droplet * R_droplet)
        {
            rho_comp[INDEX(i, j, k, RED)] = rho_0_RED;
            rho_comp[INDEX(i, j, k, BLUE)] = 0.0;
        }
        else
        {
            rho_comp[INDEX(i, j, k, RED)] = 0.0;
            rho_comp[INDEX(i, j, k, BLUE)] = rho_0_BLUE;
        }
    }
#endif

    FOR_DOMAIN
    {
        rho[INDEX_GLOB(i, j, k)] = rho_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)];
        pressure[INDEX_GLOB(i, j, k)] = rho_comp[INDEX(i, j, k, RED)] * zeta * (1.0 - alpha_RED) + rho_comp[INDEX(i, j, k, BLUE)] * zeta * (1.0 - alpha_BLUE);
    }
}

void initialize_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_i, rho_RED_i, rho_BLUE_i, uhat, vhat, what, uc, u2, eq;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;
    double **phi_eq = stencil->phi_eq;

    double *rho = glob_fields->rho;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *f1 = dists->f1;

    FOR_DOMAIN
    {
        rho_i = rho[INDEX_GLOB(i, j, k)];
        rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
        rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];
        uhat = u[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fx[INDEX(i, j, k, RED)] + Fx[INDEX(i, j, k, BLUE)]);
        vhat = v[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fy[INDEX(i, j, k, RED)] + Fy[INDEX(i, j, k, BLUE)]);
        what = w[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fz[INDEX(i, j, k, RED)] + Fz[INDEX(i, j, k, BLUE)]);
        u2 = uhat * uhat + vhat * vhat + what * what;

        for (int p = 0; p < NP; p++)
        {
            uc = uhat * (double)cx[p] + vhat * (double)cy[p] + what * (double)cz[p];
            eq = wp[p] * (uc / cs2 + uc * uc / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
            f1[INDEX_F(i, j, k, p, RED)] = rho_RED_i * (eq + phi_eq[RED][p]);
            f1[INDEX_F(i, j, k, p, BLUE)] = rho_BLUE_i * (eq + phi_eq[BLUE][p]);
        }
    }
}