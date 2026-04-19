#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <hdf5.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"
#include "../include/initialize.h"

void initialize_MPI(ParamBag *params)
{
    MPI_Comm_size(MPI_COMM_WORLD, &params->number_of_processes);

    params->number_of_processes = min(params->number_of_processes, params->NX / 2);

    int dims[3] = {params->number_of_processes, 1, 1};
    int xperiodic = 0;
#ifdef XPERIODIC
    xperiodic = 1;
#endif
    int periods[3] = {xperiodic, 0, 0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &params->comm_xslices);

    if (params->comm_xslices == MPI_COMM_NULL)
    {
        params->i_start = 0;
        params->i_end = 0;
        params->NX_proc = 0;
        params->process_rank = -1;
        printf("test\n");
        fflush(stdout);
    }
    else
    {
        MPI_Comm_rank(params->comm_xslices, &params->process_rank);
        MPI_Cart_get(params->comm_xslices, 3, dims, periods, params->process_coords);
        MPI_Cart_shift(params->comm_xslices, 0, 1, &params->process_neighbors[0], &params->process_neighbors[1]);

        params->i_start = (float)(params->process_coords[0]) / (float)params->number_of_processes * params->NX;
        params->i_end = (float)(params->process_coords[0] + 1) / (float)params->number_of_processes * params->NX;
        params->NX_proc = params->i_end - params->i_start;
    }
}

void initialize_HDF5(ParamBag *params)
{
    params->fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(params->fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    params->scalar_space = H5Screate(H5S_SCALAR);

    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;
    int NX_proc = params->NX_proc;

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    params->filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc_glob[3] = {NX_proc + 4, NY, NZ};
    params->memspace_glob = H5Screate_simple(3, dims_proc_glob, NULL);

    hsize_t dims_proc_comp[4] = {NX_proc + 4, NY, NZ, NCOMP};
    params->memspace_comp = H5Screate_simple(4, dims_proc_comp, NULL);

    params->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    params->dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(params->dxpl_id, H5FD_MPIO_COLLECTIVE);
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

    (void)NX;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double rho_0_RED = params->rho_0_RED;
    double rho_0_BLUE = params->rho_0_BLUE;

    (void)rho_0_BLUE;

    double zeta = stencil->zeta;
    double alpha_RED = params->alpha_RED;
    double alpha_BLUE = params->alpha_BLUE;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;
    double *rho_N = glob_fields->rho_N;

    double *rho_comp = comp_fields->rho_comp;

#ifdef INI_POISEUILLE
    FOR_DOMAIN
    {
        u[INDEX_GLOB(i, j, k)] = 0.0;
        v[INDEX_GLOB(i, j, k)] = 0.0;
        w[INDEX_GLOB(i, j, k)] = 0.0;

        rho_comp[INDEX(i, j, k, RED)] = rho_0_RED;
        rho_comp[INDEX(i, j, k, BLUE)] = 0.0;
    }
#endif

#ifdef INI_DROPLET
    double x, y, z, r;
    FOR_DOMAIN
    {
        x = physx(i) - INI_DROPLET_X;
        y = physy(j) - INI_DROPLET_Y;
        z = physz(k) - INI_DROPLET_Z;
        r = sqrt(x * x + y * y + z * z);

        rho_comp[INDEX(i, j, k, RED)] = 0.5 * rho_0_RED * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        rho_comp[INDEX(i, j, k, BLUE)] = 0.5 * rho_0_BLUE * (1.0 + tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        u[INDEX_GLOB(i, j, k)] = 0.5 * INI_DROPLET_U * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        v[INDEX_GLOB(i, j, k)] = 0.5 * INI_DROPLET_V * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        w[INDEX_GLOB(i, j, k)] = 0.5 * INI_DROPLET_W * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
    }
#endif

#ifdef INI_TWOCOMPONENT_POISEUILLE
    double y, r;
    FOR_DOMAIN
    {
        y = physy(j) - INI_DROPLET_Y;
        r = fabs(y);

        rho_comp[INDEX(i, j, k, RED)] = 0.5 * rho_0_RED * (1.0 - tanh((r - INI_TWOCOMPONENT_POISEUILLE_A) / INI_DROPLET_SF));
        rho_comp[INDEX(i, j, k, BLUE)] = 0.5 * rho_0_BLUE * (1.0 + tanh((r - INI_TWOCOMPONENT_POISEUILLE_A) / INI_DROPLET_SF));
        u[INDEX_GLOB(i, j, k)] = 0.0;
        v[INDEX_GLOB(i, j, k)] = 0.0;
        w[INDEX_GLOB(i, j, k)] = 0.0;
    }
#endif

#ifdef INI_TWOCOMPONENT_COUETTE
    double x;
    double width = physlx(params);
    FOR_DOMAIN
    {
        x = physy(i);

        rho_comp[INDEX(i, j, k, RED)] = 0.5 * rho_0_RED * (1.0 - tanh((x - 0.5*width) / INI_TWOCOMPONENT_COUETTE_SF));
        rho_comp[INDEX(i, j, k, BLUE)] = 0.5 * rho_0_BLUE * (1.0 + tanh((x - 0.5*width) / INI_TWOCOMPONENT_COUETTE_SF));

        u[INDEX_GLOB(i, j, k)] = 0.0;
        v[INDEX_GLOB(i, j, k)] = INI_TWOCOMPONENT_COUETTE_V_LEFT;
        w[INDEX_GLOB(i, j, k)] = 0.0;
    }
#endif

    FOR_DOMAIN
    {
        rho[INDEX_GLOB(i, j, k)] = rho_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)];
        rho_N[INDEX_GLOB(i, j, k)] = (rho_comp[INDEX(i, j, k, RED)] / params->rho_0_RED - rho_comp[INDEX(i, j, k, BLUE)] / params->rho_0_BLUE) / (rho_comp[INDEX(i, j, k, RED)] / params->rho_0_RED + rho_comp[INDEX(i, j, k, BLUE)] / params->rho_0_BLUE);
        pressure[INDEX_GLOB(i, j, k)] = rho_comp[INDEX(i, j, k, RED)] * zeta * (1.0 - alpha_RED) + rho_comp[INDEX(i, j, k, BLUE)] * zeta * (1.0 - alpha_BLUE);
    }
}

void initialize_flags(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    int *flag = glob_fields->flag;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = -2; j < NY + 2; j++)
        {
            for (int k = -2; k < NZ + 2; k++)
            {
                flag[INDEX_FLAG(i, j, k)] = FLUID;

#ifdef LEFT_HWBB_NOSLIP
                if (i < 0)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef RIGHT_HWBB_NOSLIP
                if (i > params->NX - 1)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef BOTTOM_HWBB_NOSLIP
                if (j < 0)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef TOP_HWBB_NOSLIP
                if (j > NY - 1)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef BACK_HWBB_NOSLIP
                if (k < 0)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef FRONT_HWBB_NOSLIP
                if (k > NZ - 1)
                    flag[INDEX_FLAG(i, j, k)] = BOUNDARY;
#endif

#ifdef LEFT_NEBB_VELOCITY
                if (i < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef RIGHT_NEBB_VELOCITY
                if (i > params->NX - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef BOTTOM_NEBB_VELOCITY
                if (j < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef TOP_NEBB_VELOCITY
                if (j > NY - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef BACK_NEBB_VELOCITY
                if (k < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef FRONT_NEBB_VELOCITY
                if (k > NZ - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef LEFT_NEBB_PRESSURE
                if (i < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef RIGHT_NEBB_PRESSURE
                if (i > params->NX - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef BOTTOM_NEBB_PRESSURE
                if (j < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef TOP_NEBB_PRESSURE
                if (j > NY - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef BACK_NEBB_PRESSURE
                if (k < 0)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif

#ifdef FRONT_NEBB_PRESSURE
                if (k > NZ - 1)
                    flag[INDEX_FLAG(i, j, k)] = WETNODE;
#endif
            }
        }
    }
}

void initialize_distributions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    DistributionBag *dists = sim->dists;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    Stencil *stencil = sim->stencil;

    double rho_i, u_i, v_i, w_i;
    double cs2_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *cs2 = stencil->cs2;

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
        u_i = u[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fx[INDEX(i, j, k, RED)] + Fx[INDEX(i, j, k, BLUE)]);
        v_i = v[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fy[INDEX(i, j, k, RED)] + Fy[INDEX(i, j, k, BLUE)]);
        w_i = w[INDEX_GLOB(i, j, k)] - 1.0 / (2.0 * rho_i) * (Fz[INDEX(i, j, k, RED)] + Fz[INDEX(i, j, k, BLUE)]);

        for (int n = 0; n < NCOMP; n++)
        {
            rho_i = rho_comp[INDEX(i, j, k, n)];
            cs2_i = cs2[n];

            compute_equilibrium(rho_i, u_i, v_i, w_i, cs2_i, &f(0), sim);
        }
    }
}

double physx(int i)
{
#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
    return (double)i;
#else
    return (double)i + 0.5;
#endif
}

double physy(int j)
{
#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    return (double)j;
#else
    return (double)j + 0.5;
#endif
}

double physz(int k)
{
#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    return (double)k;
#else
    return (double)k + 0.5;
#endif
}

double physlx(ParamBag *params)
{
    double result = (double)params->NX;
#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}

double physly(ParamBag *params)
{
    double result = (double)params->NY;
#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(TOP_NEBB_VELOCITY) || defined(TOP_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}

double physlz(ParamBag *params)
{
    double result = (double)params->NZ;
#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    result -= 0.5;
#endif

#if defined(FRONT_NEBB_VELOCITY) || defined(FRONT_NEBB_PRESSURE)
    result -= 0.5;
#endif
    return result;
}