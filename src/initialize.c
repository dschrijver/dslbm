#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <hdf5.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"
#include "../include/fields.h"
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
    hsize_t dims_proc[3] = {NX_proc + 4, NY, NZ};
    params->memspace = H5Screate_simple(3, dims_proc, NULL);

    params->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    params->dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(params->dxpl_id, H5FD_MPIO_COLLECTIVE);
}

void initialize_fields(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    WRITE_FIELD(rho)
    WRITE_FIELD(rho_RED)
    WRITE_FIELD(rho_BLUE)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    WRITE_FIELD(rho_N)

    PARAM(rho_0_RED)
    PARAM(rho_0_BLUE)

#ifdef INI_POISEUILLE
    FOR_DOMAIN
    {
        rho_RED[INDEX(i, j, k)] = rho_0_RED;
        rho_BLUE[INDEX(i, j, k)] = 0.0;

        u[INDEX(i, j, k)] = 0.0;
        v[INDEX(i, j, k)] = 0.0;
        w[INDEX(i, j, k)] = 0.0;
    }
#endif

#ifdef INI_DROPLET
    FOR_DOMAIN
    {
        const double x = physx(i) - INI_DROPLET_X;
        const double y = physy(j) - INI_DROPLET_Y;
        const double z = physz(k) - INI_DROPLET_Z;
        const double r = sqrt(x * x + y * y + z * z);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        rho_BLUE[INDEX(i, j, k)] = 0.5 * rho_0_BLUE * (1.0 + tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));

        u[INDEX(i, j, k)] = 0.5 * INI_DROPLET_U * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        v[INDEX(i, j, k)] = 0.5 * INI_DROPLET_V * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
        w[INDEX(i, j, k)] = 0.5 * INI_DROPLET_W * (1.0 - tanh((r - INI_DROPLET_R) / INI_DROPLET_SF));
    }
#endif

#ifdef INI_BUBBLE
    FOR_DOMAIN
    {
        const double x = physx(i) - INI_BUBBLE_X;
        const double y = physy(j) - INI_BUBBLE_Y;
        const double z = physz(k) - INI_BUBBLE_Z;
        const double r = sqrt(x * x + y * y + z * z);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (1.0 + tanh((r - INI_BUBBLE_R) / INI_BUBBLE_SF));
        rho_BLUE[INDEX(i, j, k)] = 0.5 * rho_0_BLUE * (1.0 - tanh((r - INI_BUBBLE_R) / INI_BUBBLE_SF));

        u[INDEX(i, j, k)] = 0.5 * INI_BUBBLE_U * (1.0 - tanh((r - INI_BUBBLE_R) / INI_BUBBLE_SF));
        v[INDEX(i, j, k)] = 0.5 * INI_BUBBLE_V * (1.0 - tanh((r - INI_BUBBLE_R) / INI_BUBBLE_SF));
        w[INDEX(i, j, k)] = 0.5 * INI_BUBBLE_W * (1.0 - tanh((r - INI_BUBBLE_R) / INI_BUBBLE_SF));
    }
#endif

#ifdef INI_TWOCOMPONENT_POISEUILLE
    const double width = physlx(params);
    FOR_DOMAIN
    {
        const double x = physx(i) - 0.5 * width;
        const double r = fabs(x);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (1.0 + tanh((r - INI_TWOCOMPONENT_POISEUILLE_A) / INI_TWOCOMPONENT_POISEUILLE_SF));
        rho_BLUE[INDEX(i, j, k)] = 0.5 * rho_0_BLUE * (1.0 - tanh((r - INI_TWOCOMPONENT_POISEUILLE_A) / INI_TWOCOMPONENT_POISEUILLE_SF));
        u[INDEX(i, j, k)] = 0.0;
        v[INDEX(i, j, k)] = 0.0;
        w[INDEX(i, j, k)] = 0.0;
    }
#endif

#ifdef INI_TWOCOMPONENT_COUETTE
    const double width = physlx(params);
    FOR_DOMAIN
    {
        const double x = physx(i);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (1.0 - tanh((x - 0.5 * width) / INI_TWOCOMPONENT_COUETTE_SF));
        rho_BLUE[INDEX(i, j, k)] = 0.5 * rho_0_BLUE * (1.0 + tanh((x - 0.5 * width) / INI_TWOCOMPONENT_COUETTE_SF));

        u[INDEX(i, j, k)] = 0.0;

        if (i == 0)
        {
            v[INDEX(i, j, k)] = LEFT_V_VELOCITY;
        }
        else if (i == NX - 1)
        {
            v[INDEX(i, j, k)] = RIGHT_V_VELOCITY;
        }
        else
        {
            v[INDEX(i, j, k)] = 0.0;
        }

        w[INDEX(i, j, k)] = 0.0;
    }
#endif

#ifdef INI_TWODROPLETS
    FOR_DOMAIN
    {
        const double x1 = physx(i) - INI_TWODROPLETS_X1;
        const double y1 = physy(j) - INI_TWODROPLETS_Y1;
        const double z1 = physz(k) - INI_TWODROPLETS_Z1;
        const double r1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

        const double x2 = physx(i) - INI_TWODROPLETS_X2;
        const double y2 = physy(j) - INI_TWODROPLETS_Y2;
        const double z2 = physz(k) - INI_TWODROPLETS_Z2;
        const double r2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (2.0 - tanh((r1 - INI_TWODROPLETS_R1) / INI_TWODROPLETS_SF)  - tanh((r2 - INI_TWODROPLETS_R2) / INI_TWODROPLETS_SF));
        rho_BLUE[INDEX(i, j, k)] = rho_0_BLUE - 0.5 * rho_0_BLUE * (2.0 - tanh((r1 - INI_TWODROPLETS_R1) / INI_TWODROPLETS_SF)  - tanh((r2 - INI_TWODROPLETS_R2) / INI_TWODROPLETS_SF));

        u[INDEX(i, j, k)] = 0.5 * INI_TWODROPLETS_U1 * (1.0 - tanh((r1 - INI_TWODROPLETS_R1) / INI_TWODROPLETS_SF)) + 
                            0.5 * INI_TWODROPLETS_U2 * (1.0 - tanh((r2 - INI_TWODROPLETS_R2) / INI_TWODROPLETS_SF));
        u[INDEX(i, j, k)] = 0.5 * INI_TWODROPLETS_V1 * (1.0 - tanh((r1 - INI_TWODROPLETS_R1) / INI_TWODROPLETS_SF)) + 
                            0.5 * INI_TWODROPLETS_V2 * (1.0 - tanh((r2 - INI_TWODROPLETS_R2) / INI_TWODROPLETS_SF));
        u[INDEX(i, j, k)] = 0.5 * INI_TWODROPLETS_W1 * (1.0 - tanh((r1 - INI_TWODROPLETS_R1) / INI_TWODROPLETS_SF)) + 
                            0.5 * INI_TWODROPLETS_W2 * (1.0 - tanh((r2 - INI_TWODROPLETS_R2) / INI_TWODROPLETS_SF));
    }
#endif

#ifdef INI_LAYERED_POISEUILLE
    FOR_DOMAIN
    {
        const double x = physx(i);

        rho_RED[INDEX(i, j, k)] = 0.5 * rho_0_RED * (1.0 - tanh((x - 0.5*(double)NX) / INI_LAYERED_POISEUILLE_SF));
        rho_BLUE[INDEX(i, j, k)] = 0.5 * rho_0_BLUE * (1.0 + tanh((x - 0.5*(double)NX)  / INI_LAYERED_POISEUILLE_SF));

        u[INDEX(i, j, k)] = 0.0;
        v[INDEX(i, j, k)] = 0.0;
        w[INDEX(i, j, k)] = 0.0;
    }
#endif

    FOR_DOMAIN
    {
        rho[INDEX(i, j, k)] = rho_RED[INDEX(i, j, k)] + rho_BLUE[INDEX(i, j, k)];
        rho_N[INDEX(i, j, k)] = (rho_RED[INDEX(i, j, k)] / params->rho_0_RED - rho_BLUE[INDEX(i, j, k)] / params->rho_0_BLUE) / (rho_RED[INDEX(i, j, k)] / params->rho_0_RED + rho_BLUE[INDEX(i, j, k)] / params->rho_0_BLUE);
    }

    evaluate_pressure(sim);
}

void initialize_flags(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    int *flag = fields->flag;

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
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)

    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    FOR_DOMAIN
    {
        const int idx = INDEX(i, j, k);

        const double rho_i = rho[idx];
        const double rho_RED_i = rho_RED[idx];
        const double rho_BLUE_i = rho_BLUE[idx];
        const double pressure_i = pressure[idx];
        const double u_i = u[idx] - 0.5 * Fx[idx] / rho_i;
        const double v_i = v[idx] - 0.5 * Fy[idx] / rho_i;
        const double w_i = w[idx] - 0.5 * Fz[idx] / rho_i;

        compute_equilibrium(rho_i, u_i, v_i, w_i, pressure_i, &f1_RED[INDEX_F(i, j, k, 0)], sim);
        memcpy(&f1_BLUE[INDEX_F(i, j, k, 0)], &f1_RED[INDEX_F(i, j, k, 0)], NP * sizeof(double));

        for (int p = 0; p < NP; p++)
        {
            const int idxf = INDEX_F(i, j, k, p);

            f1_RED[idxf] *= rho_RED_i / rho_i;
            f1_BLUE[idxf] *= rho_BLUE_i / rho_i;
        }
    }
}