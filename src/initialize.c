#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void dims_create(int dims[3], ParamBag *params)
{
    const int x_max = params->NX;
    const int y_max = params->NY;
    const int z_max = params->NZ;

    const int N_proc = params->number_of_processes;
    int A = 1000000000;

    for (int x = 1; x <= x_max; x++)
    {
        for (int y = 1; y <= y_max; y++)
        {
            for (int z = 1; z <= z_max; z++)
            {
                const int N_proc_local = x*y*z;

                if (N_proc_local != N_proc)
                    continue;
                
                const int x_length = (params->NX + x - 1) / x;
                const int y_length = (params->NY + y - 1) / y;
                const int z_length = (params->NZ + z - 1) / z;

                const int A_local = 2*(x_length*y_length + x_length*z_length + y_length*z_length);

                if (A_local < A)
                {
                    A = A_local;
                    dims[0] = x;
                    dims[1] = y;
                    dims[2] = z;
                }
            }
        }
    }
}

int block_start(int coord, int n, int p) 
{
    return (int)(((long)coord * n) / p);
}

void initialize_MPI(ParamBag *params)
{
    MPI_Comm_size(MPI_COMM_WORLD, &params->number_of_processes);

    int dims[3] = {0};

    dims_create(dims, params);

    MPI_Comm_rank(MPI_COMM_WORLD, &params->process_rank);

    if (dims[0] == 0)
    {
        if (params->process_rank == 0)
        {
            printf("Couldn't find decomposition!\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (params->process_rank == 0)
    {
        printf("Domain decomposition: [%d][%d][%d]\n", dims[0], dims[1], dims[2]);
    }

    int xperiodic = 0;
    int yperiodic = 0;
    int zperiodic = 0;
#ifdef XPERIODIC
    xperiodic = 1;
#endif
#ifdef YPERIODIC
    yperiodic = 1;
#endif
#ifdef ZPERIODIC
    zperiodic = 1;
#endif
    int periods[3] = {xperiodic, yperiodic, zperiodic};
    int reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &params->comm_cart);
    MPI_Comm_rank(params->comm_cart, &params->process_rank);
    MPI_Cart_get(params->comm_cart, 3, dims, periods, params->process_coords);
    MPI_Cart_shift(params->comm_cart, 0, 1, &params->process_left, &params->process_right);
    MPI_Cart_shift(params->comm_cart, 1, 1, &params->process_bottom, &params->process_top);
    MPI_Cart_shift(params->comm_cart, 2, 1, &params->process_back, &params->process_front);

    params->i_start = block_start(params->process_coords[0], params->NX, dims[0]);
    params->i_end = block_start((params->process_coords[0] + 1), params->NX, dims[0]);

    params->j_start = block_start(params->process_coords[1], params->NY, dims[1]);
    params->j_end = block_start((params->process_coords[1] + 1), params->NY, dims[1]);

    params->k_start = block_start(params->process_coords[2], params->NZ, dims[2]);
    params->k_end = block_start((params->process_coords[2] + 1), params->NZ, dims[2]);

    params->NX_proc = params->i_end - params->i_start;
    params->NY_proc = params->j_end - params->j_start;
    params->NZ_proc = params->k_end - params->k_start;

    MPI_Type_contiguous(9, MPI_DOUBLE, &params->MPI_MATRIX);
    MPI_Type_commit(&params->MPI_MATRIX);

    srand(time(NULL) + params->process_rank);
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
    int NY_proc = params->NY_proc;
    int NZ_proc = params->NZ_proc;

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    params->filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[3] = {NX_proc + 2, NY_proc + 2, NZ_proc + 2};
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
    WRITE_FIELD(phi)
    WRITE_FIELD(rho_a)
    WRITE_FIELD(rho_p)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    WRITE_FIELD(u_a)
    WRITE_FIELD(v_a)
    WRITE_FIELD(w_a)
    WRITE_FIELD(u_p)
    WRITE_FIELD(v_p)
    WRITE_FIELD(w_p)
    matrix * restrict const Q = matrices->Q;

    PARAM(rho_0_a)
    PARAM(rho_0_p)

#ifdef   DIM_2D
    const double d = 2.0;
#else // DIM_3D
    const double d = 3.0;
#endif

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

#ifdef INI_FLUID_1
        rho_a[idx0] = rho_0_a;
        rho_p[idx0] = rho_0_p;

        u[idx0] = 0.0;
        v[idx0] = 0.0;
        w[idx0] = 0.0;

#ifdef INI_ACTIVE_1
        const double phi_i = rho_0_a / (rho_0_a + rho_0_p);

        double x = (double)rand() / ((double)RAND_MAX + 1.0);
        const double phi_n = x*2.0*DS_PI;
#ifdef   DIM_3D
        x = (double)rand() / ((double)RAND_MAX + 1.0);
        const double theta_n = x*0.5*DS_PI;
        double n[3] = {sin(theta_n)*cos(phi_n), sin(theta_n)*sin(phi_n), cos(theta_n)};
#else // DIM_2D
        double n[3] = {cos(phi_n), sin(phi_n), 0.0};
#endif
#endif

        for (int a = 0; a < 3; a++)
        {
            for (int b = 0; b < 3; b++)
            {
                const double S = params->S_nem * sqrt(phi_i);
                Q[idx0][a][b] = d/(d-1.0) * S * (n[a]*n[b] - delta(a, b)/d);
            }
        }
#endif

        rho[idx0] = rho_a[idx0] + rho_p[idx0];
        phi[idx0] = rho_a[idx0] / rho[idx0];

        u_a[idx0] = u[idx0];
        v_a[idx0] = v[idx0];
        w_a[idx0] = w[idx0];
        u_p[idx0] = u[idx0];
        v_p[idx0] = v[idx0];
        w_p[idx0] = w[idx0];
    }
}

void initialize_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho_a)
    READ_FIELD(rho_p)
    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)
    READ_FIELD(u_p)
    READ_FIELD(v_p)
    READ_FIELD(w_p)
    READ_FIELD(Fx_a)
    READ_FIELD(Fy_a)
    READ_FIELD(Fz_a)
    READ_FIELD(Fx_p)
    READ_FIELD(Fy_p)
    READ_FIELD(Fz_p)

    WRITE_DIST(f1_a)
    WRITE_DIST(f1_p)

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double rho_a_i = rho_a[idx0];
        const double rho_p_i = rho_p[idx0];

        const double Fx_a_i = Fx_a[idx0];
        const double Fy_a_i = Fy_a[idx0];
        const double Fz_a_i = Fz_a[idx0];
        const double Fx_p_i = Fx_p[idx0];
        const double Fy_p_i = Fy_p[idx0];
        const double Fz_p_i = Fz_p[idx0];

        const double u_a_i = u_a[idx0] - 0.5 / rho_a_i * Fx_a_i;
        const double v_a_i = v_a[idx0] - 0.5 / rho_a_i * Fy_a_i;
        const double w_a_i = w_a[idx0] - 0.5 / rho_a_i * Fz_a_i;
        const double u_p_i = u_p[idx0] - 0.5 / rho_p_i * Fx_p_i;
        const double v_p_i = v_p[idx0] - 0.5 / rho_p_i * Fy_p_i;
        const double w_p_i = w_p[idx0] - 0.5 / rho_p_i * Fz_p_i;

        compute_equilibrium(rho_a_i, u_a_i, v_a_i, w_a_i, &f1_a[INDEX_F(i, j, k, 0)], sim);
        compute_equilibrium(rho_p_i, u_p_i, v_p_i, w_p_i, &f1_p[INDEX_F(i, j, k, 0)], sim);
    }
}