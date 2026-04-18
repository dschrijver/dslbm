#include <string.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/communicate.h"

void communicate_fields(SimulationBag *sim)
{
    ParamBag *params = sim->params;

    if (params->comm_xslices == MPI_COMM_NULL)
        return;

    GlobalFieldBag *glob_fields = sim->glob_fields;

    double *send_buffer = glob_fields->send_buffer;
    double *recv_buffer = glob_fields->recv_buffer;

    communicate_field(glob_fields->rho_N, send_buffer, recv_buffer, 0, sim);
}

void communicate_field(double *field, double *send_buffer, double *recv_buffer, int tag, SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int NY = params->NY;
    int NZ = params->NZ;
    int i_start = params->i_start;
    int i_end = params->i_end;
    int *process_neighbors = params->process_neighbors;
    MPI_Comm comm_xslices = params->comm_xslices;

    MPI_Status status_first;

    int buffer_size = 2 * NY * NZ * sizeof(double);
    int buffer_number = 2 * NY * NZ;

    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &field[INDEX_GLOB(i_start, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], tag, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], tag, comm_xslices, &status_first);
    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(&field[INDEX_GLOB(i_end, 0, 0)], recv_buffer, buffer_size);
    }

    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &field[INDEX_GLOB(i_end - 2, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], tag + 1, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], tag + 1, comm_xslices, &status_first);
    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(&field[INDEX_GLOB(i_start - 2, 0, 0)], recv_buffer, buffer_size);
    }
}

void communicate_dists(SimulationBag *sim)
{
    ParamBag *params = sim->params;

    if (params->comm_xslices == MPI_COMM_NULL)
        return;

    DistributionBag *dists = sim->dists;

    double *send_buffer = dists->send_buffer;
    double *recv_buffer = dists->recv_buffer;

    communicate_dist(dists->f2, send_buffer, recv_buffer, 2, sim);
}

void communicate_dist(double *dist, double *send_buffer, double *recv_buffer, int tag, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int i_start = params->i_start;
    int i_end = params->i_end;
    int *process_neighbors = params->process_neighbors;
    MPI_Comm comm_xslices = params->comm_xslices;

    MPI_Status status_first;

    int buffer_size = NCOMP * NY * NZ * NP * sizeof(double);
    int buffer_number = NCOMP * NY * NZ * NP;

    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &dist[INDEX_F(i_start, 0, 0, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], tag, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], tag, comm_xslices, &status_first);
    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(&dist[INDEX_F(i_end, 0, 0, 0, 0)], recv_buffer, buffer_size);
    }

    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &dist[INDEX_F(i_end - 1, 0, 0, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], tag + 1, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], tag + 1, comm_xslices, &status_first);
    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(&dist[INDEX_F(i_start - 1, 0, 0, 0, 0)], recv_buffer, buffer_size);
    }
}