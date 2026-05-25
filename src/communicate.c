#include <string.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/communicate.h"

void communicate_field(double *field, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    double *send_buffer = fields->send_buffer;
    double *recv_buffer = fields->recv_buffer;

    int *process_neighbors = params->process_neighbors;
    MPI_Comm comm_xslices = params->comm_xslices;

    MPI_Status status_first;

    int buffer_size = 2 * NY * NZ * sizeof(double);
    int buffer_number = 2 * NY * NZ;

    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &field[INDEX(i_start, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], 0, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], MPI_ANY_TAG, comm_xslices, &status_first);
    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(&field[INDEX(i_end, 0, 0)], recv_buffer, buffer_size);
    }

    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &field[INDEX(i_end - 2, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], 0, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], MPI_ANY_TAG, comm_xslices, &status_first);
    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(&field[INDEX(i_start - 2, 0, 0)], recv_buffer, buffer_size);
    }
}

void communicate_dist(double *dist, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    double *send_buffer = dists->send_buffer;
    double *recv_buffer = dists->recv_buffer;

    int *process_neighbors = params->process_neighbors;
    MPI_Comm comm_xslices = params->comm_xslices;

    MPI_Status status_first;

    int buffer_size = NY * NZ * NP * sizeof(double);
    int buffer_number = NY * NZ * NP;

    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &dist[INDEX_F(i_start, 0, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], 0, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], MPI_ANY_TAG, comm_xslices, &status_first);
    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(&dist[INDEX_F(i_end, 0, 0, 0)], recv_buffer, buffer_size);
    }

    if (process_neighbors[1] != MPI_PROC_NULL)
    {
        memcpy(send_buffer, &dist[INDEX_F(i_end - 1, 0, 0, 0)], buffer_size);
    }
    MPI_Sendrecv(send_buffer, buffer_number, MPI_DOUBLE, process_neighbors[1], 0, recv_buffer, buffer_number, MPI_DOUBLE, process_neighbors[0], MPI_ANY_TAG, comm_xslices, &status_first);
    if (process_neighbors[0] != MPI_PROC_NULL)
    {
        memcpy(&dist[INDEX_F(i_start - 1, 0, 0, 0)], recv_buffer, buffer_size);
    }
}