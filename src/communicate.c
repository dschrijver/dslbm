#include <string.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/communicate.h"

void communicate_field(double *field, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    MPI_Comm comm_cart = params->comm_cart;

    MPI_Status status_first;

    double *send_buffer_x = fields->send_buffer_x;
    double *send_buffer_y = fields->send_buffer_y;
    double *send_buffer_z = fields->send_buffer_z;

    double *recv_buffer_x = fields->recv_buffer_x;
    double *recv_buffer_y = fields->recv_buffer_y;
    double *recv_buffer_z = fields->recv_buffer_z;

    int buffer_number_x = 2 * NY_proc * NZ_proc;
    int buffer_number_y = 2 * (NX_proc + 4) * NZ_proc;
    int buffer_number_z = 2 * (NX_proc + 4) * (NY_proc + 4);
 
    // Communicate to left
    for (int i = i_start; i < i_start + 2; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                send_buffer_x[NY_proc*NZ_proc*(i - i_start) + NZ_proc*(j - j_start) + (k - k_start)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_end; i < i_end + 2; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_x[NY_proc*NZ_proc*(i - i_end) + NZ_proc*(j - j_start) + (k - k_start)];
            }
        }
    }

    // Communicate to right
    for (int i = i_end - 2; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                send_buffer_x[NY_proc*NZ_proc*(i - i_end + 2) + NZ_proc*(j - j_start) + (k - k_start)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 2; i < i_start; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_x[NY_proc*NZ_proc*(i - i_start + 2) + NZ_proc*(j - j_start) + (k - k_start)];
            }
        }
    }

    // Communicate to bottom
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start; j < j_start + 2; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                send_buffer_y[2*NZ_proc*(i - i_start + 2) + NZ_proc*(j - j_start) + (k - k_start)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_end; j < j_end + 2; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_y[2*NZ_proc*(i - i_start + 2) + NZ_proc*(j - j_end) + (k - k_start)];
            }
        }
    }

    // Communicate to top
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_end - 2; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                send_buffer_y[2*NZ_proc*(i - i_start + 2) + NZ_proc*(j - j_end + 2) + (k - k_start)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start - 2; j < j_start; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_y[2*NZ_proc*(i - i_start + 2) + NZ_proc*(j - j_start + 2) + (k - k_start)];
            }
        }
    }

    // Communicate to back
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start - 2; j < j_end + 2; j++)
        {
            for (int k = k_start; k < k_start + 2; k++)
            {
                send_buffer_z[(NY_proc + 4)*2*(i - i_start + 2) + 2*(j - j_start + 2) + (k - k_start)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start - 2; j < j_end + 2; j++)
        {
            for (int k = k_end; k < k_end + 2; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_z[(NY_proc + 4)*2*(i - i_start + 2) + 2*(j - j_start + 2) + (k - k_end)];
            }
        }
    }

    // Communicate to front
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start - 2; j < j_end + 2; j++)
        {
            for (int k = k_end - 2; k < k_end; k++)
            {
                send_buffer_z[(NY_proc + 4)*2*(i - i_start + 2) + 2*(j - j_start + 2) + (k - k_end + 2)] = field[INDEX(i, j, k)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 2; i < i_end + 2; i++)
    {
        for (int j = j_start - 2; j < j_end + 2; j++)
        {
            for (int k = k_start - 2; k < k_start; k++)
            {
                field[INDEX(i, j, k)] = recv_buffer_z[(NY_proc + 4)*2*(i - i_start + 2) + 2*(j - j_start + 2) + (k - k_start + 2)];
            }
        }
    }
}

void communicate_dist(double *dist, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    double *send_buffer_x = dists->send_buffer_x;
    double *send_buffer_y = dists->send_buffer_y;
    double *send_buffer_z = dists->send_buffer_z;

    double *recv_buffer_x = dists->recv_buffer_x;
    double *recv_buffer_y = dists->recv_buffer_y;
    double *recv_buffer_z = dists->recv_buffer_z;

    MPI_Comm comm_cart = params->comm_cart;

    MPI_Status status_first;

    int buffer_number_x = NY_proc * NZ_proc * NP;
    int buffer_number_y = (NX_proc + 2) * NZ_proc * NP;
    int buffer_number_z = (NX_proc + 2) * (NY_proc + 2) * NP;

    // Communicate to left
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_x[NZ_proc*NP*(j - j_start) + NP*(k - k_start) + p] = dist[INDEX_F(i_start, j, k, p)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, MPI_ANY_TAG, comm_cart, &status_first);
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i_end, j, k, p)] = recv_buffer_x[NZ_proc*NP*(j - j_start) + NP*(k - k_start) + p];
            }
        }
    }

    // Communicate to right
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_x[NZ_proc*NP*(j - j_start) + NP*(k - k_start) + p] = dist[INDEX_F(i_end - 1, j, k, p)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, MPI_ANY_TAG, comm_cart, &status_first);
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i_start - 1, j, k, p)] = recv_buffer_x[NZ_proc*NP*(j - j_start) + NP*(k - k_start) + p];
            }
        }
    }

    // Communicate to bottom
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_y[NZ_proc*NP*(i - i_start + 1) + NP*(k - k_start) + p] = dist[INDEX_F(i, j_start, k, p)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i, j_end, k, p)] = recv_buffer_y[NZ_proc*NP*(i - i_start + 1) + NP*(k - k_start) + p];
            }
        }
    }

    // Communicate to top
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_y[NZ_proc*NP*(i - i_start + 1) + NP*(k - k_start) + p] = dist[INDEX_F(i, j_end - 1, k, p)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i, j_start - 1, k, p)] = recv_buffer_y[NZ_proc*NP*(i - i_start + 1) + NP*(k - k_start) + p];
            }
        }
    }

    // Communicate to back
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int j = j_start - 1; j < j_end + 1; j++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_z[(NY_proc + 2)*NP*(i - i_start + 1) + NP*(j - j_start + 1) + p] = dist[INDEX_F(i, j, k_start, p)];
            }
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int j = j_start - 1; j < j_end + 1; j++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i, j, k_end, p)] = recv_buffer_z[(NY_proc + 2)*NP*(i - i_start + 1) + NP*(j - j_start + 1) + p];
            }
        }
    }

    // Communicate to front
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int j = j_start - 1; j < j_end + 1; j++)
        {
            for (int p = 0; p < NP; p++)
            {
                send_buffer_z[(NY_proc + 2)*NP*(i - i_start + 1) + NP*(j - j_start + 1) + p] = dist[INDEX_F(i, j, k_end - 1, p)];
            }   
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start - 1; i < i_end + 1; i++)
    {
        for (int j = j_start - 1; j < j_end + 1; j++)
        {
            for (int p = 0; p < NP; p++)
            {
                dist[INDEX_F(i, j, k_start - 1, p)] = recv_buffer_z[(NY_proc + 2)*NP*(i - i_start + 1) + NP*(j - j_start + 1) + p];
            }
        }
    }
}
