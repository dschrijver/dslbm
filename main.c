#include <mpi.h>
#include <stdio.h>

#include "include/datatypes.h"
#include "include/memory.h"
#include "params.h"
#include "include/initialize.h"
#include "include/stencil.h"
#include "include/forcing.h"
#include "include/output.h"
#include "include/communicate.h"
#include "include/collide.h"
#include "include/stream.h"
#include "include/wetnode.h"
#include "include/fields.h"
#include "include/wetnode.h"
#include "definitions.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    SimulationBag *sim;
    DistributionBag *dists;
    FieldBag *fields;
    ParamBag *params;
    Stencil *stencil;

    allocate_bags(&sim, &dists, &fields, &params, &stencil);

    set_params(params);

    initialize_MPI(params);

    initialize_HDF5(params);

    initialize_stencil(sim);

    allocate_distributions(sim);

    allocate_fields(sim);

    params->t = 0;
    params->t_output = 0;
    params->t_log = 0;
    params->n_output = 0;

    initialize_flags(sim);
    
    initialize_fields(sim);

    communicate_field(fields->rho_N, sim);

    evaluate_color_gradients(sim);

    communicate_field(fields->nx, sim);
    communicate_field(fields->ny, sim);
    communicate_field(fields->nz, sim);

    evaluate_forces(sim);

    communicate_field(fields->rho, sim);
    communicate_field(fields->pressure, sim);
    communicate_field(fields->u, sim);
    communicate_field(fields->v, sim);
    communicate_field(fields->w, sim);

    compute_Q_corrections(sim);

    initialize_distributions(sim);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    double start_timestep, duration_timestep;
    double start_substep, duration_substep;
    char output_info[128];
    double M_RED, M_BLUE, M_total;

    while (params->t < params->NTIME)
    {

        MPI_Barrier(MPI_COMM_WORLD);
        start_timestep = MPI_Wtime();

        // LOGGING TIME
        if ((params->process_rank == 0) && (params->t_log == params->t)) 
        {
            printf("================================================================================\n");
            printf("Time: %d\n", params->t);
            printf("--------------------------------------------------------------------------------\n");
        }

        // OUTPUT
        if (params->t_output == params->t)
        {
            TIME_OUTPUT(
                output_data(sim);
                params->t_output += params->NSTORE;
            )
        }

        TIME("> Collision...",
            collide(sim);
        )

        TIME("> Communicate distributions...",
            communicate_dist(dists->f2_RED, sim);
            communicate_dist(dists->f2_BLUE, sim);
        )

        TIME("> Streaming...",
            stream_distributions(sim);
        )

        TIME("> Computing macroscopic fields...",
            extract_moments(sim);
            evaluate_pressure(sim);

            wetnode_macroscopic_fields(sim);

            communicate_field(fields->rho_N, sim);
            evaluate_color_gradients(sim);

            communicate_field(fields->nx, sim);
            communicate_field(fields->ny, sim);
            communicate_field(fields->nz, sim);
            evaluate_forces(sim);

            wetnode_distributions(sim);

            extract_moments(sim);
            update_final_velocity(sim);

            communicate_field(fields->rho, sim);
            communicate_field(fields->pressure, sim);
            communicate_field(fields->u, sim);
            communicate_field(fields->v, sim);
            communicate_field(fields->w, sim);

            compute_Q_corrections(sim);
        )

        duration_timestep = MPI_Wtime() - start_timestep;

        // LOGGING INFORMATION
        if ((params->process_rank == 0) && (params->t_log == params->t))   
        {
            printf("--------------------------------------------------------------------------------\n");
            printf("Step completed!\n");
            printf("    RED mass: %.5e, BLUE mass: %.5e, Total mass: %.5e\n", M_RED, M_BLUE, M_total);
            printf("    Duration of time step: %.4fs\n", duration_timestep);
            printf("    Total simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0);
            printf("    Expected remaining simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0 / (double)(params->t + 1) * (double)(params->NTIME - params->t - 1));
            params->t_log += params->NLOG;
        }

        params->t++;
    }

    output_data(sim);

    if (params->process_rank == 0)
    {
        printf("\nSimulation done!\n");
    }

    free_all(sim);

    MPI_Finalize();
    return 0;
}