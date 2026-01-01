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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    SimulationBag *sim;
    DistributionBag *dists;
    GlobalFieldBag *glob_fields;
    ComponentFieldBag *comp_fields;
    ParamBag *params;
    Stencil *stencil;

    allocate_bags(&sim, &dists, &glob_fields, &comp_fields, &params, &stencil);

    set_params(params);

    initialize_MPI(params);

    initialize_stencil(sim);

    allocate_distributions(sim);

    allocate_fields(sim);

    params->t = 0;
    params->t_output = 0;
    params->n_output = 0;

    initialize_fields(sim);

    evaluate_forces(sim);

    initialize_distributions(sim);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    double start_timestep, duration_timestep;
    double start_substep, duration_substep;
    char output_info[128];

    while (params->t < params->NTIME)
    {

        MPI_Barrier(MPI_COMM_WORLD);
        start_timestep = MPI_Wtime();

        // LOGGING TIME
        if (params->process_rank == 0)
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

        TIME("> Communicate fields...",
            communicate_fields(sim);
        )

        TIME("> Collision...",
            collide_distributions_MRT(sim);
            collide_distributions_CGM(sim);
        )

        TIME("> Communicate distributions...",
            communicate_dists(sim);
        )

        TIME("> Streaming...",
            stream_distributions(sim);
        )

        TIME("> Wetnode boundary conditions...",
            wetnode_boundary_conditions(sim);
        )

        TIME("> Computing macroscopic fields...",
            extract_moments(sim);
            evaluate_forces(sim);
            update_final_velocity(sim);
        )

        duration_timestep = MPI_Wtime() - start_timestep;

        // LOGGING INFORMATION
        if (params->process_rank == 0)
        {
            printf("--------------------------------------------------------------------------------\n");
            printf("Step completed!\n");
            printf("    Duration of time step: %.4fs\n", duration_timestep);
            printf("    Total simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0);
            printf("    Expected remaining simulation time: %.2fh\n", (MPI_Wtime() - start_time) / 3600.0 / (double)(params->t + 1) * (double)(params->NTIME - params->t - 1));
        }

        params->t++;
    }

    output_data(sim);

    if (params->process_rank == 0)
    {
        printf("\nSimulation done!\n");
    }

    MPI_Finalize();
    return 0;
}