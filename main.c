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
#include "definitions.h"

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

    communicate_fields(sim);

    compute_Q_corrections(sim);

    evaluate_color_gradients(sim);

    communicate_surface_vector(sim);

    evaluate_forces(sim);

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
            communicate_dists(sim);
        )

        TIME("> Streaming...",
            stream_distributions(sim);
        )

        // TIME("> Wetnode boundary conditions...",
        //     wetnode_boundary_conditions(sim);
        // )

        TIME("> Computing macroscopic fields...",
            extract_moments(sim);
            evaluate_pressure(sim);
            communicate_fields(sim);
            compute_Q_corrections(sim);
            evaluate_color_gradients(sim);
            communicate_surface_vector(sim);
            evaluate_forces(sim);
            update_final_velocity(sim);
            M_RED = evaluate_mass(RED, sim);
            M_BLUE = evaluate_mass(BLUE, sim);
            M_total = M_RED + M_BLUE;
            if (M_total != M_total)
            {
                if (params->process_rank == 0)
                {
                    printf("\n--------------------------------------------------------------------------------\n");
                    printf("Step failed, mass is NaN!\n");
                }
                break;
            }
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