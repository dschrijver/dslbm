#include <mpi.h>
#include <stdio.h>
#include <math.h>

#include "include/memory.h"
#include "params.h"
#include "include/initialize.h"
#include "include/stencil.h"
#include "include/fields.h"
#include "include/forcing.h"
#include "include/output.h"
#include "include/collide.h"
#include "include/stream.h"
#include "include/communicate.h"
#include "include/active.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    SimulationBag *sim;
    DistributionBag *dists;
    FieldBag *fields;
    ParamBag *params;
    Stencil *stencil;
    MatrixBag *matrices;

    allocate_bags(&sim, &dists, &fields, &params, &stencil, &matrices);

    set_params(params);

    initialize_MPI(params);

    initialize_HDF5(params);

    initialize_stencil(sim);

    allocate_distributions(sim);

    allocate_fields(sim);

    allocate_matrices(sim);

    params->t = 0;
    params->t_output = 0;
    params->t_log = 0;
    params->n_output = 0;
    
    initialize_fields(sim);

    compute_chemical_potentials(sim);
    compute_molecular_field(sim);
    compute_nematic_stress(sim);
    compute_preliminary_forces(sim);
    compute_final_forces(sim);

    initialize_distributions(sim);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    double start_substep = 0.0, duration_substep = 0.0;
    char output_info[128];

    while (params->t < params->NTIME)
    {

        MPI_Barrier(MPI_COMM_WORLD);
        double start_timestep = MPI_Wtime();

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
            communicate_dist(dists->f2_a, sim);
            communicate_dist(dists->f2_p, sim);
        )

        TIME("> Streaming...",
            stream_distributions(sim);
        )

        TIME("> Computing macroscopic fields...",
            extract_moments(sim);
            compute_chemical_potentials(sim);
            compute_molecular_field(sim);
            compute_nematic_stress(sim);
            compute_preliminary_forces(sim);
            compute_final_velocity(sim);
            compute_final_forces(sim);
        )

        double duration_timestep = MPI_Wtime() - start_timestep;

        // LOGGING INFORMATION
        if (params->t_log == params->t)   
        {
            int current_total_VmRSS = get_VmRSS(params);
            if (params->process_rank == 0)
            {
                printf("--------------------------------------------------------------------------------\n");
                printf("Step completed!\n");
                printf("    Total memory in use (VmRSS): %'d kB\n", current_total_VmRSS);
                printf("    Duration of time step: %.4fs\n", duration_timestep);
                double total_time = MPI_Wtime() - start_time;
                double hours = floor(total_time / 3600.0);
                total_time -= hours*3600.0;
                double minutes = floor(total_time / 60.0);
                total_time -= minutes*60.0;
                double seconds = floor(total_time);
                printf("    Total simulation time: %02d:%02d:%02d\n", (int)hours, (int)minutes, (int)seconds);
                double remaining_time = (MPI_Wtime() - start_time) / (double)(params->t + 1) * (double)(params->NTIME - params->t - 1);
                hours = floor(remaining_time / 3600.0);
                remaining_time -= hours*3600.0;
                minutes = floor(remaining_time / 60.0);
                remaining_time -= minutes*60.0;
                seconds = floor(remaining_time);
                printf("    Remaining simulation time: %02d:%02d:%02d\n", (int)hours, (int)minutes, (int)seconds);
            }
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