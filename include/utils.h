#ifndef UTILS_H
#define UTILS_H

enum flags
{
    FLUID,
    BOUNDARY,
    WETNODE
};

#define DS_PI 3.14159265358979323846
#define INDEX(i, j, k) (NY * NZ * (i + 2 - i_start) + NZ * (j) + (k))
#define INDEX_FLAG(i, j, k) ((NY + 4)*(NZ + 4)*(i + 2 - i_start) + (NZ + 4)*(j + 2) + (k + 2))
#define INDEX_F(i, j, k, p) (NY * NZ * NP * (i + 1 - i_start) + NZ * NP * (j) + NP * (k) + (p))

inline int mod(int x, int n)
{
    if (x < 0)
        return x + n;
    else if (x > n - 1)
        return x - n;
    else
        return x;
}

inline int min(int x, int y)
{
    if (x < y) return x;
    else return y;
}

inline double delta(int a, int b)
{
    if (a == b) return 1.0;
    else return 0.0;
}

#define FOR_DOMAIN                        \
    for (int i = i_start; i < i_end; i++) \
        for (int j = 0; j < NY; j++)      \
            for (int k = 0; k < NZ; k++)

#define UNPACK_BAGS \
    ParamBag * const params = sim->params; \
    DistributionBag * const dists = sim->dists; \
    FieldBag * const fields = sim->fields; \
    Stencil * const stencil = sim->stencil; \
    (void)params; (void)dists; (void)fields; (void)stencil;

#define UNPACK_STENCIL \
    const int NP = stencil->NP; \
    const int * restrict const cx = stencil->cx; \
    const int * restrict const cy = stencil->cy; \
    const int * restrict const cz = stencil->cz; \
    const int * restrict const p_bounceback = stencil->p_bounceback;\
    const double * restrict const wp = stencil->wp; \
    (void)NP; (void)cx; (void)cy; (void)cz; (void)wp; (void)p_bounceback;

#define UNPACK_GRID \
    const int i_start = params->i_start; \
    const int i_end = params->i_end; \
    const int NX = params->NX; \
    const int NY = params->NY; \
    const int NZ = params->NZ; \
    (void)i_start; (void)i_end; (void)NX; (void)NY; (void)NZ;

#define READ_DIST(name) const double * restrict const name = dists->name; (void)name;
#define WRITE_DIST(name) double * restrict const name = dists->name;(void)name;
#define READ_FIELD(name) const double * restrict const name = fields->name;(void)name;
#define WRITE_FIELD(name) double * restrict const name = fields->name;(void)name;
#define PARAM(name) const double name = params->name;(void)name;

#define TIME(name, functions)                                        \
    if (params->t_log == params->t)                                  \
    {                                                                \
        if (params->process_rank == 0)                               \
        {                                                            \
            printf("%-70s", name);                                   \
        }                                                            \
        MPI_Barrier(MPI_COMM_WORLD);                                 \
        start_substep = MPI_Wtime();                                 \
    }                                                                \
    functions                                                        \
    if (params->t_log == params->t)                                  \
    {                                                                \
        MPI_Barrier(MPI_COMM_WORLD);                                 \
        duration_substep = MPI_Wtime() - start_substep;              \
        if (params->process_rank == 0)                               \
        {                                                            \
            printf("[%7.4fs]\n", duration_substep);                  \
        }                                                            \
    }

#define TIME_OUTPUT(functions)                                                                \
    if ((params->process_rank == 0) && (params->t_log == params->t))                          \
    {                                                                                         \
        sprintf(output_info, "> \033[0;32mOutput to data_%d.h5...\033[0m", params->n_output); \
        printf("%-81s", output_info);                                                         \
    }                                                                                         \
    MPI_Barrier(MPI_COMM_WORLD);                                                              \
    start_substep = MPI_Wtime();                                                              \
    functions                                                                                 \
        MPI_Barrier(MPI_COMM_WORLD);                                                          \
    duration_substep = MPI_Wtime() - start_substep;                                           \
    if ((params->process_rank == 0) && (params->t_log == params->t))                                                            \
    {                                                                                         \
        printf("[%7.4fs]\n", duration_substep);                                               \
    }

#endif