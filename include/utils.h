#ifndef UTILS_H
#define UTILS_H

enum components
{
    RED,
    BLUE,
    NCOMP
};

enum flags
{
    FLUID,
    BOUNDARY,
    WETNODE
};

#define DS_PI 3.14159265358979323846
#define INDEX_GLOB(i, j, k) (NY * NZ * (i + 2 - i_start) + NZ * (j) + (k))
#define INDEX(i, j, k, n) (NY * NZ * NCOMP * (i + 2 - i_start) + NZ * NCOMP * (j) + NCOMP * (k) + (n))
#define INDEX_F(i, j, k, p, n) (NY * NZ * NP * NCOMP * (i + 1 - i_start) + NZ * NP * NCOMP * (j) + NP * NCOMP * (k) + NP * (n) + (p))
#define INDEX_FLAG(i, j, k) ((NY + 4)*(NZ + 4)*(i + 2 - i_start) + (NZ + 4)*(j + 2) + (k + 2))
#define f(p) f1[INDEX_F(i, j, k, p, n)]

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

#define TIME(name, functions)                                        \
    if ((params->process_rank == 0) && (params->t_log == params->t)) \
    {                                                                \
        printf("%-70s", name);                                       \
    }                                                                \
    MPI_Barrier(MPI_COMM_WORLD);                                     \
    start_substep = MPI_Wtime();                                     \
    functions                                                        \
        MPI_Barrier(MPI_COMM_WORLD);                                 \
    duration_substep = MPI_Wtime() - start_substep;                  \
    if ((params->process_rank == 0) && (params->t_log == params->t)) \
    {                                                                \
        printf("[%7.4fs]\n", duration_substep);                      \
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