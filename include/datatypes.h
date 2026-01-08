#ifndef DATATYPES_H
#define DATATYPES_H

#include <mpi.h>

#include "utils.h"

typedef struct ParamBag
{
    // Simulation variables
    int t;        ///< Current time in simulation
    int t_output; ///< Time at which to output next h5 files.
    int t_log;    ///< Time at which to log progress.
    int n_output; ///< Index of next h5 output file.

    // General parameters
    int NTIME;  ///< Final time of simulation.
    int NSTORE; ///< Store data in data_*.h5 files every NSTORE steps.
    int NLOG;   ///< Log progress every NSTORE steps.
    int NX;     ///< Number of fluid points in x (i) direction. Note: Max_number_of_processors <= NX/2.
    int NY;     ///< Number of fluid points in y (j) direction.
    int NZ;     ///< Number of fluid points in z (k) direction.

    // Relaxation times
    double tau_RED;  ///< Viscous relaxation time of RED fluid.
    double tau_BLUE; ///< Viscous relaxation time of BLUE fluid.

    // MRT parameters
    double s_rho_RED;  ///< MRT parameter.
    double s_rho_BLUE; ///< MRT parameter.
    double s_j_RED;    ///< MRT parameter.
    double s_j_BLUE;   ///< MRT parameter.
    double s_e_RED;    ///< MRT parameter.
    double s_e_BLUE;   ///< MRT parameter.
    double s_eps_RED;  ///< MRT parameter.
    double s_eps_BLUE; ///< MRT parameter.
    double s_q_RED;    ///< MRT parameter.
    double s_q_BLUE;   ///< MRT parameter.
    double s_pi_RED;   ///< MRT parameter.
    double s_pi_BLUE;  ///< MRT parameter.
    double s_m_RED;    ///< MRT parameter.
    double s_m_BLUE;   ///< MRT parameter.

    double rho_0_RED;
    double rho_0_BLUE;

    double sigma;
    double beta;
    double alpha_RED;
    double alpha_BLUE;

    double gx;
    double gy;
    double gz;

    double G_SC;

    // MPI
    int number_of_processes;  ///< Stores number of processes
    int process_rank;         ///< Process rank
    int process_coords[3];    ///< Coords of processor in virtual MPI topology
    int process_neighbors[2]; ///< Left and right processor neighbors in virtual MPI topology.
    int i_start;              ///< Starting index of slab owned by current processor.
    int i_end;                ///< Ending index of slab owned by current processor.
    int NX_proc;              ///< Number of nodes owned by current processor, NX_proc = i_end - i_start.
    MPI_Comm comm_xslices;    ///< Communicator of slab decomposition.
} ParamBag;

typedef struct DistributionBag
{
    double *f1;
    double *f2;
} DistributionBag;

typedef struct GlobalFieldBag
{
    // Densities and pressure
    double *rho;
    double *pressure;

    // Velocities
    double *u;
    double *v;
    double *w;
} GlobalFieldBag;

typedef struct ComponentFieldBag
{
    // Densities of components
    double *rho_comp;

    // Component velocities
    double *u_comp;
    double *v_comp;
    double *w_comp;

    // Forces on components
    double *Fx;
    double *Fy;
    double *Fz;
} ComponentFieldBag;

typedef struct Stencil
{
    double cs2; // 1/3

    // Stencil
    int NP;
    int *cx;
    int *cy;
    int *cz;
    double *wp;
    int *p_bounceback;

    // MRT
    double *M;
    double *M_inv;
    double *omega[2];

    // Color-Gradient
    double zeta;
    double *B;
    double *phi_eq[2];
    double *wp_D3Q41;
} Stencil;

typedef struct SimulationBag
{
    struct ParamBag *params;
    struct DistributionBag *dists;
    struct GlobalFieldBag *glob_fields;
    struct ComponentFieldBag *comp_fields;
    struct Stencil *stencil;
} SimulationBag;

#endif