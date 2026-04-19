#ifndef DATATYPES_H
#define DATATYPES_H

#include <mpi.h>
#include <hdf5.h>

#include "utils.h"

typedef struct ParamBag
{
    // Simulation variables
    int t;       
    int t_output;
    int t_log;  
    int n_output; 

    // General parameters
    int NTIME; 
    int NSTORE; 
    int NLOG; 
    int NX;    
    int NY;  
    int NZ;  

    // Initial denities
    double rho_0_RED;
    double rho_0_BLUE;

    // Kinematic viscosities
    double nu_RED;
    double nu_BLUE;

    // Speed of sound parameters
    double alpha_RED;
    double alpha_BLUE;

    // Surface tension
    double sigma;

    // Recoloring parameter
    double beta;

    // Gravitational accelerations
    double gx;
    double gy;
    double gz;

    // MPI
    int number_of_processes;  ///< Stores number of processes
    int process_rank;         ///< Process rank
    int process_coords[3];    ///< Coords of processor in virtual MPI topology
    int process_neighbors[2]; ///< Left and right processor neighbors in virtual MPI topology.
    int i_start;              ///< Starting index of slab owned by current processor.
    int i_end;                ///< Ending index of slab owned by current processor.
    int NX_proc;              ///< Number of nodes owned by current processor, NX_proc = i_end - i_start.
    MPI_Comm comm_xslices;    ///< Communicator of slab decomposition.

    // HDF5
    hid_t fapl_id;
    hid_t scalar_space;
    hid_t filespace;
    hid_t memspace_glob;
    hid_t memspace_comp;
    hid_t dcpl_id;
    hid_t dxpl_id;
} ParamBag;

typedef struct DistributionBag
{
    double *f1;
    double *f2;
    double *raw;
    double *k_pert;
    double *k_star;
    double *feq;

    // Communication
    double *send_buffer;
    double *recv_buffer;
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

    // Color-Gradient
    double *rho_N;
    double *G_norm;
    double *Gx;
    double *Gy;
    double *Gz;

    // Boundaries
    int *flag;

    // Communication
    double *send_buffer;
    double *recv_buffer;
} GlobalFieldBag;

typedef struct ComponentFieldBag
{
    // Densities of components
    double *rho_comp;

    // Forces on components
    double *Fx;
    double *Fy;
    double *Fz;

    // Component velocities
    double *u_comp;
    double *v_comp;
    double *w_comp;
} ComponentFieldBag;

typedef struct Stencil
{
    double cs2[2];
    double tau[2];

    // Stencil
    int NP;
    int *cx;
    int *cy;
    int *cz;
    double *wp;
    int *p_bounceback;

    // Color-Gradient
    double zeta;
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