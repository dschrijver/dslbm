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

    // Dynamic viscosities
    double mu_RED;
    double mu_BLUE;

    // Speeds of sound squared
    double cs2_RED;
    double cs2_BLUE;

    // Surface tension
    double sigma;

    // Recoloring parameter
    double beta;

    // Gravitational accelerations
    double gx;
    double gy;
    double gz;

    // Constant body forces
    double Fb_x;
    double Fb_y;
    double Fb_z;

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
    hid_t memspace;
    hid_t dcpl_id;
    hid_t dxpl_id;
} ParamBag;

typedef struct DistributionBag
{
    double *f1_RED;
    double *f2_RED;

    double *f1_BLUE;
    double *f2_BLUE;

    double *meq;
    double *t_star;
    double *m_star;
    double *f_star;

    // Communication
    double *send_buffer;
    double *recv_buffer;
} DistributionBag;

typedef struct FieldBag
{
    // Densities and pressure
    double *rho;

    double *rho_RED;
    double *rho_BLUE;

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
    double *nx;
    double *ny;
    double *nz;

    // Wen's correction
    double *Qx;
    double *Qy;
    double *Qz;

    // Forces
    double *Fx;
    double *Fy;
    double *Fz;

    // Boundaries
    int *flag;

    // Communication
    double *send_buffer;
    double *recv_buffer;
} FieldBag;

typedef struct Stencil
{
    // Stencil
    int NP;
    int *cx;
    int *cy;
    int *cz;
    double *wp;
    int *p_bounceback;
} Stencil;

typedef struct SimulationBag
{
    struct ParamBag *params;
    struct DistributionBag *dists;
    struct FieldBag *fields;
    struct Stencil *stencil;
} SimulationBag;

#endif