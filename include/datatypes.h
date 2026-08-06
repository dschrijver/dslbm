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
    int number_of_processes; 
    int process_rank;     
    int process_coords[3]; 
    int process_left;
    int process_right;
    int process_bottom;
    int process_top;
    int process_back;
    int process_front;
    int i_start;              
    int i_end;                
    int NX_proc;         
    int j_start;           
    int j_end;              
    int NY_proc;            
    int k_start;             
    int k_end;              
    int NZ_proc;            
    MPI_Comm comm_cart;    

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
    double *send_buffer_x;
    double *recv_buffer_x;
    double *send_buffer_y;
    double *recv_buffer_y;
    double *send_buffer_z;
    double *recv_buffer_z;
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
    double *send_buffer_x;
    double *recv_buffer_x;
    double *send_buffer_y;
    double *recv_buffer_y;
    double *send_buffer_z;
    double *recv_buffer_z;
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