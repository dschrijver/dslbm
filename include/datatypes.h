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
    double *f1;
    double *f2;
    double *g1;
    double *g2;

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
    double *rho;
    double *phi;

    double *mu;
    double *pressure;

    // Velocities
    double *u;
    double *v;
    double *w;

    // Forces
    double *Fx;
    double *Fy;
    double *Fz;

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