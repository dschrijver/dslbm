#ifndef DATATYPES_H
#define DATATYPES_H

#include <mpi.h>
#include <hdf5.h>

#include "utils.h"

typedef double matrix[3][3];

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
    double rho_0_a;
    double rho_0_p;

    // Kinematic viscosities
    double nu_a;
    double nu_p;

    double kappa;
    double a;
    double gamma;

    double S_nem;
    double C;
    double K_LC;
    double lambda;
    double zeta;
    double Gamma_Q;

    int N_active_substeps;

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
    MPI_Datatype MPI_MATRIX;

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
    double *f1_a;
    double *f2_a;

    double *f1_p;
    double *f2_p;

    // Distribution buffers
    double *send_buffer_x;
    double *recv_buffer_x;
    double *send_buffer_y;
    double *recv_buffer_y;
    double *send_buffer_z;
    double *recv_buffer_z;
} DistributionBag;

typedef struct FieldBag
{
    // Densities and composition
    double *rho;
    double *phi;
    double *rho_a;
    double *rho_p;

    // Velocities
    double *u;
    double *v;
    double *w;
    double *u_a;
    double *v_a;
    double *w_a;
    double *u_p;
    double *v_p;
    double *w_p;

    // Forces
    double *Fx_a;
    double *Fy_a;
    double *Fz_a;
    double *Fx_p;
    double *Fy_p;
    double *Fz_p;

    // Chemical potentials
    double *mu_a;
    double *mu_p;

    // Helper fields
    double *gradient_input;
    double *gradient_output_x;
    double *gradient_output_y;
    double *gradient_output_z;
    double *laplacian_input;
    double *laplacian_output;

    // Scalar buffers
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

typedef struct MatrixBag
{
    matrix *Q;
    matrix *Q_new;
    matrix *H;
    matrix *sigma_nem;
    matrix *nabla_u;
    matrix *E;
    matrix *Omega;

    matrix *matrix_gradient_output_x;
    matrix *matrix_gradient_output_y;
    matrix *matrix_gradient_output_z;
    double *matrix_divergence_output_x;
    double *matrix_divergence_output_y;
    double *matrix_divergence_output_z;

    // Matrix buffers
    matrix *send_buffer_x;
    matrix *recv_buffer_x;
    matrix *send_buffer_y;
    matrix *recv_buffer_y;
    matrix *send_buffer_z;
    matrix *recv_buffer_z;
} MatrixBag;

typedef struct SimulationBag
{
    struct ParamBag *params;
    struct DistributionBag *dists;
    struct FieldBag *fields;
    struct Stencil *stencil;
    struct MatrixBag *matrices;
} SimulationBag;

#endif