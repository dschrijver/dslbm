#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 500000;
    params->NSTORE = 5000;
    params->NLOG = 500;
    params->NX = 160;
    params->NY = 160;
    params->NZ = 1;

    // Initial denities
    params->rho_0_RED   = 1000.0;
    params->rho_0_BLUE  = 1.0;

    // Kinematic viscosities
    params->mu_RED  = 0.0005*params->rho_0_RED;
    params->mu_BLUE = 0.0005*20.0*params->rho_0_BLUE;

    // Speeds of sound squared
    params->cs2_BLUE = 1.0/3.0;
    params->cs2_RED = 1.0/3.0/(params->rho_0_RED / params->rho_0_BLUE);

    // Surface tension
    params->sigma = 1e-4;

    // Recoloring parameter
    params->beta = 0.7;

    // Gravitational accelerations
    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;

    // Constant body forces
    params->Fb_x = 0.0;
    params->Fb_y = 0.0;
    params->Fb_z = 0.0;
}

#endif