#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 1000000;
    params->NSTORE = 10000;
    params->NLOG = 1000;
    params->NX = 100;
    params->NY = 2;
    params->NZ = 2;

    // Initial denities
    params->rho_0_RED   = 8.0;
    params->rho_0_BLUE  = 0.008;

    // Kinematic viscosities
    params->mu_RED  = 0.05826;
    params->mu_BLUE = params->mu_RED * 1.0/40.0;

    // Speeds of sound squared
    double alpha_RED = 0.9992;
    params->cs2_RED = 9.0/19.0 * (1.0 - alpha_RED);
    params->cs2_BLUE = params->cs2_RED * (params->rho_0_RED / params->rho_0_BLUE);

    // Surface tension
    params->sigma = 0.002/4.5;

    // Recoloring parameter
    params->beta = 0.7;

    // Gravitational accelerations
    params->gx = 0.0;
    params->gy = 0.0;
    params->gz = 0.0;

    // Constant body forces
    params->Fb_x = 0.0;
    params->Fb_y = 1.5e-8;
    params->Fb_z = 0.0;
}

#endif
