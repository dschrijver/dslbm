#ifndef PARAMS_H
#define PARAMS_H

#include "include/datatypes.h"

static inline void set_params(ParamBag *params)
{
    // General parameters
    params->NTIME = 2000000;
    params->NSTORE = 20000;
    params->NLOG = 2000;
    params->NX = 200;
    params->NY = 1;
    params->NZ = 1;

    // Initial denities
    params->rho_0_BLUE  = 1.0;
    params->rho_0_RED   = params->rho_0_BLUE * 1000.0;

    // Kinematic viscosities
    double nu_RED  = 0.5;
    double nu_BLUE = nu_RED / 100.0;

    // Kinematic viscosities
    params->mu_RED  = params->rho_0_RED * nu_RED;
    params->mu_BLUE = params->rho_0_BLUE * nu_BLUE;

    // Speeds of sound squared
    params->cs2_BLUE = 1.0/3.0;
    params->cs2_RED = 1.0/3.0 / (params->rho_0_RED / params->rho_0_BLUE);

    // Surface tension
    params->sigma = 0.0;

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