#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/forcing.h"
#include "../include/wetnode.h"

void wetnode_boundary_conditions(SimulationBag *sim)
{
#ifdef WETNODE
    ParamBag *params = sim->params;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;
#else
    (void)sim;
#endif

    // Mass conservation after streaming
#if defined(LEFT_NEBB_NOSLIP) || defined(LEFT_NEBB_PRESSURE)
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_mass_conservation_streaming(0, j, k, 1, 0, 0, sim);
            }
        }
    }
#endif

#if defined(RIGHT_NEBB_NOSLIP) || defined(RIGHT_NEBB_PRESSURE)
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_mass_conservation_streaming(params->NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#if defined(BOTTOM_NEBB_NOSLIP) || defined(BOTTOM_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#if defined(TOP_NEBB_NOSLIP) || defined(TOP_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#if defined(BACK_NEBB_NOSLIP) || defined(BACK_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#if defined(FRONT_NEBB_NOSLIP) || defined(FRONT_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

    // Set velocities
#ifdef LEFT_NEBB_NOSLIP
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->glob_fields->u[INDEX_GLOB(0, j, k)] = 0.0;
                sim->glob_fields->v[INDEX_GLOB(0, j, k)] = 0.0;
                sim->glob_fields->w[INDEX_GLOB(0, j, k)] = 0.0;
                sim->comp_fields->u_comp[INDEX(0, j, k, RED)] = 0.0;
                sim->comp_fields->v_comp[INDEX(0, j, k, RED)] = 0.0;
                sim->comp_fields->w_comp[INDEX(0, j, k, RED)] = 0.0;
                sim->comp_fields->u_comp[INDEX(0, j, k, BLUE)] = 0.0;
                sim->comp_fields->v_comp[INDEX(0, j, k, BLUE)] = 0.0;
                sim->comp_fields->w_comp[INDEX(0, j, k, BLUE)] = 0.0;
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_NOSLIP
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->glob_fields->u[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
                sim->glob_fields->v[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
                sim->glob_fields->w[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
                sim->comp_fields->u_comp[INDEX(params->NX - 1, j, k, RED)] = 0.0;
                sim->comp_fields->v_comp[INDEX(params->NX - 1, j, k, RED)] = 0.0;
                sim->comp_fields->w_comp[INDEX(params->NX - 1, j, k, RED)] = 0.0;
                sim->comp_fields->u_comp[INDEX(params->NX - 1, j, k, BLUE)] = 0.0;
                sim->comp_fields->v_comp[INDEX(params->NX - 1, j, k, BLUE)] = 0.0;
                sim->comp_fields->w_comp[INDEX(params->NX - 1, j, k, BLUE)] = 0.0;
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->glob_fields->u[INDEX_GLOB(i, 0, k)] = 0.0;
            sim->glob_fields->v[INDEX_GLOB(i, 0, k)] = 0.0;
            sim->glob_fields->w[INDEX_GLOB(i, 0, k)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, 0, k, RED)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, 0, k, RED)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, 0, k, RED)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, 0, k, BLUE)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, 0, k, BLUE)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, 0, k, BLUE)] = 0.0;
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->glob_fields->u[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            sim->glob_fields->v[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            sim->glob_fields->w[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, NY - 1, k, RED)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, NY - 1, k, RED)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, NY - 1, k, RED)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, NY - 1, k, BLUE)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, NY - 1, k, BLUE)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, NY - 1, k, BLUE)] = 0.0;
        }
    }
#endif

#ifdef BACK_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->glob_fields->u[INDEX_GLOB(i, j, 0)] = 0.0;
            sim->glob_fields->v[INDEX_GLOB(i, j, 0)] = 0.0;
            sim->glob_fields->w[INDEX_GLOB(i, j, 0)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, j, 0, RED)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, j, 0, RED)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, j, 0, RED)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, j, 0, BLUE)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, j, 0, BLUE)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, j, 0, BLUE)] = 0.0;
        }
    }
#endif

#ifdef FRONT_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->glob_fields->u[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
            sim->glob_fields->v[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
            sim->glob_fields->w[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, j, NZ - 1, RED)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, j, NZ - 1, RED)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, j, NZ - 1, RED)] = 0.0;
            sim->comp_fields->u_comp[INDEX(i, j, NZ - 1, BLUE)] = 0.0;
            sim->comp_fields->v_comp[INDEX(i, j, NZ - 1, BLUE)] = 0.0;
            sim->comp_fields->w_comp[INDEX(i, j, NZ - 1, BLUE)] = 0.0;
        }
    }
#endif

    // Set densities
#ifdef LEFT_NEBB_PRESSURE
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->comp_fields->rho_comp[INDEX(0, j, k, RED)] = LEFT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                sim->comp_fields->rho_comp[INDEX(0, j, k, BLUE)] = LEFT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_PRESSURE
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                sim->comp_fields->rho_comp[INDEX(params->NX - 1, j, k, RED)] = RIGHT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                sim->comp_fields->rho_comp[INDEX(params->NX - 1, j, k, BLUE)] = RIGHT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->comp_fields->rho_comp[INDEX(i, 0, k, RED)] = BOTTOM_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            sim->comp_fields->rho_comp[INDEX(i, 0, k, BLUE)] = BOTTOM_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef TOP_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            sim->comp_fields->rho_comp[INDEX(i, NY - 1, k, RED)] = TOP_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            sim->comp_fields->rho_comp[INDEX(i, NY - 1, k, BLUE)] = TOP_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef BACK_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->comp_fields->rho_comp[INDEX(i, j, 0, RED)] = BACK_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            sim->comp_fields->rho_comp[INDEX(i, j, 0, BLUE)] = BACK_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef FRONT_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            sim->comp_fields->rho_comp[INDEX(i, j, NZ - 1, RED)] = FRONT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            sim->comp_fields->rho_comp[INDEX(i, j, NZ - 1, BLUE)] = FRONT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

    // Compute densities
#ifdef LEFT_NEBB_NOSLIP
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_density(0, j, k, 1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_NOSLIP
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_density(params->NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#ifdef BACK_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_density(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_density(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

    // Compute velocities
#ifdef LEFT_NEBB_PRESSURE
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_velocity(0, j, k, 1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_PRESSURE
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_velocity(params->NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_velocity(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#ifdef TOP_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_velocity(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#ifdef BACK_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_velocity(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_velocity(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

    // Evaluate forces
#if defined(LEFT_NEBB_NOSLIP) || defined(LEFT_NEBB_PRESSURE)
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                evaluate_force(0, j, k, sim);
            }
        }
    }
#endif

#if defined(RIGHT_NEBB_NOSLIP) || defined(RIGHT_NEBB_PRESSURE)
    if (i_end == params->NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                evaluate_force(params->NX - 1, j, k, sim);
            }
        }
    }
#endif

#if defined(BOTTOM_NEBB_NOSLIP) || defined(BOTTOM_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            evaluate_force(i, 0, k, sim);
        }
    }
#endif

#if defined(TOP_NEBB_NOSLIP) || defined(TOP_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            evaluate_force(i, NY - 1, k, sim);
        }
    }
#endif

#if defined(BACK_NEBB_NOSLIP) || defined(BACK_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            evaluate_force(i, j, 0, sim);
        }
    }
#endif

#if defined(FRONT_NEBB_NOSLIP) || defined(FRONT_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            evaluate_force(i, j, NZ - 1, sim);
        }
    }
#endif

    // Non-equilibrium bounce-back
#if defined(LEFT_NEBB_NOSLIP) || defined(LEFT_NEBB_PRESSURE)
    if (i_start == 0)
    {
        non_equilibrium_bounce_back_x(0, 1, sim);
    }
#endif

#if defined(RIGHT_NEBB_NOSLIP) || defined(RIGHT_NEBB_PRESSURE)
    if (i_end == params->NX)
    {
        non_equilibrium_bounce_back_x(params->NX - 1, -1, sim);
    }
#endif

#if defined(BOTTOM_NEBB_NOSLIP) || defined(BOTTOM_NEBB_PRESSURE)
    non_equilibrium_bounce_back_y(0, 1, sim);
#endif

#if defined(TOP_NEBB_NOSLIP) || defined(TOP_NEBB_PRESSURE)
    non_equilibrium_bounce_back_y(NY - 1, -1, sim);
#endif

#if defined(BACK_NEBB_NOSLIP) || defined(BACK_NEBB_PRESSURE)
    non_equilibrium_bounce_back_z(0, 1, sim);
#endif

#if defined(FRONT_NEBB_NOSLIP) || defined(FRONT_NEBB_PRESSURE)
    non_equilibrium_bounce_back_z(NZ - 1, -1, sim);
#endif
}

void wetnode_mass_conservation_streaming(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim)
{
    DistributionBag *dists = sim->dists;
    ParamBag *params = sim->params;
    Stencil *stencil = sim->stencil;

    double cn;

    int i_start = params->i_start;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    for (int p = 0; p < NP; p++)
    {
        cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;
        if (cn < 0)
        {
            f1[INDEX_F(i, j, k, 0, RED)] += f2[INDEX_F(i, j, k, p, RED)] - f1[INDEX_F(i, j, k, p, RED)];
            f1[INDEX_F(i, j, k, 0, BLUE)] += f2[INDEX_F(i, j, k, p, BLUE)] - f1[INDEX_F(i, j, k, p, BLUE)];
        }
    }
}

void wetnode_compute_density(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i;
    double un;
    int cn;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    int i_start = params->i_start;

    double *rho_comp = comp_fields->rho_comp;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    rho_RED_i = f2[INDEX_F(i, j, k, 0, RED)];
    rho_BLUE_i = f2[INDEX_F(i, j, k, 0, BLUE)];

    for (int p = 1; p < NP; p++)
    {
        cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

        if (cn < 0)
        {
            rho_RED_i += f1[INDEX_F(i, j, k, p, RED)] + f2[INDEX_F(i, j, k, p, RED)];
            rho_BLUE_i += f1[INDEX_F(i, j, k, p, BLUE)] + f2[INDEX_F(i, j, k, p, BLUE)];
        }
        else if (cn == 0)
        {
            rho_RED_i += f1[INDEX_F(i, j, k, p, RED)];
            rho_BLUE_i += f1[INDEX_F(i, j, k, p, BLUE)];
        }
    }

    un = u[INDEX_GLOB(i, j, k)] * (double)nx + v[INDEX_GLOB(i, j, k)] * (double)ny + w[INDEX_GLOB(i, j, k)] * (double)nz;

    rho_comp[INDEX(i, j, k, RED)] = rho_RED_i / (1.0 - un);
    rho_comp[INDEX(i, j, k, BLUE)] = rho_BLUE_i / (1.0 - un);
}

void wetnode_compute_velocity(int i, int j, int k, int nx, int ny, int nz, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double un, rho_i;
    int cn;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;

    int i_start = params->i_start;

    double *rho_comp = comp_fields->rho_comp;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *u_comp = comp_fields->u_comp;
    double *v_comp = comp_fields->v_comp;
    double *w_comp = comp_fields->w_comp;

    double *f1 = dists->f1;
    double *f2 = dists->f2;

    for (int n = 0; n < NCOMP; n++)
    {
        if (rho_comp[INDEX(i, j, k, n)] == 0.0)
        {
            u_comp[INDEX(i, j, k, n)] = 0.0;
            v_comp[INDEX(i, j, k, n)] = 0.0;
            w_comp[INDEX(i, j, k, n)] = 0.0;
            continue;
        }

        un = f2[INDEX_F(i, j, k, 0, n)];

        for (int p = 1; p < NP; p++)
        {
            cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

            if (cn < 0)
            {
                un += f1[INDEX_F(i, j, k, p, n)] + f2[INDEX_F(i, j, k, p, n)];
            }
            else if (cn == 0)
            {
                un += f1[INDEX_F(i, j, k, p, n)];
            }
        }

        un = 1.0 - un / rho_comp[INDEX(i, j, k, n)];
        u_comp[INDEX(i, j, k, n)] = un * (double)nx;
        v_comp[INDEX(i, j, k, n)] = un * (double)ny;
        w_comp[INDEX(i, j, k, n)] = un * (double)nz;
    }
    rho_i = rho_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)];
    u[INDEX_GLOB(i, j, k)] = (rho_comp[INDEX(i, j, k, RED)] * u_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)] * u_comp[INDEX(i, j, k, BLUE)])/ rho_i;
    v[INDEX_GLOB(i, j, k)] = (rho_comp[INDEX(i, j, k, RED)] * v_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)] * v_comp[INDEX(i, j, k, BLUE)])/ rho_i;
    w[INDEX_GLOB(i, j, k)] = (rho_comp[INDEX(i, j, k, RED)] * w_comp[INDEX(i, j, k, RED)] + rho_comp[INDEX(i, j, k, BLUE)] * w_comp[INDEX(i, j, k, BLUE)])/ rho_i;
}

void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, vc_i, wc_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *v_comp = comp_fields->v_comp;
    double *w_comp = comp_fields->w_comp;

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int j = 0; j < NY; j++)
    {
        for (int k = 0; k < NZ; k++)
        {
            for (int n = 0; n < NCOMP; n++)
            {
                rho_i = rho_comp[INDEX(i, j, k, n)];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1[INDEX_F(i, j, k, p, n)] = 0.0;
                    }
                    continue;
                }

                u_i = u[INDEX_GLOB(i, j, k)];
                v_i = v[INDEX_GLOB(i, j, k)];
                w_i = w[INDEX_GLOB(i, j, k)];

                vc_i = v_comp[INDEX_GLOB(i, j, k)];
                wc_i = w_comp[INDEX_GLOB(i, j, k)];

                Fx_i = Fx[INDEX(i, j, k, n)];
                Fy_i = Fy[INDEX(i, j, k, n)];
                Fz_i = Fz[INDEX(i, j, k, n)];

                Nx = 0.5 * Fx_i / C_norm;
                Ny = 0.5 * Fy_i - rho_i * vc_i;
                Nz = 0.5 * Fz_i - rho_i * wc_i;

                for (int p = 1; p < NP; p++)
                {
                    if (cx[p] == 0)
                    {
                        Ny += f1[INDEX_F(i, j, k, p, n)] * (double)cy[p];
                        Nz += f1[INDEX_F(i, j, k, p, n)] * (double)cz[p];
                    }

                    else if (cx[p] * nx > 0)
                    {
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                        Ny += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cy[p];
                        Nz += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cz[p];
                    }
                }

                Ny /= C_par;
                Nz /= C_par;

                for (int p = 1; p < NP; p++)
                {
                    if (cx[p] * nx > 0)
                    {
                        c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                        f1[INDEX_F(i, j, k, p, n)] = f1[INDEX_F(i, j, k, p_bounceback[p], n)] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                    }
                }

                f1[INDEX_F(i, j, k, 0, n)] += 0.5 * Fx_i * nx;
            }
        }
    }
}

void non_equilibrium_bounce_back_y(int j, int ny, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, uc_i, wc_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *u_comp = comp_fields->u_comp;
    double *w_comp = comp_fields->w_comp;

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            for (int n = 0; n < NCOMP; n++)
            {
                rho_i = rho_comp[INDEX(i, j, k, n)];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1[INDEX_F(i, j, k, p, n)] = 0.0;
                    }
                    continue;
                }

                u_i = u[INDEX_GLOB(i, j, k)];
                v_i = v[INDEX_GLOB(i, j, k)];
                w_i = w[INDEX_GLOB(i, j, k)];

                uc_i = u_comp[INDEX(i, j, k, n)];
                wc_i = w_comp[INDEX(i, j, k, n)];

                Fx_i = Fx[INDEX(i, j, k, n)];
                Fy_i = Fy[INDEX(i, j, k, n)];
                Fz_i = Fz[INDEX(i, j, k, n)];

                Nx = 0.5 * Fx_i - rho_i * uc_i;
                Ny = 0.5 * Fy_i / C_norm;
                Nz = 0.5 * Fz_i - rho_i * wc_i;

                for (int p = 1; p < NP; p++)
                {
                    if (cy[p] == 0)
                    {
                        Nx += f1[INDEX_F(i, j, k, p, n)] * (double)cx[p];
                        Nz += f1[INDEX_F(i, j, k, p, n)] * (double)cz[p];
                    }

                    else if (cy[p] * ny > 0)
                    {
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                        Nx += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cx[p];
                        Nz += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cz[p];
                    }
                }

                Nx /= C_par;
                Nz /= C_par;

                for (int p = 1; p < NP; p++)
                {
                    if (cy[p] * ny > 0)
                    {
                        c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                        f1[INDEX_F(i, j, k, p, n)] = f1[INDEX_F(i, j, k, p_bounceback[p], n)] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                    }
                }

                f1[INDEX_F(i, j, k, 0, n)] += 0.5 * Fy_i * ny;
            }
        }
    }
}

void non_equilibrium_bounce_back_z(int k, int nz, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_i, Fx_i, Fy_i, Fz_i, Nx, Ny, Nz;
    double u_i, v_i, w_i, uc_i, vc_i, uc, c_i;

    int NY = params->NY;
    int NZ = params->NZ;
    int NP = stencil->NP;

    double cs2 = stencil->cs2;
    int *cx = stencil->cx;
    int *cy = stencil->cy;
    int *cz = stencil->cz;
    double *wp = stencil->wp;

    int *p_bounceback = stencil->p_bounceback;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    double *u_comp = comp_fields->u_comp;
    double *v_comp = comp_fields->v_comp;

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int n = 0; n < NCOMP; n++)
            {
                rho_i = rho_comp[INDEX(i, j, k, n)];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1[INDEX_F(i, j, k, p, n)] = 0.0;
                    }
                    continue;
                }

                u_i = u[INDEX_GLOB(i, j, k)];
                v_i = v[INDEX_GLOB(i, j, k)];
                w_i = w[INDEX_GLOB(i, j, k)];

                uc_i = u_comp[INDEX(i, j, k, n)];
                vc_i = v_comp[INDEX(i, j, k, n)];

                Fx_i = Fx[INDEX(i, j, k, n)];
                Fy_i = Fy[INDEX(i, j, k, n)];
                Fz_i = Fz[INDEX(i, j, k, n)];

                Nx = 0.5 * Fx_i - rho_i * uc_i;
                Ny = 0.5 * Fy_i - rho_i * vc_i;
                Nz = 0.5 * Fz_i / C_norm;

                for (int p = 1; p < NP; p++)
                {
                    if (cz[p] == 0)
                    {
                        Nx += f1[INDEX_F(i, j, k, p, n)] * (double)cx[p];
                        Ny += f1[INDEX_F(i, j, k, p, n)] * (double)cy[p];
                    }

                    else if (cz[p] * nz > 0)
                    {
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                        Nx += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cx[p];
                        Ny += 2.0 * wp[p] * rho_i * uc / cs2 * (double)cy[p];
                    }
                }

                Nx /= C_par;
                Ny /= C_par;

                for (int p = 1; p < NP; p++)
                {
                    if (cz[p] * nz > 0)
                    {
                        c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                        uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];
                        f1[INDEX_F(i, j, k, p, n)] = f1[INDEX_F(i, j, k, p_bounceback[p], n)] + 2.0 * wp[p] * rho_i * uc / cs2 - (double)cx[p] / c_i * Nx - (double)cy[p] / c_i * Ny - (double)cz[p] / c_i * Nz;
                    }
                }

                f1[INDEX_F(i, j, k, 0, n)] += 0.5 * Fz_i * nz;
            }
        }
    }
}