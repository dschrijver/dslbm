#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/forcing.h"
#include "../include/wetnode.h"

void wetnode_boundary_conditions(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;

    int NY = params->NY;
    int NZ = params->NZ;

    int i_start = params->i_start;
    int i_end = params->i_end;

    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;

    // Mass conservation after streaming
#ifdef LEFT_NEBB_NOSLIP
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

#ifdef RIGHT_NEBB_NOSLIP
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

#ifdef BOTTOM_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#ifdef BACK_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_NOSLIP
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
                u[INDEX_GLOB(0, j, k)] = 0.0;
                v[INDEX_GLOB(0, j, k)] = 0.0;
                w[INDEX_GLOB(0, j, k)] = 0.0;
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
                u[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
                v[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
                w[INDEX_GLOB(params->NX - 1, j, k)] = 0.0;
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX_GLOB(i, 0, k)] = 0.0;
            v[INDEX_GLOB(i, 0, k)] = 0.0;
            w[INDEX_GLOB(i, 0, k)] = 0.0;
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            v[INDEX_GLOB(i, NY - 1, k)] = 0.0;
            w[INDEX_GLOB(i, NY - 1, k)] = 0.0;
        }
    }
#endif

#ifdef BACK_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[INDEX_GLOB(i, j, 0)] = 0.0;
            v[INDEX_GLOB(i, j, 0)] = 0.0;
            w[INDEX_GLOB(i, j, 0)] = 0.0;
        }
    }
#endif

#ifdef FRONT_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
            v[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
            w[INDEX_GLOB(i, j, NZ - 1)] = 0.0;
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

    // Evaluate forces
#ifdef LEFT_NEBB_NOSLIP
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

#ifdef RIGHT_NEBB_NOSLIP
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

#ifdef BOTTOM_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            evaluate_force(i, 0, k, sim);
        }
    }
#endif

#ifdef TOP_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            evaluate_force(i, NY - 1, k, sim);
        }
    }
#endif

#ifdef BACK_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            evaluate_force(i, j, 0, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_NOSLIP
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            evaluate_force(i, j, NZ - 1, sim);
        }
    }
#endif

    // Non-equilibrium bounce-back
#ifdef LEFT_NEBB_NOSLIP
    if (i_start == 0)
    {
        non_equilibrium_bounce_back_x(0, 1, sim);
    }
#endif

#ifdef RIGHT_NEBB_NOSLIP
    if (i_end == params->NX)
    {
        non_equilibrium_bounce_back_x(params->NX - 1, -1, sim);
    }
#endif

#ifdef BOTTOM_NEBB_NOSLIP
    non_equilibrium_bounce_back_y(0, 1, sim);
#endif

#ifdef TOP_NEBB_NOSLIP
    non_equilibrium_bounce_back_y(NY - 1, -1, sim);
#endif

#ifdef BACK_NEBB_NOSLIP
    non_equilibrium_bounce_back_z(0, 1, sim);
#endif

#ifdef FRONT_NEBB_NOSLIP
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
    int cn, un;

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

    un = u[INDEX_GLOB(i, j, k)] * nx + v[INDEX_GLOB(i, j, k)] * ny + w[INDEX_GLOB(i, j, k)] * nz;

    rho_comp[INDEX(i, j, k, RED)] = rho_RED_i / (1.0 - un);
    rho_comp[INDEX(i, j, k, BLUE)] = rho_BLUE_i / (1.0 - un);
}

void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;
    DistributionBag *dists = sim->dists;
    Stencil *stencil = sim->stencil;

    double rho_RED_i, rho_BLUE_i, Fx_RED_i, Fy_RED_i, Fz_RED_i, Fx_BLUE_i, Fy_BLUE_i, Fz_BLUE_i;
    double Nx_RED, Ny_RED, Nz_RED, Nx_BLUE, Ny_BLUE, Nz_BLUE;
    double u_i, v_i, w_i, uc, c_i;

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

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int j = 0; j < NY; j++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
            rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

            u_i = u[INDEX_GLOB(i, j, k)];
            v_i = v[INDEX_GLOB(i, j, k)];
            w_i = w[INDEX_GLOB(i, j, k)];

            Fx_RED_i = Fx[INDEX(i, j, k, RED)];
            Fy_RED_i = Fy[INDEX(i, j, k, RED)];
            Fz_RED_i = Fz[INDEX(i, j, k, RED)];

            Fx_BLUE_i = Fx[INDEX(i, j, k, BLUE)];
            Fy_BLUE_i = Fy[INDEX(i, j, k, BLUE)];
            Fz_BLUE_i = Fz[INDEX(i, j, k, BLUE)];

            Nx_RED = 0.5 * Fx_RED_i / C_norm;
            Ny_RED = 0.5 * Fy_RED_i - rho_RED_i * v_i;
            Nz_RED = 0.5 * Fz_RED_i - rho_RED_i * w_i;

            Nx_BLUE = 0.5 * Fx_BLUE_i / C_norm;
            Ny_BLUE = 0.5 * Fy_BLUE_i - rho_BLUE_i * v_i;
            Nz_BLUE = 0.5 * Fz_BLUE_i - rho_BLUE_i * w_i;

            for (int p = 1; p < NP; p++)
            {
                if (cx[p] == 0)
                {
                    Ny_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cy[p];
                    Nz_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cz[p];

                    Ny_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cy[p];
                    Nz_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cz[p];
                }

                else if (cx[p] * nx > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Ny_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cy[p];
                    Nz_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cz[p];

                    Ny_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cy[p];
                    Nz_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cz[p];
                }
            }

            Ny_RED /= C_par;
            Nz_RED /= C_par;

            Ny_BLUE /= C_par;
            Nz_BLUE /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cx[p] * nx > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    f1[INDEX_F(i, j, k, p, RED)] = f1[INDEX_F(i, j, k, p_bounceback[p], RED)] - (double)cx[p] / c_i * Nx_RED - (double)cy[p] / c_i * Ny_RED - (double)cz[p] / c_i * Nz_RED;
                    f1[INDEX_F(i, j, k, p, BLUE)] = f1[INDEX_F(i, j, k, p_bounceback[p], BLUE)] - (double)cx[p] / c_i * Nx_BLUE - (double)cy[p] / c_i * Ny_BLUE - (double)cz[p] / c_i * Nz_BLUE;
                }
            }

            f1[INDEX_F(i, j, k, 0, RED)] += 0.5 * Fx_RED_i * nx;
            f1[INDEX_F(i, j, k, 0, BLUE)] += 0.5 * Fx_BLUE_i * nx;
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

    double rho_RED_i, rho_BLUE_i, Fx_RED_i, Fy_RED_i, Fz_RED_i, Fx_BLUE_i, Fy_BLUE_i, Fz_BLUE_i;
    double Nx_RED, Ny_RED, Nz_RED, Nx_BLUE, Ny_BLUE, Nz_BLUE;
    double u_i, v_i, w_i, uc, c_i;

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

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
            rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

            u_i = u[INDEX_GLOB(i, j, k)];
            v_i = v[INDEX_GLOB(i, j, k)];
            w_i = w[INDEX_GLOB(i, j, k)];

            Fx_RED_i = Fx[INDEX(i, j, k, RED)];
            Fy_RED_i = Fy[INDEX(i, j, k, RED)];
            Fz_RED_i = Fz[INDEX(i, j, k, RED)];

            Fx_BLUE_i = Fx[INDEX(i, j, k, BLUE)];
            Fy_BLUE_i = Fy[INDEX(i, j, k, BLUE)];
            Fz_BLUE_i = Fz[INDEX(i, j, k, BLUE)];

            Nx_RED = 0.5 * Fx_RED_i - rho_RED_i * u_i;
            Ny_RED = 0.5 * Fy_RED_i / C_norm;
            Nz_RED = 0.5 * Fz_RED_i - rho_RED_i * w_i;

            Nx_BLUE = 0.5 * Fx_BLUE_i - rho_BLUE_i * u_i;
            Ny_BLUE = 0.5 * Fy_BLUE_i / C_norm;
            Nz_BLUE = 0.5 * Fz_BLUE_i - rho_BLUE_i * w_i;

            for (int p = 1; p < NP; p++)
            {
                if (cy[p] == 0)
                {
                    Nx_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cx[p];
                    Nz_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cz[p];

                    Nx_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cx[p];
                    Nz_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cz[p];
                }

                else if (cy[p] * ny > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Nx_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cx[p];
                    Nz_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cz[p];

                    Nx_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cx[p];
                    Nz_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cz[p];
                }
            }

            Nx_RED /= C_par;
            Nz_RED /= C_par;

            Nx_BLUE /= C_par;
            Nz_BLUE /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cy[p] * ny > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    f1[INDEX_F(i, j, k, p, RED)] = f1[INDEX_F(i, j, k, p_bounceback[p], RED)] - (double)cx[p] / c_i * Nx_RED - (double)cy[p] / c_i * Ny_RED - (double)cz[p] / c_i * Nz_RED;
                    f1[INDEX_F(i, j, k, p, BLUE)] = f1[INDEX_F(i, j, k, p_bounceback[p], BLUE)] - (double)cx[p] / c_i * Nx_BLUE - (double)cy[p] / c_i * Ny_BLUE - (double)cz[p] / c_i * Nz_BLUE;
                }
            }

            f1[INDEX_F(i, j, k, 0, RED)] += 0.5 * Fy_RED_i * ny;
            f1[INDEX_F(i, j, k, 0, BLUE)] += 0.5 * Fy_BLUE_i * ny;
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

    double rho_RED_i, rho_BLUE_i, Fx_RED_i, Fy_RED_i, Fz_RED_i, Fx_BLUE_i, Fy_BLUE_i, Fz_BLUE_i;
    double Nx_RED, Ny_RED, Nz_RED, Nx_BLUE, Ny_BLUE, Nz_BLUE;
    double u_i, v_i, w_i, uc, c_i;

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

    double *f1 = dists->f1;

    double C_norm = 1.0 + 2.0 * sqrt(2.0);
    double C_par = sqrt(2.0);

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            rho_RED_i = rho_comp[INDEX(i, j, k, RED)];
            rho_BLUE_i = rho_comp[INDEX(i, j, k, BLUE)];

            u_i = u[INDEX_GLOB(i, j, k)];
            v_i = v[INDEX_GLOB(i, j, k)];
            w_i = w[INDEX_GLOB(i, j, k)];

            Fx_RED_i = Fx[INDEX(i, j, k, RED)];
            Fy_RED_i = Fy[INDEX(i, j, k, RED)];
            Fz_RED_i = Fz[INDEX(i, j, k, RED)];

            Fx_BLUE_i = Fx[INDEX(i, j, k, BLUE)];
            Fy_BLUE_i = Fy[INDEX(i, j, k, BLUE)];
            Fz_BLUE_i = Fz[INDEX(i, j, k, BLUE)];

            Nx_RED = 0.5 * Fx_RED_i - rho_RED_i * u_i;
            Ny_RED = 0.5 * Fy_RED_i - rho_RED_i * v_i;
            Nz_RED = 0.5 * Fz_RED_i / C_norm;

            Nx_BLUE = 0.5 * Fx_BLUE_i - rho_BLUE_i * u_i;
            Ny_BLUE = 0.5 * Fy_BLUE_i - rho_BLUE_i * v_i;
            Nz_BLUE = 0.5 * Fz_BLUE_i / C_norm;

            for (int p = 1; p < NP; p++)
            {
                if (cz[p] == 0)
                {
                    Nx_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cx[p];
                    Ny_RED += f1[INDEX_F(i, j, k, p, RED)] * (double)cy[p];

                    Nx_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cx[p];
                    Ny_BLUE += f1[INDEX_F(i, j, k, p, BLUE)] * (double)cy[p];
                }

                else if (cz[p] * nz > 0)
                {
                    uc = u_i * (double)cx[p] + v_i * (double)cy[p] + w_i * (double)cz[p];

                    Nx_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cx[p];
                    Ny_RED += 2.0 * wp[p] * rho_RED_i * uc / cs2 * (double)cy[p];

                    Nx_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cx[p];
                    Ny_BLUE += 2.0 * wp[p] * rho_BLUE_i * uc / cs2 * (double)cy[p];
                }
            }

            Nx_RED /= C_par;
            Ny_RED /= C_par;

            Nx_BLUE /= C_par;
            Ny_BLUE /= C_par;

            for (int p = 1; p < NP; p++)
            {
                if (cz[p] * nz > 0)
                {
                    c_i = sqrt((double)cx[p] * (double)cx[p] + (double)cy[p] * (double)cy[p] + (double)cz[p] * (double)cz[p]);
                    f1[INDEX_F(i, j, k, p, RED)] = f1[INDEX_F(i, j, k, p_bounceback[p], RED)] - (double)cx[p] / c_i * Nx_RED - (double)cy[p] / c_i * Ny_RED - (double)cz[p] / c_i * Nz_RED;
                    f1[INDEX_F(i, j, k, p, BLUE)] = f1[INDEX_F(i, j, k, p_bounceback[p], BLUE)] - (double)cx[p] / c_i * Nx_BLUE - (double)cy[p] / c_i * Ny_BLUE - (double)cz[p] / c_i * Nz_BLUE;
                }
            }

            f1[INDEX_F(i, j, k, 0, RED)] += 0.5 * Fz_RED_i * nz;
            f1[INDEX_F(i, j, k, 0, BLUE)] += 0.5 * Fz_BLUE_i * nz;
        }
    }
}