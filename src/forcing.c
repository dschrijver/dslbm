#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"
#include "../include/forcing.h"

void evaluate_forces(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    FOR_DOMAIN
    {
        evaluate_force(i, j, k, sim);
    }
}

void evaluate_force(int i, int j, int k, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(rho)

    WRITE_FIELD(Fx)
    WRITE_FIELD(Fy)
    WRITE_FIELD(Fz)

    PARAM(gx)
    PARAM(gy)
    PARAM(gz)

    PARAM(Fb_x)
    PARAM(Fb_y)
    PARAM(Fb_z)

    const double rho_i = rho[INDEX(i, j, k)];

    const int idx = INDEX(i, j, k);

    Fx[idx] = rho_i * gx + Fb_x;
    Fy[idx] = rho_i * gy + Fb_y;
    Fz[idx] = rho_i * gz + Fb_z;

    evaluate_surface_force(i, j, k, sim);
}

void evaluate_surface_force(int i, int j, int k, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(nx)
    READ_FIELD(ny)
    READ_FIELD(nz)
    READ_FIELD(Gx)
    READ_FIELD(Gy)
    READ_FIELD(Gz)

    WRITE_FIELD(Fx)
    WRITE_FIELD(Fy)
    WRITE_FIELD(Fz)

    PARAM(sigma)

    int *c_vec[3] = {cx, cy, cz};

    const int * restrict const flag = fields->flag;

    double n_i[3];
    double n_local[3];
    double dn[3][3];

    for (int alpha = 0; alpha < 3; alpha++)
    {
        for (int beta = 0; beta < 3; beta++)
        {
            dn[alpha][beta] = 0.0;
        }
    }

    for (int p = 0; p < NP; p++)
    {
        int ic = i + cx[p];
        int jc = j + cy[p];
        int kc = k + cz[p];

#ifdef YPERIODIC
        jc = mod(jc, NY);
#endif
#ifdef ZPERIODIC
        kc = mod(kc, NZ);
#endif

        if (flag[INDEX(ic, jc, kc)] > 0)
        {
            n_local[0] = extrapolate_wall_n(ic, jc, kc, 0, sim);
            n_local[1] = extrapolate_wall_n(ic, jc, kc, 1, sim);
            n_local[2] = extrapolate_wall_n(ic, jc, kc, 2, sim);
            goto skip;
        }

        n_local[0] = nx[INDEX(ic, jc, kc)];
        n_local[1] = ny[INDEX(ic, jc, kc)];
        n_local[2] = nz[INDEX(ic, jc, kc)];

    skip:

        for (int alpha = 0; alpha < 3; alpha++)
        {
            for (int beta = 0; beta < 3; beta++)
            {
                dn[alpha][beta] += 3.0 * wp[p] * n_local[beta] * (double)c_vec[alpha][p];
            }
        }
    }

    const int idx = INDEX(i, j, k);

    n_i[0] = nx[idx];
    n_i[1] = ny[idx];
    n_i[2] = nz[idx];

    double kappa = 0.0;
    for (int alpha = 0; alpha < 3; alpha++)
    {
        kappa += -dn[alpha][alpha];
        for (int beta = 0; beta < 3; beta++)
        {
            kappa += n_i[alpha] * n_i[beta] * dn[beta][alpha];
        }
    }

    Fx[idx] += 0.5 * sigma * kappa * Gx[idx];
    Fy[idx] += 0.5 * sigma * kappa * Gy[idx];
    Fz[idx] += 0.5 * sigma * kappa * Gz[idx];
}

double extrapolate_wall_n(const int i, const int j, const int k, const int alpha, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(nx)
    READ_FIELD(ny)
    READ_FIELD(nz)

    double *n_vec[3] = {nx, ny, nz};

    double sum_n = 0.0;
    double sum_wp = 0.0;

    for (int p = 1; p < NP; p++)
    {
        int ic = i + cx[p];
        int jc = j + cy[p];
        int kc = k + cz[p];

#ifndef XPERIODIC
        if ((ic < 0) || (ic > params->NX - 1))
            continue;
#endif
#ifndef YPERIODIC
        if ((jc < 0) || (jc > NY - 1))
            continue;
#else
        jc = mod(jc, NY);
#endif
#ifndef ZPERIODIC
        if ((kc < 0) || (kc > NZ - 1))
            continue;
#else
        kc = mod(kc, NZ);
#endif

        sum_n += wp[p] * n_vec[alpha][INDEX(ic, jc, kc)];
        sum_wp += wp[p];
    }

    if (sum_wp == 0.0)
        return 0.0;
    else
        return sum_n / sum_wp;
}