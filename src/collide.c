#include <math.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/collide.h"

void compute_equilibrium(double rho, double u, double v, double w, double *feq, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_STENCIL

    const double u2 = u*u + v*v + w*w;

    for (int p = 0; p < NP; p++)
    {
        const double uc = u*(double)cx[p] + v*(double)cy[p] + w*(double)cz[p];

        feq[p] = wp[p]*rho * (1.0 + 3.0*uc + 4.5*uc*uc - 1.5*u2);
    }
}

void compute_source(double Fx, double Fy, double Fz, double u, double v, double w, double *S, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_STENCIL

    for (int p = 0; p < NP; p++)
    {
        const double uc = u*(double)cx[p] + v*(double)cy[p] + w*(double)cz[p];

        const double sx = (3.0 * ((double)cx[p] - u) + 9.0*uc*(double)cx[p]) * Fx;
        const double sy = (3.0 * ((double)cy[p] - v) + 9.0*uc*(double)cy[p]) * Fy;
        const double sz = (3.0 * ((double)cz[p] - w) + 9.0*uc*(double)cz[p]) * Fz;

        S[p] = wp[p] * (sx + sy + sz);
    }
}

void collide(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho_a)
    READ_FIELD(rho_p)
    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)
    READ_FIELD(u_p)
    READ_FIELD(v_p)
    READ_FIELD(w_p)
    READ_FIELD(Fx_a)
    READ_FIELD(Fy_a)
    READ_FIELD(Fz_a)
    READ_FIELD(Fx_p)
    READ_FIELD(Fy_p)
    READ_FIELD(Fz_p)

    READ_DIST(f1_a)
    READ_DIST(f1_p)

    WRITE_DIST(f2_a)
    WRITE_DIST(f2_p)

    PARAM(eta_a)
    PARAM(eta_p)

    double feq_a[NP], feq_p[NP], S_a[NP], S_p[NP];

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double rho_a_i = rho_a[idx0];
        const double rho_p_i = rho_p[idx0];
        const double u_a_i = u_a[idx0];
        const double v_a_i = v_a[idx0];
        const double w_a_i = w_a[idx0];
        const double u_p_i = u_p[idx0];
        const double v_p_i = v_p[idx0];
        const double w_p_i = w_p[idx0];
        const double Fx_a_i = Fx_a[idx0];
        const double Fy_a_i = Fy_a[idx0];
        const double Fz_a_i = Fz_a[idx0];
        const double Fx_p_i = Fx_p[idx0];
        const double Fy_p_i = Fy_p[idx0];
        const double Fz_p_i = Fz_p[idx0];

        const double tau_a = 3.0*eta_a/rho_a_i + 0.5;
        const double tau_p = 3.0*eta_p/rho_p_i + 0.5;

        compute_equilibrium(rho_a_i, u_a_i, v_a_i, w_a_i, feq_a, sim);
        compute_equilibrium(rho_p_i, u_p_i, v_p_i, w_p_i, feq_p, sim);
        compute_source(Fx_a_i, Fy_a_i, Fz_a_i, u_a_i, v_a_i, w_a_i, S_a, sim);
        compute_source(Fx_p_i, Fy_p_i, Fz_p_i, u_p_i, v_p_i, w_p_i, S_p, sim);

        for (int p = 0; p < NP; p++)
        {
            f2_a[INDEX_F(i, j, k, p)] = (1.0 - 1.0 / tau_a)*f1_a[INDEX_F(i, j, k, p)] + 1.0 / tau_a * feq_a[p] + (1.0 - 0.5 / tau_a)*S_a[p];
            f2_p[INDEX_F(i, j, k, p)] = (1.0 - 1.0 / tau_p)*f1_p[INDEX_F(i, j, k, p)] + 1.0 / tau_p * feq_p[p] + (1.0 - 0.5 / tau_p)*S_p[p];
        }
    }
}