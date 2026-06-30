#include <math.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/forcing.h"
#include "../include/fields.h"

void extract_moments(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f1_RED)
    READ_DIST(f1_BLUE)

    WRITE_FIELD(rho)
    WRITE_FIELD(rho_RED)
    WRITE_FIELD(rho_BLUE)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    WRITE_FIELD(rho_N)

    PARAM(rho_0_RED)
    PARAM(rho_0_BLUE)

    FOR_DOMAIN
    {
        double rho_RED_i = 0.0;
        double rho_BLUE_i = 0.0;
        double u_i = 0.0;
        double v_i = 0.0;
        double w_i = 0.0;

        for (int p = 0; p < NP; p++)
        {
            const int idfx = INDEX_F(i, j, k, p);

            const double f1_RED_i = f1_RED[idfx];
            const double f1_BLUE_i = f1_BLUE[idfx];
            const double f1_i = f1_RED_i + f1_BLUE_i;

            rho_RED_i += f1_RED_i;
            rho_BLUE_i += f1_BLUE_i;
            u_i += f1_i * (double)cx[p];
            v_i += f1_i * (double)cy[p];
            w_i += f1_i * (double)cz[p];
        }

        const int idx = INDEX(i, j, k);

        const double rho_i = rho_RED_i + rho_BLUE_i;
        rho[idx] = rho_i;
        rho_N[idx] = (rho_RED_i / rho_0_RED - rho_BLUE_i / rho_0_BLUE) / (rho_RED_i / rho_0_RED + rho_BLUE_i / rho_0_BLUE);

        rho_RED[idx] = rho_RED_i;
        rho_BLUE[idx] = rho_BLUE_i;

        u[idx] = u_i / rho_i;
        v[idx] = v_i / rho_i;
        w[idx] = w_i / rho_i;
    }
}

void update_final_velocity(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    
    READ_FIELD(rho)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)

    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)

    FOR_DOMAIN
    {
        const int idx = INDEX(i, j, k);

        const double rho_i = rho[idx];

        const double Fx_i = Fx[idx];
        const double Fy_i = Fy[idx];
        const double Fz_i = Fz[idx];

        u[idx] += 0.5 * Fx_i / rho_i;
        v[idx] += 0.5 * Fy_i / rho_i;
        w[idx] += 0.5 * Fz_i / rho_i;
    }
}

void evaluate_pressure(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)

    WRITE_FIELD(pressure)

    PARAM(cs2_RED)
    PARAM(cs2_BLUE)
    PARAM(rho_0_RED)
    PARAM(rho_0_BLUE)

    const double p_star_RED = rho_0_RED * cs2_RED - rho_0_BLUE * cs2_BLUE;

    FOR_DOMAIN
    {
        const int idx = INDEX(i, j, k);

        const double rho_RED_i = rho_RED[idx];
        const double rho_BLUE_i = rho_BLUE[idx];

        const double b = rho_RED_i * cs2_RED + rho_BLUE_i * cs2_BLUE - p_star_RED;

        pressure[idx] = 0.5 * (b + sqrt(b * b + 4.0 * rho_BLUE_i * cs2_BLUE * p_star_RED));
    }
}