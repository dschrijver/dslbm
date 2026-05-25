#include "../definitions.h"
#include "../include/collide.h"
#include "../include/wetnode.h"

void wetnode_macroscopic_fields(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    WRITE_FIELD(rho_RED)
    WRITE_FIELD(rho_BLUE)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)

#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
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

#if defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_NEBB_PRESSURE)
    if (i_end == NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_mass_conservation_streaming(NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#if defined(TOP_NEBB_VELOCITY) || defined(TOP_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_mass_conservation_streaming(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#if defined(FRONT_NEBB_VELOCITY) || defined(FRONT_NEBB_PRESSURE)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_mass_conservation_streaming(i, j, NZ - 1, 0, 0, -1, sim);
        }
    }
#endif

#if defined(LEFT_NEBB_VELOCITY)
    if (i_start == 0)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                u[INDEX(0, j, k)] = LEFT_U_VELOCITY;
                v[INDEX(0, j, k)] = LEFT_V_VELOCITY;
                w[INDEX(0, j, k)] = LEFT_W_VELOCITY;
            }
        }
    }
#endif

#if defined(RIGHT_NEBB_VELOCITY)
    if (i_end == NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                u[INDEX(NX - 1, j, k)] = RIGHT_U_VELOCITY;
                v[INDEX(NX - 1, j, k)] = RIGHT_V_VELOCITY;
                w[INDEX(NX - 1, j, k)] = RIGHT_W_VELOCITY;
            }
        }
    }
#endif

#if defined(BOTTOM_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX(i, 0, k)] = BOTTOM_U_VELOCITY;
            v[INDEX(i, 0, k)] = BOTTOM_V_VELOCITY;
            w[INDEX(i, 0, k)] = BOTTOM_W_VELOCITY;
        }
    }
#endif

#if defined(TOP_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            u[INDEX(i, NY - 1, k)] = TOP_U_VELOCITY;
            v[INDEX(i, NY - 1, k)] = TOP_V_VELOCITY;
            w[INDEX(i, NY - 1, k)] = TOP_W_VELOCITY;
        }
    }
#endif

#if defined(BACK_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[INDEX(i, j, 0)] = BACK_U_VELOCITY;
            v[INDEX(i, j, 0)] = BACK_V_VELOCITY;
            w[INDEX(i, j, 0)] = BACK_W_VELOCITY;
        }
    }
#endif

#if defined(FRONT_NEBB_VELOCITY)
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            u[INDEX(i, j, NZ - 1)] = FRONT_U_VELOCITY;
            v[INDEX(i, j, NZ - 1)] = FRONT_V_VELOCITY;
            w[INDEX(i, j, NZ - 1)] = FRONT_W_VELOCITY;
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
                rho_RED[INDEX(0, j, k)] = LEFT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                rho_BLUE[INDEX(0, j, k)] = LEFT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
            }
        }
    }
#endif

#ifdef RIGHT_NEBB_PRESSURE
    if (i_end == NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                rho_RED[INDEX(NX - 1, j, k)] = RIGHT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
                rho_BLUE[INDEX(NX - 1, j, k)] = RIGHT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_RED[INDEX(i, 0, k)] = BOTTOM_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            rho_BLUE[INDEX(i, 0, k)] = BOTTOM_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef TOP_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            rho_RED[INDEX(i, NY - 1, k)] = TOP_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            rho_BLUE[INDEX(i, NY - 1, k)] = TOP_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef BACK_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            rho_RED[INDEX(i, j, 0)] = BACK_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            rho_BLUE[INDEX(i, j, 0)] = BACK_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

#ifdef FRONT_NEBB_PRESSURE
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            rho_RED[INDEX(i, j, NZ - 1)] = FRONT_PRESSURE_RED / (sim->stencil->zeta * (1.0 - params->alpha_RED));
            rho_BLUE[INDEX(i, j, NZ - 1)] = FRONT_PRESSURE_BLUE / (sim->stencil->zeta * (1.0 - params->alpha_BLUE));
        }
    }
#endif

    // Compute densities
#ifdef LEFT_NEBB_VELOCITY
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

#ifdef RIGHT_NEBB_VELOCITY
    if (i_end == NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_density(NX - 1, j, k, -1, 0, 0, sim);
            }
        }
    }
#endif

#ifdef BOTTOM_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, 0, k, 0, 1, 0, sim);
        }
    }
#endif

#ifdef TOP_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            wetnode_compute_density(i, NY - 1, k, 0, -1, 0, sim);
        }
    }
#endif

#ifdef BACK_NEBB_VELOCITY
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            wetnode_compute_density(i, j, 0, 0, 0, 1, sim);
        }
    }
#endif

#ifdef FRONT_NEBB_VELOCITY
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
    if (i_end == NX)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                wetnode_compute_velocity(NX - 1, j, k, -1, 0, 0, sim);
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
}

void wetnode_distributions(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    // Non-equilibrium bounce-back
#if defined(LEFT_NEBB_VELOCITY) || defined(LEFT_NEBB_PRESSURE)
    if (i_start == 0)
    {
        non_equilibrium_bounce_back_x(0, 1, sim);
    }
#endif

#if defined(RIGHT_NEBB_VELOCITY) || defined(RIGHT_NEBB_PRESSURE)
    if (i_end == NX)
    {
        non_equilibrium_bounce_back_x(NX - 1, -1, sim);
    }
#endif

#if defined(BOTTOM_NEBB_VELOCITY) || defined(BOTTOM_NEBB_PRESSURE)
    non_equilibrium_bounce_back_y(0, 1, sim);
#endif

#if defined(TOP_NEBB_VELOCITY) || defined(TOP_NEBB_PRESSURE)
    non_equilibrium_bounce_back_y(NY - 1, -1, sim);
#endif

#if defined(BACK_NEBB_VELOCITY) || defined(BACK_NEBB_PRESSURE)
    non_equilibrium_bounce_back_z(0, 1, sim);
#endif

#if defined(FRONT_NEBB_VELOCITY) || defined(FRONT_NEBB_PRESSURE)
    non_equilibrium_bounce_back_z(NZ - 1, -1, sim);
#endif
}

void wetnode_mass_conservation_streaming(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f2_RED)
    READ_DIST(f2_BLUE)

    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    for (int p = 0; p < NP; p++)
    {
        const int cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;
        if (cn < 0)
        {
            f1_RED[INDEX_F(i, j, k, 0)] += f2_RED[INDEX_F(i, j, k, p)] - f1_RED[INDEX_F(i, j, k, p)];
            f1_BLUE[INDEX_F(i, j, k, 0)] += f2_BLUE[INDEX_F(i, j, k, p)] - f1_BLUE[INDEX_F(i, j, k, p)];
        }
    }
}

void wetnode_compute_density(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)
    READ_DIST(f1_RED)
    READ_DIST(f1_BLUE)
    READ_DIST(f2_RED)
    READ_DIST(f2_BLUE)

    WRITE_FIELD(rho_RED)
    WRITE_FIELD(rho_BLUE)

    double rho_RED_i = f2_RED[INDEX_F(i, j, k, 0)];
    double rho_BLUE_i = f2_BLUE[INDEX_F(i, j, k, 0)];

    for (int p = 1; p < NP; p++)
    {
        const int cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

        if (cn < 0)
        {
            rho_RED_i += f1_RED[INDEX_F(i, j, k, p)] + f2_RED[INDEX_F(i, j, k, p)];
            rho_BLUE_i += f1_BLUE[INDEX_F(i, j, k, p)] + f2_BLUE[INDEX_F(i, j, k, p)];
        }
        else if (cn == 0)
        {
            rho_RED_i += f1_RED[INDEX_F(i, j, k, p)];
            rho_BLUE_i += f1_BLUE[INDEX_F(i, j, k, p)];
        }
    }

    const int idx = INDEX(i, j, k);

    const double un = u[idx] * (double)nx + v[idx] * (double)ny + w[idx] * (double)nz;

    rho_RED[idx] = rho_RED_i / (1.0 - un);
    rho_BLUE[idx] = rho_BLUE_i / (1.0 - un);
}

void wetnode_compute_velocity(const int i, const int j, const int k, const int nx, const int ny, const int nz, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_DIST(f1_RED)
    READ_DIST(f1_BLUE)
    READ_DIST(f2_RED)
    READ_DIST(f2_BLUE)

    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)

    const int idx = INDEX(i, j, k);
    const double rho_RED_i = rho_RED[idx];
    const double rho_BLUE_i = rho_BLUE[idx];

    if (rho_BLUE_i < 1e-15)
    {
        double un = f2_RED[INDEX_F(i, j, k, 0)];
        for (int p = 1; p < NP; p++)
        {
            const int cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

            if (cn < 0)
            {
                un += f1_RED[INDEX_F(i, j, k, p)] + f2_RED[INDEX_F(i, j, k, p)];
            }
            else if (cn == 0)
            {
                un += f1_RED[INDEX_F(i, j, k, p)];
            }
        }
        un = 1.0 - un / rho_RED_i;
        u[idx] = un * (double)nx;
        v[idx] = un * (double)ny;
        w[idx] = un * (double)nz;
    }   
    else
    {
        double un = f2_BLUE[INDEX_F(i, j, k, 0)];
        for (int p = 1; p < NP; p++)
        {
            const int cn = cx[p] * nx + cy[p] * ny + cz[p] * nz;

            if (cn < 0)
            {
                un += f1_BLUE[INDEX_F(i, j, k, p)] + f2_BLUE[INDEX_F(i, j, k, p)];
            }
            else if (cn == 0)
            {
                un += f1_BLUE[INDEX_F(i, j, k, p)];
            }
        }
        un = 1.0 - un / rho_BLUE_i;
        u[idx] = un * (double)nx;
        v[idx] = un * (double)ny;
        w[idx] = un * (double)nz;
    }
}

void non_equilibrium_bounce_back_x(int i, int nx, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)
    
    WRITE_DIST(f_star)
    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    double *rho_comp[2] = {rho_RED, rho_BLUE};
    double *f1_comp[2] = {f1_RED, f1_BLUE};

    for (int j = 0; j < NY; j++)
    {
        for (int k = 0; k < NZ; k++)
        {
            const int idx = INDEX(i, j, k);
            const double pressure_i = pressure[idx];
            const double rho_tot_i = rho[idx];
            const double u_i = u[idx];
            const double v_i = v[idx];
            const double w_i = w[idx];
            const double Fx_tot_i = Fx[idx];
            const double Fy_tot_i = Fy[idx];
            const double Fz_tot_i = Fz[idx];

            for (int n = 0; n < 2; n++)
            {
                const double rho_i = rho_comp[n][idx];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1_comp[n][INDEX_F(i, j, k, p)] = 0.0;
                    }
                    continue;
                }

                compute_equilibrium_comp(rho_i, rho_tot_i, u_i, v_i, w_i, pressure_i, f_star, sim);

                const double Fx_i = rho_i/rho_tot_i*Fx_tot_i;
                const double Fy_i = rho_i/rho_tot_i*Fy_tot_i;
                const double Fz_i = rho_i/rho_tot_i*Fz_tot_i;

                double Nx = 3.0 * Fx_i;
                double Ny = -rho_i*v_i + 0.5 * Fy_i;
                double Nz = -rho_i*w_i + 0.5 * Fz_i;
                for (int p = 1; p < NP; p++)
                {
                    if (cx[p] == 0)
                    {
                        const int idxf = INDEX_F(i, j, k, p);
                        Ny += f1_comp[n][idxf] * (double)cy[p];
                        Nz += f1_comp[n][idxf] * (double)cz[p];
                    }
                    else if (cx[p] * nx > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        Ny += dfeq*(double)cy[p];
                        Nz += dfeq*(double)cz[p];
                    }
                }
                Ny *= 18.0;
                Nz *= 18.0;

                for (int p = 1; p < NP; p++)
                {
                    if (cx[p] * nx > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        f1_comp[n][INDEX_F(i, j, k, p)] = f1_comp[n][INDEX_F(i, j, k, p_bounceback[p])] + dfeq - wp[p] * ((double)cx[p] * Nx + (double)cy[p] * Ny + (double)cz[p] * Nz);
                    }
                }

                f1_comp[n][INDEX_F(i, j, k, 0)] += 0.5 * Fx_i * nx;
            }

            u[idx] = u_i - 0.5*Fx_tot_i/rho_tot_i;
            v[idx] = v_i - 0.5*Fy_tot_i/rho_tot_i;
            w[idx] = w_i - 0.5*Fz_tot_i/rho_tot_i;
        }
    }
}

void non_equilibrium_bounce_back_y(int j, int ny, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)
    
    WRITE_DIST(f_star)
    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    double *rho_comp[2] = {rho_RED, rho_BLUE};
    double *f1_comp[2] = {f1_RED, f1_BLUE};

    for (int i = i_start; i < i_end; i++)
    {
        for (int k = 0; k < NZ; k++)
        {
            const int idx = INDEX(i, j, k);
            const double pressure_i = pressure[idx];
            const double rho_tot_i = rho[idx];
            const double u_i = u[idx];
            const double v_i = v[idx];
            const double w_i = w[idx];
            const double Fx_tot_i = Fx[idx];
            const double Fy_tot_i = Fy[idx];
            const double Fz_tot_i = Fz[idx];
            for (int n = 0; n < 2; n++)
            {
                const double rho_i = rho_comp[n][idx];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1_comp[n][INDEX_F(i, j, k, p)] = 0.0;
                    }
                    continue;
                }

                compute_equilibrium_comp(rho_i, rho_tot_i, u_i, v_i, w_i, pressure_i, f_star, sim);

                const double Fx_i = rho_i/rho_tot_i*Fx_tot_i;
                const double Fy_i = rho_i/rho_tot_i*Fy_tot_i;
                const double Fz_i = rho_i/rho_tot_i*Fz_tot_i;

                double Nx = -rho_i*u_i + 0.5 * Fx_i;
                double Ny = 3.0 * Fy_i;
                double Nz = -rho_i*w_i + 0.5 * Fz_i;

                for (int p = 1; p < NP; p++)
                {
                    if (cy[p] == 0)
                    {
                        const int idxf = INDEX_F(i, j, k, p);
                        Nx += f1_comp[n][idxf] * (double)cx[p];
                        Nz += f1_comp[n][idxf] * (double)cz[p];
                    }
                    else if (cy[p] * ny > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        Nx += dfeq*(double)cx[p];
                        Nz += dfeq*(double)cz[p];
                    }
                }
                Nx *= 18.0;
                Nz *= 18.0;

                for (int p = 1; p < NP; p++)
                {
                    if (cy[p] * ny > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        f1_comp[n][INDEX_F(i, j, k, p)] = f1_comp[n][INDEX_F(i, j, k, p_bounceback[p])] + dfeq - wp[p] * ((double)cx[p] * Nx + (double)cy[p] * Ny + (double)cz[p] * Nz);
                    }
                }

                f1_comp[n][INDEX_F(i, j, k, 0)] += 0.5 * Fy_i * ny;
            }

            u[idx] = u_i - 0.5*Fx_tot_i/rho_tot_i;
            v[idx] = v_i - 0.5*Fy_tot_i/rho_tot_i;
            w[idx] = w_i - 0.5*Fz_tot_i/rho_tot_i;
        }
    }
}

void non_equilibrium_bounce_back_z(int k, int nz, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    WRITE_FIELD(u)
    WRITE_FIELD(v)
    WRITE_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)
    
    WRITE_DIST(f_star)
    WRITE_DIST(f1_RED)
    WRITE_DIST(f1_BLUE)

    double *rho_comp[2] = {rho_RED, rho_BLUE};
    double *f1_comp[2] = {f1_RED, f1_BLUE};

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            const int idx = INDEX(i, j, k);
            const double pressure_i = pressure[idx];
            const double rho_tot_i = rho[idx];
            const double u_i = u[idx];
            const double v_i = v[idx];
            const double w_i = w[idx];
            const double Fx_tot_i = Fx[idx];
            const double Fy_tot_i = Fy[idx];
            const double Fz_tot_i = Fz[idx];
            for (int n = 0; n < 2; n++)
            {
                const double rho_i = rho_comp[n][idx];

                if (rho_i == 0.0)
                {
                    for (int p = 0; p < NP; p++)
                    {
                        f1_comp[n][INDEX_F(i, j, k, p)] = 0.0;
                    }
                    continue;
                }

                compute_equilibrium_comp(rho_i, rho_tot_i, u_i, v_i, w_i, pressure_i, f_star, sim);

                const double Fx_i = rho_i/rho_tot_i*Fx_tot_i;
                const double Fy_i = rho_i/rho_tot_i*Fy_tot_i;
                const double Fz_i = rho_i/rho_tot_i*Fz_tot_i;

                double Nx = -rho_i*u_i + 0.5 * Fx_i;
                double Ny = -rho_i*v_i + 0.5 * Fy_i;
                double Nz = 3.0 * Fz_i;

                for (int p = 1; p < NP; p++)
                {
                    if (cz[p] == 0)
                    {
                        const int idxf = INDEX_F(i, j, k, p);
                        Nx += f1_comp[n][idxf] * (double)cx[p];
                        Ny += f1_comp[n][idxf] * (double)cy[p];
                    }
                    else if (cz[p] * nz > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        Nx += dfeq*(double)cx[p];
                        Ny += dfeq*(double)cy[p];
                    }
                }
                Nx *= 18.0;
                Ny *= 18.0;

                for (int p = 1; p < NP; p++)
                {
                    if (cz[p] * nz > 0)
                    {
                        const double dfeq = f_star[p]-f_star[p_bounceback[p]];
                        f1_comp[n][INDEX_F(i, j, k, p)] = f1_comp[n][INDEX_F(i, j, k, p_bounceback[p])] + dfeq - wp[p] * ((double)cx[p] * Nx + (double)cy[p] * Ny + (double)cz[p] * Nz);
                    }
                }

                f1_comp[n][INDEX_F(i, j, k, 0)] += 0.5 * Fz_i * nz;
            }

            u[idx] = u_i - 0.5*Fx_tot_i/rho_tot_i;
            v[idx] = v_i - 0.5*Fy_tot_i/rho_tot_i;
            w[idx] = w_i - 0.5*Fz_tot_i/rho_tot_i;
        }
    }
}