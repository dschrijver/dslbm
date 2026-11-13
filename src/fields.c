#include <math.h>
#include <string.h>

#include "../include/datatypes.h"
#include "../include/communicate.h"

void compute_gradient(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(gradient_input)
    WRITE_FIELD(gradient_output_x)
    WRITE_FIELD(gradient_output_y)
    WRITE_FIELD(gradient_output_z)

    communicate_field(sim->fields->gradient_input, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        double gradient_x_i = 0.0;
        double gradient_y_i = 0.0;
        double gradient_z_i = 0.0;
        for (int p = 0; p < NP; p++)
        {
            // TODO: boundary conditions
            const double intput_local = gradient_input[INDEX(i + cx[p], j + cy[p], k + cz[p])];

            gradient_x_i += wp[p] * intput_local * (double)cx[p];
            gradient_y_i += wp[p] * intput_local * (double)cy[p];
            gradient_z_i += wp[p] * intput_local * (double)cz[p];
        }
        gradient_output_x[idx0] = 3.0 * gradient_x_i;
        gradient_output_y[idx0] = 3.0 * gradient_y_i;
        gradient_output_z[idx0] = 3.0 * gradient_z_i;
    }
}

void compute_laplacian(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_FIELD(laplacian_input)
    WRITE_FIELD(laplacian_output)

    communicate_field(sim->fields->laplacian_input, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        double laplacian_i = 0.0;
        for (int p = 0; p < NP; p++)
        {
            // TODO: boundary conditions
            const double intput_local = laplacian_input[INDEX(i + cx[p], j + cy[p], k + cz[p])];

            laplacian_i += wp[p] * intput_local;
        }
        laplacian_output[idx0] = 6.0 * (laplacian_i - laplacian_input[idx0]);
    }
}

void compute_chemical_potentials(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)
    READ_FIELD(rho)
    READ_FIELD(laplacian_output)

    WRITE_FIELD(laplacian_input)
    WRITE_FIELD(mu_a)
    WRITE_FIELD(mu_p)
    

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        laplacian_input[idx0] = phi[idx0];
    }

    compute_laplacian(sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        const double rho_i = rho[idx0];
        const double phi_i = phi[idx0];

        const double term = phi_i - 0.5;

        const double mu_rho_i = 1.0/3.0*(log(rho_i) + 1.0) + params->a*term*term;
        const double mu_phi_i = 2.0*params->a*term - params->kappa*laplacian_output[idx0];

        mu_a[idx0] = mu_rho_i - phi_i/rho_i*mu_phi_i;
        mu_p[idx0] = mu_rho_i + (1.0 - phi_i)/rho_i*mu_phi_i;
    }
}

void extract_moments(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    READ_DIST(f1_a)
    READ_DIST(f1_p)

    WRITE_FIELD(rho)
    WRITE_FIELD(phi)
    WRITE_FIELD(rho_a)
    WRITE_FIELD(rho_p)
    WRITE_FIELD(u_a)
    WRITE_FIELD(v_a)
    WRITE_FIELD(w_a)
    WRITE_FIELD(u_p)
    WRITE_FIELD(v_p)
    WRITE_FIELD(w_p)

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        double rho_a_i = 0.0;
        double rho_p_i = 0.0;
        double u_a_i = 0.0;
        double v_a_i = 0.0;
        double w_a_i = 0.0;
        double u_p_i = 0.0;
        double v_p_i = 0.0;
        double w_p_i = 0.0;
        
        for (int p = 0; p < NP; p++)
        {
            const int idxf0 = INDEX_F(i, j, k, p);
            const double f_a = f1_a[idxf0];
            const double f_p = f1_p[idxf0];

            rho_a_i += f_a;
            rho_p_i += f_p;

            u_a_i += f_a * (double)cx[p];
            v_a_i += f_a * (double)cy[p];
            w_a_i += f_a * (double)cz[p];
            u_p_i += f_p * (double)cx[p];
            v_p_i += f_p * (double)cy[p];
            w_p_i += f_p * (double)cz[p];
        }

        rho[idx0] = rho_a_i + rho_p_i;
        phi[idx0] = rho_a_i / (rho_a_i + rho_p_i);

        rho_a[idx0] = rho_a_i;
        rho_p[idx0] = rho_p_i;

        u_a[idx0] = u_a_i / rho_a_i;
        v_a[idx0] = v_a_i / rho_a_i;
        w_a[idx0] = w_a_i / rho_a_i;
        u_p[idx0] = u_p_i / rho_p_i;
        v_p[idx0] = v_p_i / rho_p_i;
        w_p[idx0] = w_p_i / rho_p_i;
    }
}

void compute_final_velocity(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)
    READ_FIELD(rho_a)
    READ_FIELD(rho_p)
    READ_FIELD(Fx_a)
    READ_FIELD(Fy_a)
    READ_FIELD(Fz_a)
    READ_FIELD(Fx_p)
    READ_FIELD(Fy_p)
    READ_FIELD(Fz_p)

    WRITE_FIELD(u_a)
    WRITE_FIELD(v_a)
    WRITE_FIELD(w_a)
    WRITE_FIELD(u_p)
    WRITE_FIELD(v_p)
    WRITE_FIELD(w_p)

    // Preliminary velocity shift
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double phi_i = phi[idx0];
        const double rho_a_i = rho_a[idx0];
        const double rho_p_i = rho_p[idx0];
        const double Fx_a_i = Fx_a[idx0];
        const double Fy_a_i = Fy_a[idx0];
        const double Fz_a_i = Fz_a[idx0];
        const double Fx_p_i = Fx_p[idx0];
        const double Fy_p_i = Fy_p[idx0];
        const double Fz_p_i = Fz_p[idx0];

        // Compute preliminary velociites
        const double u_a_i = u_a[idx0] + 0.5 / rho_a_i * Fx_a_i;
        const double v_a_i = v_a[idx0] + 0.5 / rho_a_i * Fy_a_i;
        const double w_a_i = w_a[idx0] + 0.5 / rho_a_i * Fz_a_i;
        const double u_p_i = u_p[idx0] + 0.5 / rho_p_i * Fx_p_i;
        const double v_p_i = v_p[idx0] + 0.5 / rho_p_i * Fy_p_i;
        const double w_p_i = w_p[idx0] + 0.5 / rho_p_i * Fz_p_i;

        const double A = params->gamma * phi_i * (1.0 - phi_i);
        const double B = 0.5 / rho_a_i * A;
        const double C = 0.5 / rho_p_i * A;
        const double D = 1.0 / (1.0 + B + C);

        u_a[idx0] = D * ((1.0 + C)*u_a_i + B*u_p_i);
        v_a[idx0] = D * ((1.0 + C)*v_a_i + B*v_p_i);
        w_a[idx0] = D * ((1.0 + C)*w_a_i + B*w_p_i);
        u_p[idx0] = D * (C*u_a_i + (1.0 + B)*u_p_i);
        v_p[idx0] = D * (C*v_a_i + (1.0 + B)*v_p_i);
        w_p[idx0] = D * (C*w_a_i + (1.0 + B)*w_p_i);
    }
}