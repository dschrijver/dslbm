#include "../include/datatypes.h"
#include "../include/fields.h"
#include "../include/active.h"

void compute_preliminary_forces(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(mu_a)
    READ_FIELD(mu_p)
    READ_FIELD(gradient_output_x)
    READ_FIELD(gradient_output_y)
    READ_FIELD(gradient_output_z)
    READ_FIELD(rho_a)
    READ_FIELD(rho_p)
    READ_FIELD(phi)
    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)
    READ_FIELD(u_p)
    READ_FIELD(v_p)
    READ_FIELD(w_p)
    
    WRITE_FIELD(gradient_input)
    WRITE_FIELD(Fx_a)
    WRITE_FIELD(Fy_a)
    WRITE_FIELD(Fz_a)
    WRITE_FIELD(Fx_p)
    WRITE_FIELD(Fy_p)
    WRITE_FIELD(Fz_p)

    // Compute gradient of mu_a
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        gradient_input[idx0] = mu_a[idx0];
    }
    compute_gradient(sim);

    // Add to active forces
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        const double rho_a_i = rho_a[idx0];

        Fx_a[idx0] = -rho_a_i * gradient_output_x[idx0];
        Fy_a[idx0] = -rho_a_i * gradient_output_y[idx0];
        Fz_a[idx0] = -rho_a_i * gradient_output_z[idx0];
    }

    // Compute gradient of mu_p
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        gradient_input[idx0] = mu_p[idx0];
    }
    compute_gradient(sim);

    // Add to passive forces
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        const double rho_p_i = rho_p[idx0];

        Fx_p[idx0] = -rho_p_i * gradient_output_x[idx0];
        Fy_p[idx0] = -rho_p_i * gradient_output_y[idx0];
        Fz_p[idx0] = -rho_p_i * gradient_output_z[idx0];
    }

    // Compute gradient of rho_a
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        gradient_input[idx0] = rho_a[idx0];
    }
    compute_gradient(sim);

    // Offset active pressure from equilibrium
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_a[idx0] += 1.0/3.0*gradient_output_x[idx0];
        Fy_a[idx0] += 1.0/3.0*gradient_output_y[idx0];
        Fz_a[idx0] += 1.0/3.0*gradient_output_z[idx0];
    }

    // Compute gradient of rho_p;
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        gradient_input[idx0] = rho_p[idx0];
    }
    compute_gradient(sim);

    // Offset passive pressure from equilibrium
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_p[idx0] += 1.0/3.0*gradient_output_x[idx0];
        Fy_p[idx0] += 1.0/3.0*gradient_output_y[idx0];
        Fz_p[idx0] += 1.0/3.0*gradient_output_z[idx0];
    }

    compute_matrix_divergence(matrices->sigma_nem, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_a[idx0] += matrices->matrix_divergence_output_x[idx0];
        Fy_a[idx0] += matrices->matrix_divergence_output_y[idx0];
        Fz_a[idx0] += matrices->matrix_divergence_output_z[idx0];
    }
}

void compute_final_forces(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)
    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)
    READ_FIELD(u_p)
    READ_FIELD(v_p)
    READ_FIELD(w_p)
    
    WRITE_FIELD(Fx_a)
    WRITE_FIELD(Fy_a)
    WRITE_FIELD(Fz_a)
    WRITE_FIELD(Fx_p)
    WRITE_FIELD(Fy_p)
    WRITE_FIELD(Fz_p)

    // Compute exchange forces
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double phi_i = phi[idx0];
        const double u_a_i = u_a[idx0];
        const double v_a_i = v_a[idx0];
        const double w_a_i = w_a[idx0];
        const double u_p_i = u_p[idx0];
        const double v_p_i = v_p[idx0];
        const double w_p_i = w_p[idx0];

        const double mom_x = params->gamma * phi_i * (1.0 - phi_i) * (u_p_i - u_a_i);
        const double mom_y = params->gamma * phi_i * (1.0 - phi_i) * (v_p_i - v_a_i);
        const double mom_z = params->gamma * phi_i * (1.0 - phi_i) * (w_p_i - w_a_i);

        Fx_a[idx0] += mom_x;
        Fy_a[idx0] += mom_y;
        Fz_a[idx0] += mom_z;

        Fx_p[idx0] -= mom_x;
        Fy_p[idx0] -= mom_y;
        Fz_p[idx0] -= mom_z;
    }
}