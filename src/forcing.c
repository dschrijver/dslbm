#include "../include/datatypes.h"
#include "../include/fields.h"
#include "../include/active.h"

void compute_preliminary_forces(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(mu_a)
    READ_FIELD(mu_p)
    READ_FIELD(dx_rho)
    READ_FIELD(dy_rho)
    READ_FIELD(dz_rho)
    READ_FIELD(rho_a)
    READ_FIELD(rho_p)
    READ_FIELD(phi)
    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)
    READ_FIELD(u_p)
    READ_FIELD(v_p)
    READ_FIELD(w_p)
    READ_FIELD(Fx_nem)
    READ_FIELD(Fy_nem)
    READ_FIELD(Fz_nem)
    
    WRITE_FIELD(Fx_a)
    WRITE_FIELD(Fy_a)
    WRITE_FIELD(Fz_a)
    WRITE_FIELD(Fx_p)
    WRITE_FIELD(Fy_p)
    WRITE_FIELD(Fz_p)

    // Compute gradient of mu_a
    compute_gradient(fields->mu_a, fields->Fx_a, fields->Fy_a, fields->Fz_a, sim);

    // Add to active forces
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        const double rho_a_i = rho_a[idx0];

        Fx_a[idx0] *= -rho_a_i;
        Fy_a[idx0] *= -rho_a_i;
        Fz_a[idx0] *= -rho_a_i;
    }

    // Compute gradient of mu_p
    compute_gradient(fields->mu_p, fields->Fx_p, fields->Fy_p, fields->Fz_p, sim);

    // Add to passive forces
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);
        const double rho_p_i = rho_p[idx0];

        Fx_p[idx0] *= -rho_p_i;
        Fy_p[idx0] *= -rho_p_i;
        Fz_p[idx0] *= -rho_p_i;
    }

    // Compute gradient of rho_a
    compute_gradient(fields->rho_a, fields->dx_rho, fields->dy_rho, fields->dz_rho, sim);

    // Offset active pressure from equilibrium
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_a[idx0] += 1.0/3.0*dx_rho[idx0];
        Fy_a[idx0] += 1.0/3.0*dy_rho[idx0];
        Fz_a[idx0] += 1.0/3.0*dz_rho[idx0];
    }

    // Compute gradient of rho_p;
    compute_gradient(fields->rho_p, fields->dx_rho, fields->dy_rho, fields->dz_rho, sim);

    // Offset passive pressure from equilibrium
    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_p[idx0] += 1.0/3.0*dx_rho[idx0];
        Fy_p[idx0] += 1.0/3.0*dy_rho[idx0];
        Fz_p[idx0] += 1.0/3.0*dz_rho[idx0];
    }

    // Add nematic force
    compute_matrix_divergence(matrices->sigma_nem, fields->Fx_nem, fields->Fy_nem, fields->Fz_nem, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        Fx_a[idx0] += Fx_nem[idx0];
        Fy_a[idx0] += Fy_nem[idx0];
        Fz_a[idx0] += Fz_nem[idx0];
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