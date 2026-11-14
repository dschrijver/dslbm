#include "../include/datatypes.h"
#include "../include/communicate.h"
#include "../include/matrix.h"

void compute_matrix_gradient(matrix *input, matrix *output_x, matrix *output_y, matrix *output_z, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        matrix gradient_x = {{0.0}};
        matrix gradient_y = {{0.0}};
        matrix gradient_z = {{0.0}};

        for (int p = 0; p < NP; p++)
        {
            const int idxp = INDEX(i + cx[p], j + cy[p], k + cz[p]);
            const double weight = wp[p];

            FOR_MAT
            {
                const double res = weight * input[idxp][a][b];

                gradient_x[a][b] += res * (double)cx[p];
                gradient_y[a][b] += res * (double)cy[p];
                gradient_z[a][b] += res * (double)cz[p];
            }
        }

        FOR_MAT
        {
            output_x[idx0][a][b] = 3.0 * gradient_x[a][b];
            output_y[idx0][a][b] = 3.0 * gradient_y[a][b];
            output_z[idx0][a][b] = 3.0 * gradient_z[a][b];
        }
    }
}

void compute_matrix_laplacian(matrix *input, matrix *output, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        matrix laplacian = {{0.0}};
        for (int p = 0; p < NP; p++)
        {
            // TODO: boundary conditions
            const int idxp = INDEX(i + cx[p], j + cy[p], k + cz[p]);
            const double weight = wp[p];

            FOR_MAT
            {
                laplacian[a][b] += weight * input[idxp][a][b];
            }
        }

        FOR_MAT
        {
            output[idx0][a][b] = 6.0 * (laplacian[a][b] - input[idx0][a][b]);
        }
    }
}

void compute_matrix_divergence(matrix *input, double *output_x, double *output_y, double *output_z, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    double *output[3] = {output_x, output_y, output_z};

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        double divergence[3] = {0.0};

        for (int p = 0; p < NP; p++)
        {
            const int idxp = INDEX(i + cx[p], j + cy[p], k + cz[p]);
            const double weight = wp[p];

            for (int a = 0; a < 3; a++)
            {
                divergence[a] += weight * (input[idxp][a][0] * (double)cx[p] + 
                                           input[idxp][a][1] * (double)cy[p] + 
                                           input[idxp][a][2] * (double)cz[p]);
            }
        }

        for (int a = 0; a < 3; a++)
        {
            output[a][idx0] = 3.0*divergence[a];
        }
    }
}

void compute_molecular_field(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)

    matrix *Q = matrices->Q;
    matrix *H = matrices->H;

    PARAM(S_nem)
    PARAM(K_LC)
    PARAM(C)

    compute_matrix_laplacian(Q, H, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double phi_i = phi[idx0];

        const double contraction = matcontraction(Q[idx0], Q[idx0]);
        const double prefac = 4.0/3.0*C*(S_nem*S_nem*phi_i - 2.0/3.0*contraction);

        FOR_MAT
        {
            H[idx0][a][b] = prefac*Q[idx0][a][b] + K_LC*H[idx0][a][b]; // H contains laplacian
        }
    }
}

void compute_nematic_stress(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)

    matrix *Q = matrices->Q;
    matrix *H = matrices->H;
    matrix *sigma_nem = matrices->sigma_nem;
    matrix *dx_Q = matrices->dx_Q;
    matrix *dy_Q = matrices->dy_Q;
    matrix *dz_Q = matrices->dz_Q;

    PARAM(K_LC)
    PARAM(lambda)
    PARAM(zeta)

    matrix *gradient[3] = {dx_Q, dy_Q, dz_Q};

    compute_matrix_gradient(Q, dx_Q, dy_Q, dz_Q, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double phi_i = phi[idx0];

        const double QH_contraction = matcontraction(Q[idx0], H[idx0]);
        matrix HQ_product; 
        matmul(H[idx0], Q[idx0], HQ_product);
        matrix QH_product;
        matmul(Q[idx0], H[idx0], QH_product);

        FOR_MAT
        {
            const double dQdQ_contraction = matcontraction(gradient[a][idx0], gradient[b][idx0]);
            const double sigma_el = 2.0*lambda*QH_contraction*(Q[idx0][a][b] + delta(a, b)/3.0) - lambda*(HQ_product[a][b] + H[idx0][a][b]/3.0) - lambda*(QH_product[a][b] + H[idx0][a][b]/3.0) - K_LC*dQdQ_contraction + QH_product[a][b] - HQ_product[a][b];
            const double sigma_act = - zeta*Q[idx0][a][b];

            sigma_nem[idx0][a][b] = phi_i * (sigma_el + sigma_act);
        }
    }
}

void evolve_Q(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(u_a)
    READ_FIELD(v_a)
    READ_FIELD(w_a)

    matrix *H = matrices->H;
    matrix *nabla_u = matrices->nabla_u;
    matrix *E = matrices->E;
    matrix *Omega = matrices->Omega;

    PARAM(lambda)
    PARAM(Gamma_Q)

    const int N_active_substeps = params->N_active_substeps;

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        // TODO: boundary conditions
        nabla_u[idx0][0][0] = 0.5 * (u_a[INDEX(i+1, j, k)] - u_a[INDEX(i-1, j, k)]);
        nabla_u[idx0][0][1] = 0.5 * (u_a[INDEX(i, j+1, k)] - u_a[INDEX(i, j-1, k)]);
        nabla_u[idx0][0][2] = 0.5 * (u_a[INDEX(i, j, k+1)] - u_a[INDEX(i, j, k-1)]);

        nabla_u[idx0][1][0] = 0.5 * (v_a[INDEX(i+1, j, k)] - v_a[INDEX(i-1, j, k)]);
        nabla_u[idx0][1][1] = 0.5 * (v_a[INDEX(i, j+1, k)] - v_a[INDEX(i, j-1, k)]);
        nabla_u[idx0][1][2] = 0.5 * (v_a[INDEX(i, j, k+1)] - v_a[INDEX(i, j, k-1)]);

        nabla_u[idx0][2][0] = 0.5 * (w_a[INDEX(i+1, j, k)] - w_a[INDEX(i-1, j, k)]);
        nabla_u[idx0][2][1] = 0.5 * (w_a[INDEX(i, j+1, k)] - w_a[INDEX(i, j-1, k)]);
        nabla_u[idx0][2][2] = 0.5 * (w_a[INDEX(i, j, k+1)] - w_a[INDEX(i, j, k-1)]);

        FOR_MAT
        {
            E[idx0][a][b] = 0.5 * (nabla_u[idx0][a][b] + nabla_u[idx0][b][a]);
            Omega[idx0][a][b] = 0.5 * (nabla_u[idx0][a][b] - nabla_u[idx0][b][a]);
        }
    }

    for (int n = 0; n < N_active_substeps; n++)
    {
        // Because they get swapped every step
        matrix *Q = matrices->Q;
        matrix *Q_new = matrices->Q_new;

        communicate_matrix(matrices->Q, sim);

        compute_molecular_field(sim);

        FOR_DOMAIN
        {
            const int idx0 = INDEX(i, j, k);

            matrix u_nabla_Q;

            FOR_MAT
            {
                // TODO: boundary conditions
                u_nabla_Q[a][b] = u_a[idx0] * 0.5 * (Q[INDEX(i+1, j, k)][a][b] - Q[INDEX(i-1, j, k)][a][b]) + 
                                  v_a[idx0] * 0.5 * (Q[INDEX(i, j+1, k)][a][b] - Q[INDEX(i, j-1, k)][a][b]) + 
                                  w_a[idx0] * 0.5 * (Q[INDEX(i, j, k+1)][a][b] - Q[INDEX(i, j, k-1)][a][b]);
            }

            matrix S1;
            matrix S2;
            matrix S3;

            FOR_MAT
            {
                S1[a][b] = lambda*E[idx0][a][b] + Omega[idx0][a][b];
                S2[a][b] = Q[idx0][a][b] + delta(a, b)/3.0;
                S3[a][b] = lambda*E[idx0][a][b] - Omega[idx0][a][b];
            }

            const double prefac = 2.0*lambda*matcontraction(Q[idx0], nabla_u[idx0]);

            matrix S4;
            matmul(S1, S2, S4);

            matrix S5;
            matmul(S2, S3, S5);

            matrix S;

            FOR_MAT
            {
                S[a][b] = S4[a][b] + S5[a][b] - prefac*S2[a][b];
            }

            FOR_MAT
            {
                Q_new[idx0][a][b] = Q[idx0][a][b] + 1.0/(double)params->N_active_substeps * (-u_nabla_Q[a][b] + S[a][b] + Gamma_Q*H[idx0][a][b]);
            }
        }
        matrix *tmp = matrices->Q;
        matrices->Q = matrices->Q_new;
        matrices->Q_new = tmp;
    }
}