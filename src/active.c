#include "../include/datatypes.h"
#include "../include/communicate.h"
#include "../include/matrix.h"

void compute_matrix_gradient(matrix *input, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    WRITE_MATRIX(matrix_gradient_output_x)
    WRITE_MATRIX(matrix_gradient_output_y)
    WRITE_MATRIX(matrix_gradient_output_z)

    communicate_matrix(input, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        FOR_MAT
        {
            double gradient_x_i = 0.0;
            double gradient_y_i = 0.0;
            double gradient_z_i = 0.0;
            for (int p = 0; p < NP; p++)
            {
                // TODO: boundary conditions
                const double input_local = input[INDEX(i + cx[p], j + cy[p], k + cz[p])][a][b];

                gradient_x_i += wp[p] * input_local * (double)cx[p];
                gradient_y_i += wp[p] * input_local * (double)cy[p];
                gradient_z_i += wp[p] * input_local * (double)cz[p];
            }
            matrix_gradient_output_x[idx0][a][b] = 3.0 * gradient_x_i;
            matrix_gradient_output_y[idx0][a][b] = 3.0 * gradient_y_i;
            matrix_gradient_output_z[idx0][a][b] = 3.0 * gradient_z_i;
        }
    }
}

void compute_matrix_laplacian(matrix *input, matrix *output, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    communicate_matrix(input, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        FOR_MAT
        {
            double laplacian_i = 0.0;
            for (int p = 0; p < NP; p++)
            {
                // TODO: boundary conditions
                const double input_local = input[INDEX(i + cx[p], j + cy[p], k + cz[p])][a][b];

                laplacian_i += wp[p] * input_local;
            }
            output[idx0][a][b] = 6.0 * (laplacian_i - input[idx0][a][b]);
        }
    }
}

void compute_matrix_divergence(matrix *input, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID
    UNPACK_STENCIL

    communicate_matrix(input, sim);

    double *output[3] = {matrices->matrix_divergence_output_x, matrices->matrix_divergence_output_y, matrices->matrix_divergence_output_z};

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        for (int a = 0; a < 3; a++)
        {
            double div = 0.0;
            
            for (int p = 0; p < NP; p++)
            {
                const int idxp = INDEX(i + cx[p], j + cy[p], k + cz[p]);

                // TODO: boundary conditions

                div += wp[p] * (input[idxp][a][0] * (double)cx[p] + 
                                input[idxp][a][1] * (double)cy[p] + 
                                input[idxp][a][2] * (double)cz[p]);
            }
            
            output[a][idx0] = 3.0*div;
        }
    }
}

void compute_molecular_field(SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    READ_FIELD(phi)

    matrix *Q = matrices->Q;
    matrix *H = matrices->Q;

    PARAM(S_nem)
    PARAM(K_LC)
    PARAM(C)

    compute_matrix_laplacian(matrices->Q, matrices->H, sim);

    FOR_DOMAIN
    {
        const int idx0 = INDEX(i, j, k);

        const double phi_i = phi[idx0];

        const double contraction = matcontraction(Q[idx0], Q[idx0]);
        const double prefac = 4.0/3.0*C*(S_nem*S_nem*phi_i - 2.0/3.0*contraction);

        FOR_MAT
        {
            H[idx0][a][b] = prefac*Q[idx0][a][b] + K_LC*H[idx0][a][b];
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

    PARAM(K_LC)
    PARAM(lambda)
    PARAM(zeta)

    matrix *gradient[3] = {matrices->matrix_gradient_output_x, matrices->matrix_gradient_output_y, matrices->matrix_gradient_output_z};

    compute_matrix_gradient(matrices->Q, sim);

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