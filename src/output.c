#include <hdf5.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/output.h"

void output_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
{
    UNPACK_BAGS
    UNPACK_GRID

    // Create dataset
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, H5T_NATIVE_DOUBLE, params->filespace, H5P_DEFAULT, params->dcpl_id, H5P_DEFAULT);

    // File hyperslab
    hsize_t start_file[3] = {i_start, j_start, k_start};
    hsize_t count[3] = {NX_proc, NY_proc, NZ_proc};
    H5Sselect_hyperslab(params->filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Process hyperslab
    hsize_t start_proc[3] = {1, 1, 1};
    H5Sselect_hyperslab(params->memspace, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, params->memspace, params->filespace, params->dxpl_id, field);

    H5Dclose(dset_id);
}

void output_data(SimulationBag *sim)
{
    UNPACK_BAGS

    char filename[32];

    int t = params->t;

    // Create file
    sprintf(filename, "data_%d.h5", params->n_output);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, params->fapl_id);

    // Write time
    hid_t dset_scalar = H5Dcreate2(file_id, "t", H5T_NATIVE_INT, params->scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    H5Dclose(dset_scalar);

    hid_t main = H5Gcreate2(file_id, "/main", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t densities = H5Gcreate2(file_id, "/densities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t velocities = H5Gcreate2(file_id, "/velocities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t forces = H5Gcreate2(file_id, "/forces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    output_field(fields->rho, "rho", main, sim);
    output_field(fields->phi, "phi", main, sim);
    output_field(fields->u, "u", main, sim);
    output_field(fields->v, "v", main, sim);
    output_field(fields->w, "w", main, sim);

    output_field(fields->rho_a, "rho_a", densities, sim);
    output_field(fields->rho_p, "rho_p", densities, sim);

    output_field(fields->u_a, "u_a", velocities, sim);
    output_field(fields->v_a, "v_a", velocities, sim);
    output_field(fields->w_a, "w_a", velocities, sim);
    output_field(fields->u_p, "u_p", velocities, sim);
    output_field(fields->v_p, "v_p", velocities, sim);
    output_field(fields->w_p, "w_p", velocities, sim);

    output_field(fields->Fx_a, "Fx_a", forces, sim);
    output_field(fields->Fy_a, "Fy_a", forces, sim);
    output_field(fields->Fz_a, "Fz_a", forces, sim);
    output_field(fields->Fx_p, "Fx_p", forces, sim);
    output_field(fields->Fy_p, "Fy_p", forces, sim);
    output_field(fields->Fz_p, "Fz_p", forces, sim);

    // Close groups
    H5Gclose(forces);
    H5Gclose(velocities);
    H5Gclose(densities);
    H5Gclose(main);

    // Close file
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    params->n_output++;
}