#include <hdf5.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/output.h"

void output_data(SimulationBag *sim)
{
    UNPACK_BAGS

    char filename[32];

    int t = params->t;

    READ_FIELD(rho)
    READ_FIELD(rho_RED)
    READ_FIELD(rho_BLUE)
    READ_FIELD(pressure)
    READ_FIELD(u)
    READ_FIELD(v)
    READ_FIELD(w)
    READ_FIELD(Fx)
    READ_FIELD(Fy)
    READ_FIELD(Fz)

    READ_FIELD(rho_N)
    READ_FIELD(Gx)
    READ_FIELD(Gy)
    READ_FIELD(Gz)
    READ_FIELD(nx)
    READ_FIELD(ny)
    READ_FIELD(nz)

    READ_FIELD(Qx)
    READ_FIELD(Qy)
    READ_FIELD(Qz)

    // Create file
    sprintf(filename, "data_%d.h5", params->n_output);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, params->fapl_id);

    // Write time
    hid_t dset_scalar = H5Dcreate2(file_id, "t", H5T_NATIVE_INT, params->scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    H5Dclose(dset_scalar);

    // Create groups
    hid_t hydro = H5Gcreate2(file_id, "/hydro", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t cg = H5Gcreate2(file_id, "/cg", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t other = H5Gcreate2(file_id, "/other", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    output_field(rho, "rho", hydro, sim);
    output_field(rho_RED, "rho_RED", hydro, sim);
    output_field(rho_BLUE, "rho_BLUE", hydro, sim);
    output_field(pressure, "pressure", hydro, sim);
    output_field(u, "u", hydro, sim);
    output_field(v, "v", hydro, sim);
    output_field(w, "w", hydro, sim);
    output_field(Fx, "Fx", hydro, sim);
    output_field(Fy, "Fy", hydro, sim);
    output_field(Fz, "Fz", hydro, sim);

    output_field(rho_N, "rho_N", cg, sim);
    output_field(Gx, "Gx", cg, sim);
    output_field(Gy, "Gy", cg, sim);
    output_field(Gz, "Gz", cg, sim);
    output_field(nx, "nx", cg, sim);
    output_field(ny, "ny", cg, sim);
    output_field(nz, "nz", cg, sim);

    output_field(Qx, "Qx", other, sim);
    output_field(Qy, "Qy", other, sim);
    output_field(Qz, "Qz", other, sim);

    // Close groups
    H5Gclose(other);
    H5Gclose(cg);
    H5Gclose(hydro);

    // Close file
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    params->n_output++;
}

void output_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int i_start = params->i_start;
    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;

    // Create dataset
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, H5T_NATIVE_DOUBLE, params->filespace, H5P_DEFAULT, params->dcpl_id, H5P_DEFAULT);

    // File hyperslab
    hsize_t start_file[3] = {i_start, 0, 0};
    hsize_t count[3] = {NX_proc, NY, NZ};
    H5Sselect_hyperslab(params->filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Process hyperslab
    hsize_t start_proc[3] = {2, 0, 0};
    H5Sselect_hyperslab(params->memspace, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, params->memspace, params->filespace, params->dxpl_id, field);

    H5Dclose(dset_id);
}