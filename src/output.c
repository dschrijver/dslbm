#include <hdf5.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/output.h"

void output_data(SimulationBag *sim)
{
    ParamBag *params = sim->params;
    GlobalFieldBag *glob_fields = sim->glob_fields;
    ComponentFieldBag *comp_fields = sim->comp_fields;

    char filename[32];

    int t = params->t;

    double *rho = glob_fields->rho;
    double *pressure = glob_fields->pressure;
    double *u = glob_fields->u;
    double *v = glob_fields->v;
    double *w = glob_fields->w;
    double *Fx = glob_fields->Fx;
    double *Fy = glob_fields->Fy;
    double *Fz = glob_fields->Fz;

    double *rho_N = glob_fields->rho_N;
    double *Gx = glob_fields->Gx;
    double *Gy = glob_fields->Gy;
    double *Gz = glob_fields->Gz;
    double *nx = glob_fields->nx;
    double *ny = glob_fields->ny;
    double *nz = glob_fields->nz;

    double *Qx = glob_fields->Qx;
    double *Qy = glob_fields->Qy;
    double *Qz = glob_fields->Qz;

    double *rho_comp = comp_fields->rho_comp;

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


    output_global_field(rho, "rho", hydro, sim);
    output_comp_field(rho_comp, "rho", hydro, sim);
    output_global_field(pressure, "pressure", hydro, sim);
    output_global_field(u, "u", hydro, sim);
    output_global_field(v, "v", hydro, sim);
    output_global_field(w, "w", hydro, sim);
    output_global_field(Fx, "Fx", hydro, sim);
    output_global_field(Fy, "Fy", hydro, sim);
    output_global_field(Fz, "Fz", hydro, sim);

    output_global_field(rho_N, "rho_N", cg, sim);
    output_global_field(Gx, "Gx", cg, sim);
    output_global_field(Gy, "Gy", cg, sim);
    output_global_field(Gz, "Gz", cg, sim);
    output_global_field(nx, "nx", cg, sim);
    output_global_field(ny, "ny", cg, sim);
    output_global_field(nz, "nz", cg, sim);

    output_global_field(Qx, "Qx", other, sim);
    output_global_field(Qy, "Qy", other, sim);
    output_global_field(Qz, "Qz", other, sim);

    // Close groups
    H5Gclose(other);
    H5Gclose(cg);
    H5Gclose(hydro);

    // Close file
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    params->n_output++;
}

void output_global_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
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
    H5Sselect_hyperslab(params->memspace_glob, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, params->memspace_glob, params->filespace, params->dxpl_id, field);

    H5Dclose(dset_id);
}

void output_comp_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int i_start = params->i_start;
    int NX_proc = params->NX_proc;
    int NY = params->NY;
    int NZ = params->NZ;

    char fieldcompname[32];

    hsize_t start_file[3] = {i_start, 0, 0};
    hsize_t count_file[3] = {NX_proc, NY, NZ};

    hsize_t start_proc[4] = {2, 0, 0, 0};
    hsize_t count_proc[4] = {NX_proc, NY, NZ, 1};

    char names[NCOMP][5] = {"RED", "BLUE"};

    hid_t dset_id;
    for (int n = 0; n < NCOMP; n++)
    {
        // Create dataset
        sprintf(fieldcompname, "%s_%s", fieldname, names[n]);

        dset_id = H5Dcreate2(loc_id, fieldcompname, H5T_NATIVE_DOUBLE, params->filespace, H5P_DEFAULT, params->dcpl_id, H5P_DEFAULT);

        // File hyperslab
        H5Sselect_hyperslab(params->filespace, H5S_SELECT_SET, start_file, NULL, count_file, NULL);

        // Process hyperslab
        start_proc[3] = n;
        H5Sselect_hyperslab(params->memspace_comp, H5S_SELECT_SET, start_proc, NULL, count_proc, NULL);

        // Write data
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, params->memspace_comp, params->filespace, params->dxpl_id, field);

        H5Dclose(dset_id);
    }
}