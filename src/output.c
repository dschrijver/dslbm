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
    double *rho_N = glob_fields->rho_N;

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    // Create file
    sprintf(filename, "data_%d.h5", params->n_output);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, params->fapl_id);

    // Write time
    hid_t dset_scalar = H5Dcreate2(file_id, "t", H5T_NATIVE_INT, params->scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    H5Dclose(dset_scalar);

    output_global_field(rho, "rho", file_id, sim);
    output_global_field(pressure, "pressure", file_id, sim);
    output_global_field(u, "u", file_id, sim);
    output_global_field(v, "v", file_id, sim);
    output_global_field(w, "w", file_id, sim);
    output_global_field(rho_N, "rho_N", file_id, sim);
    output_global_field(glob_fields->Gx, "Gx", file_id, sim);
    output_global_field(glob_fields->Gy, "Gy", file_id, sim);
    output_global_field(glob_fields->Gz, "Gz", file_id, sim);

    output_comp_field(rho_comp, "rho", file_id, sim);
    output_comp_field(Fx, "Fx", file_id, sim);
    output_comp_field(Fy, "Fy", file_id, sim);
    output_comp_field(Fz, "Fz", file_id, sim);
    output_comp_field(comp_fields->Qx, "Qx", file_id, sim);
    output_comp_field(comp_fields->Qy, "Qy", file_id, sim);
    output_comp_field(comp_fields->Qz, "Qz", file_id, sim);

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