#include <hdf5.h>

#include "../include/datatypes.h"
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

    double *rho_comp = comp_fields->rho_comp;
    double *Fx = comp_fields->Fx;
    double *Fy = comp_fields->Fy;
    double *Fz = comp_fields->Fz;

    // Create file
    sprintf(filename, "data_%d.h5", params->n_output);
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    // Write time
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t dset_scalar = H5Dcreate2(file_id, "t", H5T_NATIVE_INT, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
    H5Dclose(dset_scalar);
    H5Sclose(scalar_space);

    output_global_field(rho, "rho", file_id, sim);
    output_global_field(pressure, "pressure", file_id, sim);
    output_global_field(u, "u", file_id, sim);
    output_global_field(v, "v", file_id, sim);
    output_global_field(w, "w", file_id, sim);

    output_comp_field(rho_comp, "rho", file_id, sim);
    output_comp_field(Fx, "Fx", file_id, sim);
    output_comp_field(Fy, "Fy", file_id, sim);
    output_comp_field(Fz, "Fz", file_id, sim);

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
    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    hid_t filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[3] = {NX_proc + 4, NY, NZ};
    hid_t memspace = H5Screate_simple(3, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(loc_id, fieldname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t start_file[3] = {i_start, 0, 0};
    hsize_t count[3] = {NX_proc, NY, NZ};
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Process hyperslab
    hsize_t start_proc[3] = {2, 0, 0};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Sclose(filespace);
}

void output_comp_field(double *field, char *fieldname, hid_t loc_id, SimulationBag *sim)
{
    ParamBag *params = sim->params;

    int i_start = params->i_start;
    int NX_proc = params->NX_proc;
    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;

    char fieldcompname[32];

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};

    // Space occupied in processor memory
    hsize_t dims_proc[4] = {NX_proc + 4, NY, NZ, NCOMP};
    hid_t memspace = H5Screate_simple(4, dims_proc, NULL);

    hsize_t start_file[3] = {i_start, 0, 0};
    hsize_t count_file[3] = {NX_proc, NY, NZ};

    hsize_t start_proc[4] = {2, 0, 0, 0};
    hsize_t count_proc[4] = {NX_proc, NY, NZ, 1};

    char names[NCOMP][5] = {"RED", "BLUE"};

    hid_t dcpl_id, dset_id, filespace, dxpl_id;
    for (int n = 0; n < NCOMP; n++)
    {
        filespace = H5Screate_simple(3, dims_file, NULL);

        // Create dataset
        sprintf(fieldcompname, "%s_%s", fieldname, names[n]);
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        dset_id = H5Dcreate2(loc_id, fieldcompname, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Pclose(dcpl_id);
        H5Sclose(filespace);

        // File hyperslab
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count_file, NULL);

        // Process hyperslab
        start_proc[3] = n;
        H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_proc, NULL, count_proc, NULL);

        // Write data
        dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, field);

        H5Dclose(dset_id);
        H5Pclose(dxpl_id);
    }

    H5Sclose(memspace);
    H5Sclose(filespace);
}