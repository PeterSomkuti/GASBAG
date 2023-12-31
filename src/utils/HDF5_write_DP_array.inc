call h5screate_simple_f(rank(array), dims, dspace_id, hdferr)
call check_hdf_error(hdferr, fname, "Error creating simple dataspace")

call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
call check_hdf_error(hdferr, fname, "Error creating dataset creation property list")

! This is for compression
call h5pset_chunk_f(dcpl, size(dims), dims, hdferr)
call h5pset_deflate_f(dcpl, 1, hdferr)

if (present(fill_value)) then
    call h5pset_fill_value_f(dcpl, H5T_NATIVE_REAL, conv_fill_value, hdferr)
    call check_hdf_error(hdferr, fname, "Error setting dataset fill value")
end if

call h5dcreate_f(file_id, dset_name, H5T_NATIVE_REAL, dspace_id, &
                 dset_id, hdferr, dcpl)
call check_hdf_error(hdferr, fname, "Error creating dataset: " // trim(dset_name))

call h5dwrite_f(dset_id, H5T_NATIVE_REAL, conv_array, dims, hdferr)
call check_hdf_error(hdferr, fname, "Error writing data into dataset: " // trim(dset_name))
