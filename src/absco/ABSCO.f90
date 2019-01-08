module absco_mod

  use file_utils_mod, only: read_DP_hdf_dataset, get_HDF5_dset_dims, check_hdf_error
  use control_mod, only: CS_gas
  use logger_mod, only: logger => master_logger
  use HDF5

  public read_absco_HDF

contains

  subroutine read_absco_HDF(filename, gas, absco_dims)

    implicit none
    character(len=*), intent(in) :: filename
    type(CS_gas), intent(inout) :: gas
    integer, intent(inout) :: absco_dims

    character(len=*), parameter :: fname = "read_absco_HDF"
    integer(hid_t) :: absco_file_id, dset_id, filetype
    integer(hsize_t), allocatable :: dset_dims(:)
    integer :: hdferr
    character(len=2) :: gas_index ! Is this always going to stay 2-characters long?
    double precision, allocatable :: absco_3D(:,:,:)
    character(len=999) :: tmp_str

    call logger%trivia(fname, "Reading in ABSCO HDF file at: " // trim(filename))

    ! Try and open the ABSCO hdf file
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, absco_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Could not open ABSCO HDF file at: " // trim(filename))

    ! And now populate the required fields of the 'gas' structure ..

    ! First, find out which gas_index this particular file has
    call h5dopen_f(absco_file_id, "/Gas_Index", dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Could not open dataset /Gas_Index")
    call h5dget_type_f(dset_id, filetype, hdferr)

    call h5dread_f(dset_id, filetype, gas_index, dset_dims, hdferr)
    if (allocated(dset_dims)) deallocate(dset_dims)
    call check_hdf_error(hdferr, fname, "Error reading in: /Gas_Index")

    ! Let the user know which gas index this ABSCO file has
    call logger%debug(fname, "This ABSCO has a Gas_Index: " // gas_index)

    ! Grab the full ABSCO table
    call logger%debug(fname, "Starting to read in cross section data..")

    call get_HDF5_dset_dims(absco_file_id, "Gas_" // gas_index // "_Absorption", dset_dims)
    absco_dims = size(dset_dims)
    deallocate(dset_dims)

    if (absco_dims == 4) then
       call logger%debug(fname, "ABSCO file has 4 dimensions")
       call read_DP_hdf_dataset(absco_file_id, "Gas_" // gas_index // "_Absorption", gas%cross_section, dset_dims)
    else if (absco_dims == 3) then
       call logger%debug(fname, "ABSCO file has 3 dimensions")
       call read_DP_hdf_dataset(absco_file_id, "Gas_" // gas_index // "_Absorption", absco_3D, dset_dims)
       ! Copy to gas structure
       allocate(gas%cross_section(dset_dims(1), dset_dims(2), dset_dims(3), 1))
       gas%cross_section(:,:,:,1) = absco_3D(:,:,:)
    else
       write(tmp_str, '(A, G0.1, A)') "ABSCO file has ", absco_dims, " dimensions. This is not (yet) supported."
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if
    call logger%debug(fname, "..done!")

    deallocate(dset_dims)
    call read_DP_hdf_dataset(absco_file_id, "Temperature", gas%T, dset_dims)
    deallocate(dset_dims)
    call read_DP_hdf_dataset(absco_file_id, "Pressure", gas%p, dset_dims)

    deallocate(dset_dims)
    call read_DP_hdf_dataset(absco_file_id, "Wavenumber", gas%wavelength, dset_dims)
    ! ABSCO comes in wavenumber, but we'd rather work in wavelengths [microns], so let's convert
    gas%wavelength = 1.0d4 / gas%wavelength

    ! Now if ABSCO is 4-dimensional, we also need the broadener gas data
    if (absco_dims == 4) then

       ! While the ABSCO files have the broadener in generic terms, we just assume
       ! that it's the vater vapour dimension that we are looking at. This can potentially
       ! change in the future (?).

       deallocate(dset_dims)
       call h5dopen_f(absco_file_id, "/Broadener_Index", dset_id, hdferr)
       call check_hdf_error(hdferr, fname, "Could not open dataset /Broadener_Index")
       call h5dget_type_f(dset_id, filetype, hdferr)

       call h5dread_f(dset_id, filetype, gas_index, dset_dims, hdferr)
       if (allocated(dset_dims)) deallocate(dset_dims)
       call check_hdf_error(hdferr, fname, "Error reading in: /Broadener_Index")

       call read_DP_hdf_dataset(absco_file_id, "Broadener_" // gas_index // "_VMR", gas%H2O, dset_dims)

    end if 


  end subroutine read_absco_HDF



end module absco_mod