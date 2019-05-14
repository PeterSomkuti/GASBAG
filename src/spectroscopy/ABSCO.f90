module absco_mod

  use file_utils_mod, only: read_DP_hdf_dataset, get_HDF5_dset_dims, check_hdf_error
  use math_utils_mod, only: searchsorted_dp
  use control_mod, only: MCS, CS_gas
  use logger_mod, only: logger => master_logger
  use HDF5

  public read_absco_HDF

contains

  subroutine read_absco_HDF(filename, gas, absco_dims, wl_min, wl_max)

    implicit none
    character(len=*), intent(in) :: filename
    type(CS_gas), intent(inout) :: gas
    integer, intent(inout) :: absco_dims
    double precision, intent(in) :: wl_min, wl_max ! Not used yet!

    character(len=*), parameter :: fname = "read_absco_HDF"
    integer(hid_t) :: absco_file_id, dset_id, filetype
    integer(hsize_t), allocatable :: dset_dims(:)
    integer :: hdferr
    character(len=3) :: gas_index ! Is this always going to stay 2-characters long?
    double precision, allocatable :: tmp_absco_3D(:,:,:), tmp_absco_4D(:,:,:,:)
    character(len=999) :: tmp_str
    integer :: wl_idx_left, wl_idx_right, N_wl, i

    call logger%trivia(fname, "Reading in ABSCO HDF file at: " // trim(filename))

    ! Try and open the ABSCO hdf file
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, absco_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Could not open ABSCO HDF file at: " // trim(filename))

    ! And now populate the required fields of the 'gas' structure ..

    ! First, find out which gas_index this particular file has
    call h5dopen_f(absco_file_id, "/Gas_Index", dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Could not open dataset /Gas_Index")
    call h5dget_type_f(dset_id, filetype, hdferr)

    allocate(dset_dims(1))
    dset_dims(1) = 3
    call h5dread_f(dset_id, filetype, gas_index, dset_dims, hdferr)

    if (allocated(dset_dims)) deallocate(dset_dims)
    call check_hdf_error(hdferr, fname, "Error reading in: /Gas_Index")

    ! Let the user know which gas index this ABSCO file has
    call logger%debug(fname, "This ABSCO has a Gas_Index: " // gas_index)

    ! Grab the full ABSCO table
    call logger%debug(fname, "Starting to read in cross section data..")

    call get_HDF5_dset_dims(absco_file_id, "Gas_" // gas_index(1:2) // "_Absorption", dset_dims)
    absco_dims = size(dset_dims)
    deallocate(dset_dims)

    if (absco_dims == 4) then
       call logger%debug(fname, "ABSCO file has 4 dimensions")
       call read_DP_hdf_dataset(absco_file_id, "Gas_" // gas_index(1:2) // "_Absorption", &
            gas%cross_section, dset_dims)
       gas%has_h2o = .true.
    else if (absco_dims == 3) then
       ! The gas%cross_section array is 4-dimensional, so we need an in-between step to
       ! copy the 3-dimensional HDF array into the 4-dimensional array. We simply copy it
       ! all into a 3-dim array first, and then allocate the gas%cross_section array, but
       ! leave the H2O-dimension (dim number 2) as 1-element wide. Then it's simply copied
       ! over, skipping the H2O dimension.

       call logger%debug(fname, "ABSCO file has 3 dimensions")
       call read_DP_hdf_dataset(absco_file_id, "Gas_" // gas_index(1:2) // "_Absorption", tmp_absco_3D, dset_dims)
       ! Copy to gas structure
       if (allocated(gas%cross_section)) deallocate(gas%cross_section)
       allocate(gas%cross_section(dset_dims(1), 1, dset_dims(2), dset_dims(3)))
       gas%cross_section(:,1,:,:) = tmp_absco_3D(:,:,:)
       gas%has_h2o = .false.
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
    ! We also have to re-arrange the arrays
    N_wl = size(gas%wavelength)

    call logger%debug(fname, "Re-arranging spectroscopy arrays.")
    gas%wavelength(:) = gas%wavelength(N_wl:1:-1)
    gas%cross_section(:,:,:,:) = gas%cross_section(N_wl:1:-1,:,:,:)
    call logger%debug(fname, "..done!")


    ! Now if ABSCO is 4-dimensional, we also need the broadener gas data
    if (absco_dims == 4) then

       ! While the ABSCO files have the broadener in generic terms, we just assume
       ! that it's the vater vapour dimension that we are looking at. This can potentially
       ! change in the future (?).

       deallocate(dset_dims)
       call h5dopen_f(absco_file_id, "/Broadener_Index", dset_id, hdferr)
       call check_hdf_error(hdferr, fname, "Could not open dataset /Broadener_Index")
       call h5dget_type_f(dset_id, filetype, hdferr)

       allocate(dset_dims(1))
       dset_dims(1) = 1
       call h5dread_f(dset_id, filetype, gas_index, dset_dims, hdferr)
       if (allocated(dset_dims)) deallocate(dset_dims)
       call check_hdf_error(hdferr, fname, "Error reading in: /Broadener_Index")

       call read_DP_hdf_dataset(absco_file_id, "Broadener_" // gas_index(1:2) // "_VMR", gas%H2O, dset_dims)

    end if 


  end subroutine read_absco_HDF



end module absco_mod
