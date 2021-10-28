!> @brief ABSCO-related spectroscopy functions
!> @file ABSCO.f90
!> @author Peter Somkuti
!>
!> @details
!> This module contains functions to handle ABSCO-type spectroscopy data




module absco_mod

  use file_utils_mod, only: read_DP_hdf_dataset, get_HDF5_dset_dims, check_hdf_error
  use math_utils_mod, only: searchsorted_dp
  use control_mod, only: CS_gas_t
  use logger_mod, only: logger => master_logger
  use HDF5

  public read_absco_HDF

contains

  !> @brief Loads ABSCO spectroscopy into memory
  !> @param filename Filename of the ABSCO file
  !> @param gas Gas object entity
  !> @param hitran_index Read this specific gas
  !> @param absco_dims The dimensions of the full ABSCO array
  !> @detail
  !> This function first inspects an ABSCO file to see whether it has
  !> three or four dimensions, and then reads the full contents accordingly.
  !> The coefficients are rearranged to be in increasing wavelength order,
  !> however are not re-sampled at all in this function.
  !> If a HITRAN index other than -1 is supplied to the function, the
  !> corresponding gas will be read, and produces an error if that gas
  !> is not present in the file.
  subroutine read_absco_HDF(filename, gas, hitran_index, absco_dims)

    implicit none
    character(len=*), intent(in) :: filename
    type(CS_gas_t), intent(inout) :: gas
    integer, intent(in) :: hitran_index
    integer, intent(inout) :: absco_dims

    character(len=*), parameter :: fname = "read_absco_HDF"
    integer(hid_t) :: absco_file_id, dset_id, filetype
    integer(hsize_t), allocatable :: dset_dims(:)
    integer :: hdferr
    character(len=3) :: gas_index ! Is this always going to stay 2-characters long?
    double precision, allocatable :: tmp_absco_3D(:,:,:)
    character(len=999) :: tmp_str, gas_str
    logical :: gas_exists
    integer :: N_wl

    integer, parameter :: max_attempts = 10
    integer :: i


    call logger%trivia(fname, "Reading in ABSCO HDF file at: " // trim(filename))

    ! Try and open the ABSCO hdf file

    ! This can occasionally fail in a processing environment with many
    ! processes trying to access the same file(s). We can try to just open
    ! the file a few times and if it still does not work after some
    ! attempts - call it a day and quit.


    hdferr = -1
    i = 0
    do while ((i < max_attempts) .and. (hdferr < 0))
       call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, absco_file_id, hdferr)
       i = i + 1

       if (hdferr < 0) then
          write(tmp_str, "(A, G0.1, A, A)") "Attempt #", i, "to read ABSCO file: ", trim(filename)
          call logger%debug(fname, trim(tmp_str))
       end if

    end do

    call check_hdf_error(hdferr, fname, "Could not open ABSCO HDF file at: " // trim(filename))

    ! And now populate the required fields of the 'gas' structure ..

    ! If the user supplied a HITRAN index, we do not need to look for which
    ! gas(es) are in this file, we only need to check if the supplied index
    ! also appears in the file. Default value is -1, which means user does
    ! not know which gas index is used.

    if (hitran_index /= -1) then

       write(gas_str, "(A, I0.2, A)") "Gas_", hitran_index, "_Absorption"
       call h5lexists_f(absco_file_id, trim(gas_str), &
            gas_exists, hdferr)
       call check_hdf_error(hdferr, fname, "Error reading in: " // trim(gas_str))

    else

       ! find out which gas_index this particular file has
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

       write(gas_str, "(A,A,A)") "Gas_" // gas_index(1:2) // "_Absorption"

    end if

    ! Grab the full ABSCO table
    call logger%debug(fname, "Starting to read in cross section data..")

    call get_HDF5_dset_dims(absco_file_id, trim(gas_str), dset_dims)
    absco_dims = size(dset_dims)
    deallocate(dset_dims)

    if (absco_dims == 4) then
       call logger%debug(fname, "ABSCO file has 4 dimensions")
       call read_DP_hdf_dataset(absco_file_id, trim(gas_str), &
            gas%cross_section, dset_dims)
       gas%has_h2o = .true.
    else if (absco_dims == 3) then
       ! The gas%cross_section array is 4-dimensional, so we need an in-between step to
       ! copy the 3-dimensional HDF array into the 4-dimensional array. We simply copy it
       ! all into a 3-dim array first, and then allocate the gas%cross_section array, but
       ! leave the H2O-dimension (dim number 2) as 1-element wide. Then it's simply copied
       ! over, skipping the H2O dimension.

       call logger%debug(fname, "ABSCO file has 3 dimensions")
       call read_DP_hdf_dataset(absco_file_id, trim(gas_str), tmp_absco_3D, dset_dims)
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
    gas%wavelength(:) = 1.0d4 / gas%wavelength(:)
    ! We also have to re-arrange the arrays
    N_wl = size(gas%wavelength)

    write(tmp_str, '(A,G0.1,A)') "There are ", N_wl, " spectral points in this ABSCO file."
    call logger%debug(fname, trim(tmp_str))

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

    call h5fclose_f(absco_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Could not close ABSCO file.")



  end subroutine read_absco_HDF



end module absco_mod
