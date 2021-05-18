!> @brief perform retrievals
!> @file perform_retrievals.f90
!> @author Peter Somkuti
!>
!> @details
!> In this subroutine, we essentially perform all the tasks for the retrievals,
!> including setting them up. All necessary information about the settings etc.,
!> is avaliable through the MCS module.
!>

module perform_retrievals_mod

contains

  !> @param my_instrument Instrument entity
  subroutine perform_retrievals(my_instrument, CS)

    ! User modules
    use control_mod, only: CS_t
    use instruments_mod, only: generic_instrument
    use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error
    use guanter_model_mod
    use physical_model_mod
    use oco2_mod

    ! Third-party modules

    ! System modules
    use HDF5

    implicit none

    class(generic_instrument) :: my_instrument
    type(CS_t), intent(inout) :: CS

    ! Local variables
    character(len=*), parameter :: fname = "perform_retrieval" ! Function name for logging
    integer(hid_t) :: l1b_file_id ! L1b HDF file handler
    logical :: SoundingGeometry_exists ! Does the SoundingGeometry group exist?
    integer :: hdferr ! HDF error variable


    ! Open the L1B_HDF file - this will generally be required
    call h5fopen_f(CS%input%L1B_filename%chars(), H5F_ACC_RDONLY_F, l1b_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening L1B HDF file: " // trim(CS%input%L1B_filename%chars()))

    ! .. and store the file ID in the MCS
    CS%input%l1b_file_id = l1b_file_id

    select type(my_instrument)
    type is (oco2_instrument)
       ! We copy the SoundingGeometry group over to the results section, for
       ! easy analysis of the results later on. (This is a fairly small amount of data)

       call h5lexists_f(l1b_file_id, "/SoundingGeometry", SoundingGeometry_exists, hdferr)
       if (SoundingGeometry_exists) then
          call h5ocopy_f(l1b_file_id, "/SoundingGeometry", &
               CS%output%output_file_id, "/SoundingGeometry", hdferr)
          call check_hdf_error(hdferr, fname, "Error copying /SoundingGeometry into output file")
       else
          call logger%debug(fname, "/SoundingGeometry not part of this L1B file. Skipping.")
       end if
    end select


    if (CS%algorithm%using_GK_SIF) then
       ! Launch Guanter Retrievals
       call guanter_retrieval(my_instrument, CS)
    end if

    if (CS%algorithm%using_physical) then
       ! Launch physical-type retrievals
       call physical_retrieval(my_instrument, CS)
    end if

    ! Close L1B HDF file after we're done.
    call h5fclose_f(CS%input%l1b_file_id, hdferr)
    call check_hdf_error(hdferr, "Main", "Error closing input HDF5 file")

  end subroutine perform_retrievals

end module perform_retrievals_mod
