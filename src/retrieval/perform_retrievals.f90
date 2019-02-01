subroutine perform_retrievals(my_instrument)
!! In this subroutine, we essentially perform all the tasks for the retrievals,
!! including setting them up. All necessary information about the settings etc.,
!! is avaliable through the MCS module.

use control_mod, only: MCS
use instruments, only: generic_instrument
use logger_mod, only: logger => master_logger
use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error

use guanter_model_mod
use physical_model_mod
use HDF5

implicit none

class(generic_instrument) :: my_instrument

character(len=*), parameter :: fname = "perform_retrieval"
integer(hid_t) :: l1b_file_id

integer :: hdferr


! Open the L1B_HDF file - this will generally be required
call h5fopen_f(MCS%input%L1B_filename%chars(), H5F_ACC_RDONLY_F, l1b_file_id, hdferr)
call check_hdf_error(hdferr, fname, "Error opening L1B HDF file: " // trim(MCS%input%L1B_filename%chars()))

! .. and store the file ID in the MCS
MCS%input%l1b_file_id = l1b_file_id

! We copy the SoundingGeometry group over to the results section, for
! easy analysis of the results later on.

! TODO: this seems to case the HDF5 library to throw an "infinte loop"
! error on exiting the program. Doesn't seem to cause real problems, but
! looks 'ugly'. The alternative is to write a seperate function that does
! the copying 'bit by bit'. Checking this problem with valgrind makes me
! believe it's a problem in the HDF5 library version 1.10.4; looks like
! there's some memory allocation problems..

call h5ocopy_f(l1b_file_id, "/SoundingGeometry", &
               MCS%output%output_file_id, "/SoundingGeometry", hdferr)
call check_hdf_error(hdferr, fname, "Error copying /SoundingGeometry into output file")


if (MCS%algorithm%using_GK_SIF) then
    ! Launch Guanter Retrievals
    call guanter_retrieval(my_instrument)
end if

if (MCS%algorithm%using_physical) then
    ! Launch physical-type retrievals
    call physical_retrieval(my_instrument)
end if


end subroutine
