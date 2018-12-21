subroutine perform_retrievals(my_instrument)
!! In this subroutine, we essentially perform all the tasks for the retrievals,
!! including setting them up. All necessary information about the settings etc.,
!! is avaliable through the MCS module.

use control, only: MCS
use instruments, only: generic_instrument
use logger_mod, only: logger => master_logger
use file_utils, only: get_HDF5_dset_dims, check_hdf_error

use guanter_model
use physical_model
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


if (MCS%algorithm%using_GK_SIF) then
    ! Launch Guanter Retrievals!!!
    call guanter_retrieval(my_instrument)
end if

if (MCS%algorithm%using_physical) then
    ! Launch physical-type retrievals
    call physical_retrieval(my_instrument)
end if


end subroutine
