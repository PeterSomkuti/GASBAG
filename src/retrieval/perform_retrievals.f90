subroutine perform_retrievals(my_instrument)
!! In this subroutine, we essentially perform all the tasks for the retrievals,
!! including setting them up. All necessary information about the settings etc.,
!! is avaliable through the MCS module.

use control, only: MCS
use instruments, only: generic_instrument
use logger_mod, only: logger => master_logger
use file_utils, only: get_HDF5_dset_dims

use guanter_model
use HDF5

implicit none

class(generic_instrument) :: my_instrument

character(len=*), parameter :: fname = "perform_retrieval"
integer(hid_t) :: l1b_file_id

integer :: hdferr

! Open the L1B_HDF file - this will generally be required
call h5fopen_f(MCS%input%L1B_filename%chars(), H5F_ACC_RDONLY_F, l1b_file_id, hdferr)
if (hdferr /= 0) then
  call logger%fatal(fname, "Error opening HDF file: " // trim(MCS%input%L1B_filename%chars()))
  stop 1
end if

MCS%input%l1b_file_id = l1b_file_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Basisfunction read-in from HDF file !!
!! These are used very often throughout the retrieval process, !!
!! hence we pre-load them (like ABSCO's) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (MCS%algorithm%using_GK_SIF) then
    ! Launch Guanter Retrievals!!!
    call guanter_retrieval(my_instrument)
end if

! Do the error analysis
! How do we do this without knowing first which forward models are being run?


! Save the results into a file


end subroutine
