subroutine perform_retrievals(my_instrument)
!! In this subroutine, we essentially perform all the tasks for the retrievals,
!! including setting them up. All necessary information about the settings etc.,
!! is avaliable through the MCS module.

use control, only: MCS
use instruments, only: generic_instrument
use logger_mod, only: logger => master_logger
use file_utils, only: get_HDF5_dset_dims
use HDF5

implicit none

class(generic_instrument) :: my_instrument

character(len=*), parameter :: fname = "perform_retrieval"
integer(hid_t) :: l1b_file_id, basisfunction_file_id
logical :: using_GK_SIF
integer :: hdferr

! SIF basis functions: (number, footprint, pixel)
double precision, allocatable :: basisfunctions(:,:,:)
double precision, allocatable :: tmp_data(:)
integer(hsize_t), allocatable :: num_pixels(:)

character(len=999) :: dset_name
integer(hid_t) :: dset_id

integer :: i,j,k


! Open the L1B_HDF file - this will generally be required
call h5fopen_f(MCS%input%L1B_filename%chars(), H5F_ACC_RDONLY_F, l1b_file_id, hdferr)
if (hdferr /= 0) then
  call logger%fatal(fname, "Error opening HDF file: " // trim(MCS%input%L1B_filename%chars()))
  stop 1
end if


! If required - open the basisfunction file and populate the basisfunction array.
! See first if we are using the Guanter-Koehler method?
using_GK_SIF = .false.
! For SIF and GK algorithm, we need radiance basis functions
do i=1, MCS%algorithm%N_algorithms
    if (MCS%algorithm%name(i) == "GK") then
        using_GK_SIF = .true.
    end if
end do

if (using_GK_SIF) then
    ! If so, let us open up the basisfunction file
    call h5fopen_f(MCS%window%basisfunction_file%chars(), H5F_ACC_RDONLY_F, &
                   basisfunction_file_id, hdferr)
    if (hdferr /= 0) then
        call logger%fatal(fname, "Error opening HDF file: " // trim(MCS%window%basisfunction_file%chars()))
        stop 1
    end if

    ! And then start reading values into our own arrays
    ! Get the array dimensions by inquiring the first one..
    call get_HDF5_dset_dims(basisfunction_file_id, "/BasisFunction_SV1_FP1", num_pixels)

    ! Allocate the basisfunction array
    allocate(basisfunctions(num_pixels(1), my_instrument%num_fp, MCS%algorithm%n_basisfunctions))
    allocate(tmp_data(num_pixels(1)))
    ! And read them in, one by one

    do i=1, my_instrument%num_fp
        do j=1, MCS%algorithm%n_basisfunctions

            write(dset_name, '(A,G0.1,A,G0.1)') "/BasisFunction_SV", j, "_FP", i
            call logger%debug(fname, "Looking for " // trim(dset_name))
            call h5dopen_f(basisfunction_file_id, dset_name, dset_id, hdferr)
            if (hdferr /= 0) then
                call logger%fatal(fname, "Error. Could not open " // trim(dset_name))
                stop 1
            end if

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_data, num_pixels, hdferr)
            if (hdferr /= 0) then
                call logger%fatal(fname, "Error. Could not read " // trim(dset_name))
                stop 1
            end if
            call logger%debug(fname, "Read in " // trim(dset_name))

            basisfunctions(:,i,j) = tmp_data(:)


        end do
    end do


    do i=1, num_pixels(1)
        write(*,*) basisfunctions(i,1,1)
    end do

end if






end subroutine
