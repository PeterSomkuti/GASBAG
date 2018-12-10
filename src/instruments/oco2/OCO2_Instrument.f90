module oco2

    use stringifor
    use instruments, only: generic_instrument
    use control, only: MCS
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims

    use HDF5
    use H5LT

    implicit none

    type, extends(generic_instrument), public :: oco2_instrument

        ! dispersion coefficients, noise model, ILS array, ...
        ! all of that should be read in immediately at the start of the
        ! program, and kept in memory for quick access

    contains
        procedure, nopass :: scan_l1b_file
    end type


contains

    subroutine scan_l1b_file(l1b_file)
        !! Checks whether this is indeed a valid OCO-2 l1b file that has all
        !! the required fields and variables..

        type(string) :: l1b_file

        character(len=*), parameter :: fname = "scan_l1b_file"
        character(len=999) :: msg
        integer(hid_t) :: file_id
        integer(kind=8) :: n_fp_frames(2)
        integer :: hdferr

        ! Open the HDF file
        call h5fopen_f(l1b_file%chars(), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then
          call logger%fatal(fname, "Error opening HDF file: " // trim(l1b_file%chars()))
          stop 1
        end if

        ! Let's start with the sounding IDs and have a look how many we actually
        ! have in this file.
        call get_HDF5_dset_dims(file_id, "/SoundingGeometry/sounding_id", n_fp_frames)

        write(msg, "(A, I1, A, I8.1)") "Number of footprints: ", n_fp_frames(1), ", number of frames: ", n_fp_frames(2)
        call logger%info(fname, trim(msg))
        write(msg, "(A, I8.1, A)") "For a total of ", n_fp_frames(1)*n_fp_frames(2), " soundings."
        call logger%info(fname, trim(msg))

        ! Store the total number of soundings to be processed in the MCS. We need
        ! that later to allocate all those big arrays.
        MCS%general%N_soundings = n_fp_frames(1)*n_fp_frames(2)

        



    end subroutine
end module
