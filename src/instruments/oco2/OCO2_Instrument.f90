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
        procedure, nopass :: read_l1b_dispersion
        procedure, nopass :: read_l1b_snr_coef
        procedure, nopass :: read_num_frames
    end type


contains

    subroutine scan_l1b_file(l1b_file)
        !! Checks whether this is indeed a valid OCO-2 l1b file that has all
        !! the required fields and variables..

        type(string) :: l1b_file

        character(len=*), parameter :: fname = "scan_l1b_file"
        character(len=999) :: msg
        integer(hid_t) :: file_id
        integer(kind=8), allocatable :: n_fp_frames(:)
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
        if (size(n_fp_frames) /= 2) then
            call logger%fatal(fname, "This array -n_fp_frames- should be of size 2. But it isn't.")
            stop 1
        end if

        write(msg, "(A, I1, A, I8.1)") "Number of footprints: ", n_fp_frames(1), ", number of frames: ", n_fp_frames(2)
        call logger%info(fname, trim(msg))
        write(msg, "(A, I8.1, A)") "For a total of ", n_fp_frames(1)*n_fp_frames(2), " soundings."
        call logger%info(fname, trim(msg))

        ! Store the total number of soundings to be processed in the MCS. We need
        ! that later to allocate all those big arrays.
        MCS%general%N_soundings = n_fp_frames(1)*n_fp_frames(2)
    end subroutine


    subroutine read_l1b_dispersion(l1b_file_id, dispersion_coeffs)

        integer(hid_t), intent(in) :: l1b_file_id
        double precision, allocatable, intent(out) :: dispersion_coeffs(:,:,:)

        ! We hard-code the dimensions of the OCO-2 dispersion coefficients here,
        ! we do not really expect them to change do we?
        integer(hsize_t) :: disp_shape(3) = [6,8,3]
        character(len=*), parameter :: fname = "read_l1b_dispersion(oco2)"
        character(len=*), parameter :: dset_name = "/InstrumentHeader/dispersion_coef_samp"
        integer :: hdferr
        integer(hid_t) :: dset_id

        ! Coeffs, footprints, bands
        allocate(dispersion_coeffs(6,8,3))

        call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error opening dispersion coeffs at: " // trim(dset_name))
            stop 1
        end if

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dispersion_coeffs, disp_shape, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error reading dispersion coeffs at: " // trim(dset_name))
            stop 1
        end if

    end subroutine

    subroutine read_l1b_snr_coef(l1b_file_id, snr_coefs)

        integer(hid_t), intent(in) :: l1b_file_id
        double precision, allocatable, intent(out) :: snr_coefs(:,:,:,:)

        integer(hsize_t) :: snr_shape(4) = [2,1016,8,3]
        character(len=*), parameter :: fname = "read_l1b_snr_coefs(oco2)"
        character(len=*), parameter :: dset_name = "/InstrumentHeader/snr_coef"
        integer :: hdferr
        integer(hid_t) :: dset_id

        ! coefficients, spectral indices, footprints, bands
        allocate(snr_coefs(2,1016,8,3))

        call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error opening SNR coeffs at: " // trim(dset_name))
            stop 1
        end if

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, snr_coefs, snr_shape, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error reading SNR coeffs at: " // trim(dset_name))
            stop 1
        end if

    end subroutine


    subroutine read_num_frames(l1b_file_id, num_frames)

        integer(hid_t), intent(in) :: l1b_file_id
        integer, intent(out) :: num_frames(1)

        integer(hid_t) :: dset_id
        integer(hsize_t) :: num_frames_shape(1) = [1]
        integer :: hdferr
        character(len=*), parameter :: fname = "read_num_frames(oco2)"
        character(len=*), parameter :: dset_name = "/Metadata/ExpectedFrames"



        call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error opening SNR coeffs at: " // trim(dset_name))
            stop 1
        end if

        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_frames, num_frames_shape, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error reading SNR coeffs at: " // trim(dset_name))
            stop 1
        end if

    end subroutine

end module
