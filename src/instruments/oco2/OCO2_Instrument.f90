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
        procedure, nopass :: calculate_dispersion
        procedure, nopass :: read_one_spectrum
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

    subroutine calculate_dispersion(disp_coef, dispersion)

        double precision, intent(in) :: disp_coef(:,:,:)
        double precision, intent(out), allocatable :: dispersion(:,:,:)

        integer(8) :: pix, fp, order, band

        allocate(dispersion(1016, 8, 3))
        dispersion(:,:,:) = 0.0d0

        do pix=1, 1016
            do band=1, 3
                do fp=1, 8
                    do order=1, 6
                        dispersion(pix,fp,band) = dispersion(pix,fp,band) + (pix-1) ** (order-1) * disp_coef(order,fp,band)
                    end do
                end do
            end do
        end do

    end subroutine

    subroutine read_one_spectrum(l1b_file_id, i_fr, i_fp, band, spectrum)

        ! Select one (!) OCO-2 sounding via HDF5 hyperslab selection and feed it
        ! into the array "spectrum". The chosen spectrum is selected by the indices
        ! i_fr (frame index), i_fp (footprint index), and band (number)

        implicit none

        integer(hid_t) :: l1b_file_id
        integer, intent(in) :: i_fr, i_fp, band
        double precision, allocatable :: spectrum(:)

        character(len=*), parameter :: fname = "read_one_spectrum(oco2)"
        character(len=999) :: dset_name
        integer(hid_t) :: dset_id, dspace_id, memspace_id
        logical :: selection_valid
        integer(hsize_t) :: hs_offset(3), hs_count(3)
        integer(hsize_t) :: dim_mem(1)
        integer :: hdferr
        logical :: extent_equal


        ! Set dataset name according to the band we want
        if (band == 1) then
            dset_name = "/SoundingMeasurements/radiance_o2"
        else if (band == 2) then
            dset_name = "/SoundingMeasurements/radiance_weak_co2"
        else if (band == 3) then
            dset_name = "/SoundingMeasurements/radiance_strong_co2"
        else
            call logger%fatal(fname, "Band number must be between 1 and 3!")
            stop 1
        end if

        call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error opening spectra at: " // trim(dset_name))
            stop 1
        end if

        !! Offset - where do we start our hyperslab? We read the full spectrum, so
        !! the first index is 0, the other two depenend on the indices.
        !! Remember the order in OCO-2 files: (spectral index, footprint, frame)
        !! .. as seen by Fortran

        hs_offset(1) = 0
        hs_offset(2) = i_fp - 1
        hs_offset(3) = i_fr - 1

        !! We step 1016 in the spectral direction to get the full measurement,
        !! and 1 each in the frame and footprint directions (convention)
        hs_count(1) = 1016
        hs_count(2) = 1
        hs_count(3) = 1

        !! This is the size of the container that we will be writing the spectral
        !! data into.
        dim_mem(1) = 1016

        allocate(spectrum(1016))
        spectrum(:) = 0.0d0

        !! So this is how a hyperslab selection in HDF5 works. First, get the
        !! dataspace corresponding to the dataset you want to grab from. Then
        !! call h5sselect_hyperslab, and using offset and count, select the
        !! data you want to grab. Now we have to create a memory space which has
        !! the exact same size (not shape necessarily) as the hyperslab selection.
        !! Now using h5dread, using both memory and dataspace id, the data can
        !! be read from the file.

        call h5dget_space_f(dset_id, dspace_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error getting dataspace id for " // trim(dset_name))
            stop 1
        end if

        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
                                   hs_offset, hs_count, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error performing hyperslab selection.")
            stop 1
        end if

        call h5screate_simple_f(1, dim_mem, memspace_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error creating simple memory space.")
            stop 1
        end if

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, spectrum, dim_mem, hdferr, memspace_id, dspace_id)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error reading spectrum data from " // trim(dset_name))
            stop 1
        end if

    end subroutine

end module
