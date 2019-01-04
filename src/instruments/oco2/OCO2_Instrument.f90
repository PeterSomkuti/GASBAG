module oco2_mod

    use stringifor
    use instruments, only: generic_instrument
    use control_mod, only: MCS
    use logger_mod, only: logger => master_logger
    use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, &
                              read_DP_hdf_dataset, read_INT_hdf_dataset
    use mod_datetime

    use HDF5

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
        procedure, nopass :: calculate_noise
        procedure, nopass :: check_radiance_valid
        procedure, nopass :: read_sounding_ids
        procedure, nopass :: read_time_strings
        procedure, nopass :: convert_time_string_to_date
        procedure, nopass :: read_ils_data
        procedure, nopass :: read_sounding_geometry
        procedure, nopass :: read_sounding_location
        procedure, nopass :: read_bad_sample_list
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
        call check_hdf_error(hdferr, fname, "Error opening HDF file: " // trim(l1b_file%chars()))

        ! Let's start with the sounding IDs and have a look how many we actually
        ! have in this file.
        call get_HDF5_dset_dims(file_id, "/SoundingGeometry/sounding_id", n_fp_frames)
        if (size(n_fp_frames) /= 2) then
            call logger%fatal(fname, "This array -n_fp_frames- should be of size 2. But it isn't.")
            stop 1
        end if

        write(msg, "(A, G0.1, A, G0.1)") "Number of footprints: ", n_fp_frames(1), &
                                         ", number of frames: ", n_fp_frames(2)
        call logger%info(fname, trim(msg))
        write(msg, "(A, G0.1, A)") "For a total of ", n_fp_frames(1)*n_fp_frames(2), " soundings."
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
        call check_hdf_error(hdferr, fname, "Error opening dispersion coeffs at: " // trim(dset_name))

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dispersion_coeffs, disp_shape, hdferr)
        call check_hdf_error(hdferr, fname, "Error reading dispersion coeffs at: " // trim(dset_name))

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
        call check_hdf_error(hdferr, fname, "Error opening SNR coeffs at: " // trim(dset_name))

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, snr_coefs, snr_shape, hdferr)
        call check_hdf_error(hdferr, fname, "Error reading SNR coeffs at: " // trim(dset_name))

    end subroutine


    subroutine read_num_frames(l1b_file_id, num_frames)

        integer(hid_t), intent(in) :: l1b_file_id
        integer, intent(out) :: num_frames
        integer, dimension(1) :: tmp_num_frames

        integer(hid_t) :: dset_id
        integer(hsize_t) :: num_frames_shape(1) = [1]
        integer :: hdferr
        character(len=*), parameter :: fname = "read_num_frames(oco2)"
        character(len=*), parameter :: dset_name = "/Metadata/ExpectedFrames"

        call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error opening SNR coeffs at: " // trim(dset_name))

        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, tmp_num_frames, num_frames_shape, hdferr)
        call check_hdf_error(hdferr, fname, "Error reading SNR coeffs at: " // trim(dset_name))

        num_frames = tmp_num_frames(1)

    end subroutine

    subroutine calculate_dispersion(disp_coef, dispersion, band, fp)

        double precision, intent(in) :: disp_coef(:)
        double precision, intent(out) :: dispersion(:)
        integer, intent(in) :: band, fp
        integer :: pix, order

        dispersion(:) = 0.0d0

        do pix=1, 1016
            do order=1, 6
                dispersion(pix) = dispersion(pix) + (dble(pix) ** (order-1) * disp_coef(order))
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
        call check_hdf_error(hdferr, fname, "Error opening spectra at: " // trim(dset_name))

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
        call check_hdf_error(hdferr, fname, "Error getting dataspace id for " // trim(dset_name))

        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
                                   hs_offset, hs_count, hdferr)
        call check_hdf_error(hdferr, fname, "Error performing hyperslab selection.")


        call h5screate_simple_f(1, dim_mem, memspace_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error creating simple memory space.")

        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, spectrum, dim_mem, hdferr, memspace_id, dspace_id)
        call check_hdf_error(hdferr, fname, "Error reading spectrum data from " // trim(dset_name))


    end subroutine


    subroutine calculate_noise(snr_coefs, radiance, noise, fp, band, idx_start, idx_end)
        double precision, intent(in) :: snr_coefs(:,:,:,:)
        double precision, intent(in) :: radiance(:)
        double precision, intent(inout) :: noise(:)
        integer, intent(in) :: fp, band, idx_start, idx_end

        double precision, allocatable :: MaxMS(:)
        integer(hsize_t), allocatable :: dims(:)
        integer(hid_t) :: dset_id
        integer :: i, hdferr

        !dims(1) = 3
        !call h5dopen_f(MCS%input%l1b_file_id, "/InstrumentHeader/measureable_signal_max_observed", dset_id, hdferr)
        !call check_hdf_error(hdferr, "calculate_noise", "Error opening: /InstrumentHeader/measureable_signal_max_observed")

        !call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, MaxMS, dims, hdferr)
        !call read_DP_hdf_dataset()


        call read_DP_hdf_dataset(MCS%input%l1b_file_id, "/InstrumentHeader/measureable_signal_max_observed", MaxMS, dims)

        do i=1, size(noise)
            noise(i) = (MaxMS(band) / 100) * sqrt(abs(100 * radiance(i) / MaxMS(band)) * &
                        (snr_coefs(1,idx_start+i-1,fp,band)**2) + (snr_coefs(2,idx_start+i-1,fp,band)**2))
        end do

    end subroutine

    subroutine check_radiance_valid(l1b_file_id, radiance, idx_start, idx_end, valid)
        integer(hid_t), intent(in) :: l1b_file_id
        double precision, dimension(:), intent(in) :: radiance
        integer, intent(in) :: idx_start, idx_end
        logical, intent(out) :: valid


        valid = .true.

        if (COUNT(radiance == -999999.0) > 0) then
            valid = .false.
        end if

    end subroutine

    subroutine read_sounding_ids(l1b_file_id, sounding_ids)

        implicit none
        integer(hid_t), intent(in) :: l1b_file_id
        integer(8), dimension(:,:), allocatable, intent(out) :: sounding_ids

        character(len=*), parameter :: fname = "read_sounding_ids"
        integer(hid_t) :: dset_id, filetype
        integer(hid_t), dimension(:), allocatable :: dset_dims
        integer :: hdferr

        call h5dopen_f(l1b_file_id, "/SoundingGeometry/sounding_id", dset_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error opening: /SoundingGeometry/sounding_id")

        call get_HDF5_dset_dims(l1b_file_id, "/SoundingGeometry/sounding_id", dset_dims)

        allocate(sounding_ids(dset_dims(1), dset_dims(2)))
        ! sounding id's are 64bit little endian, but let's figure it out via this
        ! function anyway
        call h5dget_type_f(dset_id, filetype, hdferr)

        call h5dread_f(dset_id, filetype, sounding_ids, dset_dims, hdferr)
        call check_hdf_error(hdferr, fname, "Error reading in: /SoundingGeometry/sounding_id")

    end subroutine


    subroutine read_time_strings(l1b_file_id, time_strings)

        implicit none
        integer(hid_t), intent(in) :: l1b_file_id
        character(len=25), dimension(:,:), allocatable, intent(out) :: time_strings
        ! OCO-2 time strings have 24 characters!
        character(len=*), parameter :: fname = "read_time_strings"
        integer :: hdferr
        integer(hid_t) :: dset_id, filetype
        integer(hid_t), dimension(:), allocatable :: dset_dims

        call h5dopen_f(l1b_file_id, "/SoundingGeometry/sounding_time_string", dset_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error opening: /SoundingGeometry/sounding_time_string")

        call get_HDF5_dset_dims(l1b_file_id, "/SoundingGeometry/sounding_time_string", dset_dims)

        ! Difficult to figure out which kind of character-type we have in the HDF
        ! file, so let's just grab it.
        call h5dget_type_f(dset_id, filetype, hdferr)

        allocate(time_strings(dset_dims(1), dset_dims(2)))
        call h5dread_f(dset_id, filetype, time_strings, dset_dims, hdferr)
        call check_hdf_error(hdferr, fname, "Error reading in: /SoundingGeometry/sounding_time_string")

    end subroutine


    subroutine convert_time_string_to_date(time_string, date)

        ! Turn a OCO-2 time string into a datetype object

        implicit none
        character(len=25), intent(in) :: time_string
        type(datetime), intent(out) :: date

        type(string) :: tmp_str
        integer :: year, month, day, hour, minute, second, millisecond

        ! Grab the various fields/positions from the string
        read(time_string(1:4), *) year
        read(time_string(6:7), *) month
        read(time_string(9:10), *) day
        read(time_string(12:13), *) hour
        read(time_string(15:16), *) minute
        read(time_string(18:19), *) second
        read(time_string(21:23), *) millisecond

        ! Create datetime objection
        date = datetime(year, month, day, hour, minute, &
                        second, millisecond)

    end subroutine


    subroutine read_sounding_geometry(l1b_file_id, band, SZA, SAA, VZA, VAA)
        integer(hid_t), intent(in) :: l1b_file_id
        integer, intent(in) :: band
        double precision, dimension(:,:), allocatable, intent(out) :: SZA, SAA, VZA, VAA

        character(len=*), parameter :: fname = "read_sounding_geometry"
        integer :: hdferr
        integer(hid_t) :: dset_id
        integer(hsize_t), dimension(:), allocatable :: dset_dims
        double precision, dimension(:,:,:), allocatable :: tmp_array


        ! FootprintGeometry fields are (Band, FP, Frame)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_solar_zenith", tmp_array, dset_dims)
        allocate(SZA(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        SZA(:,:) = tmp_array(band,:,:)

        deallocate(tmp_array)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_solar_azimuth", tmp_array, dset_dims)
        allocate(SAA(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        SAA(:,:) = tmp_array(band,:,:)

        deallocate(tmp_array)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_zenith", tmp_array, dset_dims)
        allocate(VZA(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        VZA(:,:) = tmp_array(band,:,:)

        deallocate(tmp_array)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_azimuth", tmp_array, dset_dims)
        allocate(VAA(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        VAA(:,:) = tmp_array(band,:,:)


    end subroutine

    subroutine read_ils_data(l1b_file_id, ils_delta_lambda, ils_relative_response)

        implicit none
        integer(hid_t), intent(in) :: l1b_file_id
        double precision, allocatable, intent(inout) :: ils_delta_lambda(:,:,:,:), &
                                                        ils_relative_response(:,:,:,:)


        character(len=*), parameter :: fname = "read_sounding_location"
        integer(hsize_t), dimension(:), allocatable :: dset_dims

        call read_DP_hdf_dataset(l1b_file_id, "InstrumentHeader/ils_delta_lambda", &
                                 ils_delta_lambda, dset_dims)
        call read_DP_hdf_dataset(l1b_file_id, "InstrumentHeader/ils_relative_response", &
                                 ils_relative_response, dset_dims)

    end subroutine


    subroutine read_sounding_location(l1b_file_id, band, lon, lat, altitude, rel_vel, rel_solar_vel)

        implicit none
        integer(hid_t), intent(in) :: l1b_file_id
        integer, intent(in) :: band
        double precision, dimension(:,:), allocatable, intent(out) :: lon, lat, &
                                                altitude, rel_vel, rel_solar_vel

        character(len=*), parameter :: fname = "read_sounding_location"
        integer :: hdferr
        integer(hid_t) :: dset_id, filetype
        integer(hsize_t), dimension(:), allocatable :: dset_dims
        double precision, dimension(:,:,:), allocatable :: tmp_array


        ! FootprintGeometry fields are (Band, FP, Frame)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_longitude", tmp_array, dset_dims)
        allocate(lon(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        lon(:,:) = tmp_array(band,:,:)

        deallocate(tmp_array)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_latitude", tmp_array, dset_dims)
        allocate(lat(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        lat(:,:) = tmp_array(band,:,:)

        deallocate(tmp_array)
        call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_altitude", tmp_array, dset_dims)
        allocate(altitude(dset_dims(2), dset_dims(3))) ! We only want FP and Frame
        altitude(:,:) = tmp_array(band,:,:)

        call read_DP_hdf_dataset(l1b_file_id, "SoundingGeometry/sounding_relative_velocity", rel_vel, dset_dims)
        call read_DP_hdf_dataset(l1b_file_id, "SoundingGeometry/sounding_solar_relative_velocity", rel_solar_vel, dset_dims)

    end subroutine

    subroutine read_bad_sample_list(l1b_file_id, bad_sample_list)
        implicit none
        integer(hid_t), intent(in) :: l1b_file_id
        integer, allocatable, intent(inout) :: bad_sample_list(:,:,:)

        integer(hsize_t), allocatable :: dset_dims(:)

        call read_INT_hdf_dataset(l1b_file_id, "InstrumentHeader/bad_sample_list", bad_sample_list, dset_dims)

    end subroutine



end module
