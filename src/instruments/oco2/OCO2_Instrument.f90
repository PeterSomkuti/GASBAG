!> @brief OCO-2(-like) specific subroutines
!> @file OCO2_Instrument.f90
!> @author Peter Somkuti
!>

module oco2_mod

  ! Modules
  use instruments_mod, only: generic_instrument
  use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, &
       read_DP_hdf_dataset, read_INT_hdf_dataset
  use control_mod, only: MCS

  ! Third-party modules
  use stringifor
  use logger_mod, only: logger => master_logger
  use mod_datetime

  ! System modules
  use HDF5
  use OMP_LIB

  implicit none

  !> @brief OCO-2 instrument type, which is an extention of the
  !> (empty) generic instrument type.
  type, extends(generic_instrument), public :: oco2_instrument
   contains
     procedure, nopass :: scan_l1b_file
     procedure, nopass :: read_l1b_dispersion
     procedure, nopass :: read_l1b_snr_coef
     procedure, nopass :: read_num_frames_and_fp
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
     procedure, nopass :: read_spike_filter
  end type oco2_instrument


contains

  subroutine scan_l1b_file(l1b_file)
    !! Checks whether this is indeed a valid OCO-2 l1b file that has all
    !! the required fields and variables..

    type(string) :: l1b_file

    character(len=*), parameter :: fname = "scan_l1b_file"
    character(len=999) :: msg
    integer(hid_t) :: file_id
    integer(8), allocatable :: n_fp_frames(:), dim_spec(:)
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
    MCS%general%N_frame = n_fp_frames(2)
    MCS%general%N_fp = n_fp_frames(1)

    ! OCO-2, we have three bands
    MCS%general%N_bands = 3

    ! And we can grab the number of pixels per band individually
    allocate(MCS%general%N_spec(MCS%general%N_bands))

    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_o2", dim_spec)
    MCS%general%N_spec(1) = dim_spec(1)
    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_weak_co2", dim_spec)
    MCS%general%N_spec(2) = dim_spec(1)
    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_strong_co2", dim_spec)
    MCS%general%N_spec(3) = dim_spec(1)

  end subroutine scan_l1b_file




  !> @brief Read the dispersion coefficients from the L1B file
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param dispersion_coeffs Unallocated array for dispersion coefficients (coef, fp, band)
  subroutine read_l1b_dispersion(l1b_file_id, dispersion_coeffs)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    double precision, allocatable, intent(out) :: dispersion_coeffs(:,:,:)

    integer(hsize_t), allocatable :: disp_shape(:)
    character(len=*), parameter :: fname = "read_l1b_dispersion(oco2)"
    character(len=*), parameter :: dset_name = "/InstrumentHeader/dispersion_coef_samp"

    ! Coeffs, footprints, bands
    call read_DP_hdf_dataset(l1b_file_id, dset_name, dispersion_coeffs, disp_shape)

  end subroutine read_l1b_dispersion

  !> @brief Read the dispersion coefficients from the L1B file
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param snr_coeffs Unallocated array for SNR coefficients (coef, pixel, fp, band)
  subroutine read_l1b_snr_coef(l1b_file_id, snr_coefs)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    double precision, allocatable, intent(out) :: snr_coefs(:,:,:,:)

    integer(hsize_t), allocatable :: snr_shape(:)
    character(len=*), parameter :: fname = "read_l1b_snr_coefs(oco2)"
    character(len=*), parameter :: dset_name = "/InstrumentHeader/snr_coef"

    ! coefficients, spectral indices, footprints, bands
    call read_DP_hdf_dataset(l1b_file_id, dset_name, snr_coefs, snr_shape)

  end subroutine read_l1b_snr_coef


  subroutine read_num_frames_and_fp(l1b_file_id, num_frames, num_fp)

    integer(hid_t), intent(in) :: l1b_file_id
    integer, intent(out) :: num_frames, num_fp

    integer(hsize_t), allocatable :: num_frames_shape(:)
    integer :: hdferr
    character(len=*), parameter :: fname = "read_num_frames(oco2)"
    character(len=*), parameter :: dset_name = "/SoundingGeometry/sounding_id"

    call get_HDF5_dset_dims(l1b_file_id, trim(dset_name), num_frames_shape)
    call check_hdf_error(hdferr, fname, "Error reading data set dimesions at: " // trim(dset_name))

    num_fp = num_frames_shape(1)
    num_frames = num_frames_shape(2)

  end subroutine read_num_frames_and_fp

  !> @brief Calculate the dispersion array, based on the dispersion coefficients
  !> @param disp_coef dispersion_coefficients (coef, fp, band)
  !> @param dispersion dispersion array (pixel)
  !> @param band Band number
  !> @param fp Footprint number
  subroutine calculate_dispersion(disp_coef, dispersion, band, fp)

    implicit none

    double precision, intent(in) :: disp_coef(:)
    integer, intent(in) :: band
    integer, intent(in) :: fp
    double precision, intent(out) :: dispersion(:)

    ! Local variables
    integer :: pix, order

    dispersion(:) = 0.0d0

    do pix=1, size(dispersion)
       do order=1, size(disp_coef)
          dispersion(pix) = dispersion(pix) + (dble(pix) ** (order-1) * disp_coef(order))
       end do
    end do

  end subroutine calculate_dispersion

  !> @brief Read one single spectrum from the L1B data
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param i_fr Frame index
  !> @param i_fp Footprint index
  !> @param band Band number
  !> @param N_spec Number of spectral points
  !> @param spectrum Output spectral array
  subroutine read_one_spectrum(l1b_file_id, i_fr, i_fp, band, N_spec, spectrum)

    ! Select one (!) OCO-2 sounding via HDF5 hyperslab selection and feed it
    ! into the array "spectrum". The chosen spectrum is selected by the indices
    ! i_fr (frame index), i_fp (footprint index), and band (number)

    implicit none

    integer(hid_t) :: l1b_file_id
    integer, intent(in) :: i_fr, i_fp, band, N_spec
    double precision, allocatable, intent(inout) :: spectrum(:)

    character(len=*), parameter :: fname = "read_one_spectrum(oco2)"
    character(len=999) :: dset_name
    integer(hid_t) :: dset_id, dspace_id, memspace_id
    integer(hsize_t) :: hs_offset(3), hs_count(3)
    integer(hsize_t) :: dim_mem(1)
    integer :: hdferr
    integer(kind = OMP_lock_kind) :: lck

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

    ! Offset - where do we start our hyperslab? We read the full spectrum, so
    ! the first index is 0, the other two depenend on the indices.
    ! Remember the order in OCO-2 files: (spectral index, footprint, frame)
    ! .. as seen by Fortran

    hs_offset(1) = 0
    hs_offset(2) = i_fp - 1
    hs_offset(3) = i_fr - 1

    ! We step 1016 in the spectral direction to get the full measurement,
    ! and 1 each in the frame and footprint directions (convention)
    hs_count(1) = N_spec
    hs_count(2) = 1
    hs_count(3) = 1

    ! This is the size of the container that we will be writing the spectral
    ! data into.
    dim_mem(1) = N_spec

    allocate(spectrum(N_spec))
    spectrum(:) = 0.0d0

    ! So this is how a hyperslab selection in HDF5 works. First, get the
    ! dataspace corresponding to the dataset you want to grab from. Then
    ! call h5sselect_hyperslab, and using offset and count, select the
    ! data you want to grab. Now we have to create a memory space which has
    ! the exact same size (not shape necessarily) as the hyperslab selection.
    ! Now using h5dread, using both memory and dataspace id, the data can
    ! be read from the file.

    ! If we have OpenMP, lock this section of the code, to make sure
    ! only one thread at a time is reading in a spectrum..

!$OMP CRITICAL
    call h5dopen_f(l1b_file_id, dset_name, dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening spectra at: " // trim(dset_name))

    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error getting dataspace id for " // trim(dset_name))

    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
         hs_offset, hs_count, hdferr)
    call check_hdf_error(hdferr, fname, "Error performing hyperslab selection.")

    call h5screate_simple_f(1, dim_mem, memspace_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error creating simple memory space.")

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, spectrum, dim_mem, &
         hdferr, memspace_id, dspace_id)
    call check_hdf_error(hdferr, fname, "Error reading spectrum data from " // trim(dset_name))
!$OMP END CRITICAL

  end subroutine read_one_spectrum


  !> @brief Calculate the noise-equivalent radiances from the SNR coefficients
  !> @param snr_coefs SNR coefficients (coef, pixel, fp, band)
  !> @param radiance L1B Radiance
  !> @param noise Noise-equivalent radiances
  !> @param fp Footprint number
  !> @param band Band number
  !> @param idx_start Starting index of radiance slice
  !> @param idx_end Last index of radiance slice (not used)
  subroutine calculate_noise(snr_coefs, radiance, noise, fp, band, idx_start, idx_end)

    implicit none

    double precision, intent(in) :: snr_coefs(:,:,:,:)
    double precision, intent(in) :: radiance(:)
    double precision, intent(inout) :: noise(:)
    integer, intent(in) :: fp, band, idx_start, idx_end

    double precision, allocatable :: MaxMS(:)
    integer(hsize_t), allocatable :: dims(:)
    integer :: i

!$OMP CRITICAL
    call read_DP_hdf_dataset(MCS%input%l1b_file_id, &
         "/InstrumentHeader/measureable_signal_max_observed", MaxMS, dims)
!$OMP END CRITICAL

    do i=1, size(noise)
       noise(i) = (MaxMS(band) / 100.0d0) * sqrt(abs(100.0d0 * radiance(i) / MaxMS(band)) * &
            (snr_coefs(1,idx_start+i-1,fp,band)**2) + (snr_coefs(2,idx_start+i-1,fp,band)**2))
    end do

  end subroutine calculate_noise

  !> @brief Check if the radiance array within the bounds that we're interested in
  !> is actually valid.
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param radiance Radiance array, extracted from the L1B
  !> @param idx_start Starting index of radiance slice to check
  !> @param idx_end Last index of radiance slice to check
  subroutine check_radiance_valid(l1b_file_id, radiance, idx_start, idx_end, valid)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    double precision, dimension(:), intent(in) :: radiance
    integer, intent(in) :: idx_start, idx_end
    logical, intent(out) :: valid

    valid = .true.

    ! For OCO-2, "bad" radiance points might be flagged by
    ! -999999 values.
    if (COUNT(radiance == -999999.0d0) > 0) then
       valid = .false.
    end if

  end subroutine check_radiance_valid

  !> @brief Read the whole sounding ID array
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param sounding_ids Sounding ID array (fp, frame)
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

  end subroutine read_sounding_ids


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

  end subroutine read_time_strings


  subroutine convert_time_string_to_date(time_string, date)

    ! Turn a OCO-2 time string into a datetype object

    implicit none
    character(len=25), intent(in) :: time_string
    type(datetime), intent(out) :: date
    integer :: year, month, day, hour, minute, second, millisecond

    ! Grab the various fields/positions from the string and stick them into
    ! the corresponding date/time variables
    read(time_string(1:4), *) year
    read(time_string(6:7), *) month
    read(time_string(9:10), *) day
    read(time_string(12:13), *) hour
    read(time_string(15:16), *) minute
    read(time_string(18:19), *) second
    read(time_string(21:23), *) millisecond

    ! Create datetime object
    date = datetime(year, month, day, hour, minute, &
         second, millisecond)

  end subroutine convert_time_string_to_date


  !> @brief Read variables related to the sounding scene geometry from L1B
  !> @param l1b_file_id HDF File ID of the L1B file
  !> @param band Band number
  !> @param SZA SZA array (fp, frame)
  !> @param SAA SAA array (fp, frame)
  !> @param VZA VZA array (fp, frame)
  !> @param VAA VAA array (fp, frame)
  subroutine read_sounding_geometry(l1b_file_id, band, SZA, SAA, VZA, VAA)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    integer, intent(in) :: band
    double precision, dimension(:,:), allocatable, intent(out) :: SZA, SAA, VZA, VAA

    character(len=*), parameter :: fname = "read_sounding_geometry(oco2)"
    integer(hsize_t), dimension(:), allocatable :: dset_dims
    double precision, dimension(:,:,:), allocatable :: tmp_array

    call logger%debug(fname, "Trying to allocate sounding location arrays.")

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

  end subroutine read_sounding_geometry

  subroutine read_ils_data(l1b_file_id, ils_delta_lambda, ils_relative_response)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    double precision, allocatable, intent(inout) :: ils_delta_lambda(:,:,:,:), &
         ils_relative_response(:,:,:,:)

    character(len=*), parameter :: fname = "read_sounding_location(oco2)"
    integer(hsize_t), dimension(:), allocatable :: dset_dims

    call read_DP_hdf_dataset(l1b_file_id, "InstrumentHeader/ils_delta_lambda", &
         ils_delta_lambda, dset_dims)
    call read_DP_hdf_dataset(l1b_file_id, "InstrumentHeader/ils_relative_response", &
         ils_relative_response, dset_dims)

  end subroutine read_ils_data


  subroutine read_sounding_location(l1b_file_id, band, lon, lat, altitude, rel_vel, rel_solar_vel)

    implicit none
    integer(hid_t), intent(in) :: l1b_file_id
    integer, intent(in) :: band
    double precision, dimension(:,:), allocatable, intent(out) :: lon, lat, &
         altitude, rel_vel, rel_solar_vel

    character(len=*), parameter :: fname = "read_sounding_location(oco2)"
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

    allocate(rel_vel(dset_dims(2), dset_dims(3)))
    allocate(rel_solar_vel(dset_dims(2), dset_dims(3)))

  end subroutine read_sounding_location

  subroutine read_bad_sample_list(l1b_file_id, bad_sample_list)
    implicit none
    integer(hid_t), intent(in) :: l1b_file_id
    integer, allocatable, intent(inout) :: bad_sample_list(:,:,:)

    integer(hsize_t), allocatable :: dset_dims(:)

    call read_INT_hdf_dataset(l1b_file_id, "InstrumentHeader/bad_sample_list", bad_sample_list, dset_dims)

  end subroutine read_bad_sample_list

  subroutine read_spike_filter(l1b_file_id, spike_list, band)
    implicit none
    integer(hid_t), intent(in) :: l1b_file_id
    integer, allocatable, intent(inout) :: spike_list(:,:,:)
    integer, intent(in) :: band

    integer(hsize_t), allocatable :: dset_dims(:)

    if (band == 1) then
       call read_INT_hdf_dataset(l1b_file_id, "SpikeEOF/spike_eof_weighted_residual_o2", &
            spike_list, dset_dims)
    elseif (band == 2) then
       call read_INT_hdf_dataset(l1b_file_id, "SpikeEOF/spike_eof_weighted_residual_weak_co2", &
            spike_list, dset_dims)
    elseif (band == 3) then
       call read_INT_hdf_dataset(l1b_file_id, "SpikeEOF/spike_eof_weighted_residual_strong_co2", &
            spike_list, dset_dims)
    else
       call logger%fatal("read_spike_filter(oco2)", "Sorry - I only know bands 1 through 3!")
       stop 1
    end if

  end subroutine read_spike_filter



end module oco2_mod
