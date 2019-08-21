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
  use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  implicit none

  !> @brief OCO-2 instrument type, which is an extention of the
  !> (empty) generic instrument type.
  type, extends(generic_instrument), public :: oco2_instrument
   contains
     procedure, nopass :: scan_l1b_file

     procedure, nopass :: calculate_dispersion
     procedure, nopass :: calculate_noise
     procedure, nopass :: check_radiance_valid
     procedure, nopass :: convert_time_string_to_date

     procedure, nopass :: read_sounding_ids
     procedure, nopass :: read_time_strings
     procedure, nopass :: read_l1b_dispersion
     procedure, nopass :: read_l1b_snr_coef
     procedure, nopass :: read_ils_data
     procedure, nopass :: read_one_spectrum
     procedure, nopass :: read_spectra_and_average_by_fp
     procedure, nopass :: read_sounding_geometry
     procedure, nopass :: read_sounding_location
     procedure, nopass :: read_bad_sample_list
     procedure, nopass :: read_spike_filter
     procedure, nopass :: read_stokes_coef
     procedure, nopass :: read_MET_data
  end type oco2_instrument


contains

  !> @brief Scans the L1B file and extracts some global information
  !> @param l1b_file Path to L1B file
  subroutine scan_l1b_file(l1b_file)
    !! Checks whether this is indeed a valid OCO-2 l1b file that has all
    !! the required fields and variables..

    type(string) :: l1b_file

    character(len=*), parameter :: fname = "oco2_scan_l1b_file"
    character(len=999) :: msg
    integer(hid_t) :: file_id
    integer(8), allocatable :: n_fp_frames(:), dim_spec(:)
    logical :: sounding_id_exists, radiance_o2_exists
    integer :: hdferr
    integer :: n_fp, n_frames

    ! Open the HDF file
    call h5fopen_f(l1b_file%chars(), H5F_ACC_RDONLY_F, file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening HDF file: " // trim(l1b_file%chars()))

    ! Let's start with the sounding IDs and have a look how many frames and footprints
    ! we actually have in this file. This only works for downlooking L1B files. Solar
    ! limb scans do not have the sounding_id field.

    call h5lexists_f(file_id, "/SoundingGeometry/sounding_id", sounding_id_exists, hdferr)
    call h5lexists_f(file_id, "/SoundingMeasurements/radiance_o2", radiance_o2_exists, hdferr)

    if (sounding_id_exists) then
       call get_HDF5_dset_dims(file_id, "/SoundingGeometry/sounding_id", n_fp_frames)
       if (size(n_fp_frames) /= 2) then
          call logger%fatal(fname, "This array -n_fp_frames- should be of size 2. But it isn't.")
          stop 1
       end if

       n_fp = n_fp_frames(1)
       n_frames = n_fp_frames(2)
    else if (radiance_o2_exists) then
       call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_o2", n_fp_frames)
       if (size(n_fp_frames) /= 3) then
          call logger%fatal(fname, "This array -n_fp_frames- should be of size 3. But it isn't.")
          stop 1
       end if

       n_fp = n_fp_frames(2)
       n_frames = n_fp_frames(3)
    else
       call logger%fatal(fname, "Error in determining the frame/footprint file structure.")
       stop 1
    end if


    ! Let the user know how many
    write(msg, "(A, G0.1, A, G0.1)") "Number of footprints: ", n_fp, &
         ", number of frames: ", n_frames
    call logger%info(fname, trim(msg))
    write(msg, "(A, G0.1, A)") "For a total of ", n_fp * n_frames, " soundings."
    call logger%info(fname, trim(msg))

    ! Store the total number of soundings to be processed in the MCS. We need
    ! that later to allocate all those big arrays.
    MCS%general%N_soundings = n_fp * n_frames
    MCS%general%N_frame = n_frames
    MCS%general%N_fp = n_fp

    call get_HDF5_dset_dims(file_id, "/InstrumentHeader/measureable_signal_max_observed", dim_spec)
    ! OCO-2, we have three bands, but we grab the value from this data instead
    MCS%general%N_bands = dim_spec(1)

    ! And we can grab the number of pixels per band individually
    allocate(MCS%general%N_spec(MCS%general%N_bands))
    MCS%general%N_spec(:) = 1

    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_o2", dim_spec)
    MCS%general%N_spec(1) = dim_spec(1)
    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_weak_co2", dim_spec)
    MCS%general%N_spec(2) = dim_spec(1)
    call get_HDF5_dset_dims(file_id, "/SoundingMeasurements/radiance_strong_co2", dim_spec)
    MCS%general%N_spec(3) = dim_spec(1)

    ! Close this instance, as it will be loaded by "perform_retrievals" again
    call h5fclose_f(file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error closing HDF file: " // trim(l1b_file%chars()))

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

    ! We step N_spec in the spectral direction to get the full measurement,
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
    ! only one thread at a time is reading in a spectrum.
    ! NOTE: HDF5 (without the parallel) is not designed for concurrent access,
    ! so we MUST restrict reading from a file to one thread/process at a time.

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

    call h5dclose_f(dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error closing dataset id for " // trim(dset_name))
!$OMP END CRITICAL

  end subroutine read_one_spectrum


  subroutine read_spectra_and_average_by_fp(l1b_file_id, i_fp, band, N_spec, radiance_l1b)

    integer(hid_t) :: l1b_file_id
    integer, intent(in) :: i_fp, band, N_spec
    double precision, allocatable, intent(inout) :: radiance_l1b(:)

    character(len=*), parameter :: fname = "read_spectra_and_average_by_fp(oco2)"
    double precision, allocatable :: tmp_radiance(:)
    integer :: i_fr, i, funit
    logical :: radiance_OK
    character(len=50) :: tmp_str

    allocate(radiance_l1b(N_spec))
    radiance_l1b(:) = 0.0d0

    do i_fr=1, MCS%general%N_frame
       call read_one_spectrum(l1b_file_id, i_fr, i_fp, band, &
            MCS%general%N_spec(band), tmp_radiance)

       call check_radiance_valid(l1b_file_id, tmp_radiance, &
            1, size(tmp_radiance), radiance_OK)

       if (i_fr == 1) then
          radiance_l1b(:) = tmp_radiance(:)
          deallocate(tmp_radiance)
          cycle
       end if

       if (radiance_OK) then
          radiance_l1b(:) = (radiance_l1b(:) * dble(i_fr - 1)) / dble(i_fr) &
               + (tmp_radiance(:) / i_fr)
       else
          cycle
       end if

       deallocate(tmp_radiance)

    end do

  end subroutine read_spectra_and_average_by_fp


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

    if (any(ieee_is_nan(radiance))) then
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
    character(len=25), dimension(:), allocatable, intent(out) :: time_strings
    ! OCO-2 time strings have 24 characters!
    character(len=*), parameter :: fname = "read_time_strings"
    integer :: hdferr
    integer(hid_t) :: dset_id, filetype
    integer(hid_t), dimension(:), allocatable :: dset_dims

    call h5dopen_f(l1b_file_id, "/FrameHeader/frame_time_string", dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening: /FrameHeader/frame_time_string")

    call get_HDF5_dset_dims(l1b_file_id, "/FrameHeader/frame_time_string", dset_dims)

    ! Difficult to figure out which kind of character-type we have in the HDF
    ! file, so let's just grab it.
    call h5dget_type_f(dset_id, filetype, hdferr)

    allocate(time_strings(dset_dims(1)))
    call h5dread_f(dset_id, filetype, time_strings, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error reading in: /FrameHeader/frame_time_string")

  end subroutine read_time_strings


  subroutine convert_time_string_to_date(time_string, date, success)

    ! Turn a OCO-2 time string into a datetype object

    implicit none
    character(len=25), intent(in) :: time_string
    type(datetime), intent(out) :: date
    logical, intent(inout) :: success
    integer :: iostat
    integer :: year, month, day, hour, minute, second, millisecond

    success = .false.
    ! Grab the various fields/positions from the string and stick them into
    ! the corresponding date/time variables
    read(time_string(1:4), *, iostat=iostat) year
    if (iostat /= 0) return
    read(time_string(6:7), *, iostat=iostat) month
    if (iostat /= 0) return
    read(time_string(9:10), *, iostat=iostat) day
    if (iostat /= 0) return
    read(time_string(12:13), *, iostat=iostat) hour
    if (iostat /= 0) return
    read(time_string(15:16), *, iostat=iostat) minute
    if (iostat /= 0) return
    read(time_string(18:19), *, iostat=iostat) second
    if (iostat /= 0) return
    read(time_string(21:23), *, iostat=iostat) millisecond
    if (iostat /= 0) return

    ! Create datetime object
    date = datetime(year, month, day, hour, minute, &
         second, millisecond)

    success = .true.

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
    double precision, dimension(:,:,:), allocatable :: tmp_array3d
    double precision, dimension(:), allocatable :: tmp_array1d
    integer :: i

    call logger%debug(fname, "Trying to allocate sounding location arrays.")
    allocate(SZA(MCS%general%N_fp, MCS%general%N_frame))
    allocate(SAA(MCS%general%N_fp, MCS%general%N_frame))
    allocate(VZA(MCS%general%N_fp, MCS%general%N_frame))
    allocate(VAA(MCS%general%N_fp, MCS%general%N_frame))

    if (MCS%algorithm%observation_mode == "downlooking") then

       ! FootprintGeometry fields are (Band, FP, Frame)
       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_solar_zenith", &
            tmp_array3d, dset_dims)
       SZA(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_solar_azimuth", &
            tmp_array3d, dset_dims)
       SAA(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_zenith", &
            tmp_array3d, dset_dims)
       VZA(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_azimuth", &
            tmp_array3d, dset_dims)
       VAA(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

    else if (MCS%algorithm%observation_mode == "space_solar") then

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/solar_zenith", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          SZA(i, :) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/boresight_zenith", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          VZA(i, :) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/solar_azimuth", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          SAA(i, :) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/boresight_azimuth", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          VAA(i, :) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

    end if

  end subroutine read_sounding_geometry

  subroutine read_stokes_coef(l1b_file_id, band, stokes_coefs)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    integer, intent(in) :: band
    double precision, dimension(:,:,:), allocatable, intent(out) :: stokes_coefs

    character(len=*), parameter :: fname = "read_stokes_coefs(oco2)"
    integer(hsize_t), dimension(:), allocatable :: dset_dims
    double precision, dimension(:,:,:,:), allocatable :: tmp_array4d
    double precision, dimension(:,:,:,:,:), allocatable :: tmp_array5d
    integer :: i

    call logger%debug(fname, "Trying to allocate stokes coef. array.")
    allocate(stokes_coefs(4, MCS%general%N_fp, MCS%general%N_frame))

    ! Hack - GeoCarb files have an extra dimension here. So let's check for that first
    call get_HDF5_dset_dims(l1b_file_id, "FootprintGeometry/footprint_stokes_coefficients", &
         dset_dims)

    if (size(dset_dims) == 4) then

       deallocate(dset_dims)
       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_stokes_coefficients", &
            tmp_array4d, dset_dims)
       stokes_coefs(:,:,:) = tmp_array4d(:,band,:,:)
       deallocate(tmp_array4d)

    else if (size(dset_dims) == 5) then

       deallocate(dset_dims)
       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_stokes_coefficients", &
            tmp_array5d, dset_dims)
       stokes_coefs(:,:,:) = tmp_array5d(1,:,band,:,:)
       deallocate(tmp_array5d)

    end if

  end subroutine read_stokes_coef

  subroutine read_ils_data(l1b_file_id, ils_delta_lambda, ils_relative_response)

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    double precision, allocatable, intent(inout) :: ils_delta_lambda(:,:,:,:), &
         ils_relative_response(:,:,:,:)

    character(len=*), parameter :: fname = "read_ils_data(oco2)"
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
    double precision, dimension(:,:,:), allocatable :: tmp_array3d
    double precision, dimension(:), allocatable :: tmp_array1d
    integer :: i

    allocate(lon(MCS%general%N_fp, MCS%general%N_frame))
    allocate(lat(MCS%general%N_fp, MCS%general%N_frame))
    allocate(altitude(MCS%general%N_fp, MCS%general%N_frame))

    if (MCS%algorithm%observation_mode == "downlooking") then

       ! FootprintGeometry fields are (Band, FP, Frame)
       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_longitude", &
            tmp_array3d, dset_dims)
       lon(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_latitude", &
            tmp_array3d, dset_dims)
       lat(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

       call read_DP_hdf_dataset(l1b_file_id, "FootprintGeometry/footprint_altitude", &
            tmp_array3d, dset_dims)
       altitude(:,:) = tmp_array3d(band,:,:)
       deallocate(tmp_array3d)

    else if (MCS%algorithm%observation_mode == "space_solar") then

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/spacecraft_lon", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          lon(i,:) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/spacecraft_lat", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          lat(i,:) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

       call read_DP_hdf_dataset(l1b_file_id, "/SpacePointingFrameGeometry/spacecraft_alt", &
            tmp_array1d, dset_dims)
       do i=1, MCS%general%N_fp
          altitude(i,:) = tmp_array1d(:)
       end do
       deallocate(tmp_array1d)

    end if

    ! Keep these guys at zero for the time being?
    allocate(rel_vel(MCS%general%N_fp, MCS%general%N_frame))
    rel_vel(:,:) = 0.0d0
    allocate(rel_solar_vel(MCS%general%N_fp, MCS%general%N_frame))
    rel_solar_vel(:,:) = 0.0d0

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
       call read_INT_hdf_dataset(l1b_file_id, "/SpikeEOF/spike_eof_weighted_residual_o2", &
            spike_list, dset_dims)
    elseif (band == 2) then
       call read_INT_hdf_dataset(l1b_file_id, "/SpikeEOF/spike_eof_weighted_residual_weak_co2", &
            spike_list, dset_dims)
    elseif (band == 3) then
       call read_INT_hdf_dataset(l1b_file_id, "/SpikeEOF/spike_eof_weighted_residual_strong_co2", &
            spike_list, dset_dims)
    else
       call logger%fatal("read_spike_filter(oco2)", "Sorry - I only know bands 1 through 3!")
       stop 1
    end if

  end subroutine read_spike_filter


  subroutine read_MET_data(met_file_id, l1b_file_id, &
       met_P_levels, met_T_profiles, met_SH_profiles, met_psurf)

    integer(hid_t), intent(in) :: met_file_id
    integer(hid_t), intent(in) :: l1b_file_id
    double precision, allocatable, intent(inout) :: met_T_profiles(:,:,:)
    double precision, allocatable, intent(inout) :: met_P_levels(:,:,:)
    double precision, allocatable, intent(inout) :: met_SH_profiles(:,:,:)
    double precision, allocatable, intent(inout) :: met_psurf(:,:)

    logical :: MET_exists, ECMWF_exists
    character(len=*), parameter :: fname = "read_MET_and_some_L1B_data"
    character(len=999) :: dset_name
    integer(hsize_t), allocatable :: dset_dims(:)
    integer :: hdferr

    ! MET data read-in is not required for space-solar observation modes
    if (MCS%algorithm%observation_mode == "downlooking") then
       ! Get the necesary MET data profiles. OCO-like MET data can have either
       ! /Meteorology or /ECMWF (at least at the time of writing this). So we check
       ! which one exists (priority given to /Meteorology) and take it from there

       MET_exists = .false.
       ECMWF_exists = .false.

       call h5lexists_f(met_file_id, "/Meteorology", MET_exists, hdferr)
       if (.not. MET_exists) then
          call h5lexists_f(met_file_id, "/ECMWF", ECMWF_exists, hdferr)
       end if

       ! Let the user know which one we picked.
       if (MET_exists) then
          call logger%info(fname, "Taking MET data from /Meteorology")
       else if (ECMWF_exists) then
          call logger%info(fname, "Taking MET data from /ECMWF")
       else if ((.not. MET_exists) .and. (.not. ECMWF_exists)) then
          ! Uh-oh, neither /Meteorology nor /ECMWF are found in the
          ! MET file. Can't really go on without MET data.
          call logger%fatal(fname, "Neither /Meteorology nor /ECMWF exist in MET file.")
          stop 1
       end if

       ! Read the complete MET arrays from the corresponding HDF5 fields
       if (MET_exists) dset_name = "/Meteorology/vector_pressure_levels_met"
       if (ECMWF_exists) dset_name = "/ECMWF/vector_pressure_levels_ecmwf"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_P_levels, dset_dims)
       call logger%trivia(fname, "Finished reading in pressure levels.")

       if (MET_exists) dset_name = "/Meteorology/temperature_profile_met"
       if (ECMWF_exists) dset_name = "/ECMWF/temperature_profile_ecmwf"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_T_profiles, dset_dims)
       call logger%trivia(fname, "Finished reading in temperature profiles.")

       if (MET_exists) dset_name = "/Meteorology/specific_humidity_profile_met"
       if (ECMWF_exists) dset_name = "/ECMWF/specific_humidity_profile_ecmwf"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_SH_profiles, dset_dims)
       call logger%trivia(fname, "Finished reading in specific humidity profiles.")

       if (MET_exists) dset_name = "/Meteorology/surface_pressure_met"
       if (ECMWF_exists) dset_name = "/ECMWF/surface_pressure_ecmwf"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_psurf, dset_dims)
       call logger%trivia(fname, "Finished reading in surface pressure.")

    end if

  end subroutine read_MET_data

end module oco2_mod
