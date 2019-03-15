!> @brief Physics-based retrieval method
!> @file physical_model.f90
!> @author Peter Somkuti
!>
!! These are the subroutines required to do a physics-based retrieval. Just like the
!! guanter_model_mod, this module sets up some module-wide variables, and then loops
!! over all specified windows to perform the retrievals for all soundings according to
!! the retrieval options from the config file. All these module-wide variables are
!! then being accessed buy the physical_fm (physical forward model) subroutine. Obviously
!! this requires all the L1b and MET data (apart from the spectra) to be read into
!! memory, so the memory footprint is going to be a few GB. However, this makes it also
!! fairly fast.

module physical_model_mod

  ! User modules
  use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
       read_DP_hdf_dataset, write_INT_hdf_dataset
  use control_mod, only: MCS, MAX_WINDOWS, MAX_GASES, MCS_find_gases
  use instruments_mod, only: generic_instrument
  use oco2_mod
  use solar_model_mod
  use math_utils_mod
  use statevector_mod
  use radiance_mod
  use absco_mod
  use gas_tau_mod
  use spectroscopy_utils_mod

  ! Third-party modules
  use logger_mod, only: logger => master_logger
  use mod_datetime

  ! System modules
  use ISO_FORTRAN_ENV
  use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  implicit none

  public :: physical_retrieval
  private :: physical_fm, calculate_dispersion_limits

  !> A simple structure to keep the atmosphere data nice and tidy
  type atmosphere
     !> The name(s) of the gas(es), (gas number)
     type(string), allocatable :: gas_names(:)
     !> Gas mixing ratios (level, gas number)
     double precision, allocatable :: gas_vmr(:,:)
     !> To which spectroscopy data does this gas correspond to? (gas number)
     integer, allocatable :: gas_index(:)
     !> Model atmosphere temperature
     double precision, allocatable :: T(:)
     !> Model atmosphere pressure
     double precision, allocatable :: p(:)
     !> Model atmosphere specific humidity
     double precision, allocatable :: sh(:)
  end type atmosphere

  !> This structure contains the result data that will be stored in the output HDF file.
  type result_container
     !> State vector names (SV number)
     type(string), allocatable :: sv_names(:)
     !> Retrieved state vector (SV number, footprint, frame)
     double precision, allocatable :: sv_retrieved(:,:,:)
     !> State vector prior (SV number, footprint, frame)
     double precision, allocatable :: sv_prior(:,:,:)
     !> State vector posterior uncertainty (SV number, footprint, frame)
     double precision, allocatable :: sv_uncertainty(:,:,:)
     !> Column-average dry air mixing ratio for retrieved gases (gas number, footprint, frame)
     double precision, allocatable :: xgas(:,:,:)
     !> Final Chi2 (footprint, frame)
     double precision, allocatable :: chi2(:,:)
     !> Final residual RMS (footprint, frame)
     double precision, allocatable :: residual_rms(:,:)
     !> Final dsigma-squared
     double precision, allocatable :: dsigma_sq(:,:)
     !> Final number of iterations
     integer, allocatable :: num_iterations(:,:)
     !> Converged or not? (1=converged, 0=not converged, -1=not properly run)
     integer, allocatable :: converged(:,:)
     !> SNR estimate (mean of per-pixel SNR)
     double precision, allocatable :: SNR(:,:)
     !> SNR standard deviation (std of per-pixel SNR)
     double precision, allocatable :: SNR_std(:,:)
     !> Continuum level radiance estimate
     double precision, allocatable :: continuum(:,:)
     !> Number of moles of dry air per m2 for various sections
     !> of the model atmosphere - corresponding to retrieved
     !> gas scale factors.
     double precision, allocatable :: ndry(:,:,:)
  end type result_container


  !> Spacing of the high-resolution wavelength grid in um
  double precision :: hires_spacing
  !> Padding (left and right) of the hires wl grid to accomodate ILS
  double precision :: hires_pad
  !> High resolution wavelength grid in um
  double precision, allocatable :: hires_grid(:)
  !> Number of gridpoints on the high resolution grid
  integer :: N_hires
  !> Dispersion/wavelength array (detector pixel, footprint, band)
  double precision, allocatable :: dispersion(:,:,:)
  !> Sounding IDs (integer type 8 to accomodate OCO-2 for now), (footprint, frame)
  integer(8), allocatable :: sounding_ids(:,:)
  !> Sounding time strings
  character(len=25), allocatable :: sounding_time_strings(:,:)
  !> Solar zenith angles
  double precision, allocatable :: SZA(:,:)
  !> Solar azimuth angles
  double precision, allocatable :: SAA(:,:)
  !> Viewing/satellite zenith angles
  double precision, allocatable :: VZA(:,:)
  !> Viewing/satellite azimuth angles
  double precision, allocatable :: VAA(:,:)
  !> Sounding footprint longitudes
  double precision, allocatable :: lon(:,:)
  !> Sounding footprint latitudes
  double precision, allocatable :: lat(:,:)
  !> Sounding footprint altitudes
  double precision, allocatable :: altitude(:,:)
  !> Relative velocities between footprint on surface and spacecraft
  double precision, allocatable :: relative_velocity(:,:)
  !> Relative velocities between footprint on surface and the fixed sun
  double precision, allocatable :: relative_solar_velocity(:,:)

  !> ILS spectral dimension (delta wavelength, um), shape: (wl, pixel, footprint, band)
  double precision, allocatable :: ils_delta_lambda(:,:,:,:)
  !> ILS amplitude, shape: (wl, pixel, footprint, band)
  double precision, allocatable :: ils_relative_response(:,:,:,:)

  !> L1B dispersion coefficients (pixel, footprint, band)
  double precision, allocatable :: dispersion_coefs(:,:,:)
  !> L1B SNR coefficients needed for noise calculation (coefficient, pixel, footprint, band)
  double precision, allocatable :: snr_coefs(:,:,:,:)
  !> List of 'bad'-flagged detector pixels
  integer, allocatable :: bad_sample_list(:,:,:)
  !> If required, an array to hold the spike value (pixel, footprint, band)
  integer, allocatable :: spike_list(:,:,:)

  !> MET data temperature profiles (level, footprint, frame)
  double precision, allocatable, dimension(:,:,:) :: met_T_profiles
  !> MET data pressure levels (level, footprint, frame)
  double precision, allocatable, dimension(:,:,:) :: met_P_levels
  !> MET data specific humidity profiles (level, footprint, frame)
  double precision, allocatable, dimension(:,:,:) :: met_SH_profiles
  !> MET data surface pressure (footprint, frame)
  double precision, allocatable :: met_psurf(:,:)

  ! The initial atmosphere provided by the config file, this is the basis
  ! for all scene-dependent atmospheres. The pressure levels of initial_atm
  ! will be used for all retrievals, but the T and SH profiles are taken from
  ! the MET files (if they exist).

  !> Initial atmosphere as read from the text file
  type(atmosphere) :: initial_atm

  !> The solar (pseudo-transmittance) spectrum as read in from file (wavelength, transmission)
  double precision, allocatable :: solar_spectrum(:,:)
  !> The solar (pseudo-transmittance) spectrum on the regular wavelength grid (wavelength, transmission)
  double precision, allocatable :: solar_spectrum_regular(:,:)
  !> The number of solar spectrum points (as read from the file)
  integer :: N_solar

  !> Final modelled radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: final_radiance
  !> Measured radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: measured_radiance
  !> Calculated noise radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: noise_radiance

  !> State vector construct
  type(statevector) :: SV
  !> Result container
  type(result_container) :: results


contains

  !> @brief The physical_retrieval routine reads in all the necessary L1b and MET
  !> data and prepares for the forward model to be run.
  !> @param my_instrument Instrument entity
  !>
  !> TODO Detailed description of the workings
  subroutine physical_retrieval(my_instrument)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument

    ! HDF5 file handlers for the L1b file and the MET file
    integer(hid_t) :: l1b_file_id, met_file_id, output_file_id
    ! Variable to hold group ID
    integer(hid_t) :: result_gid
    ! Do MET or ECMWF groups exist in the MET file?
    logical :: MET_exists, ECMWF_exists
    ! Does a spike/bad_sample data field exist?
    logical :: spike_exists, bad_sample_exists
    ! Fixed dimensions for the arrays to be saved into the output HDF file
    integer(hsize_t) :: out_dims2d(2), out_dims3d(3)
    ! Dimensions for the read-in of MET/L1b data
    integer(hsize_t), allocatable :: dset_dims(:)
    ! HDF error variable
    integer :: hdferr
    ! Holders for temporary strings and dataset names to grab L1b/MET data
    character(len=999) :: dset_name, tmp_str, group_name
    ! Name of this function for debugging
    character(len=*), parameter :: fname = "physical_retrieval"
    ! What are the dimensions of the ABSCO file used?
    integer :: absco_dims
    ! What are the smallest and largest values of all ILSs in the L1b file
    ! for this given band? Needed to calculate hires_pad
    double precision :: ils_min_wl, ils_max_wl
    ! Number of frames, footprints, detector pixels and bands
    ! required vor various loops
    integer :: num_frames, num_fp, num_pixel, num_band
    ! Loop variables
    integer :: i, j, i_fp, i_fr, i_win, band
    ! Retrieval count
    integer :: retr_count
    ! Return value of physical_fm, tells us whether this one has converged or not
    logical :: this_converged
    ! File unit for debugging only..
    integer :: funit
    ! CPU time stamps and mean duration for performance analysis
    double precision :: cpu_time_start, cpu_time_stop, mean_duration

    ! Open up the MET file
    call h5fopen_f(MCS%input%met_filename%chars(), &
         H5F_ACC_RDONLY_F, MCS%input%met_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening MET file: " &
         // trim(MCS%input%met_filename%chars()))

    ! Store HDF file handler for more convenient access
    met_file_id = MCS%input%met_file_id
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id

    ! Grab number of frames and footprints
    num_frames = MCS%general%N_frame
    num_fp = MCS%general%N_fp
    ! Grab number of bands
    num_band = MCS%general%N_bands

    ! Read in MET and L1B data, this is instrument-dependent, so we need to
    ! make our instrument select here as well.
    select type(my_instrument)
    type is (oco2_instrument)

       ! Get the necesary MET data profiles. OCO-like MET data can have either
       ! /Meteorology or /ECMWF (at least at the time of writing this). So we check
       ! which one exists (priority given to /Meteorology) and take it from there

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

       ! Grab the SNR coefficients for noise calculations
       call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
       ! Conveniently grab the number of spectral pixels and number of bands


       ! Read dispersion coefficients and create dispersion array
       call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
       allocate(dispersion(num_pixel, num_fp, num_band))

       do band=1, num_band
          do i_fp=1, num_fp
             call my_instrument%calculate_dispersion(dispersion_coefs(:, i_fp, band), &
                  dispersion(:, i_fp, band), band, i_fp)
          end do
       end do

       ! Read in the sounding id's
       call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)
       ! Read the time strings
       call my_instrument%read_time_strings(l1b_file_id, sounding_time_strings)
       ! Read in the instrument ILS data
       call my_instrument%read_ils_data(l1b_file_id, ils_delta_lambda, &
            ils_relative_response)

       ! Read in Spike filter data, if it exists in this file
       call h5lexists_f(l1b_file_id, "/SpikeEOF", spike_exists, hdferr)
       !if (spike_exists) then
       !   call my_instrument%read_spike_filter(l1b_file_id, spike_list, band)
       !end if

       ! Read in bad sample list (if it exists)
       call h5lexists_f(met_file_id, "/InstrumentHeader/bad_sample_list", bad_sample_exists, hdferr)
       if (bad_sample_exists) then
          call my_instrument%read_bad_sample_list(l1b_file_id, bad_sample_list)
       end if

    end select



    ! Create the HDF group in which all the results go in the end
    call h5gcreate_f(output_file_id, "physical_retrieval_results", result_gid, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not create group: physical_retrieval_results")

    ! Loop over all potential user-defined retrieval windows.
    do i_win=1, MAX_WINDOWS

       ! Just skip unused windows
       if (.not. MCS%window(i_win)%used) cycle

       ! At the beginning, we check which gases were defined for this window,
       ! and see if a gas with the corresponding name has been defined.
       call MCS_find_gases(MCS%window, MCS%gas, i_win)

       ! If we have gases, we want to read in the corresponding spectroscopy data
       if (MCS%window(i_win)%num_gases > 0) then
          ! Read in the spectroscopy data, depending on the type
          do i=1, size(MCS%window(i_win)%gases)
             ! Which gas are we reading in? As in 'index of MCS%gases'
             j = MCS%window(i_win)%gas_index(i)

             if (MCS%gas(j)%type%lower() == "absco") then
                call logger%trivia(fname, "Reading in ABSCO-type gas: " // MCS%window(i_win)%gases(i))
                call read_absco_HDF(MCS%gas(j)%filename%chars(), MCS%gas(j), absco_dims, &
                     MCS%window(i_win)%wl_min + hires_pad, & ! This is not yet used
                     MCS%window(i_win)%wl_max + hires_pad) ! Also this
             else
                call logger%fatal(fname, "Spectroscopy type: " // MCS%gas(j)%type &
                     // " not implemented!")
                stop 1
             end if
          end do

          ! Read the atmosphere file (if present) and populate the initial_atm structure
          ! with the contents of said file. It is cross-referenced against the gases in the
          ! list of gases in this window. This is non-optional obviously, so the program is
          ! going to stop if it cannot find the specified file.

          call logger%debug(fname, "Looking into atmosphere file.")
          call read_atmosphere_file(&
               MCS%window(i_win)%atmosphere_file%chars(), &
               MCS%window(i_win)%gases, &
               initial_atm)

       end if

       ! The currently used band / spectrometer number
       band = MCS%window(i_win)%band

       ! We are reading in the sounding geometry and location on a
       ! per-band basis. This might not be necessary for all instruments,
       ! but why not make use of per-band data (like for OCO-2).
       select type(my_instrument)
       type is (oco2_instrument)
          ! Read in the measurement geometries - per band
          call my_instrument%read_sounding_geometry(l1b_file_id, band, SZA, SAA, VZA, VAA)
          ! Read in the measurement location
          call my_instrument%read_sounding_location(l1b_file_id, band, lon, lat, &
               altitude, relative_velocity, relative_solar_velocity)
       end select

       ! Grab the number of spectral pixels in this band
       num_pixel = MCS%general%N_spec(band)

       ! Find the smallest and largest delta-lambda values for all ILSs
       ! in this given band to calculate hires_pad.
       ils_min_wl = minval(ils_delta_lambda(1,:,:,band))
       ils_max_wl = maxval(ils_delta_lambda(size(ils_delta_lambda, 1),:,:,band))

       ! This is the amount by which we have to pad the hi-resolution grid in order to
       ! allow for ILS protrusion. (and add small percentage to be on the safe side)
       hires_pad = (ils_max_wl - ils_min_wl) * 1.10d0

       ! Grab the desired high-resolution wavelength grid spacing
       hires_spacing = MCS%window(i_win)%wl_spacing

       ! .. and construct the high-resolution grid from the supplied
       ! microwindow range.
       N_hires = ceiling((MCS%window(i_win)%wl_max - MCS%window(i_win)%wl_min + 2*hires_pad) &
            / hires_spacing)

       write(tmp_str, '(A,G0.1)') "Number of hires spectral points: ", N_hires
       call logger%info(fname, trim(tmp_str))

       allocate(hires_grid(N_hires))
       do i=1, N_hires
          hires_grid(i) = MCS%window(i_win)%wl_min - hires_pad + dble(i-1) * hires_spacing
       end do

       ! For a faster gas-OD calculation, we re-grid the spectroscopy data
       ! as well, such that we do not have to interpolate in the wavelength
       ! dimension every single time.

       if (MCS%window(i_win)%num_gases > 0) then
          ! Read in the spectroscopy data, depending on the type
          do i=1, size(MCS%window(i_win)%gases)

             j = MCS%window(i_win)%gas_index(i)
             call regrid_spectroscopy(MCS%gas(j), hires_grid)

          end do
       end if


       ! Read in the solar model - we do this for every band, the reason being
       ! the following. Since the re-gridding procedure is fairly costly, we want
       ! to keep the solar spectrum data as small as possible.

       if (MCS%algorithm%solar_type == "toon") then
          call read_toon_spectrum(MCS%algorithm%solar_file%chars(), &
               solar_spectrum, &
               MCS%window(i_win)%wl_min - hires_pad, &
               MCS%window(i_win)%wl_max + hires_pad)

          N_solar = size(solar_spectrum, 1)
       else
          call logger%fatal(fname, "Sorry, solar model type " &
               // MCS%algorithm%solar_type%chars() &
               // " is not known.")
       end if

       ! And we also need the solar spectrum on our user-defined
       ! high-resolution wavelength grid.
       allocate(solar_spectrum_regular(N_hires, 2))
       solar_spectrum_regular(:,1) = hires_grid

       call logger%debug(fname, "Re-gridding solar spectrum")
       call pwl_value_1d( &
            N_solar, &
            solar_spectrum(:,1), solar_spectrum(:,2), &
            N_hires, &
            solar_spectrum_regular(:,1), solar_spectrum_regular(:,2))
       call logger%debug(fname, "Finished re-gridding solar spectrum.")
       ! Note that at this point, the solar spectrum is still normalised


       ! Allocate containers to hold the radiances and noise values - if requested!
       if (MCS%output%save_radiances) then
          allocate(final_radiance(num_pixel, num_fp, num_frames))
          allocate(measured_radiance(num_pixel, num_fp, num_frames))
          allocate(noise_radiance(num_pixel, num_fp, num_frames))

          ! We fill these with NaNs. This makes plotting a bit easier, since
          ! most plotting routines (at least Matplotlib does it) just skip NaNs,
          ! so you will only see the values that are populated / used.
          final_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          measured_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          noise_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       end if

       ! Set up state vector structure here

       ! We can do this once for every window, and simply clear the contents
       ! of the state vectors inside the loop later on. Might not save an awful
       ! lot, but every little helps, I guess.

       ! At the moment, we can do following state vectors:
       !
       ! Lambertian surface albedo, arbitrary (?) polynomial order
       ! Zero-level offset / SIF (constant)
       ! Instrument disperison
       ! Surface pressure
       ! Solar shift and stretch
       ! Gas scalar factor defined on level ranges

       ! Parsing the statevector string, that was passed in the window
       ! section and initialize the state vector SV accordingly. This subroutine
       ! needs to access plenty of things in the MCS, so we only pass the
       ! window index, and the routine takes care of arranging the rest.

       ! For the beginning, we start by initialising it with the number of levels
       ! as obtained from the initial_atm

       call parse_and_initialize_SV(i_win, size(initial_atm%p), SV)
       call logger%info(fname, "Initialised SV structure")

       ! And now set up the result container with the appropriate sizes for arrays
       call create_result_container(results, num_frames, num_fp, &
            size(SV%svap), SV%num_gas)

       ! Create the SV names corresponding to the SV indices
       call assign_SV_names_to_result(results, SV, i_win)

       ! retr_count keeps track of the number of retrievals processed
       ! so far, and the mean_duration keeps track of the average
       ! processing time.
       retr_count = 0
       mean_duration = 0.0d0

       call logger%info(fname, "Starting main retrieval loop!")

       do i_fr=1, num_frames, MCS%window(i_win)%frame_skip
          do i_fp=1, num_fp, MCS%window(i_win)%footprint_skip

             call cpu_time(cpu_time_start)
             ! Do the retrieval for this particular sounding
             this_converged = physical_FM(my_instrument, i_fp, i_fr, i_win, band)
             ! After the very first call, we set first_band_call to false
             call cpu_time(cpu_time_stop)

             ! Increase the rerival count tracker and compute the average processing
             ! time per retrieval.
             retr_count = retr_count + 1
             mean_duration = mean_duration * (retr_count)/(retr_count+1) + &
                  (cpu_time_stop - cpu_time_start) / (retr_count+1)

             if (mod(retr_count, 1) == 0) then
                write(tmp_str, '(A, G0.1, A, G0.1, A, F5.2, A, F10.5, A, L1)') &
                     "Frame/FP: ", i_fr, "/", i_fp, " ( ", &
                     dble(100 * dble(retr_count) / dble(num_fp * num_frames)), "%) - ", &
                     (cpu_time_stop - cpu_time_start), ' sec. - Converged: ', this_converged
                call logger%debug(fname, trim(tmp_str))
             end if

          end do
       end do

       ! Create an HDF group for all windows separately
       group_name = "/physical_retrieval_results/" // trim(MCS%window(i_win)%name%chars())
       call h5gcreate_f(output_file_id, trim(group_name), result_gid, hdferr)
       call check_hdf_error(hdferr, fname, "Error. Could not create group: " &
            // trim(group_name))

       ! Set the dimensions of the arrays for saving them into the HDF file
       out_dims2d(1) = num_fp
       out_dims2d(2) = num_frames

       ! Writing out the prior surface pressure
       call logger%info(fname, "Writing out: " // trim(group_name) // "/prior_psurf")
       write(tmp_str, '(A,A,A)') trim(group_name) // "/prior_psurf"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), met_psurf(:,:), out_dims2d, -9999.99d0)

       ! Save the retrieved state vectors
       do i=1, size(SV%svsv)
          write(tmp_str, '(A,A,A,A)') trim(group_name) // "/" // results%sv_names(i)
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_retrieved(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the retrieved state vector uncertainties
       do i=1, size(SV%svsv)
          write(tmp_str, '(A,A,A,A,A)') trim(group_name) // "/" &
               // results%sv_names(i) // "_uncertainty"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_uncertainty(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the dry air mass (moles / cm2)
       do i=1, SV%num_gas
          if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then

             ! This is a scale-type jacobian
             do j=1, size(MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(i),:))

                if (SV%gas_idx_lookup(i) == j) then
                   write(tmp_str, "(A, A, A)") trim(group_name) // "/", &
                        results%SV_names(SV%idx_gas(i,1))%chars(), "_ndry"
                   call logger%info(fname, "Writing out: " // trim(tmp_str))
                   call write_DP_hdf_dataset(output_file_id, &
                        trim(tmp_str), &
                        results%ndry(:,:,i), out_dims2d, -9999.99d0)
                end if
             end do
          end if
       end do


       do i=1,MCS%window(i_win)%num_gases
          if (MCS%window(i_win)%gas_retrieved(i)) then
             ! Save XGAS for each gas
             write(tmp_str, '(A,A,A,A)') trim(group_name) // "/X" // MCS%window(i_win)%gases(i)
             call logger%info(fname, "Writing out: " // trim(tmp_str))
             call write_DP_hdf_dataset(output_file_id, &
                  trim(tmp_str), &
                  results%xgas(:,:,i), out_dims2d, -9999.99d0)
          end if
       end do

       ! Save number of iterations
       call logger%info(fname, "Writing out: " // trim(group_name) // "/num_iterations")
       write(tmp_str, '(A,A,A)') trim(group_name) // "/num_iterations"
       call write_INT_hdf_dataset(output_file_id, &
            trim(tmp_str), results%num_iterations(:,:), out_dims2d, -9999)

       ! Save converged status
       call logger%info(fname, "Writing out: " // trim(group_name) // "/converged")
       write(tmp_str, '(A,A,A)') trim(group_name) // "/converged"
       call write_INT_hdf_dataset(output_file_id, &
            trim(tmp_str), results%converged(:,:), out_dims2d, -9999)

       ! Retrieved CHI2
       call logger%info(fname, "Writing out: " // trim(group_name) // "/retrieved_chi2")
       write(tmp_str, '(A,A,A)') trim(group_name) // "/retrieved_chi2"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%chi2(:,:), out_dims2d, -9999.99d0)

       ! Dsigma_sq
       call logger%info(fname, "Writing out: " // trim(group_name) // "/final_dsigma_sq")
       write(tmp_str, '(A,A,A)') trim(group_name) // "/final_dsigma_sq"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%dsigma_sq(:,:), out_dims2d, -9999.99d0)

       call logger%info(fname, "Writing out: " // trim(group_name) // "/SNR")
       write(tmp_str, '(A,A)') trim(group_name) // "/SNR"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%SNR, out_dims2d)

       call logger%info(fname, "Writing out: " // trim(group_name) // "/SNR_std")
       write(tmp_str, '(A,A)') trim(group_name) // "/SNR_std"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%SNR_std, out_dims2d)

       call logger%info(fname, "Writing out: " // trim(group_name) // "/continuum")
       write(tmp_str, '(A,A)') trim(group_name) // "/continuum"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%continuum, out_dims2d)

       ! Save the radiances (this should be made optional!)
       if (MCS%output%save_radiances) then
          out_dims3d = shape(final_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/modelled_radiance")
          write(tmp_str, '(A,A)') trim(group_name) // "/modelled_radiance"
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               final_radiance, out_dims3d)

          out_dims3d = shape(measured_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/measured_radiance")
          write(tmp_str, '(A,A)') trim(group_name) // "/measured_radiance"
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               measured_radiance, out_dims3d)

          out_dims3d = shape(noise_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/noise_radiance")
          write(tmp_str, '(A,A)') trim(group_name) // "/noise_radiance"
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               noise_radiance, out_dims3d)

          ! Also deallocate containers holding the radiances
          deallocate(final_radiance)
          deallocate(measured_radiance)
          deallocate(noise_radiance)
       end if

       ! Clean-up phase! We need to deallocate / destroy a few arrays here, so
       ! they can be freshly reallocated for the next retrieval microwindow.

       ! Deallocate arrays that are allocated on per-window basis
       deallocate(solar_spectrum_regular, solar_spectrum, hires_grid)

       ! Clear and deallocate the SV structure to be ready for the next window.
       ! We really shouldn't need to check if this is deallocated already or not,
       ! since parse_and_initialize_SV should have successfully allocated
       ! all fields in SV.
       call clear_SV(SV)
       call logger%info(fname, "Clearing up SV structure")

       deallocate(SZA, SAA, VZA, VAA)
       deallocate(lon, lat, altitude, relative_velocity, relative_solar_velocity)

       if (allocated(bad_sample_list)) deallocate(bad_sample_list)
       if (allocated(spike_list)) deallocate(spike_list)

       ! Clear and deallocate the result container
       call logger%info(fname, "Clearing up results container.")
       call destroy_result_container(results)

       ! End of loop over windows
    end do

  end subroutine physical_retrieval

  !> @brief Physical-type forward model / retrieval
  !> @param my_instrument Instrument instance
  !> @param i_fp Footprint index
  !> @param i_fr Frame index
  !> @param i_win Window index for MCS
  !> @param band Band number
  !> @param converged Whether this retrieval has converged or not
  !>
  !> This function performs the full physical retrieval, and returns whether
  !> it converged or not. It accesses all the L1B/MET arrays defined in the module
  !> for fast readout and processing. The OE scheme is based on Rodgers (2000),
  !> and so far we are doing the LM-modification to the Gauss-Newton scheme.
  function physical_FM(my_instrument, i_fp, i_fr, i_win, band) result(converged)

    implicit none

    class(generic_instrument), intent(in) :: my_instrument
    integer, intent(in) :: i_fr, i_fp, i_win, band
    logical :: converged

    ! HDF file id handlers for L1B and output file
    integer(hid_t) :: l1b_file_id, output_file_id

    ! Radiances and noise arrays. We need a few of these to hold
    ! the TOA radiances for high and low-res calculations. And then
    ! some 'temp' arrays to hold radiances from e.g. perturbation or
    ! linear prediction.
    double precision, allocatable :: radiance_l1b(:), &
         radiance_tmp_work(:), &
         radiance_meas_work(:), &
         radiance_calc_work(:), &
         radiance_tmp_work_hi(:), &
         radiance_calc_work_hi(:), &
         noise_work(:)

    ! The current atmosphere. This will be a copy from the initial_atm
    ! constructed in the physical_retrieval subroutine, but T and SH
    ! profiles are taken on a scene-by-scene basis.
    type(atmosphere) :: this_atm

    ! Dispersion/wavelength indices - these tell us, where in the L1B radiance array
    ! we can extract the spectra to match the user-defined wavelength ranges.
    integer :: l1b_wl_idx_min, l1b_wl_idx_max

    ! Sounding time stuff
    type(datetime) :: date ! Datetime object for sounding date/time
    double precision :: doy_dp ! Day of year as double precision
    ! Sounding location variables
    double precision :: mu0, mu ! cos(sza) and cos(vza)

    ! Instrument stuff
    ! Instrument doppler shift based on relative motion between ground footprint
    ! and spacecraft velocity.
    double precision :: instrument_doppler
    ! Dispersion arrays to hold the wavelength-per-pixel values, as well as the
    ! arrays to hold the polynomial coefficents to create them. Also needs "tmp"
    ! and "pert" versions for perturbed values if dispersion is retrieved.
    double precision, allocatable :: this_dispersion(:), &
         this_dispersion_tmp(:), &
         this_dispersion_coefs(:), &
         this_dispersion_coefs_pert(:)

    ! For some SV elements, Jacobians are calculated separately from the
    ! main Jacobian convolution loop. This varaible indicates whether this
    ! loop iteration should be skipped.
    logical :: skip_jacobian

    ! Number of spectral points for low-res and hires grid,
    ! number of state vector elements.
    integer :: N_spec, N_spec_hi, N_sv

    ! Solar stuff
    ! Solar distance, solar relative velocity, earth relative velocity, and
    ! solar doppler shift - which happens because of the relative motion between
    ! the moving and rotating Earth and the fixed sun.
    double precision :: solar_dist, solar_rv, earth_rv, solar_doppler
    ! This solar is the full solar spectrum, i.e. the pseudo-transmittance multiplied
    ! by the solar continuum / irradiance
    double precision, dimension(:,:), allocatable :: this_solar
    ! Solar irradiance / continuum, derivative of solar spectrum with respect to wavelength
    ! (needed for the solar shift and stretch Jacobians)
    double precision, allocatable :: solar_irrad(:), dsolar_dlambda(:), solar_tmp(:)
    ! Per-iteration values for the solar shift value and the solar stretch factor
    double precision :: this_solar_shift, this_solar_stretch

    ! Atmosphere stuff
    ! Number of gases, total levels, and number of active levels (changes with surface pressure)
    integer :: num_gases, num_levels, num_active_levels
    ! Various arrays for gas optical depths
    double precision, allocatable :: gas_tau(:,:,:) ! Gas optical depth (spectral, layer, gas number)
    double precision, allocatable :: gas_tau_dvmr(:,:,:) ! dtau / dvmr (spectral, layer, gas number)
    double precision, allocatable :: gas_tau_dpsurf(:,:,:) ! dtau / dpsurf (spectral, layer, gas_number)
    double precision, allocatable :: gas_tau_pert(:,:,:,:) ! Perturbed gas optical depth (spectral, layer, gas number)

    ! Perturbed VMR profile and per-iteration-and-per-gas VMR profile for OD calculation (level)
    double precision, allocatable :: vmr_pert(:), this_vmr_profile(:)
    ! Rayleigh extinction optical depth (spectral, layer)
    double precision, allocatable :: ray_tau(:,:)
    ! Total column optical depth (spectral)
    double precision, allocatable :: total_tau(:)
    ! Start and end positions in the atmosphere of the gas scalar
    integer :: s_start(SV%num_gas), s_stop(SV%num_gas)
    ! Is this gas H2O?
    logical :: is_H2O
    ! Number of molecules of dry air per m^2
    double precision, allocatable :: ndry(:), ndry_tmp(:)

    ! Albedo
    ! Prior albedo value estimated from the radiances
    double precision :: albedo_apriori
    ! Per-wavelength albedo for hires and low-res spectra
    double precision, allocatable :: albedo(:), albedo_low(:)

    !! Surface pressure
    ! The surface pressure per-iteration
    double precision :: this_psurf
    ! Do we need/want to calculate surface pressure Jacobians?
    logical :: do_psurf_jac

    ! SIF
    ! Per-iteration SIF radiance
    double precision :: this_sif_radiance
    ! ZLO - same as SIF essentially
    double precision :: this_zlo_radiance

    ! Gases
    ! Do we want to calculate gas Jacobians
    logical :: do_gas_jac
    ! Was the calculation of gas ODs successful?
    logical :: success_gas

    ! Retrieval quantities
    ! User-defined value to scale the dsigma_sq value that essentially
    ! controls convergence. A higher value will lead to faster convergence.
    double precision :: dsigma_scale
    ! disgma_sq: see Rodgers Eq. 5.29 (change in SV as fraction of Shat_inv)
    double precision :: dsigma_sq
    ! Jacobian matrices for high-res and lowres spectral grids
    ! they have the shape (N_spec_hi, N_sv) and (N_spec, N_sv)
    double precision, allocatable :: K_hi(:,:), K(:,:)
    ! Noise covariance, prior covariance (and its inverse), shape: (N_sv, N_sv)
    double precision, allocatable :: Se_inv(:,:), Sa(:,:), Sa_inv(:,:)
    ! Posterior covariance matrix and its inverse (N_sv, N_sv)
    double precision, allocatable :: Shat(:,:), Shat_inv(:,:)

    ! Temporary matrices and vectors for computation
    double precision, allocatable :: tmp_m1(:,:), tmp_m2(:,:)
    double precision, allocatable :: tmp_v1(:), tmp_v2(:)
    double precision, allocatable :: KtSeK(:,:), gain_matrix(:,:), AK(:,:)
    ! Was the matrix inversion operation successful?
    logical :: success_inv_mat
    ! Was the convolution operation successful?
    logical :: ILS_success

    ! The SV from last iteration, as well as last successful (non-divergent) iteration
    ! shape (N_sv)
    double precision, allocatable :: old_sv(:), last_successful_sv(:)
    ! Pressure weighting function (N_sv)
    double precision, allocatable :: pwgts(:)
    ! Levenberg-Marquart gamma, per-iteration
    double precision :: lm_gamma
    ! Chi2-related variables. Chi2 of last and this current iteration,
    ! chi2 calculated from a linear prediction, and the chi2 ratio needed
    ! to adjust lm_gamma and determine a divergent step.
    double precision :: old_chi2, this_chi2, linear_prediction_chi2, chi2_ratio

    ! Iteration-related
    ! Current iteration number (starts at 1), number of divergent steps allowed.
    integer :: iteration, num_divergent_steps
    ! Variable to determine if we keep doing the iteration loop
    logical :: keep_iterating
    ! Has last iteration been a divergent one?
    logical :: divergent_step

    ! Miscellaneous stuff
    ! String to hold various names etc.
    character(len=999) :: tmp_str
    ! Function name
    character(len=*), parameter :: fname = "physical_FM"
    ! Loop variables
    integer :: i, j, l
    ! File unit for debugging
    integer :: funit


    ! Take a local copy of the HDF file ID handlers
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id

    ! What is the total number of gases in this window,
    ! regardless of whether they are retrieved or not.
    if (allocated(MCS%window(i_win)%gases)) then
       num_gases = MCS%window(i_win)%num_gases
    else
       num_gases = 0
    end if

    ! Use user-supplied value for convergence critertion
    dsigma_scale = MCS%window(i_win)%dsigma_scale
    if (dsigma_scale < 0.0d0) dsigma_scale = 1.0d0

    select type(my_instrument)
    type is (oco2_instrument)
       ! Read the L1B spectrum for this one measurement!
       call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, band, &
            MCS%general%N_spec(band), radiance_l1b)
       ! Convert the date-time-string object in the L1B to a date-time-object "date"
       call my_instrument%convert_time_string_to_date(sounding_time_strings(i_fp, i_fr), date)

    end select

    doy_dp = dble(date%yearday()) + dble(date%getHour()) / 24.0d0

    write(tmp_str, "(A, A)") "Date: ", date%isoformat()
    call logger%debug(fname, trim(tmp_str))

    write(tmp_str, "(A, F5.1)") "Day of the year: ", doy_dp
    call logger%debug(fname, trim(tmp_str))

    !write(*,*) solar_rv * 1000.0d0 + earth_rv, solar_doppler
    write(tmp_str, "(A, F6.2, A, F6.2, A, F6.2, A, F6.2)") &
         "SZA: ", SZA(i_fp, i_fr), " - SAA: ", SAA(i_fp, i_fr), &
         " - VZA: ", VZA(i_fp, i_fr), " - VAA: ", VAA(i_fp, i_fr)
    call logger%debug(fname, trim(tmp_str))

    write(tmp_str, "(A, F6.2, A, F6.2, A, F6.2)") &
         "Longitude: ", lon(i_fp, i_fr), &
         " Latitude: ", lat(i_fp, i_fr), &
         " Altitude: ", altitude(i_fp, i_fr)

    ! Dispersion array that contains the wavelenghts per pixel
    allocate(this_dispersion(size(radiance_l1b)))
    allocate(this_dispersion_tmp(size(radiance_l1b)))

    ! The dispersion coefficients use to generate this soundings' dispersion
    ! has the same size as the L1b dispersion coefficient array
    allocate(this_dispersion_coefs(size(dispersion_coefs, 1)))
    allocate(this_dispersion_coefs_pert(size(dispersion_coefs, 1)))

    ! Allocate the micro-window bounded solar arrays
    allocate(this_solar(N_hires, 2))
    allocate(dsolar_dlambda(N_hires))

    ! If solar doppler-shift is needed, calculate the distance and relative
    ! velocity between point of measurement and the (fixed) sun
    call calculate_solar_distance_and_rv(doy_dp, solar_dist, solar_rv)

    call calculate_rel_velocity_earth_sun(lat(i_fp, i_fr), SZA(i_fp, i_fr), &
         SAA(i_fp, i_fr), altitude(i_fp, i_fr), earth_rv)

    solar_doppler = 0.0d0
    instrument_doppler = 0.0d0
    !solar_doppler = (earth_rv + solar_rv * 1000.0d0) / SPEED_OF_LIGHT
    !instrument_doppler = relative_velocity(i_fp, i_fr) / SPEED_OF_LIGHT


    ! Set up retrieval quantities:
    N_sv = size(SV%svap)
    N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1
    N_spec_hi = size(this_solar, 1)
    ! Set the initial LM-Gamma parameter
    lm_gamma = MCS%window(i_win)%lm_gamma

    ! For convenience, calculate cos(sza) and cos(vza) here
    mu0 = cos(DEG2RAD * SZA(i_fp, i_fr))
    mu = cos(DEG2RAD * VZA(i_fp, i_fr))

    ! Output-resolution K is allocated within the loop, as the
    ! number of pixels might change, while the hi-res K stays the same
    allocate(K_hi(N_spec_hi, N_sv))
    K_hi(:,:) = 0.0d0

    ! Allocate OE-related matrices here
    allocate(Sa(N_sv, N_sv))
    allocate(Sa_inv(N_sv, N_sv))
    allocate(Shat_inv(N_sv, N_sv))
    allocate(Shat(N_sv, N_sv))
    allocate(tmp_m1(N_sv, N_sv), tmp_m2(N_sv, N_sv))
    allocate(tmp_v1(N_sv), tmp_v2(N_sv))
    allocate(KtSeK(N_sv, N_sv))
    allocate(AK(N_sv, N_sv))

    ! Allocate Solar continuum (irradiance) array. We do this here already,
    ! since it's handy to have it for estimating the albedo.
    allocate(solar_irrad(N_hires))

    allocate(radiance_calc_work_hi(size(this_solar, 1)))
    allocate(radiance_tmp_work_hi(size(this_solar, 1)))
    allocate(albedo(size(radiance_calc_work_hi)))


    Sa_inv(:,:) = 0.0d0
    Sa(:,:) = 0.0d0
    call populate_prior_covariance(SV, percentile(radiance_l1b, 98.0d0), &
         i_win, Sa, Sa_inv)

    ! Get the albedo prior
    select type(my_instrument)
    type is (oco2_instrument)

       ! OCO-2 has Stokes coefficient 0.5 for intensity, so we need to
       ! take that into account for the incoming solar irradiance
       call calculate_solar_planck_function(6500.0d0, solar_dist * 1000.0d0, &
            solar_spectrum_regular(:,1), solar_irrad)

       albedo_apriori = 1.0d0 * PI * maxval(radiance_l1b) / &
            (1.0d0 * maxval(solar_irrad) * mu0)

    end select


    ! We can now populate the prior state vector

    SV%svap(SV%idx_albedo(1)) = albedo_apriori
    ! Set slope etc. to zero always (why would we ever want to have a prior slope?)
    if (SV%num_albedo > 1) then
       do i=2, SV%num_albedo
          SV%svap(SV%idx_albedo(i)) = 0.0d0
       end do
    end if

    ! Solar shift is set to zero
    if (SV%num_solar_shift == 1) then
       SV%svap(SV%idx_solar_shift(1)) = 0.0d0
    end if

    ! Solar shift factor is set to one
    if (SV%num_solar_stretch == 1) then
       SV%svap(SV%idx_solar_stretch(1)) = 1.0d0
    end if

    ! Start SIF with zero
    if (SV%num_sif > 0) then
       SV%svap(SV%idx_sif(1)) = 0.0d0
    end if

    ! ZLO starts with zero too
    if (SV%num_zlo > 0) then
       SV%svap(SV%idx_zlo(1)) = 0.0d0
    end if

    ! Dispersion
    if (SV%num_dispersion > 0) then
       ! Start with the L1b dispersion values as priors
       do i=1, SV%num_dispersion
          SV%svap(SV%idx_dispersion(i)) = dispersion_coefs(i, i_fp, band)
       end do
    end if

    ! Surface pressure is taken from the MET data
    if (SV%num_psurf == 1) then
       SV%svap(SV%idx_psurf(1)) = met_psurf(i_fp, i_fr)
    end if

    ! Gases. Scale factors are set to one in the beginning.
    if (SV%num_gas > 0) then
       do i=1, SV%num_gas
          if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
             SV%svap(SV%idx_gas(i,1)) = 1.0d0
          end if
       end do
    end if



    ! Regardless of whether they are retrieved or not, the solar shift and stretch
    ! are set to 'default' values here. If we retrieve them, this value is then just
    ! updated from the state vector.
    this_solar_shift = 0.0d0
    this_solar_stretch = 1.0d0

    ! Retrival iteration loop
    iteration = 0
    num_divergent_steps = 0
    keep_iterating = .true.
    divergent_step = .false.
    allocate(old_sv(size(SV%svsv)))
    allocate(last_successful_sv(size(SV%svsv)))

    converged = .false.
    ! Main iteration loop for the retrieval process.
    do while (keep_iterating)

       if (.not. divergent_step) then
          iteration = iteration + 1
       end if

       if (iteration == 1) then

          divergent_step = .false.
          ! Initialise Chi2 with an insanely large value
          this_chi2 = 9.9d9
          ! For the first iteration, we want to use the prior albedo
          albedo(:) = albedo_apriori
          ! and the first guess state vector is the prior
          SV%svsv = SV%svap
          !
          last_successful_sv(:) = SV%svap
          old_sv(:) = SV%svap

          if (num_gases > 0) then

             ! An atmosphere is only required if there are gases present in the
             ! microwindow.
             this_psurf = met_psurf(i_fp, i_fr)
             this_atm = initial_atm
             num_levels = size(this_atm%p)

             ! And get the T and SH MET profiles onto our new atmosphere grid. We are
             ! sampling it on a log(p) grid.

             call pwl_value_1d( &
                  size(met_P_levels, 1), &
                  log(met_P_levels(:,i_fp,i_fr)), met_T_profiles(:,i_fp,i_fr), &
                  size(this_atm%p), &
                  log(this_atm%p), this_atm%T)

             call pwl_value_1d( &
                  size(met_P_levels, 1), &
                  log(met_P_levels(:,i_fp,i_fr)), met_SH_profiles(:,i_fp,i_fr), &
                  size(this_atm%p), &
                  log(this_atm%p), this_atm%sh)

             ! Obtain the number of active levels. Loop from the TOA down to the bottom,
             ! and see where the surface pressure falls into.
             do j=1, num_levels
                if (this_psurf > this_atm%p(j)) then
                   num_active_levels = j+1
                end if

                ! Should SH drop below 0 for whatever reason, shift it back
                ! some tiny value.
                if (this_atm%sh(j) < 0.0d0) this_atm%sh(j) = 1.0d-10
             end do

             ! If psurf > BOA p level, we have a problem and thus can't go on.
             if (num_active_levels > num_levels) then
                write(tmp_str, '(A, F10.3, A, F10.3)') "Psurf at ", this_psurf, &
                     " is larger than p(BOA) at ", this_atm%p(size(this_atm%p))
                call logger%error(fname, trim(tmp_str))
                return
             end if

             do i=1, num_gases
                if (MCS%window(i_win)%gases(i) == "H2O") then
                   ! If H2O is needed to be retrieved, take it from the MET atmosphere
                   ! specific humidty directly, rather than the H2O column of the
                   ! atmosphere text file.
                   ! this_atm%gas_vmr(:,i) = this_atm%sh / (1.0d0 - this_atm%sh) * SH_H2O_CONV
                end if
             end do

          end if

       else
          ! This is not the very first iteration!

          ! Save the old state vector (iteration - 1'th state vector)
          if (divergent_step) then
             old_sv(:) = last_successful_sv(:)
             SV%svsv(:) = old_sv(:)
          else
             old_sv = SV%svsv
          end if

          ! If this is not the first iteration, we grab forward model values from the !
          ! current state vector.

          ! Albedo coefficients grabbed from the SV and used to construct
          ! new albedo(:) array for high-resolution grid.
          if (SV%num_albedo > 0) then
             albedo = 0.0d0
             do i=1, SV%num_albedo
                albedo = albedo + SV%svsv(SV%idx_albedo(i)) &
                     * ((hires_grid(:) - hires_grid(1)) ** (dble(i-1)))
             end do
          endif

          ! If solar parameters are retrieved, update the solar shift and stretch from the
          ! state vector. Otherwise the values stay at 0/1 respectively.
          if (SV%num_solar_shift == 1) this_solar_shift = SV%svsv(SV%idx_solar_shift(1))
          if (SV%num_solar_stretch == 1) this_solar_stretch = SV%svsv(SV%idx_solar_stretch(1))

          ! If we are retrieving surface pressure, we need to re-calculate the
          ! "last" layer if the surface pressure jumps to the next layer.
          if (SV%num_psurf == 1) then
             this_psurf = SV%svsv(SV%idx_psurf(1))

             ! The number of active levels has to be inferred for every
             ! iteration if we are retrieving surface pressure
             do j=1, num_levels
                if (this_psurf > this_atm%p(j)) then
                   num_active_levels = j+1
                end if
             end do

             ! If psurf > BOA p level, we have a problem and thus can't go on.
             if (num_active_levels > num_levels) then
                write(tmp_str, '(A, F12.3, A, F12.3)') "Psurf at ", this_psurf, &
                     " is larger than p(BOA) at ", this_atm%p(size(this_atm%p))
                call logger%error(fname, trim(tmp_str))
                return
             end if

          end if

          do i=1, num_gases
             ! Replace SH if H2O is used in atmosphere
             if (MCS%window(i_win)%gases(i) == "H2O") then
                !    this_atm%sh = this_atm%gas_vmr(:,i) / (SH_H2O_CONV + this_atm%gas_vmr(:,i))
             end if
          end do

       endif

       ! SIF is a radiance, and we want to keep the value for this given
       ! iteration handy for calculations.
       if (SV%num_sif > 0) then
          this_sif_radiance = SV%svsv(SV%idx_sif(1))
       else
          this_sif_radiance = 0.0d0
       end if

       ! Same for ZLO
       if (SV%num_zlo > 0) then
          this_zlo_radiance = SV%svsv(SV%idx_zlo(1))
       else
          this_zlo_radiance = 0.0d0
       end if


       ! Heavy bit - calculate the optical properties given an atmosphere with gases
       ! and their VMRs. This branch of the code will only be entered if we have at least
       ! one gas present. Otherwise, gas_tau will stay unallocated.

       if (num_gases > 0) then

          allocate(gas_tau(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_dpsurf(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_dvmr(N_hires, num_levels, num_gases))
          !allocate(gas_tau_pert(N_hires, num_levels-1, num_levels-1, num_gases))
          allocate(ray_tau(N_hires, num_levels-1))
          allocate(vmr_pert(num_levels))
          allocate(this_vmr_profile(num_levels))
          allocate(ndry(num_levels), ndry_tmp(num_levels))

          ! If we retrieve surface pressure, grab it from the state vector,
          ! otherwise, grab it from the MET data.

          if (SV%num_psurf == 1) then
             this_psurf = SV%svsv(SV%idx_psurf(1))
             do_psurf_jac = .true.
          else
             this_psurf = met_psurf(i_fp, i_fr)
             do_psurf_jac = .false.
          end if

          if (this_psurf < 0) then
             call logger%error(fname, "Psurf negative!")
             return
          end if

          do j=1, num_gases

             ! Copy over this gases' VMR profile
             this_vmr_profile(:) = this_atm%gas_vmr(:,j)
             do_gas_jac = .false.

             ! Enter this branch if we have at least one retrieved gas
             if (SV%num_gas > 0) then

                if ((iteration == 1) .or. (SV%num_psurf == 1)) then

                   ! We need to 'reverse-lookup' to see which SV index belongs to this
                   ! gas to grab the right scaling factor. This is done only on the first
                   ! iteration - or if we retrieve surface pressure.
                   do i=1, SV%num_gas

                      if (MCS%window(i_win)%gas_retrieve_scale(j)) then

                         do_gas_jac = .true.
                         if (SV%gas_idx_lookup(i) == j) then

                            ! This bit here figures out which level/layer range a
                            ! certain scaling factor corresponds to. They are fractions
                            ! of surface pressure, so we use 'searchsorted' to find
                            ! where they would belong to. We also make sure it can't
                            ! go below or above the first/last level.

                            s_start(i) = searchsorted_dp((this_atm%p), &
                                 SV%gas_retrieve_scale_start(i) * (this_psurf), .true.)
                            s_start(i) = max(1, s_start(i))
                            s_stop(i) = searchsorted_dp((this_atm%p), &
                                 SV%gas_retrieve_scale_stop(i) * (this_psurf), .true.) + 1
                            s_stop(i) = min(num_active_levels, s_stop(i))

                         end if
                      end if

                   end do

                   ! We need to make sure that we are not "doubling up" on a specific
                   ! gas VMR level when retrieving scale factors. E.g. 0:0.5 0.5:1.0 will
                   ! produce overlapping s_start/s_stop.

                   do i=1, SV%num_gas
                      do l=1, SV%num_gas
                         ! Skip gases if index does not match
                         if (SV%gas_idx_lookup(i) /= j) cycle
                         if (SV%gas_idx_lookup(l) /= j) cycle

                         if (s_start(i) == s_stop(l)) then
                            s_start(i) = s_start(i) + 1
                         end if
                      end do
                   end do

                   do i=1, SV%num_gas

                      ! Skip gases if index does not match
                      if (SV%gas_idx_lookup(i) /= j) cycle

                      if (s_stop(i) - s_start(i) == 1) then
                         write(tmp_str, '(A,A)') "Scale factor index error for ", results%SV_names(SV%idx_gas(i,1))%chars()
                         call logger%error(fname, trim(tmp_str))
                         write(*,*) s_start(i), s_stop(i), SV%gas_retrieve_scale_start(i), SV%gas_retrieve_scale_stop(i)

                         do l=1, num_levels
                            write(*,*) l, this_atm%p(l), this_atm%p(l) / this_atm%p(num_levels)
                         end do


                         return
                      end if
                   end do

                end if

                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then

                      !write(*,*) "GAS: ", MCS%window(i_win)%gases(j)%chars()
                      !write(*,*) "s_start: ", s_start(i)
                      !write(*,*) "s_stop: ", s_stop(i)

                      ! Finally, apply the scaling factor to the corresponding
                      ! sections of the VMR profile.
                      this_vmr_profile(s_start(i):s_stop(i)) = &
                           this_vmr_profile(s_start(i):s_stop(i)) &
                           * SV%svsv(SV%idx_gas(i,1))
                   end if
                end do

             end if

             ! H2O is treated differently in the gas absorption calculations,
             ! hence we figure out whether this gas is water or not.
             ! CAUTION! This obviously requires that water is actually labelled
             ! H2O, and nothing else. "Standard names" like this need to be
             ! mentioned in the documentation / user guide.
             if (MCS%window(i_win)%gases(j) == "H2O") then
                is_H2O = .true.
             else
                is_H2O = .false.
             end if

             ! Call the function that calculates the gas optical depths
             call calculate_gas_tau( &
                  .true., & ! We are using pre-gridded spectroscopy!
                  is_H2O, & ! Is this gas H2O?
                  hires_grid, & ! The high-resolution wavelength grid
                  this_vmr_profile, & ! The gas VMR profile for this gas with index j
                  this_psurf, & ! Surface pressure
                  this_atm%p(:), & ! Atmospheric profile pressures
                  this_atm%T(:), & ! Atmospheric profile temperature
                  this_atm%sh(:), & ! Atmospheric profile humidity
                  MCS%gas(MCS%window(i_win)%gas_index(j)), & ! MCS%gas object for this given gas
                  MCS%window(i_win)%N_sublayers, & ! Number of sublayers for numeric integration
                  do_psurf_jac, & ! Do we require surface pressure jacobians?
                  do_gas_jac, & ! Do we require gas OD jacobians?
                  gas_tau(:,:,j), & ! Output: Gas ODs
                  gas_tau_dpsurf(:,:,j), & ! Output: dTau/dPsurf
                  gas_tau_dvmr(:,:,j), & ! Output: dTau/dVMR
                  ndry_tmp, & ! Output: ndry molecules per m2
                  success_gas) ! Output: Was the calculation successful?

             ! Only copy ndry over if it's a non-H2O gas
             if (.not. is_H2O) then
                ndry(:) = ndry_tmp(:)
             end if

             ! If the calculation goes wrong, we exit as we can't go on
             if (.not. success_gas) then
                call logger%error(fname, "Error calculating gas optical depths.")
                return
             end if
          end do

          ! Total optical depth is calculated as sum of all gas ODs
          allocate(total_tau(N_hires))
          total_tau(:) = 0.0d0
          total_tau(:) = sum(sum(gas_tau, dim=2), dim=2)

!!$          do i=1, num_active_levels
!!$             write(tmp_str,'(A,G0.1,A,G0.1,A)') "gas_od_iter", iteration, "layer", i, ".dat"
!!$             open(newunit=funit, file=trim(tmp_str))
!!$             do j=1, size(gas_tau, 1)
!!$                write(funit, *) (gas_tau(j,i,l), l=1,size(gas_tau,3))
!!$             end do
!!$             close(funit)
!!$          end do
!!$

!!$          open(newunit=funit, file="gas_od_full.dat")
!!$          do j=1, size(gas_tau, 1)
!!$             write(funit, *) (sum(gas_tau(j,:,l)), l=1,size(gas_tau,3))
!!$          end do
!!$          close(funit)

!!$          write(tmp_str,'(A,G0.1,A)') "gas_dvmr_iter", iteration, ".dat"
!!$          open(newunit=funit, file=trim(tmp_str))
!!$          do j=1, size(gas_tau, 1)
!!$             write(funit, *) (gas_tau_dvmr(j,i,2), i=1,size(gas_tau_dvmr,2))
!!$          end do
!!$          close(funit)

       end if


       ! Calculate the sun-normalized TOA radiances and store them in
       ! 'radiance_calc_work_hi'.
       call calculate_radiance(hires_grid, SZA(i_fp, i_fr), &
            VZA(i_fp, i_fr), albedo, total_tau, &
            radiance_calc_work_hi)

       ! Take a copy of the solar spectrum and re-adjust the solar spectrum wavelength grid
       ! According to both pre-computed solar shift as well as potentially retrieved
       ! solar stretch and shift.
       this_solar(:,1) = this_solar_shift + &
            this_solar_stretch * solar_spectrum_regular(:, 1) / (1.0d0 - solar_doppler)

       call calculate_solar_planck_function(6500.0d0, solar_dist * 1000.0d0, &
            this_solar(:,1), solar_irrad)

       ! And multiply to get the full solar irradiance in physical units
       this_solar(:,2) = solar_spectrum_regular(:, 2) * solar_irrad(:)

       ! If we retrieved either solar shift or stretch (or both), then we
       ! need to know the partial derivative of the solar spectrum w.r.t.
       ! wavelength: dsolar_dlambda
       if ((SV%num_solar_shift == 1) .or. (SV%num_solar_stretch == 1)) then

          ! Sample doppler-shifted and scaled solar spectrum at hires grid
          allocate(solar_tmp(N_hires))
          call pwl_value_1d( &
               N_hires, &
               this_solar(:,1), &
               solar_spectrum_regular(:, 2), &
               N_hires, &
               hires_grid, solar_tmp(:))

          this_solar(:,2) = solar_tmp(:) * solar_irrad(:)
          deallocate(solar_tmp)

          ! If needed, calculate dSolar / dlambda using central differences
          ! apart from the first and last point obviously. This is required for
          ! the solar parameter Jacobians.

          dsolar_dlambda(1) = (this_solar(2,2) - this_solar(1,2)) / &
               (this_solar(2,1) - this_solar(1,1))
          dsolar_dlambda(N_hires) = (this_solar(N_hires-1,2) - this_solar(N_hires,2)) / &
               (this_solar(N_hires-1,1) - this_solar(N_hires,1))

          do i=2, N_hires-1
             dsolar_dlambda(i) = (this_solar(i+1,2) - this_solar(i-1,2)) / &
                  (this_solar(i+1,1) - this_solar(i-1,1))
          end do

       end if


       ! Multiply with the solar spectrum for physical units and add SIF contributions. For
       ! various Jacobian calculations we also need the radiance MINUS the additive contributions,
       ! so we store them in a separate array.
       radiance_calc_work_hi(:) = this_solar(:,2) * radiance_calc_work_hi(:)
       radiance_tmp_work_hi(:) = radiance_calc_work_hi(:) !* (1.0d0 / mu0 + 1.0d0 / mu)
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:) + this_sif_radiance + this_zlo_radiance


       ! JACOBIAN CALCULATIONS

       ! Surface pressure Jacobian
       if (SV%num_psurf == 1) then
          ! This equation requires the TOA radiance before SIF is added, so if we have
          ! SIF in it, take it out beforehand.
          K_hi(:, SV%idx_psurf(1)) = radiance_tmp_work_hi(:) &
               * (1.0d0 / mu0 + 1.0d0 / mu) &
               * (sum(sum(gas_tau_dpsurf, dim=2), dim=2))
       end if

       ! Solar shift jacobian
       if (SV%num_solar_shift == 1) then
          K_hi(:, SV%idx_solar_shift(1)) = -radiance_tmp_work_hi(:) &
               / this_solar(:,2) * dsolar_dlambda(:)
       end if

       ! Solar stretch jacobian
       if (SV%num_solar_stretch == 1) then
          K_hi(:, SV%idx_solar_stretch(1)) = -radiance_tmp_work_hi(:) &
               / this_solar(:,2) * dsolar_dlambda(:) * solar_spectrum_regular(:, 1) / (1.0d0 - solar_doppler)
       end if


       ! Gas jacobians
       if (SV%num_gas > 0) then

          do i=1, SV%num_gas
             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then

                ! This is a scale-type jacobian
                do j=1, size(MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(i),:))

                   ! Loop through all potential profile "sections", but skip the unused ones
                   if (MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(i), j) == -1.0d0) cycle

                   K_hi(:, SV%idx_gas(i,1)) = -radiance_tmp_work_hi(:) * (1.0d0 / mu0 + 1.0d0 / mu) &
                        * sum(gas_tau(:, s_start(i):s_stop(i)-1, SV%gas_idx_lookup(i)), dim=2) / SV%svsv(SV%idx_gas(i,1))

                end do
             end if

          end do
       end if


       ! Stokes coefficients
       ! TODO: apply Stokes coefficients here?
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:)
       K_hi(:,:) = K_hi(:,:)

       ! Grab a copy of the dispersion coefficients from the L1B
       this_dispersion_coefs(:) = dispersion_coefs(:, i_fp, band)

       ! .. and if required, replace L1B dispersion coefficients by state vector
       ! elements from the retrieval process
       if (SV%num_dispersion > 0) then
          do i=1, SV%num_dispersion
             this_dispersion_coefs(i) = SV%svsv(SV%idx_dispersion(i))
          end do
       end if

       ! Calculate the wavelength grid using those dispersion coefficients
       select type(my_instrument)
       type is (oco2_instrument)
          call my_instrument%calculate_dispersion(this_dispersion_coefs, &
               this_dispersion(:), band, i_fp)
       end select

       ! Stretch the (total) dispersion according to the instrument doppler shift
       this_dispersion = this_dispersion / (1.0d0 - instrument_doppler)

       ! Here we grab the index limits for the radiances for
       ! the choice of our microwindow and the given dispersion relation

       call calculate_dispersion_limits(this_dispersion, i_win, l1b_wl_idx_min, l1b_wl_idx_max)

       ! Number of spectral points in the output resolution
       N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1

       ! Allocate various arrays that depend on N_spec
       allocate(K(N_spec, N_sv))
       allocate(gain_matrix(N_sv, N_spec))
       allocate(radiance_meas_work(N_spec))
       allocate(radiance_calc_work(N_spec))
       allocate(radiance_tmp_work(N_spec))
       allocate(noise_work(N_spec))

       ! Grab a copy of the L1b radiances
       radiance_meas_work(:) = radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max)

       ! We add the SIF jacobian AFTER the K matrix is allocated
       if (SV%num_sif > 0) then
          ! Plug in the Jacobians (SIF is easy)
          K(:, SV%idx_sif(1)) = 1.0d0
       end if

       if (SV%num_zlo > 0) then
          K(:, SV%idx_zlo(1)) = 1.0d0
       end if


       ! Now calculate the noise-equivalent radiances
       select type(my_instrument)
       type is (oco2_instrument)

          call my_instrument%calculate_noise( &
               snr_coefs, radiance_meas_work, &
               noise_work, i_fp, band, &
               l1b_wl_idx_min, l1b_wl_idx_max)

          if (any(ieee_is_nan(noise_work))) then
             call logger%error(fname, "NaNs in noise-equivalent-radiance.")
             return
          end if

          ! Pixels flagged with a spike need noise inflation, so
          ! that they are not really considered in the fit. This should
          ! save otherwise good spectra with just a few distorted
          ! radiance values.

          if (allocated(spike_list)) then
             do i=1, N_spec
                ! Remember: spike_list is in full l1b size, but noise work
                ! is bounded by our window choice, thus needs to be offset
                if (spike_list(i + l1b_wl_idx_min - 1, i_fp, i_fr) >= 5) then
                   noise_work(i) = noise_work(i) * 10000.0d0
                end if
             end do
          end if

       end select

       allocate(Se_inv(N_spec, N_spec))
       ! Inverse noise covariance, we keep it diagonal, as usual
       Se_inv(:,:) = 0.0d0
       do i=1, N_spec
          Se_inv(i,i) = 1 / (noise_work(i) ** 2)
       end do

!!$       write(*,*) "Albedo positions: ", SV%idx_albedo
!!$       write(*,*) "Dispersion positions: ", SV%idx_dispersion
!!$       write(*,*) "SIF position: ", SV%idx_sif
!!$       write(*,*) "psurf position: ", SV%idx_psurf
!!$
!!$       open(file="hires_jacs.dat", newunit=funit)
!!$       do i=1, size(this_solar, 1)
!!$          write(funit,*) solar_spectrum_regular(i, 1), (K_hi(i, j), j=1, N_sv)
!!$       end do
!!$       close(funit)
!!$
!!$       open(file="hires_spec.dat", newunit=funit)
!!$       do i=1, N_hires
!!$          write(funit,*) this_solar(i,1), this_solar(i,2), dsolar_dlambda(i), &
!!$               solar_spectrum_regular(i,1), &
!!$               solar_spectrum_regular(i,2), radiance_calc_work_hi(i)
!!$       end do
!!$       close(funit)

       ! Convolution with the instrument line shape function(s)
       ! Note: we are only passing the ILS arrays that correspond to the
       ! actual pixel boundaries of the chosen microwindow.

       select type(my_instrument)
       type is (oco2_instrument)

          ! Convolution of the TOA radiances
          call oco_type_convolution2(hires_grid, radiance_calc_work_hi, &
               ils_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
               ils_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
               this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), radiance_calc_work, &
               ILS_success)

          if (.not. ILS_success) then
             call logger%error(fname, "ILS convolution error.")
             return
          end if

          if (any(ieee_is_nan(radiance_calc_work))) then
             call logger%error(fname, "NaNs in convolved radiance.")
             return
          end if

          ! Convolution of any high-res Jacobians
          do i=1, N_sv

             ! We decide whether we want to go ahead to call the convolution routine
             ! on the various Jacobians. Some Jacobians (SIF, ZLO, albedo, dispersion)
             ! do not require or work with high-res spectra.
             skip_jacobian = .false.


             do j=1, SV%num_dispersion
                ! This is a dispersion Jacobian! Maybe there's a smart way of doing
                ! this analytically, but for now we just perform finite
                ! differencing in a separate loop. So skip this Jacobian.
                if (i == SV%idx_dispersion(j)) then
                   skip_jacobian = .true.
                end if
             end do

             ! SIF and ZLO Jacobians does not need convolution, since it's just 1.0
             if (i == SV%idx_sif(1)) skip_jacobian = .true.
             if (i == SV%idx_zlo(1)) skip_jacobian = .true.

             ! Albedo jacobians can be calculated later as well
             do j=1, SV%num_albedo
                if (i == SV%idx_albedo(j)) then
                   skip_jacobian = .true.
                end if
             end do

             ! If we have any reason to skip this Jacobian index for convolution,
             ! do it now.
             if (skip_jacobian) cycle

             ! Otherwise just convolve the other Jacobians and save the result in
             ! the low-resolution Jacobian matrix 'K'
             call oco_type_convolution2(hires_grid, K_hi(:,i), &
                  ils_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  ils_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), K(:,i), &
                  ILS_success)

             ! Again - check whether convolution was successful.
             if (.not. ILS_success) then
                call logger%error(fname, "ILS convolution error (Jacobians).")
                return
             end if

          end do
       end select

       ! Calculate albedo Jacobians
       ! We could just convolve the high-res albedo jacobian, however the convolution
       ! is probably a good chunk slower than this little section here.
       if (SV%num_albedo > 0) then

          allocate(albedo_low(N_spec))

          ! In order to calculate the low-resolution albedo Jacobian, we also
          ! need the low-resolution albedo, which we calculate here.
          albedo_low(:) = 0.0d0
          do i=1, SV%num_albedo
             albedo_low(:) = albedo_low(:) + SV%svsv(SV%idx_albedo(i)) * &
                  ((this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max) - &
                  this_dispersion(l1b_wl_idx_min)) ** (dble(i-1)))
          end do

          ! And calculate the derivative ..
          do i=1, SV%num_albedo
             K(:, SV%idx_albedo(i)) = (radiance_calc_work(:) - this_sif_radiance - this_zlo_radiance) / albedo_low(:) * &
                  ((this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max) - &
                  this_dispersion(l1b_wl_idx_min)) ** (dble(i-1)))
          end do

          deallocate(albedo_low)
       end if

       ! Disperion Jacobians are produced via finite differencing
       if (SV%num_dispersion > 0) then
          do i=1, SV%num_dispersion

             ! Perturb dispersion coefficient by user-supplied value
             this_dispersion_coefs_pert(:) = this_dispersion_coefs(:)
             this_dispersion_coefs_pert(i) = this_dispersion_coefs_pert(i) &
                  + MCS%window(i_win)%dispersion_pert(i)

             ! And finally do the convolution and create the Jacobian by
             ! differencing / dividing.
             select type(my_instrument)
             type is (oco2_instrument)

                ! Calculate new dispersion grid using perturbed coefficients
                call my_instrument%calculate_dispersion(this_dispersion_coefs_pert, &
                     this_dispersion_tmp(:), band, i_fp)

                ! Apply instrument Doppler shift
                this_dispersion_tmp = this_dispersion_tmp / (1.0d0 - instrument_doppler)

                ! Convolve the perturbed TOA radiance
                call oco_type_convolution2(hires_grid, radiance_calc_work_hi, &
                     ils_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                     ils_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                     this_dispersion_tmp(l1b_wl_idx_min:l1b_wl_idx_max), radiance_tmp_work, &
                     ILS_success)

                if (.not. ILS_success) then
                   call logger%error(fname, "ILS convolution error.")
                   return
                end if

                ! Compute and store the Jacobian into the right place in the matrix
                K(:, SV%idx_dispersion(i)) = (radiance_tmp_work - radiance_calc_work) &
                     / MCS%window(i_win)%dispersion_pert(i)
             end select

          end do
       end if

       ! See Rodgers (2000) equation 5.36: calculating x_i+1 from x_i

       ! K^T Se K
       KtSeK(:,:) = matmul(matmul(transpose(K), Se_inv), K)
       ! (1+gamma) * Sa^-1 + (K^T Se K) 
       tmp_m1 = (1.0d0 + lm_gamma) * Sa_inv + KtSeK

       call invert_matrix(tmp_m1, tmp_m2, success_inv_mat)
       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert K^T Se K")
          return
       end if

       tmp_v1 = matmul(matmul(transpose(K), Se_inv), radiance_meas_work - radiance_calc_work)
       tmp_v2 = matmul(Sa_inv, SV%svsv - SV%svap)

       ! Update state vector
       SV%svsv = SV%svsv + matmul(tmp_m2, tmp_v1 - tmp_v2)

       ! Calculate Shat_inv
       Shat_inv = KtSeK + Sa_inv

       ! Check delta sigma square for this iteration
       dsigma_sq = dot_product(old_sv - SV%svsv, matmul(Shat_inv, old_sv - SV%svsv))

       ! In the case of retrieving gas - we have to adjust the retrieved state vector
       ! if the retrieval wants to push it below 0.

       do i=1, SV%num_gas
          do j=1, size(sv%idx_gas(i,:))
             if (SV%idx_gas(i,j) /= -1) then
                if (SV%svsv(sv%idx_gas(i, j)) < 1.0d-10) then
                   SV%svsv(sv%idx_gas(i, j)) = 1.0d-10
                end if
             end if
          end do
       end do

       if ( &
            (dsigma_sq < dble(N_sv) * dsigma_scale) .or. &
            (iteration > MCS%window(i_win)%max_iterations) .or. &
            (num_divergent_steps > 1) &
            ) then

          ! Stop iterating - we've either coverged or exeeded the max. number of
          ! iterations or max. number of divergences.
          keep_iterating = .false.

          if (iteration <= MCS%window(i_win)%max_iterations) then
             converged = .true.
             results%converged(i_fp, i_fr) = 1
          else
             converged = .false.
             results%converged(i_fp, i_fr) = 0
          end if


          allocate(pwgts(num_active_levels))

          ! Calculate the XGAS for every retreived gas only, same as above with the
          ! gas OD calculation, we loop through all gases, apply the gas scaling
          ! factor from the state vector, and calculate the pressure weighting function
          ! as well as the XGAS.

          if (SV%num_gas > 0) then
             do j=1, num_gases


                ! Here we need to do the same thing as before when calculating
                ! gas OD's. Take local copy of VMR profile, and re-scale the portions
                ! of the profile which, according to the retrieval, have changed.
                this_vmr_profile = this_atm%gas_vmr(:,j)
                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then

                      ! Finally, apply the scaling factor to the corresponding
                      ! sections of the VMR profile.
                      this_vmr_profile(s_start(i):s_stop(i)) = this_vmr_profile(s_start(i):s_stop(i)) &
                           * SV%svsv(SV%idx_gas(i,1))

                      ! We also want to have the corresponding number of molecules of dry air
                      ! for the various sections of the atmopshere.
                      results%ndry(i_fp, i_fr, i) = sum(ndry(s_start(i):s_stop(i)))
                   end if
                end do

                ! Based on this 'current' VMR profile
                call pressure_weighting_function( &
                     this_atm%p(1:num_active_levels), &
                     this_psurf, &
                     this_vmr_profile(1:num_active_levels), &
                     pwgts)

                ! Compute XGAS as the sum of pgwts times GAS VMRs.
                results%xgas(i_fp, i_fr, j) =  sum( &
                     pwgts(:) &
                     * this_vmr_profile(1:num_active_levels) &
                     )

             end do
          end if

          ! Save the final dSigma-squared value (in case anyone needs it)
          results%dsigma_sq(i_fp, i_fr) = dsigma_sq

          ! Calculate Shat from Shat_inverse
          call invert_matrix(Shat_inv, Shat, success_inv_mat)

          if (.not. success_inv_mat) then
             call logger%error(fname, "Failed to invert Shat^-1")
             return
          end if

          ! Calculate the Gain matrix
          gain_matrix(:,:) = matmul(matmul(Shat(:,:), transpose(K)), Se_inv)

          ! Calculate the averaging kernel
          AK(:,:) = matmul(Shat, KtSeK) !matmul(gain_matrix, K)

          ! Calculate state vector element uncertainties from Shat
          do i=1, N_sv
             SV%sver(i) = sqrt(Shat(i,i))
          end do

          do i=1, N_sv
             write(*,*) i, AK(i,i), SV%svsv(i), SV%sver(i), 100.0d0 * (SV%sver(i) / sqrt(Sa(i,i))), results%sv_names(i)%chars()
          end do

          ! Put the SV uncertainty into the result container
          results%sv_uncertainty(i_fp, i_fr, :) = SV%sver(:)

          ! Save retrieved CHI2 - this is the predicted chi2 from the 'last+1'
          ! linear step that is not evaluated.
          results%chi2(i_fp, i_fr) = linear_prediction_chi2

          ! Get an SNR (mean and std) estimate
          results%SNR(i_fp, i_fr) = mean(radiance_meas_work / noise_work)
          results%SNR_std(i_fp, i_fr) = std(radiance_meas_work / noise_work)

          ! Save also the continuum level radiance
          results%continuum(i_fp, i_fr) = percentile(radiance_meas_work, 99.0d0)

          ! Save statevector at last iteration
          do i=1, size(SV%svsv)
             results%sv_retrieved(i_fp, i_fr,i) = SV%svsv(i)
          end do

          ! Save the number of iterations
          results%num_iterations(i_fp, i_fr) = iteration

          if (MCS%output%save_radiances) then
             ! Save the radiances and noises - if required by the user
             final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_calc_work(:)
             measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_meas_work(:)
             noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)
          end if

          ! Write a debug message
          write(tmp_str, '(A,G3.1,A,F6.1,A,G2.1,A,F10.2,A,E10.3,A,F10.2)') "Iteration: ", iteration ,&
               ", Chi2: ",  results%chi2(i_fp, i_fr), &
               ", Active Levels: ", num_active_levels, &
               ", Psurf: ", this_psurf, &
               ", LM-Gamma: ", lm_gamma, &
               ", SNR: ", results%SNR(i_fp, i_fr)
          call logger%debug(fname, trim(tmp_str))

       else
          ! Not converged yet - change LM-gamma according to the chi2_ratio.

          !! See Rogers (2000) Section 5.7 - if Chi2 increases as a result of the iteration,
          !! we revert to the last state vector, and increase lm_gamma

          if (.not. divergent_step) then
             ! Make sure we keep the 'old' Chi2 only from a valid iteration
             old_chi2 = this_chi2
          end if

          this_chi2 = calculate_chi2(radiance_meas_work, radiance_calc_work, noise_work, N_spec - N_sv)

          if (iteration == 1) then
             chi2_ratio = 0.5d0
             ! For very first iteration, we set R in the trust region, so
             ! lm_gamma will not be changed.
          else
             ! Otherwise R is ratio between linear forecast value and actual computed Chi2
             chi2_ratio = (old_chi2 - this_chi2) / (old_chi2 - linear_prediction_chi2)
          end if

          ! Make linear prediction of CHI2 for the next iteration
          radiance_tmp_work(:) = radiance_calc_work(:) + matmul(K, SV%svsv(:) - old_sv(:))
          linear_prediction_chi2 = calculate_chi2(radiance_meas_work, radiance_tmp_work, &
               noise_work, N_spec - N_sv)

          if (MCS%window(i_win)%allow_divergences) then
             ! Do we allow for divergent steps and the per-iteration adjustment of
             ! the LM-gamma parameter?

             if (chi2_ratio > 0.75d0) then
                ! If fit is much better than predicted / fairly linear - decrese gamma by a bit
                lm_gamma = lm_gamma * 0.5d0
                divergent_step = .false.
                ! If we have a successful iteration, we want to store the last
                ! successful one here as well.
                last_successful_sv(:) = SV%svsv(:)
             else if ((chi2_ratio < 0.25d0) .and. (chi2_ratio > 1.0d-4)) then
                ! Fit is not as good as predicted so increase gamma
                lm_gamma = lm_gamma * 10.0d0
                divergent_step = .false.
                ! If we have a successful iteration, we want to store the last
                ! successful one here as well.
                last_successful_sv(:) = SV%svsv(:)
             else if (chi2_ratio <= 1.0d-4) then
                ! We consider this a divergence!
                divergent_step = .true.
                num_divergent_steps = num_divergent_steps + 1

                lm_gamma = lm_gamma * 10.0d0

                call logger%debug(fname, "Divergent step!")

             else
                ! Otherweise, we are in the trust region, do nothing with gamma
                divergent_step = .false.
             end if
          end if

       end if

       open(file="jacobian.dat", newunit=funit)
       do i=1, N_spec
          write(funit,*) (K(i, j), j=1, N_sv)
       end do
       close(funit)
!!$
!!$
!!$       open(file="l1b_spec.dat", newunit=funit)
!!$       do i=1, N_spec
!!$          write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), radiance_calc_work(i), &
!!$               noise_work(i)!, solar_low(i)
!!$
!!$       end do
!!$       close(funit)
!!$
!!$
!!$       write(*,*) num_active_levels
!!$       do i=1, num_active_levels
!!$          write(*,*) this_atm%p(i), (this_atm%gas_vmr(i,j), j=1, size(this_atm%gas_vmr, 2)), ndry(i)
!!$       end do

!!$       write(*,*) "old, current and delta state vector, and errors"
!!$
!!$       write(*,*) "Iteration: ", iteration
!!$       do i=1, N_sv
!!$          write(*, '(I3.1,A40,ES15.6,ES15.6,ES15.6,ES15.6)') &
!!$               i, results%sv_names(i)%chars(), old_sv(i), SV%svsv(i), &
!!$               SV%svsv(i) - old_sv(i), sqrt(Shat(i,i))
!!$       end do
!!$       write(*,*) "Chi2:    ", SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)
!!$       write(*,*) "Dsigma2: ", dsigma_sq, '/', dble(N_sv) * dsigma_scale
!!$       write(*,*) "LM-Gamma: ", lm_gamma
!!$       write(*,*) "Ratio R: ", chi2_ratio
!!$       write(*,*) this_sif_radiance
!!$       write(*,*) this_zlo_radiance


       ! These quantities are all allocated within the iteration loop, and
       ! hence need explicit de-allocation.
       deallocate(radiance_meas_work, radiance_calc_work, radiance_tmp_work, &
            noise_work, Se_inv, K, gain_matrix, solar_irrad)

       if (allocated(gas_tau)) deallocate(gas_tau)
       if (allocated(gas_tau_dpsurf)) deallocate(gas_tau_dpsurf)
       if (allocated(gas_tau_dvmr)) deallocate(gas_tau_dvmr)
       if (allocated(gas_tau_pert)) deallocate(gas_tau_pert)
       if (allocated(ray_tau)) deallocate(ray_tau)
       if (allocated(total_tau)) deallocate(total_tau)
       if (allocated(vmr_pert)) deallocate(vmr_pert)
       if (allocated(this_vmr_profile)) deallocate(this_vmr_profile)
       if (allocated(ndry)) deallocate(ndry)
       if (allocated(ndry_tmp)) deallocate(ndry_tmp)

       !read(*,*)
    end do

  end function physical_FM


  !> @brief Set the entries of the prior covariance matrix (and inverse)
  !>
  !> @param SV State vector object
  !> @param continuum Continuum level radiance (needed for SIF prior cov)
  !> @param i_win Window index for MCS
  !> @param Sa Prior covariance matrix
  !> @param Sa_inv Inverse of prior covariance matrix

  subroutine populate_prior_covariance(SV, continuum, i_win, Sa, Sa_inv)

    implicit none

    type(statevector), intent(in) :: SV
    double precision, intent(in) :: continuum
    integer, intent(in) :: i_win
    double precision, intent(inout) :: Sa(:,:)
    double precision, intent(inout) :: Sa_inv(:,:)

    ! Local variables
    integer :: i

    ! Make the albedo essentially unconstrained. We keep this larger than 1, simply
    ! because we can end up with effective albedos much larger than 1.0 due to
    ! non-Lambertian BRDFs.
    if (SV%num_albedo > 0) then
       do i=1, SV%num_albedo
          if (i==1) then
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 10.0d0
          else
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 10.0d0 ** dble(i)
          end if
       end do
    end if

    if (SV%num_solar_shift == 1) then
       Sa(SV%idx_solar_shift(1), SV%idx_solar_shift(1)) = 1.0d-2
    end if

    if (SV%num_solar_stretch == 1) then
       Sa(SV%idx_solar_stretch(1), SV%idx_solar_stretch(1)) = 1.0d0
    end if

    if (SV%num_sif > 0) then
       ! Put SIF prior covariance at the continuum level of the band
       Sa(SV%idx_sif(1), SV%idx_sif(1)) = continuum * continuum
    end if

    if (SV%num_zlo > 0) then
       ! Put ZLO prior covariance at the continuum level of the band
       Sa(SV%idx_zlo(1), SV%idx_zlo(1)) = continuum * continuum
    end if

    if (SV%num_psurf == 1) Sa(SV%idx_psurf(1), SV%idx_psurf(1)) = 1.0d6

    ! Make the dispersion prior covariance about ten times the perturbation value
    if (SV%num_dispersion > 0) then
       do i=1, SV%num_dispersion
          Sa(SV%idx_dispersion(i), SV%idx_dispersion(i)) = &
               (10.0d0 * MCS%window(i_win)%dispersion_pert(i))**2
       end do
    end if

    ! The scale factor covariances are taken from the configuration file / MCS
    do i=1, SV%num_gas
       if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
          Sa(SV%idx_gas(i,1), SV%idx_gas(i,1)) = &
               (SV%gas_retrieve_scale_cov(i) * SV%gas_retrieve_scale_cov(i))
       end if
    end do

    ! At the moment, Sa will be diagonal, so inverting is trivial
    do i=1, size(Sa, 1)
       Sa_inv(i,i) = 1.0d0 / Sa(i,i)
    end do

    !! If at some point we want to have non-diagonal elements in
    !! Sa, just uncomment the lines below and a proper inversion is being done.

    !call invert_matrix(Sa, Sa_inv, success_inv_mat)
    !if (.not. success_inv_mat) then
    !   call logger%error(fname, "Failed to invert Sa")
    !   return
    !end if

  end subroutine populate_prior_covariance



  !> @brief Given a dispersion array "this_dispersion", in window "i_win", this
  !> function calculates the first and last pixel indices as the boundaries in
  !> the detector.
  !> @param this_dispersion Wavelength per detector index
  !> @param i_win Retrieval Window index for MCS
  !> @param l1b_wl_idx_min Lower wavelength pixel index corresponding to
  !> user-defined "wl_min"
  !> @param l1b_wl_idx_min Upper wavelength pixel index corresponding to
  !> user-defined "wl_max"
  subroutine calculate_dispersion_limits(this_dispersion, i_win, &
       l1b_wl_idx_min, l1b_wl_idx_max)

    implicit none
    double precision, intent(in) :: this_dispersion(:)
    integer, intent(in) :: i_win
    integer, intent(inout) :: l1b_wl_idx_min, l1b_wl_idx_max

    integer :: i

    l1b_wl_idx_min = 0
    l1b_wl_idx_max = 0

    ! This here grabs the boundaries of the L1b data by simply looping over
    ! the dispersion array and setting the l1b_wl_idx_* accordingly.
    do i=1, size(this_dispersion)
       if (this_dispersion(i) < MCS%window(i_win)%wl_min) then
          l1b_wl_idx_min = i
       end if
       if (this_dispersion(i) < MCS%window(i_win)%wl_max) then
          l1b_wl_idx_max = i
       end if
       if (this_dispersion(i) > MCS%window(i_win)%wl_max) then
          exit
       end if
    end do

    ! If window lower limit is below the first wavelength value, set it
    ! to the beginning (index 1)
    if (l1b_wl_idx_min == 0) then
       l1b_wl_idx_min = 1
    end if

    ! Have to increase the higher-wavelength index by one so that we can
    ! use the full stride l1b_wl_idx_min: l1b_wl_idx_min to access the L1b
    ! radiance corresponding to the user-defined values.
    if (l1b_wl_idx_max < size(this_dispersion)) then
       l1b_wl_idx_max = l1b_wl_idx_max + 1
    end if

    ! If the index goes past the maximal size of the array, simply
    ! set it back to the boundary.
    if (l1b_wl_idx_max > size(this_dispersion)) then
       l1b_wl_idx_max = size(this_dispersion)
    end if

  end subroutine calculate_dispersion_limits



  !> @brief Creates / allocates the "results" container to hold all retrieval results
  !>
  !> @param results Result container
  !> @param num_frames Number of frames
  !> @param num_fp Number of footprints
  !> @param num_SV Number of state vector elements
  !> @param num_gas Number of retrieved (!) gases
  subroutine create_result_container(results, num_frames, num_fp, num_SV, num_gas)
    implicit none
    type(result_container), intent(inout) :: results
    integer, intent(in) :: num_frames, num_fp, num_SV, num_gas

    allocate(results%sv_names(num_SV))

    allocate(results%sv_retrieved(num_fp, num_frames, num_SV))
    allocate(results%sv_prior(num_fp, num_frames, num_SV))
    allocate(results%sv_uncertainty(num_fp, num_frames, num_SV))
    allocate(results%xgas(num_fp, num_frames, num_gas))

    allocate(results%chi2(num_fp, num_frames))
    allocate(results%residual_rms(num_fp, num_frames))
    allocate(results%dsigma_sq(num_fp, num_frames))
    allocate(results%SNR(num_fp, num_frames))
    allocate(results%SNR_std(num_fp, num_frames))
    allocate(results%continuum(num_fp, num_frames))

    allocate(results%num_iterations(num_fp, num_frames))
    allocate(results%converged(num_fp, num_frames))

    allocate(results%ndry(num_fp, num_frames, num_gas))

    results%sv_names = "NONE"

    ! This might cause problems for some, but I find it convenient
    ! to set 'unused' fields to NaNs. Remember this only works for
    ! reals/double precision, but not for integers.
    results%sv_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_uncertainty = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%chi2 = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%residual_rms = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%dsigma_sq = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR_std = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%continuum = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

    results%num_iterations = -1
    results%converged = -1

  end subroutine create_result_container


  !> @brief Destroys the "results" container allocated in "create_result_container"
  !>
  !> @param results Result container
  subroutine destroy_result_container(results)
    implicit none
    type(result_container), intent(inout) :: results

    deallocate(results%sv_names)
    deallocate(results%sv_retrieved)
    deallocate(results%sv_prior)
    deallocate(results%sv_uncertainty)
    deallocate(results%xgas)
    deallocate(results%chi2)
    deallocate(results%residual_rms)
    deallocate(results%dsigma_sq)
    deallocate(results%SNR)
    deallocate(results%SNR_std)
    deallocate(results%continuum)
    deallocate(results%num_iterations)
    deallocate(results%converged)
    deallocate(results%ndry)

  end subroutine destroy_result_container



  !> @brief Creates human-readable names for state vector elements
  !>
  !> The idea is faily simple: we loop through all the state vector elements,
  !> and then check for each one if there is a corresponding SV\%idx_* associated with
  !> that element position. Based on that, we create a name for the state vector
  !> element, which usually has the parameter number (e.g. albedo order) baked in.
  !> @param results Result container
  !> @param SV State vector object
  !> @param i_win Retrieval window index for MCS
  subroutine assign_SV_names_to_result(results, SV, i_win)
    implicit none
    type(result_container), intent(inout) :: results
    type(statevector), intent(in) :: SV
    integer, intent(in) :: i_win

    character(len=999) :: tmp_str
    integer :: i,j,k

    i = 1
    do while (i <= size(SV%svsv))

       ! Albedo names
       do j=1, SV%num_albedo
          if (SV%idx_albedo(j) == i) then
             write(tmp_str, '(A,G0.1)') "albedo_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! SIF name (really only one at this point)
       do j=1, SV%num_sif
          if (SV%idx_sif(j) == i) then
             write(tmp_str, '(A)') "SIF_absolute"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       do j=1, SV%num_zlo
          if (SV%idx_zlo(j) == i) then
             write(tmp_str, '(A)') "ZLO"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Solar shift name
       if (SV%idx_solar_shift(1) == i) then
          write(tmp_str, '(A)') "solar_shift"
          results%sv_names(i) = trim(tmp_str)
       end if

       ! Solar stretch name
       if (SV%idx_solar_stretch(1) == i) then
          write(tmp_str, '(A)') "solar_stretch"
          results%sv_names(i) = trim(tmp_str)
       end if

       ! Surface pressure name
       do j=1, SV%num_psurf
          if (SV%idx_psurf(j) == i) then
             write(tmp_str, '(A)') "psurf"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Dispersion parameter names
       do j=1, SV%num_dispersion
          if (SV%idx_dispersion(j) == i) then
             write(tmp_str, '(A,G0.1)') "dispersion_coef_", j
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved gas names (scalar retrieval only so far)
       k = 1
       do j=1, SV%num_gas

          ! Check if this SV element is a scalar retrieval
          if (SV%idx_gas(j,1) == i) then

             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(j))) then
                if (MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(j), k) == -1.0) cycle

                write(tmp_str, '(A,A)') trim(MCS%window(i_win)%gases(sv%gas_idx_lookup(j))%chars() // "_scale_")
                write(tmp_str, '(A, F4.2)') trim(tmp_str), &
                     SV%gas_retrieve_scale_start(j)
                write(tmp_str, '(A,A,F4.2)') trim(tmp_str), "_" , &
                     SV%gas_retrieve_scale_stop(j)
                results%sv_names(i) = trim(tmp_str)

                k = k + 1
             end if

          end if

       end do

       i = i+1
    end do

  end subroutine assign_SV_names_to_result

  !> @brief Reads in the atmosphere file, which contains column-based profiles.
  !>
  !> @param filename Path to the atmosphere file
  !> @param Array of 'strings' containing names of required gases
  !> @param Atmosphere object that will be populated by this function
  !>
  !> @detail We supply here a filename, a list of gas strings and the
  !> atmosphere-object. The function checks whether all required gases
  !> are present in the atmosphere file, and will throw a fatal error if
  !> that is not the case.
  subroutine read_atmosphere_file(filename, gas_strings, atm)

    implicit none
    character(len=*), intent(in) :: filename
    type(string), intent(in) :: gas_strings(:)
    type(atmosphere), intent(inout) :: atm

    ! Function name
    character(len=*), parameter :: fname = "read_atmosphere_file"
    ! File handler and IO stat variable
    integer :: funit, iostat
    ! Whether file exists
    logical :: file_exist
    ! Various counting variables to figure out where the
    ! contents of the atmosphere file start.
    integer :: line_count, nonempty_line_count, file_start, level_count
    ! Indices which tell us which columns are pressure and temperature
    integer :: idx_p, idx_t
    ! Character variables to store lines
    character(len=999) :: dummy, tmp_str
    ! This dummy is for reading in numerical values
    double precision :: dummy_dp
    ! This dummy is for reading in strings
    type(string) :: dummy_string
    ! Lines will be split using this dummy variable
    type(string), allocatable :: split_string(:)
    ! Various loop counters
    integer :: i, j, cnt
    ! Too keep track of whether we have found the required gases, and
    ! at which position/column they are.
    integer, allocatable :: this_gas_index(:)

    integer :: num_gases

    ! Check whether the file exists or not.
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Atmosphere file does not exist: " // filename)
       stop 1
    else
       call logger%debug(fname, "File does exist.")
    end if

    ! First pass: we scan the file to see how many levels our atmosphere has
    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    level_count = 0
    nonempty_line_count = 0
    file_start = -1

    ! Loop through the file until we have reached the end
    do
       read(funit, '(A)', iostat=iostat) dummy

       if (iostat == iostat_end) then
          ! End of file?
          exit
       end if

       ! Keep track of how many lines we have traversed
       line_count = line_count + 1

       if (scan(dummy, "!#;") > 0) then
          ! Skip this line, as it's commented
          cycle
       else if (trim(dummy) == "") then
          ! Skip empty lines
          cycle
       end if

       ! And keep track of now many lines are non-empty
       nonempty_line_count = nonempty_line_count + 1

       ! First non-empty line should be saved here.
       ! file_start is where the real file contents begin
       if (file_start == -1) then
          file_start = line_count
       else
          if (nonempty_line_count > 1) then
             level_count = level_count + 1
          end if
       end if

    end do

    ! Go back to the top of the file.
    rewind(unit=funit, iostat=iostat)

    ! Reset the line counter since we're starting from the top again.
    line_count = 0
    do
       ! Read the line
       read(funit, '(A)', iostat=iostat) dummy
       line_count = line_count + 1

       ! .. and immediately skip until we are at the
       ! beginning of the contents of the file.
       if (line_count < file_start) cycle

       if (line_count == file_start) then
          ! This is the proper atmosphere header that should contain
          ! the information about the gases. So first, let's check if
          ! the numbers match

          ! Read the line into a string object and split it by whitespaces
          dummy_string = dummy
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now that we know both the number of levels and gases, we can allocate the
          ! arrays in the atmosphere structure.
          num_gases = 0
          idx_p = -1
          idx_t = -1
          do j=1, size(split_string)
             ! Skip temp or pressure - this requires the pressure and temperature
             ! columns to be explicitly labeled p/P and t/T
             if (split_string(j)%lower() == "p") then
                idx_p = j
                cycle
             end if
             if (split_string(j)%lower() == "t") then
                idx_t = j
                cycle
             end if
             num_gases = num_gases + 1
          end do

          ! Let the user know how many gases and levels we have found
          write(tmp_str, '(A,G0.1,A,A)') "There seem to be ", num_gases, " gases in ", filename
          call logger%info(fname, trim(tmp_str))
          write(tmp_str, '(A, G0.1)') "The number of atmospheric levels is: ", level_count
          call logger%info(fname, trim(tmp_str))

          ! If the atmosphere structure was already allocated, deallocate it first
          if (allocated(atm%p)) deallocate(atm%p)
          if (allocated(atm%T)) deallocate(atm%T)
          if (allocated(atm%sh)) deallocate(atm%sh)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)

          ! Allocate according to the file structure
          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))
          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))

          allocate(this_gas_index(size(gas_strings)))

          ! But we also want to know what gas index to use for storage
          do i=1, size(gas_strings)
             this_gas_index(i) = -1
             cnt = 1
             do j=1, size(split_string)
                ! Skip temp or pressure
                if (split_string(j)%lower() == "p") cycle
                if (split_string(j)%lower() == "t") cycle
                ! Check if gas description matches gases we know
                if (split_string(j) == gas_strings(i)) then
                   this_gas_index(i) = j
                   write(tmp_str, '(A,A,A,G0.1)') "Index for atmosphere gas ", &
                        split_string(j)%chars(), ": ", j
                   call logger%debug(fname, trim(tmp_str))
                   exit
                end if
                cnt = cnt + 1
             end do
          end do

          ! And last check - do all required gases have a VMR column in the
          ! atmosphere file.

          do i=1, size(gas_strings)
             if (this_gas_index(i) == -1) then
                ! Uh-oh, gas that was speficied in the "gases" option of the window
                ! could not be found in the atmosphere! Hard exit.
                write(tmp_str, "(A,A)") "The following gas was not found in the atmosphere file: " &
                     // gas_strings(i)%chars()
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if
          end do

       end if


       if (line_count > file_start) then
          ! Right after the header, we should have the data in rows. We
          ! use the string split option to split the row string into substrings,
          ! and then convert each into double precision values and feed them into
          ! the atmosphere structure - at the right position!

          dummy_string = dummy
          ! Need to deallocate the split_string object first
          if (allocated(split_string)) deallocate(split_string)
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now here we need to check again whether a certain line has more than
          ! num_gases+1 columns.
          if (size(split_string) /= (num_gases + 2)) then
             write(tmp_str, '(A, G0.1)') "Too many values in line ", line_count
             call logger%fatal(fname, trim(tmp_str))
             stop 1
          end if

          ! Get the pressure value
          tmp_str = split_string(idx_p)%chars()
          read(tmp_str, *) dummy_dp
          atm%p(line_count - file_start) = dummy_dp

          ! Get the temperature value
          tmp_str = split_string(idx_t)%chars()
          read(tmp_str, *) dummy_dp
          atm%T(line_count - file_start) = dummy_dp

          ! And the gas value(s) from the other column(s)
          do i=1, size(gas_strings)
             tmp_str = split_string(this_gas_index(i))%chars()
             read(tmp_str, *) dummy_dp
             atm%gas_vmr(line_count - file_start, i) = dummy_dp
          end do

       end if

       if (line_count == (file_start + level_count)) exit

    end do

    close(funit)

  end subroutine read_atmosphere_file

end module physical_model_mod
