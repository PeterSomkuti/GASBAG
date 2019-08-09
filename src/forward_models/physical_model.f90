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
  use control_mod, only: MCS, MAX_WINDOWS, MAX_GASES, MAX_AEROSOLS, MCS_find_gases
  use instruments_mod, only: generic_instrument
  use oco2_mod
  use solar_model_mod
  use math_utils_mod
  use statevector_mod
  use Beer_Lambert_mod
  use XRTM_mod
  use absco_mod
  use gas_tau_mod
  use spectroscopy_utils_mod
  use smart_first_guess_mod
  use aerosols_mod
  use jacobians_mod
  use Rayleigh_mod

  ! Third-party modules
  use logger_mod, only: logger => master_logger
  use mod_datetime
  use doppler_solar_module
  use stringifor

  ! XRTM
  use xrtm_int_f90

  ! System modules
  use ISO_FORTRAN_ENV
  USE OMP_LIB
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

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
  !> Sounding IDs from MET file (integer type 8 to accomodate OCO-2 for now), (footprint, frame)
  integer(8), allocatable :: sounding_ids_met(:,:)
  !> Sounding time strings
  character(len=25), allocatable :: frame_time_strings(:)
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
  !> L1B Stokes coefficients
  double precision, allocatable :: stokes_coef(:,:,:)
  !> List of 'bad'-flagged detector pixels
  integer, allocatable :: bad_sample_list(:,:,:)
  !> If required, an array to hold the spike value (pixel, footprint, band)
  integer, allocatable :: spike_list(:,:,:)

  !> MET data temperature profiles (level, footprint, frame)
  double precision, allocatable :: met_T_profiles(:,:,:)
  !> MET data pressure levels (level, footprint, frame)
  double precision, allocatable :: met_P_levels(:,:,:)
  !> MET data specific humidity profiles (level, footprint, frame)
  double precision, allocatable :: met_SH_profiles(:,:,:)
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
  double precision, allocatable :: solar_continuum_from_hdf(:,:)
  integer :: N_solar

  !> Final modelled radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: final_radiance
  !> Measured radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: measured_radiance
  !> Calculated noise radiances (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: noise_radiance
  !> Wavelength array (pixel, footprint, frame)
  double precision, dimension(:,:,:), allocatable :: wavelength_radiance

  !> Global state vector construct
  type(statevector) :: global_SV
  !> Result container
  type(result_container) :: results

  !> RT model Beer-Lambert
  integer, parameter :: RT_BEER_LAMBERT = 0
  !> RT model XRTM
  integer, parameter :: RT_XRTM = 1


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
    ! Does a spike/bad_sample data field exist?
    logical :: spike_exists, bad_sample_exists, result_exists, met_sounding_exists
    ! Fixed dimensions for the arrays to be saved into the output HDF file
    integer(hsize_t) :: out_dims2d(2), out_dims3d(3)
    ! Dimensions for the read-in of MET/L1b data
    integer(hsize_t), allocatable :: dset_dims(:)
    ! HDF error variable
    integer :: hdferr
    ! Holders for temporary strings and dataset names to grab L1b/MET data
    character(len=999) :: dset_name, tmp_str, group_name
    ! String needed for converting to lower-case
    type(string) :: lower_str
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
    ! Thread number for openMP
    integer :: this_thread
    ! Retrieval count
    integer :: retr_count, total_number_todo
    ! Return value of physical_fm, tells us whether this one has converged or not
    logical :: this_converged
    ! File unit for debugging only..
    integer :: funit
    ! CPU time stamps and mean duration for performance analysis
    double precision :: cpu_time_start, cpu_time_stop, mean_duration
    integer :: frame_start, frame_skip, frame_stop

    double precision, allocatable :: land_fraction(:,:)

    ! Open up the MET file
    if (MCS%algorithm%observation_mode == "downlooking") then
       call h5fopen_f(MCS%input%met_filename%chars(), &
            H5F_ACC_RDONLY_F, MCS%input%met_file_id, hdferr)
       call check_hdf_error(hdferr, fname, "Error opening MET file: " &
            // trim(MCS%input%met_filename%chars()))
       ! Store HDF file handler for more convenient access
       met_file_id = MCS%input%met_file_id
    end if

    ! Store HDF file handler for more convenient access
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id
    this_thread = 0

    ! Grab number of frames and footprints
    num_frames = MCS%general%N_frame
    num_fp = MCS%general%N_fp
    ! Grab number of bands
    num_band = MCS%general%N_bands

    !---------------------------------------------------------------------
    ! INSTRUMENT-DEPENDENT SET UP OF L1B AND MET DATA
    !---------------------------------------------------------------------


    ! Read in MET and L1B data, this is instrument-dependent, so we need to
    ! make our instrument select here as well.
    select type(my_instrument)
    type is (oco2_instrument)

       call my_instrument%read_MET_data(met_file_id, l1b_file_id, &
            met_P_levels, met_T_profiles, met_SH_profiles, met_psurf)

       ! Pressure levels of p=0 don't agree well with the rest of the code,
       ! especially when logarithms are calculated. So I'm replacing them
       ! right here with a small value.
       ! (This really shouldn't occur anyway)
       where (met_P_levels < 1d-10) met_P_levels = 1d-10

       if (MCS%algorithm%observation_mode == "downlooking") then
          ! These here also only make sense in a downlooking position.
          ! Read in the sounding id's
          call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)

          ! We can check if the MET and L1B match up by comparing the sounding_id
          ! field in both files - which should be identical!
          call h5lexists_f(met_file_id, "/SoundingGeometry/sounding_id", &
               met_sounding_exists, hdferr)
          if (met_sounding_exists) then
             call my_instrument%read_sounding_ids(met_file_id, sounding_ids_met)

             if (all(sounding_ids == sounding_ids_met)) then
                call logger%debug(fname, "L1B and MET have same Sounding IDs - good.")
             else
                call logger%error(fname, "L1B and MET have differing Sounding IDs!")
                stop 1
             end if
          else
             call logger%debug(fname, "MET does not have Sounding IDs - cannot verify if files match.")
          end if

          if (allocated(land_fraction)) deallocate(land_fraction)
          dset_name = "/SoundingGeometry/sounding_land_fraction"
          call read_DP_hdf_dataset(l1b_file_id, dset_name, land_fraction, dset_dims)
          call logger%trivia(fname, "Finished reading in land fraction.")
       end if

       ! Grab the SNR coefficients for noise calculations
       call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
       ! Read the time strings
       call my_instrument%read_time_strings(l1b_file_id, frame_time_strings)
       ! Read in the instrument ILS data
       call my_instrument%read_ils_data(l1b_file_id, ils_delta_lambda, &
            ils_relative_response)
       ! Read in bad sample list (if it exists)
       call h5lexists_f(l1b_file_id, "/InstrumentHeader/bad_sample_list", &
            bad_sample_exists, hdferr)
       if (bad_sample_exists) then
          call my_instrument%read_bad_sample_list(l1b_file_id, bad_sample_list)
       end if

       ! Read dispersion coefficients and create dispersion array
       call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)

       ! Grab the number of spectral pixels in this band - we need to grab the largest value
       ! of all bands here. TODO: This will cause problems if we have different number of
       ! pixels for the different bands!
       num_pixel = 0
       do band=1, num_band
          if (MCS%general%N_spec(band) > num_pixel) then
             num_pixel = MCS%general%N_spec(band)
          end if
       end do

       allocate(dispersion(num_pixel, num_fp, num_band))

       do band=1, num_band
          do i_fp=1, num_fp
             call my_instrument%calculate_dispersion(dispersion_coefs(:, i_fp,&
                  & band), dispersion(:, i_fp, band), band, i_fp)
          end do
       end do

    end select

    ! We can also read in the aerosol data right here, as they are the same for
    ! all retrieval windows.
    do i=1, MAX_AEROSOLS
       if (.not. MCS%aerosol(i)%used) cycle
       call ingest_aerosol_files(MCS%aerosol(i))
    end do

    ! Create the HDF group in which all the results go in the end
    call h5lexists_f(l1b_file_id, "/RetrievalResults", result_exists, hdferr)
    if (.not. result_exists) then
       call h5gcreate_f(output_file_id, "RetrievalResults", result_gid, hdferr)
    end if
    call check_hdf_error(hdferr, fname, "Error. Could not create group: RetrievalResults")
    call h5gcreate_f(output_file_id, "/RetrievalResults/physical", result_gid, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not create group: RetrievalResults/physical")



    !---------------------------------------------------------------------
    ! MAIN RETRIEVAL WINDOW LOOP
    !---------------------------------------------------------------------

    ! Loop over all potential user-defined retrieval windows.
    do i_win=1, MAX_WINDOWS

       ! Just skip unused windows
       if (.not. MCS%window(i_win)%used) cycle

       ! At the beginning, we check which gases were defined for this window,
       ! and see if a gas with the corresponding name has been defined.
       call MCS_find_gases(MCS%window, MCS%gas, i_win)

       ! If we have gases, we want to read in the corresponding spectroscopy data
       ! We have to do this for every microwindow since the cross sections are being
       ! re-gridded for every microwindow.
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
       ! Grab the number of spectral pixels in this band
       num_pixel = MCS%general%N_spec(band)

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
          ! Grab the L1B stokes coefficients
          call my_instrument%read_stokes_coef(l1b_file_id, band, stokes_coef)
          ! Read in Spike filter data, if it exists in this file
          call h5lexists_f(l1b_file_id, "/SpikeEOF", spike_exists, hdferr)
          if (spike_exists) then
             if (allocated(spike_list)) deallocate(spike_list)
             call my_instrument%read_spike_filter(l1b_file_id, spike_list, band)
          end if

       end select

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
          call read_toon_solar_spectrum(MCS%algorithm%solar_file%chars(), &
               solar_spectrum, &
               MCS%window(i_win)%wl_min - hires_pad, &
               MCS%window(i_win)%wl_max + hires_pad)

       else if (MCS%algorithm%solar_type == "oco_hdf") then
          call read_oco_hdf_solar_spectrum(MCS%algorithm%solar_file%chars(), &
               band, &
               solar_spectrum, &
               solar_continuum_from_hdf, &
               MCS%window(i_win)%wl_min - hires_pad, &
               MCS%window(i_win)%wl_max + hires_pad)
       else
          call logger%fatal(fname, "Sorry, solar model type " &
               // MCS%algorithm%solar_type%chars() &
               // " is not known.")
       end if

       ! Need this
       N_solar = size(solar_spectrum, 1)

       ! And we also need the solar spectrum on our user-defined
       ! high-resolution wavelength grid.
       allocate(solar_spectrum_regular(N_hires, 2))
       solar_spectrum_regular(:,1) = hires_grid

       call logger%debug(fname, "Re-gridding solar transmission spectrum")
       call pwl_value_1d( &
            N_solar, &
            solar_spectrum(:,1), solar_spectrum(:,2), &
            N_hires, &
            solar_spectrum_regular(:,1), solar_spectrum_regular(:,2))

       call logger%debug(fname, "Finished re-gridding solar transmission spectrum.")
       ! Note that at this point, the solar spectrum is still normalised


       ! Allocate containers to hold the radiances and noise values - if requested!
       if (MCS%output%save_radiances) then
          allocate(final_radiance(num_pixel, num_fp, num_frames))
          allocate(measured_radiance(num_pixel, num_fp, num_frames))
          allocate(noise_radiance(num_pixel, num_fp, num_frames))
          allocate(wavelength_radiance(num_pixel, num_fp, num_frames))

          ! We fill these with NaNs. This makes plotting a bit easier, since
          ! most plotting routines (at least Matplotlib does it) just skip NaNs,
          ! so you will only see the values that are populated / used.
          final_radiance(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          measured_radiance(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          noise_radiance(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          wavelength_radiance(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       end if

       ! Set up state vector structure here

       ! We can do this once for every window, and simply clear the contents
       ! of the state vectors inside the loop later on. Might not save an awful
       ! lot, but every little helps, I guess.

       ! At the moment, we can do following state vectors:
       !
       ! Lambertian surface albedo, arbitrary (?) polynomial order
       ! Zero-level offset / SIF (constant)
       ! Instrument disperison, arbitrary (?) polynomial order
       ! Surface pressure
       ! Solar shift and stretch
       ! Gas scalar factor defined on level ranges
       ! ILS stretch/squeeze, arbitrary (?) polynomial order
       ! Temperature offset

       ! Parsing the statevector string, that was passed in the window
       ! section and initialize the state vector SV accordingly. This subroutine
       ! needs to access plenty of things in the MCS, so we only pass the
       ! window index, and the routine takes care of arranging the rest.

       ! For the beginning, we start by initialising it with the number of levels
       ! as obtained from the initial_atm

       call parse_and_initialize_SV(i_win, size(initial_atm%p), global_SV)
       call logger%info(fname, "Initialised SV structure")

       ! And now set up the result container with the appropriate sizes for arrays
       call create_result_container(results, num_frames, num_fp, &
            size(global_SV%svap), global_SV%num_gas)

       ! Create the SV names corresponding to the SV indices
       call assign_SV_names_to_result(results, global_SV, i_win)


       call logger%info(fname, "Starting main retrieval loop!")

       ! If we use averaging for solar spectra, we really only have to retrieve
       ! one frame.
       if ((MCS%algorithm%solar_footprint_averaging) .and. &
            (MCS%algorithm%observation_mode%lower() == "space_solar")) then
          frame_start = 1
          frame_skip = 1
          frame_stop = 1
       else
          frame_start = 1
          frame_skip = MCS%window(i_win)%frame_skip
          frame_stop = num_frames
       end if

       ! retr_count keeps track of the number of retrievals processed
       ! so far, and the mean_duration keeps track of the average
       ! processing time.
       retr_count = 0
       total_number_todo = (num_fp * frame_stop / frame_skip) / &
            (MCS%window(i_win)%frame_skip * MCS%window(i_win)%footprint_skip)
       mean_duration = 0.0d0

       ! At the moment, OpenMP is implemented to spread the loop over many threads,
       ! and should be maxed out at the number of available cores on the same machine.
       ! This could be further expanded to use MPI, but at the moment, this seems fast
       ! enough for my humble purposes.
       ! As you can see, it does not require much, whereas MPI would be more effort.

       !$OMP PARALLEL DO SHARED(retr_count, mean_duration) SCHEDULE(static) &
       !$OMP PRIVATE(i_fr, i_fp, cpu_time_start, cpu_time_stop, this_thread, this_converged)
       do i_fr=frame_start, frame_stop, frame_skip
          do i_fp=1, num_fp !, MCS%window(i_win)%footprint_skip

#ifdef _OPENMP
             this_thread = OMP_GET_THREAD_NUM()
             cpu_time_start = omp_get_wtime()
#else
             this_thread = 0
             call cpu_time(cpu_time_start)
#endif

             !write(*,*) land_fraction(i_fp, i_fr)
             !if (land_fraction(i_fp, i_fr) < 100) then
             !   cycle
             !end if

             ! ---------------------------------------------------------------------
             ! Do the retrieval for this particular sounding -----------------------
             this_converged = physical_FM(my_instrument, i_fp, i_fr, i_win, band)
             ! ---------------------------------------------------------------------

#ifdef _OPENMP
             cpu_time_stop = omp_get_wtime()
#else
             call cpu_time(cpu_time_stop)
#endif
             ! Increase the rerival count tracker and compute the average processing
             ! time per retrieval.
             retr_count = retr_count + 1
             mean_duration = mean_duration * (retr_count)/(retr_count+1) + &
                  (cpu_time_stop - cpu_time_start) / (retr_count+1)

             if (mod(retr_count, 1) == 0) then
                write(tmp_str, '(A, G0.1, A, G0.1, A, F6.2, A, F10.5, A, L1, A, G0.1)') &
                     "Frame/FP: ", i_fr, "/", i_fp, " ( ", &
                     dble(100 * dble(retr_count) / dble(total_number_todo)), "%) - ", &
                     (cpu_time_stop - cpu_time_start), ' sec. - Converged: ', &
                     this_converged, ", Thread #: ", this_thread
                call logger%debug(fname, trim(tmp_str))
             end if

          end do
       end do
       !$OMP END PARALLEL DO



       !---------------------------------------------------------------------
       ! HDF OUTPUT
       ! Here, we write out the various arrays into HDF datasets
       !---------------------------------------------------------------------

       ! Create an HDF group for all windows separately
       group_name = "RetrievalResults/physical/" // trim(MCS%window(i_win)%name%chars())
       call h5gcreate_f(output_file_id, trim(group_name), result_gid, hdferr)
       call check_hdf_error(hdferr, fname, "Error. Could not create group: " &
            // trim(group_name))

       ! Set the dimensions of the arrays for saving them into the HDF file
       out_dims2d(1) = num_fp
       out_dims2d(2) = num_frames

       ! Writing out the prior surface pressure, but obviously only if allocated
       if (allocated(met_psurf)) then
          call logger%info(fname, "Writing out: " // trim(group_name) // &
               "/surface_pressure_apriori_" // MCS%general%code_name)
          write(tmp_str, '(A,A,A)') trim(group_name), "/surface_pressure_apriori_", MCS%general%code_name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), met_psurf(:,:), out_dims2d, -9999.99d0)
       end if

       ! Save the prior state vectors
       do i=1, size(global_SV%svsv)
          write(tmp_str, '(A,A,A,A,A)') trim(group_name) , "/", &
               results%sv_names(i)%chars() , "_apriori_", MCS%general%code_name
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_prior(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the retrieved state vectors
       do i=1, size(global_SV%svsv)
          write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", results%sv_names(i)%chars(), &
               "_", MCS%general%code_name
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_retrieved(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the retrieved state vector uncertainties
       do i=1, size(global_SV%svsv)
          write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", &
               results%sv_names(i)%chars() , "_uncertainty_", MCS%general%code_name
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_uncertainty(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the dry air mass (molecules / cm2)
       do i=1, global_SV%num_gas
          if (MCS%window(i_win)%gas_retrieve_scale(global_SV%gas_idx_lookup(i))) then

             ! This is a scale-type jacobian
             do j=1, size(MCS%window(i_win)%gas_retrieve_scale_start(global_SV%gas_idx_lookup(i),:))

                if (global_SV%gas_idx_lookup(i) == j) then
                   write(tmp_str, "(A,A,A,A)") trim(group_name) // "/", &
                        results%SV_names(global_SV%idx_gas(i,1))%chars(), &
                        "_ndry_", MCS%general%code_name
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
             lower_str = MCS%window(i_win)%gases(i)%lower()
             write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/x", &
                  lower_str%chars(), "_", MCS%general%code_name
             call logger%info(fname, "Writing out: " // trim(tmp_str))
             call write_DP_hdf_dataset(output_file_id, &
                  trim(tmp_str), &
                  results%xgas(:,:,i), out_dims2d, -9999.99d0)
          end if
       end do

       ! Save number of iterations
       call logger%info(fname, "Writing out: " // trim(group_name) // "/num_iterations_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/num_iterations_", MCS%general%code_name
       call write_INT_hdf_dataset(output_file_id, &
            trim(tmp_str), results%num_iterations(:,:), out_dims2d, -9999)

       ! Save converged status
       call logger%info(fname, "Writing out: " // trim(group_name) // "/converged_flag_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/converged_flag_", MCS%general%code_name
       call write_INT_hdf_dataset(output_file_id, &
            trim(tmp_str), results%converged(:,:), out_dims2d, -9999)

       ! Retrieved CHI2
       call logger%info(fname, "Writing out: " // trim(group_name) // "/reduced_chi_squared_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/reduced_chi_squared_", MCS%general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%chi2(:,:), out_dims2d, -9999.99d0)

       ! Residual RMS
       call logger%info(fname, "Writing out: " // trim(group_name) // "/residual_rms_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/residual_rms_", MCS%general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%residual_rms(:,:), out_dims2d, -9999.99d0)

       ! Dsigma_sq
       call logger%info(fname, "Writing out: " // trim(group_name) // "/final_dsigma_sq_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/final_dsigma_sq_", MCS%general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%dsigma_sq(:,:), out_dims2d, -9999.99d0)

       ! Signal-to-noise ratio (mean)
       call logger%info(fname, "Writing out: " // trim(group_name) // "/snr_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/snr_", MCS%general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%SNR, out_dims2d)

       !call logger%info(fname, "Writing out: " // trim(group_name) // "/snr_std_" &
       !     // MCS%general%code_name)
       !write(tmp_str, '(A,A,A)') trim(group_name), "/snr_std_", MCS%general%code_name
       !call write_DP_hdf_dataset(output_file_id, &
       !     trim(tmp_str), results%SNR_std, out_dims2d)

       call logger%info(fname, "Writing out: " // trim(group_name) // "/continuum_level_radiance_" &
            // MCS%general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/continuum_level_radiance_", MCS%general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), results%continuum, out_dims2d)

       ! Save the radiances, only on user request (non-default)
       if (MCS%output%save_radiances) then

          out_dims3d = shape(final_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/modelled_radiance_" &
               // MCS%general%code_name)
          write(tmp_str, '(A,A,A)') trim(group_name), "/modelled_radiance_", MCS%general%code_name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), final_radiance, out_dims3d)

          out_dims3d = shape(measured_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/measured_radiance_" &
               // MCS%general%code_name)
          write(tmp_str, '(A,A,A)') trim(group_name), "/measured_radiance_", MCS%general%code_name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), measured_radiance, out_dims3d)

          out_dims3d = shape(noise_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/noise_radiance_" &
               // MCS%general%code_name)
          write(tmp_str, '(A,A,A)') trim(group_name), "/noise_radiance_", MCS%general%code_name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), noise_radiance, out_dims3d)

          out_dims3d = shape(wavelength_radiance)
          call logger%info(fname, "Writing out: " // trim(group_name) // "/wavelength_" &
               // MCS%general%code_name)
          write(tmp_str, '(A,A,A)') trim(group_name), "/wavelength_", MCS%general%code_name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), wavelength_radiance, out_dims3d)

          ! Also deallocate containers holding the radiances
          deallocate(final_radiance)
          deallocate(measured_radiance)
          deallocate(noise_radiance)
          deallocate(wavelength_radiance)

       end if

       !---------------------------------------------------------------------
       ! CLEAN UP PHASE
       ! We need to deallocate / destroy a few arrays here, so
       ! they can be freshly reallocated for the next retrieval microwindow.
       !---------------------------------------------------------------------


       ! Deallocate arrays that are allocated on per-window basis
       deallocate(solar_spectrum_regular, solar_spectrum, hires_grid)
       if (allocated(solar_continuum_from_hdf)) deallocate(solar_continuum_from_hdf)

       ! Clear and deallocate the SV structure to be ready for the next window.
       ! We really shouldn't need to check if this is deallocated already or not,
       ! since parse_and_initialize_SV should have successfully allocated
       ! all fields in SV.
       call clear_SV(global_SV)
       call logger%info(fname, "Clearing up SV structure")

       ! Clear up sounding geometry
       deallocate(SZA, SAA, VZA, VAA)
       deallocate(lon, lat, altitude, relative_velocity, relative_solar_velocity)

       ! If present, deallocate a bad sample list as well as the spike filter list
       if (allocated(bad_sample_list)) deallocate(bad_sample_list)
       if (allocated(spike_list)) deallocate(spike_list)
       if (allocated(stokes_coef)) deallocate(stokes_coef)

       ! Clear and deallocate the result container
       call logger%info(fname, "Clearing up results container.")
       call destroy_result_container(results)

       ! We also de-allocate ABSCO-related fields for safety. The ABSCO read routine
       ! will overwrite the data anyway, but you know .. bad practice.
       do i=1, size(MCS%gas)
          ! Skip over unused
          if (.not. MCS%gas(i)%used) cycle

          if (allocated(MCS%gas(i)%cross_section)) deallocate(MCS%gas(i)%cross_section)
          if (allocated(MCS%gas(i)%wavelength)) deallocate(MCS%gas(i)%wavelength)
          if (allocated(MCS%gas(i)%T)) deallocate(MCS%gas(i)%T)
          if (allocated(MCS%gas(i)%p)) deallocate(MCS%gas(i)%p)
          if (allocated(MCS%gas(i)%H2O)) deallocate(MCS%gas(i)%H2O)

       end do

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
         radiance_tmp_hi_nosif_nozlo(:), &
         radiance_calc_work_hi(:), &
         noise_work(:)

    ! Radiative transfer models
    integer :: RT_model

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
    ! Epoch, which is just the date split into an integer array
    integer, dimension(7) :: epoch

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
    integer :: center_pixel, center_pixel_hi

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
    double precision, allocatable :: gas_tau_dtemp(:,:,:) ! dtau / dvmr (spectral, layer, gas number)
    double precision, allocatable :: gas_tau_dpsurf(:,:,:) ! dtau / dpsurf (spectral, layer, gas_number)
    double precision, allocatable :: gas_tau_pert(:,:,:,:) ! Perturbed gas optical depth (spectral, layer, gas number)
    double precision, allocatable :: rayleigh_tau(:,:) ! Rayleigh optical depth (spectral, layer)
    double precision, allocatable :: rayleigh_depolf(:) ! Rayleigh depolarization factor (spectral)

    ! Perturbed VMR profile and per-iteration-and-per-gas VMR profile for OD calculation (level)
    double precision, allocatable :: vmr_pert(:), this_vmr_profile(:,:)
    ! Rayleigh extinction optical depth (spectral, layer)
    double precision, allocatable :: ray_tau(:,:)
    ! Total column optical depth (spectral)
    double precision, allocatable :: total_tau(:)
    ! Start and end positions in the atmosphere of the gas scalar
    integer :: s_start(global_SV%num_gas), s_stop(global_SV%num_gas)
    ! Is this gas H2O?
    logical :: is_H2O
    ! Number of molecules of dry air per m^2
    double precision, allocatable :: ndry(:), ndry_tmp(:)

    ! Albedo
    ! Prior albedo value estimated from the radiances
    double precision :: albedo_apriori
    ! Per-wavelength albedo for hires and low-res spectra
    double precision, allocatable :: albedo(:)

    ! Surface pressure
    ! The surface pressure per-iteration
    double precision :: this_psurf
    ! Do we need/want to calculate surface pressure Jacobians?
    logical :: do_psurf_jac

    ! SIF
    ! Per-iteration SIF radiance
    double precision :: this_sif_radiance
    ! ZLO - same as SIF essentially
    double precision :: this_zlo_radiance

    ! Temperature
    double precision :: this_temp_offset

    ! Gases
    ! Do we want to calculate gas Jacobians
    logical :: do_gas_jac
    ! Was the calculation of gas ODs successful?
    logical :: success_gas
    ! Was the subcolumn boundary calculation successful?
    logical :: success_scale_levels
    double precision, allocatable :: scale_first_guess(:)

    ! Retrieval quantities
    type(statevector) :: SV
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
    ! Correlation matrix of Shat
    double precision, allocatable :: Shat_corr(:,:)


    ! Temporary matrices and vectors for computation
    double precision, allocatable :: tmp_m1(:,:), tmp_m2(:,:)
    double precision, allocatable :: tmp_v1(:), tmp_v2(:)
    double precision, allocatable :: KtSeK(:,:), gain_matrix(:,:), AK(:,:)
    ! Was the matrix inversion operation successful?
    logical :: success_inv_mat
    ! Was the convolution operation successful?
    logical :: ILS_success
    ! Was the time string conversion successful?
    logical :: success_time_convert
    ! Are the radiances proper and valid?
    logical :: radiance_OK

    ! ILS stuff
    double precision, allocatable :: this_ILS_stretch(:)
    double precision, allocatable :: this_ILS_delta_lambda(:,:)

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
    double precision :: old_chi2, this_chi2, &
         linear_prediction_chi2, chi2_ratio, chi2_rel_change

    ! Iteration-related
    ! Current iteration number (starts at 1), number of divergent steps allowed.
    integer :: iteration, num_divergent_steps
    ! Variable to determine if we keep doing the iteration loop
    logical :: keep_iterating
    ! Has last iteration been a divergent one?
    logical :: divergent_step

    ! XRTM Radiative Transfer model handler
    type(xrtm_type) :: xrtm
    type(string) :: xrtm_options_string(1)
    type(string) :: xrtm_solvers_string(1)
    logical :: xrtm_success

    ! Miscellaneous stuff
    ! String to hold various names etc.
    character(len=999) :: tmp_str
    ! Function name
    character(len=*), parameter :: fname = "physical_FM"
    ! Loop variables
    integer :: i, j, l
    ! File unit for debugging
    integer :: funit


    ! DEBUG STUFF
    double precision :: cpu_start, cpu_stop

    ! Grab a copy of the state vector for local use
    SV = global_SV

    ! Take a local copy of the HDF file ID handlers
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id

    ! Set the used radiative transfer model
    if (MCS%window(i_win)%RT_model%lower() == "beer-lambert") then
       RT_model = RT_BEER_LAMBERT
    else if (MCS%window(i_win)%RT_model%lower() == "xrtm") then
       RT_model = RT_XRTM
    else
       call logger%error(fname, "RT Method: " // MCS%window(i_win)%RT_model%chars() // " unknown.")
       stop 1
    end if

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
       ! READ THE MEASUREMENT

       if ((MCS%algorithm%solar_footprint_averaging) .and. &
            (MCS%algorithm%observation_mode%lower() == "space_solar")) then
          ! If the user chooses to, we can average all frames of a certain footprint
          ! To one single measurement.
          call my_instrument%read_spectra_and_average_by_fp(l1b_file_id, i_fp, band, &
               MCS%general%N_spec(band), radiance_l1b)
       else
          ! Read the L1B spectrum for this one measurement in normal mode!
          call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, band, &
               MCS%general%N_spec(band), radiance_l1b)
       end if
       ! Convert the date-time-string object in the L1B to a date-time-object "date"
       call my_instrument%convert_time_string_to_date(frame_time_strings(i_fr), &
            date, success_time_convert)
    end select

    ! If extraction of the date-time fails, abort and return. Date and time is
    ! needed for the solar doppler calculation.
    if (.not. success_time_convert) then
       call logger%error(fname, "Time string conversion error!")
       return
    end if


    ! Estimate a smart first guess for the gas scale factor, if the user supplied
    ! values for expected delta tau etc.
    if (allocated(MCS%window(i_win)%smart_scale_first_guess_wl_in)) then

       allocate(scale_first_guess(size(MCS%window(i_win)%smart_scale_first_guess_wl_in)))

       call estimate_first_guess_scale_factor(dispersion(:, i_fp, band), &
            radiance_l1b, &
            MCS%window(i_win)%smart_scale_first_guess_wl_in(:), &!expected_wavelengths_in, &
            MCS%window(i_win)%smart_scale_first_guess_wl_out(:), &!expected_wavelengths_out, &
            MCS%window(i_win)%smart_scale_first_guess_delta_tau(:), &!expected_delta_tau, &
            SZA(i_fp, i_fr), VZA(i_fp, i_fr), scale_first_guess)
    else
       ! Otherwise just start with 1.0
       allocate(scale_first_guess(1))
       scale_first_guess(1) = 1.0d0
    end if


    ! Calculate the day of the year into a full fractional value
    doy_dp = dble(date%yearday()) + dble(date%getHour()) / 24.0d0

    ! Epoch is needed by the MS3 solar doppler code
    epoch(1) = date%getYear()
    epoch(2) = date%getMonth()
    epoch(3) = date%getDay()
    epoch(4) = date%getHour()
    epoch(5) = date%getMinute()
    epoch(6) = date%getSecond()


    ! ----------------------------------------------
    ! Print out some debug information for the scene
    ! ----------------------------------------------

    write(tmp_str, "(A, A)") "Date: ", date%isoformat()
    call logger%debug(fname, trim(tmp_str))

    write(tmp_str, "(A, F5.1)") "Day of the year: ", doy_dp
    call logger%debug(fname, trim(tmp_str))

    write(tmp_str, "(A, F6.2, A, F6.2, A, F6.2, A, F6.2)") &
         "SZA: ", SZA(i_fp, i_fr), " - SAA: ", SAA(i_fp, i_fr), &
         " - VZA: ", VZA(i_fp, i_fr), " - VAA: ", VAA(i_fp, i_fr)
    call logger%debug(fname, trim(tmp_str))

    write(tmp_str, "(A, F8.2, A, F8.2, A, ES15.3)") &
         "Lon.: ", lon(i_fp, i_fr), &
         " Lat.: ", lat(i_fp, i_fr), &
         " Alt.: ", altitude(i_fp, i_fr)
    call logger%debug(fname, trim(tmp_str))

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

    ! The "instrument doppler shift" is caused by the relative velocity
    ! between the point on the surface and the spacecraft. Obviously, this
    ! contribution is zero for space-solar geometries.
    if (MCS%algorithm%observation_mode == "downlooking") then
       instrument_doppler = relative_velocity(i_fp, i_fr) / SPEED_OF_LIGHT

       ! THIS IS "BORROWED" FROM THE MS3 CODE
       call solar_doppler_velocity(SZA(i_fp, i_fr), SAA(i_fp, i_fr), &
            epoch, lat(i_fp, i_fr), altitude(i_fp, i_fr), solar_rv, solar_dist)
       solar_doppler = solar_rv / SPEED_OF_LIGHT

    else if (MCS%algorithm%observation_mode == "space_solar") then

       ! For space-solar observation mode, the doppler is obviously different
       ! so we need to change the calculation slightly.
       call solar_distance_and_velocity_v2(epoch, solar_dist, solar_rv)

       instrument_doppler = 0.0d0
       solar_doppler = solar_rv / SPEED_OF_LIGHT
    end if

    ! Set up retrieval quantities:
    N_sv = size(SV%svap)
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
    allocate(Shat_corr(N_sv, N_sv))
    allocate(tmp_m1(N_sv, N_sv), tmp_m2(N_sv, N_sv))
    allocate(tmp_v1(N_sv), tmp_v2(N_sv))
    allocate(KtSeK(N_sv, N_sv))
    allocate(AK(N_sv, N_sv))

    allocate(radiance_calc_work_hi(size(this_solar, 1)))
    allocate(radiance_tmp_hi_nosif_nozlo(size(this_solar, 1)))
    allocate(albedo(size(radiance_calc_work_hi)))


    Sa_inv(:,:) = 0.0d0
    Sa(:,:) = 0.0d0

    ! Separate function to populate the prior covariance - this contains
    ! a good number of hard-coded values, which in future should be
    ! given through the config file.
    call populate_prior_covariance(SV, percentile(radiance_l1b, 98.0d0), &
         i_win, Sa, Sa_inv)

    ! -----------------------------------------------------------------------
    ! Get the albedo prior estimated through the L1B radiances using a simple
    ! Lambertian model.
    select type(my_instrument)
    type is (oco2_instrument)

       allocate(solar_irrad(N_hires))
       if (MCS%algorithm%solar_type == "toon") then
          ! Allocate Solar continuum (irradiance) array. We do this here already,
          ! since it's handy to have it for estimating the albedo

          call calculate_solar_planck_function(6500.0d0, solar_dist, &
               solar_spectrum_regular(:,1), solar_irrad)

       else if (MCS%algorithm%solar_type == "oco_hdf") then
          ! The OCO-HDF-type spectrum needs to be first re-gridded here
          call pwl_value_1d( &
               N_solar, &
               solar_continuum_from_hdf(:,1), solar_continuum_from_hdf(:,2), &
               N_hires, &
               solar_spectrum_regular(:,1), solar_irrad(:))

          solar_irrad(:) = solar_irrad(:) / ((solar_dist / AU_UNIT )** 2)

       else
          call logger%error(fname, "Solar type: " // MCS%algorithm%solar_type%chars() // "is unknown.")
       end if

       ! Otherwise, if we use an OCO/HDF-like solar spectrum, that already
       ! comes with its own irradiance

       albedo_apriori = 1.0d0 * PI * maxval(radiance_l1b) / &
            (1.0d0 * maxval(solar_irrad) * mu0)
       ! Take that into account Stokes coef of instrument for the incoming solar irradiance
       albedo_apriori = albedo_apriori / stokes_coef(1, i_fp, i_fr)
    end select
    ! -----------------------------------------------------------------------

    ! We can now populate the prior state vector

    ! Set slope etc. to zero always (why would we ever want to have a prior slope?)
    if (SV%num_albedo > 0) then
       SV%svap(SV%idx_albedo(1)) = albedo_apriori
       do i=2, SV%num_albedo
          SV%svap(SV%idx_albedo(i)) = 0.0d0
       end do
    end if

    ! Solar shift is set to zero
    if (SV%num_solar_shift == 1) then
       SV%svap(SV%idx_solar_shift(1)) = 0.0d0
    end if

    ! Solar stretch factor is set to one
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

    ! Temperature offset starts with zero
    if (SV%num_temp > 0) then
       SV%svap(SV%idx_temp(1)) = 0.0d0
    end if

    ! Dispersion prior is taken from L1B
    if (SV%num_dispersion > 0) then
       ! Start with the L1b dispersion values as priors
       do i=1, SV%num_dispersion
          SV%svap(SV%idx_dispersion(i)) = dispersion_coefs(i, i_fp, band)
       end do
    end if

    ! ILS stretch - we assume the L1B is unstretched
    if (SV%num_ils_stretch > 0) then
       ! Set the first coefficient to 1.0d0, i.e. no stretch
       SV%svap(SV%idx_ils_stretch(1)) = 1.0d0
       do i=2, SV%num_ils_stretch
          ! And set all other coefficients to zero at first
          SV%svap(SV%idx_ils_stretch(i)) = 0.0d0
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
             SV%svap(SV%idx_gas(i,1)) = mean(scale_first_guess(:))
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

          ! Here we set up quantities that need to be done
          ! for the very first iteration.

          divergent_step = .false.
          ! Initialise Chi2 with an insanely large value
          this_chi2 = 9.9d9
          old_chi2 = 9.9d10
          linear_prediction_chi2 = 9.9d9
          ! For the first iteration, we want to use the prior albedo
          albedo(:) = albedo_apriori
          ! and the first guess state vector is the prior
          SV%svsv = SV%svap
          !
          results%sv_prior(i_fp, i_fr, :) = SV%svap
          last_successful_sv(:) = SV%svap
          old_sv(:) = SV%svap

          if (num_gases > 0) then
             ! An atmosphere is only required if there are gases present in the
             ! microwindow.
             this_psurf = met_psurf(i_fp, i_fr)
             this_atm = initial_atm

             if (this_psurf < 1.0d-10) then
                call logger%error(fname, "MET surface pressure is almost 0.")
                return
             end if

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
                   ! If H2O needs to be retrieved, take it from the MET atmosphere
                   ! specific humidty directly, rather than the H2O column of the
                   ! atmosphere text file.
                   this_atm%gas_vmr(:,i) = this_atm%sh / (1.0d0 - this_atm%sh) * SH_H2O_CONV

                end if
             end do

          else
             ! If we don't do gases, just set this variable to zero, mainly to avoid
             ! potential problems .. nasty segfaults etc.
             num_active_levels = 0
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
          ! new albedo(:) array for high-resolution grid. Reminder: albedo slope is
          ! defined with respect to the center pixel (and higher orders).
          if (SV%num_albedo > 0) then
             albedo = 0.0d0
             do i=1, SV%num_albedo
                albedo = albedo + SV%svsv(SV%idx_albedo(i)) &
                     * ((hires_grid(:) - hires_grid(center_pixel_hi)) ** (dble(i-1)))
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
       endif

       ! NOTE
       ! SIF and ZLO are (right now) exactly the same, i.e. an additive radiance
       ! contribution that is constant w.r.t. wavelength.

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


       ! For XRTM, we need to initialize the model with the current number of
       ! active levels / layers. This probably eats up a few cycles, but nothing
       ! we should be worried about. For the rare cases where the number of layers
       ! changes per iterations, this MUST be within the iteration loop

       if (RT_model == RT_XRTM) then
          xrtm_solvers_string(1) = "EIG_BVP"
          xrtm_options_string(1) = ""

          call setup_XRTM(xrtm, xrtm_options_string, &
               xrtm_solvers_string, &
               xrtm_success)
          if (.not. xrtm_success) return

          call create_XRTM(xrtm, 7, 2, 1, 1, num_active_levels - 1, 1, 50, xrtm_success)
          if (.not. xrtm_success) return
       end if



       ! Heavy bit - calculate the optical properties given an atmosphere with gases
       ! and their VMRs. This branch of the code will only be entered if we have at least
       ! one gas present. Otherwise, gas_tau will stay unallocated.

       if ((num_gases > 0) .and. (MCS%algorithm%observation_mode == "downlooking")) then

          allocate(gas_tau(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_pert(N_hires, num_levels-1, num_gases, SV%num_gas))
          allocate(gas_tau_dpsurf(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_dtemp(N_hires, num_levels-1, num_gases))
          allocate(ray_tau(N_hires, num_levels-1))
          allocate(vmr_pert(num_levels))
          allocate(this_vmr_profile(num_levels, num_gases))
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
             ! Main gas loop. This loops over all present gases in the
             ! microwindow, and depending on whether we want to retrieve
             ! this gas or not, the code determines the partial column
             ! segment, and applies the retrieved scale factor to it.

             ! Copy over this gases' VMR profile
             this_vmr_profile(:,j) = this_atm%gas_vmr(:,j)
             do_gas_jac = .false.

             ! Enter this branch if we have at least one retrieved gas
             if (SV%num_gas > 0) then

                if ((iteration == 1) .or. (SV%num_psurf == 1)) then

                   ! From the user SV input, determine which levels belong to the
                   ! gas subcolumn retrieval.
                   call set_gas_scale_levels(SV, j, i_win, this_atm, this_psurf, &
                        s_start, s_stop, do_gas_jac, success_scale_levels)

                   if (.not. success_scale_levels) then
                      call logger%error(fname, "Error calculating subcolumn boundaries.")
                      return
                   end if
                end if

                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then
                      ! Finally, apply the scaling factor to the corresponding
                      ! sections of the VMR profile.
                      this_vmr_profile(s_start(i):s_stop(i), j) = &
                           this_vmr_profile(s_start(i):s_stop(i), j) &
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

             if (SV%num_temp == 1) then
                this_temp_offset = SV%svsv(SV%idx_temp(1))
             else
                this_temp_offset = 0.0d0
             end if

             if (SV%num_temp == 1) then
                ! First, we calculate gas OD's with a 1K temperature perturbation, but we only need
                ! to do this if we retrieve the T offset.

                call calculate_gas_tau( &
                     .true., & ! We are using pre-gridded spectroscopy!
                     is_H2O, & ! Is this gas H2O?
                     hires_grid, & ! The high-resolution wavelength grid
                     this_vmr_profile(:,j), & ! The gas VMR profile for this gas with index j
                     this_psurf, & ! Surface pressure
                     this_atm%p(:), & ! Atmospheric profile pressures
                     this_atm%T(:) + this_temp_offset + 1.0d0, & ! Atmospheric profile temperature plus 1K perturbation
                     this_atm%sh(:), & ! Atmospheric profile humidity
                     MCS%gas(MCS%window(i_win)%gas_index(j)), & ! MCS%gas object for this given gas
                     MCS%window(i_win)%N_sublayers, & ! Number of sublayers for numeric integration
                     do_psurf_jac, & ! Do we require surface pressure jacobians?
                     gas_tau_dtemp(:,:,j), & ! Output: Gas ODs
                     gas_tau_dpsurf(:,:,j), & ! Output: dTau/dPsurf
                     ndry_tmp, & ! Output: ndry molecules per m2
                     success_gas) ! Output: Was the calculation successful?
             end if

             ! Call the function that calculates the gas optical depths
             call calculate_gas_tau( &
                  .true., & ! We are using pre-gridded spectroscopy!
                  is_H2O, & ! Is this gas H2O?
                  hires_grid, & ! The high-resolution wavelength grid
                  this_vmr_profile(:,j), & ! The gas VMR profile for this gas with index j
                  this_psurf, & ! Surface pressure
                  this_atm%p(:), & ! Atmospheric profile pressures
                  this_atm%T(:) + this_temp_offset, & ! Atmospheric profile temperature
                  this_atm%sh(:), & ! Atmospheric profile humidity
                  MCS%gas(MCS%window(i_win)%gas_index(j)), & ! MCS%gas object for this given gas
                  MCS%window(i_win)%N_sublayers, & ! Number of sublayers for numeric integration
                  do_psurf_jac, & ! Do we require surface pressure jacobians?
                  gas_tau(:,:,j), & ! Output: Gas ODs
                  gas_tau_dpsurf(:,:,j), & ! Output: dTau/dPsurf
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

          ! ----------------------------------------------------------
          ! Rayleigh optical depth calculations are done here as well,
          ! and reside in a matrix for all wavelengths and layers
          ! The depolarization factors depend on wavelength only
          ! ----------------------------------------------------------

          allocate(rayleigh_tau(N_hires, size(this_atm%p) - 1))
          allocate(rayleigh_depolf(N_hires))
          call calculate_rayleigh_tau(hires_grid, this_atm%p, &
               rayleigh_tau, rayleigh_depolf)

          ! ---------------------------------------------------------------
          ! If there are aerosols in the scene, calculate the optical depth
          ! profiles here.
          ! ---------------------------------------------------------------


          ! ---------------------------------------------------------------
          ! Optical depth cleanup
          ! ---------------------------------------------------------------

          ! Set tiny gas OD values to some lower threshold. Some RT solvers
          ! do not like gas OD = 0
          where(gas_tau < 1d-10) gas_tau = 1d-10
          where(rayleigh_tau < 1d-10) rayleigh_tau = 1d-10
          ! Total optical depth is calculated as sum of all gas ODs
          allocate(total_tau(N_hires))
          total_tau(:) = sum(sum(gas_tau, dim=2) + rayleigh_tau, dim=2)



          ! ----------------------------------
          ! END SECTION FOR GASES / ATMOSPHERE
          ! ----------------------------------
       end if


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
       call calculate_dispersion_limits(this_dispersion, i_win, &
            l1b_wl_idx_min, l1b_wl_idx_max)

       ! Number of spectral points in the output resolution
       N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1

       center_pixel = N_spec / 2
       center_pixel_hi = N_spec_hi / 2

       ! Allocate various arrays that depend on N_spec
       allocate(K(N_spec, N_sv))
       allocate(gain_matrix(N_sv, N_spec))
       allocate(radiance_meas_work(N_spec))
       allocate(radiance_calc_work(N_spec))
       allocate(radiance_tmp_work(N_spec))
       allocate(noise_work(N_spec))

       !---------------------------------------------------------------------
       ! RT CALCULATIONS
       !---------------------------------------------------------------------

       ! Calculate the sun-normalized TOA radiances and store them in
       ! 'radiance_calc_work_hi'.
       ! One last check if the radiances (measured) are valid

       select type(my_instrument)
       type is (oco2_instrument)
          call my_instrument%check_radiance_valid(l1b_file_id, radiance_l1b, &
               l1b_wl_idx_min, l1b_wl_idx_max, radiance_OK)
       end select

       if (radiance_OK .eqv. .false.) then
          write(tmp_str, '(A,G0.1,A,G0.1,A)') "Sounding (", i_fr, ",", i_fp, &
               ") has invalid L1B radiances. Skipping."
          call logger%error(fname, trim(tmp_str))
          return
       end if

       select case (RT_model)
       case (RT_BEER_LAMBERT)
          call calculate_BL_radiance(hires_grid, mu0, mu, &
               albedo, total_tau, &
               radiance_calc_work_hi)
       case (RT_XRTM)
          call calculate_XRTM_radiance(xrtm, hires_grid, &
               SZA(i_fp, i_fr), VZA(i_fp, i_fr), SAA(i_fp, i_fr), VAA(i_fp, i_fr), &
               albedo, sum(gas_tau, dim=3), 1, 1, num_active_levels - 1, &
               radiance_calc_work_hi)
       end select

       ! Take a copy of the solar spectrum and re-adjust the solar spectrum wavelength grid
       ! According to both pre-computed solar shift as well as potentially retrieved
       ! solar stretch and shift.
       this_solar(:,1) = this_solar_shift + &
            this_solar_stretch * solar_spectrum_regular(:, 1) / (1.0d0 - solar_doppler)

       if (MCS%algorithm%solar_type == "toon") then
          call calculate_solar_planck_function(6500.0d0, solar_dist, &
               this_solar(:,1), solar_irrad)
          ! And multiply to get the full solar irradiance in physical units
          this_solar(:,2) = solar_spectrum_regular(:, 2) * solar_irrad(:)
       else if (MCS%algorithm%solar_type == "oco_hdf") then

          ! Have to manually re-allocate solar_irrad, which is done in the
          ! calculate_solar_planck_function subroutine.
          if (allocated(solar_irrad)) deallocate(solar_irrad)
          allocate(solar_irrad(size(this_solar, 1)))

          call pwl_value_1d_v2( &
               N_solar, &
               solar_continuum_from_hdf(:,1), solar_continuum_from_hdf(:,2), &
               size(this_solar, 1), &
               this_solar(:,1), solar_irrad(:))

          solar_irrad(:) = solar_irrad(:) / ((solar_dist / AU_UNIT )** 2)
          this_solar(:,2) = solar_spectrum_regular(:,2) * solar_irrad(:)
       end if

       ! If we retrieve either solar shift or stretch (or both), then we
       ! need to know the partial derivative of the solar spectrum w.r.t.
       ! wavelength: dsolar_dlambda
       if ((SV%num_solar_shift == 1) .or. (SV%num_solar_stretch == 1)) then

          allocate(solar_tmp(N_hires))
          call calculate_solar_jacobian(this_solar(:,1), this_solar(:,2), &
               solar_irrad, hires_grid, solar_tmp, dsolar_dlambda)

          this_solar(:,2) = solar_tmp(:)
          deallocate(solar_tmp)

       end if


       ! Multiply with the solar spectrum for physical units and add SIF and ZLO contributions.
       ! For various Jacobian calculations we also need the radiance MINUS the additive
       ! contributions, so we store them in a separate array.
       radiance_calc_work_hi(:) = this_solar(:,2) * radiance_calc_work_hi(:)
       radiance_tmp_hi_nosif_nozlo(:) = radiance_calc_work_hi(:)
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:) + this_sif_radiance + this_zlo_radiance

       !---------------------------------------------------------------------
       ! JACOBIAN CALCULATIONS - FOR HIGH-RES SPECTRA BEFORE CONVOLUTION
       !---------------------------------------------------------------------

       ! Solar jacobians - these are independent of the RT model, since the solar
       ! transmittance is just multiplied onto the TOA spectrum after RT calls

       ! Solar shift jacobian
       if (SV%num_solar_shift == 1) then
          K_hi(:, SV%idx_solar_shift(1)) = -radiance_tmp_hi_nosif_nozlo(:) &
               / this_solar(:,2) * dsolar_dlambda(:)
       end if

       ! Solar stretch jacobian
       if (SV%num_solar_stretch == 1) then
          K_hi(:, SV%idx_solar_stretch(1)) = -radiance_tmp_hi_nosif_nozlo(:) &
               / this_solar(:,2) * dsolar_dlambda(:) * solar_spectrum_regular(:, 1) &
               / (1.0d0 - solar_doppler)
       end if

       ! Surface pressure Jacobian
       if (SV%num_psurf == 1) then
          ! This equation requires the TOA radiance before SIF is added, so if we have
          ! SIF in it, take it out beforehand.
          select case (RT_model)
          case (RT_BEER_LAMBERT)
             call calculate_BL_psurf_jacobian(radiance_tmp_hi_nosif_nozlo, &
                  gas_tau_dpsurf, mu0, mu, K_hi(:, SV%idx_psurf(1)))
          case default
             call logger%error(fname, "Surface pressure Jacobian not implemented " &
                  // "for RT Model: " // MCS%window(i_win)%RT_model%chars())
             stop 1
          end select
       end if

       ! Temperature offset Jacobian
       if (SV%num_temp == 1) then
          select case (RT_model)
          case (RT_BEER_LAMBERT)
             call calculate_BL_temp_jacobian(radiance_tmp_hi_nosif_nozlo, &
                  gas_tau, gas_tau_dtemp, mu0, mu,  K_hi(:, SV%idx_temp(1)))
          case default
             call logger%error(fname, "Temperature offset Jacobian not implemented " &
                  // "for RT Model: " // MCS%window(i_win)%RT_model%chars())
             stop 1
          end select
       end if

       ! Gas jacobians
       if (SV%num_gas > 0) then
          do i=1, SV%num_gas
             ! Jacobian for a scalar-type gas retrieval
             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
                select case (RT_model)
                case (RT_BEER_LAMBERT)
                   call calculate_BL_gas_subcolumn_jacobian(radiance_tmp_hi_nosif_nozlo, &
                        mu0, mu, gas_tau(:, s_start(i):s_stop(i)-1, SV%gas_idx_lookup(i)), &
                        SV%svsv(SV%idx_gas(i,1)), K_hi(:, SV%idx_gas(i,1)))
                case default
                   call logger%error(fname, "Gas sub-column Jacobian not implemented " &
                        // "for RT Model: " // MCS%window(i_win)%RT_model%chars())
                end select
             end if
          end do
       end if

       ! Albedo Jacobians
       if (SV%num_albedo > 0) then
          select case (RT_model)
          case (RT_BEER_LAMBERT)
             do i=1, SV%num_albedo
                call calculate_BL_albedo_jacobian(radiance_tmp_hi_nosif_nozlo, albedo, &
                     hires_grid, center_pixel_hi, i, K_hi(:, SV%idx_albedo(i)))
             end do
          case default
             call logger%error(fname, "Albedo Jacobian not implemented " &
                  // "for RT Model: " // MCS%window(i_win)%RT_model%chars())
             stop 1
          end select
       end if


       ! Stokes coefficients
       ! TODO: apply Stokes coefficients here via instrument parameters from L1b?
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:) * stokes_coef(1, i_fp, i_fr)
       K_hi(:,:) = K_hi(:,:) * stokes_coef(1, i_fp, i_fr)

       ! Grab a copy of the L1b radiances
       radiance_meas_work(:) = radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max)

       ! These are the 'trivial' Jacobians that don't really need their own functions
       ! We add the SIF jacobian AFTER the K matrix is allocated
       if (SV%num_sif > 0) then
          ! Plug in the Jacobians (SIF is easy)
          K(:, SV%idx_sif(1)) = 1.0d0
       end if
       ! Same for ZLO
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
                ! TODO: This threshold value should be user-supplied, as well
                ! as the noise inflation factor.
                ! -127 is a fill value of some sort, so ignore that one
                if ((spike_list(i + l1b_wl_idx_min - 1, i_fp, i_fr) >= 5) .and. &
                     (spike_list(i + l1b_wl_idx_min - 1, i_fp, i_fr) /= -127)) then
                   noise_work(i) = noise_work(i) * 10000.0d0
                end if
             end do
          end if

       end select

       ! If we are doing solar measurements and average spectra,
       ! we also need to re-scale the noise to get adjusted CHI2
       if ((MCS%algorithm%solar_footprint_averaging) .and. &
            (MCS%algorithm%observation_mode%lower() == "space_solar")) then

          noise_work(:) = noise_work(:) / sqrt(dble(MCS%general%N_frame))

       end if

       allocate(Se_inv(N_spec, N_spec))
       ! Inverse noise covariance, we keep it diagonal, as usual
       Se_inv(:,:) = 0.0d0
       do i=1, N_spec
          Se_inv(i,i) = 1.0d0 / (noise_work(i) ** 2)
       end do


       ! Allocate ILS stretch factor (pixel dependent)
       allocate(this_ILS_stretch(N_spec))
       this_ILS_stretch(:) = 0.0d0
       ! Build the ILS stretch polynomial
       if (SV%num_ils_stretch > 0) then
          do i=1, N_spec
             do j=1, SV%num_ils_stretch
                this_ILS_stretch(i) = this_ILS_stretch(i) &
                     + (dble(i - center_pixel) ** (j-1) * SV%svsv(SV%idx_ils_stretch(j)))
             end do
          end do
       end if

       ! Allocate the CURRENTLY USED ILS wavelength array
       allocate(this_ILS_delta_lambda(size(ils_delta_lambda, 1), N_spec))
       this_ILS_delta_lambda(:,:) = ils_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band)

       ! If we retrieve an ILS stretch, then apply the stretch factor to the ILS
       ! delta lambda data here.

       if (SV%num_ils_stretch > 0) then
          do i=1, N_spec
             this_ILS_delta_lambda(:,i) = this_ILS_delta_lambda(:,i) * this_ILS_stretch(i)
          end do
       end if

       ! Convolution with the instrument line shape function(s)
       ! Note: we are only passing the ILS arrays that correspond to the
       ! actual pixel boundaries of the chosen microwindow.

       select type(my_instrument)
       type is (oco2_instrument)

          ! Convolution of the TOA radiances
          call oco_type_convolution(hires_grid, radiance_calc_work_hi, &
               this_ILS_delta_lambda(:,:), &
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

             ! If we have any reason to skip this Jacobian index for convolution,
             ! do it now.
             if (skip_jacobian) cycle

             ! Otherwise just convolve the other Jacobians and save the result in
             ! the low-resolution Jacobian matrix 'K'
             call oco_type_convolution(hires_grid, K_hi(:,i), &
                  this_ILS_delta_lambda(:,:), &
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


       !---------------------------------------------------------------------
       ! JACOBIAN CALCULATIONS - FOR SPECTRA AFTER CONVOLUTION
       !---------------------------------------------------------------------

       ! Calculate disperion Jacobians
       if (SV%num_dispersion > 0) then
          do i=1, SV%num_dispersion
             call calculate_dispersion_jacobian(my_instrument, i, &
                  this_dispersion_coefs, band, i_win, i_fp, N_spec, &
                  this_ILS_delta_lambda, ILS_relative_response, &
                  l1b_wl_idx_min, l1b_wl_idx_max, &
                  instrument_doppler, hires_grid, radiance_calc_work, &
                  radiance_calc_work_hi, K(:, SV%idx_dispersion(i)))
          end do
       end if


       ! ILS Jacobians are produced via finite differencing
       if (SV%num_ILS_stretch > 0) then
          ! Loop over all required ILS stretch orders
          do i=1, SV%num_ils_stretch
             call calculate_ILS_stretch_jacobian(my_instrument, i, SV, &
                  band, i_win, i_fp, N_spec, &
                  ils_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  ILS_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  this_ILS_stretch, this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), &
                  l1b_wl_idx_min, l1b_wl_idx_max, &
                  hires_grid, center_pixel, radiance_calc_work, &
                  radiance_calc_work_hi, K(:, SV%idx_ils_stretch(i)))
          end do
       end if


       !---------------------------------------------------------------------
       ! INVERSE SOLVER
       !---------------------------------------------------------------------

       if (MCS%window(i_win)%inverse_method%lower() == "imap") then
          ! Use iterative maximum a-posteriori solution (linear retrieval)
          ! (IMAP) as done by Christian Frankenberg

          ! K^T Se^-1 K
          KtSeK(:,:) = matmul(matmul(transpose(K), Se_inv), K)

          tmp_m1 = Sa_inv + KtSeK
          call invert_matrix(tmp_m1, tmp_m2, success_inv_mat)
          if (.not. success_inv_mat) then
             call logger%error(fname, "Failed to invert K^T Se K")
             return
          end if

          !
          ! tmp_v2 = matmul(matmul(matmul(tmp_m2, transpose(K)), Se_inv), &
          !          radiance_meas_work - radiance_calc_work + matmul(K, SV%svsv - SV%svap))

          ! This here updates the AP essentially, so the next step launches from the
          ! last iteration's result.
          tmp_v2 = matmul(matmul(matmul(tmp_m2, transpose(K)), Se_inv), &
               radiance_meas_work - radiance_calc_work + matmul(K, SV%svsv - old_sv))

          ! Update state vector
          !SV%svsv = SV%svap + tmp_v2
          SV%svsv = old_sv + tmp_v2
       else if (MCS%window(i_win)%inverse_method%lower() == "lm") then

          ! (1+gamma) * Sa^-1 + (K^T Se^-1 K)
          tmp_m1 = (1.0d0 + lm_gamma) * Sa_inv + KtSeK

          ! K^T Se^-1 K
          KtSeK(:,:) = matmul(matmul(transpose(K), Se_inv), K)

          tmp_m1 = Sa_inv + KtSeK
          call invert_matrix(tmp_m1, tmp_m2, success_inv_mat)
          if (.not. success_inv_mat) then
             call logger%error(fname, "Failed to invert K^T Se K")
             return
          end if

          tmp_v1 = matmul(matmul(transpose(K), Se_inv), radiance_meas_work - radiance_calc_work)
          tmp_v2 = matmul(Sa_inv, SV%svsv - SV%svap)
          SV%svsv = SV%svsv + matmul(tmp_m2, tmp_v1 - tmp_v2)
       else
          call logger%error(fname, "Inverse method: " &
               // MCS%window(i_win)%inverse_method%lower() // ", not known!")
          stop 1
       end if

       ! Statevector clamping - for select elements of the state vector, we
       ! want them to NOT change until a certain iteration is reached.

       ! In the case of retrieving gases - we have to adjust the retrieved state vector
       ! if the retrieval wants to push it below 0. Also, we do not allow the gas scale
       ! factors to drop below 50% of the last iteration's value.
       do i=1, SV%num_gas
          do j=1, size(sv%idx_gas(i,:))
             if (SV%idx_gas(i,j) /= -1) then
                if (SV%svsv(sv%idx_gas(i, j)) < (0.5d0 * old_sv(sv%idx_gas(i, j)))) then
                   SV%svsv(sv%idx_gas(i, j)) = (0.5d0 * old_sv(sv%idx_gas(i, j)))
                end if
                if (SV%svsv(sv%idx_gas(i, j)) < 0.01d0) then
                   SV%svsv(sv%idx_gas(i, j)) = 0.01d0
                end if
             end if
          end do
       end do

       ! We limit the retrieved albedo, but only in downlooking mode
       if (MCS%algorithm%observation_mode == "downlooking") then
          if (SV%num_albedo > 0) then
             if (SV%svsv(SV%idx_albedo(1)) < 0.0d0) then
                SV%svsv(SV%idx_albedo(1)) = 1.0d-4
             end if
          end if
       end if

       ! Calculate Shat_inv
       Shat_inv = KtSeK + Sa_inv

       ! Calculate Shat from Shat_inverse
       call invert_matrix(Shat_inv, Shat, success_inv_mat)

       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert Shat^-1")
          return
       end if

       ! Check delta sigma square for this iteration
       dsigma_sq = dot_product(old_sv - SV%svsv, matmul(Shat_inv, old_sv - SV%svsv))

       ! Calculate the chi2 of this iteration
       this_chi2 = calculate_chi2(radiance_meas_work, radiance_calc_work, &
            noise_work, N_spec - N_sv)

       chi2_rel_change = abs(this_chi2 - old_chi2) / old_chi2
       ! Now we check for convergence!
       if ( &
            (dsigma_sq < dble(N_sv) * dsigma_scale) .or. &
            (iteration > MCS%window(i_win)%max_iterations) .or. &
            (chi2_rel_change < 0.01) .or. &
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

          ! Calculate the XGAS for every retreived gas only, same as above with the
          ! gas OD calculation, we loop through all gases, apply the gas scaling
          ! factor from the state vector, and calculate the pressure weighting function
          ! as well as the XGAS.

          if (SV%num_gas > 0) then

             ! Allocate array for pressure weights
             allocate(pwgts(num_active_levels))

             do j=1, num_gases

                ! Here we need to do the same thing as before when calculating
                ! gas OD's. Take local copy of VMR profile, and re-scale the portions
                ! of the profile which, according to the retrieval, have changed.
                this_vmr_profile(:,j) = this_atm%gas_vmr(:,j)

                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then

                      ! Finally, apply the scaling factor to the corresponding
                      ! sections of the VMR profile.
                      this_vmr_profile(s_start(i):s_stop(i),j) = &
                           this_vmr_profile(s_start(i):s_stop(i),j) &
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
                     this_vmr_profile(1:num_active_levels,j), &
                     pwgts)

                ! Compute XGAS as the sum of pgwts times GAS VMRs.
                results%xgas(i_fp, i_fr, j) =  sum( &
                     pwgts(:) &
                     * this_vmr_profile(1:num_active_levels,j) &
                     )

             end do

             deallocate(pwgts)
          end if

          ! Save the final dSigma-squared value (in case anyone needs it)
          results%dsigma_sq(i_fp, i_fr) = dsigma_sq

          ! Calculate the Gain matrix
          gain_matrix(:,:) = matmul(matmul(Shat(:,:), transpose(K)), Se_inv)

          ! Calculate the averaging kernel
          AK(:,:) = matmul(Shat, KtSeK) !matmul(gain_matrix, K) ! - these should be the same?

          ! Calculate state vector element uncertainties from Shat
          do i=1, N_sv
             SV%sver(i) = sqrt(Shat(i,i))
          end do

          ! Put the SV uncertainty into the result container
          results%sv_uncertainty(i_fp, i_fr, :) = SV%sver(:)

          ! Save retrieved CHI2 (before the last update)
          results%chi2(i_fp, i_fr) = this_chi2

          ! Save the residual RMS
          results%residual_rms(i_fp, i_fr) = sqrt(mean((radiance_meas_work - radiance_calc_work)**2))

          ! Get an SNR (mean and std) estimate
          results%SNR(i_fp, i_fr) = mean(radiance_meas_work / noise_work)
          results%SNR_std(i_fp, i_fr) = std(radiance_meas_work / noise_work)

          ! Save also the continuum level radiance
          results%continuum(i_fp, i_fr) = percentile(radiance_meas_work, 99.0d0)

          ! Save statevector at last iteration
          do i=1, size(SV%svsv)
             results%sv_retrieved(i_fp, i_fr,i) = SV%svsv(i)
          end do

          ! Save the number of iterations, remember we start counting at 1,
          ! so the first update is the 2nd iteration etc.
          results%num_iterations(i_fp, i_fr) = iteration

          if (MCS%output%save_radiances) then
             ! Save the radiances and noises - if required by the user
             final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_calc_work(:)
             measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_meas_work(:)
             noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)
             wavelength_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) =  this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max)
          end if

          ! Write a debug message
          write(tmp_str, '(A,G3.1,A,F6.1,A,G2.1,A,F10.2,A,E10.3,A,F10.2)') "Iteration: ", iteration , &
               ", Chi2: ",  results%chi2(i_fp, i_fr), &
               ", Active Levels: ", num_active_levels, &
               ", Psurf: ", this_psurf, &
               ", LM-Gamma: ", lm_gamma, &
               ", SNR: ", results%SNR(i_fp, i_fr)
          call logger%debug(fname, trim(tmp_str))

       else
          ! Not converged yet - change LM-gamma according to the chi2_ratio.

          ! See Rogers (2000) Section 5.7 - if Chi2 increases as a result of the iteration,
          ! we revert to the last state vector, and increase lm_gamma

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


       !---------------------------------------------------------------------
       ! STEP-THROUGH MODE
       !---------------------------------------------------------------------

       ! If the user requests a step-through, then we print a bunch of debug information about
       ! the retrieval for each iteration and wait for a return by the user.
       ! Also, various variables (SV, spectra, AK, gain matrix etc.) are written out.
       if (MCS%algorithm%step_through) then

          call logger%debug(fname, "---------------------------------")
          write(tmp_str, '(A,G0.1,A,G0.1)') "Frame: ", i_fr, ", Footprint: ", i_fp
          call logger%debug(fname, trim(tmp_str))

          ! Gain matrix and AK are normally just calculated after convergence, so
          ! we need to compute them here per-iteration

          ! Calculate the Gain matrix
          gain_matrix(:,:) = matmul(matmul(Shat(:,:), transpose(K)), Se_inv)

          ! Calculate the averaging kernel
          AK(:,:) = matmul(Shat, KtSeK) !matmul(gain_matrix, K) ! - these should be the same?

          do i=1, N_sv
             do j=1, N_sv
                Shat_corr(i,j) = Shat(i,j) / sqrt(Shat(i,i) * Shat(j,j))
             end do
          end do

          open(file="shat_corr.dat", newunit=funit)
          do i=1, N_sv
             write(funit,*) (Shat_corr(i, j), j=1, N_sv)
          end do
          close(funit)
          call logger%debug(fname, "Written file: shat_corr.dat (statevector x statevector)")

          open(file='sv_names.dat', newunit=funit)
          do i=1, N_sv
             write(funit, *) results%sv_names(i)%chars()
          end do
          close(funit)
          call logger%debug(fname, "Written file: sv_names.dat (statevector)")

          open(file="jacobian.dat", newunit=funit)
          do i=1, N_spec
             write(funit,*) (K(i, j), j=1, N_sv)
          end do
          close(funit)
          call logger%debug(fname, "Written file: jacobian.dat (spectral x statevector)")

          open(file="ak_matrix.dat", newunit=funit)
          do i=1, N_sv
             write(funit,*) (AK(i, j), j=1, N_sv)
          end do
          close(funit)
          call logger%debug(fname, "Written file: ak_matrix.dat (statevector x statevector)")

          open(file="gain_matrix.dat", newunit=funit)
          do i=1, N_spec
             write(funit,*) (gain_matrix(j, i), j=1, N_sv)
          end do
          close(funit)
          call logger%debug(fname, "Written file: gain_matrix.dat (spectral x statevector)")

          open(file="spectra.dat", newunit=funit)
          do i=1, N_spec
             write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), &
                  radiance_calc_work(i), noise_work(i)
          end do
          close(funit)
          call logger%debug(fname, "Written file: spectra.dat (wavelength, measured, modelled, noise)")

          open(file="hires_spectra.dat", newunit=funit)
          do i=1, N_hires
             write(funit, *) radiance_calc_work_hi(i)
          end do
          close(funit)
          call logger%debug(fname, "Written file: hires_spectra.dat (modelled)")

          call logger%debug(fname, "---------------------------------")

          write(tmp_str, '(A,G0.1)') "Iteration: ", iteration
          call logger%debug(fname, trim(tmp_str))


          if ((iteration == 1) .and. allocated(scale_first_guess)) then
             do i=1, size(scale_first_guess)
                write(tmp_str, '(A,G0.1,A,F15.5)') "Scale factor first guess #", &
                     i, ": ", scale_first_guess(i)
                call logger%debug(fname, trim(tmp_str))
             end do
          end if

          write(tmp_str,*) "SV index, SV name, SV prior, SV last iteration, " &
               // "SV current, SV delta, prior cov., posterior cov., diagonal of AK"
          call logger%debug(fname, trim(tmp_str))

          do i=1, N_sv
             write(tmp_str, '(I3.1,A25,ES15.6,ES15.6,ES15.6,ES15.6,ES15.6,ES15.6,ES15.6)') &
                  i, results%sv_names(i)%chars(), SV%svap(i), old_sv(i), SV%svsv(i), &
                  SV%svsv(i) - old_sv(i), sqrt(Sa(i,i)), sqrt(Shat(i,i)), AK(i,i)
             call logger%debug(fname, trim(tmp_str))
          end do

          call logger%debug(fname, "---------------------------------")

          write(tmp_str, '(A,F10.3)') "Chi2 (linear prediction): ", linear_prediction_chi2
          call logger%debug(fname, trim(tmp_str))
          write(tmp_str, '(A,F10.3)') "Chi2 (this iteration):    ", this_chi2
          call logger%debug(fname, trim(tmp_str))
          write(tmp_str, '(A,F10.3)') "Chi2 (last iteration):    ", old_chi2
          call logger%debug(fname, trim(tmp_str))
          write(tmp_str, '(A,F10.3)') "Chi2 (relative change):   ", abs(this_chi2 - old_chi2) / old_chi2
          call logger%debug(fname, trim(tmp_str))

          write(tmp_str, '(A,ES15.6,A,ES15.6)') "Dsigma-squared / convergence threshold: ", &
               dsigma_sq, ' / ', dble(N_sv) * dsigma_scale

          call logger%debug(fname, trim(tmp_str))
          if (MCS%window(i_win)%allow_divergences) then
             write(tmp_str, '(A,F10.3)') "LM-gamma: ", lm_gamma
             call logger%debug(fname, trim(tmp_str))

             write(tmp_str, '(A,F10.3)') "Chi2 ratio (R): ", chi2_ratio
             call logger%debug(fname, trim(tmp_str))
          end if

          call logger%debug(fname, "---------------------------------")
          call logger%debug(fname, "Press ENTER to continue")
          read(*,*)

       end if

       ! These quantities are all allocated within the iteration loop, and
       ! hence need explicit de-allocation.
       deallocate(radiance_meas_work, radiance_calc_work, radiance_tmp_work, &
            noise_work, Se_inv, K, gain_matrix, solar_irrad, &
            this_ILS_stretch, this_ILS_delta_lambda)

       if (allocated(gas_tau)) deallocate(gas_tau)
       if (allocated(gas_tau_dpsurf)) deallocate(gas_tau_dpsurf)
       if (allocated(gas_tau_dtemp)) deallocate(gas_tau_dtemp)
       if (allocated(gas_tau_pert)) deallocate(gas_tau_pert)
       if (allocated(ray_tau)) deallocate(ray_tau)
       if (allocated(total_tau)) deallocate(total_tau)
       if (allocated(vmr_pert)) deallocate(vmr_pert)
       if (allocated(this_vmr_profile)) deallocate(this_vmr_profile)
       if (allocated(ndry)) deallocate(ndry)
       if (allocated(ndry_tmp)) deallocate(ndry_tmp)

       if (.not. divergent_step) then
          ! Make sure we keep the 'old' Chi2 only from a valid iteration
          old_chi2 = this_chi2
       end if

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
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 1.0d0 !* SV%svap(SV%idx_albedo(i))
          else
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 1.0d0 !* SV%svap(SV%idx_albedo(1)) ** dble(i)
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

    if (SV%num_temp > 0) then
       ! Put ZLO prior covariance at the continuum level of the band
       Sa(SV%idx_temp(1), SV%idx_temp(1)) = sqrt(5.0d0)
    end if

    if (SV%num_psurf == 1) Sa(SV%idx_psurf(1), SV%idx_psurf(1)) = 1.0d6

    ! Make the dispersion prior covariance about ten times the perturbation value
    if (SV%num_dispersion > 0) then
       do i=1, SV%num_dispersion
          Sa(SV%idx_dispersion(i), SV%idx_dispersion(i)) = &
               MCS%window(i_win)%dispersion_cov(i)**2
               !(10.0d0 * MCS%window(i_win)%dispersion_pert(i))**2
       end do
    end if

    ! Make the ILS stretch prior covariance about ten times the perturbation value
    if (SV%num_ils_stretch > 0) then
       do i=1, SV%num_ils_stretch
          Sa(SV%idx_ils_stretch(i), SV%idx_ils_stretch(i)) = &
               MCS%window(i_win)%ils_stretch_cov(i)**2
               !(10.0d0 * MCS%window(i_win)%ils_stretch_pert(i))**2
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

  !> @brief Calculates the boundaries of the subcolumns
  !>
  !> @param SV State vector object
  !> @param gas_idx gas index (from MCS%gas)

  subroutine set_gas_scale_levels(SV, gas_idx, i_win, atm, psurf, &
       s_start, s_stop, do_gas_jac, success)

    type(statevector), intent(in) :: SV
    integer, intent(in) :: gas_idx
    integer, intent(in) :: i_win
    type(atmosphere), intent(in) :: atm
    double precision, intent(in) :: psurf
    integer, intent(inout) :: s_start(:)
    integer, intent(inout) :: s_stop(:)
    logical, intent(inout) :: do_gas_jac
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "set_gas_scale_levels"
    character(len=999) :: tmp_str
    integer :: i, l


    success = .false.

    ! We need to 'reverse-lookup' to see which SV index belongs to this
    ! gas to grab the right scaling factor. This is done only on the first
    ! iteration - or if we retrieve surface pressure.
    do i=1, SV%num_gas

       if (MCS%window(i_win)%gas_retrieve_scale(gas_idx)) then

          do_gas_jac = .true.
          if (SV%gas_idx_lookup(i) == gas_idx) then

             ! This bit here figures out which level/layer range a
             ! certain scaling factor corresponds to. They are fractions
             ! of surface pressure, so we use 'searchsorted' to find
             ! where they would belong to. We also make sure it can't
             ! go below or above the first/last level.

             s_start(i) = searchsorted_dp(log(atm%p), &
                  SV%gas_retrieve_scale_start(i) * log(psurf), .true.)
             s_start(i) = max(1, s_start(i))
             s_stop(i) = searchsorted_dp(log(atm%p), &
                  SV%gas_retrieve_scale_stop(i) * log(psurf), .true.) + 1
             s_stop(i) = min(size(atm%p), s_stop(i))

          end if
       end if

    end do

    ! We need to make sure that we are not "doubling up" on a specific
    ! gas VMR level when retrieving scale factors. E.g. 0:0.5 0.5:1.0 will
    ! produce overlapping s_start/s_stop.

    do i=1, SV%num_gas
       do l=1, SV%num_gas
          ! Skip gases if index does not match
          if (SV%gas_idx_lookup(i) /= gas_idx) cycle
          if (SV%gas_idx_lookup(l) /= gas_idx) cycle

          if (s_start(i) == s_stop(l)) then
             s_start(i) = s_start(i) + 1
          end if
       end do
    end do

    ! Last check - we run through all gas statevectors and
    ! check if they are at least 2 apart - meaning you can't
    ! (as of now) retrieve a single gas layer.
    do i=1, SV%num_gas

       ! Skip gases if index does not match
       if (SV%gas_idx_lookup(i) /= gas_idx) cycle

       if (s_stop(i) - s_start(i) == 1) then
          write(tmp_str, '(A,A)') "Scale factor index error for ", &
               results%SV_names(SV%idx_gas(i,1))%chars()
          call logger%error(fname, trim(tmp_str))
          return
       end if
    end do

    success = .true.

  end subroutine set_gas_scale_levels


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
