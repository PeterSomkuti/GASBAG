module physical_model_mod

  use ISO_FORTRAN_ENV

  use logger_mod, only: logger => master_logger
  use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
       read_DP_hdf_dataset, write_INT_hdf_dataset
  use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  use control_mod, only: MCS, MAX_WINDOWS, MAX_GASES, MCS_find_gases
  use instruments, only: generic_instrument
  use oco2_mod
  use solar_model_mod
  use math_utils_mod
  use statevector_mod
  use radiance_mod
  use absco_mod
  use Rayleigh_mod
  use gas_tau_mod
  use spectroscopy_utils_mod

  use mod_datetime

  implicit none

  public physical_retrieval

  private

  ! A simple structure to keep the atmosphere data nice
  ! and tidy.
  type atmosphere
     ! The name(s) of the gas(es)
     type(string), allocatable :: gas_names(:)
     ! Gas mixing ratios
     double precision, allocatable :: gas_vmr(:,:)
     ! To which spectroscopy data does this gas correspond to?
     integer, allocatable :: gas_index(:)
     ! Temperature, pressure, and specific humidity
     double precision, allocatable :: T(:), p(:), sh(:)
  end type atmosphere

  ! This structure contains the result data that will be stored in the
  ! output HDF file.

  type result_container
     type(string), allocatable :: sv_names(:)
     double precision, allocatable :: sv_retrieved(:,:,:), &
          sv_prior(:,:,:), sv_uncertainty(:,:,:)
     double precision, allocatable :: xgas(:,:,:)
     double precision, allocatable, dimension(:,:) :: &
          chi2, residual_rms, dsigma_sq
     integer, allocatable, dimension(:,:) :: num_iterations
     logical, allocatable :: converged(:,:)
  end type result_container


  ! High resolution wavelength grid spacing (hires_spacing), and the
  ! padding required make sure that the ILS does not protrude
  ! outside of the grid.
  double precision :: hires_spacing, hires_pad
  ! High resolution wavelength grid
  double precision, allocatable :: hires_grid(:)
  ! Number of gridpoints on the high resolution grid
  integer :: N_hires
  ! Dispersion/wavelength array
  double precision, allocatable :: dispersion(:,:,:)
  ! Sounding_ids
  integer(8), allocatable :: sounding_ids(:,:)
  ! Sounding time strings
  character(len=25), allocatable :: sounding_time_strings(:,:)
  ! Sounding scene geometry (solar and viewing zenith and azimuth)
  double precision, dimension(:,:), allocatable :: SZA, SAA, VZA, VAA
  ! Sounding scene location (lon, lat, altitude, rel-vel instrument-earth,
  !                          rel-vel instrument-sun)
  double precision, dimension(:,:), allocatable :: lon, lat, altitude, &
       relative_velocity, &
       relative_solar_velocity

  ! ILS data: shape: (wl, pixel, fp, band)
  double precision, dimension(:,:,:,:), allocatable :: ils_delta_lambda, &
       ils_relative_response
  ! ILS mapped onto the high resolution grid for fast
  ! ILS "convolution" calculations
  double precision, dimension(:,:,:,:), allocatable :: ils_hires_delta_lambda, &
       ils_hires_relative_response
  double precision :: ils_hires_min_wl, ils_hires_max_wl, ils_hires_spacing
  double precision, allocatable :: ils_hires_grid(:)
  integer :: num_ils_hires
  double precision :: ils_norm_factor

  ! L1B dispersion coefficients
  double precision, allocatable :: dispersion_coefs(:,:,:)
  ! L1B SNR coefficients needed for noise calculation
  double precision, allocatable :: snr_coefs(:,:,:,:)

  ! For OCO-2: bad sample list. Bad samples should not be counted towards
  ! fitting results (noise inflation)
  integer, allocatable :: bad_sample_list(:,:,:)
  integer, allocatable :: spike_list(:,:,:)

  ! We grab the full MET profiles from the MET data file, for quick access
  double precision, allocatable, dimension(:,:,:) :: met_T_profiles, &
       met_P_levels, &
       met_SH_profiles
  ! MET surface pressure
  double precision, allocatable :: met_psurf(:,:)

  ! The initial atmosphere provided by the config file, this will stay the same
  ! for all retrievals. Within the physical_FM function, this atmosphere is
  ! re-gridded, depending on the surface pressure of the particular scene, but
  ! saved into a new variable
  type(atmosphere) :: initial_atm

  ! The solar spectrum (wavelength, transmission)
  double precision, allocatable :: solar_spectrum(:,:), &
       solar_spectrum_regular(:,:), solar_tmp(:)
  integer :: N_solar

  ! Radiances
  double precision, dimension(:,:,:), allocatable :: final_radiance, &
       measured_radiance, &
       noise_radiance

  !
  logical :: first_band_call
  ! State vector construct!
  type(statevector) :: SV
  ! Result container
  type(result_container) :: results


contains

  subroutine physical_retrieval(my_instrument)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument

    integer(hid_t) :: l1b_file_id, met_file_id, output_file_id, dset_id, result_gid
    logical :: MET_exists, ECMWF_exists, spike_exists
    integer(hsize_t) :: out_dims2d(2), out_dims3d(3)
    integer(hsize_t), dimension(:), allocatable :: dset_dims
    integer :: hdferr
    character(len=999) :: dset_name, tmp_str
    character(len=*), parameter :: fname = "physical_retrieval"

    logical :: gas_found, all_gases_found
    integer :: absco_dims

    integer :: num_frames, num_fp, num_spec, num_band
    integer :: i_fp, i_fr, i_pix, i_win, band
    integer :: i, j, retr_count
    logical :: this_converged

    integer :: funit


    double precision :: cpu_time_start, cpu_time_stop, mean_duration

    !! Open up the MET file
    call h5fopen_f(MCS%input%met_filename%chars(), &
         H5F_ACC_RDONLY_F, MCS%input%met_file_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error opening MET file: " &
         // trim(MCS%input%met_filename%chars()))

    !! Store HDF file handler for more convenient access
    met_file_id = MCS%input%met_file_id
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id


    select type(my_instrument)
    type is (oco2_instrument)

       ! How many frames do we have in this file again?
       call my_instrument%read_num_frames_and_fp(l1b_file_id, num_frames, num_fp)

       !! Get the necesary MET data profiles. OCO-like MET data can have either
       !! /Meteorology or /ECMWF (at least at the time of writing this). So we check
       !! which one exists (priority given to /Meteorology) and take it from there

       call h5lexists_f(met_file_id, "/Meteorology", MET_exists, hdferr)
       if (.not. MET_exists) then
          call h5lexists_f(met_file_id, "/ECMWF", ECMWF_exists, hdferr)
       end if

       if ((.not. MET_exists) .and. (.not. ECMWF_exists)) then
          ! Uh-oh, neither /Meteorology nor /ECMWF are found in the
          ! MET file. Can't really go on.
          call logger%fatal(fname, "Neither /Meteorology nor /ECMWF exist in MET file.")
          stop 1
       end if

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
       num_spec = size(snr_coefs, 2)
       num_band = size(snr_coefs, 4)

       ! Read dispersion coefficients and create dispersion array
       call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
       allocate(dispersion(num_spec, num_fp, num_band))

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

       ! Read in bad sample list
       !call my_instrument%read_bad_sample_list(l1b_file_id, bad_sample_list)

    end select

    ! Main loop over all soundings to perform actual retrieval
    ! footprint, frame, microwindow #, band

    ! Create the HDF group in which all the results go in the end
    call h5gcreate_f(output_file_id, "physical_retrieval_results", &
         result_gid, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not create group: physical_retrieval_results")


    do i_win=1, MAX_WINDOWS

       ! Just skip unused windows
       if (.not. MCS%window(i_win)%used) cycle

       ! At the beginning, we check which gases were defined for this window,
       ! and see if a gas with the corresponding name has been defined.
       call MCS_find_gases(MCS%window, MCS%gas, i_win)

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
                call logger%fatal(fname, "Spectroscopy type: " // MCS%gas(j)%type // " not implemented!")
                stop 1
             end if
          end do

          ! Read the atmosphere file (if present) and populate the initial_atm structure
          ! with the contents of said file. It is cross-referenced against the gases in the
          ! list of gases in this window.

          call logger%debug(fname, "Looking into atmosphere file.")
          call read_atmosphere_file(&
               MCS%window(i_win)%atmosphere_file%chars(), &
               MCS%window(i_win)%gas_index, &
               MCS%window(i_win)%gases, &
               initial_atm)

       end if

       ! To make the convolution operation faster, we will here interpolate the
       ! ILS's to higher resolution for all pixels on the very first occassion,
       ! and then we won't have to do it again for every retrieval and every
       ! detector pixel. This obviously only works for a regular solar grid, and
       ! the new hires ILS will be sitting on a regular grid with the same spacing.

       ! Of course every pixel has a different delta_lambda range, so we are going
       ! to create the new ILS grid in such a way that it fits the one ILS with the
       ! largest range.

       band = MCS%window(i_win)%band

       ils_hires_min_wl = minval(ils_delta_lambda(1,:,:,band))
       ils_hires_max_wl = maxval(ils_delta_lambda(size(ils_delta_lambda, 1),:,:,band))

       ! This is the amount by which we have to pad the hi-resolution grid in order to
       ! allow for ILS protrusion. (and add small percentage to be on the safe side)
       hires_pad = (ils_hires_max_wl - ils_hires_min_wl) * 1.10d0

       ! Grab the desired high-resolution wavelength grid spacing
       hires_spacing = MCS%window(i_win)%wl_spacing

       ! .. and construct the high-resolution grid from the supplied
       ! microwindow range.
       N_hires = ceiling((MCS%window(i_win)%wl_max - MCS%window(i_win)%wl_min + 2*hires_pad) / hires_spacing)

       write(tmp_str, '(A,G0.1)') "Number of hires spectral points: ", N_hires
       call logger%debug(fname, trim(tmp_str))

       allocate(hires_grid(N_hires))
       do i=1, N_hires
          hires_grid(i) = MCS%window(i_win)%wl_min - hires_pad + dble(i-1) * hires_spacing
       end do

       ! Now that we have the appropriate high-resolution spacing, we can finish the ILS business
       num_ils_hires = ceiling((ils_hires_max_wl - ils_hires_min_wl) / hires_spacing)

       write(tmp_str, '(A,G0.1)') "Number of hires spectral points for ILS: ", num_ils_hires
       call logger%debug(fname, trim(tmp_str))

       allocate(ils_hires_grid(num_ils_hires))
       allocate(ils_hires_delta_lambda(num_ils_hires, &
            size(ils_delta_lambda, 2), &
            size(ils_delta_lambda, 3), &
            size(ils_delta_lambda, 4)))
       allocate(ils_hires_relative_response, mold=ils_hires_delta_lambda)

       ! We need to ensure that the hi-res grid and the hi-res ILS grid line
       ! up fully, so ils_hires_min_wl needs to be shifted to the closest multiple
       ! of solar_grid_spacing.
       ils_hires_min_wl = hires_spacing * ceiling(ils_hires_min_wl / hires_spacing)

       ! And construct the ILS hires grid using the solar grid spacing
       do i=1, num_ils_hires
          ils_hires_grid(i) = ils_hires_min_wl + dble(i-1) * hires_spacing
       end do

       call logger%debug(fname, "Re-gridding ILS to hi-res grid.")

       do i_fp=1, num_fp
          do i_pix=1, size(dispersion, 1)

             ils_hires_delta_lambda(:, i_pix, i_fp, band) = ils_hires_grid(:)

             ! Interpolate the ILS onto the high-resolution grid
             call pwl_value_1d( &
                  size(ils_delta_lambda, 1), &
                  ils_delta_lambda(:, i_pix, i_fp, band), &
                  ils_relative_response(:, i_pix, i_fp, band), &
                  num_ils_hires, &
                  ils_hires_grid, ils_hires_relative_response(:, i_pix, i_fp, band) )

             ! Since we have resampled the ILS, it also needs re-normalising, a simple
             ! trapezoidal integration should do.

             !ils_norm_factor = 0.0d0
             !do i=1, size(ils_hires_relative_response, 1)-1
             !   ils_norm_factor = ils_norm_factor + &
             !        (ils_hires_relative_response(i, i_pix, i_fp, band) + &
             !        ils_hires_relative_response(i+1, i_pix, i_fp, band)) * &
             !        (ils_hires_delta_lambda(i+1, i_pix, i_fp, band) - &
             !        ils_hires_delta_lambda(i, i_pix, i_fp, band)) / 2.0d0
             !end do

             !ils_hires_relative_response(:, i_pix, i_fp, band) = &
             !     ils_hires_relative_response(:, i_pix, i_fp, band) / ils_norm_factor

             ! And also divide by the sum here, such that we don't have to
             ! re-do the sum every time during the convolution.

             ils_hires_relative_response(:, i_pix, i_fp, band) = &
                  ils_hires_relative_response(:, i_pix, i_fp, band) / &
                  sum(ils_hires_relative_response(:, i_pix, i_fp, band))

          end do
       end do
       call logger%debug(fname, "Done re-gridding ILS.")


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

       ! To make life easier, we want to keep the solar model on the regular,
       ! evenly-spaced highres-grid in wavelength space.

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

       select type(my_instrument)
       type is (oco2_instrument)

          ! Read in the measurement geometries - these are the same for all three bands?
          ! TODO: this should be taking the footprint data, rather than sounding data..
          call my_instrument%read_sounding_geometry(l1b_file_id, band, SZA, SAA, VZA, VAA)
          ! Read in the measurement location
          call my_instrument%read_sounding_location(l1b_file_id, band, lon, lat, &
               altitude, relative_velocity, &
               relative_solar_velocity)

          ! Read in Spike filter data, if it exists in this file
          call h5lexists_f(l1b_file_id, "/SpikeEOF", spike_exists, hdferr)
          if (spike_exists) then
             call my_instrument%read_spike_filter(l1b_file_id, spike_list, band)
          end if

       end select

       ! Allocate containers to hold the retrieval results
       allocate(final_radiance(size(dispersion, 1), num_fp, num_frames))
       allocate(measured_radiance(size(dispersion, 1), num_fp, num_frames))
       allocate(noise_radiance(size(dispersion, 1), num_fp, num_frames))

       final_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       measured_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       noise_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

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
       call create_result_container(results, num_frames, num_fp, size(SV%svap), SV%num_gas)

       ! Create the SV names corresponding to the SV indices
       call assign_SV_names_to_result(results, SV, i_win)

       ! retr_count keeps track of the number of retrievals processed
       ! so far, and the mean_duration keeps track of the average
       ! processing time.
       retr_count = 0
       mean_duration = 0.0d0

       ! If we want to store the pre-calculated cross sections,
       ! we first need to tell the gas_tau function that this is indeed
       ! the very first call, and thus the cross sections need to be
       ! computed.
       first_band_call = .true.

       do i_fr=1, num_frames
          do i_fp=1, num_fp

             call cpu_time(cpu_time_start)
             ! Do the retrieval for this particular sounding
             this_converged = physical_FM(my_instrument, i_fp, i_fr, i_win, band)
             ! After the very first call, we set first_band_call to false, 
             first_band_call = .false.
             call cpu_time(cpu_time_stop)

             ! Increase the rerival count tracker and compute the average processing
             ! time per retrieval.
             retr_count = retr_count + 1
             mean_duration = mean_duration * (retr_count)/(retr_count+1) + &
                  (cpu_time_stop - cpu_time_start) / (retr_count+1)

             if (mod(retr_count, 50) == 0) then
                write(tmp_str, '(A, G0.1, A, G0.1, A, F10.5, A, L1)') &
                     "Frame/FP: ", i_fr, "/", i_fp, " - ", mean_duration, ' - ', this_converged
                call logger%debug(fname, trim(tmp_str))
             end if

          end do
       end do

       ! Deallocate arrays that are allocated on per-window basis
       deallocate(solar_spectrum_regular, solar_spectrum, hires_grid, ils_hires_grid, &
            ils_hires_delta_lambda, ils_hires_relative_response)


       ! Set the dimensions of the arrays for saving them into the HDF file
       out_dims2d(1) = num_fp
       out_dims2d(2) = num_frames

       ! Save the retrieved state vectors
       do i=1, size(SV%svsv)
          write(tmp_str, '(A,A,A,A)') "/physical_retrieval_results/" &
               // MCS%window(i_win)%name // "_" // results%sv_names(i)
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_retrieved(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save the retrieved state vector uncertainties
       do i=1, size(SV%svsv)
          write(tmp_str, '(A,A,A,A,A)') "/physical_retrieval_results/" &
               // MCS%window(i_win)%name // "_" // results%sv_names(i) // "_uncertainty"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%sv_uncertainty(:,:,i), out_dims2d, -9999.99d0)
       end do

       ! Save XGAS for each gas
       do i=1,MCS%window(i_win)%num_gases
          if (MCS%window(i_win)%gas_retrieved(i)) then
             write(tmp_str, '(A,A,A,A)') "/physical_retrieval_results/" &
                  // MCS%window(i_win)%name // "_X" // MCS%window(i_win)%gases(i)
             call write_DP_hdf_dataset(output_file_id, &
                  trim(tmp_str), &
                  results%xgas(:,:,i), out_dims2d, -9999.99d0)
          end if
       end do

       ! Save number of iterations
       write(tmp_str, '(A,A,A)') "/physical_retrieval_results/" &
            // MCS%window(i_win)%name // "_num_iterations"
       call write_INT_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%num_iterations(:,:), out_dims2d, -9999)

       ! Retrieved CHI2
       write(tmp_str, '(A,A,A)') "/physical_retrieval_results/" &
            // MCS%window(i_win)%name // "_retrieved_chi2"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%chi2(:,:), out_dims2d, -9999.99d0)

       ! Dsigma_sq
       write(tmp_str, '(A,A,A)') "/physical_retrieval_results/" &
            // MCS%window(i_win)%name // "_final_dsigma_sq"
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%dsigma_sq(:,:), out_dims2d, -9999.99d0)

       out_dims3d = shape(final_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/modelled_radiance_" &
            // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            final_radiance, out_dims3d)

       out_dims3d = shape(measured_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/measured_radiance_" &
            // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            measured_radiance, out_dims3d)

       out_dims3d = shape(noise_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/noise_radiance_" &
            // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            noise_radiance, out_dims3d)

       ! Clear and deallocate the SV structure to be ready for the next window.
       ! We really shouldn't need to check if this is deallocated already or not,
       ! since parse_and_initialize_SV should have successfully allocated
       ! all fields in SV.
       call clear_SV(SV)
       call logger%info(fname, "Clearing up SV structure")

       ! Also deallocate containers holding the radiances
       deallocate(final_radiance)
       deallocate(measured_radiance)
       deallocate(noise_radiance)

       !deallocate(bad_sample_list)
       if (allocated(spike_list)) deallocate(spike_list)

       ! Clear and deallocate the result container
       call destroy_result_container(results)

       ! End the loop over windows
    end do

  end subroutine physical_retrieval


  function physical_FM(my_instrument, i_fp, i_fr, i_win, band) result(converged)

    implicit none

    class(generic_instrument), intent(in) :: my_instrument
    integer, intent(in) :: i_fr, i_fp, i_win, band
    logical :: converged

    !!
    integer(hid_t) :: l1b_file_id, output_file_id

    !! Radiances and noise arrays
    double precision, dimension(:), allocatable :: radiance_l1b, &
         radiance_tmp_work, &
         radiance_meas_work, &
         radiance_calc_work, &
         radiance_tmp_work_hi, &
         radiance_calc_work_hi, &
         noise_work

    ! The current atmosphere, dependent on surface pressure. This can change
    type(atmosphere) :: this_atm

    !! Dispersion indices
    double precision :: high_res_wl_min, high_res_wl_max
    integer :: l1b_wl_idx_min, l1b_wl_idx_max, solar_idx_min, solar_idx_max

    !! Sounding location /time stuff
    type(datetime) :: date ! Datetime object for sounding date/time
    double precision :: doy_dp ! Day of year as double precision

    !! Instrument stuff
    double precision :: instrument_doppler
    double precision, dimension(:), allocatable :: this_dispersion, &
         this_dispersion_tmp, &
         this_dispersion_coefs, &
         this_dispersion_coefs_pert

    !! For some SV elements, Jacobians are calculated
    !! separately from the main Jacobian convolution loop
    logical :: skip_jacobian

    !! Solar stuff
    double precision :: solar_dist, solar_rv, earth_rv, solar_doppler
    double precision, dimension(:,:), allocatable :: this_solar
    double precision, allocatable :: solar_irrad(:), dsolar_dlambda(:), solar_low(:)
    double precision :: this_solar_shift, this_solar_stretch

    !! Atmosphere
    integer :: num_gases, num_levels, num_active_levels
    double precision, allocatable :: gas_tau(:,:,:), gas_tau_dvmr(:,:,:), &
         gas_tau_dpsurf(:,:,:), gas_tau_pert(:,:,:,:), gas_tau_dsh(:,:,:), &
         vmr_pert(:), this_vmr_profile(:), ray_tau(:,:), total_tau(:)
    integer :: s_start(SV%num_gas), s_stop(SV%num_gas)
    double precision :: gas_scaling_factor
    logical :: is_H2O
    !! Albedo
    double precision :: albedo_apriori
    double precision, allocatable :: albedo(:)

    !! Surface pressure
    double precision :: psurf, this_psurf
    logical :: do_psurf_jac

    !! SIF?
    double precision :: this_sif_radiance

    !! GASES
    logical :: do_gas_jac, success_gas, precompute_CS

    double precision :: continuum

    !! Retrieval quantities
    double precision :: dsigma_scale
    double precision, allocatable :: K_hi(:,:), K(:,:) ! Jacobian matrices
    ! Noise covariance, prior covariance (and its inverse)
    double precision, allocatable :: Se_inv(:,:), Sa(:,:), Sa_inv(:,:)
    double precision, allocatable :: Shat(:,:), Shat_inv(:,:)
    double precision :: dsigma_sq
    ! Temporary matrices and vectors for computation
    double precision, allocatable :: tmp_m1(:,:), tmp_m2(:,:), tmp_m3(:,:)
    double precision, allocatable :: tmp_v1(:), tmp_v2(:), old_sv(:)
    ! Pressure weighting function
    double precision, allocatable :: pwgts(:)
    ! Levenberg-Marquart gamma
    double precision :: lm_gamma
    double precision :: old_chi2, this_chi2
    ! Iteration-related
    integer :: iteration, max_iterations
    logical :: keep_iterating
    logical :: log_retrieval
    logical :: success_inv_mat


    !! Misc
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "physical_FM"
    integer :: N_spec, N_spec_hi, N_spec_tmp, N_sv
    integer :: i, j, l, cnt_j, cnt_l
    integer :: funit
    double precision :: cpu_start, cpu_end
    logical :: ILS_success

    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id

    log_retrieval = .false.

    if (allocated(MCS%window(i_win)%gases)) then
       num_gases = MCS%window(i_win)%num_gases
    else
       num_gases = 0
    end if

    ! Create a datetime object for the current measurement
    select type(my_instrument)
    type is (oco2_instrument)
       call my_instrument%convert_time_string_to_date(sounding_time_strings(i_fp, i_fr), date)
    end select
    doy_dp = dble(date%yearday()) + dble(date%getHour()) / 24.0d0


    ! Use user-supplied value for convergence critertion
    dsigma_scale = MCS%window(i_win)%dsigma_scale
    if (dsigma_scale < 0.0d0) dsigma_scale = 1.0d0

    !write(*,*) "Isoformat: ", date%isoformat()
    !write(*,*) "Day of the year: ", doy_dp

    !write(*,*) solar_rv * 1000.0d0 + earth_rv, solar_doppler
    !write(*,*) "SZA: ", SZA(i_fp, i_fr)
    !write(*,*) "SAA: ", SAA(i_fp, i_fr)
    !write(*,*) "VZA: ", VZA(i_fp, i_fr)
    !write(*,*) "VAA: ", VAA(i_fp, i_fr)
    !write(*,*) "Lon: ", lon(i_fp, i_fr), " Lat: ", lat(i_fp, i_fr), "Altitude: ", altitude(i_fp, i_fr)

    ! Grab the radiance for this particular sounding
    select type(my_instrument)
    type is (oco2_instrument)
       call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, band, &
            size(snr_coefs, 2), radiance_l1b)
    end select

    ! Dispersion array that contains the wavelenghts per pixel
    allocate(this_dispersion(size(radiance_l1b)))
    allocate(this_dispersion_tmp(size(radiance_l1b)))

    ! The dispersion coefficients use to generate this soundings' dispersion
    ! has the same size as the L1b dispersion coefficient array
    allocate(this_dispersion_coefs(size(dispersion_coefs, 1)))
    allocate(this_dispersion_coefs_pert(size(dispersion_coefs, 1)))

    ! Allocate the micro-window bounded solar and radiance arrays
    allocate(this_solar(N_hires, 2))
    allocate(dsolar_dlambda(N_hires))

    ! If solar doppler-shift is needed, calculate the distance and relative
    ! velocity between point of measurement and the (fixed) sun
    call calculate_solar_distance_and_rv(doy_dp, solar_dist, solar_rv)

    call calculate_rel_velocity_earth_sun(lat(i_fp, i_fr), SZA(i_fp, i_fr), &
         SAA(i_fp, i_fr), altitude(i_fp, i_fr), earth_rv)

    solar_doppler = (earth_rv + solar_rv * 1000.0d0) / SPEED_OF_LIGHT
    instrument_doppler = relative_velocity(i_fp, i_fr) / SPEED_OF_LIGHT


    ! Set up retrieval quantities:
    N_sv = size(SV%svap)
    N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1
    N_spec_hi = size(this_solar, 1)
    lm_gamma = 5.0d0

    ! Output-resolution K is allocated within the loop, as the
    ! number of pixels might change, while the hi-res K stays the same
    allocate(K_hi(N_spec_hi, N_sv))
    K_hi(:,:) = 0.0d0

    allocate(Sa(N_sv, N_sv))
    allocate(Sa_inv(N_sv, N_sv))
    allocate(Shat_inv(N_sv, N_sv))
    allocate(Shat(N_sv, N_sv))
    allocate(tmp_m1(N_sv, N_sv), tmp_m2(N_sv, N_sv))
    allocate(tmp_v1(N_sv), tmp_v2(N_sv))

    ! (Inverse) prior covariance matrix Sa(_inv)
    ! Most SV element Sa's are hardcoded for now, as
    ! it would just blow up the config file to specify them
    ! via the user. 

    Sa_inv(:,:) = 0.0d0
    Sa(:,:) = 0.0d0

    if (SV%num_albedo > 0) then
       do i=1, SV%num_albedo
          if (i==1) then
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 100.0d0
          else
             Sa(SV%idx_albedo(i), SV%idx_albedo(i)) = 100.0d0 ** dble(i)
          end if
       end do
    end if

    if (SV%num_solar_shift == 1) then
       Sa(SV%idx_solar_shift(1), SV%idx_solar_shift(1)) = 1.0d0
    end if

    if (SV%num_solar_stretch == 1) then
       Sa(SV%idx_solar_stretch(1), SV%idx_solar_stretch(1)) = 1.0d0
    end if

    if (SV%num_sif > 0) then
       ! Put SIF prior covariance at the continuum level of the band
       Sa(SV%idx_sif(1), SV%idx_sif(1)) = percentile(radiance_l1b, 98.0d0)**2
    end if

    if (SV%num_psurf == 1) Sa(SV%idx_psurf(1), SV%idx_psurf(1)) = 1000.0d0 ** 2

    if (SV%num_dispersion > 0) then
       do i=1, SV%num_dispersion
          Sa(SV%idx_dispersion(i), SV%idx_dispersion(i)) = MCS%window(i_win)%dispersion_pert(i)**2
       end do
    end if

    do i=1, SV%num_gas
       if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
          Sa(SV%idx_gas(i,1), SV%idx_gas(i,1)) = (SV%gas_retrieve_scale_cov(i)) ** 2
       end if
    end do

    ! At the moment, Sa will be diagonal, so inverting is trivial
    do i=1, size(Sa, 1)
       Sa_inv(i,i) = 1.0d0 / Sa(i,i)
    end do

    !! If at some point we want to have non-diagonal elements in
    !! Sa, just uncomment the lines 
    !call invert_matrix(Sa, Sa_inv, success_inv_mat)
    !if (.not. success_inv_mat) then
    !   call logger%error(fname, "Failed to invert Sa")
    !   return
    !end if


    ! Allocate Solar continuum (irradiance) array. We do this here already,
    ! since it's handy to have it for estimating the albedo.
    allocate(solar_irrad(N_hires))

    ! Get the albedo prior

    select type(my_instrument)
    type is (oco2_instrument)

       ! OCO-2 has Stokes coefficient 0.5 for intensity, so we need to
       ! take that into account for the incoming solar irradiance
       call calculate_solar_planck_function(6500.0d0, solar_dist * 1000.0d0, &
            solar_spectrum_regular(:,1), solar_irrad)

       albedo_apriori = 1.0d0 * PI * maxval(radiance_l1b) / &
            (1.0d0 * maxval(solar_irrad) * cos(DEG2RAD * SZA(i_fp, i_fr)))

    end select

    ! Clearly, we could have an albedo that falls outside of (0,1). That does
    ! not mean that it's all screwed up, especially if the true surface is far
    ! away from the Lambertian assumption. Nevertheless, just spit out a quick
    ! warning. (or not)

    if (albedo_apriori > 1) then
       write(tmp_str, '(A, F8.5)') "Albedo too large: ", albedo_apriori
       !call logger%warning(fname, trim(tmp_str))
    end if


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

    ! Dispersion
    if (SV%num_dispersion > 0) then
       ! Start with the L1b dispersion values as priors
       do i=1, SV%num_dispersion
          SV%svap(SV%idx_dispersion(i)) = dispersion_coefs(i,i_fp,band)
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

    ! Retrival iteration loop
    iteration = 1
    keep_iterating = .true.
    allocate(old_sv(size(SV%svsv)))


    ! Calculate the high-resolution modelled spectrum
    ! High-res spectral resolution determined by solar model
    allocate(radiance_calc_work_hi(size(this_solar, 1)))
    allocate(radiance_tmp_work_hi(size(this_solar, 1)))
    allocate(albedo(size(radiance_calc_work_hi)))

    ! Regardless of whether they are retrieved or not, the solar shift and stretch
    ! are set to 'default' values here. If we retrieve them, this value is then just
    ! updated from the state vector.
    this_solar_shift = 0.0d0
    this_solar_stretch = 1.0d0

    converged = .false.
    ! Main iteration loop for the retrieval process.
    do while (keep_iterating)

       if (iteration == 1) then
          !
          this_chi2 = 9.9d9
          ! For the first iteration, we want to use the prior albedo
          albedo(:) = albedo_apriori
          ! and the state vector is the prior (maybe we want to change that later..)
          SV%svsv = SV%svap

          if (num_gases > 0) then

             ! An atmosphere is only required if there are gases present in the
             ! microwindow.
             this_psurf = met_psurf(i_fp, i_fr)
             this_atm = initial_atm
             num_levels = size(this_atm%p)

             ! And get the T and SH profiles onto our new atmosphere grid. We are
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
                call logger%error(fname, "Psurf > p(BOA)!")
                return
             end if

             do i=1, num_gases
                if (MCS%window(i_win)%gases(i) == "H2O") then
                   ! If H2O is needed to be retrieved, take it from the MET atmosphere
                   ! specific humidty directly, rather than the H2O column of the
                   ! atmosphere text file.
                   this_atm%gas_vmr(:,i) = this_atm%sh / (1.0d0 - this_atm%sh) * SH_H2O_CONV
                end if
             end do

          end if
       else

          ! If this is not the first iteration, we grab values from the current
          ! state vector.

          ! Albedo coefficients grabbed from the SV and used to construct
          ! new albedo(:) array for high-resolution grid.
          if (SV%num_albedo > 0) then
             albedo = 0.0d0
             do i=1, SV%num_albedo
                albedo = albedo + SV%svsv(SV%idx_albedo(i)) * ((hires_grid(:) - hires_grid(1)) ** (dble(i-1)))
             end do
          endif

          ! If solar parameters are retrieved, update the solar shift and stretch from the
          ! state vector.
          if (SV%num_solar_shift == 1) this_solar_shift = SV%svsv(SV%idx_solar_shift(1))
          if (SV%num_solar_stretch == 1) this_solar_stretch = SV%svsv(SV%idx_solar_stretch(1))

          ! If we are retrieving surface pressure, we possibly need
          ! to re-grid the atmosphere if the surface pressure jumps to the
          ! next layer.
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
                call logger%error(fname, "Psurf > p(BOA)!")
                return
             end if

          end if

          do i=1, num_gases
             ! Replace SH if H2O is used in atmosphere
             if (MCS%window(i_win)%gases(i) == "H2O") then
                !   this_atm%sh = this_atm%gas_vmr(:,i) / (SH_H2O_CONV + this_atm%gas_vmr(:,i))
             end if
          end do


       endif


       ! Save the old state vector (iteration - 1'th state vector)
       old_sv = SV%svsv

       ! SIF is a radiance, and we want to keep the value for this given
       ! iteration handy for calculations.
       if (SV%num_sif > 0) then
          this_sif_radiance = SV%svsv(SV%idx_sif(1))
       else
          this_sif_radiance = 0.0d0
       end if


       ! Heavy bit - calculate the optical properties given an atmosphere with gases
       ! and their VMRs. This branch of the code will only be entered if we have at least
       ! one gas present. Otherwise, gas_tau will stay unallocated.

       if (num_gases > 0) then

          allocate(gas_tau(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_dpsurf(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_dvmr(N_hires, num_levels, num_gases))
          allocate(gas_tau_dsh(N_hires, num_levels-1, num_gases))
          allocate(gas_tau_pert(N_hires, num_levels-1, num_levels-1, num_gases))
          allocate(ray_tau(N_hires, num_levels-1))
          allocate(vmr_pert(num_levels))
          allocate(this_vmr_profile(num_levels))

          ! If we retrieve surface pressure, grab it from the state vector,
          ! otherwise, grab it from the MET data.

          if (SV%num_psurf == 1) then
             psurf = SV%svsv(SV%idx_psurf(1))
             do_psurf_jac = .true.
          else
             psurf = met_psurf(i_fp, i_fr)
             do_psurf_jac = .false.
          end if

          if (psurf < 0) then
             call logger%error(fname, "Psurf negative!")
             return
          end if

          call calculate_Rayleigh_tau(hires_grid, this_atm%p, ray_tau, first_band_call)

          do j=1, num_gases

             ! Copy over this gases' VMR profile
             this_vmr_profile(:) = this_atm%gas_vmr(:,j)
             do_gas_jac = .false.

             ! Enter this branch if we have at least one retrieved gas
             if (SV%num_gas > 0) then

                ! We need to 'reverse-lookup' to see which SV index belongs to this
                ! gas to grab the right scaling factor.
                do i=1, SV%num_gas

                   if (MCS%window(i_win)%gas_retrieve_scale(j)) then

                      do_gas_jac = .true.
                      if (SV%gas_idx_lookup(i) == j) then

                         ! This bit here figures out which level/layer range a
                         ! certain scaling factor corresponds to. They are fractions
                         ! of surface pressure, so we use 'searchsorted' to find
                         ! where they would belong to. We also make sure it can't
                         ! go below or above the first/last level.

                         s_start(i) = searchsorted_dp(this_atm%p, &
                              SV%gas_retrieve_scale_start(i) * psurf, .false.)
                         s_start(i) = max(1, s_start(i))
                         s_stop(i) = searchsorted_dp(this_atm%p, &
                              SV%gas_retrieve_scale_stop(i) * psurf, .false.)
                         s_stop(i) = min(num_active_levels, s_stop(i))

                      end if
                   end if

                end do

                ! We need to make sure that we are not "doubling up" on a specific
                ! gas VMR level when retrieving scale factors. E.g. 0:0.5 0.5:1.0 will
                ! produce overlapping s_start/s_stop.

                do i=1, SV%num_gas
                   do l=1, SV%num_gas
                      if (SV%gas_idx_lookup(i) /= j) cycle
                      if (SV%gas_idx_lookup(l) /= j) cycle

                      if (s_start(i) == s_stop(l)) then
                         s_start(i) = s_start(i) + 1
                      end if

                   end do
                end do

                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then

                      ! Finally, apply the scaling factor to the corresponding
                      ! sections of the VMR profile.
                      this_vmr_profile(s_start(i):s_stop(i)) = this_vmr_profile(s_start(i):s_stop(i)) &
                           * SV%svsv(SV%idx_gas(i,1))
                   end if
                end do

             end if

             if (MCS%window(i_win)%gases(j) == "H2O") then
                is_H2O = .true.
             else
                is_H2O = .false.
             end if

             call calculate_gas_tau( &
                  .true., & ! We are using pre-gridded spectroscopy!
                  is_H2O, & ! Is this gas H2O?
                  hires_grid, &
                  this_vmr_profile, &
                  psurf, &
                  this_atm%p(:), &
                  this_atm%T(:), &
                  this_atm%sh(:), &
                  MCS%gas(MCS%window(i_win)%gas_index(j)), &
                  MCS%window(i_win)%N_sublayers, &
                  do_psurf_jac, &
                  do_gas_jac, &
                  gas_tau(:,:,j), &
                  gas_tau_dpsurf(:,:,j), &
                  gas_tau_dvmr(:,:,j), &
                  gas_tau_dsh(:,:,j), &
                  .true., &
                  first_band_call, &
                  num_gases, & ! Total number of gases (for precomputed CS)
                  j, & ! Which gas? (for precomputed CS)
                  success_gas)

             if (.not. success_gas) then
                call logger%error(fname, "Error calculating gas optical depths.")
                return
             end if
          end do

          ! Total optical depth is calculated as sum of all gas ODs
          ! as well as the Rayleigh extinction OD, and summing over
          ! all layers as well.
          allocate(total_tau(N_hires))
          total_tau(:) = 0.0d0
          total_tau(:) = sum(sum(gas_tau, dim=2), dim=2) + sum(ray_tau, dim=2)

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

       ! And multiply to get the solar irradiance in physical units
       this_solar(:,2) = solar_spectrum_regular(:, 2) * solar_irrad(:)

       if ((SV%num_solar_shift == 1) .or. (SV%num_solar_stretch == 1)) then

          ! Sample doppler-shifted and scaled solar spectrum at hires grid
          ! TODO: this is quite slow!! I think we need to write a fast
          ! interpolation routine that works on sorted arrays.
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


       ! Multiply with the solar spectrum for physical units and add SIF contributions
       radiance_calc_work_hi(:) = this_solar(:,2) * radiance_calc_work_hi(:) &
            + this_sif_radiance


       !!!!!!!!!!!!!!!!!!!!!!!!!!!! JACOBIAN CALCULATIONS
       ! This probably should go into the radiance module?
       if (SV%num_albedo > 0) then
          do i=1, SV%num_albedo
             K_hi(:, SV%idx_albedo(i)) = (radiance_calc_work_hi(:) - this_sif_radiance) / albedo * &
                  ((hires_grid(:) - hires_grid(1)) ** (dble(i-1)))
          end do
       end if

       ! Surface pressure Jacobian
       if (SV%num_psurf == 1) then
          ! This equation requires the TOA radiance before SIF is added, so if we have
          ! SIF in it, take it out beforehand.
          K_hi(:, SV%idx_psurf(1)) = (radiance_calc_work_hi(:) - this_sif_radiance) &
               * (1.0d0 / cos(DEG2RAD * SZA(i_fp, i_fr)) + 1.0d0 / cos(DEG2RAD * VZA(i_fp, i_fr))) &
               * (sum(sum(gas_tau_dpsurf, dim=2), dim=2))
       end if

       ! Solar shift jacobian
       if (SV%num_solar_shift == 1) then
          K_hi(:, SV%idx_solar_shift(1)) = -(radiance_calc_work_hi(:) - this_sif_radiance) &
               / this_solar(:,2) * dsolar_dlambda(:)
       end if

       ! Solar stretch jacobian
       if (SV%num_solar_stretch == 1) then
          K_hi(:, SV%idx_solar_stretch(1)) = -(radiance_calc_work_hi(:) - this_sif_radiance) &
               / this_solar(:,2) * dsolar_dlambda(:) * solar_spectrum_regular(:, 1) / (1.0d0 - solar_doppler)
       end if


       ! Gas jacobians
       if (SV%num_gas > 0) then
          do i=1, SV%num_gas
             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
                ! This is a scale-type jacobian

                do j=1, size(MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(i),:))

                   ! Loop through all potential profile "sections", but skip the unused ones
                   if (MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(i), j) == -1) cycle

                   K_hi(:, SV%idx_gas(i,1)) = -(radiance_calc_work_hi(:) - this_sif_radiance) &
                        * (1.0d0 / cos(DEG2RAD * SZA(i_fp, i_fr)) + 1.0d0 / cos(DEG2RAD * VZA(i_fp, i_fr))) &
                        * sum(gas_tau(:, s_start(i):s_stop(i)-1, SV%gas_idx_lookup(i)), dim=2) / SV%svsv(SV%idx_gas(i,1))

                   ! If this is a H2O Jacobian, we need to add the derivative
                   ! dtau / dsh * dsh / dh2o
                   !if (MCS%window(i_win)%gases(sv%gas_idx_lookup(i)) == "H2O") then

                      !   K_hi(:, SV%idx_gas(i,1)) = K_hi(:, SV%idx_gas(i,1)) &
                      !        + sum(gas_tau_dsh(:,:,SV%gas_idx_lookup(i)) &
                      !        * (SH_H2O_CONV / (SH_H2O_CONV + gas_tau(:,:,SV%gas_idx_lookup(i)))) ** 2, &
                      !        dim=2)
                   !end if

                end do
             end if

          end do
       end if


       ! Stokes coefficients
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:)
       K_hi(:,:) = K_hi(:,:)

       ! Grab the dispersion coefficients from the L1B .. 
       this_dispersion_coefs(:) = dispersion_coefs(:, i_fp, band)
       ! .. and if required, replace L1b dispersion coefficients by state vector
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

       ! Stretch the dispersion
       this_dispersion = this_dispersion / (1.0d0 - instrument_doppler)

       ! Here we grab the index limits for the radiances for
       ! the choice of our microwindow and the given dispersion relation
       call calculate_dispersion_limits(this_dispersion, i_win, l1b_wl_idx_min, l1b_wl_idx_max)

       ! Number of spectral points in the output resolution
       N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1

       ! Allocate various arrays that depend on N_spec
       allocate(K(N_spec, N_sv))
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


       ! Now calculate the noise-equivalent radiances
       select type(my_instrument)
       type is (oco2_instrument)
          call my_instrument%calculate_noise( &
               snr_coefs, radiance_meas_work, &
               noise_work, i_fp, band, &
               l1b_wl_idx_min, l1b_wl_idx_max)

          ! Pixels flagged with a spike need noise inflation, so
          ! that they are not really considered in the fit. This should
          ! save otherwise good spectra with just a few distorted
          ! radiance values.
          if (allocated(spike_list)) then
             do i=1, N_spec
                if (spike_list(i + l1b_wl_idx_min - 1, i_fp, i_fr) >= 5) then
                   noise_work(i) = noise_work(i) * 10000.0d0
                end if
             end do
          end if

       end select

       allocate(Se_inv(N_spec, N_spec))
       ! Inverse noice covariance, we keep it diagonal, as usual
       Se_inv(:,:) = 0.0d0
       do i=1, N_spec
          Se_inv(i,i) = 1 / (noise_work(i) ** 2)
       end do

       ! Allocation of the Jacobian matrix at the resolution of the
       ! instrument.


       ! Convolution with the instrument line shape function(s)
       ! Note: we are only passing the ILS arrays that correspond to the
       ! actual pixel boundaries of the chosen microwindow.

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


       select type(my_instrument)
       type is (oco2_instrument)

          ! Convolution of the TOA radiances

          call oco_type_convolution(hires_grid, radiance_calc_work_hi, &
               ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
               ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
               this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), radiance_calc_work, &
               ILS_success)

          if (.not. ILS_success) then
             call logger%error(fname, "ILS convolution error.")
             return
          end if


          do i=1, N_sv

             skip_jacobian = .false.

             do j=1, SV%num_dispersion
                ! This is a dispersion jacobian! Maybe there's a way of doing
                ! this analytically, but for now we just perform finite
                ! differencing in a separate loop. So skip this Jacobian.
                if (i == SV%idx_dispersion(j)) then
                   skip_jacobian = .true.
                end if
             end do

             ! SIF Jacobian does not need convolution, since it's just 1.0
             if (i == SV%idx_sif(1)) skip_jacobian = .true.

             ! If we have any reason to skip this Jacobian index for convolution,
             ! do it.
             if (skip_jacobian) cycle

             ! Otherwise just convolve the other Jacobians and save the result in
             ! the low-resolution Jacobian matrix 'K'

             call oco_type_convolution(hires_grid, K_hi(:,i), &
                  ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), K(:,i), &
                  ILS_success)

             if (.not. ILS_success) then
                call logger%error(fname, "ILS convolution error.")
                return
             end if

          end do
       end select

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

                call my_instrument%calculate_dispersion(this_dispersion_coefs_pert, &
                     this_dispersion_tmp(:), band, i_fp)

                this_dispersion_tmp = this_dispersion_tmp / (1.0d0 - instrument_doppler)

                ! Convolve the perturbed TOA radiance
                call oco_type_convolution(hires_grid, radiance_calc_work_hi, &
                     ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                     ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
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
       !tmp_m1 = (1.0d0 + lm_gamma) * Sa_inv + matmul(matmul(transpose(K), Se_inv), K)
       tmp_m1 = matmul(matmul(transpose(K), Se_inv), K)

       call invert_matrix(tmp_m1, tmp_m2, success_inv_mat)
       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert K^T Se K")
          !return
       end if

       tmp_v1 = matmul(matmul(transpose(K), Se_inv), radiance_meas_work - radiance_calc_work)
       !tmp_v2 = matmul(Sa_inv, SV%svsv - SV%svap)

       ! Update state vector
       SV%svsv = SV%svsv + matmul(tmp_m2, tmp_v1)! - tmp_v2)

       ! Calculate Shat_inv
       Shat_inv = matmul(matmul(transpose(K), Se_inv), K) !+ Sa_inv
       call invert_matrix(Shat_inv, Shat, success_inv_mat)

       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert Shat^-1")
          !return
       end if

       ! Check delta sigma square for this iteration
       dsigma_sq = dot_product(old_sv - SV%svsv, matmul(Shat_inv, old_sv - SV%svsv))

       ! In the case of retrieving gas - we have to adjust the retrieved state vector
       ! if the retrieval wants to push it below 0.

       do i=1, SV%num_gas
          do j=1, size(sv%idx_gas(i,:))
             if (SV%idx_gas(i,j) /= -1) then
                if (SV%svsv(sv%idx_gas(i, j)) < 0.0d0) then
                   SV%svsv(sv%idx_gas(i, j)) = 1.0d-10
                end if
             end if
          end do
       end do

!!$       open(file="jacobian.dat", newunit=funit)
!!$       do i=1, N_spec
!!$          write(funit,*) (K(i, j), j=1, N_sv)
!!$       end do
!!$       close(funit)
!!$
!!$       if (allocated(solar_low)) deallocate(solar_low)
!!$       allocate(solar_low(N_spec))
!!$
!!$       call oco_type_convolution(hires_grid, this_solar(:,2), &
!!$            ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
!!$            ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
!!$            this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), solar_low(:), &
!!$            ILS_success)

!!$       open(file="l1b_spec.dat", newunit=funit)
!!$       do i=1, N_spec
!!$          write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), radiance_calc_work(i), &
!!$               noise_work(i)!, solar_low(i)
!!$
!!$       end do
!!$       close(funit)
!!$
!!$       deallocate(solar_low)


!!$       do i=1, num_active_levels
!!$          write(*,*) this_atm%p(i), (this_atm%gas_vmr(i,j), j=1, size(this_atm%gas_vmr, 2))
!!$       end do
!!$
!!$       write(*,*) "old, current and delta state vector, and errors"
!!$       write(*,*) "Iteration: ", iteration
!!$       do i=1, N_sv
!!$          write(*, '(I3.1,A40,ES15.6,ES15.6,ES15.6,ES15.6)') &
!!$               i, results%sv_names(i)%chars(), old_sv(i), SV%svsv(i), &
!!$               SV%svsv(i) - old_sv(i), sqrt(Shat(i,i))
!!$       end do
!!$       write(*,*) "Chi2:    ", SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)
!!$       write(*,*) "Dsigma2: ", dsigma_sq, '/', dble(N_sv) * dsigma_scale

       old_chi2 = this_chi2
       this_chi2 = SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)

       if ((dsigma_sq < dble(N_sv) * dsigma_scale) .or. &
            (iteration > MCS%window(i_win)%max_iterations)) then

       !if ((abs(sqrt(old_chi2)-sqrt(this_chi2)) <= 0.01) .or. &
       !     (iteration > MCS%window(i_win)%max_iterations)) then

          ! Stop iterating - we've either coverged to exeeded the max. number of
          ! iterations.
          keep_iterating = .false.
          if (iteration <= MCS%window(i_win)%max_iterations) then
             converged = .true.
          else
             converged = .false.
          end if

          allocate(pwgts(num_active_levels))

          ! Calculate the XGAS for every retreived gas only, same as above with the
          ! gas OD calculation, we loop through all gases, apply the gas scaling
          ! factor from the state vector, and calculate the pressure weighting function
          ! as well as the XGAS.

          do j=1, num_gases
             if (SV%num_gas > 0) then

                gas_scaling_factor = 1.0d0

                ! We need to 'reverse-lookup' to see which SV index belongs to this
                ! gas to grab the right scaling factor.
                do i=1, SV%num_gas

                   if (SV%gas_idx_lookup(i) == j) then
                      if (MCS%window(i_win)%gas_retrieve_scale(j)) then
                         gas_scaling_factor = SV%svsv(SV%idx_gas(j,1))
                      end if
                   end if

                end do

                call pressure_weighting_function(this_atm%p(1:num_active_levels), &
                     psurf, this_atm%gas_vmr(1:num_active_levels, j) * gas_scaling_factor, pwgts)

                results%xgas(i_fp, i_fr, j) =  sum(pwgts(:) * this_atm%gas_vmr(1:num_active_levels, j) * gas_scaling_factor)

!!$                do l=1, num_active_levels
!!$                   write(*,*) l, this_atm%sh(l), this_atm%gas_vmr(l, j) * gas_scaling_factor, pwgts(l), sum(pwgts(:))
!!$                end do
!!$                write(*,*) "xgas:", sum(pwgts(:) * this_atm%gas_vmr(1:num_active_levels, j) * gas_scaling_factor)

             end if
          end do


          results%dsigma_sq(i_fp, i_fr) = dsigma_sq

          ! Calculate state vector element uncertainties from Shat
          do i=1, N_sv
             SV%sver(i) = sqrt(Shat(i,i))
          end do

          results%sv_uncertainty(i_fp, i_fr, :) = SV%sver(:)

          ! Save retrieved CHI2
          results%chi2(i_fp, i_fr) = SUM(((radiance_meas_work - radiance_calc_work) ** 2) / &
               (noise_work ** 2)) / (N_spec - N_sv)

          ! Save statevector at last iteration
          do i=1, size(SV%svsv)
             results%sv_retrieved(i_fp, i_fr,i) = SV%svsv(i)
          end do

          results%num_iterations(i_fp, i_fr) = iteration

          final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_calc_work(:)
          measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_meas_work(:)
          noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)

          write(tmp_str, '(A,G3.1,A,F6.3,A,G2.1,A,F10.2)') "Iteration: ", iteration ,&
               ", Chi2: ",  results%chi2(i_fp, i_fr), &
               ", Active Levels: ", num_active_levels, &
               ", Psurf: ", this_psurf
          !call logger%debug(fname, trim(tmp_str))

       end if

       ! These quantities are all allocated within the iteration loop, and
       ! hence need explicit de-allocation.
       deallocate(radiance_meas_work, radiance_calc_work, radiance_tmp_work, &
            noise_work, Se_inv, K, solar_irrad)

       if (allocated(gas_tau)) deallocate(gas_tau)
       if (allocated(gas_tau_dpsurf)) deallocate(gas_tau_dpsurf)
       if (allocated(gas_tau_dvmr)) deallocate(gas_tau_dvmr)
       if (allocated(gas_tau_dsh)) deallocate(gas_tau_dsh)
       if (allocated(gas_tau_pert)) deallocate(gas_tau_pert)
       if (allocated(ray_tau)) deallocate(ray_tau)
       if (allocated(total_tau)) deallocate(total_tau)
       if (allocated(vmr_pert)) deallocate(vmr_pert)
       if (allocated(this_vmr_profile)) deallocate(this_vmr_profile)

       iteration = iteration + 1
       !read(*,*)
    end do

  end function physical_FM



  subroutine calculate_dispersion_limits(this_dispersion, i_win, &
       l1b_wl_idx_min, l1b_wl_idx_max)

    implicit none
    double precision, intent(in) :: this_dispersion(:)
    integer, intent(in) :: i_win
    integer, intent(inout) :: l1b_wl_idx_min, l1b_wl_idx_max

    integer :: i

    l1b_wl_idx_min = 0
    l1b_wl_idx_max = 0

    ! This here grabs the boundaries of the L1b data
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

    if (l1b_wl_idx_max < size(this_dispersion)) then
       l1b_wl_idx_max = l1b_wl_idx_max + 1
    end if

    if (l1b_wl_idx_max > size(this_dispersion)) then
       l1b_wl_idx_max = size(this_dispersion)
    end if


  end subroutine calculate_dispersion_limits

  subroutine read_atmosphere_file(filename, gas_indices, gas_strings, atm)
    ! This function is a little misplaced here, I would have like to have it
    ! in the File_utils module. However that creates a circular dependence with
    ! the Control module, so we just stick it here instead..

    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: gas_indices(:)
    type(string), intent(in) :: gas_strings(:)
    type(atmosphere), intent(inout) :: atm

    character(len=*), parameter :: fname = "read_atmosphere_file"
    integer :: funit, iostat
    logical :: file_exist
    integer :: line_count, nonempty_line_count, file_start, level_count
    integer :: idx_p, idx_t
    character(len=999) :: dummy, tmp_str
    double precision :: dummy_dp
    type(string) :: dummy_string
    type(string), allocatable :: split_string(:)

    integer :: i,j,cnt
    integer, allocatable :: this_gas_index(:)

    integer :: num_gases

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

    ! Loop through the file until we have reached EOF
    do
       read(funit, '(A)', iostat=iostat) dummy

       if (iostat == iostat_end) then
          ! End of file?
          exit
       end if

       line_count = line_count + 1

       if (scan(dummy, "!#;") > 0) then
          ! Skip this line, as it's commented
          cycle
       else if (trim(dummy) == "") then
          ! Skip empty lines
          cycle
       end if

       nonempty_line_count = nonempty_line_count + 1

       ! First non-empty line should be saved here
       if (file_start == -1) then
          file_start = line_count
       else
          if (nonempty_line_count > 1) then
             level_count = level_count + 1
          end if
       end if

    end do

    !close(funit)
    !open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')

    ! Go back to the top of the file.
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    do
       read(funit, '(A)', iostat=iostat) dummy

       line_count = line_count + 1

       if (line_count < file_start) cycle

       if (line_count == file_start) then
          ! This is the proper atmosphere header that should contain
          ! the information about the gases. So first, let's check if
          ! the numbers match
          dummy_string = dummy
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now that we know both the number of levels and gases, we can allocate the
          ! arrays in the atmosphere structure.

          num_gases = 0
          idx_p = -1
          idx_t = -1
          do j=1, size(split_string)
             ! Skip temp or pressure
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

          write(tmp_str, '(A,G0.1,A,A)') "There seem to be ", num_gases, " gases in ", filename
          call logger%info(fname, trim(tmp_str))

          if (allocated(atm%p)) deallocate(atm%p)
          if (allocated(atm%T)) deallocate(atm%T)
          if (allocated(atm%sh)) deallocate(atm%sh)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)

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


  function regrid_atmosphere(old_atm, psurf) result(new_atm)
    implicit none
    type(atmosphere), intent(in) :: old_atm
    double precision, intent(in) :: psurf
    type(atmosphere) :: new_atm

    integer :: i, N_lev, N_lev_new, N_gases

    N_lev = size(old_atm%p)
    ! If the surface pressure lies between the last and penultimate
    ! level, there's nothing to be done here.
    if ((psurf < old_atm%p(N_lev)) .and. (psurf > old_atm%p(N_lev-1))) then
       new_atm = old_atm
       return
    end if

    ! Otherwise, loop through the atmospheric levels to see where we need to cut off
    do i=1, N_lev-1
       if ((psurf > old_atm%p(i)) .and. (psurf < old_atm%p(i+1))) then
          ! Found
          exit
       end if
    end do

    ! Found the new number of levels
    N_lev_new = i+1
    N_gases = size(old_atm%gas_names)

    ! Allocate new atmosphere data, and in case we have an existing object,
    ! deallocate first

    if (allocated(new_atm%p)) deallocate(new_atm%p)
    if (allocated(new_atm%T)) deallocate(new_atm%T)
    if (allocated(new_atm%sh)) deallocate(new_atm%sh)
    if (allocated(new_atm%gas_names)) deallocate(new_atm%gas_names)
    if (allocated(new_atm%gas_index)) deallocate(new_atm%gas_index)
    if (allocated(new_atm%gas_vmr)) deallocate(new_atm%gas_vmr)

    allocate(new_atm%p(N_lev_new))
    allocate(new_atm%T(N_lev_new))
    allocate(new_atm%sh(N_lev_new))
    allocate(new_atm%gas_names(N_gases))
    allocate(new_atm%gas_index(N_gases))
    allocate(new_atm%gas_vmr(N_lev_new, N_gases))

    ! And copy over old data
    new_atm%p(:) = old_atm%p(1:N_lev_new)
    new_atm%T(:) = old_atm%T(1:N_lev_new)
    new_atm%sh(:) = old_atm%sh(1:N_lev_new)
    new_atm%gas_names(:) = old_atm%gas_names(:)
    new_atm%gas_index(:) = old_atm%gas_index(:)
    new_atm%gas_vmr(:,:) = old_atm%gas_vmr(1:N_lev_new, :)

  end function regrid_atmosphere



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

    allocate(results%num_iterations(num_fp, num_frames))
    allocate(results%converged(num_fp, num_frames))

    results%sv_names = "NONE"

    results%sv_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%chi2 = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%residual_rms = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%dsigma_sq = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

    results%num_iterations = -1
    results%converged = .false.

  end subroutine create_result_container



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
    deallocate(results%num_iterations)
    deallocate(results%converged)

  end subroutine destroy_result_container


  subroutine assign_SV_names_to_result(results, SV, i_win)
    implicit none
    type(result_container), intent(inout) :: results
    type(statevector), intent(in) :: SV
    integer, intent(in) :: i_win

    character(len=999) :: tmp_str
    integer :: i,j,k,l

    i = 1
    do while (i <= size(SV%svsv))

       ! Albedo names
       do j=1, SV%num_albedo
          if (SV%idx_albedo(j) == i) then
             write(tmp_str, '(A,G0.1)') "albedo_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       do j=1, SV%num_sif
          if (SV%idx_sif(j) == i) then
             write(tmp_str, '(A)') "SIF_absolute"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       if (SV%idx_solar_shift(1) == i) then
          write(tmp_str, '(A)') "solar_shift"
          results%sv_names(i) = trim(tmp_str)
       end if

       if (SV%idx_solar_stretch(1) == i) then
          write(tmp_str, '(A)') "solar_stretch"
          results%sv_names(i) = trim(tmp_str)
       end if

       do j=1, SV%num_psurf
          if (SV%idx_psurf(j) == i) then
             write(tmp_str, '(A)') "psurf"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       do j=1, SV%num_dispersion
          if (SV%idx_dispersion(j) == i) then
             write(tmp_str, '(A,G0.1)') "dispersion_coef_", j
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       k = 1
       do j=1, SV%num_gas


          ! Check if this SV element is a scalar retrieval
          if (SV%idx_gas(j,1) == i) then

             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(j))) then
                if (MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(j), k) == -1) cycle

                write(tmp_str, '(A,A)') trim(MCS%window(i_win)%gases(sv%gas_idx_lookup(j))%chars() // "_scale_")
                write(tmp_str, '(A, G0.3)') trim(tmp_str), &
                     SV%gas_retrieve_scale_start(j)
                write(tmp_str, '(A,A,G0.3)') trim(tmp_str), ":" , &
                     SV%gas_retrieve_scale_stop(j)
                results%sv_names(i) = trim(tmp_str)

                k = k + 1
             end if

          end if

          ! Check if it's a profile retrieval, and then loop through the profile
          ! elements to get each name.
          if (MCS%window(i_win)%gas_retrieve_profile(sv%gas_idx_lookup(j))) then
             do k=1, size(SV%idx_gas, 2)
                if (SV%idx_gas(j,k) == i) then

                   write(tmp_str,'(A,A,G0.1)') MCS%window(i_win)%gases(sv%gas_idx_lookup(j))%chars(), &
                        "_profile_", k
                   results%sv_names(i) = trim(tmp_str)

                end if
             end do
          end if

       end do

       i = i+1
    end do


  end subroutine assign_SV_names_to_result



end module physical_model_mod
