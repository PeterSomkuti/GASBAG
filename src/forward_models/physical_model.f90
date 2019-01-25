module physical_model_mod

  use control_mod, only: MCS, MAX_WINDOWS, MAX_GASES, MCS_find_gases
  use instruments, only: generic_instrument
  use oco2_mod
  use logger_mod, only: logger => master_logger
  use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
       read_DP_hdf_dataset
  use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  use solar_model_mod
  use math_utils_mod
  use statevector_mod
  use radiance_mod
  use absco_mod
  use gas_tau_mod

  use mod_datetime

  implicit none

  public physical_retrieval

  private

  ! A simple structure to keep the atmosphere data nice
  ! and tidy - one per sounding.
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
     double precision, allocatable :: sv_retrieved(:,:,:), sv_prior(:,:,:)
     double precision, allocatable, dimension(:,:) :: &
          chi2, residual_rms, dsigma_sq
     double precision, allocatable, dimension(:,:) :: num_iterations
     logical, allocatable :: converged(:,:)
  end type result_container


  double precision :: hires_spacing, hires_pad
  double precision, allocatable :: hires_grid(:)
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
  double precision, dimension(:,:,:,:), allocatable :: ils_hires_delta_lambda, &
       ils_hires_relative_response
  ! Averaged ILS in case of FFT convolution (wl, fp, band )
  double precision, dimension(:,:,:), allocatable :: ils_hires_avg_relative_response
  double precision :: ils_hires_min_wl, ils_hires_max_wl, ils_hires_spacing
  double precision, allocatable :: ils_hires_grid(:)
  integer :: num_ils_hires
  double precision :: ils_norm_factor

  double precision, allocatable :: dispersion_coefs(:,:,:)
  double precision, allocatable :: snr_coefs(:,:,:,:)
  double precision, allocatable :: land_fraction(:,:)


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
  double precision, allocatable :: solar_spectrum(:,:), solar_spectrum_regular(:,:)
  integer :: N_solar

  ! Final SIF value in physical radiance units, chi2
  double precision, dimension(:,:), allocatable :: retrieved_SIF_abs, &
       retrieved_SIF_rel, chi2, &
       num_iterations

  ! Radiances
  double precision, dimension(:,:,:), allocatable :: final_radiance, &
       measured_radiance, &
       noise_radiance

  ! State vector construct!
  type(statevector) :: SV
  ! Result container
  type(result_container) :: results


contains

  subroutine physical_retrieval(my_instrument)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument

    integer(hid_t) :: l1b_file_id, met_file_id, output_file_id, dset_id, result_gid
    integer(hsize_t) :: out_dims2d(2), out_dims3d(3)
    integer(hsize_t), dimension(:), allocatable :: dset_dims
    integer :: hdferr
    character(len=999) :: dset_name, tmp_str
    character(len=*), parameter :: fname = "physical_retrieval"

    logical :: gas_found, all_gases_found
    integer :: absco_dims

    integer :: num_frames
    integer :: i_fp, i_fr, i_pix, i_win, band
    integer :: i, j, retr_count
    logical :: this_converged


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

       !! Get the necesary MET data profiles. This probably needs expanding for
       !! other instrument types / file structures. (will Geocarb be different?)
       dset_name = "/Meteorology/vector_pressure_levels_met"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_P_levels, dset_dims)
       call logger%trivia(fname, "Finished reading in pressure levels.")

       dset_name = "/Meteorology/temperature_profile_met"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_T_profiles, dset_dims)
       call logger%trivia(fname, "Finished reading in temperature profiles.")

       dset_name = "/Meteorology/specific_humidity_profile_met"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_SH_profiles, dset_dims)
       call logger%trivia(fname, "Finished reading in specific humidity profiles.")

       dset_name = "/Meteorology/surface_pressure_met"
       call read_DP_hdf_dataset(met_file_id, dset_name, met_psurf, dset_dims)
       call logger%trivia(fname, "Finished reading in surface pressure.")

       dset_name = "/SoundingGeometry/sounding_land_fraction"
       call read_DP_hdf_dataset(l1b_file_id, dset_name, land_fraction, dset_dims)
       call logger%trivia(fname, "Finished reading in land fraction.")

       ! Read dispersion coefficients and create dispersion array
       call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
       allocate(dispersion(1016,8,3))

       do band=1,3
          do i_fp=1,8
             call my_instrument%calculate_dispersion(dispersion_coefs(:,i_fp,band), &
                  dispersion(:,i_fp,band), band, i_fp)
          end do
       end do

       ! Grab the SNR coefficients for noise calculations
       call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
       ! How many frames do we have in this file again?
       call my_instrument%read_num_frames(l1b_file_id, num_frames)
       ! Read in the sounding id's
       call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)
       ! Read the time strings
       call my_instrument%read_time_strings(l1b_file_id, sounding_time_strings)
       ! Read in the instrument ILS data
       call my_instrument%read_ils_data(l1b_file_id, ils_delta_lambda, &
            ils_relative_response)

       ! Read in bad sample list
       call my_instrument%read_bad_sample_list(l1b_file_id, bad_sample_list)

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
                     MCS%window(i_win)%wl_min - hires_pad, &
                     MCS%window(i_win)%wl_max + hires_pad)
             else
                call logger%fatal(fname, "Spectroscopy type: " // MCS%gas(j)%type // " not implemented!")
                stop 1
             end if
          end do

          ! Read the atmosphere file (if present) and populate the initial_atm structure
          ! with the contents of said file. It is cross-referenced against the gases in the
          ! list of gases in this window.

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
       hires_pad = (ils_hires_max_wl - ils_hires_min_wl) * 1.1d0

       ! Grab the desired high-resolution wavelength grid spacing
       hires_spacing = MCS%window(i_win)%wl_spacing

       ! .. and construct the high-resolution grid from the supplied
       ! microwindow range.
       N_hires = ceiling((MCS%window(i_win)%wl_max - MCS%window(i_win)%wl_min + 2*hires_pad) / hires_spacing)
       allocate(hires_grid(N_hires))
       do i=1, N_hires
          hires_grid(i) = MCS%window(i_win)%wl_min - hires_pad + dble(i-1) * hires_spacing
       end do

       ! Now that we have the appropriate high-resolution spacing, we can finish the ILS business
       num_ils_hires = ceiling((ils_hires_max_wl - ils_hires_min_wl) / hires_spacing)

       allocate(ils_hires_grid(num_ils_hires))
       allocate(ils_hires_delta_lambda(num_ils_hires, &
            size(ils_delta_lambda, 2), &
            size(ils_delta_lambda, 3), &
            size(ils_delta_lambda, 4)))
       allocate(ils_hires_relative_response, mold=ils_hires_delta_lambda)

       if (MCS%window(i_win)%fft_convolution) then
          allocate(ils_hires_avg_relative_response(num_ils_hires, &
               size(ils_relative_response, 3), &
               size(ils_relative_response, 4)))

       end if


       ! We need to ensure that the hi-res grid and the hi-res ILS grid line
       ! up fully, so ils_hires_min_wl needs to be shifted to the closest multiple
       ! of solar_grid_spacing.
       ils_hires_min_wl = hires_spacing * ceiling(ils_hires_min_wl / hires_spacing)

       ! And construct the ILS hires grid using the solar grid spacing
       do i=1, num_ils_hires
          ils_hires_grid(i) = ils_hires_min_wl + dble(i-1) * hires_spacing
       end do

       call logger%debug(fname, "Re-gridding ILS to hi-res grid.")

       do i_fp=1, my_instrument%num_fp
          do i_pix=1, size(dispersion, 1)

             ils_hires_delta_lambda(:, i_pix, i_fp, band) = ils_hires_grid(:)

             ! Interpolate the ILS onto the high-resolution grid
             call pwl_value_1d(&
                  size(ils_delta_lambda, 1), &
                  ils_delta_lambda(:, i_pix, i_fp, band), &
                  ils_relative_response(:, i_pix, i_fp, band), &
                  num_ils_hires, &
                  ils_hires_grid, ils_hires_relative_response(:, i_pix, i_fp, band) )

             ! Since we have resampled the ILS, it also needs re-normalising, a simple
             ! trapezoidal integration should do.

             ils_norm_factor = 0.0d0
             do i=1, size(ils_hires_relative_response, 1)-1
                ils_norm_factor = ils_norm_factor + &
                     (ils_hires_relative_response(i, i_pix, i_fp, band) + &
                     ils_hires_relative_response(i+1, i_pix, i_fp, band)) * &
                     (ils_hires_delta_lambda(i+1, i_pix, i_fp, band) - &
                     ils_hires_delta_lambda(i, i_pix, i_fp, band)) / 2.0d0
             end do

             ils_hires_relative_response(:, i_pix, i_fp, band) = &
                  ils_hires_relative_response(:, i_pix, i_fp, band) / ils_norm_factor

          end do

          ! If FFT convolution is performed, we have to assign one effective
          ! ILS for the entire band
          if (MCS%window(i_win)%fft_convolution) then

             ! The easiest approach is to simply average the ILS over all
             ! pixels for a given footprint and band.

             ! Dimensions: (delta lambda, fp, band)
             do i=1, size(ils_hires_relative_response, 1)
                ils_hires_avg_relative_response(i, i_fp, band) = &
                     sum(ils_hires_relative_response(i, :, i_fp, band)) / &
                     size(ils_hires_relative_response, 2)
             end do

          end if

       end do
       call logger%debug(fname, "Done re-gridding ILS.")


       ! Read in the solar model - we do this for every band, the reason being
       ! the following. Since the re-gridding procedure is fairly costly, we want
       ! to keep the solar spectrum data as small as possible. Hence we call all of this
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

       ! To make life easier, we want to keep the solar model on a regular,
       ! evenly-spaced grid in wavelength space. So we are going to interpolate
       ! it onto a slightly changed wavelength grid.
       allocate(solar_spectrum_regular(N_hires, 2))
       solar_spectrum_regular(:,1) = hires_grid

       call logger%debug(fname, "Re-gridding solar spectrum")
       call pwl_value_1d(size(solar_spectrum, 1), &
            solar_spectrum(:,1), solar_spectrum(:,2), &
            N_hires, &
            solar_spectrum_regular(:,1), solar_spectrum_regular(:,2))
       call logger%debug(fname, "Finished re-gridding solar spectrum.")
       ! Note that at this point, the solar spectrum is still normalised

       select type(my_instrument)
       type is (oco2_instrument)

          ! Read in the measurement geometries - these are the same for all three bands?
          call my_instrument%read_sounding_geometry(l1b_file_id, band, SZA, SAA, VZA, VAA)
          ! Read in the measurement location
          call my_instrument%read_sounding_location(l1b_file_id, band, lon, lat, &
               altitude, relative_velocity, &
               relative_solar_velocity)

          ! Read in Spike filter data
          call my_instrument%read_spike_filter(l1b_file_id, spike_list, band)
       end select

       ! Allocate containers to hold the retrieval results
       allocate(retrieved_SIF_abs(my_instrument%num_fp, num_frames))
       allocate(retrieved_SIF_rel(my_instrument%num_fp, num_frames))
       allocate(num_iterations(my_instrument%num_fp, num_frames))
       allocate(chi2(my_instrument%num_fp, num_frames))

       allocate(final_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
       allocate(measured_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
       allocate(noise_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))

       num_iterations = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       chi2 = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       retrieved_SIF_abs = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       retrieved_SIF_rel = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

       final_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       measured_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       noise_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)


       ! Set up state vector structure here

       ! We can do this once for every window, and simply clear the contents
       ! of the state vectors inside the loop later on. Might not save an awful
       ! lot, but every little helps, I guess..

       ! At the moment, we can do following state vectors:
       !
       ! Lambertian surface albedo, arbitrary (?) polynomial order
       ! Zero-level offset / SIF (constant)
       ! Instrument disperison
       ! Surface pressure

       ! Parsing the statevector string, that was passed in the window
       ! section and initialize the state vector SV accordingly. This subroutine
       ! needs to access plenty of things in the MCS, so we only pass the
       ! window index, and the routine takes care of arranging the rest.

       ! For the beginning, we start by initialising it with the number of levels
       ! as obtained from the initial_atm

       call parse_and_initialize_SV(i_win, size(initial_atm%p), SV)
       call logger%info(fname, "Initialised SV structure")


       ! And now set up the result container with the appropriate sizes for arrays
       call create_result_container(results, num_frames, my_instrument%num_fp, size(SV%svap))

       retr_count = 0
       mean_duration = 0.0d0

       do i_fr=1, 10 !num_frames
          do i_fp=1, my_instrument%num_fp

             !if (land_fraction(i_fp, i_fr) < 0.95) then
             !   cycle
             !end if

             call cpu_time(cpu_time_start)
             this_converged = physical_FM(my_instrument, i_fp, i_fr, i_win, band)
             call cpu_time(cpu_time_stop)

             retr_count = retr_count + 1
             !mean_duration = ((mean_duration * retr_count) + (cpu_time_stop - cpu_time_start)) / (retr_count)

             if (mod(retr_count, 1) == 0) then
                !mean_duration = (cpu_time_stop - cpu_time_start)
                write(tmp_str, '(A, G0.1, A, G0.1, A, F10.5, A, L1)') &
                     "Frame/FP: ", i_fr, "/", i_fp, " - ", cpu_time_stop-cpu_time_start, ' - ', this_converged
                call logger%debug(fname, trim(tmp_str))
                !read(*,*)
             end if
          end do
       end do

       deallocate(solar_spectrum_regular, solar_spectrum, hires_grid, ils_hires_grid, &
            ils_hires_delta_lambda, ils_hires_relative_response)

       if (MCS%window(i_win)%fft_convolution) deallocate(ils_hires_avg_relative_response)


       out_dims2d = shape(chi2)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/retrieved_chi2_" // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            chi2, out_dims2d, -9999.99d0)

       out_dims2d = shape(num_iterations)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/retrieved_num_iterations_" // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            num_iterations, out_dims2d, -9999.99d0)

       if (SV%num_sif == 1) then
          out_dims2d = shape(retrieved_sif_abs)
          write(tmp_str, '(A,A)') "/physical_retrieval_results/retrieved_sif_abs_" // MCS%window(i_win)%name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               retrieved_sif_abs, out_dims2d, -9999.99d0)

          out_dims2d = shape(retrieved_sif_rel)
          write(tmp_str, '(A,A)') "/physical_retrieval_results/retrieved_sif_rel_" // MCS%window(i_win)%name
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               retrieved_sif_rel, out_dims2d, -9999.99d0)
       end if

       out_dims3d = shape(final_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/modelled_radiance_" // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            final_radiance, out_dims3d)

       out_dims3d = shape(measured_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/measured_radiance_" // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            measured_radiance, out_dims3d)

       out_dims3d = shape(noise_radiance)
       write(tmp_str, '(A,A)') "/physical_retrieval_results/noise_radiance_" // MCS%window(i_win)%name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            noise_radiance, out_dims3d)

       ! Clear and deallocate the SV structure to be ready for the next window.
       ! We really shouldn't need to check if this is deallocated already or not,
       ! since parse_and_initialize_SV should have successfully allocated
       ! all fields in SV.
       call clear_SV(SV)
       call logger%info(fname, "Clearing up SV structure")

       ! Also deallocate containers holding the retrieval results
       deallocate(retrieved_SIF_abs)
       deallocate(retrieved_SIF_rel)
       deallocate(num_iterations)
       deallocate(chi2)

       deallocate(final_radiance)
       deallocate(measured_radiance)
       deallocate(noise_radiance)

       !deallocate(bad_sample_list)
       !deallocate(spike_list)

       call destroy_result_container(results)

    end do


  end subroutine physical_retrieval


  function physical_FM(my_instrument, i_fp, i_fr, i_win, band) result(converged)

    implicit none

    class(generic_instrument), intent(in) :: my_instrument
    integer, intent(in) :: i_fr, i_fp, i_win, band
    logical :: converged

    !!
    integer(hid_t) :: l1b_file_id, output_file_id

    !! Radiances
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

    logical :: skip_jacobian

    !! Solar stuff
    double precision :: solar_dist, solar_rv, earth_rv, solar_doppler
    double precision, dimension(:,:), allocatable :: this_solar
    double precision, allocatable :: solar_irrad(:)

    !! Atmosphere
    integer :: num_gases, num_levels, num_active_levels
    double precision, allocatable :: gas_tau(:,:,:), gas_tau_dvmr(:,:,:), &
         gas_tau_dpsurf(:,:,:), gas_tau_pert(:,:,:), gas_tau_dpsurf2(:,:,:)
    double precision :: gas_scaling_factor

    !! Albedo
    double precision :: albedo_apriori
    double precision, allocatable :: albedo(:)

    !! Surface pressure
    double precision :: psurf, this_psurf
    logical :: do_psurf_jac

    !! SIF?
    double precision :: this_sif_radiance

    !! GASES
    logical :: do_gas_jac

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
    ! Levenberg-Marquart gamma
    double precision :: lm_gamma
    ! Iteration-related
    integer :: iteration, max_iterations
    logical :: keep_iterating
    logical :: log_retrieval
    logical :: success_inv_mat


    !! Misc
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "physical_FM"
    integer :: N_spec, N_spec_hi, N_spec_tmp, N_sv
    integer :: i, j
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
    doy_dp = dble(date%yearday()) + dble(date%getHour()) / 24.0


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
       call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, band, radiance_l1b)
    end select

    allocate(this_dispersion(size(radiance_l1b)))
    allocate(this_dispersion_tmp(size(radiance_l1b)))
    ! The dispersion coefficients use to generate this soundings' dispersion
    ! has the same size as the L1b dispersion coefficient array
    allocate(this_dispersion_coefs(size(dispersion_coefs, 1)))
    allocate(this_dispersion_coefs_pert(size(dispersion_coefs, 1)))


    ! Allocate the micro-window bounded solar and radiance arrays
    allocate(this_solar(N_hires, 2))
    ! If solar doppler-shift is needed, calculate the distance and relative
    ! velocity between point of measurement and the (fixed) sun

    call calculate_solar_distance_and_rv(doy_dp, solar_dist, solar_rv)
    call calculate_rel_velocity_earth_sun(lat(i_fp, i_fr), SZA(i_fp, i_fr), &
         SAA(i_fp, i_fr), altitude(i_fp, i_fr), earth_rv)

    ! solar_doppler =  (solar_rv * 1000.0d0 + earth_rv) / SPEED_OF_LIGHT
    solar_doppler = (earth_rv + solar_rv * 1000.0d0) / SPEED_OF_LIGHT
    !solar_doppler = (relative_solar_velocity(i_fp, i_fr)) / SPEED_OF_LIGHT
    ! Take a copy of the solar spectrum and re-adjust the solar spectrum wavelength grid
    ! According to the Doppler shift
    this_solar(:,1) = solar_spectrum_regular(:, 1) / (1.0d0 - solar_doppler)

    ! Since we now have the solar distance, we can calculate the proper solar
    ! continuum for the given day of the year.
    allocate(solar_irrad(size(this_solar, 1)))
    call calculate_solar_planck_function(6500.d0, solar_dist * 1000.d0, &
         this_solar(:,1), solar_irrad)
    ! And multiply to get the solar irradiance in physical units
    this_solar(:,2) = solar_spectrum_regular(:, 2) * solar_irrad(:)

    ! Correct for instrument doppler shift
    instrument_doppler = relative_velocity(i_fp, i_fr) / SPEED_OF_LIGHT

    ! write(*,*) "Earth RV: ", earth_rv
    ! write(*,*) "Solar RV: ", solar_rv * 1000
    ! write(*,*) "Relative velocity solar: ", relative_solar_velocity(i_fp, i_fr)
    ! write(*,*) "Relative velocity instrument: ", relative_velocity(i_fp, i_fr)
    !
    !
    ! write(*,*) "Solar doppler: ", solar_doppler
    ! write(*,*) "Instrument doppler: ", instrument_doppler


    ! Set up retrieval quantities:
    N_sv = size(SV%svap)
    N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1
    N_spec_hi = size(this_solar, 1)
    lm_gamma = 5.0d0

    ! Output-resolution K is allocated within the loop, as the
    ! number of pixels might change?
    allocate(K_hi(N_spec_hi, N_sv))
    K_hi(:,:) = 0.0d0

    allocate(Sa_inv(N_sv, N_sv))
    allocate(Shat_inv(N_sv, N_sv))
    allocate(Shat(N_sv, N_sv))
    allocate(tmp_m1(N_sv, N_sv), tmp_m2(N_sv, N_sv))
    allocate(tmp_v1(N_sv), tmp_v2(N_sv))

    ! Inverse prior covariance
    Sa_inv(:,:) = 0.0d0
    if (SV%num_albedo > 0) then
       do i=1, SV%num_albedo
          Sa_inv(SV%idx_albedo(i), SV%idx_albedo(i)) = 1.0d0
       end do
    end if

    if (SV%num_sif > 0) then
       Sa_inv(SV%idx_sif(1), SV%idx_sif(1)) = 1.0d0 / (1.0d18 ** 2)
    end if

    ! Populate state vector priors (this should also be read from config)

    ! Albedo
    select type(my_instrument)
    type is (oco2_instrument)
       ! OCO-2 has Stokes coefficient 0.5 for intensity, so we need to
       ! take that into account for the incoming solar irradiance

       albedo_apriori = PI * maxval(radiance_l1b) / &
            (1.0d0 * maxval(this_solar(:,2)) * cos(DEG2RAD * SZA(i_fp, i_fr)))

    end select

    if (albedo_apriori > 1) then
       write(tmp_str, '(A, F8.5)') "Albedo too large: ", albedo_apriori
       call logger%warning(fname, trim(tmp_str))
    end if


    ! TODO: What do we do in case of invalid albedo (<0, >1)?

    SV%svap(SV%idx_albedo(1)) = albedo_apriori
    ! Set slope etc. to zero always (why would we ever want to have a prior slope?)
    if (SV%num_albedo > 1) then
       do i=2, SV%num_albedo
          SV%svap(SV%idx_albedo(i)) = 0.0d0
       end do
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

    ! Surface pressure
    if (SV%num_psurf == 1) then
       SV%svap(SV%idx_psurf(1)) = met_psurf(i_fp, i_fr)
    end if

    ! Gases
    if (SV%num_gas > 0) then
       do i=1, SV%num_gas

          if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
             SV%svap(SV%idx_gas(i,1)) = 1.0d0
          end if
          if (MCS%window(i_win)%gas_retrieve_profile(sv%gas_idx_lookup(i))) then
             SV%svap(SV%idx_gas(i,:)) = initial_atm%gas_vmr(:, sv%gas_idx_lookup(i))
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




    converged = .false.
    ! Main iteration loop for the retrieval process.
    do while (keep_iterating)

       if (iteration == 1) then
          ! For the first iteration, we want to use the prior albedo
          albedo(:) = albedo_apriori
          ! and the state vector is the prior (maybe we want to change that later..)
          SV%svsv = SV%svap

          ! On the first iteration, we populate the this_atm structure by grabbing the
          ! initial atmosphere, and potentially re-gridding it - depending on where the
          ! surface pressure lies. We use the initial surface pressure coming from MET.
          !this_atm = regrid_atmosphere(initial_atm, met_psurf(i_fp, i_fr))
          if (num_gases > 0) then

             this_psurf = met_psurf(i_fp, i_fr)
             this_atm = initial_atm
             num_levels = size(this_atm%p)

             ! And get the T and SH profiles onto our new atmosphere grid
             call pwl_value_1d(size(met_P_levels, 1), &
                  log(met_P_levels(:,i_fp,i_fr)), met_T_profiles(:,i_fp,i_fr), &
                  size(this_atm%p), log(this_atm%p), this_atm%T)

             call pwl_value_1d(size(met_P_levels, 1), &
                  log(met_P_levels(:,i_fp,i_fr)), met_SH_profiles(:,i_fp,i_fr), &
                  size(this_atm%p), log(this_atm%p), this_atm%sh)

             do j=1, num_levels
                if (this_psurf > this_atm%p(j)) then
                   num_active_levels = j+1
                end if

                if (this_atm%sh(j) < 0.0d0) this_atm%sh(j) = 1.0d-10
             end do

             !this_atm%sh = 0.0d0
          end if
       else
          ! Otherwise, calculate it from the state vector
          ! Obviously ONLY if we want to retrieve it, otherwise the albedo
          ! will just stay at the prior value.
          if (SV%num_albedo > 0) then
             albedo = 0.0d0
             do i=1, SV%num_albedo
                albedo = albedo + SV%svsv(SV%idx_albedo(i)) * ((this_solar(:, 1) - this_solar(1,1)) ** (dble(i-1)))
             end do
          endif

          ! If we are retrieving surface pressure, we possibly need
          ! to re-grid the atmosphere if the surface pressure jumps to the
          ! next layer.

          if (SV%num_psurf == 1) then
             this_psurf = SV%svsv(SV%idx_psurf(1))

             if (this_psurf < this_atm%p(size(this_atm%p) - 1)) then
                call pwl_value_1d(size(met_P_levels, 1), &
                     log(met_P_levels(:,i_fp,i_fr)), met_T_profiles(:,i_fp,i_fr), &
                     size(this_atm%p), log(this_atm%p), this_atm%T)
                call pwl_value_1d(size(met_P_levels, 1), &
                     log(met_P_levels(:,i_fp,i_fr)), met_SH_profiles(:,i_fp,i_fr), &
                     size(this_atm%p), log(this_atm%p), this_atm%sh)
             end if

             if (this_psurf > this_atm%p(size(this_atm%p))) then
                write(*,*) "Psurf: ", SV%svsv(SV%idx_psurf(1))
                write(*,*) "Lowest p: ", this_atm%p(size(this_atm%p))
                return
             end if

             ! The number of active levels has to be inferred for every
             ! iteration if we are retrieving surface pressure 
             do j=1, num_levels
                if (this_psurf > this_atm%p(j)) then
                   num_active_levels = j+1
                end if
             end do

          end if
       endif


       !this_atm%sh = 1d-5
       !this_atm%T = 220.0d0

       ! Save the old state vector (iteration - 1'th state vector)
       old_sv = SV%svsv

       if (SV%num_sif > 0) then
          this_sif_radiance = SV%svsv(SV%idx_sif(1))
       else
          this_sif_radiance = 0.0d0
       end if

       !write(*,*) num_active_levels
       !write(*,*) this_atm%p

       do j=1, num_active_levels
       !   write(*,*) j, this_atm%p(j), this_atm%T(j), this_atm%sh(j), &
       !        (this_atm%gas_vmr(j,1), i=1, size(this_atm%gas_vmr, 2))
       end do
 
       ! Heavy bit - calculate the optical properties given an atmosphere with gases
       ! and their VMRs. This branch of the code will only be entered if we have at least
       ! one gas present. Otherwise, gas_tau will stay unallocated.

       if (num_gases > 0) then

          allocate(gas_tau(size(this_solar, 1), num_levels-1, num_gases))
          allocate(gas_tau_dpsurf(size(this_solar, 1), num_levels-1, num_gases))
          allocate(gas_tau_dpsurf2(size(this_solar, 1), num_levels-1, num_gases))
          allocate(gas_tau_dvmr(size(this_solar, 1), num_levels-1, num_gases))
          allocate(gas_tau_pert(size(this_solar, 1), num_levels-1, num_gases))

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

          do j=1, num_gases

             do_gas_jac = .false.

             if (SV%num_gas > 0) then

                gas_scaling_factor = 1.0d0

                ! We need to 'reverse-lookup' to see which SV index belongs to this
                ! gas to grab the right scaling factor.
                do i=1, SV%num_gas
                   if (SV%gas_idx_lookup(i) == j) then
                      do_gas_jac = .true.
                      if (MCS%window(i_win)%gas_retrieve_scale(j)) then
                         gas_scaling_factor = SV%svsv(SV%idx_gas(j,1))
                      end if
                   end if
                end do
             else
                ! If no gases are retrieved, we keep this at 1!
                gas_scaling_factor = 1.0d0
             end if

             call calculate_gas_tau( &
                  this_solar(:,1), &
                  this_atm%gas_vmr(:,j) * gas_scaling_factor, &
                  psurf, &
                  this_atm%p(:), &
                  this_atm%T(:), &
                  this_atm%sh(:), &
                  MCS%gas(MCS%window(i_win)%gas_index(j)), &
                  MCS%window(i_win)%N_sublayers, &
                  do_psurf_jac, &
                  do_gas_jac, & ! In case we ever want per-layer gas jacobians
                  gas_tau(:,:,j), &
                  gas_tau_dpsurf(:,:,j), &
                  gas_tau_dvmr(:,:,j))
          end do

!!$          write(tmp_str,'(A,G0.1,A)') "gas_od_iter", iteration, ".dat"
!!$          open(newunit=funit, file=trim(tmp_str))
!!$          do j=1, size(gas_tau, 1)
!!$             write(funit, *) (gas_tau_dpsurf(j,num_active_levels-1,i), i=1,size(gas_tau,3))
!!$          end do
!!$          close(funit)
!!$
!!$          write(tmp_str,'(A,G0.1,A)') "gas_dvmr_iter", iteration, ".dat"
!!$          open(newunit=funit, file=trim(tmp_str))
!!$          do j=1, size(gas_tau, 1)
!!$             write(funit, *) (gas_tau_dvmr(j,num_active_levels-1,i), i=1,size(gas_tau_dvmr,3))
!!$          end do
!!$          close(funit)

       end if


       ! Calculate the sun-normalized TOA radiances and store them in
       ! 'radiance_calc_work_hi'
       call calculate_radiance(this_solar(:,1), SZA(i_fp, i_fr), &
            VZA(i_fp, i_fr), albedo, gas_tau, &
            radiance_calc_work_hi)

       ! And multiply with the solar spectrum for physical units and add SIF contributions
       radiance_calc_work_hi(:) = this_solar(:,2) * radiance_calc_work_hi(:)
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:) + this_sif_radiance

       ! this probably should go into the radiance module?
       if (SV%num_albedo > 0) then
          do i=1, SV%num_albedo
             K_hi(:, SV%idx_albedo(i)) = (radiance_calc_work_hi(:) - this_sif_radiance) / albedo * &
                  ((this_solar(:,1) - this_solar(1,1)) ** (dble(i-1)))
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


       ! Gas jacobians
       if (SV%num_gas > 0) then
          do i=1, SV%num_gas

             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(i))) then
                ! This is a scale-type jacobian

                K_hi(:, SV%idx_gas(i,1)) = -(radiance_calc_work_hi(:) - this_sif_radiance) &
                     * (1.0d0 / cos(DEG2RAD * SZA(i_fp, i_fr)) + 1.0d0 / cos(DEG2RAD * VZA(i_fp, i_fr))) &
                     * sum(gas_tau(:,:,SV%gas_idx_lookup(i)), dim=2) / SV%svsv(SV%idx_gas(i,1))

             end if

             if (MCS%window(i_win)%gas_retrieve_profile(sv%gas_idx_lookup(i))) then

                do j=1, num_active_levels-1

                   K_hi(:, SV%idx_gas(i,j)) = -(radiance_calc_work_hi(:) - this_sif_radiance) &
                        * (1.0d0 / cos(DEG2RAD * SZA(i_fp, i_fr)) + 1.0d0 / cos(DEG2RAD * VZA(i_fp, i_fr))) &
                        * gas_tau_dvmr(:,j,SV%gas_idx_lookup(i)) 

                end do

             end if

          end do
       end if


       ! Stokes coefficients
       radiance_calc_work_hi(:) = radiance_calc_work_hi(:)
       K_hi(:,:) = K_hi(:,:)

       ! Dispersion
       this_dispersion_coefs(:) = dispersion_coefs(:, i_fp, band)
       ! If required, replace L1b dispersion coefficients by state vector
       ! elements from the retrieval process

       if (SV%num_dispersion > 0) then
          do i=1, SV%num_dispersion
             this_dispersion_coefs(i) = SV%svsv(SV%idx_dispersion(i))
          end do
       end if

       ! Calculate the wavelength grid
       select type(my_instrument)
       type is (oco2_instrument)
          call my_instrument%calculate_dispersion(this_dispersion_coefs, &
               this_dispersion(:), band, i_fp)
       end select

       ! Grab a copy of the dispersion relation that's shifted according to
       ! spacecraft velocity.
       this_dispersion = this_dispersion / (1.0d0 - instrument_doppler)

       ! Here we grab the index limits for the radiances for
       ! the choice of our microwindow and the given dispersion relation
       call calculate_dispersion_limits(this_dispersion, i_win, l1b_wl_idx_min, l1b_wl_idx_max)

       N_spec = l1b_wl_idx_max - l1b_wl_idx_min + 1

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
          ! that they are not considered in the fit. This should
          ! save otherwise good spectra with just a few distorted
          ! radiance values.

          do i=1, N_spec
             if (spike_list(i + l1b_wl_idx_min - 1, i_fp, i_fr) >= 5) then
                noise_work(i) = noise_work(i) * 10000.0d0
             end if
          end do

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
!!$          write(funit,*) (K_hi(i, j), j=1, N_sv)
!!$       end do
!!$       close(funit)

!!$       open(file="hires_spec.dat", newunit=funit)
!!$       do i=1, size(this_solar, 1)
!!$          write(funit,*) this_solar(i,1), this_solar(i,2), solar_spectrum_regular(i,1), &
!!$               solar_spectrum_regular(i,2), radiance_calc_work_hi(i) !, K_hi(i, SV%idx_gas(1)), K_hi(i, SV%idx_gas(2))
!!$       end do
!!$       close(funit)
       ! ! !


       select type(my_instrument)
       type is (oco2_instrument)

          ! Convolution of the TOA radiances

          if (MCS%window(i_win)%fft_convolution) then
             call fft_convolution(radiance_calc_work_hi, ils_hires_avg_relative_response(:,i_fp,band), &
                  this_solar(:,1), this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), radiance_calc_work)
          else
             call oco_type_convolution(this_solar(:,1), radiance_calc_work_hi, &
                  ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                  this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), radiance_calc_work, &
                  ILS_success)

             if (.not. ILS_success) then
                call logger%error(fname, "ILS convolution error.")
                return
             end if
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

             ! SIF does not need convolution
             if (i == SV%idx_sif(1)) skip_jacobian = .true.

             if (skip_jacobian) cycle

             ! Otherwise just convolve the other Jacobians and save the result in
             ! the low-resolution Jacobian matrix 'K'
             if (MCS%window(i_win)%fft_convolution) then

                call fft_convolution(K_hi(:,i), ils_hires_avg_relative_response(:,i_fp,band), &
                     this_solar(:,1), this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), K(:,i))

             else

                call oco_type_convolution(this_solar(:,1), K_hi(:,i), &
                     ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                     ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                     this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), K(:,i), &
                     ILS_success)

                if (.not. ILS_success) then
                   call logger%error(fname, "ILS convolution error.")
                   return
                end if
             end if

          end do
       end select

       ! Disperion Jacobians are produced via finite differencing
       if (SV%num_dispersion > 0) then
          do i=1, SV%num_dispersion

             ! Set prior inverse covariance (again all diagonal for now)
             ! Sa_inv(SV%idx_dispersion(i), SV%idx_dispersion(i)) = 1.0d0 / MCS%window(i_win)%dispersion_cov(i)

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

                if (MCS%window(i_win)%fft_convolution) then

                   call fft_convolution(radiance_calc_work_hi, &
                        ils_hires_avg_relative_response(:,i_fp,band), &
                        this_solar(:,1), this_dispersion_tmp(l1b_wl_idx_min:l1b_wl_idx_max), &
                        radiance_tmp_work)

                else

                   call oco_type_convolution(this_solar(:,1), radiance_calc_work_hi, &
                        ils_hires_delta_lambda(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                        ils_hires_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
                        this_dispersion_tmp(l1b_wl_idx_min:l1b_wl_idx_max), radiance_tmp_work, &
                        ILS_success)

                   if (.not. ILS_success) then
                      call logger%error(fname, "ILS convolution error.")
                      return
                   end if

                end if

                ! Compute and store the Jacobian into the right place in the matrix
                K(:, SV%idx_dispersion(i)) = (radiance_tmp_work - radiance_calc_work) &
                     / MCS%window(i_win)%dispersion_pert(i)
             end select

          end do
       end if

!!$       open(file="jacobian.dat", newunit=funit)
!!$       do i=1, N_spec
!!$          write(funit,*) (K(i, j), j=1, N_sv)
!!$       end do
!!$       close(funit)



       ! See Rodgers (2000) equation 5.36: calculating x_i+1 from x_i
       !tmp_m1 = (1.0d0 + lm_gamma) * Sa_inv + matmul(matmul(transpose(K), Se_inv), K)
       tmp_m1 = matmul(matmul(transpose(K), Se_inv), K)

       call invert_matrix(tmp_m1, tmp_m2, success_inv_mat)
       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert K^T Se K")
          return
       end if


       tmp_v1 = matmul(matmul(transpose(K), Se_inv), radiance_meas_work - radiance_calc_work)
       !tmp_v2 = matmul(Sa_inv, SV%svsv - SV%svap)

       ! Update state vector
       SV%svsv = SV%svsv + matmul(tmp_m2, tmp_v1) !  - tmp_v2)

       ! Calculate Shat_inv
       Shat_inv = matmul(matmul(transpose(K), Se_inv), K)! + Sa_inv
       call invert_matrix(Shat_inv, Shat, success_inv_mat)

       if (.not. success_inv_mat) then
          call logger%error(fname, "Failed to invert Shat^-1")
          return
       end if

       ! Check delta sigma square for this iteration
       dsigma_sq = dot_product(old_sv - SV%svsv, matmul(Shat_inv, old_sv - SV%svsv))

       do i=1, N_sv
          SV%sver(i) = sqrt(Shat(i,i))
       end do


       open(file="l1b_spec.dat", newunit=funit)
       do i=1, N_spec
          write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), radiance_calc_work(i), &
               noise_work(i) !, K(i, SV%idx_gas(1)), K(i, SV%idx_gas(2)) !, K(i, SV%idx_psurf(1)), K(i, SV%idx_dispersion(1)), K(i, SV%idx_dispersion(2))
       end do
       close(funit)

!!$
!!$
!!$       write(*,*) "old, current and delta state vector, and errors"
!!$       write(*,*) "Iteration: ", iteration-1
!!$       do i=1, N_sv
!!$          write(*,*) i, old_sv(i), SV%svsv(i), SV%svsv(i) - old_sv(i), sqrt(Shat(i,i))
!!$       end do
!!$       write(*,*) "Chi2: ", SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)
!!$       write(*,*) "Dsigma2: ", dsigma_sq

!!$

       if (dsigma_sq < dble(N_sv) * dsigma_scale) then
          keep_iterating = .false.
          !write(*,*) "Converged!"
          !write(*,*) "Chi2: ", SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)

          do i=1, N_sv
             SV%sver(i) = sqrt(Shat(i,i))
          end do

          ! Save values
          chi2(i_fp, i_fr) = SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)

          if (SV%num_sif == 1) then
             retrieved_sif_abs(i_fp, i_fr) = SV%svsv(SV%idx_sif(1))
             retrieved_sif_rel(i_fp, i_fr) = retrieved_sif_abs(i_fp, i_fr) / maxval(radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max))
          end if
          num_iterations(i_fp, i_fr) = iteration

          final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_calc_work(:)
          measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_meas_work(:)
          noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)

          write(tmp_str, '(A,G0.1,A,F10.3)') "Iteration: ", iteration ,", Chi2: ",  chi2(i_fp, i_fr)
          call logger%debug(fname, trim(tmp_str))

          converged = .true.

!!$       open(file="l1b_spec.dat", newunit=funit)
!!$       do i=1, N_spec
!!$          write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), radiance_calc_work(i), &
!!$               noise_work(i), K(i, SV%idx_psurf(1))
!!$       end do
!!$       close(funit)

          !if (i_fr == 51) then
          ! end if
          !read(*,*)
       end if

       if (iteration == 5) then
          call logger%warning(fname, "Not converged!")
          keep_iterating = .false.

          return

          ! Print out SV for visual inspection!

!!$          do i=1, N_sv
!!$             write(*,*) i, old_sv(i), SV%svsv(i), SV%svsv(i) - old_sv(i), sqrt(Shat(i,i))
!!$          end do
!!$          write(*,*) "Chi2: ", SUM(((radiance_meas_work - radiance_calc_work) ** 2) / (noise_work ** 2)) / (N_spec - N_sv)
!!$
!!$          write(*,*) "Dsigma_sq: ", dsigma_sq

!!$          open(file="l1b_spec.dat", newunit=funit)
!!$          do i=1, N_spec
!!$             write(funit,*) this_dispersion(i+l1b_wl_idx_min-1), radiance_meas_work(i), radiance_calc_work(i), &
!!$                  noise_work(i), bad_sample_list(i+l1b_wl_idx_min-1,i_fp,band)
!!$          end do
!!$          close(funit)
!!$
!!$          read(*,*)
       end if

       ! These quantities are all allocated within the iteration loop, and
       ! hence need explicit de-allocation.
       deallocate(radiance_meas_work, radiance_calc_work, radiance_tmp_work, &
            noise_work,Se_inv, K)

       if (allocated(gas_tau)) deallocate(gas_tau)
       if (allocated(gas_tau_dpsurf)) deallocate(gas_tau_dpsurf)
       if (allocated(gas_tau_dvmr)) deallocate(gas_tau_dvmr)
       if (allocated(gas_tau_dpsurf2)) deallocate(gas_tau_dpsurf2)
       if (allocated(gas_tau_pert)) deallocate(gas_tau_pert)

       iteration = iteration + 1

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
    end if

    ! First pass: we scan the file to see how many levels our atmosphere has
    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
    line_count = 0
    level_count = 0
    nonempty_line_count = 0
    file_start = -1
    do
       read(funit, '(A)', iostat=iostat) dummy
       if (iostat < 0) then
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
    close(funit)

    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
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

          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))
          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))

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
          deallocate(split_string)
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



  subroutine create_result_container(results, num_frames, num_fp, num_SV)
    implicit none
    type(result_container), intent(inout) :: results
    integer, intent(in) :: num_frames, num_fp, num_SV

    allocate(results%sv_names(num_SV))

    allocate(results%sv_retrieved(num_frames, num_fp, num_SV))
    allocate(results%sv_prior(num_frames, num_fp, num_SV))
    allocate(results%chi2(num_frames, num_fp))
    allocate(results%residual_rms(num_frames, num_fp))
    allocate(results%dsigma_sq(num_frames, num_fp))

    allocate(results%num_iterations(num_frames, num_fp))
    allocate(results%converged(num_frames, num_fp))

    results%sv_names = "NONE"

    results%sv_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
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
    deallocate(results%chi2)
    deallocate(results%residual_rms)
    deallocate(results%dsigma_sq)
    deallocate(results%num_iterations)
    deallocate(results%converged)

  end subroutine destroy_result_container


  subroutine assign_SV_names_to_result(results, SV)
    implicit none
    type(result_container), intent(inout) :: results
    type(statevector), intent(in) :: SV

    character(len=999) :: tmp_str
    integer :: i,j

    do i=1, size(SV%svsv)

       do j=1, SV%num_albedo
          if (SV%idx_albedo(j) == i) then
             write(tmp_str, '(A,G0.1)') "albedo_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

    end do


  end subroutine assign_SV_names_to_result



end module physical_model_mod
