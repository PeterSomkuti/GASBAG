!> @brief Guanter-type empirical SIF retrieval method
!> @file guanter_model.f90
!> @author Peter Somkuti
!>
!! This module contains the subroutines to execute the empirical, data-driven SIF
!! retrieval scheme, outlined in Guanter et al. (2012). There will be a good chunk
!! in here, which is redundant when compared to the physical-based retrieval module,
!! however since the two retrieval methods are so different - I decided to put them
!! into separate modules at the cost of some code doubling-up.

module guanter_model_mod

  ! User modules
  use control_mod, only: MCS, MAX_WINDOWS
  use instruments_mod, only: generic_instrument
  use oco2_mod
  use math_utils_mod

  ! Third party modules
  use logger_mod, only: logger => master_logger
  use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset

  ! System modules
  use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  use OMP_LIB
  use HDF5

  implicit none

  public :: guanter_retrieval
  private :: guanter_fm, slope_correction


  ! In the Guanter-scheme, dispersion is not really touched, hence we can just
  ! keep it as a module-wide fixed set of numbers (pixel, footprint, band)
  !> Dispersion array to hold the wavelength-pixel mapping (detector pixels, footprint, band).
  double precision, allocatable :: dispersion(:,:,:)
  !> Dispersion coefficient array that holds the polynomial coefficients
  !> required to construct the dispersion/wavelength arrays.
  double precision, allocatable :: dispersion_coefs(:,:,:)
  !> SNR coefficient array required to calculate instrument
  !> noise equivalent radiances.
  double precision, allocatable :: snr_coefs(:,:,:,:)

  !> Basisfunctions are stored as (window, spectral pixel, footprint, #SV)
  double precision, allocatable :: basisfunctions(:,:,:,:)
  !> Sounding IDs (if available)
  integer(8), allocatable :: sounding_ids(:,:)

  !> Final statevector and post-uncert for each measurement
  double precision, allocatable :: final_SV(:,:,:)
  double precision, allocatable :: final_SV_uncert(:,:,:)
  !> Final SIF value in physical radiance units
  double precision, allocatable :: retrieved_SIF_abs(:,:)
  !> Final SIF uncertainty in physical radiance units
  double precision, allocatable :: retrieved_SIF_abs_uncertainty(:,:)
  !> Final SIF value as a fraction of continuum level radiance
  double precision, allocatable :: retrieved_SIF_rel(:,:)
  !> Final SIF uncertainty as a fraction of continuum level radiance
  double precision, allocatable :: retrieved_SIF_rel_uncertainty(:,:)
  !> Final retrieved CHI2
  double precision, allocatable :: chi2(:,:)

  !> Final radiance at last iteration
  double precision, allocatable :: final_radiance(:,:,:)
  !> Measured radiance
  double precision, allocatable :: measured_radiance(:,:,:)
  !> Noise-equivalent radiance
  double precision, allocatable :: noise_radiance(:,:,:)

contains

  !> @brief Function to handle all the retrieval logistics, mainly
  !> reading in the needed arrays, and in the big loop call the forward
  !> model / retrieval function.
  !> @param my_instrument Instrument entity
  !>
  subroutine guanter_retrieval(my_instrument)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument

    !> Various HDF file IDs for L1B, basisfunction file, and output file
    integer(hid_t) :: l1b_file_id, basisfunction_file_id, output_file_id

    !> Function name for logging
    character(len=*), parameter :: fname = "guanter_retrieval"
    !> Temporary string for various messages
    character(len=999) :: tmp_str

    !> Placeholder to grab a basisfunction
    double precision, allocatable :: tmp_basisfunctions(:)
    !> 2D HDF dimension array
    integer(hsize_t) :: out_dims2d(2)
    !> 3D HDF dimension array
    integer(hsize_t) :: out_dims3d(3)
    !> Number of pixels read from the basisfunction HDF file
    integer(hsize_t), allocatable :: num_pixels(:)

    !> CPU time stamps and mean duration for performance analysis
    double precision :: cpu_time_start, cpu_time_stop, mean_duration

    !> Indices for frames, footprints, microwindows
    integer :: i_fr, i_fp, i_win
    !> Variables to hold total number of frames, footprints, and the product
    integer :: N_fr, N_fp, N_bands, num_total_soundings
    !> Loop counters
    integer :: cnt, band, fp
    integer :: i,j

    !> Character for HDF dataset name
    character(len=999) :: dset_name
    !> HDF dataset and group IDs
    integer(hid_t) :: dset_id, sif_result_gid
    !> HDF group name
    character(len=999) :: group_name
    !> HDF error variable
    integer :: hdferr


    ! Grab the HDF IDs for easy access
    l1b_file_id = MCS%input%l1b_file_id
    output_file_id = MCS%output%output_file_id

    !
    N_fp = MCS%general%N_fp
    N_fr = MCS%general%N_frame
    N_bands = MCS%general%N_bands

    ! And also create the result group in the output file
    call h5gcreate_f(output_file_id, "linear_fluorescence_results", &
         sif_result_gid, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not create group: linear_fluorescence_results")

    ! Loop over all potential user-defined retrieval windows. We are reading in all
    ! basisfunctions for all retrieval windows already here.
    do i_win=1, MAX_WINDOWS

       if (MCS%window(i_win)%used .eqv. .false.) then
          cycle ! This window is not used!
       end if

       ! Open up the basisfunction file
       call h5fopen_f(MCS%window(i_win)%basisfunction_file%chars(), &
            H5F_ACC_RDONLY_F, basisfunction_file_id, hdferr)
       call check_hdf_error(hdferr, fname, "Error opening HDF file: " &
            // trim(MCS%window(i_win)%basisfunction_file%chars()))

       if (.not. allocated(basisfunctions)) then
          ! And then start reading values into our own arrays
          ! Get the array dimensions by inquiring the first one
          call get_HDF5_dset_dims(basisfunction_file_id, "/BasisFunction_SV1_FP1", num_pixels)
          ! Allocate the basisfunction array
          allocate(basisfunctions(MAX_WINDOWS, num_pixels(1), &
               MCS%general%N_fp, MCS%algorithm%n_basisfunctions))
          ! Allocate the temporary array to hold one single basisfunction
          allocate(tmp_basisfunctions(num_pixels(1)))
       end if

       ! And read them in, one by one
       do i=1, N_fp
          do j=1, MCS%algorithm%n_basisfunctions

             ! For now, I've decided the naming pattern for the basisfunctions
             ! are BasisFunction_SV{X}_FP{Y}, for basisfunction number X and footprint Y
             write(dset_name, '(A,G0.1,A,G0.1)') "/BasisFunction_SV", j, "_FP", i

             call logger%debug(fname, "Looking for " // trim(dset_name))
             call h5dopen_f(basisfunction_file_id, dset_name, dset_id, hdferr)
             call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

             call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_basisfunctions, num_pixels, hdferr)
             call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

             call logger%debug(fname, "Read in " // trim(dset_name))

             ! Copy over data to basisfunctions array. Note!! This might have a
             ! good number of NaNs in there. Does not mean the file is bad!
             basisfunctions(i_win, :, i, j) = tmp_basisfunctions(:)

          end do
       end do

       call h5fclose_f(basisfunction_file_id, hdferr)
       call check_hdf_error(hdferr, fname, "Could not close basisfunction file " &
            // trim(MCS%window(i_win)%basisfunction_file%chars()))

    end do

    ! Read-in of dispersion and noise coefficients (THIS IS INSTRUMENT SPECIFIC!)
    select type(my_instrument)
    type is (oco2_instrument)

       ! Read dispersion coefficients and create dispersion array
       call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
       ! Grab the SNR coefficients for noise calculations
       call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
       ! Read in the sounding id's
       call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)

       allocate(dispersion(maxval(MCS%general%N_spec),MCS%general%N_fp, MCS%general%N_bands))
       do band=1, N_bands
          do fp=1, N_fp
             call my_instrument%calculate_dispersion(dispersion_coefs(:,fp,band), dispersion(:,fp,band), band, fp)
          end do
       end do

    end select

    ! Allocate the (module-wide) result arrays here
    num_total_soundings = MCS%general%N_soundings

    allocate(retrieved_SIF_abs(N_fp, N_fr))
    allocate(retrieved_SIF_rel(N_fp, N_fr))

    allocate(retrieved_SIF_abs_uncertainty(N_fp, N_fr))
    allocate(retrieved_SIF_rel_uncertainty(N_fp, N_fr))

    allocate(chi2(N_fp, N_fr))
    allocate(final_SV(N_fp, N_fr, MCS%algorithm%n_basisfunctions + 1))
    allocate(final_SV_uncert(N_fp, N_fr, MCS%algorithm%n_basisfunctions + 1))


    ! Allocate radiance containers if requested by user
    if (MCS%output%save_radiances) then
       allocate(final_radiance(size(dispersion, 1), MCS%general%N_fp, N_fr))
       allocate(measured_radiance(size(dispersion, 1), MCS%general%N_fp, N_fr))
       allocate(noise_radiance(size(dispersion, 1), MCS%general%N_fp, N_fr))
    end if


    !! Loop through all frames and footprints, and perform the retrieval
    !! The outermose loop simply does it for the various microwindows

    do i_win=1, MAX_WINDOWS

       if (MCS%window(i_win)%used .eqv. .false.) then
          cycle ! This window is not used!
       end if


       call logger%info(fname, "Processing window: " // trim(MCS%window(i_win)%name))

       ! Fill result containers with NaNs
       retrieved_SIF_abs(:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       retrieved_SIF_rel(:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       retrieved_SIF_abs_uncertainty(:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       retrieved_SIF_rel_uncertainty(:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

       final_SV(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       final_SV_uncert(:,:,:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       chi2(:,:) =IEEE_VALUE(1D0, IEEE_QUIET_NAN)

       if (MCS%output%save_radiances) then
          final_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          measured_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
          noise_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
       end if

       cnt = 0
       do i_fr=1, N_fr
          do i_fp=1, N_fp

             write(tmp_str,'(A,I7.1,A,I1.1)') "Frame: ", i_fr, ", FP: ", i_fp
             call logger%trivia(fname, trim(tmp_str))

             ! Retrieval time!
             call cpu_time(cpu_time_start)
             call guanter_FM(my_instrument, i_fr, i_fp, i_win, MCS%window(i_win)%band)
             call cpu_time(cpu_time_stop)

             if (mod(cnt, 1) == 0) then
                write(tmp_str, '(A, G0.1, A, G0.1, A, F6.2, A, F12.6, A)') &
                     "Frame/FP: ", i_fr, "/", i_fp, " ( ", &
                     dble(100 * dble(cnt) / dble(MCS%general%N_frame * MCS%general%N_fp)), "%) - ", &
                     (cpu_time_stop - cpu_time_start), ' sec.'
                call logger%debug(fname, trim(tmp_str))
             end if

             cnt = cnt + 1
          end do
       end do

       !! Write the results into the output HDF5 file

       call logger%info(fname, "Writing out results.")
       ! And then just write out the datasets / arrays

       group_name = "linear_fluorescence_results/" // trim(MCS%window(i_win)%name)

       call logger%info(fname, "Writing group: " // trim(group_name))
       call h5gcreate_f(output_file_id, group_name, sif_result_gid, hdferr)
       call check_hdf_error(hdferr, fname, "Error. Could not create group: " // trim(group_name))

       do i=1, MCS%algorithm%n_basisfunctions

          out_dims2d(1) = N_fp
          out_dims2d(2) = N_fr

          write(tmp_str, '(A,A,G0.1)') trim(group_name), "/basisfunction_coefficient_", i
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               final_SV(:,:,i), out_dims2d)

          write(tmp_str, '(A,A,G0.1,A)') trim(group_name), &
               "/basisfunction_coefficient_", i, "_uncertainty"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               final_SV_uncert(:,:,i), out_dims2d)
       end do

       out_dims2d = shape(retrieved_SIF_rel)
       write(tmp_str, '(A,A)') trim(group_name), "/SIF_relative"
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            retrieved_SIF_rel, out_dims2d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

       out_dims2d = shape(retrieved_SIF_rel)
       write(tmp_str, '(A,A)') trim(group_name), "/SIF_relative_uncertainty"
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            retrieved_SIF_rel_uncertainty(:,:), &
            out_dims2d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

       out_dims2d = shape(retrieved_SIF_abs)
       write(tmp_str, '(A,A)') trim(group_name), "/SIF_absolute"
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            retrieved_SIF_abs, out_dims2d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

       out_dims2d = shape(retrieved_SIF_rel)
       write(tmp_str, '(A,A)') trim(group_name), "/SIF_absolute_uncertainty"
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            retrieved_SIF_abs_uncertainty(:,:), &
            out_dims2d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

       out_dims2d = shape(chi2)
       write(tmp_str, '(A,A)') trim(group_name), "/retrieved_chi2"
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            chi2, out_dims2d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

       if (MCS%output%save_radiances) then
          out_dims3d = shape(final_radiance)
          write(tmp_str, '(A,A)') trim(group_name), "/modelled_radiance"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               final_radiance, out_dims3d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

          out_dims3d = shape(measured_radiance)
          write(tmp_str, '(A,A)') trim(group_name), "/measured_radiance"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               measured_radiance, out_dims3d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))

          out_dims3d = shape(noise_radiance)
          write(tmp_str, '(A,A)') trim(group_name), "/noise_radiance"
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               noise_radiance, out_dims3d, IEEE_VALUE(1D0, IEEE_QUIET_NAN))
       end if
       call logger%info(fname, "Finished writing out results.")
    end do ! End retrieval window loop
  end subroutine guanter_retrieval


  !> @brief Performs the retrival using the Guanter method
  !> @param my_instrument Instrument instance
  !> @param i_fr Frame counter
  !> @param i_fp Footprint counter
  !> @param i_win Retrieval window counter
  !> @param band Current band
  subroutine guanter_FM(my_instrument, i_fr, i_fp, i_win, band)

    implicit none

    class(generic_instrument), intent(in) :: my_instrument
    integer, intent(in) :: i_fr
    integer, intent(in) :: i_fp
    integer, intent(in) :: i_win
    integer, intent(in) :: band

    ! Lowest and highest L1B pixel index
    integer :: l1b_wl_idx_min, l1b_wl_idx_max
    ! Loop counters
    integer :: i, j, l, cnt
    ! L1B HDF file ID
    integer(hid_t) :: l1b_file_id

    ! Function name for logging
    character(len=*), parameter :: fname = "guanter_FM"
    ! Temporary string for misc. stuff
    character(len=999) :: tmp_str

    ! This is where most of the juice is. The retrieval concept is fairly
    ! simple, but there are still some tasks to do to get it done. i_fr and
    ! i_fp tell us the location of the sounding we will process here.

    ! The measured radiance coming straight from the L1B file, and then to
    ! be further processed (slope-corrected)
    double precision, allocatable :: radiance_work(:), noise_work(:), rad_conv(:)
    double precision, allocatable :: radiance_l1b(:)

    ! Continuum level radiance
    double precision :: cont_rad
    ! The intercept from the fit used to normalize the radiances
    ! which is also needed to re-scale the retrieved sif back up.
    double precision :: fit_intercept
    ! Are the radiances that we process sensible?
    logical :: radiance_OK
    ! Was matrix inversion successful?
    logical :: success_inv_mat
    ! Was slope normalization successful?
    logical :: success_slope
    ! Temp. placeholder for CHI2 calculation
    double precision :: tmp_chi2
    ! Variables to hold the Bayesian information content
    double precision :: BIC_all, this_BIC, min_BIC
    ! Variables to hold the position of the lowest BIC
    integer :: BIC_pos
    ! BIC's for the current or last iteration
    double precision, allocatable :: BIC_last_iteration(:), BIC_this_iteration(:)
    ! Number of eliminated PCs
    integer :: N_eliminated
    ! Arrays for holding the indices of which basisfunctions are used or dropped
    logical, allocatable :: used_basisfunctions(:), dropped_basisfunctions(:)
    ! Number of statevector elements
    integer :: N_sv
    ! Number of spectral points
    integer :: N_spec
    ! Number of available basisfunctions
    integer :: N_basisfunctions
    ! Have we finished the backwards elimination?
    logical :: finished_elimination

    ! File unit
    integer :: funit

    ! Retrieval quantities

    ! Solution state vector
    double precision, allocatable :: xhat(:)
    ! Jacobian, inverse noise covariance matrix
    double precision, allocatable :: K(:,:), Se_inv(:,:)
    ! Regular and inverse of the posterior covariance matrix
    double precision, allocatable :: Shat(:,:), Shat_inv(:,:)
    ! Temporary matrices for calculations
    double precision, allocatable :: m_tmp1(:,:), m_tmp2(:,:)

    ! Grab the L1B HDF file handler
    l1b_file_id = MCS%input%l1b_file_id

    ! Read the spectrum from the L1B HDF file
    select type(my_instrument)
    type is (oco2_instrument)
       call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, &
            band, MCS%general%N_spec(band), radiance_l1b)
    end select

    ! Calculate the continuum radiance for this spectrum (the whole band)
    cont_rad = percentile(radiance_l1b, 99.9d0)

    ! We now have the full L1b spectrum of a given index, we need to grab
    ! the relevant portion of the spectrum.
    l1b_wl_idx_min = 0
    l1b_wl_idx_max = 0

    do i=1, size(dispersion, 1)
       if (dispersion(i, i_fp, band) < MCS%window(i_win)%wl_min) then
          l1b_wl_idx_min = i
       end if
       if (dispersion(i, i_fp, band) < MCS%window(i_win)%wl_max) then
          l1b_wl_idx_max = i
       end if
    end do

    ! If window lower limit is below the first wavelength value, set it
    ! to the beginning (index 1)
    if (l1b_wl_idx_min == 0) then
       l1b_wl_idx_min = 1
    end if

    if (l1b_wl_idx_max < size(dispersion, 1)) then
       l1b_wl_idx_max = l1b_wl_idx_max + 1
    end if

    if (l1b_wl_idx_max > size(dispersion, 1)) then
       l1b_wl_idx_max = size(dispersion, 1)
    end if

    ! Allocate some space for the new work array which we can do the
    ! retrieval with, and copy over the corrsponding part of the spectrum
    allocate(radiance_work(l1b_wl_idx_max - l1b_wl_idx_min + 1))
    allocate(rad_conv(l1b_wl_idx_max - l1b_wl_idx_min + 1))
    allocate(noise_work(l1b_wl_idx_max - l1b_wl_idx_min + 1))

    ! Copy the relevant part of the spectrum to radiance_work
    radiance_work(:) = radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max)
    N_spec = size(radiance_work)

    ! Here we check the radiance
    select type(my_instrument)
    type is (oco2_instrument)
       call my_instrument%check_radiance_valid(l1b_file_id, radiance_work, &
            l1b_wl_idx_min, l1b_wl_idx_max, &
            radiance_OK)
    end select

    if (radiance_OK .eqv. .false.) then
       write(tmp_str, '(A,G0.1,A,G0.1,A)') "Sounding (", i_fr, ",", i_fp, &
            ") has invalid radiances. Skipping."
       call logger%warning(fname, trim(tmp_str))
       return
    end if

    ! Now calculate the noise-equivalent radiances
    select type(my_instrument)
    type is (oco2_instrument)
       call my_instrument%calculate_noise(snr_coefs, radiance_work, &
            noise_work, i_fp, band, &
            l1b_wl_idx_min, l1b_wl_idx_max)
    end select

    ! Do the slope correction for the selected spectrum. The second
    ! parameter sets from which percentage onwards the radiances should
    ! be considered for the fit of the slope. This should remove points from
    ! absorption/solar lines so the fit really goes for the continuum level signal.
    !
    ! 40% seems to work and there is no good reason now to change this.
    call slope_correction(radiance_work, 40.0d0, fit_intercept, success_slope)
    if (.not. success_slope) then
       call logger%error(fname, "Slope normalization failed!")
       return
    end if
    ! And we have our "y" mesurement vector ready!!

    ! Calculate the inverse(!) noise matrix Se_inv
    allocate(Se_inv(N_spec, N_spec))
    ! Since we work in slope-normalized space, we also need to scale our
    ! noise levels.
    noise_work(:) = noise_work(:) / fit_intercept

    ! Construct the inverse noise covariance matrix.
    Se_inv(:,:) = 0.0d0
    do i=1, size(radiance_work)
       Se_inv(i,i) = 1 / (noise_work(i) ** 2)
    end do

    ! NOTE: HERE, RADIANCE_L1B IS IS IN PHYSICAL UNITS, AS IS NOISE_WORK.
    ! RADIANCE_WORK, HOWEVER, HAS BEEN SLOPE-NORMALIZED AND IS OF UNIT 1.
    ! WE THUS NORMALIZE THE NOISE DATA AS WELL (see lines below)

    ! We use P. Koehler's backward elimination algorithm here. Retrieval
    ! solutions (xhat) are calculated for all cases of removing one SV
    ! element (basisfunction, not the SIF), and then we assess the BIC to
    ! see if any of the cases with a removed one is lowest. This is then
    ! the new baseline, and the process is repeated for removing another
    ! basisfunction, until no further improvement is seen.

    N_basisfunctions = MCS%algorithm%n_basisfunctions
    N_eliminated = 0

    allocate(BIC_this_iteration(N_basisfunctions))
    allocate(dropped_basisfunctions(N_basisfunctions))
    dropped_basisfunctions(:) = .false.

    finished_elimination = .true.

    ! Do initial retrieval, using all supplied basisfunctions
    !call calculate_xhat(basisfunctions(i_win, l1b_wl_idx_min:l1b_wl_idx_max, i_fp, :), &
    !     radiance_work, noise_work, Se_inv, &
    !     xhat, Shat, tmp_chi2, BIC_all)
    !deallocate(xhat)
    !deallocate(Shat)

    min_BIC = BIC_all
    !write(*,*) "Initial BIC: ", BIC_all
    !write(*,*) "Initial CHI2: ", tmp_chi2

    N_eliminated = 0
    do while (.not. finished_elimination)

       if (count(dropped_basisfunctions .eqv. .true.) == (N_basisfunctions - 1)) then
          exit
       end if

       BIC_this_iteration(:) = 0.0d0
       do i=2, N_basisfunctions

          allocate(used_basisfunctions(N_basisfunctions))

          used_basisfunctions(:) = .true.
          used_basisfunctions(i) = .false.

          do l=2, size(used_basisfunctions)
             if (dropped_basisfunctions(l)) used_basisfunctions(l) = .false.
          end do

          if (dropped_basisfunctions(i)) then
             deallocate(used_basisfunctions)
             cycle
          end if

!!$       N_sv = count(used_basisfunctions .eqv. .true.) + 1
!!$
!!$       write(*,*) "Iteration", i
!!$
          call calculate_xhat(basisfunctions(i_win, l1b_wl_idx_min:l1b_wl_idx_max, i_fp, :), &
               radiance_work, noise_work, Se_inv, &
               xhat, rad_conv, &
               Shat, tmp_chi2, this_BIC, used_basisfunctions)
          deallocate(Shat)

          BIC_this_iteration(i) = this_BIC
!!$
          !write(*,*) used_basisfunctions
          !write(*,*) "BIC: ", this_BIC,  min_BIC
          !write(*,*) "CHI2: ", tmp_chi2
!!$
          deallocate(used_basisfunctions, xhat)
       end do

       !write(*,*) "Iterations: ", i
       BIC_pos = minloc(BIC_this_iteration, dim=1)
       !write(*,*) "MIN BIC POS: ", BIC_pos

       ! We want the BIC to improve by some value
       if ((BIC_this_iteration(BIC_pos) - min_BIC) > -5.0d0) then
          finished_elimination = .true.
       else
          !write(*, *) "BIC is better by: ", BIC_this_iteration(BIC_pos) - min_BIC
          min_BIC = BIC_this_iteration(BIC_pos)
          dropped_basisfunctions(BIC_pos) = .true.
       end if

       !if (count(dropped_basisfunctions .eqv. .true.) > 5) finished_elimination = .true.

    end do

    allocate(used_basisfunctions(N_basisfunctions))

    used_basisfunctions(:) = .true.

    do l=1, size(used_basisfunctions)
       if (dropped_basisfunctions(l)) used_basisfunctions(l) = .false.
    end do

    ! Calculate the maximum posterior statevector for the purely linear case
    call calculate_xhat(basisfunctions(i_win, l1b_wl_idx_min:l1b_wl_idx_max, i_fp, :), &
         radiance_work, noise_work, Se_inv, &
         xhat, rad_conv, &
         Shat, tmp_chi2, this_BIC, used_basisfunctions)

    ! Store the calculated CHI2 value
    chi2(i_fp, i_fr) = tmp_chi2
    !write(*,*) "FINAL BIC: ", this_BIC
    !write(*,*) "FINAL CHI2: ", tmp_chi2

    N_sv = size(xhat)
    ! Report the state vector posteriori uncertainty
    do i=1, N_sv
       final_SV_uncert(i_fp, i_fr, i) = sqrt(Shat(i,i))
    end do

    ! Plug the final value into the container
    final_SV(i_fp, i_fr, :) = xhat(:)

    ! Re-scale the retrieved relative SIF to physical units
    retrieved_SIF_abs(i_fp, i_fr) = xhat(N_sv) * fit_intercept
    !write(tmp_str, '(A, E12.5)') "SIF abs: ", xhat(N_sv) * fit_intercept
    !call logger%trivia(fname, trim(tmp_str))

    retrieved_SIF_rel(i_fp, i_fr) = xhat(N_sv)
    retrieved_SIF_rel_uncertainty(i_fp, i_fr) = sqrt(Shat(N_sv, N_sv))
    retrieved_SIF_abs_uncertainty(i_fp, i_fr) = retrieved_SIF_rel_uncertainty(i_fp, i_fr) * fit_intercept
    !write(tmp_str, '(A, F12.3,A)') "SIF rel: ", xhat(N_sv) * 100, "%"
    !call logger%trivia(fname, trim(tmp_str))

    if (MCS%output%save_radiances) then
       ! If requested, store the radiances and the noise values
       final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = rad_conv(:)
       measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_work(:)
       noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)
    end if

    write(tmp_str, '(A, I3.1)') "N_sv: ", N_sv
    !call logger%trivia(fname, trim(tmp_str))
    write(tmp_str, '(A,F7.3)') "Chi2: ", tmp_chi2
    !call logger%trivia(fname, trim(tmp_str))

    !deallocate(Shat)
  end subroutine guanter_FM


  !> @brief Function for performing the actual retrieval and calculating xhat solution
  !> @param basisfunctions Array that contains the (pixel, j) basisfunctions
  !> @param radiance_work The normalized and slope-corrected measured radiances
  !> @param noise_work The normalized per-pixel noise-equivalent radiance
  !> @param Se_inv Inverse noise error covariance (must correspond to noise_work)
  !> @param x_hat Linear solution (maximum a-posteriori)
  !> @param rad_conv Modelled radiances
  !> @param Shat Posterior convariance matrix
  !> @param chi2 Retrieved CHI2 (no prior)
  !> @param BIC Bayesian information content of this retrieval
  !> @param used_basisfunctions Whicih basisfunctions to use in the retrieval?
  subroutine calculate_xhat(basisfunctions, radiance_work, noise_work, Se_inv, &
       xhat, rad_conv, Shat, chi2, BIC, used_basisfunctions)

    double precision, intent(in) :: basisfunctions(:,:)
    double precision, intent(in) :: radiance_work(:)
    double precision, intent(in) :: noise_work(:)
    double precision, intent(in) :: Se_inv(:,:)
    double precision, intent(inout), allocatable :: xhat(:)
    double precision, intent(inout) :: rad_conv(:)
    double precision, allocatable, intent(inout) :: Shat(:,:)
    double precision, intent(inout) :: chi2
    double precision, intent(inout) :: BIC
    logical, intent(in), optional :: used_basisfunctions(:)

    ! Function name
    character(len=*), parameter :: fname = "calculate_xhat"

    ! Array sizes
    integer :: N_sv
    integer :: N_spec
    ! Jacobian, gain matrix, averaging kernel
    double precision, allocatable :: K(:,:), G(:,:), AK(:,:)
    ! Inverse of posterior covariance matrix
    double precision, allocatable :: Shat_inv(:,:)
    ! Temp. matrices for calculations
    double precision, allocatable :: m_tmp1(:,:), m_tmp2(:,:)
    ! Was matrix inversion successful?
    logical :: success_inv_mat
    ! Loop and counter variables
    integer :: i, j, l, cnt

    ! If requested by the user, we can choose to only use
    ! certain basisfunctions (required by the backward elimination algorithm)
    if (present(used_basisfunctions)) then
       N_sv = count(used_basisfunctions .eqv. .true.) + 1
    else
       N_sv = size(basisfunctions, 2) + 1
    end if

    ! Determine the number of pixels
    N_spec = size(basisfunctions, 1)

    allocate(xhat(N_sv))
    xhat(:) = 0.0d0

    allocate(Shat(N_sv, N_sv))
    Shat(:,:) = 0.0d0

    allocate(Shat_inv(N_sv, N_sv))
    Shat_inv(:,:) = 0.0d0

    allocate(K(N_spec, N_sv))
    K(:,:) = 0.0d0

    allocate(G(N_sv, N_spec))
    G(:,:) = 0.0d0

    allocate(AK(N_sv, N_sv))
    AK(:,:) = 0.0d0

    ! Fill Jacobian with the basisfunctions. This needs to be populated
    ! here as the used basisfunctions might change and thus K will look
    ! different with every retrieval.
    if (present(used_basisfunctions)) then
       cnt = 1
       do i=1, size(basisfunctions, 2)
          if (used_basisfunctions(i)) then
             K(:, cnt) = basisfunctions(:, i)
             cnt = cnt + 1
          end if
       end do
    else
       do i=1, size(basisfunctions, 2)
          K(:, i) = basisfunctions(:, i)
       end do
    end if

    ! SIF Jacobian - REPLACE THIS BY SHAPE IF YOU REALLY WANT TO
    K(:, N_sv) = 1.0d0

    ! Start by calculating posterior covariance for linear case
    ! Shat = (K^T S_e^-1 K)^-1 (Rodgers 2.27 when S_a is large)
    allocate(m_tmp1, mold=Shat)
    m_tmp1 = matmul(matmul(transpose(K), Se_inv), K) ! this is (K^T Se_inv K)^-1
    Shat_inv(:,:) = m_tmp1

    allocate(m_tmp2, mold=m_tmp1)
    call invert_matrix(m_tmp1, m_tmp2, success_inv_mat)
    if (.not. success_inv_mat) then
       call logger%error(fname, "Failed to invert K^T Se^-1 K")
       return
    end if

    ! m_tmp2 is know (K^T Se_inv K)^-1 = Shat
    Shat(:,:) = m_tmp2(:,:)

    ! Calculate xhat!
    xhat = matmul(matmul(matmul(Shat, transpose(K)), Se_inv), radiance_work)

    rad_conv(:) = 0.0d0
    !! Calculate the modeled radiance via xhat
    do l=1, N_sv
       !write(*,*) l, xhat(l), sqrt(Shat(l,l))
       rad_conv(:) = rad_conv(:) + K(:,l) * xhat(l)
    end do

    !G = matmul(matmul(Shat, transpose(K)), Se_inv)
    !AK = matmul(G, K)

    chi2 = calculate_chi2(radiance_work, rad_conv, &
         noise_work, N_spec - N_sv)

    BIC = -2.0d0 * (&
         dot_product(radiance_work - rad_conv, &
         matmul(Se_inv, radiance_work - rad_conv))) &
         + N_sv * log(real(N_spec))

  end subroutine calculate_xhat


  subroutine slope_correction(radiance, perc, intercept, success)

    implicit none
    double precision, intent(inout) :: radiance(:)
    double precision, intent(in) :: perc
    double precision, intent(inout) :: intercept
    logical, intent(inout) :: success

    double precision :: perc_value
    integer :: num_used
    double precision, allocatable :: adata(:,:), bdata(:,:)
    double precision, allocatable :: work(:)
    integer :: dgels_info
    character(len=*), parameter :: fname = "slope_correction"

    integer :: i, cnt
    ! Slope normalization works like this: first, we grab the top X percent
    ! all radiance points (to make sure we mostly have continuum). X ~ 60
    ! After that, we fit a linear function through these points to get
    ! slope and intersect - by which function we then divide the
    ! radiance.

    success = .false.

    ! This is the percentile value from the radiance array
    perc_value = percentile(radiance, perc)
    num_used = count(radiance > perc_value)

    if (num_used <= 2) then
       ! Could not find more than 2 points to be used for slope correction?
       return
    end if


    ! We need a new container for values ABOVE the percentile one
    allocate(adata(num_used, 2))
    allocate(bdata(num_used, 1))
    allocate(work(2 * num_used))

    ! Take out all the radiance values (with pixel positions) that are
    ! larger than the percentile value, and stick them into tmp_rad and tmp_coord
    adata(:, 1) = 1.0 ! For calculating intercepts
    cnt = 1
    do i=1, size(radiance)
       if (radiance(i) > perc_value) then
          bdata(cnt, 1) = radiance(i)
          adata(cnt, 2)= i-1
          cnt = cnt + 1
       end if
    end do

    !! Now we need to perform a linear regression y = ax + b.
    !! Use the handy LAPACK routine DGELS here.

    call DGELS('N', &          ! TRANS
         num_used, &     ! M
         2, &            ! N
         1, &            ! NRHS
         adata, &        ! A
         num_used, &     ! LDA
         bdata, &        ! B
         num_used, &     ! LDB
         work, &         ! WORK
         2*num_used, &   ! LWORK
         dgels_info)

    if (dgels_info /= 0) then
       call logger%fatal(fname, "Error from DGELS.")
       stop 1
    end if

    ! We now have the intercept in bdata(1,1), and the slope in bdata(2,1)
    ! and can thus slope-correct the spectrum.
    do i=1, size(radiance)
       radiance(i) = radiance(i) / (bdata(1,1) + ((i-1) * bdata(2,1)))
    end do

    intercept = bdata(1,1)
    success = .true.

  end subroutine slope_correction



end module guanter_model_mod
