

module physical_model_addon_mod

  ! User modules
  use math_utils_mod
  use statevector_mod
  use Rayleigh_mod
  use scene_mod
  use aerosols_mod
  use file_utils_mod
  use control_mod

  ! Third-party modules
  use stringifor
  use mod_datetime
  use logger_mod, only: logger => master_logger

  ! System modules
  use HDF5
  use ISO_FORTRAN_ENV
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

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
     !> Column-average dry air mixing ratio for retrieved gases (footprint, frame, gas_number)
     double precision, allocatable :: xgas(:,:,:)
     !> Column-average dry air mixing ratio for prior gases! (footprint, frame, gas_number)
     double precision, allocatable :: xgas_prior(:,:,:)
     !> Pressure levels (footprint, frame, pressure level)
     double precision, allocatable :: pressure_levels(:,:,:)
     !> Prior gas VMRs per level (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: vmr_prior(:,:,:,:)
     !> Retrieved gas VMRs per level (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: vmr_retrieved(:,:,:,:)
     !> Pressure weighting functions (footprint, frame, pressure level)
     double precision, allocatable :: pwgts(:,:,:)
     !> Column averaging kernels (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: col_ak(:,:,:,:)
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
     !> Single-sounding retrieval processing time
     double precision, allocatable :: processing_time(:,:)
     !> Number of moles of dry air per m2
     double precision, allocatable :: ndry(:,:)
  end type result_container


  public scene_altitude

contains


  subroutine compute_coef_at_wl( &
       scn, &
       CS_aerosol, &
       SV, &
       wl, &
       n_mom, &
       n_derivs, &
       ray_tau, &
       aer_sca_tau, &
       coef, &
       lcoef, &
       success)

    type(scene), intent(in) :: scn
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    type(statevector), intent(in) :: SV
    double precision, intent(in) :: wl
    integer, intent(in) :: n_mom ! 1 for scalar, 6 for vector
    integer, intent(in) :: n_derivs
    double precision, intent(in) :: ray_tau(:)
    double precision, intent(in) :: aer_sca_tau(:,:)
    double precision, allocatable, intent(inout) :: coef(:,:,:)
    double precision, allocatable, intent(inout) :: lcoef(:,:,:,:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "compute_coef_at_wl"

    integer :: n_aer
    integer :: n_layer
    integer :: n_pfmom
    integer :: a, l, p, i
    integer :: aer_idx
    integer :: aer_sv_idx
    integer :: l_aero_idx

    double precision, allocatable :: aer_fac(:)

    double precision :: fac
    double precision :: denom
    double precision :: this_aero_height
    double precision, allocatable :: ray_coef(:,:)
    double precision, allocatable :: aerpmom(:,:,:,:)

    success = .false.

    n_layer = scn%num_active_levels - 1
    n_aer = scn%num_aerosols
    ! Set number of phase function coefficients to either 3 (Rayleigh only)
    ! or whatever number we need if we have aerosols
    n_pfmom = max(3, scn%max_pfmom)

    allocate(ray_coef(3, n_mom))
    allocate(aerpmom(n_aer, n_pfmom, n_mom, n_layer))
    allocate(coef(n_pfmom, n_mom, n_layer))
    allocate(lcoef(n_pfmom, n_mom, n_derivs, n_layer))

    coef(:,:,:) = 0.0d0
    lcoef(:,:,:,:) = 0.0d0
    ray_coef(:,:) = 0.0d0
    aerpmom(:,:,:,:) = 0.0d0

    call calculate_rayleigh_scatt_matrix(&
         calculate_rayleigh_depolf(wl), ray_coef(:,:))

    ! Remember, the total phasefunction coefficients are calculated as follows
    ! beta = (beta_aer * tau_aer_sca + beta_ray * tau_ray) / (tau_aer_sca + tau_ray)

    ! Loop over every type of aerosol used
    do l = 1, n_layer

       ! The denominator is the sum of aerosol and Rayleigh scattering optical depth
       denom = ray_tau(l)
       if (n_aer > 0) denom = denom + sum(aer_sca_tau(l, :))

       do a = 1, n_aer

          aer_idx = scn%op%aer_mcs_map(a)

          ! Calculate the "interpolation" factor needed for coef interpolation
          ! Using this here makes the subroutine a bit more general and useful
          fac = (wl - CS_aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_l(a))) &
               / (CS_aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_r(a)) &
               - CS_aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_l(a)) )

          ! This gives us the phase function moments from the mom file, but interpolated
          ! AT the wavelength requested by the user. If this is used to calculate the
          ! total moments at the wavelengths that are given by the mom file itself, "fac"
          ! should be either 0 or 1, and essentially just use the numbers as they are
          ! stored in the file.

          do p = 1, n_mom
             aerpmom(a,1:CS_aerosol(aer_idx)%max_n_coef,p,l) = &
               (1.0d0 - fac) * CS_aerosol(aer_idx)%coef(:,p,scn%op%aer_wl_idx_l(a)) &
               + fac * CS_aerosol(aer_idx)%coef(:,p,scn%op%aer_wl_idx_r(a))

             ! Add aerosol contributions here
             coef(:, p, l) = coef(:, p, l) + aerpmom(a, :, p, l) * aer_sca_tau(l, a)

          end do

       end do ! end aerosol type loop
       !write(*,*) "debug info"
       !write(*,*) "l, fac, ray, denom: ", l, fac, ray_tau(l), denom
       !write(*,*) "aer sca: ", aer_sca_tau(l, :)
       !write(*,*) "coef: ", coef(1:3, :, l) !, aerpmom(1, 1:3, 1, l)


       ! And add Rayleigh contributions here (after partial aerosol sum)
       coef(1:3, :, l) = coef(1:3, :, l) + ray_coef(:, :) * ray_tau(l)
       ! and divide the entire layer-moments by the denominator
       coef(:, :, l) = coef(:, :, l) / denom

       ! Now that we have beta, we can 'simply' calculate the derivative inputs
       ! needed by the RT model(s).

       ! We start counting the position of aerosol indicies here
       l_aero_idx = 2 * n_layer
       if (SV%num_albedo > 0) then
          l_aero_idx = l_aero_idx + 1
       end if

       ! Aerosol AODs
       do i = 1, SV%num_aerosol_aod

          ! The position of the AOD derivative inputs must essentially
          ! match the SV structure, and we are generally using this ordering
          ! TODO: this is a bit hacky, would be nice to have some global
          ! dictionary where one can look these up

          l_aero_idx = l_aero_idx + i
          ! Which aerosol belongs to SV index 'i'?
          aer_sv_idx = SV%aerosol_aod_idx_lookup(i)
          ! What is the corresponding aerosol in the MCS?
          aer_idx = scn%op%aer_mcs_map(aer_sv_idx)

          ! And calculate dBeta/dAOD
          lcoef(:,:,l_aero_idx,l) = aer_sca_tau(l, aer_sv_idx) / scn%op%reference_aod(aer_sv_idx) * &
               (aerpmom(aer_sv_idx,:,:,l) - coef(:,:,l)) / (aer_sca_tau(l, aer_sv_idx) + ray_tau(l))
       end do

       ! Aerosol heights
       do i = 1, SV%num_aerosol_height

          allocate(aer_fac(n_layer))

          ! The position of the AOD derivative inputs must essentially
          ! match the SV structure, and we are generally using this ordering
          ! TODO: this is a bit hacky, would be nice to have some global
          ! dictionary where one can look these up

          ! Keep using the l_aero_idx
          l_aero_idx =  l_aero_idx + i
          ! Which aerosol belongs to SV index 'i'?
          aer_sv_idx = SV%aerosol_height_idx_lookup(i)
          ! What is the corresponding aerosol in the MCS?
          aer_idx = scn%op%aer_mcs_map(aer_sv_idx)

          this_aero_height = exp(SV%svsv(SV%idx_aerosol_height(i))) * scn%atm%psurf

          call calculate_aero_height_factors( &
               scn%atm%p_layers(1:n_layer), &
               this_aero_height, &
               CS_aerosol(scn%op%aer_mcs_map(aer_sv_idx))%default_width, &
               aer_fac)

          ! And calculate dBeta/dAerosolHeight for layer l
          lcoef(:,:,l_aero_idx,l) = aer_sca_tau(l, aer_sv_idx) * aer_fac(l) * &
               (aerpmom(aer_sv_idx,:,:,l) - coef(:,:,l)) / (aer_sca_tau(l, aer_sv_idx) + ray_tau(l))

          deallocate(aer_fac)

       end do

    end do ! Layer loop

    if (any(ieee_is_nan(coef))) then
       call logger%error(fname, "I found NaN(s) in the COEF calculation.")

       read(*,*)
       do a = 1, size(coef, 1)
          do l = 1, size(coef, 2)
             do p = 1, size(coef, 3)

                if (ieee_is_nan(coef(a,l,p))) then
                   write(*,*) "NaN at: ",  a, l, p
                end if
             end do
          end do
       end do

       read(*,*)
       return
    end if

    if (any(ieee_is_nan(lcoef))) then
       call logger%error(fname, "I found NaN(s) in the LCOEF calculation.")
       return
    end if

    success = .true.

  end subroutine compute_coef_at_wl



  subroutine calculate_active_levels(scn)

    type(scene), intent(inout) :: scn

    ! Value to reduce surface pressure if too close to layer boundary
    double precision, parameter :: psurf_bump = 1e-2
    character(len=99) :: tmp_str
    integer :: i

    ! If surface pressure exceeds atmospheric set-up,
    ! set it to -1 and return.
    if (scn%atm%psurf > scn%atm%p(size(scn%atm%p))) then
       scn%num_active_levels = -1
       return
    end if

    ! Counting from TOA down to BOA, the layer for which
    ! psurf is located in, is the last layer - hence the
    ! lower boundary of that layer is considered the last
    ! level and thus sets the number of active levels.

    do i = 1, scn%num_levels - 1
       if ((scn%atm%psurf > scn%atm%p(i)) .and. &
            (scn%atm%psurf < scn%atm%p(i+1))) then

          scn%num_active_levels = i + 1
          exit
       end if

       ! Fringe case:
       ! In the event that the surface pressure is too close to the
       ! as a certain pressure level, it will cause problems later on.
       ! So if we find the surface pressure is too close to a pressure
       ! level, we reduce the surface pressure by a small amount to
       ! avoid the problem.

       if (abs(scn%atm%psurf - scn%atm%p(i+1)) < 1d-5) then

          write(tmp_str, '(A, F14.4, A, F14.4)') "Surface pressure bumped up: ", &
          scn%atm%psurf, " -> ", scn%atm%psurf - psurf_bump

          call logger%debug("calculate_active_levels", trim(tmp_str))

          scn%atm%psurf = scn%atm%psurf - psurf_bump
          scn%num_active_levels = i + 1

          exit
       end if

    end do

    ! Allocate other objects in the scene which depend on the number
    ! of active levels / layers
    if (allocated(scn%atm%pwgts)) deallocate(scn%atm%pwgts)
    allocate(scn%atm%pwgts(scn%num_active_levels))

    if (allocated(scn%atm%pwgts_layers)) deallocate(scn%atm%pwgts_layers)
    allocate(scn%atm%pwgts_layers(scn%num_active_levels + 1))

  end subroutine calculate_active_levels



  !> @brief Reads in the atmosphere file, which contains column-based profiles.
  !>
  !> @param filename Path to the atmosphere file
  !> @param Array of 'strings' containing names of required gases
  !> @param Atmosphere object that will be populated by this function
  !>
  !> @details We supply here a filename, a list of gas strings and the
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

    idx_p = -1
    idx_t = -1

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
          if (allocated(atm%sh_layers)) deallocate(atm%sh_layers)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)
          if (allocated(atm%altitude_levels)) deallocate(atm%altitude_levels)
          if (allocated(atm%altitude_layers)) deallocate(atm%altitude_layers)
          if (allocated(atm%grav)) deallocate(atm%grav)
          if (allocated(atm%grav_layers)) deallocate(atm%grav_layers)
          if (allocated(atm%ndry)) deallocate(atm%ndry)

          ! Allocate according to the file structure
          atm%num_levels = level_count
          atm%num_gases = num_gases

          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))
          allocate(atm%sh_layers(level_count - 1))
          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))
          allocate(atm%altitude_levels(level_count))
          allocate(atm%altitude_layers(level_count - 1))
          allocate(atm%grav(level_count))
          allocate(atm%grav_layers(level_count - 1))
          allocate(atm%ndry(level_count))

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
                   atm%gas_index(i) = j
                   atm%gas_names(i) = split_string(j)%lower()

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

  !> @brief Creates / allocates the "results" container to hold all retrieval results
  !>
  !> @param results Result container
  !> @param num_frames Number of frames
  !> @param num_fp Number of footprints
  !> @param num_SV Number of state vector elements
  !> @param num_gas Number of retrieved (!) gases
  subroutine create_result_container(results, num_frames, num_fp, num_SV, num_gas, num_level)
    implicit none
    type(result_container), intent(inout) :: results
    integer, intent(in) :: num_frames, num_fp, num_SV, num_gas, num_level

    allocate(results%sv_names(num_SV))

    allocate(results%sv_retrieved(num_fp, num_frames, num_SV))
    allocate(results%sv_prior(num_fp, num_frames, num_SV))
    allocate(results%sv_uncertainty(num_fp, num_frames, num_SV))
    allocate(results%xgas(num_fp, num_frames, num_gas))
    allocate(results%xgas_prior(num_fp, num_frames, num_gas))
    allocate(results%pwgts(num_fp, num_frames, num_level))
    allocate(results%col_ak(num_fp, num_frames, num_gas, num_level))
    allocate(results%pressure_levels(num_fp, num_frames, num_level))
    allocate(results%vmr_prior(num_fp, num_frames, num_gas, num_level))
    allocate(results%vmr_retrieved(num_fp, num_frames, num_gas, num_level))
    allocate(results%chi2(num_fp, num_frames))
    allocate(results%residual_rms(num_fp, num_frames))
    allocate(results%dsigma_sq(num_fp, num_frames))
    allocate(results%SNR(num_fp, num_frames))
    allocate(results%SNR_std(num_fp, num_frames))
    allocate(results%continuum(num_fp, num_frames))
    allocate(results%processing_time(num_fp, num_frames))

    allocate(results%num_iterations(num_fp, num_frames))
    allocate(results%converged(num_fp, num_frames))

    allocate(results%ndry(num_fp, num_frames))

    results%sv_names = "NONE"

    ! This might cause problems for some, but I find it convenient
    ! to set 'unused' fields to NaNs. Remember this only works for
    ! reals/double precision, but not for integers.
    results%sv_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_uncertainty = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%pwgts = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%vmr_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%vmr_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%pressure_levels = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%col_ak = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%chi2 = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%residual_rms = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%dsigma_sq = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR_std = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%continuum = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%processing_time = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%ndry = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

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
    deallocate(results%xgas_prior)
    deallocate(results%pwgts)
    deallocate(results%vmr_prior)
    deallocate(results%vmr_retrieved)
    deallocate(results%pressure_levels)
    deallocate(results%col_ak)
    deallocate(results%chi2)
    deallocate(results%residual_rms)
    deallocate(results%dsigma_sq)
    deallocate(results%SNR)
    deallocate(results%SNR_std)
    deallocate(results%continuum)
    deallocate(results%num_iterations)
    deallocate(results%converged)
    deallocate(results%ndry)
    deallocate(results%processing_time)

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
  subroutine assign_SV_names_to_result(results, SV, CS_win)
    implicit none

    type(result_container), intent(inout) :: results
    type(statevector), intent(in) :: SV
    type(CS_window_t), intent(in) :: CS_win

    type(string) :: lower_str
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

       do j=1, SV%num_solar_irrad_scale
          if (SV%idx_solar_irrad_scale(j) == i) then
             write(tmp_str, '(A,G0.1)') "solar_irrad_scale_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! SIF name (really only one at this point)
       do j=1, SV%num_sif
          if (SV%idx_sif(j) == i) then
             write(tmp_str, '(A)') "sif_radiance"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! ZLO name
       do j=1, SV%num_zlo
          if (SV%idx_zlo(j) == i) then
             write(tmp_str, '(A)') "zero_level_offset"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       do j=1, SV%num_temp
          if (SV%idx_temp(j) == i) then
             write(tmp_str, '(A)') "temperature_offset"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Solar shift name
       if (SV%idx_solar_shift(1) == i) then
          write(tmp_str, '(A,A)') "solar_shift"
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
             write(tmp_str, '(A)') "surface_pressure"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Dispersion parameter names
       do j=1, SV%num_dispersion
          if (SV%idx_dispersion(j) == i) then
             write(tmp_str, '(A,G0.1)') "dispersion_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! ILS parameter names
       do j=1, SV%num_ils_stretch
          if (SV%idx_ils_stretch(j) == i) then
             write(tmp_str, '(A,G0.1)') "ils_stretch_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved aerosol AOD names
       do j=1, SV%num_aerosol_aod
          if (SV%idx_aerosol_aod(j) == i) then
             lower_str = CS_win%aerosol(sv%aerosol_aod_idx_lookup(j))%lower()
             write(tmp_str, '(A,A)') lower_str%chars(), "_aod"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved aerosol AOD names
       do j=1, SV%num_aerosol_height
          if (SV%idx_aerosol_height(j) == i) then
             lower_str = CS_win%aerosol(sv%aerosol_height_idx_lookup(j))%lower()
             write(tmp_str, '(A,A)') lower_str%chars(), "_height"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved gas names (scalar retrieval only so far)
       k = 1
       do j=1, SV%num_gas
          ! Check if this SV element is a scalar retrieval
          if (SV%idx_gas(j,1) == i) then
             if (CS_win%gas_retrieve_scale(sv%gas_idx_lookup(j))) then
                if (abs(CS_win%gas_retrieve_scale_start(sv%gas_idx_lookup(j), k) + 1.0) < 1e-10) cycle

                lower_str = CS_win%gases(sv%gas_idx_lookup(j))%lower()
                write(tmp_str, '(A)') trim(lower_str%chars() // "_scale_")
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


  !> @brief Writes global result arrays into the output HDF file
  subroutine write_results_into_hdf_output( &
       CS_win, &
       CS_general, &
       CS_output, &
       output_file_id, &
       results, &
       met_psurf, &
       final_radiance, &
       measured_radiance, &
       noise_radiance, &
       wavelength_radiance)

    type(CS_window_t), intent(in) :: CS_win
    type(CS_general_t), intent(in) :: CS_general
    type(CS_output_t), intent(in) :: CS_output

    integer(hid_t), intent(in) :: output_file_id
    type(result_container), intent(in) :: results
    double precision, allocatable, intent(in) :: met_psurf(:,:)
    double precision, allocatable, intent(in) :: final_radiance(:,:,:)
    double precision, allocatable, intent(in) :: measured_radiance(:,:,:)
    double precision, allocatable, intent(in) :: noise_radiance(:,:,:)
    double precision, allocatable, intent(in) :: wavelength_radiance(:,:,:)

    ! Function name
    character(len=*), parameter :: fname = "write_results_into_hdf_output"
    ! HDF error variable
    integer :: hdferr
    ! Group ID for result group
    integer(hid_t) :: result_gid
    ! Group name for results
    character(len=999) :: group_name
    ! Temporary string
    character(len=999) :: tmp_str
    ! Another temporary string
    type(string) :: lower_str
    ! Containers for 1d, 2d, 3d, and 4d output array dimensions
    integer(hsize_t) :: out_dims2d(2)
    integer(hsize_t) :: out_dims3d(3)
    ! How many SV elements do we have?
    integer :: N_SV
    ! Loop variable
    integer :: i
    ! Does this key exist?
    logical :: hdfexists

    N_SV = size(results%sv_prior, 3)

    !---------------------------------------------------------------------
    ! HDF OUTPUT
    ! Here, we write out the various arrays into HDF datasets
    !---------------------------------------------------------------------
    !
    ! The general idea is to
    ! a) Create a string that corresponds to the output HDF dataset, like
    !    RetrievalResults/physical/ch4/reduced_chi_squared_gbg
    ! b) (optional) re-assign the dataset dimensions needed to write the
    !    arrays into the HDF file using write_*_hdf_dataset
    ! c) Call the write_*_hdf_dataset, usually using the results%** array
    !    to store those values into the HDF file
    !
    !---------------------------------------------------------------------

    ! Create an HDF group for all windows separately
    group_name = "RetrievalResults/physical/" // trim(CS_win%name%chars())
    call h5gcreate_f(output_file_id, trim(group_name), result_gid, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not create group: " &
         // trim(group_name))

    ! Set the dimensions of the arrays for saving them into the HDF file
    out_dims2d(1) = CS_general%n_fp
    out_dims2d(2) = CS_general%n_frame

    ! Save the processing time
    call logger%info(fname, "Writing out: " // trim(group_name) // &
         "/processing_time_" // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/processing_time_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%processing_time(:,:), out_dims2d, -9999.99d0)

    ! Writing out the prior surface pressure, but obviously only if allocated
    if (allocated(met_psurf)) then
       call logger%info(fname, "Writing out: " // trim(group_name) // &
            "/surface_pressure_apriori_" // CS_general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/met_surface_pressure_", &
            CS_general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), met_psurf(:,:), out_dims2d, -9999.99d0)
    end if

    ! Save the prior state vectors
    do i=1, N_SV
       write(tmp_str, '(A,A,A,A,A)') trim(group_name) , "/", &
            results%sv_names(i)%chars() , "_apriori_", CS_general%code_name
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%sv_prior(:,:,i), out_dims2d, -9999.99d0)
    end do

    ! Save the retrieved state vectors
    do i=1, N_SV
       write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", results%sv_names(i)%chars(), &
            "_", CS_general%code_name
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%sv_retrieved(:,:,i), out_dims2d, -9999.99d0)
    end do

    ! Save the retrieved state vector uncertainties
    do i=1, N_SV
       write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", &
            results%sv_names(i)%chars() , "_uncertainty_", CS_general%code_name
       call logger%info(fname, "Writing out: " // trim(tmp_str))
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), &
            results%sv_uncertainty(:,:,i), out_dims2d, -9999.99d0)
    end do


    ! Here we re-set the output dimensions for 3D fields, needed
    ! for column AKs

    out_dims3d(1) = CS_general%n_fp
    out_dims3d(2) = CS_general%n_frame
    out_dims3d(3) = size(results%col_AK, 4)

    do i=1, CS_win%num_gases
       if (CS_win%gas_retrieved(i)) then

          ! Save XGAS for each gas
          lower_str = CS_win%gases(i)%lower()

          write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/x", &
               lower_str%chars(), "_", CS_general%code_name
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%xgas(:,:,i), out_dims2d, -9999.99d0)

          write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/x", &
               lower_str%chars(), "_apriori_", CS_general%code_name
          call logger%info(fname, "Writing out: " // trim(tmp_str))
          call write_DP_hdf_dataset(output_file_id, &
               trim(tmp_str), &
               results%xgas_prior(:,:,i), out_dims2d, -9999.99d0)

          if (CS_output%pressure_weights) then
             write(tmp_str, '(A,A,A,A)') trim(group_name), &
                  "/pressure_weights_", CS_general%code_name
             call h5lexists_f(output_file_id, trim(tmp_str), hdfexists, hdferr)

             if (.not. hdfexists) then
                call logger%info(fname, "Writing out: " // trim(tmp_str))
                call write_DP_hdf_dataset(output_file_id, &
                     trim(tmp_str), &
                     results%pwgts(:,:,:), out_dims3d, -9999.99d0)
             end if
          end if

          if (CS_output%gas_averaging_kernels) then
             ! These variables only need to be saved if the user requested
             ! column averaging kernels. They add quite a bit to the output
             ! file size, hence you want to make sure that you really need
             ! them.

             ! Write pressure levels only once!
             write(tmp_str, '(A,A,A,A)') trim(group_name), "/pressure_levels", &
                  "_", CS_general%code_name
             call h5lexists_f(output_file_id, trim(tmp_str), hdfexists, hdferr)

             if (.not. hdfexists) then
                call logger%info(fname, "Writing out: " // trim(tmp_str))
                call write_DP_hdf_dataset(output_file_id, &
                     trim(tmp_str), &
                     results%pressure_levels(:,:,:), out_dims3d, -9999.99d0)
             end if

             if (CS_output%gas_averaging_kernels) then
                write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/x", &
                     lower_str%chars(), "_column_ak_", CS_general%code_name
                call logger%info(fname, "Writing out: " // trim(tmp_str))
                call write_DP_hdf_dataset(output_file_id, &
                     trim(tmp_str), &
                     results%col_AK(:,:,i,:), out_dims3d, -9999.99d0)
             end if

             write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", &
                  lower_str%chars(), "_profile_apriori_", CS_general%code_name
             call logger%info(fname, "Writing out: " // trim(tmp_str))
             call write_DP_hdf_dataset(output_file_id, &
                  trim(tmp_str), &
                  results%vmr_prior(:,:,i,:), out_dims3d, -9999.99d0)

             write(tmp_str, '(A,A,A,A,A)') trim(group_name), "/", &
                  lower_str%chars(), "_profile_retrieved_", CS_general%code_name
             call logger%info(fname, "Writing out: " // trim(tmp_str))
             call write_DP_hdf_dataset(output_file_id, &
                  trim(tmp_str), &
                  results%vmr_retrieved(:,:,i,:), out_dims3d, -9999.99d0)

          end if

       end if
    end do

    ! Save number of iterations
    call logger%info(fname, "Writing out: " // trim(group_name) // "/num_iterations_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/num_iterations_", CS_general%code_name
    call write_INT_hdf_dataset(output_file_id, &
         trim(tmp_str), results%num_iterations(:,:), out_dims2d, -9999)

    ! Save converged status
    call logger%info(fname, "Writing out: " // trim(group_name) // "/converged_flag_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/converged_flag_", CS_general%code_name
    call write_INT_hdf_dataset(output_file_id, &
         trim(tmp_str), results%converged(:,:), out_dims2d, -9999)

    ! Retrieved CHI2
    call logger%info(fname, "Writing out: " // trim(group_name) // "/reduced_chi_squared_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/reduced_chi_squared_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%chi2(:,:), out_dims2d, -9999.99d0)

    ! Residual RMS
    call logger%info(fname, "Writing out: " // trim(group_name) // "/residual_rms_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/residual_rms_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%residual_rms(:,:), out_dims2d, -9999.99d0)

    ! Dsigma_sq
    call logger%info(fname, "Writing out: " // trim(group_name) // "/final_dsigma_sq_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/final_dsigma_sq_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%dsigma_sq(:,:), out_dims2d, -9999.99d0)

    ! Signal-to-noise ratio (mean)
    call logger%info(fname, "Writing out: " // trim(group_name) // "/snr_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/snr_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%SNR, out_dims2d)

    call logger%info(fname, "Writing out: " // trim(group_name) // "/continuum_level_radiance_" &
         // CS_general%code_name)
    write(tmp_str, '(A,A,A)') trim(group_name), "/continuum_level_radiance_", CS_general%code_name
    call write_DP_hdf_dataset(output_file_id, &
         trim(tmp_str), results%continuum, out_dims2d)

    ! Save the radiances, only on user request (non-default)
    if (CS_output%save_radiances) then

       out_dims3d = shape(final_radiance)

       call logger%info(fname, "Writing out: " // trim(group_name) // "/modelled_radiance_" &
            // CS_general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/modelled_radiance_", CS_general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), final_radiance, out_dims3d)

       out_dims3d = shape(measured_radiance)
       call logger%info(fname, "Writing out: " // trim(group_name) // "/measured_radiance_" &
            // CS_general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/measured_radiance_", CS_general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), measured_radiance, out_dims3d)

       out_dims3d = shape(noise_radiance)
       call logger%info(fname, "Writing out: " // trim(group_name) // "/noise_radiance_" &
            // CS_general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/noise_radiance_", CS_general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), noise_radiance, out_dims3d)

       out_dims3d = shape(wavelength_radiance)
       call logger%info(fname, "Writing out: " // trim(group_name) // "/wavelength_" &
            // CS_general%code_name)
       write(tmp_str, '(A,A,A)') trim(group_name), "/wavelength_", CS_general%code_name
       call write_DP_hdf_dataset(output_file_id, &
            trim(tmp_str), wavelength_radiance, out_dims3d)

    end if

  end subroutine write_results_into_hdf_output

  !> @brief Given a dispersion array "this_dispersion", in window "i_win", this
  !> function calculates the first and last pixel indices as the boundaries in
  !> the detector.
  !> @param this_dispersion Wavelength per detector index
  !> @param i_win Retrieval Window index for MCS
  !> @param l1b_wl_idx_min Lower wavelength pixel index corresponding to
  !> user-defined "wl_min"
  !> @param l1b_wl_idx_min Upper wavelength pixel index corresponding to
  !> user-defined "wl_max"
  subroutine calculate_dispersion_limits(this_dispersion, CS_win, &
       l1b_wl_idx_min, l1b_wl_idx_max)

    implicit none
    double precision, intent(in) :: this_dispersion(:)
    type(CS_window_t), intent(in) :: CS_win
    integer, intent(inout) :: l1b_wl_idx_min, l1b_wl_idx_max

    integer :: i

    l1b_wl_idx_min = 0
    l1b_wl_idx_max = 0

    ! This here grabs the boundaries of the L1b data by simply looping over
    ! the dispersion array and setting the l1b_wl_idx_* accordingly.
    do i=1, size(this_dispersion)
       if (this_dispersion(i) < CS_win%wl_min) then
          l1b_wl_idx_min = i
       end if
       if (this_dispersion(i) < CS_win%wl_max) then
          l1b_wl_idx_max = i
       end if
       if (this_dispersion(i) > CS_win%wl_max) then
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

  !> @brief Calculates the boundaries of the subcolumns
  !> @param SV State vector object
  !> @param gas_idx gas index (from MCS%gas)

  subroutine set_gas_scale_levels( &
       SV, &
       gas_idx, &
       CS_win, &
       atm, &
       psurf, &
       s_start, &
       s_stop, &
       do_gas_jac, &
       success)

    type(statevector), intent(inout) :: SV
    integer, intent(in) :: gas_idx
    type(CS_window_t), intent(in) :: CS_win
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

       if (CS_win%gas_retrieve_scale(gas_idx)) then

          do_gas_jac = .true.
          if (SV%gas_idx_lookup(i) == gas_idx) then

             ! This bit here figures out which level/layer range a
             ! certain scaling factor corresponds to. They are fractions
             ! of surface pressure, so we use 'searchsorted' to find
             ! where they would belong to. We also make sure it can't
             ! go below or above the first/last level.

             s_start(i) = searchsorted_dp((atm%p), &
                  SV%gas_retrieve_scale_start(i) * (psurf), .true.)
             s_start(i) = max(1, s_start(i))

             s_stop(i) = searchsorted_dp((atm%p), &
                  SV%gas_retrieve_scale_stop(i) * (psurf), .true.) + 1
             s_stop(i) = min(size(atm%p), s_stop(i))

             SV%s_start(i) = s_start(i)
             SV%s_stop(i) = s_stop(i)

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
             SV%s_start(i) = SV%s_start(i) + 1
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
          write(tmp_str, '(A,A)') "Scale factor index error for gas ", &
               CS_win%gases(SV%gas_idx_lookup(i))%chars()
          call logger%error(fname, trim(tmp_str))
          return
       end if
    end do

    success = .true.

  end subroutine set_gas_scale_levels



  !> @begin Wrapper to replace prior VMRs with special functions
  subroutine replace_prior_VMR(scn, prior_types)

    type(scene), intent(inout) :: scn
    type(string), intent(in) :: prior_types(:)

    ! Function name
    character(len=*), parameter :: fname = "replace_prior_VMR"
    character(len=999) :: tmp_str

    integer :: i

    ! The prior_types are in order of scn%atm%gas_vmr, so
    ! we can calculate a new prior using prior_type(i) and stick
    ! it into scn%atm%gas_vmr(:,i).

    do i=1, size(prior_types)

       ! Nothing to do if this string is empty
       if (trim(prior_types(i)) == "") cycle

       if (prior_types(i) == "SC4C2018") then
          ! Max Reuter / Oliver Schneising
          ! CH4 profiles as derived from a climatology file
          call logger%debug(fname, "Using Reuter/Schneising SC4C2018 model for CH4.")

       else
          ! If the prior type is not implemented, the user has made a mistake.
          ! Again, we are terminating here immediately, since falling back to
          ! some default behavior is not a good option..

          write(tmp_str, '(A,A)') "Sorry, the following prior VMR function " &
              // "is not implemented: ", prior_types(i)%chars()
          call logger%fatal(fname, trim(tmp_str))
          !stop 1
       end if

    end do

  end subroutine replace_prior_VMR



end module physical_model_addon_mod
