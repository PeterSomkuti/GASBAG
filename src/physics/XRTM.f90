!> @brief Module to house the interfaces and calculations for the XRTM RT model
!> @file XRTM.f90
!> @author Peter Somkuti

module XRTM_mod

  ! User modules
  use control_mod, only: CS_window_t
  use misc_utils_mod
  use math_utils_mod
  use Rayleigh_mod, only : calculate_rayleigh_scatt_matrix, calculate_rayleigh_depolf
  use statevector_mod
  use aerosols_mod
  use scene_mod
  use physical_model_addon_mod
  use PCART_mod

  ! Third-party modules
  use xrtm_int_f90
  use stringifor
  use logger_mod, only: logger => master_logger

  ! System modules
  use OMP_LIB
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_is_finite


  implicit none

  public :: solve_RT_problem_XRTM, calculate_XRTM_radiance

contains


  subroutine calculate_XRTM_ltaugas_inputs(SV, scn, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)


    integer :: i, j, l
    integer :: n_lay

    n_lay = scn%num_active_levels - 1

    do i = 1, n_lay
       do j = 1, n_lay

          if (i /= j) cycle

          do l = 1, size(scn%op%wl)

             ltau(l, SV%idx_wf_taugas(j), j) = 1.0d0
             lomega(l, SV%idx_wf_taugas(j), j) = -scn%op%layer_omega(l, j) / scn%op%layer_tau(l, j)

          end do
       end do
    end do

  end subroutine calculate_XRTM_ltaugas_inputs

  subroutine calculate_XRTM_lgas_inputs(SV, scn, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)

    integer :: idx_start, idx_stop
    integer :: j
    integer :: l
    integer :: n_lay

    n_lay = scn%num_active_levels - 1

    do j = 1, SV%num_gas

       idx_start = SV%s_start(j)
       idx_stop = SV%s_stop(j) - 1

       do l = 1, size(scn%op%wl)

          ltau(l, SV%idx_wf_gas(j), idx_start:idx_stop) = &
               scn%op%gas_tau(l, idx_start:idx_stop, SV%gas_idx_lookup(j)) / &
               SV%svsv(SV%idx_gas(j, 1))

          lomega(l, SV%idx_wf_gas(j), idx_start:idx_stop) = &
               (-scn%op%layer_omega(l, idx_start:idx_stop)) &
               / scn%op%layer_tau(l, idx_start:idx_stop) &
               * ltau(l, SV%idx_wf_gas(j), idx_start:idx_stop)
       end do
    end do

  end subroutine calculate_XRTM_lgas_inputs

  subroutine calculate_XRTM_ltemp_inputs(SV, scn, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)

    integer :: l
    integer :: n_lay

    if (SV%num_temp > 0) then

       n_lay = scn%num_active_levels - 1
       do l = 1, size(scn%op%wl)

          ltau(l, SV%idx_wf_temp(1), 1:n_lay) = &
               sum(scn%op%gas_tau_dtemp(l,1:n_lay,:), dim=2)

          lomega(l, SV%idx_wf_temp(1), 1:n_lay) = &
               (-scn%op%layer_omega(l, 1:n_lay)) / scn%op%layer_tau(l, 1:n_lay) * &
               ltau(l, SV%idx_wf_temp(1), 1:n_lay)

       end do
    end if

  end subroutine calculate_XRTM_ltemp_inputs

  subroutine calculate_XRTM_lpsurf_inputs(SV, scn, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)

    integer :: l
    integer :: n_lay

    if (SV%num_psurf > 0) then

       n_lay = scn%num_active_levels - 1

       do l = 1, size(scn%op%wl)

          ltau(l, SV%idx_wf_psurf(1), 1:n_lay) = &
               (-sum(scn%op%gas_tau_dpsurf(l, 1:n_lay, :), dim=2))

          lomega(l, SV%idx_wf_psurf(1), 1:n_lay) = &
               (-scn%op%layer_omega(l, 1:n_lay)) / scn%op%layer_tau(l, 1:n_lay) * &
               ltau(l, SV%idx_wf_psurf(1), 1:n_lay)

       end do
    end if

  end subroutine calculate_XRTM_lpsurf_inputs

  subroutine calculate_XRTM_lalbedo_inputs(SV, scn, lsurf)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: lsurf(:,:)

    integer :: j

    do j = 1, SV%num_albedo
       lsurf(:, SV%idx_wf_albedo(j)) = &
            (scn%op%wl(:) - scn%op%wl(int(size(scn%op%wl) / 2)))**(dble(j-1))
    end do

  end subroutine calculate_XRTM_lalbedo_inputs

  subroutine calculate_XRTM_laerosol_aod_inputs(SV, scn, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)

    integer :: l
    integer :: j
    integer :: n_lay
    integer :: aer_idx

    n_lay = scn%num_active_levels - 1

    do j = 1, SV%num_aerosol_aod

       aer_idx = SV%aerosol_aod_idx_lookup(j)

       do l = 1, size(scn%op%wl)

          ltau(l, SV%idx_wf_aerosol_aod(j), 1:n_lay) = &
               scn%op%aer_ext_tau(l, 1:n_lay, aer_idx) / scn%op%reference_aod(aer_idx)

          lomega(l, SV%idx_wf_aerosol_aod(j), 1:n_lay) = &
               ( scn%op%aer_sca_tau(l, 1:n_lay, aer_idx) / scn%op%aer_ext_tau(l, 1:n_lay, aer_idx) &
               - scn%op%layer_omega(l, 1:n_lay) &
               ) / scn%op%layer_tau(l, 1:n_lay) &
               * ltau(l, SV%idx_wf_aerosol_aod(j), 1:n_lay)

       end do
    end do

  end subroutine calculate_XRTM_laerosol_aod_inputs

  subroutine calculate_XRTM_laerosol_height_inputs(SV, scn, CS_aerosol, ltau, lomega)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    double precision, intent(inout) :: ltau(:,:,:)
    double precision, intent(inout) :: lomega(:,:,:)

    double precision, allocatable :: aero_fac(:)
    integer :: l
    integer :: j
    integer :: n_lay
    integer :: aer_idx

    n_lay = scn%num_active_levels - 1
    allocate(aero_fac(n_lay))


    do j = 1, SV%num_aerosol_height

       aer_idx = SV%aerosol_height_idx_lookup(j)

       call calculate_aero_height_factors_Gauss( &
            scn%atm%p_layers(:), &
            exp(SV%svsv(SV%idx_aerosol_height(j))) * scn%atm%psurf, &
            CS_aerosol(scn%op%aer_mcs_map(aer_idx))%default_width, &
            aero_fac)

       do l = 1, size(scn%op%wl)

          ltau(l, SV%idx_wf_aerosol_height(j), 1:n_lay) = &
               scn%op%aer_ext_tau(l, 1:n_lay, aer_idx) * aero_fac(:) * scn%atm%psurf

          lomega(l, SV%idx_wf_aerosol_height(j), 1:n_lay) = &
               ( &
               scn%op%aer_sca_tau(l, 1:n_lay, aer_idx) / scn%op%aer_ext_tau(l, 1:n_lay, aer_idx)  &
               - scn%op%layer_omega(l, 1:n_lay) &
               ) / scn%op%layer_tau(l, 1:n_lay) &
               * ltau(l, SV%idx_wf_aerosol_height(j), 1:n_lay)

       end do
    end do

  end subroutine calculate_XRTM_laerosol_height_inputs


  !> @brief "Translates" configuration values into XRTM settings
  !> @param xrtm_options_string Array of strings containing XRTM options
  !> @param xrtm_solvers_string Array of strings containing XRTM solvers
  !> @param do_polarization Run RT in vector mode?
  !> @param xrtm_options Integer bitmask for XRTM options
  !> @param xrtm_solvers Integer bitmask for XRTM solvers
  !> @param xrtm_kernels Integer array specifying which surface kernels
  !> @param success Did any error occur?
  subroutine setup_XRTM( &
       xrtm_options_string, &
       xrtm_solvers_string, &
       do_polarization, &
       keep_coef_constant, &
       max_pfmom, &
       n_streams, &
       xrtm_options, &
       xrtm_solvers, &
       xrtm_kernels, &
       success)

    implicit none
    type(string), allocatable, intent(in) :: xrtm_options_string(:)
    type(string), allocatable, intent(in) :: xrtm_solvers_string(:)
    logical, intent(in) :: do_polarization
    logical, intent(in) :: keep_coef_constant
    integer, intent(in) :: max_pfmom
    integer, intent(in) :: n_streams
    integer, allocatable, intent(inout) :: xrtm_options(:)
    integer, allocatable, intent(inout) :: xrtm_solvers(:)
    integer, intent(inout) :: xrtm_kernels(:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "setup_XRTM"

    type(string) :: tmp_str
    integer :: i
    integer :: num_solvers

    success = .false.

    if (.not. allocated(xrtm_solvers_string)) then
       call logger%fatal(fname, "XRTM solver strings not allocated! Aborting.")
       return
    end if

    num_solvers = size(xrtm_solvers_string)

    allocate(xrtm_options(num_solvers))
    allocate(xrtm_solvers(num_solvers))

    ! Initialize XRTM options and solvers with 0
    xrtm_options(:) = 0
    xrtm_solvers(:) = 0

    ! Surface kernels, for the time being make Lambert only
    xrtm_kernels(1) = XRTM_KERNEL_LAMBERTIAN

    ! Populate and option bit fields, depending on requested options,
    ! for each solver separately.

    do i=1, num_solvers

       ! Do weighting functions and pseudo-spherical geometry
       ! by default. Also, only return upwelling radiances.
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_CALC_DERIVS)
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_PSA)
       !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_NO_AZIMUTHAL)
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_UPWELLING_OUTPUT)
       !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DOWNWELLING_OUTPUT)
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_OUTPUT_AT_LEVELS)
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SOURCE_SOLAR)

       if (keep_coef_constant) then
          ! If we use spectrally constant phase function expansion
          ! coeffs, we can make use of this XRTM option which saves
          ! the phase matrix between XRTM calls.
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SAVE_PHASE_MATS)
       end if

       ! If polarization is requested, we run
       ! the RT models in vector mode
       if (do_polarization) then
          call logger%debug(fname, "Using XRTM in vector mode.")
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_VECTOR)
       else
          call logger%debug(fname, "Using XRTM in scalar mode.")
       end if

       ! Grab local copy of string
       tmp_str = xrtm_solvers_string(i)
       if (tmp_str == "TWO_OS") then

          ! Vijay-like 2OS, gives you SECOND ORDER ONLY
          call logger%debug(fname, "Using XRTM in 2-orders-of-scattering mode")
          xrtm_solvers(i) = XRTM_SOLVER_TWO_OS

       else if (tmp_str == "SINGLE") then

          ! Single scattering only
          call logger%debug(fname, "Using XRTM in single-scattering mode")
          xrtm_solvers(i) = XRTM_SOLVER_SINGLE
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SFI)
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)

       else if (tmp_str == "TWO_STREAM") then

          ! Dedicated 2-Stream solver
          call logger%debug(fname, "Using XRTM in two-stream mode")
          !call logger%debug(fname, "Using N-T TMS correction and Delta-M scaling.")
          xrtm_solvers(i) = XRTM_SOLVER_TWO_STREAM
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SFI)

          if (max_pfmom >= (2 * n_streams + 1)) then
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)
          end if

       else if (tmp_str == "EIG_BVP") then

          ! Quadrature (LIDORT-like)
          call logger%debug(fname, "Using XRTM in discrete-ordinate (EIG_BVP) mode")
          !call logger%debug(fname, "Using N-T TMS correction and Delta-M scaling.")
          xrtm_solvers(i) = XRTM_SOLVER_EIG_BVP
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SFI)

          if (max_pfmom >= (2 * n_streams + 1)) then
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)
          end if

       else if (tmp_str == "MEM_BVP") then

          ! Quadrature (LIDORT-like)
          call logger%debug(fname, "Using XRTM in matrix exponential discrete-ordinate mode")
          xrtm_solvers(i) = XRTM_SOLVER_MEM_BVP

          if (max_pfmom >= (2 * n_streams + 1)) then
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
             xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)
          end if

       else if (tmp_str == "SOS") then

          ! Pade approximation and adding
          call logger%debug(fname, "Using XRTM in successive-orders-of-scattering mode")
          call logger%debug(fname, "-rt_streams- is used to infer scattering order")
          xrtm_solvers(i) = XRTM_SOLVER_SOS

       else if (tmp_str == "PADE_ADD") then

          ! Pade approximation and adding
          call logger%debug(fname, "Using XRTM in Pade-adding mode")
          xrtm_solvers(i) = XRTM_SOLVER_PADE_ADD

       else


          call logger%error(fname, "XRTM solver option is not implemented, and will be ignored: " &
               // tmp_str%chars())

       end if

    end do

    success = .true.

  end subroutine setup_XRTM


  subroutine create_XRTM( &
       xrtm, &
       xrtm_options, &
       xrtm_solvers, &
       max_coef, &
       n_quad, &
       n_stokes, &
       n_derivs, &
       n_layers, &
       n_kernels, &
       n_kernel_quad, &
       xrtm_kernels, &
       success)

    implicit none

    type(xrtm_type), intent(inout) :: xrtm
    integer, intent(in) :: xrtm_options
    integer, intent(in) :: xrtm_solvers
    integer, intent(in) :: max_coef
    integer, intent(in) :: n_quad
    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers
    integer, intent(in) :: n_kernels
    integer, intent(in) :: n_kernel_quad
    integer, intent(in) :: xrtm_kernels(:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "create_XRTM"
    character(len=999) :: tmp_str
    integer :: xrtm_error

    success = .false.

    write(tmp_str, '(A, G0.1)') "XRTM Option bit mask: ", xrtm_options
    call logger%debug(fname, trim(tmp_str))
    write(tmp_str, '(A, G0.1)') "XRTM Solver bit mask: ", xrtm_solvers
    call logger%debug(fname, trim(tmp_str))

    call xrtm_create_f90( &
         xrtm, &  ! XRTM object
         xrtm_options, & ! XRTM option bitmask
         xrtm_solvers, & ! XRTM solver bitmask
         max_coef, & ! Number of phase matrix Legendre coefficients
         n_quad, & ! Number of quadrature streams for some solvers
         n_stokes, & ! Number of stokes coefficients
         n_derivs, & ! Number of derivatives to calculate
         n_layers, & ! Number of layers in model atmosphere
         1, & ! Number of solar zeniths
         n_kernels, & ! Number of surface kernels (needs to match up with xrtm_kernels)
         n_kernel_quad, & ! Number of surface kernel quadratures to use
         xrtm_kernels, & ! Which surface kernels to use
         1, & ! Number of outgoing levels
         1, & ! Number of outgoing zeniths
         xrtm_error)

    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_create_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    endif

    success = .true.

  end subroutine create_XRTM


  subroutine solve_RT_problem_XRTM( &
       CS_win, &
       CS_general, &
       CS_aerosol, &
       SV, &
       scn, &
       n_stokes, &
       radiance_calc_work_hi_stokes, &
       dI_dgas, &
       dI_dsurf, &
       dI_dTemp, &
       dI_dAOD, &
       dI_dAHeight, &
       dI_dpsurf, &
       dI_dtaugas, &
       need_gas_averaging_kernels, &
       xrtm_success)

    type(CS_window_t), intent(in) :: CS_win
    type(CS_general_t), intent(in) :: CS_general
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    type(statevector), intent(inout) :: SV
    type(scene), intent(in) :: scn
    integer, intent(in) :: n_stokes

    double precision, intent(inout) :: radiance_calc_work_hi_stokes(:,:)

    double precision, allocatable, intent(inout) :: dI_dgas(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dsurf(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dTemp(:,:)
    double precision, allocatable, intent(inout) :: dI_dAOD(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dAHeight(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dpsurf(:,:)
    double precision, allocatable, intent(inout) :: dI_dtaugas(:,:,:)

    logical, intent(in) :: need_gas_averaging_kernels

    logical, intent(inout) :: xrtm_success

    ! Local variables

    ! Function name for logger
    character(len=*), parameter :: fname = "solve_RT_problem_XRTM"
    character(len=999) :: tmp_str

    double precision, allocatable :: ltau(:,:,:)
    double precision, allocatable :: lomega(:,:,:)
    double precision, allocatable :: lsurf(:,:)

    double precision, allocatable :: ltau_pca(:,:,:)
    double precision, allocatable :: lomega_pca(:,:,:)
    double precision, allocatable :: lsurf_pca(:,:)

    double precision, allocatable :: monochr_radiance(:,:)
    double precision, allocatable :: monochr_weighting_functions(:,:,:)

    double precision, allocatable :: binned_results_lo(:,:,:,:)
    double precision, allocatable :: binned_weighting_functions_lo(:,:,:,:,:)
    double precision, allocatable :: binned_results_hi(:,:,:,:)
    double precision, allocatable :: binned_weighting_functions_hi(:,:,:,:,:)

    type(PCA_handler_t) :: PCA_handler
    type(string), allocatable :: PCA_xrtm_solvers_lo(:), PCA_xrtm_solvers_hi(:)
    type(scene), allocatable :: PCA_scn(:,:) ! (bin, eof)

    integer :: bin, eof

    integer :: num_act_lev
    integer :: num_lay
    integer :: num_wl
    integer :: n_mom
    integer :: max_pfmom
    ! XRTM Radiative Transfer model handler
    type(xrtm_type) :: xrtm
    ! Actual options variable
    integer, allocatable :: xrtm_options(:)
    ! Actual solvers variable
    integer, allocatable :: xrtm_solvers(:)
    ! XRTM surface kernels
    integer :: xrtm_kernels(1)
    ! XRTM error variable
    integer :: xrtm_error
    ! How many derivatives does XRTM need to calculate?
    integer :: xrtm_n_derivs
    !
    integer :: xrtm_streams
    !
    integer, allocatable :: n_coef(:)
    !
    integer :: i, j, q
    !
    logical :: PCA_success
    double precision :: cpu_start_pca, cpu_stop_pca
    double precision :: pca_overhead_time

    ! Initialize success variable
    xrtm_success = .false.

    ! Grab the number of active levels for easy access
    num_act_lev = scn%num_active_levels
    num_lay = num_act_lev - 1
    num_wl = size(scn%op%wl)
    max_pfmom = max(scn%max_pfmom, 3)

    allocate(n_coef(num_lay))

    if (n_stokes == 1) then
       n_mom = 1
    else if (n_stokes == 3) then
       n_mom = 6
    else
       call logger%fatal(fname, "N_Stokes is neither 1 or 3 - not sure what to do!")
       stop 1
    end if

    ! Assign RT weighting function <-> SV Jacobian mapping
    call assign_RT_jacobian_indices(SV, need_gas_averaging_kernels, num_lay)


    ! -----------------------------------------------------------------------
    ! RT strategy formulation
    ! -----------------------------------------------------------------------
    !
    ! Unfortunately, there is a large number of possible RT solvers + fast RT
    ! methods that can be used in and with XRTM. For example, a monochromatic
    ! run can use different solvers whose results are all added up. For other
    ! fast RT methods, one might need to split the solvers up into different
    ! groups (e.g. low- and high-stream paths). And they all need to harmonize
    ! with and without weighting functions as well.
    !
    !
    ! -----------------------------------------------------------------------

    if (CS_win%RT_strategy == "") then
       call logger%fatal(fname, "Must have an RT strategy for XRTM!")
       stop 1
    end if

    if (CS_win%RT_strategy%lower() == "monochromatic") then

       call logger%debug(fname, "Using monochromatic RT strategy.")

       ! --------------------------------------------------------------------
       ! Monochromatic RT
       ! --------------------------------------------------------------------
       !
       ! This is the most straight-forward option. Using the user-supplied
       ! RT model and options, we run through the entire band, point for point.
       ! --------------------------------------------------------------------

       ! This function "translates" our verbose configuration file options
       ! into proper XRTM language and sets the corresponding options and solvers.

       if (.not. allocated(CS_win%RT_streams)) then
          xrtm_streams = 1
       else

          if (mod(CS_win%RT_streams(1), 2) /= 0) then
             call logger%error(fname, "Numnber of RT streams must be even!")
             stop 1
          end if

          xrtm_streams = CS_win%RT_streams(1) / 2
       end if

       call setup_XRTM( &
            CS_win%xrtm_options, &
            CS_win%xrtm_solvers, &
            CS_win%do_polarization, &
            CS_win%constant_coef, &
            max_pfmom, &
            xrtm_streams, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_kernels, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to setup_XRTM unsuccessful.")
          return
       end if

       xrtm_n_derivs = SV%num_rt_wf

       ! Allocate the arrays which hold the radiances and weighting functions
       allocate(monochr_radiance(num_wl, scn%num_stokes))
       allocate(monochr_weighting_functions(num_wl, scn%num_stokes, xrtm_n_derivs))

       monochr_radiance(:,:) = 0.0d0
       monochr_weighting_functions(:,:,:) = 0.0d0

       allocate(ltau(num_wl, xrtm_n_derivs, num_lay))
       allocate(lomega(num_wl, xrtm_n_derivs, num_lay))
       allocate(lsurf(num_wl, xrtm_n_derivs))

       ltau(:,:,:) = 0.0d0
       lomega(:,:,:) = 0.0d0
       lsurf(:,:) = 0.0d0

       call calculate_XRTM_lgas_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_ltaugas_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_ltemp_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_lpsurf_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_laerosol_aod_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_laerosol_height_inputs(SV, scn, CS_aerosol, ltau, lomega)
       call calculate_XRTM_lalbedo_inputs(SV, scn, lsurf)

       ! -------------------------------------------
       !
       ! This outer loop runs over all supplied
       ! XRTM solvers, and the results for jacobians
       ! AND radiances will be ADDED up. The user
       ! must know what they are doing.
       !
       ! -------------------------------------------

       do i = 1, size(xrtm_solvers)

          ! This is where the XRTM handler is created

          call create_XRTM( &
               xrtm, & ! XRTM handler
               xrtm_options(i), & ! XRTM options bitmask
               xrtm_solvers(i), & ! XRTM solvers combined bitmask
               max_pfmom, & ! Max coef - either 3 for Rayleigh, or whatever we have for aerosols
               xrtm_streams, & ! Quadrature points
               n_stokes, & ! Number of stokes coeffs
               xrtm_n_derivs, & ! Number of derivatives
               num_lay, & ! Number of layers
               1, & ! Number of surface kernels
               50, & ! Number of kernel quadrature points
               xrtm_kernels, & ! XRTM surface kernels
               xrtm_success) ! Call successful?

          if (.not. xrtm_success) then
             call logger%error(fname, "Call to create_XRTM unsuccessful.")
             call xrtm_destroy_f90(xrtm, xrtm_error)
             return
          end if

          ! This calculates all radiances and weighting functions. It DOES NOT map the weighting function
          ! array back to partial derivatives. This is done right at the end of this subroutine.
          ! NOTE
          ! The results of this particular XRTM pass are ADDED to the radiance and derivative arrays,
          ! so all contributions added up here.

          call XRTM_monochromatic( &
               xrtm, &
               xrtm_n_derivs, &
               scn, &
               SV, &
               CS_win, &
               CS_general, &
               CS_aerosol, &
               ltau, &
               lomega, &
               lsurf, &
               monochr_radiance, &
               monochr_weighting_functions, &
               xrtm_success &
               )

          if (.not. xrtm_success) then
             call logger%error(fname, "Call to XRTM_monochromatic unsuccessful.")
             call xrtm_destroy_f90(xrtm, xrtm_error)
             return
          end if

          ! After calculations are done, we really don't need the XRTM object anymore,
          ! so independent of the success, we can safely destroy it.
          !
          ! XRTM must be destroyed at some point, otherwise it will just keep
          ! creating arrays and filling up memory (quickly!). Why not destroy
          ! it right here - all the results are transferred to various arrays
          ! so we do not really need the XRTM instance anymore.
          ! If the calculation of radiances failed for whatever reason, the program
          ! returns here anyway such that XRTM is destroyed as well.
          call xrtm_destroy_f90(xrtm, xrtm_error)

       end do

       radiance_calc_work_hi_stokes(:,:) = monochr_radiance(:,:)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to calculate_XRTM_radiance unsuccessful.")
          return
       end if

       deallocate(ltau, lomega, lsurf)

       call logger%debug(fname, "Mapping weighting functions back to Jacobians.")

    else if (CS_win%RT_strategy%lower() == "pca") then

       ! --------------------------------------------------------------------
       ! PCA-based fast RT
       ! --------------------------------------------------------------------
       !
       ! Look at Somkuti et al. 2017 for details
       ! --------------------------------------------------------------------

       if (.not. allocated(CS_win%RT_streams)) then
          call logger%error(fname, "Sorry - must specify RT streams for PCA method.")
          stop 1
       end if

       if (size(CS_win%RT_streams) /= 2) then
          call logger%error(fname, "Sorry - must specify two different "&
               // "RT streams for PCA method")
          stop 1
       end if



       pca_overhead_time = 0.0d0

       cpu_start_pca = get_cpu_time()

       call logger%debug(fname, "Using PCA-based fast RT strategy.")

       ! Initialize the PCA object
       call initialize_PCA(scn, PCA_handler)
       ! Perform spectral binning
       call assign_pushkar_bins(scn, PCA_handler)
       ! Bail-out: if unassigned points remain, exit
       if (count(PCA_handler%assigned_bin(:) == -1) > 0) then
          call logger%error(fname, "Unassigned spectral points remain.")
          return
       end if

       ! Allocate some of the PCA object matrices
       call allocate_first_PCA(PCA_handler)

       ! Construct optical property matrices
       call create_F_matrix(scn, PCA_handler)

       ! Perform optical property matrix forward transform
       call forward_transform_F_matrix(PCA_handler)

       ! Perform optical property matrix mean removal
       ! (also computes mean optical states)
       call F_matrix_mean_removal(PCA_handler)

       ! Perform principal component analysis on mean-removed
       ! optical state matrices. Also allocates the necessary
       ! matrices within PCA_handler
       call perform_PCA(PCA_handler, PCA_success)

       if (.not. PCA_success) then
          call logger%error(fname, "Could not successfully perform PCA.")
          return
       end if

       ! Check here if any of the perturbed optical properties is a non-finite
       ! value..
       do j = 1, PCA_handler%N_bin
          if (.not. all(ieee_is_finite(PCA_handler%PCA_bin(j)%pert_opt_states(:,:)))) then
             call logger%error(fname, "Sorry - non-finite value encountered after PCA.")
             return
          end if
       end do
       ! Create scene object array from PCA data
       ! (mean and perturbed states, which we can then stick into XRTM)
       call create_scenes_from_PCA(PCA_handler, CS_aerosol, scn, PCA_scn)

       ! --------------------------------------------------------------------
       ! PCA set-up done!
       !
       ! All required values are now available inside of the object
       ! "PCA_handler"
       ! --------------------------------------------------------------------

       cpu_stop_pca = get_cpu_time()

       pca_overhead_time = pca_overhead_time + cpu_stop_pca - cpu_start_pca

       ! For the time being, we want to keep the PCA method as:
       ! Low accuracy: 2S
       ! High accuracy: EIG_BVP
       ! (maybe change later on..)
       ! NOTE
       ! In XRTM, TWO_STREAM computes both single scattering AND diffuse 2S contributions.
       ! This is different in LIDORT, 2STREAM etc.

       allocate(PCA_xrtm_solvers_lo(1))
       allocate(PCA_xrtm_solvers_hi(1))

       PCA_xrtm_solvers_lo(1) = "TWO_STREAM"
       PCA_xrtm_solvers_hi(1) = "EIG_BVP"

       call setup_XRTM( &
            CS_win%xrtm_options, &
            PCA_xrtm_solvers_lo, &
            .false., & ! Polarization? Not for low accuracy run!
            CS_win%constant_coef, &
            max_pfmom, &
            CS_win%RT_streams(1) / 2, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_kernels, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to setup_XRTM unsuccessful.")
          return
       end if

       ! Set up number of derivatives needed by XRTM
       xrtm_n_derivs = SV%num_rt_wf

       ! --------------------------------------------------------------------
       ! Perform low-accuracy line-by-line run for TWO_STREAM
       ! contributions.
       !
       ! All required values are now available inside of the object
       ! "PCA_handler"
       ! --------------------------------------------------------------------

       call create_XRTM( &
            xrtm, & ! XRTM handler
            xrtm_options(1), & ! XRTM options bitmask
            xrtm_solvers(1), & ! XRTM solvers combined bitmask
            max_pfmom, & ! Max coef - either 3 for Rayleigh or max_pfmom for aerosols
            CS_win%RT_streams(1) / 2, & ! Quadrature points per hemisphere
            n_stokes, & ! Number of stokes coeffs
            xrtm_n_derivs, & ! Number of derivatives
            num_lay, & ! Number of layers
            1, & ! Number of surface kernels
            50, & ! Number of kernel quadrature points
            xrtm_kernels, & ! XRTM surface kernels
            xrtm_success) ! Call successful?

       ! Allocate the arrays which hold the radiances and weighting functions
       ! for the line-by-line low accuracy run.
       allocate(monochr_radiance(size(scn%op%wl), scn%num_stokes))
       allocate(monochr_weighting_functions(size(scn%op%wl), scn%num_stokes, xrtm_n_derivs))

       monochr_radiance(:,:) = 0.0d0
       monochr_weighting_functions(:,:,:) = 0.0d0

       allocate(ltau(num_wl, xrtm_n_derivs, num_lay))
       allocate(lomega(num_wl, xrtm_n_derivs, num_lay))
       allocate(lsurf(num_wl, xrtm_n_derivs))

       allocate(ltau_pca(1, xrtm_n_derivs, num_lay))
       allocate(lomega_pca(1, xrtm_n_derivs, num_lay))
       allocate(lsurf_pca(1, xrtm_n_derivs))

       ltau(:,:,:) = 0.0d0
       lomega(:,:,:) = 0.0d0
       lsurf(:,:) = 0.0d0

       ltau_pca(:,:,:) = 0.0d0
       lomega_pca(:,:,:) = 0.0d0
       lsurf_pca(:,:) = 0.0d0

       call calculate_XRTM_lgas_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_ltaugas_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_ltemp_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_lpsurf_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_laerosol_aod_inputs(SV, scn, ltau, lomega)
       call calculate_XRTM_laerosol_height_inputs(SV, scn, CS_aerosol, ltau, lomega)
       call calculate_XRTM_lalbedo_inputs(SV, scn, lsurf)

       call XRTM_monochromatic( &
            xrtm, &
            xrtm_n_derivs, &
            scn, &
            SV, &
            CS_win, &
            CS_general, &
            CS_aerosol, &
            ltau, &
            lomega, &
            lsurf, &
            monochr_radiance, &
            monochr_weighting_functions, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to XRTM_monochromatic unsuccessful.")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! --------------------------------------------------------------------
       !
       !
       ! Perform binned calculations: low-accuracy
       ! (here we use the same XRTM options / solvers as before)
       !
       ! --------------------------------------------------------------------

       allocate(binned_results_lo(PCA_handler%N_bin, &
            -maxval(PCA_handler%PCA_bin(:)%N_EOF):maxval(PCA_handler%PCA_bin(:)%N_EOF), &
            1, scn%num_stokes))

       allocate(binned_weighting_functions_lo(PCA_handler%N_bin, &
            -maxval(PCA_handler%PCA_bin(:)%N_EOF):maxval(PCA_handler%PCA_bin(:)%N_EOF), &
            1, scn%num_stokes, xrtm_n_derivs))

       do bin = 1, PCA_handler%N_bin

          call extract_linputs_for_PCA(PCA_handler, bin, &
               ltau, lomega, lsurf, &
               ltau_pca, lomega_pca, lsurf_pca)

          do eof = -PCA_handler%PCA_bin(bin)%N_EOF, PCA_handler%PCA_bin(bin)%N_EOF

             binned_results_lo(bin, eof, :, :) = 0.0d0
             binned_weighting_functions_lo(bin, eof, :, :, :) = 0.0d0

             call XRTM_monochromatic( &
                  xrtm, &
                  xrtm_n_derivs, &
                  PCA_scn(bin, eof), &
                  SV, &
                  CS_win, &
                  CS_general, &
                  CS_aerosol, &
                  ltau_pca, &
                  lomega_pca, &
                  lsurf_pca, &
                  binned_results_lo(bin, eof, :, :), &
                  binned_weighting_functions_lo(bin, eof, :, :, :), &
                  xrtm_success)

             if (.not. xrtm_success) then
                call logger%error(fname, "Error while calculating low-accuracy bins.")
                call xrtm_destroy_f90(xrtm, xrtm_error)
                return
             end if

          end do
       end do

       ! --------------------------------------------------------------------
       !
       !
       ! Perform binned calculations: high-accuracy
       !
       !
       ! --------------------------------------------------------------------


       allocate(binned_results_hi(PCA_handler%N_bin, &
            -maxval(PCA_handler%PCA_bin(:)%N_EOF):maxval(PCA_handler%PCA_bin(:)%N_EOF), &
            1, scn%num_stokes))

       allocate(binned_weighting_functions_hi(PCA_handler%N_bin, &
            -maxval(PCA_handler%PCA_bin(:)%N_EOF):maxval(PCA_handler%PCA_bin(:)%N_EOF), &
            1, scn%num_stokes, xrtm_n_derivs))

       ! We need to create new XRTM solvers, options
       call xrtm_destroy_f90(xrtm, xrtm_error)

       if (xrtm_error /= 0) then
          call logger%error(fname, "Call to destroy XRTM unsuccessful.")
          return
       end if


       ! deallocate options and solvers
       deallocate(xrtm_options, xrtm_solvers)

       call setup_XRTM( &
            CS_win%xrtm_options, &
            PCA_xrtm_solvers_hi, &
            CS_win%do_polarization, &
            CS_win%constant_coef, &
            max_pfmom, &
            CS_win%RT_streams(2) / 2, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_kernels, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to setup_XRTM unsuccessful.")
          return
       end if

       call create_XRTM( &
            xrtm, & ! XRTM handler
            xrtm_options(1), & ! XRTM options bitmask
            xrtm_solvers(1), & ! XRTM solvers combined bitmask
            max_pfmom, & ! Max coef - either 3 for Rayleigh or max_pfmom for aerosols
            CS_win%RT_streams(2) / 2, & ! Quadrature points per hemisphere
            n_stokes, & ! Number of stokes coeffs
            xrtm_n_derivs, & ! Number of derivatives
            num_lay, & ! Number of layers
            1, & ! Number of surface kernels
            50, & ! Number of kernel quadrature points
            xrtm_kernels, & ! XRTM surface kernels
            xrtm_success) ! Call successful?

       do bin = 1, PCA_handler%N_bin

          call extract_linputs_for_PCA(PCA_handler, bin, &
               ltau, lomega, lsurf, &
               ltau_pca, lomega_pca, lsurf_pca)

          do eof = -PCA_handler%PCA_bin(bin)%N_EOF, PCA_handler%PCA_bin(bin)%N_EOF

             binned_results_hi(bin, eof, :, :) = 0.0d0
             binned_weighting_functions_hi(bin, eof, :, :, :) = 0.0d0

             call XRTM_monochromatic( &
                  xrtm, &
                  xrtm_n_derivs, &
                  PCA_scn(bin, eof), &
                  SV, &
                  CS_win, &
                  CS_general, &
                  CS_aerosol, &
                  ltau_pca, &
                  lomega_pca, &
                  lsurf_pca, &
                  binned_results_hi(bin, eof, :, :), &
                  binned_weighting_functions_hi(bin, eof, :, :, :), &
                  xrtm_success)

             if (.not. xrtm_success) then
                call xrtm_destroy_f90(xrtm, xrtm_error)
                call logger%error(fname, "Error while calculating high-accuracy bins.")
                return
             end if

          end do
       end do


       call xrtm_destroy_f90(xrtm, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Call to destroy XRTM unsuccessful.")
          return
       end if

       deallocate(ltau, lomega, lsurf)
       deallocate(ltau_pca, lomega_pca, lsurf_pca)

       ! --------------------------------------------------------------------
       !
       !
       ! Map the results from binned calculations back to wavelength space
       ! and add the results to the monochromatic calculations
       !
       !
       ! --------------------------------------------------------------------

       cpu_start_pca = get_cpu_time()

       call map_PCA_radiances(PCA_handler, &
            binned_results_lo(:,:,1,:), &
            binned_results_hi(:,:,1,:), &
            monochr_radiance)

       call map_PCA_jacobians(PCA_handler, &
            binned_weighting_functions_lo(:,:,1,:,:), &
            binned_weighting_functions_hi(:,:,1,:,:), &
            monochr_weighting_functions)

       cpu_stop_pca = get_cpu_time()

       pca_overhead_time = pca_overhead_time + cpu_stop_pca - cpu_start_pca

       write(tmp_str, '(A, F10.5)') "Total PCA overhead time: ",  pca_overhead_time
       call logger%debug(fname, trim(tmp_str))

       radiance_calc_work_hi_stokes(:,:) = monochr_radiance(:,:)

       do j=1, SV%num_aerosol_aod
          ! If this AOD retrieval is in log-space, we need to multiply by
          ! the (linear-space) AOD itself. df/d(ln(x)) = df/dx * x
          if (CS_win%aerosol_retrieve_aod_log(SV%aerosol_aod_idx_lookup(j))) then
             dI_dAOD(:,j,:) = dI_dAOD(:,j,:) * exp(SV%svsv(SV%idx_aerosol_aod(j)))
          end if
       end do

       ! Store aerosol layer height Jacobians
       do j=1, SV%num_aerosol_height
          if (CS_win%aerosol_retrieve_height_log(SV%aerosol_height_idx_lookup(j))) then
             dI_dAHeight(:,j,:) = dI_dAHeight(:,j,:) * exp(SV%svsv(SV%idx_aerosol_height(j)))
          end if
       end do



       deallocate(PCA_handler%PCA_bin)


    else

       call logger%fatal(fname, "RT strategy: '" // CS_win%RT_strategy%chars() &
            // "' not recognized. Valid options are..")

       call logger%fatal(fname, "(1) monochromatic")
       !call logger%fatal(fname, "(2) PCA (work in progress)")
       stop 1

    end if

    ! --------------------------------------------------------------------
    ! Atmospheric and surface RT weighting functions are now
    ! converted into proper SV-related Jacobians
    ! --------------------------------------------------------------------

    do i = 1, SV%num_gas
       dI_dgas(:,i,:) = monochr_weighting_functions(:,:,SV%idx_wf_gas(i))
    end do

    if (allocated(SV%idx_wf_taugas)) then
       do i = 1, num_lay
          dI_dtaugas(:,i,:) = monochr_weighting_functions(:,:,SV%idx_wf_taugas(i))
       end do
    end if

    do i = 1, SV%num_albedo
       dI_dsurf(:,i,:) = monochr_weighting_functions(:,:,SV%idx_wf_albedo(i))
    end do

    if (SV%num_temp > 0) then
       dI_dtemp(:,:) = monochr_weighting_functions(:,:,SV%idx_wf_temp(1))
    end if

    if (SV%num_psurf > 0) then
       dI_dpsurf(:,:) = monochr_weighting_functions(:,:,SV%idx_wf_psurf(1))
    end if

    do i = 1, SV%num_aerosol_aod
       dI_dAOD(:,i,:) = monochr_weighting_functions(:,:,SV%idx_wf_aerosol_aod(i))
    end do

    do i = 1, SV%num_aerosol_height
       dI_dAHeight(:,i,:) = monochr_weighting_functions(:,:,SV%idx_wf_aerosol_height(i))
    end do

    ! --------------------------------------------------------------------
    ! Aerosol jacobians might need to be converted from logspace
    ! --------------------------------------------------------------------

    do j=1, SV%num_aerosol_aod
       ! If this AOD retrieval is in log-space, we need to multiply by
       ! the (linear-space) AOD itself. df/d(ln(x)) = df/dx * x
       if (CS_win%aerosol_retrieve_aod_log(SV%aerosol_aod_idx_lookup(j))) then
          dI_dAOD(:,j,:) = dI_dAOD(:,j,:) * exp(SV%svsv(SV%idx_aerosol_aod(j)))
       end if
    end do

    ! Store aerosol layer height Jacobians
    do j=1, SV%num_aerosol_height
       if (CS_win%aerosol_retrieve_height_log(SV%aerosol_height_idx_lookup(j))) then
          dI_dAHeight(:,j,:) = dI_dAHeight(:,j,:) * exp(SV%svsv(SV%idx_aerosol_height(j)))
       end if
    end do



  end subroutine solve_RT_problem_XRTM




  !> @brief Calculates radiances and weighting functions for all wavelengths in the band
  subroutine XRTM_monochromatic( &
       xrtm, &
       n_derivs, &
       scn, &
       SV, &
       CS_win, &
       CS_general, &
       CS_aerosol, &
       ltau, &
       lomega, &
       lsurf, &
       radiance, &
       weighting_functions, &
       success )

    type(xrtm_type), intent(inout) :: xrtm
    integer, intent(in) :: n_derivs
    type(scene), intent(in) :: scn
    type(statevector), intent(in) :: SV
    type(CS_window_t), intent(in) :: CS_win
    type(CS_general_t), intent(in) :: CS_general
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    ! Linearized optical depth input (spectral, n_derivs, n_layers)
    double precision, allocatable, intent(in) :: ltau(:,:,:)
    ! Linearized single scatter alebdos (spectral, n_derivs, n_layers)
    double precision, allocatable, intent(in) :: lomega(:,:,:)
    ! Linearized surface parameters (spectral, n_derivs)
    double precision, allocatable, intent(in) :: lsurf(:,:)
    ! Container for radiance result (specetral, n_stokes)
    double precision, intent(inout) :: radiance(:,:)
    ! Container for weighting function results (spectral, n_stokes, n_derivs)
    double precision, intent(inout) :: weighting_functions(:,:,:)
    logical, intent(inout) :: success
    !double precision, optional :: ltau_input(:,:,:) ! (spectral, n_derivs, n_layers)
    !double precision, optional :: lomega_input(:,:,:) ! (spectral, n_derivs, n_layers)

    character(len=*), parameter :: fname = "XRTM_monochromatic"
    character(len=999) :: tmp_str


    ! Were the calls to XRTM successful?
    logical :: xrtm_success
    ! Were the coef calculations successful?
    logical :: coef_success
    ! Number of active levels
    integer :: n_levels
    ! Number of (active) layers
    integer :: n_layers
    ! Number of wavelengths
    integer :: n_wl
    ! Aerosol index to map SV element to aerosol number
    integer :: aer_idx
    ! Number of phase matrix elements
    integer :: n_mom
    ! Number of used phasefunction expansion coefficients, per layer
    integer, allocatable :: n_coef(:)
    ! Linearized surface inputs (spectral, n_derivs)
    !double precision, allocatable :: lsurf(:,:)
    ! Aerosol shape factor required for Jacobians (n_layers)
    double precision, allocatable :: aero_fac(:)
    ! Phase function coefficients for left hand side (lower-wavelength edge of band)
    ! (n_pfmom, n_mom, n_layer)
    double precision, allocatable :: coef_left(:,:,:)
    ! Phase function coefficients for right hand side (lower-wavelength edge of band)
    ! (n_pfmom, n_mom, n_layer)
    double precision, allocatable :: coef_right(:,:,:)
    ! Linearized phase function coefficients for left hand side (lower-wavelength edge of band)
    ! (n_pfmom, n_mom, n_derivs, n_layer)
    double precision, allocatable :: lcoef_left(:,:,:,:)
    ! Linearized phase function coefficients for right hand side (lower-wavelength edge of band)
    ! (n_pfmom, n_mom, n_derivs, n_layer)
    double precision, allocatable :: lcoef_right(:,:,:,:)

    ! Counter for weighting functions (debugging mostly)
    integer :: deriv_counter
    ! Loop variables
    integer :: i, j, l

    success = .false.
    coef_success = .false.
    xrtm_success = .false.

    ! Derive some useful values here
    n_levels = scn%num_active_levels
    n_layers = n_levels - 1
    n_wl = size(scn%op%wl)

    ! Allocate the RT inputs, and zero them
    !allocate(ltau(n_wl, n_derivs, n_layers))
    !allocate(lomega(n_wl, n_derivs, n_layers))
    !allocate(lsurf(n_wl, n_derivs))

    !ltau(:,:,:) = 0.0d0
    !lomega(:,:,:) = 0.0d0
    !lsurf(:,:) = 0.0d0

    if (scn%num_stokes == 1) then
       n_mom = 1
    else if (scn%num_stokes == 3) then
       n_mom = 6
    else
       call logger%fatal(fname, "n_stokes is neither 1 or 3 - not sure what to do!")
       stop 1
    end if

    ! -------------------------------------------------------------------------------
    ! Calculate all required linearized inputs, and do some bookkeeping and debugging
    ! to keep track of the weighting function positions.
    ! -------------------------------------------------------------------------------

    call logger%debug(fname, "Calculating linearized inputs for RT.")

    call compute_coef_at_wl(scn, CS_aerosol, SV, scn%op%wl(1), n_mom, n_derivs, &
         scn%op%ray_tau(1,:), scn%op%aer_sca_tau(1,:,:), &
         coef_left, lcoef_left, coef_success)

    if (.not. coef_success) then
       return
    end if

    ! Only need right-side coeffs if we have more than one spectral point in this
    ! scene collection.
    if (n_wl > 1) then

       call compute_coef_at_wl(scn, CS_aerosol, SV, scn%op%wl(n_wl), n_mom, n_derivs, &
            scn%op%ray_tau(n_wl,:), scn%op%aer_sca_tau(n_wl,:,:), &
            coef_right, lcoef_right, coef_success)

       if (.not. coef_success) then
          return
       end if

    else
       allocate(coef_right, source=coef_left)
       allocate(lcoef_right, source=lcoef_left)

       coef_right = coef_left
       lcoef_right = lcoef_left
    end if

    call logger%debug(fname, "Determining number of phase function coefficients per layer used in XRTM.")
    ! Figure out how many coefficients we are giving to XRTM
    ! The default is 3 for Rayleigh-only

    allocate(n_coef(n_layers))
    n_coef(:) = max(3, scn%max_pfmom)

    !n_coef(:) = 3
    !if (scn%num_aerosols > 0) then
    ! Loop through all layers and look at how much aerosol
    ! extinction is present if below some threshold value,
    ! we just assume that 3 coefficients is enough for a
    ! (at best) Rayleigh-only layer
    !   do l = 1, n_layers
    !      if (sum(scn%op%aer_ext_tau(1,l,:)) < 1.0d-6) then
    !         n_coef(l) = 3
    !      else
    !         n_coef(l) = scn%max_pfmom
    !      end if
    !   end do
    !end if

    ! ------------------------------------------
    !
    ! Run the monochromatic radiance calculation
    !
    ! ------------------------------------------

    ! We pass all per-wavelength quantities to the function which then
    ! calculates all per-wavelength radiances as well as per-wavelength
    ! Jacobians.

    call calculate_XRTM_radiance( &
         xrtm, & ! XRTM handler
         SV, & ! State vector structure - needed for decoding SV positions
         CS_general, & ! Control structure general section
         CS_aerosol, & ! Control structure aerosol array
         scn%op%wl, & ! per pixel wavelength
         scn%SZA, & ! solar zenith angle
         scn%VZA, & ! viewing zenith angle
         scn%SAA, & ! solar azimuth angle
         scn%VAA, & ! viewing azimuth angle
         scn%atm%altitude_levels, & ! per-level altitude
         scn%op%albedo, & ! Surface albedo
         scn%op%layer_tau, & ! per-wl per-layer optical depths
         scn%op%layer_omega, & ! per-wl per-layer SSA
         ltau, & ! per-wl per layer optical depth derivatives
         lomega, & ! per-wl per layer SSA derivatives
         lsurf, & ! per-wl surface derivatives
         scn%num_stokes, & ! Number of Stokes parameters to calculate
         n_derivs, & ! Number of derivatives to be calculated
         n_layers, & ! Number of atmospheric layers
         n_coef, & ! Number of per-layer coefficients to be used by XRTM
         CS_win%constant_coef, & ! Keep coefficients constant along window?
         radiance, & ! Output radiances
         weighting_functions, & ! Output weighting functions
         xrtm_success, & ! Success indicator
         coef_left=coef_left, & ! per-wl phasefunction coefficients
         lcoef_left=lcoef_left, & ! per-wl phasefunction derivatives
         coef_right=coef_right, & ! per-wl phasefunction coefficients
         lcoef_right=lcoef_right) ! per-wl phasefunction derivatives

    if (.not. xrtm_success) then
       call logger%error(fname, "Call to calculate_XRTM_radiance unsuccessful.")
       return
    end if

    success = .true.

  end subroutine XRTM_monochromatic



  !> @details Phasefunction expansion coefficients can be provided in two ways: either supply directly the
  !> (coef, lcoef) arrays which have a wavelength dimension. Alternatively, supply the (coef_left, lcoef_left)
  !> and (coef_right, lcoef_right) arrays without a wavelength dimension, and linear interpolation will be used
  !> to obtain the values at each wavelength.
  !> NOTE that radiance and derivative values are added to "radiance" and "derivs" arrays, so you MUST
  !> make sure that they are properly initialized/zeroed before they are passed into this function.
  subroutine calculate_XRTM_radiance( &
       xrtm, &
       SV, &
       CS_general, &
       CS_aerosol, &
       wavelengths, SZA, VZA, SAA, VAA, &
       altitude_levels, albedo, &
       layer_tau, &
       layer_omega, &
       ltau, &
       lomega, &
       lsurf, &
       n_stokes, n_derivs, n_layers, n_coef, &
       keep_coef_constant, &
       radiance, &
       derivs, &
       success, &
       coef, &
       lcoef, &
       coef_left, &
       lcoef_left, &
       coef_right, &
       lcoef_right)

    implicit none
    type(xrtm_type), intent(inout) :: xrtm
    type(statevector), intent(in) :: SV
    type(CS_general_t), intent(in) :: CS_general
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)

    double precision, intent(in) :: wavelengths(:)
    double precision, intent(in) :: SZA
    double precision, intent(in) :: VZA
    double precision, intent(in) :: SAA
    double precision, intent(in) :: VAA
    double precision, intent(in) :: albedo(:)
    double precision, intent(in) :: altitude_levels(:)
    double precision, intent(in) :: layer_tau(:,:)
    double precision, intent(in) :: layer_omega(:,:)
    double precision, intent(in) :: ltau(:,:,:)
    double precision, intent(in) :: lomega(:,:,:)
    double precision, intent(in) :: lsurf(:,:)

    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers
    integer, intent(in) :: n_coef(:)
    logical, intent(in) :: keep_coef_constant

    double precision, intent(inout) :: radiance(:,:)
    double precision, intent(inout) :: derivs(:,:,:)

    double precision, allocatable, intent(in), optional :: coef(:,:,:,:) !(pfmom, elem, n_layers, spectral)
    double precision, allocatable, intent(in), optional :: lcoef(:,:,:,:,:) !(pfmom, elem, n_derivs, n_layers, spectral)

    double precision, allocatable, intent(in), optional :: coef_left(:,:,:) !(pfmom, elem, n_layers)
    double precision, allocatable, intent(in), optional :: lcoef_left(:,:,:,:) !(pfmom, elem, n_derivs, n_layers)
    double precision, allocatable, intent(in), optional :: coef_right(:,:,:) !(pfmom, elem, n_layers)
    double precision, allocatable, intent(in), optional :: lcoef_right(:,:,:,:) !(pfmom, elem, n_derivs, n_layers)

    logical, intent(inout) :: success

    ! Local
    integer :: xrtm_options
    integer :: xrtm_solver

    double precision, allocatable :: I_p(:,:,:,:)
    double precision, allocatable :: I_m(:,:,:,:)
    double precision, allocatable :: K_p(:,:,:,:,:)
    double precision, allocatable :: K_m(:,:,:,:,:)

    character(len=*), parameter :: fname = "calculate_XRTM_radiance"
    character(len=999) :: tmp_str
    character(len=999) :: fmt_str
    integer :: xrtm_error

    double precision :: out_thetas(1)
    double precision :: out_phis(1,1)
    integer :: out_levels(1)
    !integer :: n_coef(n_layers)
    double precision :: single_ltau(n_derivs, n_layers)
    double precision :: single_lomega(n_derivs, n_layers)
    double precision :: single_lsurf(n_derivs)

    double precision :: tau(n_layers), omega(n_layers)
    ! (pfmom, elem, n_layers)
    double precision, allocatable :: this_coef(:,:,:), this_lcoef(:,:,:,:)

    double precision :: wl_fac
    integer :: N_spec
    integer :: i, q, m, k

    double precision :: cpu_start, cpu_end
    double precision :: cpu_start2, cpu_end2
    double precision :: pure_XRTM_duration
    success = .false.

    ! Make some sanity checks for correct usage of this function.

    if (present(coef) .and. .not. present(lcoef)) then
       call logger%fatal(fname, "Argument coef and lcoef must both be present!")
       stop 1
    end if

    ! If coef and lcoef were not supplied, we MUST have the other four arrays
    if (.not. (present(coef) .and. present(lcoef))) then
       if (.not. (present(coef_left) .and. present(lcoef_left) &
            .and. present(coef_right) .and. present(lcoef_right))) then
          call logger%fatal(fname, "Arguments coef_left, lcoef_left, coef_right, lcoef_right " &
               // "must ALL be present!")
          stop 1
       end if
    end if

    ! Also, we can't mix the two
    if ((present(coef) .or. present(lcoef)) .and. &
         (present(coef_left) .or. present(lcoef_left) .or. &
         present(coef_right) .or. present(lcoef_right))) then

       call logger%fatal(fname, "Cannot mix the two coef supply methods!")
       stop 1
    end if

    xrtm_options = xrtm_get_options_f90(xrtm)
    xrtm_solver = xrtm_get_solvers_f90(xrtm)

    N_spec = size(wavelengths)

    ! Viewing geometry
    out_thetas(1) = VZA
    out_phis(1,1) = VAA
    out_levels(1) = 0

    single_ltau(:,:) = 0.0d0
    single_lomega(:,:) = 0.0d0
    single_lsurf(:) = 0.0d0

    ! Depending on which coefficient strategy we use, we allocate the coefficients
    ! according to the input files.
    ! There are two cases, and two cases within those two:
    !
    ! a) The user supplied per-wavelength coeffs
    !   In this case, we allocate the "this_coef" and "this_lcoef" fields, according to
    !   the shape of the supplied data.
    !
    !   a1) If the user wants to keep the scattering coefficients constant along
    !       the band, then we just grab the values from the band center
    !   a2) If not, there is nothing to do, but to just take the per-wavelength
    !       values and stick them into XRTM
    !
    ! b) The user supplied coefs at the two band edges "coef_left" and "coef_right"
    !   In this case, the "this_coef" and "this_lcoef" fields are allocated
    !   according to the shape of "coef_left".
    !
    !   b1) If the user wants to keep the scattering coefficients constant along
    !       the band, we grab the band-center values using the "left" and "right"
    !       coefficient arrays.
    !   b2) If not, we make the same computation in the spectral loop, using a
    !       wavelength interpolation factor derived from the current spectral
    !       point (wl_fac)

    if (present(coef)) then
       allocate(this_coef(size(coef, 1), size(coef, 2), size(coef, 3)))
       allocate(this_lcoef(size(lcoef, 1), size(lcoef, 2), size(lcoef, 3), size(lcoef, 4)))

       ! Spectrally flat scattering properties, this operation then only
       ! needs to be done once for the entire RT procedure

       if (keep_coef_constant) then
          ! Grab the value from the center
          wl_fac = 0.5d0
          this_coef(:,:,:) = coef(:,:,:,int(size(coef, 4) / 2))
          this_lcoef(:,:,:,:) = lcoef(:,:,:,:,int(size(lcoef, 5) / 2))
       end if

    else
       allocate(this_coef(size(coef_left, 1), size(coef_left, 2), size(coef_left, 3)))
       allocate(this_lcoef(size(lcoef_left, 1), size(lcoef_left, 2), &
            size(lcoef_left, 3), size(lcoef_left, 4)))

       if (keep_coef_constant) then
          ! Choose band center, grab a mixture of left and right l/coefs
          wl_fac = 0.5d0
          this_coef(:,:,:) = (1.0d0 - wl_fac) * coef_left(:,:,:) + wl_fac * coef_right(:,:,:)
          this_lcoef(:,:,:,:) = (1.0d0 - wl_fac) * lcoef_left(:,:,:,:) + wl_fac * lcoef_right(:,:,:,:)
       end if

    end if




    ! Allocate output containers - these will not change with wavelength, hence
    ! we only need to do this once per microwindow
    allocate(I_p(n_stokes, 1, 1, 1))
    allocate(I_m(n_stokes, 1, 1, 1))
    allocate(K_p(n_stokes, 1, 1, n_derivs, 1))
    allocate(K_m(n_stokes, 1, 1, n_derivs, 1))

    ! Some general set-up's that are not expected to change, and can
    ! be safely hard-coded.
    call xrtm_set_fourier_tol_f90(xrtm, .0001d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_fourier_tol_f90")
       return
    endif

    ! Set isotropic emission at TOA and BOA to zero.
    call xrtm_set_F_iso_top_f90(xrtm, 0.0d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_F_iso_top_f90")
       return
    endif

    call xrtm_set_F_iso_bot_f90(xrtm, 0.0d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_F_iso_bot_f90")
       return
    endif

    ! ----------------------------------------------------------------------------
    ! First, plug in the values that XRTM needs for its calculations, which DO NOT
    ! depend on wavelength (i.e. do not need to be called repeatedly), such as
    ! viewing geometry.
    ! ----------------------------------------------------------------------------

    if (iand(xrtm_options, XRTM_OPTION_PSA) /= 0) then
       ! The following two calls are only valid if XRTM has been created
       ! with the flag to use pseudo-spherical approximation for a curved
       ! atmosphere.

       ! Set Earth's radius
       call xrtm_set_planet_r_f90(xrtm, EARTH_EQUATORIAL_RADIUS, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_planet_r_f90")
          return
       end if

       ! Set model atmosphere altitude levels
       call xrtm_set_levels_z_f90(xrtm, altitude_levels(1:n_layers+1), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_levels_z_f90")
          return
       end if

    end if

    if (iand(xrtm_solver, XRTM_SOLVER_SOS) /= 0) then

       call xrtm_set_sos_params_f90(xrtm, 10, 1.0d0, 0.01d0, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_sos_params_f90")
          return
       endif
    end if

    if (iand(xrtm_solver, XRTM_SOLVER_PADE_ADD) /= 0) then
       call xrtm_set_pade_params_f90(xrtm, 0, 2, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_pade_params_f90")
          return
       endif
    end if

    ! Set the output level to TOA
    call xrtm_set_out_levels_f90(xrtm, out_levels, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_out_levels_f90")
       return
    endif

    ! We want to calculate sun-normalized values, so set F_0 to 1.0
    call xrtm_set_F_0_f90(xrtm, 1.0d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_F_0_f90")
       return
    end if

    ! Plug in the solar zenith angle
    call xrtm_set_theta_0_f90(xrtm, SZA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_theta_0_f90")
       return
    end if

    ! Plug in the viewing zenith angle
    call xrtm_set_out_thetas_f90(xrtm, out_thetas, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_out_thetas_f90")
       return
    end if

    ! Plug in the solar azimuth angle
    call xrtm_set_phi_0_f90(xrtm, SAA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_phi_0_f90")
       return
    end if


    ! -----------------------------------------------------------------------------------
    ! Then, loop through the wavelengths and stick in wavelength-dependent quantities,
    ! i.e. coefficients, optical depths and single-scattering albedos etc., and call
    ! the XRTM solver to get radiances - this is obviously a performance-critical section
    ! and is easily the most costly portion of the entire forward model.
    !
    ! NOTE
    ! This subroutine is intrinsically monochromatic, so the only thing you really need
    ! to make sure is that the array positions between optical paramaters make sense.
    !
    ! NOTE
    ! Most XRTM functions to set optical inputs return an error variable. Ideally we want
    ! to check every one of those variables whether the call was successful or not. It
    ! turns out, however, that these checks within the monochromatic loop have a fairly
    ! significant performance penalty, and can blow up the calculation times by a factor
    ! of 2x or so. I have thus decided to make these checks only in DEBUG mode, and use
    ! preprocessor directives to include/exclude them. This is generally bad practice,
    ! since an error here could only be discovered only when compiling in DEBUG mode.
    !
    ! -----------------------------------------------------------------------------------

#ifdef DEBUG
    call logger%debug(fname, "GASBAG was compiled in DEBUG mode. " &
         // "Additional checks for XRTM will be performed.")
#endif


    pure_XRTM_duration = 0.0d0
    cpu_start = get_cpu_time()

    do i=1, N_spec

       ! TOTAL atmospheric optical properties - these go into the
       ! RT code for the forward calculation.
       tau(:) = layer_tau(i, 1:n_layers)
       omega(:) = layer_omega(i, 1:n_layers)

       single_ltau(:,:) = ltau(i,:,1:n_layers)
       single_lomega(:,:) = lomega(i,:,1:n_layers)
       single_lsurf(:) = lsurf(i,:)

#ifdef DEBUG
       if (any(ieee_is_nan(tau))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for TAU at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(omega))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for OMEGA at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(single_ltau))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for L_TAU at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(single_lomega))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for L_OMEGA at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(single_lsurf))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for L_SURF at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if
#endif


       ! Interpolate phasefunction coefficients at wavelength, but only if
       ! requested by the user.
       ! NOTE WARNING
       ! This is a fairly slow section of the code and can increase processing
       ! time by a couple of factors easily!


       if (.not. keep_coef_constant) then
          if (present(coef_left)) then

             wl_fac = (wavelengths(i) - wavelengths(1)) / (wavelengths(N_spec) - wavelengths(1))
             this_coef(:,:,:) = (1.0d0 - wl_fac) * coef_left(:,:,:) + wl_fac * coef_right(:,:,:)
             this_lcoef(:,:,:,:) = (1.0d0 - wl_fac) * lcoef_left(:,:,:,:) + wl_fac * lcoef_right(:,:,:,:)
          else
             ! Plug in the layer and wavelength-dependent scattering matrix
             ! Spectrally varying scattering properties
             this_coef(:,:,:) = coef(:,:,:,i)
             this_lcoef(:,:,:,:) = lcoef(:,:,:,:,i)
          end if
       end if

#ifdef DEBUG
       if (any(ieee_is_nan(this_coef))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for COEF at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(this_lcoef))) then
          write(tmp_str, '(A, G0.1)') "NaN(s) found for L_COEF at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if
#endif


       ! Plug in optical depth
       call xrtm_set_ltau_n_f90(xrtm, tau, xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_n_f90")
          return
       end if
#endif

       ! Plug in the layer optical depth derivatives
       call xrtm_set_ltau_l_nn_f90(xrtm, single_ltau, xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_nn_f90")
          return
       end if
#endif

       ! Plug in the single-scatter albedo
       call xrtm_set_omega_n_f90(xrtm, omega, xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_n_f90")
          return
       end if
#endif

       ! Plug in the single-scatter albedo derivatives
       call xrtm_set_omega_l_nn_f90(xrtm, single_lomega, xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_l_nn_f90")
          return
       end if
#endif

       if ((.not. keep_coef_constant) .or. (i == 1)) then
          ! Plug in the phase function expansion moments

          ! XRTM might complain if the first phase function
          ! expansion coefficient is not EXACTLY 1.0d0
          ! If the actual value is < 1d-10 away from 1.0d0,
          ! we can somewhat safely assume that the calculation
          ! was done correctly, but floating point business
          ! ends up giving us a value away from 1.0

          where(abs(this_coef(1,1,:) - 1.0d0) < 1d-10) this_coef(1,1,:) = 1.0d0


          call xrtm_set_coef_n_f90(xrtm, n_coef, this_coef, xrtm_error)
#ifdef DEBUG
          if (xrtm_error /= 0) then
             call logger%error(fname, "Error calling xrtm_set_coef_n_f90")
             return
          end if
#endif
          ! Plug in the phase function expansion moment derivatives
          call xrtm_set_coef_l_nn_f90(xrtm, this_lcoef, xrtm_error)
#ifdef DEBUG
          if (xrtm_error /= 0) then
             call logger%error(fname, "Error calling xrtm_set_coef_l_nn_f90")
             return
          end if
#endif
       end if

       ! Plug in surface property
       call xrtm_set_kernel_ampfac_f90(xrtm, 0, albedo(i), xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampfac_f90")
          return
       end if
#endif
       ! Plug in the surface derivatives
       call xrtm_set_kernel_ampfac_l_n_f90(xrtm, 0, single_lsurf, xrtm_error)
#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampflac_l_n_f90")
          return
       end if
#endif
       ! Call this function to update the layer structure in XRTM,
       ! according to the manual it must only change if a linearized value
       ! changed from non-zero to zero (or the other way round). Hence,
       ! we only call it at the first wavelength.
       if (i == 1) then
          call xrtm_update_varied_layers_f90(xrtm, xrtm_error)
#ifdef DEBUG
          if (xrtm_error /= 0) then
             call logger%error(fname, "Error calling xrtm_update_varied_layers_f90")
             return
          end if
#endif
       endif

       ! NOTE
       ! Absolutely NEVER use "call cpu_time()" in an OpenMP environment
       ! as your performance will completely tank with more than a few cores.
       ! I suspect it's some scheduling issue where a thread has to wait
       ! until the no other thread uses this function.
       ! So the solution here is to compile with "omp_get_wtime" whenever
       ! we compile with OpenMP, and use "cpu_time" otherwise.

       cpu_start2 = get_cpu_time()
       ! Calculate TOA radiance!
       call xrtm_radiance_f90( &
            xrtm, & ! XRTM object
            xrtm_solver, & ! XRTM solver bitmask
            1 , & ! Number of output azimuths
            out_phis, & ! Output azimuths
            I_p, & ! Upward radiances,
            I_m, & ! Downward radiances
            K_p, & ! Upward jacobians
            K_m, & ! Downward jacobians
            xrtm_error)

       cpu_end2 = get_cpu_time()

       pure_XRTM_duration = pure_XRTM_duration + cpu_end2 - cpu_start2

#ifdef DEBUG
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_radiance_f90")
          return
       end if

       if (any(ieee_is_nan(I_p))) then
          write(tmp_str, '(A, G0.1)') "XRTM produced a NaN for radiances " &
               // "at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if

       if (any(ieee_is_nan(K_p))) then
          write(tmp_str, '(A, G0.1)') "XRTM produced a NaN for weighting functions " &
               // "at wavelength index: ", i
          call logger%error(fname, trim(tmp_str))
          return
       end if
#endif


       ! Store radiances into wavelength-indexed array
       radiance(i,:) = radiance(i,:) + I_p(:,1,1,1)

       ! Store derivatives into wavelength-indexed array
       derivs(i,:,:) = derivs(i,:,:) + K_p(:,1,1,:,1)


       ! Format string for derivative debug output
       ! .. but only go into this branch when debug is requested
       !if ((i == 1) .and. (CS_general%loglevel <= 10)) then
       !call logger%debug(fname, "--- First wavelength ----------------------------------------------")

       ! Write XRTM inputs
       !call logger%debug(fname, "-------------------------------------------------------------------")
       !call logger%debug(fname, "Layer, input optical depth, input single scatter albedo, # of coefs")
       !do m = 1, n_layers
       !   write(tmp_str, '(I5, ES15.5, ES15.5, I6)') m, tau(m), omega(m), n_coef(m)
       !   call logger%debug(fname, trim(tmp_str))
       !end do

       !call logger%debug(fname, "-------------------------------------------------------")
       !call logger%debug(fname, "Layer, input linearized optical depth inputs")
       !do m = 1, n_layers
       !   write(fmt_str, '(A, G0.1, A)') "(I5, ", size(derivs, 3), "ES15.5)"
       !   write(tmp_str, trim(fmt_str)) m, single_ltau(:,m)
       !   call logger%debug(fname, trim(tmp_str))
       !end do

       !call logger%debug(fname, "----------------------------------------------------")
       !call logger%debug(fname, "Layer, input linearized single scatter albedo inputs")
       !do m = 1, n_layers
       !   write(fmt_str, '(A, G0.1, A)') "(I5, ", size(derivs, 3), "ES15.5)"
       !   write(tmp_str, trim(fmt_str)) m, single_lomega(:,m)
       !   call logger%debug(fname, trim(tmp_str))
       !end do

       !call logger%debug(fname, "------------------------------------------------")
       !call logger%debug(fname, "Layer, first few phasefunction expansion coeffs.")
       !do m = 1, n_layers
       !   do q = 1, size(this_coef, 2)
       !      ! (pfmom, elem, n_layers)
       !      write(tmp_str, '(A, G0.1, I5, 3ES15.5)') "Element #", q, m, this_coef(1:3, q, m)
       !      call logger%debug(fname, trim(tmp_str))
       !   end do
       !end do

       !call logger%debug(fname, "-----------------------------------------------------------------------")
       !call logger%debug(fname, "Layer, derivative, first few linearized phasefunction expansion coeffs.")
       !do m = 1, n_layers
       !   do k = 1, n_derivs
       !      ! (pfmom, elem, n_derivs, n_layers)
       !      write(tmp_str, '(I5, I5, 3ES15.5)') m, k, this_lcoef(1:3, 1, k, m)
       !      call logger%debug(fname, trim(tmp_str))
       !   end do
       !end do
       !call logger%debug(fname, "-----------------------------------------------------------------------")

       !write(fmt_str, '(A,G0.1,A)') "(A, G0.1, A, ", size(derivs, 3), "ES15.5)"
       !do q = 1, n_stokes
       !   write(tmp_str, '(A, G0.1, A, ES20.10)') "XRTM radiance, Stokes # ", q, "          : ", I_p(q,1,1,1)
       !   call logger%debug(fname, trim(tmp_str))
       !   write(tmp_str, '(A, G0.1, A, ES20.10)') "Cumulative XRTM radiance, Stokes #", q, ": ", radiance(i,q)
       !   call logger%debug(fname, trim(tmp_str))
       !end do

       !do q = 1, n_stokes
       !   write(tmp_str, trim(fmt_str)) "XRTM weighting functions, Stokes #", q, "           : ", K_p(q,1,1,:,1)
       !   call logger%debug(fname, trim(tmp_str))
       !   write(tmp_str, trim(fmt_str)) "Cumulative XRTM weighting functions, Stokes #", q, ": ", derivs(i,q,:)
       !   call logger%debug(fname, trim(tmp_str))
       !end do


       !end if

#ifdef DEBUG
       if (N_spec > 10) then
          if (mod(i, N_spec/10) == 0) then
             write(tmp_str, '(A, F6.2, A, G0.1, A, G0.1, A)') &
                  "XRTM calls (", 100.0 * i/N_spec, "%, ", i, "/", N_spec, ")"
             call logger%debug(fname, trim(tmp_str))
          end if
       end if
#endif

    end do

    cpu_end = get_cpu_time()

    write(tmp_str, '(A, F7.3, A)') "Pure XRTM radiance calculations: ", pure_XRTM_duration, " sec."
    call logger%debug(fname, trim(tmp_str))
    write(tmp_str, '(A, F7.3, A)') "XRTM monochr. calc. with overhead: ", cpu_end - cpu_start, " sec."
    call logger%debug(fname, trim(tmp_str))

    success = .true.

  end subroutine calculate_XRTM_radiance

  subroutine calculate_XRTM_scale_AK_corr( &
       solar_hires, &
       dI_dtaugas, &
       stokes_coeffs, &
       scn, &
       SV, &
       gain_matrix, &
       ILS_delta_lambda, &
       ILS_relative_response, &
       dispersion, &
       psurf, &
       num_active_levels, &
       N_stokes, &
       N_spec, &
       idx_start, &
       idx_stop, &
       col_AK &
       )

    implicit none
    double precision, intent(in) :: solar_hires(:)
    double precision, intent(in) :: dI_dtaugas(:,:,:)
    double precision, intent(in) :: stokes_coeffs(:)
    type(scene), intent(in) :: scn
    type(statevector), intent(in) :: SV
    double precision, intent(in) :: gain_matrix(:,:)
    real, intent(in) :: ILS_delta_lambda(:,:)
    real, intent(in) :: ILS_relative_response(:,:)
    double precision, intent(in) :: dispersion(:)
    double precision, intent(in) :: psurf
    integer, intent(in) :: num_active_levels
    integer, intent(in) :: N_stokes
    integer, intent(in) :: N_spec
    integer, intent(in) :: idx_start(:)
    integer, intent(in) :: idx_stop(:)
    double precision, intent(inout) :: col_AK(:,:)


    logical :: ILS_success

    double precision, allocatable :: dI_dVMR_stokes(:,:)
    double precision, allocatable :: dI_dVMR(:)
    double precision, allocatable :: dI_dVMR_conv(:,:)
    double precision :: pwgts(num_active_levels)
    double precision :: prior_VMR(num_active_levels), this_VMR(num_active_levels)
    double precision :: s_bar(num_active_levels)

    double precision, allocatable :: AK_profile(:,:), AK_profile_total(:)
    double precision, allocatable :: tmp_v1(:), tmp_v2(:)

    integer :: idx1, idx2

    integer :: N_gas, nlay, N_hires
    integer :: i_gas, i_SV
    integer :: i, j, q

    N_hires = size(scn%op%wl)
    N_gas = size(scn%op%gas_tau, 3)
    nlay = scn%num_active_levels - 1

    allocate(dI_dVMR(N_hires))
    allocate(dI_dVMR_stokes(N_hires, N_stokes))
    allocate(dI_dVMR_conv(N_spec, num_active_levels))
    allocate(tmp_v1(size(SV%svap)))
    allocate(tmp_v2(num_active_levels))
    allocate(AK_profile(SV%num_gas, num_active_levels))
    allocate(AK_profile_total(num_active_levels))

    ! dI_dtaugas: (N_hires, N_gas, N_stokes)

    ! Calculate this for every retrieved gas
    do i_gas=1, N_gas

       ! Is this gas found in ANY of the state vector elements?
       ! If not - skip this gas
       if (count(SV%gas_idx_lookup(:) == i_gas) == 0) cycle

       dI_dVMR_stokes(:,:) = 0.0d0
       dI_dVMR_conv(:,:) = 0.0d0
       pwgts(:) = 0.0d0
       prior_VMR(:) = 0.0d0
       this_VMR(:) = 0.0d0

       do i=1, num_active_levels

          ! THIS IS WRONG
          ! dI_dgas must be replaced by dI_dTau
          ! which unfortunately we do not get from XRTM
          ! at the moment.
          ! The main reason is that the number of weighting
          ! functions to be calculated would increase so much
          ! and slow down the entire algorithm.

          if (i == 1) then
             ! TOA layer
             do j = 1, N_hires
                dI_dVMR_stokes(j,:) = solar_hires(j) * dI_dtaugas(j, i, :) &
                     * scn%op%gas_tau_dvmr(1,j,i,i_gas)
             end do
          else if (i == num_active_levels) then
             ! BOA layer
             do j = 1, N_hires
                dI_dVMR_stokes(j,:) = solar_hires(j) * dI_dtaugas(j, i-1, :) &
                     * scn%op%gas_tau_dvmr(2,j,i-1,i_gas)
             end do
          else
             ! everything in between
             do j = 1, N_hires
                dI_dVMR_stokes(j,:) = solar_hires(j) * dI_dtaugas(j, i, :) &
                     * (scn%op%gas_tau_dvmr(2,j,i-1,i_gas) &
                     + scn%op%gas_tau_dvmr(1,j,i,i_gas))
             end do
          end if

          ! Sum up the dI_dVMR stokes components to get
          dI_dVMR(:) = 0.0d0
          do q=1, N_stokes
             dI_dVMR(:) = dI_dVMR(:) + dI_dVMR_stokes(:, q) * stokes_coeffs(q)
          end do

          call oco_type_convolution(scn%op%wl, &
               dI_dVMR(:), &
               ILS_delta_lambda(:,:), &
               ILS_relative_response(:,:), &
               dispersion(:), &
               dI_dVMR_conv(:,i), &
               ILS_success)

       end do ! End level loop

       ! Compute / adjust the VMRs given the prior scale factor (not necessarily 1.0)
       ! and the region of the column in which the sub-column is retrieved
       prior_VMR = scn%atm%gas_vmr(1:num_active_levels, i_gas)
       this_VMR = scn%atm%gas_vmr(1:num_active_levels, i_gas)

       do i_SV=1, SV%num_gas
          if (i_gas /= SV%gas_idx_lookup(i_SV)) cycle

          idx1 = max(idx_start(i_SV), 1)
          idx2 = min(idx_stop(i_SV), num_active_levels)

          prior_VMR(idx1:idx2) = prior_VMR(idx1:idx2) &
               * SV%svap(SV%idx_gas(i_SV, 1))
          this_VMR(idx1:idx2) = this_VMR(idx1:idx2) &
               * SV%svsv(SV%idx_gas(i_SV, 1))
       end do

       AK_profile_total(:) = 0.0d0
       do i_SV=1, SV%num_gas
          if (i_gas /= SV%gas_idx_lookup(i_SV)) cycle

          idx1 = max(idx_start(i_SV), 1)
          idx2 = min(idx_stop(i_SV), num_active_levels)

          s_bar(:) = 0.0d0
          s_bar(idx1:idx2) = 1.0d0

          tmp_v2(:) = 0.0d0
          do i=1, num_active_levels
             ! G dot dI/dVMR is a vector of length(SV)
             tmp_v1(:) = matmul(gain_matrix, dI_dVMR_conv(:, i))
             ! Now multiply by (h^T dot [prior_VMR * s_bar])
             AK_profile(i_SV, i) = tmp_v1(SV%idx_gas(i_SV, 1)) * dot_product(scn%atm%pwgts, prior_VMR * s_bar)
             AK_profile_total(i) = AK_profile_total(i) + AK_profile(i_SV, i)

          end do

       end do

       col_AK(i_gas, 1:num_active_levels) = AK_profile_total(:)

    end do ! End gas loop

  end subroutine calculate_XRTM_scale_AK_corr

end module XRTM_mod
