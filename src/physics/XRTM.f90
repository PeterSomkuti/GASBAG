!> @brief Module to house the interfaces and calculations for the XRTM RT model
!> @file XRTM.f90
!> @author Peter Somkuti

module XRTM_mod

  ! User modules
  use control_mod, only: CS_window
  use math_utils_mod
  use Rayleigh_mod, only : calculate_rayleigh_scatt_matrix, calculate_rayleigh_depolf
  use statevector_mod
  use aerosols_mod
  use physical_model_addon_mod, only : scene
  
  ! Third-party modules
  use xrtm_int_f90
  use stringifor
  use logger_mod, only: logger => master_logger

  implicit none

  public :: setup_XRTM, create_XRTM, solve_RT_problem_XRTM, calculate_XRTM_radiance

contains


  !> @brief "Translates" configuration values into XRTM settings
  !> @param xrtm_options_string Array of strings containing XRTM options
  !> @param xrtm_solvers_string Array of strings containing XRTM solvers
  !> @param do_polarization Run RT in vector mode?
  !> @param xrtm_options Integer bitmask for XRTM options
  !> @param xrtm_solvers Integer bitmask for XRTM solvers
  !> @param xrtm_kernels Integer array specifying which surface kernels
  !> @param success Did any error occur?
  subroutine setup_XRTM(xrtm_options_string, xrtm_solvers_string, &
       do_polarization, &
       xrtm_options, xrtm_solvers, &
       xrtm_separate_solvers, xrtm_kernels, success)

    implicit none
    type(string), allocatable, intent(in) :: xrtm_options_string(:)
    type(string), allocatable, intent(in) :: xrtm_solvers_string(:)
    logical, intent(in) :: do_polarization
    integer, intent(inout) :: xrtm_options
    integer, intent(inout) :: xrtm_solvers
    integer, allocatable, intent(inout) :: xrtm_separate_solvers(:)
    integer, intent(inout) :: xrtm_kernels(:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "setup_XRTM"

    type(string) :: tmp_str
    integer :: N, i, j
    integer :: xrtm_error
    integer :: num_solvers

    success = .false.

    ! Initialize XRTM options and solvers with 0
    xrtm_options = 0
    xrtm_solvers = 0


    ! Surface kernels
    xrtm_kernels(1) = XRTM_KERNEL_LAMBERTIAN

    num_solvers = 0
    ! Populate solvers bit field, depending on requested options
    if (allocated(xrtm_solvers_string)) then
       do i=1, size(xrtm_solvers_string)
          ! Grab local copy of string
          tmp_str = xrtm_solvers_string(i)

          if (tmp_str == "TWO_OS") then
             ! Vijay-like 2OS, gives you SECOND ORDER ONLY
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_TWO_OS)
             num_solvers = num_solvers + 1
             call logger%debug(fname, "Using XRTM in 2-orders-of-scattering mode")
          else if (tmp_str == "SINGLE") then
             ! Single scattering only
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_SINGLE)
             num_solvers = num_solvers + 1
             call logger%debug(fname, "Using XRTM in single-scattering mode")
          else if (tmp_str == "EIG_BVP") then
             ! Quadrature (LIDORT-like)
             call logger%debug(fname, "Using XRTM in discrete-ordinate mode")
             num_solvers = num_solvers + 1
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_EIG_BVP)
          else if (tmp_str == "SOS") then
             ! Pade approximation and adding
             call logger%debug(fname, "Using XRTM in successive-orders-of-scattering mode")
             call logger%debug(fname, "-rt_streams- is used to infer scattering order")
             num_solvers = num_solvers + 1
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_SOS)
          else if (tmp_str == "PADE_ADD") then
             ! Pade approximation and adding
             call logger%debug(fname, "Using XRTM in Pade-adding mode")
             num_solvers = num_solvers + 1
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_PADE_ADD)
          else
             call logger%error(fname, "XRTM solver option is not implemented, and will be ignored: " &
                  // tmp_str%chars())
          end if
       end do
    end if

    ! XRTM solvers come as a 32-bit mask, and during the set-up stage,
    ! XRTM ingests the inclusive-or combined one such that we can
    ! then call every solver individually. This little piece of
    ! code simply creates an array with integers, where each
    ! array element is one specific solver that we can hand
    ! off to "xrtm_radiance".
    ! Then feed each separate solver (integer value) into a separate
    ! element of the "separate_solvers" array.
    allocate(xrtm_separate_solvers(num_solvers))
    i = 1
    do j=0, 31
       if (iand(xrtm_solvers, int(2**j)) /= 0) then
          xrtm_separate_solvers(i) = int(2**j)
          i = i + 1
       end if
    end do

    ! Do weighting functions and pseudo-spherical geometry
    ! by default. Also, only return upwelling radiances.
    xrtm_options = ior(XRTM_OPTION_CALC_DERIVS, XRTM_OPTION_PSA)
    xrtm_options = ior(xrtm_options, XRTM_OPTION_UPWELLING_OUTPUT)

    ! If polarization is requested, we run
    ! the RT models in vector mode
    if (do_polarization) then
       call logger%debug(fname, "Using XRTM in vector mode.")
       xrtm_options = ior(xrtm_options, XRTM_OPTION_VECTOR)
    else
       call logger%debug(fname, "Using XRTM in scalar mode.")
    end if

    xrtm_options = ior(xrtm_options, XRTM_OPTION_N_T_TMS)
    xrtm_options = ior(xrtm_options, XRTM_OPTION_DELTA_M)

    success = .true.

  end subroutine setup_XRTM


  subroutine create_XRTM(xrtm, xrtm_options, xrtm_solvers, &
       max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, &
       n_kernel_quad, xrtm_kernels, success)
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
    integer :: xrtm_error

    success = .false.

    call xrtm_create_f90(xrtm, &  ! XRTM object
         xrtm_options, & ! XRTM option bitmask
         xrtm_solvers, & ! XRTM solver bitmask
         max_coef, & ! Number of phase matrix Legendre coefficients
         n_quad, & ! Number of quadrature streams for some solvers
         n_stokes, & ! Number of stokes coefficients
         n_derivs, & ! Number of derivatives to calculate
         n_layers, & ! Number of layers in model atmosphere
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


  subroutine solve_RT_problem_XRTM(window, SV, scn, n_stokes, s_start, s_stop, &
       radiance_calc_work_hi_stokes, dI_dgas, dI_dsurf, dI_dTemp, dI_dAOD, &
       xrtm_success)

    type(CS_window), intent(in) :: window
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    integer, intent(in) :: n_stokes
    integer, intent(in) :: s_start(:)
    integer, intent(in) :: s_stop(:)
    double precision, intent(inout) :: radiance_calc_work_hi_stokes(:,:)
    double precision, allocatable, intent(inout) :: dI_dgas(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dsurf(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dTemp(:,:)
    double precision, allocatable, intent(inout) :: dI_dAOD(:,:,:)

    logical, intent(inout) :: xrtm_success

    ! Local variables

    ! Function name for logger
    character(len=*), parameter :: fname = "solve_RT_problem_XRTM"
    character(len=999) :: tmp_str

    double precision, allocatable :: weighting_functions(:,:,:)

    integer :: num_act_lev
    integer :: num_lay
    ! XRTM Radiative Transfer model handler
    type(xrtm_type) :: xrtm
    ! Actual options variable
    integer :: xrtm_options
    ! Actual solvers variable
    integer :: xrtm_solvers
    ! Solvers separated into individual one's
    integer, allocatable :: xrtm_separate_solvers(:)
    ! XRTM surface kernels
    integer :: xrtm_kernels(1)
    ! XRTM error variable
    integer :: xrtm_error
    ! How many derivatives does XRTM need to calculate?
    integer :: xrtm_n_derivs
    !
    integer :: xrtm_streams

    double precision, allocatable :: coef(:,:,:,:)
    double precision, allocatable :: lcoef(:,:,:,:,:)

    double precision, allocatable :: ltau(:,:,:)
    double precision, allocatable :: lomega(:,:,:)
    double precision, allocatable :: lsurf(:,:)
    !
    double precision :: cpu_start, cpu_stop
    integer :: i,j
    integer :: aer_idx

    ! Initialize success variable
    xrtm_success = .false.

    ! Grab the number of active levels for easy access
    num_act_lev = scn%num_active_levels
    num_lay = num_act_lev - 1

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

    if (window%RT_strategy == "") then
       call logger%fatal(fname, "Must have an RT strategy for XRTM!")
       stop 1
    end if

    if (window%RT_strategy%lower() == "monochromatic") then

       call logger%debug(fname, "Using monochromatic RT strategy.")

       ! --------------------------------------------------------------------
       ! Monochromatic RT
       ! --------------------------------------------------------------------
       !
       ! This is the most straight-forward option. Using the user-supplied
       ! RT model and options, we run through the entire band, point for point.
       ! --------------------------------------------------------------------

       
       ! This function "translates" our verbose configuration file options
       ! into proper XRTM language and sets the corresponding options.

       call setup_XRTM( &
            window%xrtm_options, &
            window%xrtm_solvers, &
            window%do_polarization, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_separate_solvers, &
            xrtm_kernels, &
            xrtm_success)

          if (.not. xrtm_success) then
             call logger%error(fname, "Call to setup_XRTM unsuccessful.")
             return
          end if

          xrtm_n_derivs = SV%num_gas + SV%num_temp + SV%num_albedo + SV%num_aerosol_aod

          ! --------------------------------------
          ! Derive the number of quadrature points
          !
          ! If more than one was supplied by the user,
          ! simply grab the first value and ignore the
          ! rest for monochromatic calculations. At least
          ! one value has to be present. Of course this
          ! option is only required for RT solvers that
          ! require some notion of "streams", such as
          ! BVP.
          !
          ! --------------------------------------

          if ( &
               (iand(xrtm_solvers, XRTM_SOLVER_EIG_BVP) /= 0) .or. &
               (iand(xrtm_solvers, XRTM_SOLVER_TWO_OS) /= 0) &
               ) then

             ! Is it even allocated?
             if (.not. allocated(window%RT_streams)) then
                call logger%fatal(fname, "You MUST supply a -rt_streams- option!")
                stop 1
             end if

             ! And if allocated, does it have the right size?
             if (size(window%RT_streams) < 1) then
                call logger%fatal(fname, "You need to supply at least one value for -rt_streams-")
                stop 1
             end if

             ! And then finally it needs to have the right value
             if (window%RT_streams(1) < 2) then
                write(tmp_str, '(A, G0.1)') "Need at least 2 streams, but you said: ", window%RT_streams(1)
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if

             xrtm_streams = window%RT_streams(1) / 2

          else

             xrtm_streams = 1
             
          end if


          ! Precompute all  phase function coefficients and derivatives
          call cpu_time(cpu_start)
          call precompute_all_coef(scn, SV, n_stokes, xrtm_n_derivs, &
               window%constant_coef, coef, lcoef)
          call cpu_time(cpu_stop)

          call create_XRTM( &
               xrtm, & ! XRTM handler
               xrtm_options, & ! XRTM options bitmask
               xrtm_solvers, & ! XRTM solvers combined bitmask
               max(3, scn%max_pfmom), & ! Max coef
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

          ! Precompute the derivatives, so we do not need to deal with it
          ! inside the radiance calculation code.
          allocate(ltau(size(scn%op%wl), xrtm_n_derivs, num_lay))
          allocate(lomega(size(scn%op%wl), xrtm_n_derivs, num_lay))
          allocate(lsurf(size(scn%op%wl), xrtm_n_derivs))

          allocate(weighting_functions(size(scn%op%wl), n_stokes, xrtm_n_derivs))


          ltau(:,:,:) = 0.0d0
          lomega(:,:,:) = 0.0d0
          lsurf(:,:) = 0.0d0

          do j=1, SV%num_gas
             ltau(:, j, s_start(j):s_stop(j)-1) = &
                  scn%op%gas_tau(:, s_start(j):s_stop(j)-1, SV%gas_idx_lookup(j)) / SV%svsv(SV%idx_gas(j,1))
             lomega(:, j, s_start(j):s_stop(j)-1) = &
                  -scn%op%layer_omega(:,s_start(j):s_stop(j)-1) / scn%op%layer_tau(:,s_start(j):s_stop(j)-1) &
                  * ltau(:, j, s_start(j):s_stop(j)-1)
          end do

          if (SV%num_temp > 0) then
             j = SV%num_gas + 1

             ltau(:,j,:) = sum(scn%op%gas_tau_dtemp(:,1:num_lay,:), dim=3)
             lomega(:,j,:) = -scn%op%layer_omega(:,:) / scn%op%layer_tau(:,:) * ltau(:,j,:)

          end if

          do j=1, SV%num_albedo
             i = SV%num_gas + SV%num_temp + j

             lsurf(:,i) = 1.0d0
          end do

          do j=1, SV%num_aerosol_aod
             i = SV%num_gas + SV%num_temp + SV%num_albedo + j
             aer_idx = SV%aerosol_idx_lookup(j)

             ltau(:,i,:) = scn%op%aer_ext_tau(:,:,aer_idx) / scn%op%reference_aod(aer_idx)
             lomega(:,i,:) = ltau(:,i,:) * (1.0d0 - scn%op%layer_omega(:,:) / scn%op%layer_tau(:,:))
          end do

          ! ------------------------------------------
          ! Run the monochromatic radiance calculation
          ! ------------------------------------------


          call calculate_XRTM_radiance( &
               xrtm, & ! XRTM handler
               xrtm_separate_solvers, & ! separate XRTM solvers
               SV, & ! State vector structure - useful for decoding SV positions
               scn%op%wl, & ! per pixel wavelength
               scn%SZA, & ! solar zenith angle
               scn%VZA, & ! viewing zenith angle
               scn%SAA, & ! solar azimuth angle
               scn%VAA, & ! viewing azimuth angle
               scn%atm%altitude_levels, & ! per-level altitude
               scn%op%albedo, & ! Surface albedo
               scn%op%layer_tau, & ! per-wl per-layer optical depths
               scn%op%layer_omega, & ! per-wl per-layer SSA
               ltau, &
               lomega, &
               lsurf, &
               coef, lcoef, &
               n_stokes, & ! Number of Stokes parameters to calculate
               xrtm_n_derivs, & ! Number of derivatives to be calculated
               num_lay, & ! Number of atmospheric layers
               radiance_calc_work_hi_stokes, &
               weighting_functions, &
               xrtm_success)

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

          if (.not. xrtm_success) then
             call logger%error(fname, "Call to calculate_XRTM_radiance unsuccessful.")
             return
          end if


          ! Recover the Jacobians from the XRTM container
          ! Store gas subcolumn derivatives
          do j=1, SV%num_gas
             dI_dgas(:,j,:) = dI_dgas(:,j,:) + weighting_functions(:,:,j)
          end do

          ! Store the temperature offset Jacobian if needed
          if (SV%num_temp > 0) then
             dI_dTemp(:,:) = dI_dTemp(:,:) + weighting_functions(:,:,SV%num_gas + 1)
          end if

          ! Store surface jacobian
          do j=1, SV%num_albedo
             dI_dsurf(:,j,:) = weighting_functions(:,:,SV%num_gas + SV%num_temp + j)
          end do

          do j=1, SV%num_aerosol_aod
             dI_dAOD(:,j,:) = weighting_functions(:,:,SV%num_gas + SV%num_temp + SV%num_albedo + j)
          end do

    else

       call logger%fatal(fname, "RT strategy: '" // window%RT_strategy%chars() &
            // "' not recognized. Valid options are..")

       call logger%fatal(fname, "(1) monochromatic")
       !call logger%fatal(fname, "(2) PCA (work in progress)")
       stop 1

    end if


  end subroutine solve_RT_problem_XRTM





  subroutine calculate_XRTM_radiance(xrtm, xrtm_separate_solvers, &
       SV, &
       wavelengths, SZA, VZA, SAA, VAA, &
       altitude_levels, albedo, &
       layer_tau, &
       layer_omega, &
       ltau, &
       lomega, &
       lsurf, &
       coef, &
       lcoef, &
       n_stokes, n_derivs, n_layers, &
       radiance, &
       derivs, &
       success)

    implicit none
    type(xrtm_type), intent(inout) :: xrtm
    integer, intent(in) :: xrtm_separate_solvers(:)
    type(statevector), intent(in) :: SV
    double precision, intent(in) :: wavelengths(:)
    double precision, intent(in) :: SZA, VZA, SAA, VAA, albedo(:)
    double precision, intent(in) :: altitude_levels(:)
    double precision, intent(in) :: layer_tau(:,:)
    double precision, intent(in) :: layer_omega(:,:)
    double precision, intent(in) :: ltau(:,:,:)
    double precision, intent(in) :: lomega(:,:,:)
    double precision, intent(in) :: lsurf(:,:)
    double precision, allocatable, intent(in) :: coef(:,:,:,:) !(spectral, pfmom, elem, n_layers)
    double precision, allocatable, intent(in) :: lcoef(:,:,:,:,:)  !(spectral, pfmom, elem, n_derivs, n_layers)

    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers

    double precision, intent(inout) :: radiance(:,:)
    double precision, intent(inout) :: derivs(:,:,:)
    
    !double precision, allocatable, intent(inout) :: dI_dgas(:,:,:)
    !double precision, allocatable, intent(inout) :: dI_dsurf(:,:,:)
    !double precision, allocatable, intent(inout) :: dI_dTemp(:,:)
    !double precision, allocatable, intent(inout) :: dI_dAOD(:,:,:)

    logical, intent(inout) :: success

    ! Local
    integer :: xrtm_options
    integer :: xrtm_solvers

    double precision, allocatable :: I_p(:,:,:,:)
    double precision, allocatable :: I_m(:,:,:,:)
    double precision, allocatable :: K_p(:,:,:,:,:)
    double precision, allocatable :: K_m(:,:,:,:,:)

    character(len=*), parameter :: fname = "calculate_XRTM_radiance"
    character(len=999) :: tmp_str
    integer :: xrtm_error

    integer :: num_solvers
    integer, allocatable :: separate_solvers(:)

    double precision :: out_thetas(1)
    double precision :: out_phis(1,1)
    double precision, allocatable :: ray_coef(:,:)
    integer :: out_levels(1)
    integer :: n_coef(n_layers)
    double precision :: single_ltau(n_derivs, n_layers)
    double precision :: single_lomega(n_derivs, n_layers)
    double precision :: single_lsurf(n_derivs)
    
    double precision :: tau(n_layers), omega(n_layers), omega_tmp(n_layers)
    double precision, allocatable :: this_coef(:,:,:), this_lcoef(:,:,:,:)
    integer :: N_spec
    integer :: max_coefs
    integer :: i,j,k,l,m

    integer :: funit

    double precision :: cpu_start, cpu_end

    success = .false.

    xrtm_options = xrtm_get_options_f90(xrtm)
    xrtm_solvers = xrtm_get_solvers_f90(xrtm)

    N_spec = size(wavelengths)

    radiance(:,:) = 0.0d0
    derivs(:,:,:) = 0.0d0
    !dI_dgas(:,:,:) = 0.0d0
    
    !if (allocated(dI_dTemp)) dI_dTemp(:,:) = 0.0d0
    !if (allocated(dI_dsurf)) dI_dsurf(:,:,:) = 0.0d0
    !if (allocated(dI_dAOD)) dI_dAOD(:,:,:) = 0.0d0

    out_thetas(1) = VZA
    out_phis(1,1) = VAA
    out_levels(1) = 0

    single_ltau(:,:) = 0.0d0
    single_lomega(:,:) = 0.0d0
    single_lsurf(:) = 0.0d0

    ! Derivative dI/dsurf - for the time being only Albedo
    !do i = 1, SV%num_albedo
    !   lsurf(SV%num_gas + SV%num_temp + i) = 1.0d0
    !end do
    
    !if (n_stokes > 1) then
    !   allocate(ray_coef(3, 6))
       !allocate(coef(3, 6, n_layers))
    !else
    !   allocate(ray_coef(3, 1))
       !allocate(coef(3, 1, n_layers))
    !endif

    n_coef(:) = size(coef, 2)



    ! Allocate output containers - these will not change with wavelength, hence
    ! we only need to do this once per microwindow
    allocate(I_p(n_stokes, 1, 1, 1))
    allocate(I_m(n_stokes, 1, 1, 1))
    allocate(K_p(n_stokes, 1, 1, n_derivs, 1))
    allocate(K_m(n_stokes, 1, 1, n_derivs, 1))

    allocate(this_coef(size(coef, 1), size(coef, 2), size(coef, 3)))
    allocate(this_lcoef(size(lcoef, 1), size(lcoef, 2), size(lcoef, 3), size(lcoef, 4)))


    ! Some general set-up's that are not expected to change, and can
    ! be safely hard-coded.
    call xrtm_set_fourier_tol_f90(xrtm, .0001d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_fourier_tol_f90")
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

    if (iand(xrtm_solvers, XRTM_SOLVER_SOS) /= 0) then

       call xrtm_set_sos_params_f90(xrtm, 4, 10.0d0, 0.01d0, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_sos_params_f90")
          return
       endif
    end if

    if (iand(xrtm_solvers, XRTM_SOLVER_PADE_ADD) /= 0) then
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
    ! -----------------------------------------------------------------------------------

    call cpu_time(cpu_start)
    do i=1, N_spec

       if (mod(i, N_spec/10) == 0) then
          write(tmp_str, '(A, F6.2, A, G0.1, A, G0.1, A)') &
               "XRTM calls (", 100.0 * i/N_spec, "%, ", i, "/", N_spec, ")"
          call logger%debug(fname, trim(tmp_str))
       end if

       ! TOTAL atmospheric optical properties - these go into the
       ! RT code for the forward calculation.
       tau(:) = layer_tau(i, :)
       omega(:) = layer_omega(i, :) !(ray_tau(i, :) + sum(aer_tau(i,:,:), dim=2))/ tau(:)

       ! Plug in surface property
       call xrtm_set_kernel_ampfac_f90(xrtm, 0, albedo(i), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampfac_f90")
          return
       end if

       ! Plug in optical depth
       call xrtm_set_ltau_n_f90(xrtm, tau, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_n_f90")
          return
       end if

       ! Plug in the layer optical depth derivatives
       single_ltau(:,:) = ltau(i,:,:)
       call xrtm_set_ltau_l_nn_f90(xrtm, single_ltau, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_nn_f90")
          return
       end if

       ! Plug in the single-scatter albedo
       call xrtm_set_omega_n_f90(xrtm, omega, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_n_f90")
          return
       end if

       ! Plug in the single-scatter albedo derivatives
       single_lomega(:,:) = lomega(i,:,:)
       call xrtm_set_omega_l_nn_f90(xrtm, single_lomega, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_l_nn_f90")
          return
       end if

       ! Plug in the layer scattering matrix
       if (size(coef, 4) > 1) then
          this_coef(:,:,:) = coef(:,:,:,i)
          this_lcoef(:,:,:,:) = lcoef(:,:,:,:,i)
       else
          this_coef(:,:,:) = coef(:,:,:,1)
          this_lcoef(:,:,:,:) = lcoef(:,:,:,:,1)
       end if

       call xrtm_set_coef_n_f90(xrtm, n_coef, this_coef, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_coef_n_f90")
          return
       end if



       !call xrtm_set_coef_l_nn_f90(xrtm, this_lcoef, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_coef_l_nn_f90")
          return
       end if

       ! Plug in the surface derivatives
       single_lsurf(:) = lsurf(i,:)
       call xrtm_set_kernel_ampfac_l_n_f90(xrtm, 0, single_lsurf(:), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampflac_l_n_f90")
          return
       end if

       ! Call this function for whatever reason
       call xrtm_update_varied_layers_f90(xrtm, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_update_varied_layers_f90")
          return
       end if


       ! XRTM has been initialized with whatever number of solvers are stored in
       ! "xrtm_solvers", however only one is executed at a time (ask Greg?).
       ! Thus, we need to loop over the possible bitmask positions.
       !


       do l=1, size(xrtm_separate_solvers)
          ! Calculate TOA radiance!
          call xrtm_radiance_f90(xrtm, & ! XRTM object
               xrtm_separate_solvers(l), & ! XRTM solver bitmask
               1 , & ! Number of output azimuths
               out_phis, & ! Output azimuths
               I_p, & ! Upward radiances,
               I_m, & ! Downward radiances
               K_p, & ! Upward jacobians
               K_m, & ! Downward jacobians
               xrtm_error)

          if (xrtm_error /= 0) then
             call logger%error(fname, "Error calling xrtm_radiance_f90")
             return
          end if

          ! Store radiances
          radiance(i,:) = radiance(i,:) + I_p(:,1,1,1)

          ! Store derivatives
          derivs(i,:,:) = derivs(i,:,:) + K_p(:,1,1,:,1)

       end do

    end do

    call cpu_time(cpu_end)
    write(tmp_str, '(A, F7.3, A)') "XRTM monochromatic calculations: ", cpu_end - cpu_start, " sec"
    call logger%debug(fname, trim(tmp_str))

    success = .true.

  end subroutine calculate_XRTM_radiance

end module XRTM_mod
