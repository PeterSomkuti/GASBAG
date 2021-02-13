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
  subroutine setup_XRTM( &
       xrtm_options_string, &
       xrtm_solvers_string, &
       do_polarization, &
       keep_coef_constant, &
       xrtm_options, &
       xrtm_solvers, &
       xrtm_kernels, &
       success)

    implicit none
    type(string), allocatable, intent(in) :: xrtm_options_string(:)
    type(string), allocatable, intent(in) :: xrtm_solvers_string(:)
    logical, intent(in) :: do_polarization
    logical, intent(in) :: keep_coef_constant
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
       xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_UPWELLING_OUTPUT)
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
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)

       else if (tmp_str == "TWO_STREAM") then

          ! Dedicated 2-Stream solver
          call logger%debug(fname, "Using XRTM in two-stream mode")
          xrtm_solvers(i) = XRTM_SOLVER_TWO_STREAM
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SFI)

       else if (tmp_str == "EIG_BVP") then

          ! Quadrature (LIDORT-like)
          xrtm_solvers(i) = XRTM_SOLVER_EIG_BVP
          xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_SFI)

          !call logger%debug(fname, "Using XRTM in discrete-ordinate mode")
          !call logger%debug(fname, "Using N-T TMS correction and Delta-M scaling.")
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_N_T_TMS)
          !xrtm_options(i) = ior(xrtm_options(i), XRTM_OPTION_DELTA_M)

       else if (tmp_str == "MEM_BVP") then

          ! Quadrature (LIDORT-like)
          call logger%debug(fname, "Using XRTM in matrix exponential discrete-ordinate mode")
          xrtm_solvers(i) = XRTM_SOLVER_MEM_BVP

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
       xrtm_success)

    type(CS_window_t), intent(in) :: CS_win
    type(CS_general_t), intent(in) :: CS_general
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    type(statevector), intent(in) :: SV
    type(scene), intent(in) :: scn
    integer, intent(in) :: n_stokes

    double precision, intent(inout) :: radiance_calc_work_hi_stokes(:,:)
    double precision, allocatable, intent(inout) :: dI_dgas(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dsurf(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dTemp(:,:)
    double precision, allocatable, intent(inout) :: dI_dAOD(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dAHeight(:,:,:)
    double precision, allocatable, intent(inout) :: dI_dpsurf(:,:)

    logical, intent(inout) :: xrtm_success

    ! Local variables

    ! Function name for logger
    character(len=*), parameter :: fname = "solve_RT_problem_XRTM"
    character(len=999) :: tmp_str

    double precision, allocatable :: monochr_radiance(:,:)
    double precision, allocatable :: monochr_weighting_functions(:,:,:)

    type(PCA_handler_t) :: PCA_handler
    type(string), allocatable :: PCA_xrtm_solvers_lo(:), PCA_xrtm_solvers_hi(:)
    type(scene), allocatable :: PCA_scn(:)

    integer :: num_act_lev
    integer :: num_lay
    integer :: num_wl
    integer :: n_mom
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
    integer :: i, j
    !
    double precision :: cpu_start_pca, cpu_stop_pca


    ! Initialize success variable
    xrtm_success = .false.

    ! Grab the number of active levels for easy access
    num_act_lev = scn%num_active_levels
    num_lay = num_act_lev - 1
    num_wl = size(scn%op%wl)

    allocate(n_coef(num_lay))

    if (n_stokes == 1) then
       n_mom = 1
    else if (n_stokes == 3) then
       n_mom = 6
    else
       call logger%fatal(fname, "N_Stokes is neither 1 or 3 - not sure what to do!")
       stop 1
    end if


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

       call setup_XRTM( &
            CS_win%xrtm_options, &
            CS_win%xrtm_solvers, &
            CS_win%do_polarization, &
            CS_win%constant_coef, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_kernels, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to setup_XRTM unsuccessful.")
          return
       end if

       xrtm_n_derivs = ( &
            + SV%num_gas &
            + SV%num_temp &
            + SV%num_albedo &
            + SV%num_aerosol_aod &
            + SV%num_aerosol_height &
            + SV%num_psurf &
            )

       ! Allocate the arrays which hold the radiances and weighting functions
       allocate(monochr_radiance(size(scn%op%wl), scn%num_stokes))
       allocate(monochr_weighting_functions(size(scn%op%wl), scn%num_stokes, xrtm_n_derivs))

       monochr_radiance(:,:) = 0.0d0
       monochr_weighting_functions(:,:,:) = 0.0d0

       ! -------------------------------------------
       !
       ! This outer loop runs over all supplied
       ! XRTM solvers, and the results for jacobians
       ! AND radiances will be ADDED up. The user
       ! must know what they are doing.
       ! 
       ! -------------------------------------------

       do i = 1, size(xrtm_solvers)

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
          ! --------------------------------------

          if ( &
               (iand(xrtm_solvers(i), XRTM_SOLVER_EIG_BVP) /= 0) .or. &
               (iand(xrtm_solvers(i), XRTM_SOLVER_TWO_OS) /= 0) &
               ) then

             ! Is it even allocated?
             if (.not. allocated(CS_win%RT_streams)) then
                call logger%fatal(fname, "You MUST supply a -rt_streams- option!")
                stop 1
             end if

             ! And if allocated, does it have the right size?
             if (size(CS_win%RT_streams) < 1) then
                call logger%fatal(fname, "You need to supply at least one value for -rt_streams-")
                stop 1
             end if

             ! And then finally it needs to have the right value
             if (CS_win%RT_streams(1) < 2) then
                write(tmp_str, '(A, G0.1)') "Need at least 2 streams, but you said: ", &
                     CS_win%RT_streams(1)
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if

             xrtm_streams = CS_win%RT_streams(1) / 2

          else

             call logger%debug(fname, "Setting RT streams to 1 per hemisphere.")
             xrtm_streams = 1

          end if


          ! This is where the XRTM handler is created

          call create_XRTM( &
               xrtm, & ! XRTM handler
               xrtm_options(i), & ! XRTM options bitmask
               xrtm_solvers(i), & ! XRTM solvers combined bitmask
               max(3, scn%max_pfmom), & ! Max coef - either 3 for Rayleigh, or whatever we have for aerosols
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
               monochr_radiance, &
               monochr_weighting_functions, &
               xrtm_success)

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

       call logger%debug(fname, "Mapping weighting functions back to Jacobians.")
       ! Recover the Jacobians from the XRTM container
       ! Store gas subcolumn derivatives
       do j=1, SV%num_gas
          dI_dgas(:,j,:) = monochr_weighting_functions(:,:,j)
       end do

       ! Store the temperature offset Jacobian if needed
       if (SV%num_temp == 1) then
          dI_dTemp(:,:) = monochr_weighting_functions(:,:,SV%num_gas + 1)
       end if

       ! Store surface Jacobian
       do j=1, SV%num_albedo
          dI_dsurf(:,j,:) = monochr_weighting_functions(:,:,SV%num_gas + SV%num_temp + j)
       end do

       ! Store aerosol optical depth Jacobians
       do j=1, SV%num_aerosol_aod
          dI_dAOD(:,j,:) = monochr_weighting_functions(:,:,SV%num_gas + SV%num_temp + SV%num_albedo + j)
          ! If this AOD retrieval is in log-space, we need to multiply by
          ! the (linear-space) AOD itself. df/d(ln(x)) = df/dx * x
          if (CS_win%aerosol_retrieve_aod_log(SV%aerosol_aod_idx_lookup(j))) then
             dI_dAOD(:,j,:) = dI_dAOD(:,j,:) * exp(SV%svsv(SV%idx_aerosol_aod(j)))
          end if
       end do

       ! Store aerosol layer height Jacobians
       do j=1, SV%num_aerosol_height
          dI_dAHeight(:,j,:) = monochr_weighting_functions(:,:, &
               SV%num_gas + SV%num_temp + SV%num_albedo + SV%num_aerosol_aod + j)
          if (CS_win%aerosol_retrieve_height_log(SV%aerosol_height_idx_lookup(j))) then
             dI_dAHeight(:,j,:) = dI_dAHeight(:,j,:) * exp(SV%svsv(SV%idx_aerosol_height(j)))
          end if
       end do

       ! Store surface pressure jacobians
       if (SV%num_psurf == 1) then
          dI_dpsurf(:,:) = monochr_weighting_functions(:,:, &
               SV%num_gas + &
               SV%num_temp + &
               SV%num_albedo + &
               SV%num_aerosol_aod + &
               SV%num_aerosol_height + 1)
       end if

    else if (CS_win%RT_strategy%lower() == "pca") then

       ! --------------------------------------------------------------------
       ! PCA-based fast RT
       ! --------------------------------------------------------------------
       !
       ! Look at Somkuti et al. 2016 for details
       ! --------------------------------------------------------------------

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
       call perform_PCA(PCA_handler)

       ! Create scene object array from PCA data
       ! (mean and perturbed states, which we can then stick into XRTM)
       call create_scenes_from_PCA(PCA_handler, scn, PCA_scn)

       ! --------------------------------------------------------------------
       ! PCA set-up done!
       !
       ! All required values are now available inside of the object
       ! "PCA_handler"
       ! --------------------------------------------------------------------

       cpu_stop_pca = get_cpu_time()

       write(*,*) "PCA cost: ", cpu_stop_pca - cpu_start_pca

       ! For the time being, we want to keep the PCA method as:
       ! Low accuracy: SS + 2S
       ! High accuracy: SS + EIG_BVP
       ! (maybe change later on..)
       ! NOTE
       ! In XRTM, TWO_STREAM computes both single scattering AND diffuse 2S contributions.
       ! This is different in LIDORT, 2STREAM etc.

       allocate(PCA_xrtm_solvers_lo(1))
       allocate(PCA_xrtm_solvers_hi(1))

       PCA_xrtm_solvers_lo(1) = "TWO_STREAM"
       PCA_xrtm_solvers_hi(1) = "SINGLE"


       call setup_XRTM( &
            CS_win%xrtm_options, &
            PCA_xrtm_solvers_lo, &
            .false., & ! Polarization? Not for low accuracy run!
            CS_win%constant_coef, &
            xrtm_options, &
            xrtm_solvers, &
            xrtm_kernels, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to setup_XRTM unsuccessful.")
          return
       end if

       ! Set up number of derivatives needed by XRTM
       xrtm_n_derivs = ( &
            + SV%num_gas &
            + SV%num_temp &
            + SV%num_albedo &
            + SV%num_aerosol_aod &
            + SV%num_aerosol_height &
            + SV%num_psurf &
            )

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
            max(3, scn%max_pfmom), & ! Max coef - either 3 for Rayleigh, or whatever we have for aerosols
            1, & ! Quadrature points per hemisphere
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

       call XRTM_monochromatic( &
            xrtm, &
            xrtm_n_derivs, &
            scn, &
            SV, &
            CS_win, &
            CS_general, &
            CS_aerosol, &
            monochr_radiance, &
            monochr_weighting_functions, &
            xrtm_success)

       if (.not. xrtm_success) then
          call logger%error(fname, "Call to XRTM_monochromatic unsuccessful.")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Need to clean up xrtm object
       call xrtm_destroy_f90(xrtm, xrtm_error)

       ! --------------------------------------------------------------------
       !
       !
       ! Perform high-accuracy binned calculations.
       !
       !
       ! --------------------------------------------------------------------

       


       stop 1

    else

       call logger%fatal(fname, "RT strategy: '" // CS_win%RT_strategy%chars() &
            // "' not recognized. Valid options are..")

       call logger%fatal(fname, "(1) monochromatic")
       !call logger%fatal(fname, "(2) PCA (work in progress)")
       stop 1

    end if


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
       radiance, &
       weighting_functions, &
       success)

    type(xrtm_type), intent(inout) :: xrtm
    integer, intent(in) :: n_derivs
    type(scene), intent(in) :: scn
    type(statevector), intent(in) :: SV
    type(CS_window_t), intent(in) :: CS_win
    type(CS_general_t), intent(in) :: CS_general
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    ! Container for radiance result (specetral, n_stokes)
    double precision, intent(inout) :: radiance(:,:)
    ! Container for weighting function results (spectral, n_stokes, n_derivs)
    double precision, intent(inout) :: weighting_functions(:,:,:)
    logical, intent(inout) :: success


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
    ! Linearized optical depth input (spectral, n_derivs, n_layers)
    double precision, allocatable :: ltau(:,:,:)
    ! Linearized single scatter alebdos (spectral, n_derivs, n_layers)
    double precision, allocatable :: lomega(:,:,:)
    ! Linearized surface inputs (specetral, n_derivs)
    double precision, allocatable :: lsurf(:,:)
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
    allocate(ltau(n_wl, n_derivs, n_layers))
    allocate(lomega(n_wl, n_derivs, n_layers))
    allocate(lsurf(n_wl, n_derivs))

    ltau(:,:,:) = 0.0d0
    lomega(:,:,:) = 0.0d0
    lsurf(:,:) = 0.0d0

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
    deriv_counter = 0

    ! --------------------------
    !
    ! Gas sub-column derivatives
    !
    ! --------------------------

    do j = 1, SV%num_gas
       ltau(:, j, SV%s_start(j):SV%s_stop(j)-1) = &
            scn%op%gas_tau(:, SV%s_start(j):SV%s_stop(j)-1, SV%gas_idx_lookup(j)) / SV%svsv(SV%idx_gas(j,1))

       lomega(:, j, SV%s_start(j):SV%s_stop(j)-1) = &
            -scn%op%layer_omega(:,SV%s_start(j):SV%s_stop(j)-1) &
            / scn%op%layer_tau(:,SV%s_start(j):SV%s_stop(j)-1) &
            * ltau(:, j, SV%s_start(j):SV%s_stop(j)-1)

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A, G0.1)') "Derivative #", deriv_counter, ": gas sub-column #", j
       call logger%debug(fname, trim(tmp_str))
    end do

    ! ----------------------
    !
    ! Temperature derivative
    !
    ! ----------------------

    if (SV%num_temp == 1) then
       j = SV%num_gas + 1

       ltau(:,j,:) = sum(scn%op%gas_tau_dtemp(:,1:n_layers,:), dim=3)
       lomega(:,j,:) = -scn%op%layer_omega(:,1:n_layers) / scn%op%layer_tau(:,1:n_layers) * ltau(:,j,:)

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A)') "Derivative #", deriv_counter, ": temperature"
       call logger%debug(fname, trim(tmp_str))
    end if

    ! --------------------------------------
    !
    ! Surface albedo coefficient derivatives
    !
    ! --------------------------------------

    do j = 1, SV%num_albedo
       i = SV%num_gas + SV%num_temp + j

       lsurf(:,i) = (scn%op%wl(:) - scn%op%wl(int(size(scn%op%wl) / 2)))**(dble(j-1))

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A, G0.1)') "Derivative #", deriv_counter, ": albedo parameter #", j
       call logger%debug(fname, trim(tmp_str))
    end do

    ! -----------------------
    !
    ! Aerosol AOD derivatives
    !
    ! -----------------------

    do j = 1, SV%num_aerosol_aod
       i = SV%num_gas + SV%num_temp + SV%num_albedo + j
       aer_idx = SV%aerosol_aod_idx_lookup(j)

       ltau(:,i,:) = scn%op%aer_ext_tau(:,1:n_layers,aer_idx) / scn%op%reference_aod(aer_idx)
       ! This should be tau_aer_ext / AOD * (omega_a/tau - omega/tau)
       lomega(:,i,:) = ltau(:,i,:) * ( &
            (scn%op%aer_sca_tau(:,1:n_layers,aer_idx) / scn%op%aer_ext_tau(:,1:n_layers,aer_idx)  &
            - scn%op%layer_omega(:,1:n_layers)) / scn%op%layer_tau(:,1:n_layers) &
            )

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A, G0.1)') "Derivative #", deriv_counter, ": aerosol OD #", j
       call logger%debug(fname, trim(tmp_str))
    end do

    ! --------------------------
    !
    ! Aerosol height derivatives
    !
    ! --------------------------
    !
    ! NOTE / TODO
    ! This is also a bit hacky since we take the aerosol width from the MCS, making this
    ! code a bit messy. It also works because right now, we can't change the aerosol
    ! layer width at this point in time. (OCO-like instruments not really sensitive to layer width)

    do j = 1, SV%num_aerosol_height

       allocate(aero_fac(n_layers))

       i = SV%num_gas + SV%num_temp + SV%num_albedo + SV%num_aerosol_aod + j
       aer_idx = SV%aerosol_height_idx_lookup(j)

       ! This factor is required for the aerosol height Jacobians
       call calculate_aero_height_factors( &
            scn%atm%p_layers(:), &
            exp(SV%svsv(SV%idx_aerosol_height(j))) * scn%atm%psurf, &
            CS_aerosol(scn%op%aer_mcs_map(aer_idx))%default_width, &
            aero_fac)

       do l = 1, n_layers
          ltau(:,i,l) = scn%op%aer_ext_tau(:,l,aer_idx) * aero_fac(l) * scn%atm%psurf
       end do

       lomega(:,i,:) = ltau(:,i,:) * ( &
            (scn%op%aer_sca_tau(:,:,aer_idx) / scn%op%aer_ext_tau(:,:,aer_idx)  &
            - scn%op%layer_omega(:,1:n_layers)) / scn%op%layer_tau(:,1:n_layers) &
            )

       deallocate(aero_fac)

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A, G0.1)') "Derivative #", deriv_counter, ": aerosol height #", j
       call logger%debug(fname, trim(tmp_str))
    end do

    ! ---------------------------
    !
    ! Surface pressure derivative
    !
    ! ---------------------------

    if (SV%num_psurf == 1) then
       i = SV%num_gas + SV%num_temp + SV%num_albedo + SV%num_aerosol_aod + SV%num_aerosol_height + 1

       ! Make sure to keep the minus sign here!
       do l = 1, n_layers
          ! Sum over contributions from all gases, since we have dTau/dPsurf per gas
          ltau(:, i, l) = -sum(scn%op%gas_tau_dpsurf(:, l, :), dim=2)
          lomega(:, i, l) = -scn%op%layer_omega(:, l) / scn%op%layer_tau(:, l) * ltau(:, i, l)
       end do

       deriv_counter = deriv_counter + 1
       write(tmp_str, '(A, G0.1, A)') "Derivative #", deriv_counter, ": surface pressure"
       call logger%debug(fname, trim(tmp_str))
    end if

    ! Calculate coef and lcoef at band edges
    call logger%debug(fname, "Computing phase function coefficients and derivatives at band edges.")
    call compute_coef_at_wl(scn, CS_aerosol, SV, scn%op%wl(1), n_mom, n_derivs, &
         scn%op%ray_tau(1,:), scn%op%aer_sca_tau(1,:,:), &
         coef_left, lcoef_left, coef_success)

    if (.not. coef_success) then
       return
    end if

    call compute_coef_at_wl(scn, CS_aerosol, SV, scn%op%wl(n_wl), n_mom, n_derivs, &
         scn%op%ray_tau(n_wl,:), scn%op%aer_sca_tau(n_wl,:,:), &
         coef_right, lcoef_right, coef_success)

    if (.not. coef_success) then
       return
    end if

    call logger%debug(fname, "Determining number of phase function coefficients per layer used in XRTM.")
    ! Figure out how many coefficients we are giving to XRTM
    ! The default is 3 for Rayleigh-only
    allocate(n_coef(n_layers))
    n_coef(:) = 3
    if (scn%num_aerosols > 0) then
       ! Loop through all layers and look at how much aerosol
       ! extinction is present if below some threshold value,
       ! we just assume that 3 coefficients is enough for a
       ! (at best) Rayleigh-only layer
       do l = 1, n_layers
          if (sum(scn%op%aer_ext_tau(1,l,:)) < 1.0d-6) then
             n_coef(l) = 3
          else
             n_coef(l) = scn%max_pfmom
          end if
       end do
    end if

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

       call xrtm_set_sos_params_f90(xrtm, 4, 10.0d0, 0.01d0, xrtm_error)
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

       single_ltau(:,:) = ltau(i,:,:)
       single_lomega(:,:) = lomega(i,:,:)
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
       if (mod(i, N_spec/10) == 0) then
          write(tmp_str, '(A, F6.2, A, G0.1, A, G0.1, A)') &
               "XRTM calls (", 100.0 * i/N_spec, "%, ", i, "/", N_spec, ")"
          call logger%debug(fname, trim(tmp_str))
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

end module XRTM_mod
