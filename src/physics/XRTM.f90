!> @brief Module to house the interfaces and calculations for the XRTM RT model
!> @file XRTM.f90
!> @author Peter Somkuti

module XRTM_mod

  ! User modules
  use math_utils_mod
  use Rayleigh_mod, only : calculate_rayleigh_scatt_matrix

  ! Third-party modules
  use xrtm_int_f90
  use stringifor
  use logger_mod, only: logger => master_logger

  implicit none

  public :: setup_XRTM, create_XRTM, calculate_XRTM_radiance

contains


  !> @brief "Translates" configuration values into XRTM settings
  !> @param xrtm_options_string Array of strings containing XRTM options
  !> @param xrtm_solvers_string Array of strings containing XRTM solvers
  !> @param xrtm_options Integer bitmask for XRTM options
  !> @param xrtm_solvers Integer bitmask for XRTM solvers
  !> @param xrtm_kernels Integer array specifying which surface kernels
  !> @param success Did any error occur?
  subroutine setup_XRTM(xrtm_options_string, xrtm_solvers_string, &
       xrtm_options, xrtm_solvers, xrtm_kernels, success)

    implicit none
    type(string), allocatable, intent(in) :: xrtm_options_string(:)
    type(string), allocatable, intent(in) :: xrtm_solvers_string(:)
    integer, intent(inout) :: xrtm_options
    integer, intent(inout) :: xrtm_solvers
    integer, intent(inout) :: xrtm_kernels(:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "setup_XRTM"

    type(string) :: tmp_str
    integer :: N, i, j
    integer :: xrtm_error

    success = .false.

    ! Initialize XRTM options and solvers with 0
    xrtm_options = 0
    xrtm_solvers = 0

    ! Surface kernels
    xrtm_kernels(1) = XRTM_KERNEL_LAMBERTIAN


    ! Populate solvers bit field, depending on requested options
    if (allocated(xrtm_solvers_string)) then
       do i=1, size(xrtm_solvers_string)
          ! Grab local copy of string
          tmp_str = xrtm_solvers_string(i)

          if (tmp_str == "TWO_OS") then
             ! Vijay-like 2OS, gives you SECOND ORDER ONLY
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_TWO_OS)
             call logger%debug(fname, "Using XRTM in 2-orders-of-scattering mode")
          else if (tmp_str == "SINGLE") then
             ! Single scattering only
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_SINGLE)
             call logger%debug(fname, "Using XRTM in single-scattering mode")
          else if (tmp_str == "EIG_BVP") then
             ! Quadrature (LIDORT-like)
             call logger%debug(fname, "Using XRTM in discrete-ordinate mode")
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_EIG_BVP)
          else if (tmp_str == "PADE_ADD") then
             ! Pade approximation and adding
             call logger%debug(fname, "Using XRTM in Pade-adding mode")
             xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_PADE_ADD)
          else
             call logger%error(fname, "XRTM solver option is not implemented, and will be ignored: " &
                  // tmp_str%chars())
          end if
       end do
    end if


    xrtm_options = ior(XRTM_OPTION_CALC_DERIVS, XRTM_OPTION_PSA)
    !xrtm_options = ior(xrtm_options, XRTM_OPTION_N_T_TMS)
    !xrtm_options = ior(xrtm_options, XRTM_OPTION_DELTA_M)
    !xrtm_options = XRTM_OPTION_CALC_DERIVS
    !xrtm_solvers = ior(XRTM_SOLVER_SINGLE, XRTM_SOLVER_TWO_OS)
    !xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_EIG_BVP)
    !xrtm_solvers = XRTM_SOLVER_SINGLE
    !write(*,*) "Solvers: ", xrtm_solvers

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



  subroutine calculate_XRTM_radiance(xrtm, wavelengths, SZA, VZA, SAA, VAA, &
       altitude_levels, albedo, gas_tau, ray_tau, ray_depolf, gas_scale, &
       n_stokes, n_derivs, n_layers, &
       s_start, s_stop, gas_lookup, &
       radiance, dI_dgas, dI_dsurf, &
       success)

    implicit none
    type(xrtm_type), intent(inout) :: xrtm
    double precision, intent(in) :: wavelengths(:)
    double precision, intent(in) :: SZA, VZA, SAA, VAA, albedo(:)
    double precision, intent(in) :: altitude_levels(:)
    double precision, intent(in) :: gas_tau(:,:,:)
    double precision, intent(in) :: ray_tau(:,:)
    double precision, intent(in) :: ray_depolf(:)
    double precision, allocatable, intent(in) :: gas_scale(:)
    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers
    integer, intent(in) :: s_start(:)
    integer, intent(in) :: s_stop(:)
    integer, allocatable, intent(in) :: gas_lookup(:)

    double precision, intent(inout) :: radiance(:)
    double precision, intent(inout) :: dI_dgas(:,:)
    double precision, intent(inout) :: dI_dsurf(:)

    logical, intent(inout) :: success


    ! Local

    integer :: xrtm_options
    integer :: xrtm_solvers

    double precision, allocatable    :: I_p(:,:,:,:)
    double precision, allocatable    :: I_m(:,:,:,:)
    double precision, allocatable    :: K_p(:,:,:,:,:)
    double precision, allocatable    :: K_m(:,:,:,:,:)

    character(len=*), parameter :: fname = "calculate_XRTM_radiance"
    character(len=999) :: tmp_str
    integer :: xrtm_error
    double precision :: out_thetas(1)
    double precision :: out_phis(1,1)
    double precision, allocatable :: ray_coef(:,:)
    double precision, allocatable :: coef(:,:,:)
    integer :: out_levels(1)
    integer :: n_coef(n_layers)
    double precision :: ltau(n_derivs, n_layers), lomega(n_derivs, n_layers), lsurf(n_derivs)
    double precision :: lcoef(3, 1, n_derivs, n_layers)
    double precision :: tau(n_layers), omega(n_layers), omega_tmp(n_layers)
    integer :: N_spec
    integer :: i,j,k,l,m

    integer :: funit

    double precision :: cpu_start, cpu_end

    success = .false.

    xrtm_options = xrtm_get_options_f90(xrtm)
    xrtm_solvers = xrtm_get_solvers_f90(xrtm)

    N_spec = size(wavelengths)

    radiance(:) = 0.0d0
    dI_dgas(:,:) = 0.0d0
    dI_dsurf(:) = 0.0d0

    out_thetas(1) = VZA
    out_phis(1,1) = VAA
    out_levels(1) = 0

    ltau(:,:) = 0.0d0
    lomega(:,:) = 0.0d0
    lcoef(:,:,:,:) = 0.0d0
    lsurf(:) = 0.0d0

    ! Third derivative dI/dsurf
    lsurf(n_derivs) = 1.0d0

    allocate(ray_coef(3, 1))
    allocate(coef(3, 1, n_layers))
    n_coef(:) = 3

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
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    endif

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, plug in the values that XRTM needs for its calculations, which DO NOT
    ! depend on wavelength (i.e. do not need to be called repeatedly), such as
    ! viewing geometry.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (iand(xrtm_options, XRTM_OPTION_PSA) /= 0) then
       ! The following two calls are only valid if XRTM has been created
       ! with the flag to use pseudo-spherical approximation for a curved
       ! atmosphere.

       ! Set Earth's radius
       call xrtm_set_planet_r_f90(xrtm, EARTH_EQUATORIAL_RADIUS, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_planet_r_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Set model atmosphere altitude levels
       call xrtm_set_levels_z_f90(xrtm, altitude_levels(1:n_layers+1), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_levels_z_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

    end if

    if (iand(xrtm_solvers, XRTM_SOLVER_SOS) /= 0) then
       call xrtm_set_sos_params_f90(xrtm, 2, 10.0d0, 0.01d0, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_sos_params_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       endif
    end if

    if (iand(xrtm_solvers, XRTM_SOLVER_PADE_ADD) /= 0) then
       call xrtm_set_pade_params_f90(xrtm, 0, 2, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_pade_params_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       endif
    end if

    ! Set the output level to TOA
    call xrtm_set_out_levels_f90(xrtm, out_levels, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_out_levels_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    endif

    ! We want to calculate sun-normalized values, so set F_0 to 1.0
    call xrtm_set_F_0_f90(xrtm, 1.0d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_F_0_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    end if

    ! Plug in the solar zenith angle
    call xrtm_set_theta_0_f90(xrtm, SZA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_theta_0_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    end if

    ! Plug in the viewing zenith angle
    call xrtm_set_out_thetas_f90(xrtm, out_thetas, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_out_thetas_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    end if

    ! Plug in the solar azimuth angle
    call xrtm_set_phi_0_f90(xrtm, SAA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_phi_0_f90")
       call xrtm_destroy_f90(xrtm, xrtm_error)
       return
    end if

    ! -----------------------------------------------------------------------------------
    ! Then, loop through the wavelengths and stick in wavelength-dependent quantities,
    ! i.e. coefficients, optical depths and single-scattering albedos etc., and call
    ! the XRTM solver to get radiances - this is obviously a performance-critical section
    ! and is easily the most costly portion of the entire forward model.
    ! -----------------------------------------------------------------------------------

    !open(file="xrtm.debug", newunit=funit)
    !write(*,*) "Starting monochromatic loop"
    call cpu_time(cpu_start)
    do i=1, N_spec

       ! if (mod(i, 100) == 0) write(*,*) i, N_spec

       ! TOTAL atmospheric optical properties - these go into the
       ! RT code for the forward calculation.

       ! Calculate total tau
       tau(:) = sum(gas_tau(i,:,:), dim=2) + ray_tau(i, :)
       ! Calculate the total single scatter albedo
       omega(:) = ray_tau(i, :) / tau(:)

       ! Linearized inputs:
       ! Gas sub-columns - NOTE that the l_tau and l_omega are WITH
       ! RESPECT TO THE gas species, and NOT the total tau / omega
       ltau(:,:) = 0.0d0
       lomega(:,:) = 0.0d0

       do j=1, size(s_start)
          ltau(j, s_start(j):s_stop(j)-1) = gas_tau(i, s_start(j):s_stop(j)-1, gas_lookup(j)) / gas_scale(j)
          lomega(j, s_start(j):s_stop(j)-1) = -omega(s_start(j):s_stop(j)-1) / tau(s_start(j):s_stop(j)-1) * ltau(j, s_start(j):s_stop(j)-1)
       end do

       ! Calculate the rayleigh scattering matrix
       call calculate_rayleigh_scatt_matrix(ray_depolf(i), ray_coef)

       do j=1, n_layers
          coef(:,:,j) = ray_coef(:,:)
       end do

       ! Plug in surface property
       call xrtm_set_kernel_ampfac_f90(xrtm, 0, albedo(i), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampfac_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in optical depth
       call xrtm_set_ltau_n_f90(xrtm, tau, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_n_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the layer optical depth derivatives
       call xrtm_set_ltau_l_nn_f90(xrtm, ltau, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_nn_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the single-scatter albedo
       call xrtm_set_omega_n_f90(xrtm, omega, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_n_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the single-scatter albedo derivatives
       call xrtm_set_omega_l_nn_f90(xrtm, lomega, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_l_nn_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the layer scattering matrix
       call xrtm_set_coef_n_f90(xrtm, n_coef, coef, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_coef_n_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the layer scattering matrix derivatives
       call xrtm_set_coef_l_nn_f90(xrtm, lcoef, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_coef_l_nn_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Plug in the surface derivatives
       call xrtm_set_kernel_ampfac_l_n_f90(xrtm, 0, lsurf(:), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampflac_l_n_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Call this function for whatever reason
       call xrtm_update_varied_layers_f90(xrtm, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_update_varied_layers_f90")
          call xrtm_destroy_f90(xrtm, xrtm_error)
          return
       end if

       ! Calculate TOA radiance!
       call xrtm_radiance_f90(xrtm, &
            xrtm_solvers, &
            1 , & ! Number of output azimuths
            out_phis, & ! Output azimuths
            I_p, & ! Upward radiances,
            I_m, & ! Downward radiances
            K_p, & ! Upward jacobians
            K_m, & ! Downward jacobians
            xrtm_error)

       if (xrtm_error /= 0) then
           call logger%error(fname, "Error calling xrtm_radiance_f90")
           call xrtm_destroy_f90(xrtm, xrtm_error)
           return
       end if

       radiance(i) = I_p(1,1,1,1)

       ! Store gas subcolumn derivatives
       do j=1, size(s_start)
          dI_dgas(i,j) = K_p(1,1,1,j,1)
       end do
       ! Store surface jacobian
       dI_dsurf(i) = K_p(1,1,1,size(s_start)+1,1)

       !write(funit, *) I_p(1,1,1,1), (K_p(1,1,1,j,1), j=1, n_derivs)

    end do

    !close(funit)

    !write(*,*) "Finished monochromatic loop"
    call cpu_time(cpu_end)
    write(tmp_str, '(A, F7.3, A)') "XRTM monochromatic calculations: ", cpu_end - cpu_start, " sec"
    call logger%debug(fname, trim(tmp_str))

    success = .true.

  end subroutine calculate_XRTM_radiance

end module XRTM_mod
