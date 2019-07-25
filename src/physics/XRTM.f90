!> @brief Module to house the interfaces and calculations for the XRTM RT model
!> @file XRTM.f90
!> @author Peter Somkuti


module XRTM_mod

  ! User modules
  use math_utils_mod

  ! Third-party modules
  use xrtm_int_f90
  use stringifor
  use logger_mod, only: logger => master_logger

  implicit none

  integer :: xrtm_options
  integer :: xrtm_solvers
  integer :: xrtm_kernels(1)

  public :: setup_XRTM, calculate_XRTM_radiance

contains



  subroutine setup_XRTM(xrtm, xrtm_options_string, xrtm_solvers_string, success)
    implicit none
    type(xrtm_type), intent(inout) :: xrtm
    type(string), intent(in) :: xrtm_options_string(:)
    type(string), intent(in) :: xrtm_solvers_string(:)
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "setup_XRTM"

    type(string) :: tmp_str
    integer :: N, i, j
    integer :: xrtm_error

    integer :: out_levels(1)


    success = .false.

    ! Initialize XRTM options and solveres with 0
    xrtm_options = 0
    xrtm_solvers = 0
    ! We always want TOA radiance
    out_levels(1) = 0
    ! Surface kernels
    xrtm_kernels(1) = XRTM_KERNEL_LAMBERTIAN


    ! Populate solvers bit field, depending on requested options
    do i=1, size(xrtm_solvers_string)
       ! Grab local copy of string
       tmp_str = xrtm_solvers_string(i)

       if (tmp_str == "TWO_OS") then
          ! Vijay-like 2OS, gives you SECOND ORDER ONLY
          xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_TWO_OS)
       else if (tmp_str == "SINGLE") then
          ! Single scattering only
          xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_SINGLE)
       else if (tmp_str == "EIG_BVP") then
          ! Quadrature (LIDORT-like)
          xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_EIG_BVP)
       else
          call logger%error(fname, "XRTM solver option is not implemented, and will be ignored: " &
               // tmp_str%chars())
       end if
    end do

    xrtm_options = XRTM_OPTION_CALC_DERIVS
    !xrtm_solvers = ior(XRTM_SOLVER_SINGLE, XRTM_SOLVER_TWO_OS)
    !xrtm_solvers = ior(xrtm_solvers, XRTM_SOLVER_EIG_BVP)
    xrtm_solvers = XRTM_SOLVER_EIG_ADD
    write(*,*) "Solvers: ", xrtm_solvers
    ! Some general set-up's that are not expected to change, and can
    ! be safely hard-coded.
    call xrtm_set_fourier_tol_f90(xrtm, .0001d0, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_fourier_tol_f90()")
       return
    endif

    success = .true.

  end subroutine setup_XRTM


  subroutine create_XRTM(xrtm, &
       max_coef, n_quad, n_stokes, n_derivs, n_layers, n_kernels, &
       n_kernel_quad, success)
    implicit none

    type(xrtm_type), intent(inout) :: xrtm
    integer, intent(in) :: max_coef
    integer, intent(in) :: n_quad
    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers
    integer, intent(in) :: n_kernels
    integer, intent(in) :: n_kernel_quad
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
       call logger%error(fname, "Error calling xrtm_create_f90)")
       return
    endif

    success = .true.

  end subroutine create_XRTM



  subroutine calculate_XRTM_radiance(xrtm, wavelengths, SZA, VZA, SAA, VAA, &
       albedo, tau, &
       n_stokes, n_derivs, n_layers, radiance)

    implicit none
    type(xrtm_type), intent(inout) :: xrtm
    double precision, intent(in) :: wavelengths(:)
    double precision, intent(in) :: SZA, VZA, SAA, VAA, albedo(:)
    double precision, intent(in) :: tau(:,:)
    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    integer, intent(in) :: n_layers
    double precision, intent(inout) :: radiance(:)

    double precision, allocatable    :: I_p(:,:,:,:)
    double precision, allocatable    :: I_m(:,:,:,:)
    double precision, allocatable    :: K_p(:,:,:,:,:)
    double precision, allocatable    :: K_m(:,:,:,:,:)

    character(len=*), parameter :: fname = "calculate_XRTM_radiance"
    integer :: xrtm_error
    double precision :: out_thetas(1)
    double precision :: out_phis(1,1)
    double precision, allocatable :: coef(:,:,:)
    integer :: out_levels(1)
    integer :: n_coef(n_layers)
    double precision :: ltau(n_layers)
    integer :: N_spec
    integer :: i

    integer :: funit

    double precision :: cpu_start, cpu_end

    double precision, allocatable :: omega_test(:)

    N_spec = size(wavelengths)
    out_thetas(1) = VZA
    out_phis(1,1) = VAA
    out_levels(1) = 0
    allocate(omega_test(size(tau, 2)))
    omega_test(:) = 1.0d-4
    omega_test(10) = 0.25d0
    ltau(:) = 1.0d0

    allocate(coef(7, 1, n_layers))
    n_coef(:) = 7

    coef(1, :, :) = 1.0d0
    coef(2:7, :, :) = 0.0d0

    coef(1, 1, 10) = 1.000000e+00
    coef(2, 1, 10) = 1.865569e+00
    coef(3, 1, 10) = 1.789985e+00
    coef(4, 1, 10) = 1.220838e+00
    coef(5, 1, 10) = 7.472409e-01
    coef(6, 1, 10) = 4.017337e-01
    coef(7, 1, 10) = 2.173326e-01


    ! Allocate output containers - these will not change with wavelength, hence
    ! we only need to do this once per microwindow
    allocate(I_p(n_stokes, 1, 1, 1))
    allocate(I_m(n_stokes, 1, 1, 1))
    allocate(K_p(n_stokes, 1, 1, n_derivs, 1))
    allocate(K_m(n_stokes, 1, 1, n_derivs, 1))

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! First, plug in the values that XRTM needs for its calculations, which DO NOT
    ! depend on wavelength (i.e. do not need to be called repeatedly), such as
    ! viewing geometry.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    endif

    ! Plug in the solar zenith angle
    call xrtm_set_theta_0_f90(xrtm, SZA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_theta_0_f90")
       return
    endif

    ! Plug in the viewing zenith angle
    call xrtm_set_out_thetas_f90(xrtm, out_thetas, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_out_thetas_f90")
       return
    endif

    ! Plug in the solar azimuth angle
    call xrtm_set_phi_0_f90(xrtm, SAA, xrtm_error)
    if (xrtm_error /= 0) then
       call logger%error(fname, "Error calling xrtm_set_phi_0_f90")
       return
    endif

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Then, loop through the wavelengths and stick in wavelength-dependent quantities,
    ! i.e. coefficients, optical depths and single-scattering albedos etc., and call
    ! the XRTM solver to get radiances - this is obviously a performance-critical section
    ! and is easily the most costly portion of the entire forward model.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) albedo(1)

    open(file="xrtm.debug", newunit=funit)
    write(*,*) "Starting monochromatic loop"
    call cpu_time(cpu_start)
    do i=1, N_spec

       ! Plug in surface property
       call xrtm_set_kernel_ampfac_f90(xrtm, 0, albedo(i), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_kernel_ampfac_f90")
          return
       endif

       ! Plug in optical depth
       call xrtm_set_ltau_n_f90(xrtm, tau(i,:), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_n_f90")
          return
       endif

       call xrtm_set_ltau_l_1n_f90(xrtm, 1, ltau(:), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_ltau_1n_f90")
          return
       endif

       ! Plug in the single-scatter albedo
       call xrtm_set_omega_n_f90(xrtm, omega_test(:), xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_omega_n_f90")
          return
       endif

       ! Plug in the layer scattering matrix
       call xrtm_set_coef_n_f90(xrtm, n_coef, coef, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_set_coef_n_f90")
          return
       endif

       ! Call this function for whatever reason
       call xrtm_update_varied_layers_f90(xrtm, xrtm_error)
       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_update_varied_layers_f90")
          return
       endif


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

       write(funit, *) I_p(1,1,1,1), K_p(1,1,1,1,1)

       if (xrtm_error /= 0) then
          call logger%error(fname, "Error calling xrtm_radiance_f90)")
          return
       endif

    end do

    close(funit)

    write(*,*) "Finished monochromatic loop"
    call cpu_time(cpu_end)
    writE(*,*) cpu_end - cpu_start, "seconds"


  end subroutine calculate_XRTM_radiance

end module XRTM_mod
