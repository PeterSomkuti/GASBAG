module Rayleigh_mod

  use math_utils_mod, only : NA
  public calculate_rayleigh_tau


contains

  !> @brief Calculates extinction and depolarization due to Rayleigh scattering
  !> @param wl Wavelenghts
  !> @param p Pressure levels
  !> @param rayleigh_tau Rayleigh optical depth
  !> @param rayleigh_depolf Rayleigh depolarization factor
  !> @param co2 Optional CO2 VMR profile
  subroutine calculate_rayleigh_tau(wl, p, rayleigh_tau, rayleigh_depolf, co2)

    double precision, intent(in) :: wl(:)
    double precision, intent(in) :: p(:)
    double precision, optional, intent(in) :: co2(:)
    double precision, intent(inout) :: rayleigh_tau(:,:)
    double precision, intent(inout) :: rayleigh_depolf(:)

    double precision, parameter :: N2_C = 78.084d0
    double precision, parameter :: O2_C = 20.946d0
    double precision, parameter :: Ar_C = 1.0d0
    double precision, parameter :: CO2_C = 4d-4

    integer :: i,j,k
    double precision :: F_N2, F_O2
    double precision, parameter :: F_Ar = 0.934d0
    double precision, parameter :: F_CO2 = 1.15d0
    double precision :: depolf
    double precision :: wl2, wl4

    double precision :: n300_1, nCO2_1, ns
    double precision :: ray_sigma

    do i=1, size(wl)

       wl2 = wl(i) * wl(i)
       wl4 = wl2 * wl2

       F_N2 = 1.034d0 + 3.17d0 * 1.0d-4 / wl2
       F_O2 = 1.096d0 + 1.385d0 * 1.0d-3 / wl2 + 1.448 * 1.0d-4 / wl4

       depolf = N2_C * F_N2 + O2_C * F_O2 + Ar_C * F_ar + CO2_C * F_CO2
       depolf = depolf / (N2_C + O2_C + Ar_C + CO2_C)

       rayleigh_depolf(i) = depolf

       n300_1 = 8060.51d0 + 2480990.0d0 / (132.274d0 - (wl(i)**(-2))) + &
            17455.7d0 / (39.32957d0 - (wl(i)**(-2)))
       n300_1 = n300_1 * 1.0d-8
       ns = 1.0d0 + n300_1

       do j=1, size(p) - 1

          ! If we have the CO2 VMRs, we can adjust the value for the
          ! refractive index that is calculated at a reference of 300ppm
          if (present(co2)) then
             nCO2_1 = n300_1 * (1.0d0 + 0.54d0 * ((co2(j) + co2(j+1)) * 0.5d0 - 3.0d-4))
             ns = 1.0d0 + nCO2_1
          end if

          ray_sigma = (ns**2 - 1)**2 / (ns**2 + 2)**2 * 1.1471954d-24 * depolf / wl(i)**4
          rayleigh_tau(i,j) = NA * ray_sigma * (p(j+1) - p(j)) / (9.81d0)

       end do
    end do

  end subroutine calculate_rayleigh_tau


  !> @brief Calculates the Rayleigh scattering matrix for a given depol factor
  !> @param depolf Depolarization factor (see calculate_rayleigh_tau)
  !> @param Phase matrix (max_coef, n_elem)
  subroutine calculate_rayleigh_scatt_matrix(depolf, coeffs)
    double precision, intent(in) :: depolf
    double precision, intent(inout) :: coeffs(:,:)
    double precision :: depol_rho

    depol_rho = (6.0d0 * depolf - 6.0d0) / (7.0d0 * depolf + 3.0d0)

    ! Reset the matrix
    coeffs(:,:) = 0.0d0

    ! Set beta (for scalar transport)
    coeffs(1, 1) = 1.0d0 ! beta_0
    coeffs(3, 1) = (1.0d0 - depol_rho) / (2.0d0 + depol_rho) ! beta_2

    ! Set other coefficients for vector transport
    if (size(coeffs, 2) == 6) then
       coeffs(3, 2) = 6.0d0 * (1.0d0 - depol_rho) / (2.0d0 + depol_rho) ! alpha_2
       coeffs(2, 4) = 3.0d0 * (1.0d0 - 2.0d0 * depol_rho) / (2.0d0 + depol_rho) ! delta_1
       coeffs(3, 5) = 2.449489742d0 * (1.0d0 - depol_rho) / (2.0d0 + depol_rho)
       ! gamma_2 (that funky number is ~sqrt(6.0)
    end if

  end subroutine calculate_rayleigh_scatt_matrix




end module Rayleigh_mod
