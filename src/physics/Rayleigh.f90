module Rayleigh_mod

  use math_utils_mod, only : NA, pi, DRY_AIR_MASS
  public calculate_rayleigh_tau, calculate_rayleigh_depolf, calculate_rayleigh_scatt_matrix


contains

  !> @brief King factor according to Tomasi et al. 2005
  !> @param wl Wavelength in microns
  !> @param King factor
  pure function King_factor(wl) result(fac)

    double precision, intent(in) :: wl

    double precision, parameter :: N2_C = 78.084d0
    double precision, parameter :: O2_C = 20.946d0
    double precision, parameter :: Ar_C = 1.0d0
    double precision, parameter :: CO2_C = 4d-4

    double precision :: F_N2, F_O2
    double precision, parameter :: F_Ar = 0.934d0
    double precision, parameter :: F_CO2 = 1.15d0

    double precision :: fac
    double precision :: wl2, wl4

    wl2 = wl * wl
    wl4 = wl2 * wl2

    F_N2 = 1.034d0 + 3.17d0 * 1.0d-4 / wl2
    F_O2 = 1.096d0 + 1.385d0 * 1.0d-3 / wl2 + 1.448 * 1.0d-4 / wl4

    fac = N2_C * F_N2 + O2_C * F_O2 + Ar_C * F_ar + CO2_C * F_CO2
    fac = fac / (N2_C + O2_C + Ar_C + CO2_C)

  end function King_factor


  !> @brief Depolarization ratio
  !> @param wl Wavelength in microns
  !> @param depolf Depolarization ratio
  pure function calculate_rayleigh_depolf(wl) result(depol_r)

    double precision, intent(in) :: wl

    double precision :: depol_r
    double precision :: F ! King factor

    F = King_factor(wl)
    depol_r = (6.0d0 * F - 6.0d0) / (7.0d0 * F + 3.0d0)

  end function calculate_rayleigh_depolf



  !> @brief Calculates extinction and depolarization due to Rayleigh scattering
  !> @param wl Wavelengths
  !> @param p Pressure levels
  !> @param rayleigh_tau Rayleigh optical depth
  !> @param rayleigh_depolf Rayleigh depolarization factor
  !> @param co2 Optional CO2 VMR profile
  subroutine calculate_rayleigh_tau(wl, p, g, T, rayleigh_tau, co2)

    double precision, intent(in) :: wl(:)
    double precision, intent(in) :: p(:)
    double precision, intent(in) :: g(:)
    double precision, intent(in) :: T(:)
    double precision, optional, intent(in) :: co2(:)
    
    double precision, intent(inout) :: rayleigh_tau(:,:)

    double precision, parameter :: Nair_ref = 2.546899d19
    double precision, parameter :: p_ref = 101325d0
    double precision, parameter :: T_ref= 288.15d0

    double precision :: F
    double precision :: depolf
    double precision :: Nair
    double precision :: wl2, wl4
    double precision :: wli2
    double precision, parameter :: pi3 = pi * pi * pi
    integer :: i, j, N

    double precision :: n300_1, nCO2_1, ns
    double precision :: ray_sigma


    N = size(p)
    Nair = Nair_ref * (p(N) / p_ref) * (T_ref / T(N))

    do i=1, size(wl)

       ! These are all in microns
       wl2 = wl(i) * wl(i)
       wli2 = 1.0d0 / wl2
       wl4 = wl2 * wl2

       F = King_factor(wl(i))
       depolf = calculate_rayleigh_depolf(wl(i))

       n300_1 = 8060.51d0 + 2480990.0d0 / (132.274d0 - wli2) + &
            17455.7d0 / (39.32957d0 - wli2)
       n300_1 = n300_1 * 1.0d-8
       ns = 1.0d0 + n300_1

       do j=1, size(p) - 1

          ! If we have the CO2 VMRs, we can adjust the value for the
          ! refractive index that is calculated at a reference of 300ppm
          if (present(co2)) then
             nCO2_1 = n300_1 * (1.0d0 + 0.54d0 * ((co2(j) + co2(j+1)) * 0.5d0 - 3.0d-4))
             ns = 1.0d0 + nCO2_1
          end if

          ! Equation 4 of Tomasi et al. 2005
          ! Nair is in [cm^-3]
          ! To convert wl to [cm^-3], 1um = 10^4 cm, so wl in [um^4] = 10^-16 [cm^4]
          ! I'm not sure why Nair is only the surface-level Nair..
          ray_sigma = ( &
               (24.0d0 * pi3 * (ns*ns - 1)**2) / &
               (wl4 * 1.0d-16 * Nair * Nair * (ns*ns + 2)**2) * &
               F &
               )

          rayleigh_tau(i,j) = 0.1d0 * ray_sigma * NA * (p(j+1) - p(j)) / &
               (0.5d0 * (g(j+1) + g(j)) * 1.0d3 * DRY_AIR_MASS)

       end do

    end do

  end subroutine calculate_rayleigh_tau


  !> @brief Calculates the Rayleigh scattering matrix for a given depol factor
  !> @param depolf Depolarization factor (see calculate_rayleigh_tau)
  !> @param Phase matrix (max_coef, n_elem)
  subroutine calculate_rayleigh_scatt_matrix(depolf, coeffs)
    double precision, intent(in) :: depolf
    double precision, intent(inout) :: coeffs(:,:)

    ! Reset the matrix
    coeffs(:,:) = 0.0d0

    ! Set beta (for scalar transport)
    coeffs(1, 1) = 1.0d0 ! beta_0
    coeffs(3, 1) = (1.0d0 - depolf) / (2.0d0 + depolf) ! beta_2

    ! Set other coefficients for vector transport
    if (size(coeffs, 2) == 6) then
       coeffs(3, 2) = 6.0d0 * (1.0d0 - depolf) / (2.0d0 + depolf) ! alpha_2
       coeffs(2, 4) = 3.0d0 * (1.0d0 - 2.0d0 * depolf) / (2.0d0 + depolf) ! delta_1
       coeffs(3, 5) = 2.449489742d0 * (1.0d0 - depolf) / (2.0d0 + depolf) ! -gamma_2
       ! gamma_2 (that funky number is ~sqrt(6.0)
    end if

  end subroutine calculate_rayleigh_scatt_matrix




end module Rayleigh_mod
