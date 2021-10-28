module Rayleigh_mod

  use math_utils_mod, only : NA, pi, DRY_AIR_MASS, EPSILON
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



  !> @brief Calculates extinction and depolarization due to Rayleigh scattering
  !> @param wl Wavelengths
  !> @param p Pressure levels
  !> @param rayleigh_tau Rayleigh optical depth
  !> @param rayleigh_depolf Rayleigh depolarization factor
  !> @param co2 Optional CO2 VMR profile
  subroutine calculate_rayleigh_tau_MS3(wl, p, g, sh, T, rayleigh_tau, co2)

    double precision, intent(in) :: wl(:)
    double precision, intent(in) :: p(:)
    double precision, intent(in) :: g(:)
    double precision, intent(in) :: sh(:)
    double precision, intent(in) :: T(:)
    double precision, optional, intent(in) :: co2(:)

    double precision, intent(inout) :: rayleigh_tau(:,:)

    double precision :: F

    double precision :: Nair
    double precision :: wl2, wl4
    double precision :: wli2
    double precision, parameter :: pi3 = pi * pi * pi
    integer :: i, j, N

    double precision :: n300_1, nCO2_1, ns
    double precision :: depolf
    double precision, allocatable :: ray_sigma(:)


    allocate(ray_sigma(size(wl)))
    N = size(p)

    depolf = ray_depol_MS3(wl(1))
    do j=1, size(p) - 1

       ray_sigma(:) = ray_tau_MS3(wl(:), depolf)
       rayleigh_tau(:,j) = ray_sigma(:) * (p(j+1) - p(j)) / &
            (0.5d0 * (g(j+1) + g(j)) * DRY_AIR_MASS) &
            * (1.0d0 + 0.5d0 * (sh(j+1) + sh(j)) * (1.0d0 - EPSILON) / EPSILON)

    end do

  end subroutine calculate_rayleigh_tau_MS3



  FUNCTION ray_depol_MS3(wl) RESULT(dpf)

    ! calculates the depolarization factor for a given wavelength
    ! currently just uses band-averaged values for OCO.
    ! finds the closest one for the input wavelength

    double precision, intent(in) :: wl
    double precision             :: dpf

    double precision, dimension(3), parameter :: &
         OCO_DEPOLF = (/ 0.027706d0, 0.027260d0, 0.027213d0 /)
    double precision, dimension(3), parameter :: &
         OCO_WVL    = (/ 0.76d0, 1.61d0, 2.03d0 /)

    integer :: wl_units = 1
    double precision :: wl_um
    integer :: i1(1)

    select case(wl_units)
    case(1) ! microns
       wl_um = wl
    case(2) ! nm
       wl_um = 1d-3 * wl
    case default ! wavenumbers
       wl_um = 1d4/wl
    end select

    i1 = minloc( abs(OCO_WVL - wl_um) )
    dpf = OCO_DEPOLF(i1(1))

  END FUNCTION ray_depol_MS3




  FUNCTION ray_tau_MS3(wl, dpf) RESULT(ray_ext)

    !
    !  USES Bodhaine et al. (1999) parameterization of rayleigh scattering
    !  Bodhaine et al., 1999, Journ. Atm. Oc. Tech., 16, 1854.
    !
    !  I use a simple constant in each wavelength band for the depolarization
    !  factor; this differs only very slightly from Bodhaine.
    !  The full formula is:
    !
    !                  24 PI^3    (ns^2-1)^2   (6+3dpf)
    !  cross-section = ------- *  ---------- * --------
    !                 wl^4 Ns^2   (ns^2+2)^2   (6-7dpf)
    !
    ! where ns is the index of refraction of standard air at some P,T, and
    ! Ns is the number of molecules per cubic meter AT THE SAME P,T.
    ! (the product of those funky terms turns out to be independent of P & T,
    !  according to Lorenz-Lorentz theory); dpf is the depolarization factor,
    !  and wl is the wavelength in meters.
    !
    ! The result is then the Rayleigh cross section in m^2 per molecule of air.
    ! Multiplying then by Avogradro's number yields the cross section in m^2/mol
    ! of (total=moist) air.

    double precision, intent(in), dimension(:) :: wl ! wavelength in units "wl_units"
    double precision, intent(in) :: dpf ! depolarization factor
    double precision, dimension(size(wl))   :: ray_ext  ! output: cross section in m^2/mol of moist air

    ! Fit to Peck & Reeder (1972) index of refraction factor, which was for wavelengths > 0.23 microns.
    ! They had a 4-term fit for the index of refraction of air at T=288.15 K and P = 1013.25 mbar.
    ! I have fitted the index of refraction factor, ((ns^2-1)/(ns^2+2))^2, which is used in the rayleigh
    ! formula.  It is well fit by a function of the form A(1) + A(2)/wl^2 + A(3)/wl^4, with wl in microns.
    double precision, dimension(3), parameter :: Icoeff = (/ 3.302565d-8, 3.7118985d-10, 4.4216716d-12 /)

    ! This is the prefactor in the rayleigh cross section equation: 24*PI^3/(Ns^2)*NA,
    ! such that you must multiply by 1/wl^4, with wl in microns, and the result will be the cross section
    ! in m^2 per mole of gas. This is also for P=1013.25 mbar, T=288.15 K air.
    double precision, parameter :: const = 1.1471954d-24 * NA
    double precision, dimension(size(wl)) :: index_factor, wl2i, wl4i
    double precision, dimension(size(wl))   :: wl_um    ! wavelengths in um
    double precision :: king_factor
    integer :: wl_units = 1 ! 0=cm-1, 1=um, 2=nm

    select case(wl_units)
    case(0) ! wavenumbers
       wl_um = 1d4/wl
    case(1) ! microns
       wl_um = wl
    case(2) ! nm
       wl_um = 1d-3 * wl
    end select

    king_factor = (6.0d0+3.0d0*dpf)/(6.0d0-7.0d0*dpf)
    wl2i = 1.0d0/(wl_um*wl_um)
    wl4i = wl2i * wl2i
    index_factor = Icoeff(1) + Icoeff(2)*wl2i + Icoeff(3)*wl4i ! ((ns^2-1)/(ns^2+2))^2
    ray_ext = const*wl4i*index_factor*king_factor ! rayleigh cross section in m^2/mole

  END FUNCTION ray_tau_MS3




end module Rayleigh_mod
