!> @brief Math Utils module
!> @file math_utils.f90
!> @author Peter Somkuti
!>
!> @details
!> Here we collect 'math-y' utils/subroutines that are used by
!> various other modules, i.e. calculation of a mean, std or the
!> chi2.

module math_utils_mod

  ! User modules
  use logger_mod, only: logger => master_logger
  use misc_utils_mod, only: get_cpu_time

  ! System modules
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use :: iso_fortran_env

  implicit none

  ! Various constants
  double precision, parameter :: SPEED_OF_LIGHT = 299792458d0
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DEG2RAD = 0.017453292519943295d0
  double precision, parameter :: RAD2DEG = 57.29577951308232d0
  double precision, parameter :: NA = 6.022140857d23
  double precision, parameter :: H2Om = 0.018015422135d0
  double precision, parameter :: DRY_AIR_MASS = 0.0289644d0
  double precision, parameter :: PSURF_PERTURB = 100.0d0 !Pa
  double precision, parameter :: SH_H2O_CONV = DRY_AIR_MASS / H2Om
  double precision, parameter :: Rd = 8.314472d0 / DRY_AIR_MASS
  double precision, parameter :: EPSILON = H2Om / DRY_AIR_MASS
  double precision, parameter :: PLANCK_CONST = 6.62607015d-34
  double precision, parameter :: BOLTZMANN_CONST = 1.38064852d-23
  double precision, parameter :: SOLAR_RADIUS = 695.660d6 ! [m]
  double precision, parameter :: AU_UNIT = 149597870d3 ! [m]
  double precision, parameter :: EARTH_EQUATORIAL_RADIUS = 6378178.0d0 ! [M]


contains


  !> @brief Calculates reduced CHI2 statistic for two spectra
  !> @param array1 One of two spectra
  !> @param array2 The other spectrum
  !> @param std_dev Standard deviation / noise per pixel
  !> @param dof Degrees of freedom (number of SV elements - 1 usually)
  function calculate_chi2(array1, array2, std_dev, dof) result(chi2)
    implicit none
    double precision, intent(in) :: array1(:), array2(:), std_dev(:)
    integer, intent(in) :: dof

    double precision :: chi2

    chi2 = sum(((array1 - array2) * (array1 - array2)) / (std_dev * std_dev)) / dble(dof)
  end function calculate_chi2


  !> @brief Simple function to calculate the mean for double-precision arrays
  !> @param array Array you want the mean of
  !> @param mean_value Arithmetic mean of array
  function mean(array) result(mean_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value

    mean_value = sum(array) / dble(size(array))
  end function mean

  !> @brief Simple function to calculate the standard deviation for double-precision arrays
  !> @param array Array you want the std of
  !> @param std_value (Unbiased) standard deviation of array (Bessel correction)
  function std(array) result(std_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value, std_value

    mean_value = mean(array)
    std_value = sqrt(sum((array - mean_value) * (array - mean_value)) / dble(size(array) - 1))
  end function std


  !> @brief Tropopause pressure calculation (C) T. Reichler
  !> @detail This code was copied from http://www.inscc.utah.edu/~reichler/research/projects/TROPO/code.txt
  !>         and only adjusted for data types (real to DP), and reversed layer ordering.
  !>         Values used in Leicester preprocessing:
  !>         plimu = 45000 Pa
  !>         pliml = 7500 Pa
  !>         gamma = -0.002 K/m
  subroutine twmo(level, t, p, plimu, pliml, gamma, trp)

    implicit none
    integer,intent(in)                  :: level
    double precision, intent(in)        :: t(:), p(:)
    double precision, intent(in)        :: plimu, pliml, gamma
    double precision, intent(out)       :: trp

    double precision, parameter         :: kap=0.286
    double precision, parameter         :: faktor = -9.81/287.0
    double precision, parameter         :: deltaz = 2000.0
    double precision, parameter         :: ka1=kap-1.

    double precision                    :: pmk, pm, a, b, tm, dtdp, dtdz
    double precision                    :: ag, bg, ptph
    double precision                    :: pm0, pmk0, dtdz0
    double precision                    :: p2km, asum, aquer
    double precision                    :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
    integer                             :: icount, jj
    integer                             :: j

    trp=-99.0                           ! negative means not valid
    do j=level, 2, -1

       ! dt/dz
       pmk= .5 * (p(j-1)**kap+p(j)**kap)
       pm = pmk**(1/kap)
       a = (t(j-1)-t(j))/(p(j-1)**kap-p(j)**kap)
       b = t(j)-(a*p(j)**kap)
       tm = a * pmk + b
       dtdp = a * kap * (pm**ka1)
       dtdz = faktor*dtdp*pm/tm

       ! dt/dz valid?
       if (j.eq.level)    go to 999     ! no, start level, initialize first
       if (dtdz.le.gamma) go to 999     ! no, dt/dz < -2 K/km
       if (pm.gt.plimu)   go to 999     ! no, too low

       ! dtdz is valid, calculate tropopause pressure
       if (dtdz0.lt.gamma) then
          ag = (dtdz-dtdz0) / (pmk-pmk0)
          bg = dtdz0 - (ag * pmk0)
          ptph = exp(log((gamma-bg)/ag)/kap)
       else
          ptph = pm
       endif

       if (ptph.lt.pliml) go to 999
       if (ptph.gt.plimu) go to 999

       ! 2nd test: dtdz above 2 km must not exceed gamma
       p2km = ptph + deltaz*(pm/tm)*faktor          ! p at ptph + 2km
       asum = 0.0                                   ! dtdz above
       icount = 0                                   ! number of levels above

       ! test until apm < p2km
       do jj=j, 2, -1

          pmk2 = .5 * (p(jj-1)**kap+p(jj)**kap)    ! p mean ^kappa
          pm2 = pmk2**(1/kap)                      ! p mean
          if(pm2.gt.ptph) go to 110                ! doesn't happen
          if(pm2.lt.p2km) go to 888                ! ptropo is valid

          a2 = (t(jj-1)-t(jj))                     ! a
          a2 = a2/(p(jj-1)**kap-p(jj)**kap)
          b2 = t(jj)-(a2*p(jj)**kap)               ! b
          tm2 = a2 * pmk2 + b2                     ! T mean
          dtdp2 = a2 * kap * (pm2**(kap-1))        ! dt/dp
          dtdz2 = faktor*dtdp2*pm2/tm2
          asum = asum+dtdz2
          icount = icount+1
          aquer = asum/float(icount)               ! dt/dz mean

          ! discard ptropo ?
          if (aquer.le.gamma) go to 999           ! dt/dz above < gamma

110       continue
       enddo                           ! test next level

888    continue                        ! ptph is valid
       trp = ptph
       return

999    continue                        ! continue search at next higher level
       pm0 = pm
       pmk0 = pmk
       dtdz0  = dtdz

    enddo

    ! no tropopouse found
    return
  end subroutine twmo


  !> @brief Calculates the O'Dell-ian pressure weighting function
  !
  !> @param N Number of active levels
  !> @param p_levels Pressure level array (levels)
  !> @param psurf Surface pressure
  !> @param q Specific humidity (layers)
  !> @param g Gravity (layers)
  !> @param pwgts Pressure weights (levels)
  !
  !> @details For details of the calculation go and read O'Dell et al. (2012 ACOS paper)
  subroutine pressure_weighting_function(N, p_levels, psurf, q, g, pwgts)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: p_levels(:)
    double precision, intent(in) :: psurf
    double precision, intent(in) :: q(:)
    double precision, intent(in) :: g(:)
    double precision, intent(inout) :: pwgts(:)

    ! Loop and size variables
    integer :: i
    ! H-prime, c-bar and delta pressure
    double precision, allocatable :: hp(:), cbar(:), dp(:)
    ! Interpolation factors
    double precision :: f, fs

    pwgts(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

    if (size(pwgts) /= N) return

    ! First, calculate h-prime (Equation A4 in O'Dell 2012)
    allocate(hp(N-1))
    allocate(cbar(N-1))
    allocate(dp(N-1))

    do i=1, N-1
       ! Loop over levels starting from top down to
       ! the surface layer.
       cbar(i) = (1.0d0 - q(i)) / (g(i) * DRY_AIR_MASS)
       dp(i) = p_levels(i+1) - p_levels(i)
       hp(i) = cbar(i) * dp(i)
    end do

    ! Normalize h-prime
    hp(:) = hp(:) / sum(hp(:))

    ! F is set to 0.5 for all levels, I guess this could be changed
    ! if someone really wanted to..
    f = 0.5d0
    fs = (psurf - p_levels(N-1)) &
         / (p_levels(N) - p_levels(N-1))

    do i=1, N

       if (i == 1) then
          pwgts(i) = (1.0d0 - f) * hp(i)
       else if (i == (N-1)) then
          pwgts(i) = f * hp(i-1) + (1.0d0 - fs * f) * hp(i)
       else if (i == N) then
          pwgts(i) = fs * f * hp(i-1)
       else if ((i > 1) .and. (i < N - 1)) then
          pwgts(i) = f * hp(i-1) + (1.0d0 - f) * hp(i)
       end if

    end do

    if (abs(sum(pwgts(:)) - 1.0d0) > 1.0d-5) then
       call logger%warning("pressure_weighting_function", &
            "Pressure weighting function does not sum up to 1.0")
    end if
  end subroutine pressure_weighting_function


  !> @brief Piecewise linear interpolation, based on the code
  !> by John Burkardt - should be roughly the same speed.
  !
  ! NOTE: THIS ASSUMES AN ORDERED ARRAY! THE ORDERING IS NOT CHECKED!
  !
  !> @param nd Number of data points given
  !> @param xd x-coordinates of data
  !> @param yd y-coordinates of data
  !> @param ni Number of points to be calculated
  !> @param xi x-coordinates where interpolation should be performed
  !> @param yi Interpolated values
  subroutine pwl_value_1d_v2(nd, xd, yd, ni, xi, yi)

    implicit none

    integer, intent(in) :: nd, ni
    double precision, intent(in) :: xd(nd)
    double precision, intent(in) :: yd(nd)
    double precision, intent(in) :: xi(ni)
    double precision, intent(inout) :: yi(ni)

    double precision :: fac
    integer :: i, d, last_d

    last_d = 1
    do i=1, ni

       if (xi(i) <= xd(1)) then
          ! Edge case - if smaller than smallest data point,
          ! set to that value.
          yi(i) = yd(1)
       else if (xi(i) >= xd(nd)) then
          ! Edge case - if larger than largest data point,
          ! set to that value.
          yi(i) = yd(nd)
       else
          do d=last_d, nd-1
             if ((xi(i) >= xd(d)) .and. (xi(i) <= xd(d+1))) then

                ! Since we are using a sorted array, we can
                ! skip searching the first d points, as we know
                ! the value is not going to be found there.
                ! This provides a significant speed-up by ~10x
                last_d = d

                fac = (xi(i) - xd(d)) / (xd(d+1) - xd(d))
                yi(i) = (1.0d0 - fac) * yd(d) + fac * yd(d+1)
                exit

             end if
          end do
       end if
    end do

  end subroutine pwl_value_1d_v2


  !> @brief See Numpy's searchsorted function, using binary search
  !> @param x Array to be searched for
  !> @param val Value, whose position in the array needs determining
  !> @param left Left insert? (or right if set to .false.)
  !> @param idx The index at which val would be inserted in x
  pure function searchsorted_dp(x, val, left, user_L) result(idx)
    implicit none
    double precision, intent(in) :: x(:), val
    logical, intent(in), optional :: left
    integer, intent(in), optional :: user_L
    integer :: idx

    logical :: from_left
    integer :: L, R, m

    ! Insert to the left of value is the standard
    if (.not. present(left)) then
       from_left = .true.
    else
       from_left = left
    end if

    if (present(user_L)) then
       L = user_L
    else
       L = 1
    end if

    R = size(x)
    idx = -1

    do while (L <= R)
       ! We want floor here! That's exactly what Fortran
       ! integer division does by default.
       m = (L + R) / 2

       if (m >= size(x) - 1) then
          idx = m
          return
       end if

       if (m == 0) then
          idx = 1
          return
       end if

       if ((x(m) <= val) .and. (x(m+1) >= val)) then
          ! Found!
          idx = m
          if (.not. from_left) idx = idx+1
          return
       else if (x(m) < val) then
          L = m + 1
       else if (x(m) > val) then
          R = m - 1
       else
          idx = m
       end if
    end do

  end function searchsorted_dp



  !> @brief Convolve function for OCO-type instruments with different
  !> ILS's for every detector pixel
  !> @param wl_input Wavelengths of high-res radiance array
  !> @param input High-res radiance array
  !> @param wl_kernels Wavelength for ILS per pixel for this footprint/band
  !> @param kernels ILS relative response per pixel for this footprint/band
  !> @param wl_output Requested wavelength grid for convolved radiance
  !> @param output Convolved radiance
  !> @param success Did it all go well?
  !> @note
  !> I've experimented a bit with making this function faster, since the convolution
  !> and the gas OD calculation are essentially the top two time-consuming portions
  !> of the physical algorithm in Beer-Lambert mode. Another way of doing this was to
  !> re-grid all ILS tables to a common high-resolution grid, so that we would not need to
  !> interpolate the ILS onto the high-res grid of the current pixel. This however seemed to
  !> create issues when the chosen high-res spacing was not fine enough, since the ILS wl spacing
  !> gets finer the closer you are at the center (at least for OCO-2/3). Hence I reverted to
  !> this particular option, where the ILS is interpolated to the high-res wavelength grid
  !> for every pixel at every calculation. While this makes it somewhat slower, it also
  !> seemed to have eliminated the issue of bad results. These bad results were non-trivially
  !> seen as stripy patterns which were related to time-of-day and doppler shift.
  !> Since the doppler shift changes the wavelength grid, one can end up with misaligned
  !> spectra if the ILS convolution (however it is done) does not work accordingly and shifts
  !> the line cores around..
  subroutine oco_type_convolution(wl_input, input, wl_kernels, kernels, &
       wl_output, output, success)

    ! This is an implementation of the OCO-type ILS application, which ins't
    ! quite a convolution since we have to treat every pixel with a different
    ! ILS. wl_output is the DESIRED output wavelength grid

    ! High-resolution bits
    double precision, contiguous, intent(in) :: wl_input(:)
    double precision, contiguous, intent(in) :: input(:)

    ! "Convolution" wavelengths and kernels
    real, contiguous, intent(in) :: wl_kernels(:,:)
    real, contiguous, intent(in) :: kernels(:,:)

    double precision, contiguous, intent(in) :: wl_output(:)
    double precision, contiguous, intent(inout) :: output(:)
    logical, intent(out) :: success

    character(len=*), parameter :: fname = "oco_type_convolution"

    ! These should be constant, and helps the compiler to optimize the loop
    integer :: N_pix
    integer :: N_ils_pix
    integer :: N_wl

    integer :: N_this_wl
    integer :: idx_pix
    integer :: idx_hires_ILS_min
    integer :: idx_hires_ILS_max

    double precision, allocatable :: tmp1(:), tmp2(:)

    real, allocatable :: wl_kernels_absolute(:,:)
    double precision, allocatable :: ILS_upsampled(:)

    double precision :: cpu_start, cpu_stop
    character(len=999) :: tmp_str

    double precision, external :: ddot
    double precision, external :: dasum

    cpu_start = get_cpu_time()

    N_wl = size(wl_input)
    N_pix = size(wl_output)
    N_ils_pix = size(wl_kernels, 1)

    success = .false.

    if (N_pix /= size(wl_kernels, 2)) then
       call logger%fatal(fname, "wl_kernels or wl_output have incompatible sizes.")
       stop 1
    end if

    if (N_pix /= size(kernels, 2)) then
       call logger%fatal(fname, "kernels or wl_output have incompatible sizes.")
       stop 1
    end if

    if (any(ieee_is_nan(wl_input))) then
       call logger%fatal(fname, "NaN(s) in wl_input")
       return
    end if

    if (any(ieee_is_nan(wl_output))) then
       call logger%fatal(fname, "NaN(s) in wl_output")
       return
    end if

    if (any(ieee_is_nan(wl_kernels))) then
       call logger%fatal(fname, "NaN(s) in wl_kernels")
       return
    end if

    ! Note the ILS boundary in wavelength space. wl_kernels spans usually
    ! some range from -lambda to +lambda. We thus bring them onto an absolute
    ! wavelength grid, dependent on the detector pixel position.
    allocate(wl_kernels_absolute, mold=wl_kernels)
    do idx_pix = 1, N_pix
        wl_kernels_absolute(:, idx_pix) = wl_kernels(:, idx_pix) + wl_output(idx_pix)
    end do

    ! Main loop over all (instrument) output pixels. Note that idx_pix is
    ! RELATIVE to the array that is supplied here, which means you need to
    ! make sure that you only pass arrays here that in wavelength (detector pixel)
    ! space are at the same positions as the output wavelengths wl_output.

    allocate(tmp1(size(wl_kernels_absolute, 1)))
    allocate(tmp2(size(kernels, 1)))

    allocate(ILS_upsampled(N_wl))
    do idx_pix=1, N_pix

       idx_hires_ILS_min = searchsorted_dp(wl_input, dble(wl_kernels_absolute(1, idx_pix)))
       ! When doing the binary search for the upper WL index, we can safely
       ! set the "left" search limit to the other index.
       idx_hires_ILS_max = searchsorted_dp(wl_input, &
            dble(wl_kernels_absolute(N_ils_pix, idx_pix)), &
            user_L=idx_hires_ILS_min)

       if ((idx_hires_ILS_min < 1) .or. (idx_hires_ILS_min > size(wl_input))) then
          call logger%error(fname, "idx_hires_ILS_min out of bounds.")
          return
       end if

       if ((idx_hires_ILS_max < 1) .or. (idx_hires_ILS_max > size(wl_input))) then
          call logger%error(fname, "idx_hires_ILS_max out of bounds.")
          return
       end if

       N_this_wl = idx_hires_ILS_max - idx_hires_ILS_min + 1


       ! We need this crutch, because ILS tables are now in real(4), but the
       ! function itself requires real(8)
       tmp1(:) = wl_kernels_absolute(:, idx_pix)
       tmp2(:) = kernels(:, idx_pix)

       ! Interpolate the ILS at the actual hi-res grid
       call pwl_value_1d_v2( &
            N_ils_pix, &
            tmp1, tmp2, &
            N_this_wl, &
            wl_input(idx_hires_ILS_min:idx_hires_ILS_max), ILS_upsampled(1:N_this_wl))

       ! Use BLAS for dot product and sum - gives a nice speed-up!
       output(idx_pix) = ddot( &
            N_this_wl, ILS_upsampled(1:N_this_wl), 1, &
            input(idx_hires_ILS_min:idx_hires_ILS_max), 1) &
            / dasum(N_this_wl, ILS_upsampled(1:N_this_wl), 1)

    end do
    deallocate(ILS_upsampled)

    success = .true.

    cpu_stop = get_cpu_time()
    write(tmp_str, '(A, F10.7)') "Convolution time (s): ", cpu_stop - cpu_start
    call logger%debug(fname, trim(tmp_str))

  end subroutine oco_type_convolution

  !> @brief Inverts a quadratic matrix using DGETRF/DGETRI
  !> @param mat_in Matrix to be inverted
  !> @param mat_out Inverted matrix
  !> @param success False if DGETRF or DGETRI fail
  subroutine invert_matrix(mat_in, mat_out, success)

    double precision, dimension(:,:), intent(in) :: mat_in
    double precision, dimension(:,:), intent(out) :: mat_out
    logical, intent(out):: success

    double precision, dimension(size(mat_in,1)) :: work  ! work array for LAPACK
    integer, dimension(size(mat_in,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    mat_out = mat_in
    n = size(mat_in, 1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, mat_out, n, ipiv, info)

    if (info /= 0) then
       call logger%fatal("invert_matrix", "Matrix is numerically singular!")
       write(*,*) "DGETRF Error Code: ", info
       success = .false.
       return
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, mat_out, n, ipiv, work, n, info)

    if (info /= 0) then
       call logger%fatal("invert_matrix", "Matrix inversion failed!")
       success = .false.
       return
    end if

    success = .true.

  end subroutine invert_matrix


  !> @brief Calculates the perc'th percentile for array x
  !> @param x Array you want the perc'th percentile of
  !> @param perc Percentile to be calculated
  function percentile(x, perc)

    implicit none
    double precision :: percentile

    double precision, dimension(:), intent(in) :: x
    double precision, intent(in) :: perc
    double precision, dimension(:), allocatable :: x_sort
    double precision :: position, position_remainder
    integer :: position_int
    !character(len=*), parameter :: fname = "percentile"

    if ((perc < 0) .or. (perc > 100)) then
       percentile = -9999
       return
    end if

    allocate(x_sort(size(x)))
    x_sort(:) = x
    call combsort(x_sort)

    position = (perc / 100) * size(x)

    position_int = floor(position)
    position_remainder = position - floor(position)

    if ((position_int < 1) .or. (position_int > size(x))) then
       percentile = -9999
       return
    end if

    percentile = (x_sort(position_int) * (1.0d0 - position_remainder)) + &
         (x_sort(position_int + 1) * position_remainder)

  end function percentile

  !> @brief Combsort routine, taken from Rosettacode
  !> @param a Array to be sorted
  subroutine combsort(a)

    double precision, intent(in out) :: a(:)
    double precision :: temp
    integer :: i, gap
    logical :: swapped = .true.

    gap = size(a)
    do while (gap > 1 .or. swapped)
       gap = int(gap / 1.3)
       if (gap < 1) gap = 1
       swapped = .false.
       do i = 1, size(a)-gap
          if (a(i) > a(i+gap)) then
             temp = a(i)
             a(i) = a(i+gap)
             a(i+gap) = temp;
             swapped = .true.
          end if
       end do
    end do

  end subroutine combsort


end module math_utils_mod
