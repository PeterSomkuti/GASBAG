!> @brief Math Utils module
!> @file math_utils.f90
!> @author Peter Somkuti
!>
!! Here we collect 'math-y' utils/subroutines that are used by
!! various other modules, i.e. calculation of a mean, std or the
!! chi2.

module math_utils_mod

  ! User modules
  use logger_mod, only: logger => master_logger

  ! System modules
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
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
  double precision, parameter :: PLANCK_CONST = 6.62607015d-34
  double precision, parameter :: BOLTZMANN_CONST = 1.38064852d-23
  double precision, parameter :: SOLAR_RADIUS = 695.660d6 ! [m]
  double precision, parameter :: AU_UNIT = 149597870d3 ! [m]

contains

  !> @brief Calculates reduced CHI2 statistic for two spectra
  !> @param array1 One of two spectra
  !> @param array2 The other spectrum
  !> @param std_dev Standard deviation / noise per pixel
  !> @param dof Degrees of freedom (number of SV elements - 1 usually)
  pure function calculate_chi2(array1, array2, std_dev, dof) result(chi2)
    implicit none
    double precision, intent(in) :: array1(:), array2(:), std_dev(:)
    integer, intent(in) :: dof

    double precision :: chi2

    chi2 = sum(((array1 - array2) * (array1 - array2)) / (std_dev * std_dev)) / dble(dof)
  end function calculate_chi2


  !> @brief Simple function to calculate the mean for double-precision arrays
  !> @param array Array you want the mean of
  !> @param mean_value Arithmetic mean of array
  pure function mean(array) result(mean_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value

    mean_value = sum(array) / dble(size(array))
  end function mean

  !> @brief Simple function to calculate the standard deviation for double-precision arrays
  !> @param array Array you want the std of
  !> @param std_value (Unbiased) standard deviation of array (Bessel correction)
  pure function std(array) result(std_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value, std_value

    mean_value = mean(array)
    std_value = sqrt(sum((array - mean_value) * (array - mean_value)) / dble(size(array) - 1))
  end function std


  !> @brief Calculates the O'Dell-ian pressure weighting function
  !
  !> @param p_levels Pressure level array
  !> @param psurf Surface pressure
  !> @param vrms Volume mixing ratio array for this gas
  !> @param pwgts Pressure weights
  !
  !> @detail For details of the calculation go and read O'Dell et al. (2012 ACOS paper)
  subroutine pressure_weighting_function(p_levels, psurf, vmrs, pwgts)
    implicit none
    double precision, intent(in) :: p_levels(:)
    double precision, intent(in) :: psurf
    double precision, intent(in) :: vmrs(:)
    double precision, intent(inout) :: pwgts(:)

    ! Loop and size variables
    integer :: i, N
    ! H-prime, c-bar and delta pressure
    double precision, allocatable :: hp(:), cbar(:), dp(:)
    ! Interpolation factors
    double precision :: f, fs

    N = size(p_levels)

    ! First, calculate h-prime (Equation A4 in O'Dell 2012)
    allocate(hp(N-1))
    allocate(cbar(N-1))
    allocate(dp(N-1))


    do i=1, N-1
       ! Loop over levels starting from top down to
       ! the surface layer.
       cbar(i) = 0.5d0 * (vmrs(i) + vmrs(i+1))
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
       else
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
    double precision, intent(in) :: xd(nd), yd(nd), xi(ni)
    double precision, intent(inout) :: yi(ni)

    double precision :: fac
    integer :: i, d, last_d

    last_d = 1
    do i=1, ni
       if (xi(i) <= xd(1)) then
          yi(i) = yd(1)
       else if (xi(i) >= xd(nd)) then
          yi(i) = yd(nd)
       else
          do d=last_d, nd-1
             if ((xi(i) >= xd(d)) .and. (xi(i) <= xd(d+1))) then

                ! Since we are using a sorted array, we can
                ! skip searching the first d points, as we know
                ! the value is not going to be found there.
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
  pure function searchsorted_dp(x, val, left) result(idx)
  implicit none
  double precision, intent(in) :: x(:), val
  logical, intent(in), optional :: left
  integer :: idx

  logical :: from_left
  integer :: L, R, m

  ! Insert to the left of value is the standard
  if (.not. present(left)) then
       from_left = .true.
    else
       from_left = left
    end if

    L = 1
    R = size(x)
    idx = -1

    do while (L <= R)
       ! We want floor here! That's exactly what Fortran
       ! integer division does by default.
       m = (L + R) / 2

       if (m == size(x)) then
          idx = m
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
  subroutine oco_type_convolution(wl_input, input, wl_kernels, kernels, &
       wl_output, output, success)

    ! This is an implementation of the OCO-type ILS application, which ins't
    ! quite a convolution since we have to treat every pixel with a different
    ! ILS. wl_output is the DESIRED output wavelength grid

    ! Notes Peter Somkuti:
    ! I've experimented a bit with making this function faster, since the convolution
    ! and the gas OD calculation are essentially the top two time-consuming portions
    ! of the physical algorithm. Another way of doing this was to re-grid all ILS tables
    ! to a common high-resolution grid, so that we would not need to interpolate the ILS
    ! onto the high-res grid of the current pixel. This however seemed to create issues
    ! when the chosen high-res spacing was not fine enough, since the ILS wl spacing gets
    ! finer the closer you are at the center (at least for OCO-2/3). Hence I reverted to
    ! this particular option, where the ILS is interpolated to the high-res wavelength grid
    ! for every pixel at every calculation. While this makes it somewhat slower, it also
    ! seemed to have eliminated the issue of bad results. These bad results were non-trivially
    ! seen as stripy patterns which were somewhat related to time-of-day and doppler shift.
    ! Since the doppler shift changes the wavelength grid, one can end up with misaligned
    ! spectra if the ILS convolution (however it is done) does work accordingly and shifts
    ! the line cores around..

    ! High-resolution bits
    double precision, intent(in) :: wl_input(:), input(:)
    ! "Convolution" wavelengths and kernels
    double precision, intent(in) :: wl_kernels(:,:), kernels(:,:)
    double precision, intent(in) :: wl_output(:)
    double precision, intent(inout) :: output(:)
    logical, intent(out) :: success

    character(len=*), parameter :: fname = "oco_type_convolution"
    integer :: N_pix, N_ils_pix, N_wl, N_this_wl
    integer :: idx_pix, idx_hires_ILS_min, idx_hires_ILS_max
    double precision :: ILS_delta_min, ILS_delta_max
    double precision :: ILS_wl_spacing

    double precision, allocatable :: ILS_upsampled(:)

    N_wl = size(wl_input)
    N_pix = size(wl_output)
    N_ils_pix = size(wl_kernels, 1)

    success = .false.

    ILS_wl_spacing = wl_input(2) - wl_input(1)

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

    ! Main loop over all (instrument) output pixels. Note that idx_pix is
    ! RELATIVE to the array that is supplied here, which means you need to
    ! make sure that you only pass arrays here that in wavelength (detector pixel)
    ! space are at the same positions as the output wavelengths wl_output.

    do idx_pix=1, N_pix

       ! Note the ILS boundary in wavelength space. wl_kernels spans usually
       ! some range from -lambda to +lambda.
       ILS_delta_min = wl_output(idx_pix) + wl_kernels(1, idx_pix)
       ILS_delta_max = wl_output(idx_pix) + wl_kernels(N_ils_pix, idx_pix)

       idx_hires_ILS_min = searchsorted_dp(wl_input, ILS_delta_min)
       idx_hires_ILS_max = searchsorted_dp(wl_input, ILS_delta_max)

       if ((idx_hires_ILS_min < 1) .or. (idx_hires_ILS_min > size(wl_input))) then
          call logger%error(fname, "idx_hires_ILS_min out of bounds.")
          return
       end if

       if ((idx_hires_ILS_max < 1) .or. (idx_hires_ILS_max > size(wl_input))) then
          call logger%error(fname, "idx_hires_ILS_max out of bounds.")
          return
       end if


       N_this_wl = idx_hires_ILS_max - idx_hires_ILS_min
       allocate(ILS_upsampled(N_this_wl + 1))
       ILS_upsampled = 0.0d0

       call pwl_value_1d_v2( &
            N_ils_pix, &
            wl_output(idx_pix) + wl_kernels(:, idx_pix), kernels(:, idx_pix), &
            N_this_wl, &
            wl_input(idx_hires_ILS_min:idx_hires_ILS_max), ILS_upsampled(:))

       output(idx_pix) = dot_product(ILS_upsampled(:), input(idx_hires_ILS_min:idx_hires_ILS_max)) &
            / sum(ILS_upsampled)

       deallocate(ILS_upsampled)

    end do

    success = .true.
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
    mat_out(:,:) = mat_in(:,:)
    n = size(mat_in,1)

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
    character(len=*), parameter :: fname = "percentile"

    if ((perc < 0) .or. (perc > 100)) then
       !call logger%fatal(fname, "Percentile must be between 0 and 100!")
       !stop 1
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
