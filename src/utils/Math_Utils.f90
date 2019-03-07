module math_utils_mod

  use logger_mod, only: logger => master_logger
  use, intrinsic:: ieee_arithmetic, only: ieee_is_nan

  implicit none

  double precision, parameter :: SPEED_OF_LIGHT = 299792458d0
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DEG2RAD = 0.017453292519943295d0
  double precision, parameter :: RAD2DEG = 57.29577951308232d0
  double precision, parameter :: NA = 6.022140857d23
  double precision, parameter :: H2Om = 0.018015422135d0
  double precision, parameter :: DRY_AIR_MASS = 0.0289644d0
  double precision, parameter :: PSURF_PERTURB = 100.0d0 !Pa
  double precision, parameter :: SH_H2O_CONV = DRY_AIR_MASS / H2Om

contains

  pure function calculate_chi2(array1, array2, std_dev, dof) result(chi2)
    implicit none
    double precision, intent(in) :: array1(:), array2(:), std_dev(:)
    integer, intent(in) :: dof

    double precision :: chi2

    chi2 = sum(((array1 - array2) * (array1 - array2)) / (std_dev * std_dev)) / dble(dof)

  end function calculate_chi2

  pure function mean(array) result(mean_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value

    mean_value = sum(array) / dble(size(array))

  end function mean

  pure function std(array) result(std_value)
    implicit none
    double precision, intent(in) :: array(:)
    double precision :: mean_value, std_value

    mean_value = mean(array)
    std_value = sqrt(sum((array - mean_value) * (array - mean_value)) / dble(size(array) - 1))

  end function std


  subroutine pressure_weighting_function(p_levels, psurf, vmrs, pwgts)
    implicit none
    double precision, intent(in) :: p_levels(:), psurf, vmrs(:)
    double precision, intent(inout) :: pwgts(:)

    integer :: i, N
    double precision, allocatable :: hp(:), cbar(:), dp(:)
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

    hp(:) = hp(:) / sum(hp(:))

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


  subroutine pwl_value_1d_v2(nd, xd, yd, ni, xi, yi)

    implicit none

    integer, intent(in) :: nd, ni
    double precision, intent(in) :: xd(:), yd(:), xi(:)
    double precision, intent(inout) :: yi(:)

    double precision :: fac
    integer :: i, j, d, last_d

    last_d = 1
    do i=1, ni
       if (xi(i) <= xd(1)) then
          yi(i) = yd(1)
       else if (xi(i) >= xd(nd)) then
          yi(i) = yd(nd)
       else
          do d=last_d, nd-1
             if ((xi(i) >= xd(d)) .and. (xi(i) <= xd(d+1))) then

                last_d = d

                fac = (xi(i) - xd(d)) / (xd(d+1) - xd(d))
                yi(i) = (1.0d0 - fac) * yd(d) + fac * yd(d+1)

             end if
          end do
       end if
    end do



  end subroutine pwl_value_1d_v2


  pure function searchsorted_dp(x, val, left) result(idx)
    implicit none
    double precision, intent(in) :: x(:), val
    logical, intent(in), optional :: left
    integer :: idx

    logical :: from_left
    integer :: i, L, R, m

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


  subroutine fft_convolution(input, kernel, wl_input, wl_output, output)

    implicit none
    double precision, intent(in) :: input(:), kernel(:), wl_input(:), wl_output(:)
    double precision, intent(inout) :: output(:)

    integer :: power
    integer :: N_input, N_kernel, N_fft
    integer :: input_shift, kernel_shift, pad_space
    integer(4) :: lensav, lenwrk
    integer :: fft_err
    double precision, allocatable :: wsave(:), work(:)
    complex(8), allocatable :: input_fft(:), kernel_fft(:), ifft(:)

    integer :: i, funit
    logical :: found_power


    ! FFT-type convolution using FFTPACK5.1, which is most efficient and
    ! quickest if our two arrays are powers of 2. We thus must zero pad
    ! to a) get the dimension of the array to 2^x (or at least some low
    ! number of prime factors), and b) to mitigate the wrap-around problem
    ! of the FFT method.

    ! Simple procedure:
    ! a) shift the input and ILS kernel such that the center is at
    !    the first quarter of the respective array
    ! b) Perform forward FFT of both input and kernel
    ! c) Multiply the result and back-transform
    ! d) Sample the convolved output at the low-resultion
    !    wavelengths (usually given by the dispersion)

    N_input = size(input)
    N_kernel = size(kernel)

    ! For most efficient FFT, choose closest
    ! power of two for input array.

    found_power = .false.
    power = 1
    do while (.not. found_power)
       if ((dble(N_input)) < (2 ** power)) then
          found_power = .true.
          exit
       end if
       power = power + 1
    end do

    N_fft = 2 ** power

    !pad_space = 1000
    !N_fft = 2 * N_input + pad_space

    ! Look up the rfft1i function - it recommends this size for he
    ! wsave and work arrays
    !lensav = N_fft + int(log(dble(N_fft)) / log(2.0d0)) + 4
    lensav = 2 * N_fft + power + 4
    lenwrk = 2 * N_fft

    allocate(work(lenwrk))
    allocate(wsave(lensav))

    ! Initialise the FFT solver
    call cfft1i(N_fft, wsave, lensav, fft_err)

    ! Transform the input (radiances)
    allocate(input_fft(N_fft))
    allocate(kernel_fft(N_fft))

    input_fft(:) = (0.0d0, 0.0d0)
    kernel_fft(:) = (0.0d0, 0.0d0)

    ! Copy over input data and kernel to complex array
    ! and convert to complex data type with zero imaginary part
    do i=1, N_input
       input_fft(i) = cmplx(input(i), 0.0d0)
    end do

    do i=1, N_kernel
       kernel_fft(i) = cmplx(kernel(i), 0.0d0)
    end do

    input_shift = (N_fft - N_input)/2
    input_shift = -N_input/2
    !kernel_shift = N_kernel / 2
    kernel_shift = maxloc(kernel, dim=1) !N_kernel/2

    input_fft = cshift(input_fft, -input_shift)
    ! The kernel has to be shifted to the center, such
    ! that it's max value is essentially at 1.
    kernel_fft = cshift(kernel_fft, kernel_shift)

    kernel_fft = kernel_fft / sum(kernel_fft)

!!$    open(newunit=funit, file='fft_input_before.dat')
!!$    do i=1, N_input
!!$       write(funit, *) wl_input(i), input(i)
!!$    end do
!!$    close(funit)

    call cfft1f(N_fft, 1, input_fft, N_fft, wsave, &
         lensav, work, lenwrk, fft_err)

    ! Transform the ILS kernel
    call cfft1f(N_fft, 1, kernel_fft, N_fft, wsave, &
         lensav, work, lenwrk, fft_err)

!!$    open(newunit=funit, file='fft_input.dat')
!!$    do i=1, N_fft
!!$       write(funit, *) input_fft(i), kernel_fft(i)
!!$    end do
!!$    close(funit)

    input_fft(:) = input_fft(:) * kernel_fft(:)

!!$    open(newunit=funit, file='fft_mult.dat')
!!$    do i=1, N_fft
!!$       write(funit, *) input_fft(i)
!!$    end do
!!$    close(funit)

    !call cfft1i(N_fft, wsave, lensav, fft_err)
    call cfft1b(N_fft, 1, input_fft, N_fft, wsave, &
         lensav, work, lenwrk, fft_err)

!!$    open(newunit=funit, file='fft_back_raw.dat')
!!$    do i=1, N_fft
!!$       write(funit, *) abs(input_fft(i))
!!$    end do
!!$    close(funit)

    ! Scale the iFFT result and shift it back
    !input_fft = input_fft / sqrt(dble(N_fft) / 5.0d0)
    !input_fft = cshift(input_fft, input_shift)

    input_fft = cshift(input_fft * dble(N_fft), input_shift)

    ! Interpolate the values at the wavelengths given by
    ! the dispersion relation. Note that we are only
    ! using the REAL part of the iFFT result. In Fortran,
    ! dble/real of a complex type returns the real part,
    ! rather than the absolute value.

    call pwl_value_1d(N_input, &
         wl_input(:), dble(input_fft(1:N_input)), &
         size(output), &
         wl_output(:), output(:))

!!$    open(newunit=funit, file='fft_back_hires.dat')
!!$    do i=1, N_input
!!$       write(funit, *) wl_input(i), abs(input_fft(i))
!!$    end do
!!$    close(funit)
!!$
!!$    open(newunit=funit, file='fft_back.dat')
!!$    do i=1, size(output)
!!$       write(funit, *) wl_output(i), output(i)
!!$    end do
!!$    close(funit)
!!$
!!$    stop 1


  end subroutine fft_convolution


  subroutine oco_type_convolution2(wl_input, input, wl_kernels, kernels, &
       wl_output, output, success)

    !! This is an implementation of the OCO-type ILS application, which ins't
    !! quite a convolution since we have to treat every pixel with a different
    !! ILS. wl_output is the DESIRED output wavelength grid

    ! High-resolution bits
    double precision, intent(in) :: wl_input(:), input(:)
    ! "Convolution" wavelengths and kernels
    double precision, intent(in) :: wl_kernels(:,:), kernels(:,:)
    double precision, intent(in) :: wl_output(:)
    double precision, intent(inout) :: output(:)
    logical, intent(out) :: success

    character(len=*), parameter :: fname = "oco_type_convolution"
    integer :: N_pix, N_ils_pix, N_wl, N_this_wl
    integer :: idx_pix, idx_hires_closest, idx_hires_ILS_min, idx_hires_ILS_max
    integer :: kernel_idx_min, kernel_idx_max
    double precision :: ILS_delta_min, ILS_delta_max
    double precision :: ILS_wl_spacing
    double precision :: wl_diff

    double precision, allocatable :: ILS_upsampled(:), input_upsampled(:)
    double precision :: time_start, time_stop
    integer :: i, funit

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

       N_this_wl = idx_hires_ILS_max - idx_hires_ILS_min
       allocate(ILS_upsampled(N_this_wl + 1))
       ILS_upsampled = 0.0d0

       call pwl_value_1d( &
            N_ils_pix, &
            wl_output(idx_pix) + wl_kernels(:, idx_pix), kernels(:, idx_pix), &
            N_this_wl, &
            wl_input(idx_hires_ILS_min:idx_hires_ILS_max), ILS_upsampled(:))

       output(idx_pix) = dot_product(ILS_upsampled(:), input(idx_hires_ILS_min:idx_hires_ILS_max)) &
            / sum(ILS_upsampled)

       deallocate(ILS_upsampled)

    end do

    success = .true.


  end subroutine oco_type_convolution2

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


  function percentile(x, perc)

    !! Calculates the percentile "perc" of a double precision array "x"

    implicit none
    double precision :: percentile

    double precision, dimension(:), intent(in) :: x
    double precision, intent(in) :: perc
    double precision, dimension(:), allocatable :: x_sort
    double precision :: position, position_remainder
    integer :: position_int
    character(len=*), parameter :: fname = "percentile"

    if ((perc < 0) .or. (perc > 100)) then
       call logger%fatal(fname, "Percentile must be between 0 and 100!")
       stop 1
    end if

    allocate(x_sort(size(x)))
    x_sort(:) = x
    call combsort(x_sort)

    position = (perc / 100) * size(x)

    position_int = floor(position)
    position_remainder = position - floor(position)

    percentile = (x_sort(position_int) * (1.0d0 - position_remainder)) + &
         (x_sort(position_int + 1) * position_remainder)

  end function percentile

  subroutine combsort(a)
    ! Taken from Rosettacode

    double precision, intent(in out) :: a(:)
    double precision :: temp
    integer :: i, gap
    logical :: swapped = .true.

    gap = size(a)
    do while (gap > 1 .or. swapped)
       gap = gap / 1.3
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
