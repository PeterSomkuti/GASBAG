module math_utils_mod

    use logger_mod, only: logger => master_logger

    implicit none

    double precision, parameter :: SPEED_OF_LIGHT = 299792458d0
    double precision, parameter :: PI = 3.141592653589793d0
    double precision, parameter :: DEG2RAD = 0.017453292519943295d0
    double precision, parameter :: RAD2DEG = 57.29577951308232d0
    double precision, parameter :: NA = 6.022140857d23
    double precision, parameter :: H2Om = 0.018015422135d0
    double precision, parameter :: DRY_AIR_MASS = 0.02897d0
    double precision, parameter :: PSURF_PERTURB = 100.0d0 !Pa

contains


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

       if ((x(m) < val) .and. (x(m+1) > val)) then
          ! Found!
          idx = m
          if (.not. from_left) idx = m+1
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


  subroutine oco_type_convolution(wl_input, input, wl_kernels, kernels, &
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
    integer :: N_pix, N_ils_pix, N_wl
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

    ILS_wl_spacing = wl_input(2) - wl_input(1)

    if (N_pix /= size(wl_kernels, 2)) then
       call logger%fatal(fname, "wl_kernels or wl_output have incompatible sizes.")
       stop 1
    end if

    if (N_pix /= size(kernels, 2)) then
       call logger%fatal(fname, "kernels or wl_output have incompatible sizes.")
       stop 1
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

       ! Again, we need to make sure that these wavelengths here ILS_delta_* are
       ! multiples of the hires grid
       ILS_delta_min = ILS_wl_spacing * ceiling(ILS_delta_min / ILS_wl_spacing)
       !ILS_delta_max = ILS_wl_spacing * ceiling(ILS_delta_max / ILS_wl_spacing)

       ! Find out, which index of the high-res input is closest to the
       ! detector pixel, and also which is the center pixel

       idx_hires_ILS_min = -1
       !idx_hires_ILS_max = -1

       if (wl_input(1) > ILS_delta_min) then
          call logger%warning(fname, "ILS produdes out of lower wavelength range!")
          kernel_idx_min = searchsorted_dp(wl_kernels(:, idx_pix) + wl_output(idx_pix), wl_input(1))
          kernel_idx_max = N_ils_pix
          idx_hires_ILS_min = 1
          idx_hires_ILS_max = 1 + kernel_idx_max - kernel_idx_min
       else if (wl_input(N_wl) < ILS_delta_max) then
          call logger%warning(fname, "ILS produdes out of higher wavelength range!")
          kernel_idx_max = searchsorted_dp(wl_kernels(:, idx_pix) + wl_output(idx_pix), wl_input(N_wl))
          kernel_idx_min = 1
          idx_hires_ILS_min = N_wl - (kernel_idx_max - kernel_idx_min) - 1
          idx_hires_ILS_max = N_wl
       else
          ! Normal, ILS bounds are within high-res wavelength array
          idx_hires_ILS_min = searchsorted_dp(wl_input(:), ILS_delta_min)
          idx_hires_ILS_max = idx_hires_ILS_min + N_ils_pix - 1
          kernel_idx_min = 1
          kernel_idx_max = N_ils_pix
       end if

!!$       do i=1, size(wl_input)-1
!!$          ! The ILS should be on a high-resolution grid at the same spacing as the
!!$          ! input radiances, so we only need to find the index at which they are
!!$          ! essentially the same wavelength
!!$          if (abs(ILS_delta_min - wl_input(i)) < (0.1d0 * ILS_wl_spacing)) then
!!$             idx_hires_ILS_min = i
!!$             exit
!!$          end if
!!$       end do


       if (idx_hires_ILS_min == -1) then
          write(*,*) "idx_hires_ILS_min is -1"
          success = .false.
          return
       end if

      

       if (idx_hires_ILS_max > size(input)) then
          success = .false.
          return
       end if

       output(idx_pix) = dot_product(input(idx_hires_ILS_min:idx_hires_ILS_max), &
            kernels(kernel_idx_min:kernel_idx_max, idx_pix)) &
            / sum(kernels(kernel_idx_min:kernel_idx_max, idx_pix))

    end do

    success = .true.

  end subroutine oco_type_convolution

    subroutine linear_upsample(x_hires, x_lowres, y_lowres, output)

        implicit none
        double precision, intent(in) :: x_hires(:), x_lowres(:), y_lowres(:)
        double precision, intent(inout) :: output(:)

        integer :: idx_closest, idx_left, idx_right
        integer :: i,j
        double precision :: frac_value

        ! Output and hires arrays should be the same size!
        if ((size(output) /= size(x_hires))) then
            call logger%fatal("linear_upsample", "Array size mismatch!")
            stop 1
        end if

        do i=1, size(x_hires)

            ! Between which two lowres values is our hires value?
            idx_closest = minloc(abs(x_lowres - x_hires(i)), dim=1) !find_closest_index_DP(x_lowres, x_hires(i))
            if (x_lowres(idx_closest) < x_hires(i)) then
                idx_right = idx_closest + 1
                idx_left = idx_closest
            else
                idx_left = idx_closest - 1
                idx_right = idx_closest
            end if

            ! Boundary cases:
            if (idx_left == 0) then
                idx_left = 1
                idx_right = 2
            end if
            if (idx_right > size(x_lowres)) then
                idx_right = size(x_lowres)
                idx_left = idx_right - 1
            end if

            frac_value = (x_hires(i) - x_lowres(idx_left)) / (x_lowres(idx_right) - x_lowres(idx_left))
            output(i) = (1.0d0 - frac_value) * y_lowres(idx_left) + frac_value * y_lowres(idx_right)

        end do

    end subroutine


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

    end function

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


end module
