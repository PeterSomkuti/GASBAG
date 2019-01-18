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


  subroutine fft_convolution(input, kernel, wl_spacing, wl_output, output)

    implicit none
    double precision, intent(in) :: input(:), kernel(:), wl_spacing, wl_output(:)
    double precision, intent(inout) :: output(:)

    integer :: N_input, N_kernel, N_input_fft, N_kernel_fft
    integer(4) :: lensav, lenwrk
    integer :: fft_err
    double precision, allocatable :: wsave(:), work(:)
    double precision, allocatable :: input_fft(:), kernel_fft(:)


    ! FFT-type convolution using FFTPACK5.1, which is most efficient and
    ! quickest if our two arrays are powers of 2. We thus must zero pad
    ! to a) get the dimension of the array to 2^x (or at least some low
    ! number of prime factors), and b) to mitigate the wrap-around problem
    ! of the FFT method.

    ! Simple procedure:
    ! a) shift the ILS kernel such that the center is at index 1 and wraps
    !    around the array
    ! b) Perform forward FFT of both input and kernel
    ! c) Multiply the result and back-transform
    ! d) Sample the convolved output at the low-resultion
    !    wavelengths (usually given by the dispersion)

    N_input = size(input)
    N_kernel = size(kernel)

    ! Look up the rfft1i function - it recommends this size for he
    ! wsave array
    lensav = N_input + int(log(dble(N_input)) / log(2.0d0)) + 4

    allocate(work(N_input))
    allocate(wsave(lensav))

    ! Initialise the FFT solver
    call rfft1i(N_input, wsave, lensav, fft_err)

    ! Transform the input (radiances)
    N_input_fft = N_input
    allocate(input_fft(N_input_fft))
    input_fft(1:N_input) = input(:)
    call rfft1f(N_input, 1, input_fft, N_input_fft, lensav, work, fft_err)

    ! Transform the ILS kernel
    N_kernel_fft = N_kernel
    allocate(kernel_fft(N_kernel_fft))
    kernel_fft(1:N_kernel) = kernel(:)
    call rfft1f(N_kernel, 1, kernel_fft, N_kernel_fft, lensav, work, fft_err)

    read(*,*)


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
    integer :: N_pix, N_ils_pix
    integer :: idx_pix, idx_hires_closest, idx_hires_ILS_min, idx_hires_ILS_max
    double precision :: ILS_delta_min, ILS_delta_max
    double precision :: wl_diff, wl_diff_old

    double precision, allocatable :: ILS_upsampled(:), input_upsampled(:)
    double precision :: time_start, time_stop
    integer :: i, funit

    N_pix = size(wl_output)
    N_ils_pix = size(wl_kernels, 1)

    if (N_pix /= size(wl_kernels, 2)) then
       call logger%fatal(fname, "wl_kernels or wl_output have incompatible sizes.")
       stop 1
    end if

    if (N_pix /= size(kernels, 2)) then
       call logger%fatal(fname, "kernels or wl_output have incompatible sizes.")
       stop 1
    end if

    ! Main loop over all (instrument) output pixels
    do idx_pix=1, N_pix

       ! Note the ILS boundary in wavelength space
       ILS_delta_min = wl_output(idx_pix) + wl_kernels(1, idx_pix)
       ILS_delta_max = wl_output(idx_pix) + wl_kernels(N_ils_pix, idx_pix)

       ! Find out, which index of the high-res input is closest to the
       ! detector pixel, and also which is the center pixel
       !idx_hires_closest = -1
       idx_hires_ILS_min = -1
       idx_hires_ILS_max = -1


       !idx_hires_ILS_min = minloc(abs(wl_input - ILS_delta_min), dim=1)
       do i=1, size(wl_input)-1
          if ((ILS_delta_min >= wl_input(i)) .and. (ILS_delta_min <= wl_input(i+1))) then
             idx_hires_ILS_min = i
             exit
          end if
       end do

       if (idx_hires_ILS_min == -1) then
          success = .false.
          return
       end if

       !idx_hires_closest = find_closest_index_DP(wl_input, wl_output(idx_pix))
       !idx_hires_ILS_min = find_closest_index_DP(wl_input, ILS_delta_min)
       idx_hires_ILS_max = idx_hires_ILS_min + N_ils_pix - 1

       if (idx_hires_ILS_max > size(input)) then
          success = .false.
          return
       end if

       !if (idx_hires_ILS_max > size(input)) idx_hires_ILS_max = size(input)
       !call find_closest_index_DP(wl_input, ILS_delta_max, idx_hires_ILS_max)

       ! Now we have to either interpolate the ILS data onto the high-res grid,
       ! or, if the ILS data is higher resolution, do the other way round.
       !if (N_ils_pix > (idx_hires_ILS_max - idx_hires_ILS_min + 1)) then
       ! ILS is higher resolution than hires radiances
       !    write(*,*) N_ils_pix
       !    write(*,*) idx_hires_ILS_min, idx_hires_ILS_max, idx_hires_ILS_max - idx_hires_ILS_min
       !    call logger%fatal(fname, "Oops - not implemented yet! Call Peter.")
       !    stop 1
       !else
       ! ILS is lower resolution than hires radiances
       !allocate(ILS_upsampled(idx_hires_ILS_max - idx_hires_ILS_min + 1))

       ! This is a fairly expensive call! 0.3ms - but we do it 1016 times..
       ! Upsample the ILS to the hires wavelength grid
       !call linear_upsample(wl_input(idx_hires_ILS_min:idx_hires_ILS_max), &
       !                     wl_output(idx_pix) + wl_kernels(:, idx_pix), &
       !                     kernels(:, idx_pix), &
       !                     ILS_upsampled)
       ! And do the 'convolution'
       !write(*,*) '--------------'
       !write(*,*) shape(kernels)
       !write(*,*) idx_hires_ILS_min, idx_hires_ILS_max, idx_hires_ILS_max-idx_hires_ILS_min

       if ((idx_hires_ILS_min < 1) .or. (idx_hires_ILS_max > size(input))) then
          write(*,*) "dimension issue in oco_convolution routine"
          write(*,*) ILS_delta_min, ILS_delta_max
          write(*,*) idx_hires_ILS_min, idx_hires_ILS_max, size(input)
          write(*,*) wl_input(1), wl_input(size(wl_input)), ILS_delta_min, wl_output(idx_pix)
          !stop 1
       end if


       output(idx_pix) = dot_product(input(idx_hires_ILS_min:idx_hires_ILS_max), kernels(:, idx_pix)) &
            / sum(kernels(:, idx_pix))

       success = .true.

       !deallocate(ILS_upsampled)
       !end if
    end do

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

    pure function find_closest_index_DP(array, value) result(index)

        implicit none
        double precision, intent(in) :: array(:), value
        integer :: index

        double precision :: old_diff, diff
        integer :: i

        old_diff = abs(array(1) - value) + 1d3
        index = -1
        do i=1, size(array)
            diff = abs(array(i) - value)
            if (diff < old_diff) then
                old_diff = diff
                index = i
            end if
        end do

    end function


    pure function single_convolve(input, kernel)

        double precision, intent(in) :: input(:), kernel(:)
        double precision :: single_convolve
        integer :: i

        !! Direct convolution for a single pixel, so input and kernel need
        !! to be the same size.

        if (size(input) /= size(kernel)) then
            single_convolve = -9999.99
            return
        end if

        single_convolve = 0.0d0
        do i=1, size(input)
            single_convolve = single_convolve + input(i) * kernel(i)
        end do

        single_convolve = single_convolve / sum(kernel)

    end function


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
