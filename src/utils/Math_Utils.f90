module math_utils

    use logger_mod, only: logger => master_logger

    public :: combsort, percentile

contains


    subroutine invert_matrix(mat_in, mat_out)

        double precision, dimension(:,:), intent(in) :: mat_in
        double precision, dimension(:,:), intent(out) :: mat_out

        double precision, dimension(size(mat_in,1)) :: work  ! work array for LAPACK
        integer, dimension(size(mat_in,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        !external DGETRF
        !external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        mat_out(:,:) = mat_in(:,:)
        n = size(mat_in,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, mat_out, n, ipiv, info)

        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, mat_out, n, ipiv, work, n, info)

        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if

    end subroutine


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
