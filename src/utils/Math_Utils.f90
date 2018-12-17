module math_utils

  use logger_mod, only: logger => master_logger

  public :: combsort, percentile

contains

  function percentile(x, perc)

    !! Calculates the percentile "perc" of a double precision array "x"

    implicit none
    double precision :: percentile

    double precision, intent(in) :: x(:)
    double precision, intent(in) :: perc
    double precision, allocatable :: x_sort(:)
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
