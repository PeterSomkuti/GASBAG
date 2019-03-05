module Rayleigh_mod

  use math_utils_mod
  public calculate_rayleigh_tau


contains

  subroutine calculate_rayleigh_tau(wl, p_levels, ray_tau)

    implicit none
    double precision, intent(in) :: wl(:), p_levels(:)
    double precision, intent(inout) :: ray_tau(:,:)
    double precision, allocatable, save :: ray_CS(:)
    integer :: i, j

    ray_tau(:,:) = 0.0d0

    if (allocated(ray_CS)) deallocate(ray_CS)

    ! Calculate the rayleigh CS
    if (.not. allocated(ray_CS)) then
       call calculate_rayleigh_CS(wl, ray_CS)
    end if

    do i=1, size(p_levels) - 1
       do j=1, size(wl)
          ray_tau(j,i) = NA * ray_CS(j) * (p_levels(i+1) - p_levels(i)) &
               / (1.0d3 * DRY_AIR_MASS * 9.80665d0)
       end do
    end do

  end subroutine calculate_rayleigh_tau

  subroutine calculate_rayleigh_CS(wl, ray_CS)

    implicit none
    double precision, intent(in) :: wl(:)
    double precision, allocatable, intent(inout) :: ray_CS(:)

    double precision, parameter :: delta = 0.02790d0
    double precision, parameter :: a = 2.871d-4
    double precision, parameter :: b = 5.670d-3
    double precision, parameter :: cnst = 1.031d-24

    double precision :: wl2i, aniso, r, r2
    integer :: i

    allocate(ray_CS(size(wl)))

    aniso = (6.0d0 + 3.0d0 * delta) / (6.0d0 - 7.0d0 * delta)

    do i=1, size(wl)

       wl2i = 1.0d0 / (wl(i) * wl(i))
       r = 1.0d0 + a * (1.0d0 + b * wl2i)
       r2 = r * r + 2.0d0

       ray_CS(i) = cnst * wl2i * wl2i * aniso * ((r*r - 1.0d0) / r2) ** 2

    end do



  end subroutine calculate_rayleigh_CS



end module Rayleigh_mod
