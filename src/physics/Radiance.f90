module radiance_mod

    use math_utils_mod

    public calculate_radiance

contains


    subroutine calculate_radiance(wavelengths, SZA, VZA, albedo, tau, radiance)

        implicit none
        double precision, dimension(:), intent(in) :: wavelengths
        double precision, intent(in) :: SZA, VZA, albedo(:)
        double precision, allocatable, intent(in) :: tau(:,:,:)
        double precision, dimension(:), intent(inout) :: radiance
        double precision :: mu0, mu

        double precision :: tau_flat(size(wavelengths))
        integer :: i, funit
        mu0 = cos(DEG2RAD * SZA)
        mu = cos(DEG2RAD * VZA)

        radiance(:) = albedo(:) / PI * mu0

        if (allocated(tau)) then
           ! Have gas absorbers in the atmosphere?
           write(*,*) "We are calculating tau.."
           tau_flat = sum(SUM(tau, dim=2), dim=2)
           open(newunit=funit, file='gas_test2.txt')
           do i=1, size(wavelengths)
              write(funit,*) wavelengths(i), tau_flat(i)
           end do
           close(funit)
           

           radiance(:) = radiance(:) * exp(-1.0d0/mu0 * sum(SUM(tau, dim=2), dim=2))
           radiance(:) = radiance(:) * exp(-1.0d0/mu * sum(SUM(tau, dim=2), dim=2))
        endif

    end subroutine



end module
