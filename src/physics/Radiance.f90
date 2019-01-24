module radiance_mod

    use math_utils_mod

    public calculate_radiance

contains


  subroutine calculate_radiance(wavelengths, SZA, VZA, albedo, tau, &
       radiance)

        implicit none
        double precision, intent(in) :: wavelengths(:)
        double precision, intent(in) :: SZA, VZA, albedo(:)
        double precision, allocatable, intent(in) :: tau(:,:,:)
        double precision, intent(inout) :: radiance(:)
        double precision :: mu0, mu

        double precision :: tau_flat(size(wavelengths))
        integer :: i, funit
        mu0 = cos(DEG2RAD * SZA)
        mu = cos(DEG2RAD * VZA)

        radiance(:) = albedo(:) / PI * mu0

        if (allocated(tau)) then
           ! Have gas absorbers in the atmosphere?
           radiance(:) = radiance(:) * exp(-1.0d0/mu0 * sum(sum(tau, dim=3), dim=2))
           radiance(:) = radiance(:) * exp(-1.0d0/mu * sum(sum(tau, dim=3), dim=2))
        endif

    end subroutine



end module
