module radiance_mod

    use math_utils_mod

    public calculate_radiance

contains


  subroutine calculate_Beer_Lambert(wavelengths, mu0, mu, albedo, tau, &
       radiance)

        implicit none
        double precision, intent(in) :: wavelengths(:)
        double precision, intent(in) :: mu0, mu, albedo(:)
        double precision, allocatable, intent(in) :: tau(:)
        double precision, intent(inout) :: radiance(:)

        radiance(:) = albedo(:) / PI * mu0

        if (allocated(tau)) then
           ! Have gas absorbers in the atmosphere?
           radiance(:) = radiance(:) * exp(-1.0d0/mu0 * tau(:)) * exp(-1.0d0/mu * tau(:))
        endif

    end subroutine



end module
