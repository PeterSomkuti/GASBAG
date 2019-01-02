module radiance_mod

    use math_utils_mod

    public calculate_radiance

contains


    subroutine calculate_radiance(wavelengths, SZA, VZA, albedo, radiance)

        implicit none
        double precision, dimension(:), intent(in) :: wavelengths
        double precision, intent(in) :: SZA, VZA, albedo(:)
        double precision, dimension(:), intent(inout) :: radiance

        radiance(:) = albedo(:) / PI * cos(DEG2RAD * SZA)

    end subroutine



end module
