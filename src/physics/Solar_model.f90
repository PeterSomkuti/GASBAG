module solar_model_mod

    use logger_mod, only: logger => master_logger
    use math_utils_mod

    public read_toon_spectrum, calculate_solar_distance_and_rv

contains

    subroutine read_toon_spectrum(filename, solar_spectrum, wl_min, wl_max)

        implicit none
        character(len=*), intent(in) :: filename
        double precision, dimension(:,:), allocatable, intent(inout) :: solar_spectrum
        double precision, optional, intent(in) :: wl_min, wl_max
        !!!!!

        character(len=*), parameter :: fname = "read_toon_spectrum"
        integer :: funit, line_count, i, io_status
        logical :: exist_file
        double precision :: dummy
        double precision, allocatable :: tmp_array(:,:)
        integer :: idx_min, idx_max

        call logger%debug(fname, "Reading Toon solar spectrum at: " // trim(filename))

        inquire(file=filename, exist=exist_file)
        if (exist_file .eqv. .false.) then
            call logger%fatal(fname, "Toon spectrum file could not be found.")
            stop 1
        end if

        ! First count how many lines we have
        open(newunit=funit, file=filename, status='old', action='read')
        line_count = 0
        do
            read(funit, *, iostat=io_status) dummy
            if (io_status < 0) then
                exit
            end if
            line_count = line_count + 1
        end do
        close(funit)

        ! And then read it all in, properly
        allocate(solar_spectrum(line_count, 2))
        open(newunit=funit, file=filename, status='old', action='read')
        do i=1, line_count
            read(funit, *) solar_spectrum(i,1), solar_spectrum(i,2)
        end do
        close(funit)
        call logger%debug(fname, "Finished reading in Toon spectrum.")

        ! Now turn them into wavelengths (microns) and re-order the array, so that
        ! we have the spectrum in units of increasing wavelength
        solar_spectrum(:,:) = solar_spectrum(line_count:1:-1,:)
        solar_spectrum(:,1) = 1e4 / solar_spectrum(:,1)

        ! If the user wants to, we can cut off the fairly large solar spectrum
        ! array, depending on the wavelength limits supplied to the function

        if(present(wl_min) .and. present(wl_max)) then
            idx_min = -1
            idx_max = -1

            idx_min = minloc(abs(solar_spectrum(:,1) - wl_min), dim=1)
            idx_max = minloc(abs(solar_spectrum(:,1) - wl_max), dim=1)

            allocate(tmp_array(idx_max - idx_min + 1, 2))
            tmp_array(:,:) = solar_spectrum(idx_min:idx_max, :)

            deallocate(solar_spectrum)
            allocate(solar_spectrum(idx_max - idx_min + 1, 2))
            solar_spectrum(:,:) = tmp_array(:,:)
        end if


    end subroutine


    subroutine calculate_solar_distance_and_rv(dayofyear, distance, rv)

        ! This function returns the distance in km and relative velocity in km/s
        ! between Earth and the fixed Sun. It's a simple polynomial fit to some
        ! ephemeris data that was constructed using sun_calculate.py. We hope this
        ! is exact enough for the purpose of correcting for Solar doppler shift,
        ! as this is somewhat faster than using a full ephemeris module.

        ! The concept is shamelessly taken from Hartmut Bösch's solar module.

        implicit none
        double precision, intent(in) :: dayofyear
        double precision, intent(out) :: distance, rv

        ! Coefficients for the distance [AU] polynomial, derived using "sun_calculate.py"
        double precision, dimension(*), parameter :: poly_dist = &
        [ &
        0.983346446687037d0, &
        -1.8104806582301097d-05, &
        2.51272811877413d-06, &
        2.5878414474217157d-09, &
        -9.823924798139483d-11, &
        2.2221419917641396d-13, &
        1.63122084729368d-16, &
        -8.109432394654854d-19, &
        5.213684249628703d-22 &
        ]

        ! .. and the same for relative velocity [km/s]
        double precision, dimension(*), parameter :: poly_rv = &
        [ &
        -0.041756570228324644d0, &
        0.009976285299598773d0, &
        -3.595896055785801d-05, &
        2.059174396382475d-07, &
        -6.6266091284662725d-09, &
        4.8784827830332416d-11, &
        -1.579046809797818d-13, &
        2.5440629508031043d-16, &
        -1.6979380980265063d-19 &
        ]

        integer :: i


        distance = 0.0d0
        rv = 0.0d0

        do i=1, size(poly_dist)
            distance = distance + (dayofyear ** (i-1)) * poly_dist(i)
            rv = rv + (dayofyear ** (i-1)) * poly_rv(i)
        end do

        ! Convert AU to km
        distance = distance * 149597870.700d0

    end subroutine


    subroutine calculate_rel_velocity_earth_sun(lat, sza, saa, altitude, rv)
        implicit none
        double precision, intent(in) :: lat, sza, saa, altitude
        double precision, intent(out) :: rv

        double precision, parameter :: EARTH_ROT_SPD = 7.2921150D-5
        double precision, parameter :: CON = 6.73951496D-3
        double precision, parameter :: EARTH_RADIUS = 6378.137 ! [km]

        double precision :: gclat, gcrad

         gclat = ATAN(TAN(DEG2RAD*lat)/(1.d0+CON))
         gcrad = 1000.d0 * ((altitude / 1000.d0) + EARTH_RADIUS/SQRT(1.d0+CON*SIN(gclat)**2))

         rv = -EARTH_ROT_SPD * gcrad * SIN(DEG2RAD*sza) * COS(DEG2RAD*(saa-90.d0)) * COS(gclat)

     end subroutine

     subroutine calculate_solar_planck_function(temp, solar_distance, wavelength, irradiance)
         implicit none
         double precision, intent(in) :: temp ![K]
         double precision, intent(in) :: solar_distance ! [m]
         double precision, dimension(:), intent(in) :: wavelength ! in microns

         double precision, dimension(:), allocatable, intent(out) :: irradiance
         double precision, parameter :: SPEED_OF_LIGHT = 299792458d+0
         double precision, parameter :: PLANCK_CONST = 6.62607015d-34
         double precision, parameter :: BOLTZMANN_CONST = 1.38064852d-23
        double precision, parameter :: SOLAR_RADIUS = 695.660d6 ! [m] 

         double precision :: solar_solid_angle, solar_angular_radius

         double precision, dimension(:), allocatable :: denom

         allocate(irradiance(size(wavelength)))
         allocate(denom(size(wavelength)))

         ! Wavelength comes in as microns, so convert them to meters

         solar_angular_radius = atan(SOLAR_RADIUS / solar_distance)
         solar_solid_angle = PI * solar_angular_radius * solar_angular_radius

         irradiance(:) = 2d18 * SPEED_OF_LIGHT / ((wavelength(:))**4)
         denom(:) = exp(1d6 * PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * BOLTZMANN_CONST * temp)) - 1.0d0
         irradiance(:) = solar_solid_angle * irradiance(:) / denom(:)

         !irradiance(:) = 2 * PLANCK_CONST * (SPEED_OF_LIGHT**2) / ((wavelength(:) * 1d-6) ** 5)
         !denom(:) = exp(PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * 1d-6 * BOLTZMANN_CONST * temp)) - 1.0d0
         !irradiance(:) = irradiance(:) / denom(:)
         ! And convert to ph/s/m2/sr/um
         !irradiance(:) = SUN_APP_SIZE * 1d-6 * irradiance(:) / (PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * 1d-6))
     end subroutine


end module
