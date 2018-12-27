module solar_model

    use logger_mod, only: logger => master_logger

    public read_toon_spectrum, calculate_solar_distance_and_rv

contains

    subroutine read_toon_spectrum(filename, solar_spectrum)

        implicit none
        character(len=*), intent(in) :: filename
        double precision, dimension(:,:), allocatable, intent(inout) :: solar_spectrum
        !!!!!
        character(len=*), parameter :: fname = "read_toon_spectrum"
        integer :: funit, line_count, i, io_status
        logical :: exist_file
        double precision :: dummy

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

    end subroutine


    subroutine calculate_solar_distance_and_rv(dayofyear, distance, rv)

        ! This function returns the distance in AU and relative velocity in km/s
        ! between Earth and the fixed Sun. It's a simple polynomial fit to some
        ! ephemeris data that was constructed using sun_calculate.py. We hope this
        ! is exact enough for the purpose of correcting for Solar doppler shift,
        ! as this is somewhat faster than using a full ephemeris module.

        implicit none
        integer, intent(in) :: dayofyear
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


    end subroutine

end module
