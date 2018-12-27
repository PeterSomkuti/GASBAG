module solar_model

    use logger_mod, only: logger => master_logger

    public read_toon_spectrum

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

        ! Now turn them into wavelengths (microns) and re-order the array
        solar_spectrum(:,:) = solar_spectrum(line_count:1:-1,:)
        solar_spectrum(:,1) = 1e4 / solar_spectrum(:,1)

    end subroutine

end module
