module startup

    implicit none

    public :: initialize_config
    private :: populate_global_state

contains
    subroutine initialize_config
        ! This routine attempts to read the command-line argument that should
        ! contain the configuration file path, and then tries to parse it. If
        ! successful, it then initialises a logger object. The settings of the
        ! logger object are stored in the logger module, so it can be loaded
        ! from every other module/subroutine with the correct settings.

        use iso_fortran_env
        ! Pick up the logging module for status messages
        use logger_mod, only: logger => master_logger
        ! Also load the config file parser module "Finer"
        use finer, only: file_ini
        ! And the CLI argument parser module
        use flap, only : command_line_interface
        use stringifor

        use version

        implicit none

        type(command_line_interface) :: cli ! Command line interface handler
        type(file_ini) :: fini ! Config file interface handler
        integer :: error ! Error variable

        character(len=*), parameter :: fname = 'initialize_config'
        type(string) :: version_string

        character(999) :: config_file ! Path to the configuration file
        logical :: config_file_exists ! Does the file exist?
        integer :: config_file_size ! What size is the file?
        logical :: config_file_OK ! Is the config OK (resonable)?


        version_string = " -- Branch: " // git_branch // ", Hash: " // git_commit_hash

        ! Initialize the CLI interface - we only really need that one option,
        ! which is the location to the configuration file. Everything else
        ! should be contained within that text file.
        call cli%init(description = 'GeoCARBSIF Retrieval Algorithm', &
                      authors='Peter Somkuti (CSU/CIRA)', &
                      version=version_string%chars())

        call cli%add(switch='--config', &
                     switch_ab='-c', &
                     help='Path to the configuration file', &
                     required=.true., &
                     act='store', &
                     error=error)

        ! Check if setting up the command line arguments worked or not.
        if (error /= 0) then
            write(*, '(A)') "Could not set up command line arguments."
            stop
        end if

        ! Check if grabbing the command line argument worked or not.
        call cli%get(switch='-c', val=config_file, error=error)
        if (error /= 0) then
            write(*, '(A)') "Could not parse command line arguments."
            stop
        end if

        call logger%info(fname, &
                         "Attempting to read configuration file at: " &
                         // trim(config_file))

        ! Check if file exists, and how large it is
        inquire(file=config_file, exist=config_file_exists, &
                size=config_file_size)

        ! End program if the file does not exist
        if (.not. config_file_exists) then
            call logger%fatal(fname, "File: " // trim(config_file) // " does not exist.")
            stop
        end if

        ! End program if the file is empty
        if (config_file_size == 0) then
            call logger%fatal(fname, "File: " // trim(config_file) // " seems to be empty.")
            stop
        end if

        ! Now try to read the config file:
        call fini%load(filename=trim(config_file), error=error)

        ! If FINER can read this fine, we are good to go
        if (error /= 0) then
            call logger%fatal(fname, "Could not parse configuration file.")
            stop
        else
            call logger%info(fname, "Parsing of configuration file successful!")
        end if

        ! Check if we have any invalid sections, or invalid keywords within
        ! those sections. Rather than allowing all keywords, we stop the
        ! program immediately if something does not match up at this stage.
        call check_config(fini, config_file_OK)

        ! If the check failed, stop the program.
        if (.not. config_file_OK) then
            call logger%fatal(fname, "Config file " // trim(config_file) // &
                              " has not passed the check. Aborting.")
            stop
        else
            call logger%info(fname, "Config file seems OK!")
        end if

        ! Now push all the configuration file contents into the global
        ! configuration structure.
        call populate_global_state(fini)

    end subroutine

    subroutine check_config(fini, config_file_OK)
        ! This checks only whether the configuration file has any unexpected
        ! items in it. The downside - we MUST update the Kewords module anytime
        ! there is an update in the code itself
        use finer, only: file_ini ! CLI argument parser module
        ! Pick up the logging module for status messages
        use logger_mod, only: logger => master_logger
        use stringifor ! For maniupulating strings
        use Keywords

        implicit none
        ! DUMMY ARGUMENTS
        type(file_ini), intent(in out) :: fini ! Config file interface handler
        logical, intent(out) :: config_file_OK ! Is the config OK?

        character(len=*), parameter :: fname = 'Check_Config'

        character(len=:), allocatable :: sections(:) ! All sections in the file
        type(string) :: tmp_str

        logical :: found_section ! Have we found the section yet?
        integer :: i, j ! Loop variables

        ! Initialise valid_sections from the keywords module here
        call initialize_valid_sections()

        ! Assume the config file is OK for now
        config_file_OK = .true.

        ! First, check if each section is a valid one
        call fini%get_sections_list(sections)

        do i=1, size(sections)
            ! Loop through every section item and check it against every known
            ! and valid keyword. If found, set found_section to true.
            found_section = .false.
            tmp_str = sections(i)
            do j=1, size(valid_sections)
                if (trim(tmp_str%lower()) == valid_sections(j)%chars()) then
                    found_section = .true.
                    call logger%info(fname, "Found section [" // trim(tmp_str%chars()) // "]")
                end if
            end do

            if (.not. found_section) then
                ! Oh boy, this is not a valid section
                call logger%error(fname, 'Sorry, section keyword [' // &
                                      tmp_str%chars() // &
                                      '] is not in the list of known sections.')
                config_file_OK = .false.
            endif
        end do



    end subroutine

    subroutine populate_global_state(fini)

        use finer, only: file_ini ! CLI argument parser module
        use stringifor ! For maniupulating strings

        implicit none

        type(file_ini), intent(in out) :: fini ! Config file interface handler


    end subroutine


end module
