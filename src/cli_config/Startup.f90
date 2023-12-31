!> @brief Startup module for initialization of the main control structure (MCS)
!> @file Startup.f90
!> @author Peter Somkuti
!>
!> @details Startup module. Routines in here are called during the startup-phase of the
!> program, where we also check if the config/inputs are sane. At this point, the
!> logger is not yet initialized, so the various messages do not have the usual
!> formatting as the messages later on.

module startup_mod

    implicit none

    public :: initialize_config
    private :: check_config

  contains
    !> This routine attempts to read the command-line argument that should
    !> contain the configuration file path, and then tries to parse it. If
    !> successful, it then initialises a logger object. The settings of the
    !> logger object are stored in the logger module, so it can be loaded
    !> from every other module/subroutine with the correct settings.
    !>
    !! @param fini FINER ini object
    subroutine initialize_config(fini)


      ! User modules
      use version_mod
      use file_utils_mod, only: fini_extract_CHAR, fini_extract_INT

      ! Third-party modules
      use logger_mod, only: logger_init, logger => master_logger
      use finer, only: file_ini
      use flap, only: command_line_interface
      use stringifor, only: string


      ! System modules
      use iso_fortran_env

      implicit none

      type(file_ini), intent(inout) :: fini ! Config file interface handler

      ! Local variables
      type(command_line_interface) :: cli ! Command line interface handler
      integer :: error ! Error variable

      character(len=*), parameter :: fname = 'initialize_config'
      type(string) :: version_string

      character(999) :: config_file ! Path to the configuration file
      logical :: config_file_exists ! Does the file exist?
      integer :: config_file_size ! What size is the file?
      logical :: config_file_OK ! Is the config OK (resonable)?

      character(999) :: logfile
      integer :: loglevel = 10

      character(999) :: tmp_str

      ! Create a version string that can be displayed at the command line
      ! if you ask it via -v or --version
      version_string = "[" // git_branch // " " // git_commit_hash //  &
           " #" // git_rev_no // "]"

      ! Initialize the CLI interface - we only really need that one option,
      ! which is the location to the configuration file. Everything else
      ! should be contained within that text file.
      call cli%init(description = 'Generic Algorithm for Single-Band Acquisition of Gases', &
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
         stop 1
      end if

      ! Check if grabbing the command line argument worked or not.
      call cli%get(switch='-c', val=config_file, error=error)
      if (error /= 0) then
         write(*, '(A)') "Could not parse command line arguments."
         stop 1
      end if

      write(*, '(A)') "Attempting to read configuration file at: " &
           // trim(config_file)

      ! Check if file exists, and how large it is
      inquire(file=config_file, exist=config_file_exists, &
           size=config_file_size)

      ! End program if the file does not exist
      if (.not. config_file_exists) then
         write(*, '(A)') "File: " // trim(config_file) // " does not exist."
         stop 1
      end if

      ! End program if the file is empty
      if (config_file_size == 0) then
         write(*, '(A)') "File: " // trim(config_file) // " seems to be empty."
         stop 1
      end if

      ! Now try to read the config file
      call fini%load(filename=trim(config_file), error=error)

      ! If FINER can read this fine, we are good to go
      if (error /= 0) then
         write(*, '(A)') "Could not parse configuration file."
         stop 1
      else
         write(*, '(A)') "Parsing of configuration file successful!"
      end if

      ! Check if we have any invalid sections, or invalid keywords within
      ! those sections. Rather than allowing all keywords, we stop the
      ! program immediately if something does not match up at this stage.
      call check_config(fini, config_file_OK)

      ! If the check failed, stop the program.
      if (.not. config_file_OK) then
         write(*, '(A)') "Config file " // trim(config_file) // &
              " has not passed the check. Aborting."
         stop 1
      else
         write(*, '(A)') "Config file seems OK! Moving on."
      end if


      ! Extract the logfile name. This is a required field.
      call fini_extract_CHAR(fini, 'logger', 'logfile', .true., logfile)

      ! Extract the loglevel (optional)
      call fini_extract_INT(fini, 'logger', 'loglevel', .false., loglevel)
      ! If not specified, just keep it at 10=debug, which shows essentially
      ! all log messages. This will be a lot, but hey, you wanted it!
      if (loglevel == -9999) loglevel = 10

      ! Initialize the logger entity with the path to the logfile and keep all
      ! loggers at the same log-level (stdout, stderr, logfile)
      call logger_init(trim(logfile), loglevel, loglevel, loglevel)

      ! From here on, NO NEED TO USE write(*,*), just use the logging commands
      ! to get debug,info etc. outputs. It'll create timestamps and everything
      ! nicely for you.

      write(tmp_str, '(A, I2.0)') "Logger initialized at [" // trim(logfile) // &
           "] with loglevel: ", loglevel
      call logger%info(fname, trim(tmp_str))

    end subroutine initialize_config

    !> Here we check whether the contents of the configuration file, as they
    !> are passed into "fini" (FINER object), are sensible and without contradictions.
    !>
    !! @param fini FINER ini object
    !! @param config_file_OK Is the configuration file deemed sane?
    subroutine check_config(fini, config_file_OK)
      ! This checks only whether the configuration file has any unexpected
      ! items in it. The downside - we MUST update the Kewords module anytime
      ! there is an update in the code itself. Note that this routine
      ! merely checks whether the config file structure makes sense. Any
      ! further checks about the contents of the config file items themselves
      ! are made down the line, when the respective item is actually used.

      use finer, only: file_ini ! CLI argument parser module
      use stringifor, only: string, trim ! For maniupulating strings
      use keywords_mod, only: initialize_valid_sections, valid_sections, valid_options

      implicit none

      type(file_ini), intent(in out) :: fini ! Config file interface handler
      logical, intent(out) :: config_file_OK ! Is the config OK?

      ! Local variables
      character(len=:), allocatable :: sections(:) ! All sections in the file
      character(len=:), allocatable :: option_pairs(:) ! Options in section
      type(string) :: tmp_section, tmp_option

      logical :: found_section ! Have we found the section yet?
      logical :: found_option ! .. or the option?
      ! The index used for valid_sections when a valid section has been found
      integer :: section_idx
      integer :: i, j, k ! Loop variables

      ! Initialise valid_sections from the keywords module here, so we can
      ! access them from outside.
      call initialize_valid_sections()

      ! Assume the config file is OK for now
      config_file_OK = .true.
      section_idx = 0

      ! First, check if each section is a valid one, and whether each option
      ! within that section is a valid option
      call fini%get_sections_list(sections)

      do i=1, size(sections)
         ! Loop through every section item in the config and check it against
         !  every known and valid keyword. If found, set found_section=true.
         found_section = .false.
         tmp_section = fini%section(i)

         do j=1, size(valid_sections)
            if (trim(tmp_section%lower()) == valid_sections(j)%chars()) then
               found_section = .true.
               write(*, '(A)') "Found section [" // &
                    trim(tmp_section%chars()) // "]"
               ! This section in the config file corresponds to section_idx
               ! in valid_sections, and valid_options(section_idx, :)
               section_idx = j
               exit
            end if
         end do

         if (.not. found_section) then
            ! Oh boy, this is not a valid section
            write(*, '(A)') 'Sorry, section keyword [' // &
                 tmp_section%chars() // &
                 '] is not in the list of known sections.'
            config_file_OK = .false.
            stop 1
         else
            ! But if successful, we now need to check if all section options
            ! are valid.
            do while(fini%loop(section_name=fini%section(i), &
                 option_pairs=option_pairs))

               found_option = .false.

               if (.not. allocated(option_pairs)) then
                  write(*, '(A)') "Could not read an option-value pair."
                  write(*, '(A)') "Make sure there are no strange characters, " &
                       // "especially equal-signs in the option values!"
                  stop 1
               end if

               ! option_pairs is (option name, whatever-is-right-of-equal-sign)

               tmp_option = option_pairs(1)

               do k=1, size(valid_options, 2)
                  ! Loop over the 'dictionary' of valid options, corresponding
                  ! to this particular section.

                  ! If that option is emtpy, just skip it
                  if (trim(valid_options(section_idx, k)) == "") then
                     cycle
                  end if

                  ! .. otherwise, compare it against the config entry
                  if (trim(tmp_option%lower()) == &
                       valid_options(section_idx, k)%chars()) then
                     ! Yes, we found a valid option!
                     write(*, '(A)') 'Option: "' &
                          // trim(option_pairs(1)) // '" found'
                     write(*, '(A)') trim(option_pairs(1)) // ' value(s): ' // option_pairs(2)

                     found_option = .true.
                  end if
               end do

               if (.not. found_option) then
                  ! This option is not in the list of recognized ones
                  write(*, '(A)') 'Sorry - option "' &
                       // option_pairs(1) &
                       // '" is not recognized.'
                  config_file_OK = .false.
                  stop 1
               end if

            end do

         endif
      end do ! end loop over sections

    end subroutine check_config

end module ! end module Startup
