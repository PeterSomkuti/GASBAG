!
! Main GeoCARB SIF Retrieval program. This is where it all starts. This main
! program handles the command-line input via FLAP as well as setting up the
! configuration through the text file via FINER.
!

program GeoCARBSIF

    ! Pick up the logging module for status messages
    use logger_mod, only: logger_init
    use Startup, only: initialize_config

    implicit none

    ! Before we can read the contents of the logfile, cast it all into /dev/null
    ! logger_init REQUIRES a logfile as an argument
    call logger_init('/dev/null')
    ! Initialize the whole thing by reading the configuration file
    call initialize_config()

end program
