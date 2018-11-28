!
! Main GeoCARB SIF Retrieval program. This is where it all starts. This main
! program handles the command-line input via FLAP as well as setting up the
! configuration through the text file via FINER.
!
program GeoCARBSIF

    use startup, only: initialize_config
    use version, only: git_branch, git_commit_hash, git_rev_no
    use control, only: populate_control_structure
    use finer, only: file_ini

    use iso_fortran_env

    implicit none


    type(file_ini) :: fini ! The config file structure

    ! Greet the user and display information about the build itself.

    write(*,'(A)') "------------------------------"
    write(*,'(A)') "╔═╗┌─┐┌─┐╔═╗╔═╗╦═╗╔╗   ╔═╗╦╔═╗"
    write(*,'(A)') "║ ╦├┤ │ │║  ╠═╣╠╦╝╠╩╗  ╚═╗║╠╣ "
    write(*,'(A)') "╚═╝└─┘└─┘╚═╝╩ ╩╩╚═╚═╝  ╚═╝╩╚  "
    write(*,'(A)') "------------------------------"

    write(*,'(A)') "Version [" // git_branch // " " // git_commit_hash // &
                   " #" // git_rev_no // "]"

    ! Initialize the whole thing by reading the configuration file
    call initialize_config(fini)

    ! Initialize the program control structure with the settings from the
    ! config file.
    call populate_control_structure(fini)

    ! Check the supplied L1B file for sanity - are all variables present that
    ! we need? Do they have the right shapes?

    ! Go and perform the retrieval process

end program
