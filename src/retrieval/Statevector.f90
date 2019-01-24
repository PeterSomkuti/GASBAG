module statevector_mod

  use logger_mod, only: logger => master_logger
  use control_mod, only: MCS, MAX_GASES
  use stringifor, only: string

    type statevector
        ! Number of state vector elements per type
        integer :: num_albedo, num_sif, num_dispersion, num_psurf, num_gas
        integer, dimension(:), allocatable :: idx_albedo, idx_sif, &
             idx_dispersion, idx_psurf, idx_gas, gas_idx_lookup
        ! State vector (current), state vector a priori
        double precision, dimension(:), allocatable :: svsv, svap, sver
        double precision, dimension(:,:), allocatable :: sv_ap_cov, sv_post_cov
    end type

    public initialize_statevector, parse_and_initialize_SV, clear_SV

contains



  subroutine clear_SV(SV)
    implicit none
    type(statevector), intent(inout) :: SV

    if (allocated(SV%idx_albedo)) deallocate(SV%idx_albedo)
    if (allocated(SV%idx_sif)) deallocate(SV%idx_sif)
    if (allocated(SV%idx_dispersion)) deallocate(SV%idx_dispersion)
    if (allocated(SV%idx_psurf)) deallocate(SV%idx_psurf)

    if (allocated(SV%svap)) deallocate(SV%svap)
    if (allocated(SV%svsv)) deallocate(SV%svsv)
    if (allocated(SV%sver)) deallocate(SV%sver)

    if (allocated(SV%sv_ap_cov)) deallocate(SV%sv_ap_cov)
    if (allocated(SV%sv_post_cov)) deallocate(SV%sv_post_cov)

    SV%num_albedo = -1
    SV%num_sif = -1
    SV%num_dispersion = -1
    SV%num_psurf = -1
    SV%num_gas = -1

  end subroutine clear_SV


  subroutine parse_and_initialize_SV(i_win, SV)

    implicit none
    integer, intent(in) :: i_win
    type(statevector), intent(inout) :: SV

    character(len=*), parameter :: fname = "parse_and_initialize_statevector"
    type(string), allocatable :: split_string(:)

    ! Which are the SV elements known by the code, and which one's are used?
    type(string) :: known_SV(99)
    logical :: is_gas_SV(99)

    integer :: num_SV, i, j, k, gas_index
    logical :: SV_found

    integer :: num_albedo_parameters, num_dispersion_parameters, &
         num_sif_parameters, num_psurf_parameters

    known_SV(:) = ""
    is_gas_SV(:) = .false.

    ! These are the known state vector elements - only these will be properly
    ! parsed. The order in which they are defined is not significant.
    known_SV(1) = "albedo"
    known_SV(2) = "dispersion"
    known_SV(3) = "sif"
    known_SV(4) = "psurf"

    ! Add gases as 'known' state vector elements. CAUTION! There is
    ! obviously a danger if someone decides to name their gas "psurf".
    ! Maybe I should add a list of protected names that can't be used.
    do i=1, MCS%window(i_win)%num_gases
       known_SV(4+i) = MCS%window(i_win)%gases(i)
       ! This flags the known SV as one of a gas type
       is_gas_SV(4+i) = .true.
    end do


    ! Split the state vector string from the config file
    call MCS%window(i_win)%SV_string%split(tokens=split_string, sep=' ')

    do i=1, size(split_string)
       SV_found = .false.

       ! Loop through the known SV elements, and exit if any requested
       ! one is not found.
       do j=1, size(known_SV)
          if (known_SV(j) /= "") then
             if (split_string(i)%lower() == known_SV(j)%lower()) then
                SV_found = .true.
                exit
             end if
          end if
       end do

       if (.not. SV_found) then
          call logger%fatal(fname, "This SV element was not recognized: " // &
               trim(split_string(i)%chars()))
          stop 1
       else
          call logger%info(fname, "SV element recognized: " // &
               trim(split_string(i)%chars()))
       end if
    end do

    ! Now that we have made sure that only SV elements known by the code are being
    ! supplied, we can loop through them, and depending on the type also check
    ! if necessary parameters are available. Example: if we want to retrieve
    ! dispersion, we also require the order and perturbation values.

    num_albedo_parameters = 0
    num_dispersion_parameters = 0
    num_sif_parameters = 0
    num_psurf_parameters = 0

    MCS%window(i_win)%gas_retrieved(:) = .false.

    do i=1, size(split_string)

       ! We are retrieving ALBEDO!
       if (split_string(i)%lower() == "albedo") then

          ! Albedo order needs to be >= 0
          if (MCS%window(i_win)%albedo_order < 0) then
             call logger%fatal(fname, "We are retrieving albedo, but the albedo order " &
                  // "needs to be >= 0. Check if you've supplied a sensible value (or at all).")
             stop 1
          else
             ! Albedo order 0 means simple factor, order 1 is with slope etc.
             num_albedo_parameters = MCS%window(i_win)%albedo_order + 1
          end if

       end if


       ! We are retrieving DISPERSION!
       if (split_string(i)%lower() == "dispersion") then
          ! Dispersion order needs to be > 0
          if (MCS%window(i_win)%dispersion_order <= 0) then
             call logger%fatal(fname, "We are retrieving dispersion, but the dispersion order " &
                  // "needs to be > 0. Check if you've supplied a sensible value (or at all).")
             stop 1
          else
             ! Dispersion order 1 means shift, 2 is stretch etc., this is not quite
             ! consistent with the albedo order notation, but whatever.
             num_dispersion_parameters = MCS%window(i_win)%dispersion_order

             ! We MUST have at least the same number of dispersion perturbation
             ! elements.
             if (.not. allocated(MCS%window(i_win)%dispersion_pert)) then
                call logger%fatal(fname, "Dispersion perturbation not in config file!")
                stop 1
             end if

             if (num_dispersion_parameters > size(MCS%window(i_win)%dispersion_pert)) then
                call logger%fatal(fname, "Not enough disperison perturbation values!")
                stop 1
             end if
          end if
       end if

       ! We are retrieving SIF!
       if (split_string(i)%lower() == "sif") then
          num_sif_parameters = 1
       end if

       ! We are retrieving surface pressure!
       if (split_string(i)%lower() == "psurf") then
          ! Surface pressure only makes sense if we have gases defined
          ! in our current microwindow. 

          if (MCS%window(i_win)%num_gases == 0) then
             call logger%fatal(fname, "Sorry, you must have at least one gas in the window " &
                  // "to retrieve surface pressure!")
             stop 1
          else
             num_psurf_parameters = 1
          end if
       end if


       ! Now check for the gases
       do j=1, size(known_SV)
          ! These are case-sensitive (not sure why?)
          ! Is this SV element a gas-type state vector?
          if ((split_string(i) == known_SV(j)) .and. (is_gas_SV(j))) then
             ! Yes, it is - now we look which particular gas this is, and set the
             ! retrieved-flag accordingly.
             do k=1, MCS%window(i_win)%num_gases
                gas_index = MCS%window(i_win)%gas_index(k)
                ! If this gas is not used in this window, skip!
                if (.not. MCS%gas(gas_index)%used) cycle
                ! Otherwise, loop through the gases we have in the window
                ! and set them as 'retrieved'
                if (MCS%window(i_win)%gases(k) == split_string(i)) then
                   MCS%window(i_win)%gas_retrieved(k) = .true.
                   call logger%info(fname, "We are retrieving gas: " &
                        // MCS%window(i_win)%gases(k)%chars())
                end if
             end do
          end if
       end do

    end do

    ! Once we have all the data, initialize the state vector with the appropriate
    ! number of elements for each type.

    call initialize_statevector( &
         i_win, &
         SV, &
         num_albedo_parameters, &
         num_sif_parameters, &
         num_dispersion_parameters, &
         num_psurf_parameters)

  end subroutine parse_and_initialize_SV




  subroutine initialize_statevector(i_win, sv, count_albedo, count_sif, &
       count_dispersion, count_psurf)

    implicit none
    integer, intent(in) :: i_win
    type(statevector), intent(inout) :: sv
    integer, intent(in) :: count_albedo, count_sif, &
         count_dispersion, count_psurf

    integer :: count_gas
    character(len=*), parameter :: fname = "initialize_statevector"
    character(len=999) :: tmp_str
    integer :: sv_count
    integer :: i, cnt

    ! Set The Number of Statevector parameters
    sv%num_albedo = count_albedo
    sv%num_sif = count_sif
    sv%num_dispersion = count_dispersion
    sv%num_psurf = count_psurf
    sv%num_gas = 0

    sv_count = 0
    ! And determine the position (indices) of the state vecto elements within
    ! the state vector

    ! Albedo: arbitrary number of parameters allowed
    if (sv%num_albedo > 0) then

       write(tmp_str, '(A, G0.1)') "Number of albedo SV elements: ", sv%num_albedo
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_albedo(sv%num_albedo))
       do i=1, sv%num_albedo
          sv_count = sv_count + 1
          sv%idx_albedo(i) = sv_count
       end do
    else
       allocate(sv%idx_albedo(1))
       sv%idx_albedo(1) = -1
    end if

    ! SIF: we can do only two things here, SIF magnitude and slope
    if (sv%num_sif > 2) then
       call logger%fatal(fname, "Sorry! Only up to 2 SIF SV elements supported!")
       stop 1
    end if

    if (sv%num_sif > 0) then

       write(tmp_str, '(A, G0.1)') "Number of SIF SV elements: ", sv%num_sif
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_sif(sv%num_sif))
       do i=1, sv%num_sif
          sv_count = sv_count + 1
          sv%idx_sif(i) = sv_count
       end do
    else
       allocate(sv%idx_sif(1))
       sv%idx_sif(1) = -1
    end if

    ! Dispersion: we do arbitrary coefficients here, and we'll have to check
    ! in the retrieval subroutines whether the number passed into this function
    ! actually makes sense.

    if (sv%num_dispersion > 0) then

       write(tmp_str, '(A, G0.1)') "Number of dispersion SV elements: ", sv%num_dispersion
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_dispersion(sv%num_dispersion))
       do i=1, sv%num_dispersion
          sv_count = sv_count + 1
          sv%idx_dispersion(i) = sv_count
       end do
    else
       allocate(sv%idx_dispersion(1))
       sv%idx_dispersion(1) = -1
    end if

    ! Surface pressure
    if (sv%num_psurf == 1) then
       write(tmp_str, '(A)') "Number of surface pressure SV elements: 1"
       call logger%info(fname, trim(tmp_str))
       allocate(sv%idx_psurf(1))
       sv_count = sv_count + 1
       sv%idx_psurf = sv_count
    end if


    ! Gases
    ! In order to find out which gases are retrieved, we acces the MCS
    ! rather than have it passed into this function. For the time being, we are
    ! retrieving a scale factor?


    count_gas = count(MCS%window(i_win)%gas_retrieved(:) .eqv. .true.)
    sv%num_gas = count_gas
    
    if (count_gas > 0) then
       allocate(sv%idx_gas(count_gas))
       allocate(sv%gas_idx_lookup(count_gas))

       write(*,*) (MCS%window(i_win)%gases(i)%chars(), i=1, MCS%window(i_win)%num_gases)
       write(*,*) MCS%window(i_win)%gas_retrieved(:)
       write(*,*) count_gas

       
       sv%idx_gas(:) = -1
       sv%gas_idx_lookup(:) = -1

       if (count_gas > 0) then
          write(tmp_str, '(A,G0.1)') "Number of gas SV elements: ", count_gas
          call logger%info(fname, trim(tmp_str))

          cnt = 1
          do i=1, MCS%window(i_win)%num_gases
             if (MCS%window(i_win)%gas_retrieved(i)) then
                sv_count = sv_count + 1
                sv%idx_gas(cnt) = sv_count
                ! This is vitally needed to make sure we can match the
                ! gas state vector element with the position in the
                ! gas optical depth arrays.
                sv%gas_idx_lookup(cnt) = i
                cnt = cnt + 1
             end if
          end do
       end if
    end if

    write(tmp_str, '(A, G0.1, A)') "We have ", sv_count, " SV elements."
    call logger%debug(fname, trim(tmp_str))

    allocate(sv%svap(sv_count))
    allocate(sv%svsv(sv_count))
    allocate(sv%sver(sv_count))

    sv%svap(:) = -9999.99
    sv%svsv(:) = -9999.99
    sv%sver(:) = -9999.99

    allocate(sv%sv_ap_cov(sv_count, sv_count))
    allocate(sv%sv_post_cov(sv_count, sv_count))

  end subroutine initialize_statevector

end module
