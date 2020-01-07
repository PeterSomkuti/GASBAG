module statevector_mod

  use logger_mod, only: logger => master_logger
  use control_mod, only: MCS, MAX_GASES
  use stringifor, only: string

  type statevector
     ! Number of state vector elements per type
     integer :: num_albedo
     integer :: num_sif
     integer :: num_dispersion
     integer :: num_psurf
     integer :: num_gas
     integer :: num_solar_shift
     integer :: num_solar_stretch
     integer :: num_zlo
     integer :: num_temp
     integer :: num_ils_stretch
     integer :: num_aerosol_aod
     integer :: num_aerosol_height


     integer, dimension(:), allocatable :: idx_albedo
     integer, dimension(:), allocatable :: idx_sif
     integer, dimension(:), allocatable :: idx_dispersion
     integer, dimension(:), allocatable :: idx_psurf
     integer, dimension(:), allocatable :: idx_solar_shift
     integer, dimension(:), allocatable :: idx_solar_stretch
     integer, dimension(:), allocatable :: idx_zlo
     integer, dimension(:), allocatable :: idx_temp
     integer, dimension(:), allocatable :: idx_ils_stretch
     integer, dimension(:), allocatable :: idx_aerosol_aod
     integer, dimension(:), allocatable :: idx_aerosol_height

     ! Lookup arrays that tell us which aerosol/gas this statevector
     ! element belongs to.
     ! E.g. aerosol_aod_idx_lookup(2) = 1
     ! means that the second AOD SV element (at position idx_aerosol_aod)
     ! refers to the aerosol found in MCS%window(i_win)%aerosol(1)
     ! .. and so on for gases.

     integer, allocatable :: aerosol_aod_idx_lookup(:)
     integer, allocatable :: aerosol_height_idx_lookup(:)
     integer, allocatable :: gas_idx_lookup(:)

     integer, allocatable :: idx_gas(:,:)

     ! These are the scale start/stop AS GIVEN BY THE CONFIG
     double precision, allocatable :: gas_retrieve_scale_start(:)
     double precision, allocatable :: gas_retrieve_scale_stop(:)
     double precision, allocatable :: gas_retrieve_scale_cov(:)

     ! These are the scale start/stop positions within the model atmosphere
     integer, allocatable :: s_start(:)
     integer, allocatable :: s_stop(:)

     ! State vector (current), state vector a priori
     double precision, dimension(:), allocatable :: svsv, svap, sver
     double precision, dimension(:,:), allocatable :: sv_ap_cov, sv_post_cov
  end type statevector


  public initialize_statevector, parse_and_initialize_SV, &
       clear_SV

contains



  subroutine clear_SV(SV)
    implicit none
    type(statevector), intent(inout) :: SV

    if (allocated(SV%idx_albedo)) deallocate(SV%idx_albedo)
    if (allocated(SV%idx_sif)) deallocate(SV%idx_sif)
    if (allocated(SV%idx_dispersion)) deallocate(SV%idx_dispersion)
    if (allocated(SV%idx_psurf)) deallocate(SV%idx_psurf)
    if (allocated(SV%idx_solar_shift)) deallocate(SV%idx_solar_shift)
    if (allocated(SV%idx_solar_stretch)) deallocate(SV%idx_solar_stretch)
    if (allocated(SV%idx_zlo)) deallocate(SV%idx_zlo)
    if (allocated(SV%idx_temp)) deallocate(SV%idx_temp)
    if (allocated(SV%idx_ils_stretch)) deallocate(SV%idx_ils_stretch)
    if (allocated(SV%idx_aerosol_aod)) deallocate(SV%idx_aerosol_aod)
    if (allocated(SV%idx_aerosol_height)) deallocate(SV%idx_aerosol_height)

    if (allocated(SV%svap)) deallocate(SV%svap)
    if (allocated(SV%svsv)) deallocate(SV%svsv)
    if (allocated(SV%sver)) deallocate(SV%sver)

    if (allocated(SV%sv_ap_cov)) deallocate(SV%sv_ap_cov)
    if (allocated(SV%sv_post_cov)) deallocate(SV%sv_post_cov)

    if (allocated(SV%aerosol_aod_idx_lookup)) deallocate(SV%aerosol_aod_idx_lookup)

    if (allocated(SV%idx_gas)) deallocate(SV%idx_gas)
    if (allocated(SV%gas_idx_lookup)) deallocate(SV%gas_idx_lookup)

    if (allocated(SV%gas_retrieve_scale_start)) deallocate(SV%gas_retrieve_scale_start)
    if (allocated(SV%gas_retrieve_scale_stop)) deallocate(SV%gas_retrieve_scale_stop)
    if (allocated(SV%gas_retrieve_scale_cov)) deallocate(SV%gas_retrieve_scale_cov)

    if (allocated(SV%s_start)) deallocate(SV%s_start)
    if (allocated(SV%s_stop)) deallocate(SV%s_stop)

    SV%num_albedo = -1
    SV%num_sif = -1
    SV%num_dispersion = -1
    SV%num_psurf = -1
    SV%num_gas = -1
    SV%num_aerosol_aod = -1
    SV%num_solar_shift = -1
    SV%num_solar_stretch = -1
    SV%num_zlo = -1
    SV%num_temp = -1
    SV%num_ils_stretch = -1
    SV%num_aerosol_aod = -1
    SV%num_aerosol_height = -1

  end subroutine clear_SV


  subroutine parse_and_initialize_SV(i_win, num_levels, SV)

    implicit none
    integer, intent(in) :: i_win, num_levels
    type(statevector), intent(inout) :: SV

    character(len=*), parameter :: fname = "parse_and_initialize_statevector"

    type(string), allocatable :: split_string(:)
    type(string), allocatable :: split_sv_string(:)
    type(string), allocatable :: split_svval_string(:)
    type(string), allocatable :: split_aero_string(:)
    type(string) :: check_sv_name
    type(string) :: check_sv_retr_type
    character(len=999) :: tmp_str

    ! Which are the SV elements known by the code, and which one's are used?
    type(string) :: known_SV(99)
    logical :: is_gas_SV(99)
    logical :: is_aerosol_SV(99)

    integer :: i, j, k
    integer :: gas_index
    integer :: aero_index
    integer :: known_SV_max
    integer :: last_known
    logical :: SV_found

    integer, allocatable :: gas_retr_count(:)

    integer :: num_albedo_parameters = 0
    integer :: num_dispersion_parameters = 0
    integer :: num_sif_parameters = 0
    integer :: num_psurf_parameters = 0
    integer :: num_solar_shift_parameters = 0
    integer :: num_solar_stretch_parameters = 0
    integer :: num_zlo_parameters = 0
    integer :: num_temp_parameters = 0
    integer :: num_ils_stretch_parameters = 0
    integer :: num_aerosol_aod_parameters = 0
    integer :: num_aerosol_height_parameters = 0

    known_SV(:) = ""
    is_gas_SV(:) = .false.
    is_aerosol_SV(:) = .false.

    ! How many SV elements do we have per gas? We can retrieve
    ! multiple sub-column VMRs for our gases.
    allocate(gas_retr_count(MCS%window(i_win)%num_gases))
    gas_retr_count(:) = 0

    ! These are the known state vector elements - only these will be properly
    ! parsed. The order in which they are defined is not significant.
    known_SV(1) = "albedo"
    known_SV(2) = "dispersion"
    known_SV(3) = "sif"
    known_SV(4) = "psurf"
    known_SV(5) = "solar_shift"
    known_SV(6) = "solar_stretch"
    known_SV(7) = "zlo"
    known_SV(8) = "temp"
    known_SV(9) = "ils_stretch"
    last_known = 9

    ! Add gases as 'known' state vector elements. CAUTION! There is
    ! obviously a danger if someone decides to name their gas "psurf".
    ! Maybe I should add a list of protected names that can't be used.
    do i=1, MCS%window(i_win)%num_gases
       known_SV(last_known+i) = MCS%window(i_win)%gases(i)
       ! This flags the known SV as one of a gas type
       is_gas_SV(last_known+i) = .true.
    end do
    last_known = last_known + MCS%window(i_win)%num_gases

    ! Do the same for aerosols
    do i=1, MCS%window(i_win)%num_aerosols
       known_SV(last_known + i) = MCS%window(i_win)%aerosols(i)
       is_aerosol_SV(last_known + i) = .true.
    end do

    ! We really only need to look for this many items in the list
    ! of known SV's
    known_SV_max = last_known+i

    ! Split the state vector string from the config file
    call MCS%window(i_win)%SV_string%split(tokens=split_string, sep=' ')

    do i=1, size(split_string)
       SV_found = .false.

       ! Loop through the known SV elements, and exit if any requested
       ! one is not found.
       do j=1, known_SV_max
          if (known_SV(j) == "") cycle

          if (split_string(i)%lower() == known_SV(j)%lower()) then
             ! Found a regular, non-gas state vector element
             SV_found = .true.
             exit
          end if

          ! If there is a (one) '|' in the state vector name, we might
          ! have a gas or aerosol SV element
          if (split_string(i)%count('|') >= 1) then

             ! Split that string, so we can grab the gas bit before the
             ! dash, and check that one against the list of known SVs
             call split_string(i)%split(tokens=split_sv_string, sep='|')
             if (split_sv_string(1)%lower() == known_SV(j)%lower()) then
                SV_found = .true.
                deallocate(split_sv_string)
                exit
             end if

             deallocate(split_sv_string)
          end if

       end do

       if (.not. SV_found) then
          call logger%fatal(fname, "This SV element was not recognized: " // &
               trim(split_string(i)%chars()))
          call logger%fatal(fname, "Make sure your SV gases and aerosols are present in the window!")
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

    ! Again, we can do these only if the arrays are allocated,
    ! which they aren't if no gases are defined..
    if (MCS%window(i_win)%num_gases > 0) then
       MCS%window(i_win)%gas_retrieved(:) = .false.
       MCS%window(i_win)%gas_retrieve_profile(:) = .false.
       MCS%window(i_win)%gas_retrieve_scale(:) = .false.
       MCS%window(i_win)%gas_retrieve_scale_start(:,:) = -1.0d0
       MCS%window(i_win)%gas_retrieve_scale_stop(:,:) = -1.0d0
    end if

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
          if (MCS%window(i_win)%dispersion_order < 0) then
             call logger%fatal(fname, "We are retrieving dispersion, but the dispersion order " &
                  // "needs to be >= 0. Check if you've supplied a sensible value (or at all).")
             stop 1
          else
             ! Dispersion order 1 means shift, 2 is stretch etc., this is not quite
             ! consistent with the albedo order notation, but whatever.
             num_dispersion_parameters = MCS%window(i_win)%dispersion_order + 1

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

             if (num_dispersion_parameters > size(MCS%window(i_win)%dispersion_cov)) then
                call logger%fatal(fname, "Not enough disperison covariance values!")
                stop 1
             end if
          end if
       end if

       ! We are retrieving ILS stretch!
       if (split_string(i)%lower() == "ils_stretch") then
          ! Stretch order needs to be > 0
          if (MCS%window(i_win)%ils_stretch_order < 0) then
             call logger%fatal(fname, "We are retrieving ILS stretch, but the ILS stretch order " &
                  // "needs to be >= 0. Check if you've supplied a sensible value (or at all).")
             stop 1
          else
             num_ils_stretch_parameters = MCS%window(i_win)%ils_stretch_order + 1

             ! We MUST have at least the same number of dispersion perturbation
             ! elements.
             if (.not. allocated(MCS%window(i_win)%ils_stretch_pert)) then
                call logger%fatal(fname, "ILS perturbation not in config file!")
                stop 1
             end if

             if (num_ils_stretch_parameters > size(MCS%window(i_win)%ils_stretch_pert)) then
                call logger%fatal(fname, "Not enough ILS perturbation values!")
                stop 1
             end if

             if (num_ils_stretch_parameters > size(MCS%window(i_win)%ils_stretch_cov)) then
                call logger%fatal(fname, "Not enough ILS covariance values!")
                stop 1
             end if

          end if
       end if

       ! We are retrieving solar shift
       if (split_string(i)%lower() == "solar_shift") then
          num_solar_shift_parameters = 1
       end if

       ! We are retrieving solar stretch
       if (split_string(i)%lower() == "solar_stretch") then
          num_solar_stretch_parameters = 1
       end if

       ! We are retrieving SIF!
       if (split_string(i)%lower() == "sif") then
          num_sif_parameters = 1
       end if

       ! We are retrieving ZLO! (this is the SAME as SIF essentially)
       if (split_string(i)%lower() == "zlo") then
          num_zlo_parameters = 1
       end if

       ! We are retrieving a temperature offset!
       if (split_string(i)%lower() == "temp") then
          num_temp_parameters = 1
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


       ! Now check for the gases and aerosols
       do j=1, known_SV_max
          ! These are case-sensitive (not sure why?)

          ! We tell the algorithm to retrieve either a scaling factor,
          ! or a full profile by attaching a "|scale" or "|profile". So when
          ! checking for gases, we must drop that part of the string.

          if (split_string(i)%count("|") >= 1) then
             call split_string(i)%split(tokens=split_sv_string, sep='|')
             check_sv_name = split_sv_string(1)
             check_sv_retr_type = split_sv_string(2)
          else
             check_sv_name = split_string(i)
             check_sv_retr_type = ""
          end if

          ! Is this SV element a gas-type state vector?
          if ((check_sv_name == known_SV(j)) .and. (is_gas_SV(j))) then
             ! Yes, it is - now we look which particular gas this is, and set the
             ! retrieved-flag accordingly.

             do k=1, MCS%window(i_win)%num_gases
                gas_index = MCS%window(i_win)%gas_index(k)
                ! If this gas is not used in this window, skip!
                if (.not. MCS%gas(gas_index)%used) cycle
                ! Otherwise, loop through the gases we have in the window
                ! and set them as 'retrieved'
                if (MCS%window(i_win)%gases(k) == check_sv_name) then
                   MCS%window(i_win)%gas_retrieved(k) = .true.

                   ! Depending on the type, we must set the profile/scale retrieval
                   ! flags accordingly. The default setting is a profile retrieval.
                   if ((check_sv_retr_type == "profile") .or. (check_sv_retr_type == "")) then
                      write(tmp_str, '(A, A, A)') "We are retrieving gas ", check_sv_name%chars(), &
                           " as a full profile."
                      call logger%debug(fname, trim(tmp_str))
                      MCS%window(i_win)%gas_retrieve_profile(k) = .true.
                   else if (check_sv_retr_type == "scale") then
                      write(tmp_str, '(A, A, A)') "We are retrieving gas ", check_sv_name%chars(), &
                           " as scale factor(s)."
                      call logger%debug(fname, trim(tmp_str))
                      ! Increment gas retrieval counter
                      gas_retr_count(k) = gas_retr_count(k) + 1

                      MCS%window(i_win)%gas_retrieve_scale(k) = .true.
                      call split_sv_string(3)%split(tokens=split_svval_string, sep=':')

                      if (size(split_svval_string) /= 3) then
                         call logger%fatal(fname, "Error in gas-scale string. Must have exactly 2 vertical bars.")
                         stop 1
                      end if

                      ! Stick the gas scalar start and end psurf factors into the MCS
                      write(tmp_str, *) split_svval_string(1)%chars()
                      read(tmp_str, *) MCS%window(i_win)%gas_retrieve_scale_start(k, gas_retr_count(k))
                      write(tmp_str, *) split_svval_string(2)%chars()
                      read(tmp_str, *) MCS%window(i_win)%gas_retrieve_scale_stop(k, gas_retr_count(k))
                      write(tmp_str, *) split_svval_string(3)%chars()
                      read(tmp_str, *) MCS%window(i_win)%gas_retrieve_scale_cov(k, gas_retr_count(k))

                      ! Of course the values cannot be outside of the interval [0,1]
                      if ((MCS%window(i_win)%gas_retrieve_scale_start(k, gas_retr_count(k)) < 0.0d0) .or. &
                           (MCS%window(i_win)%gas_retrieve_scale_start(k, gas_retr_count(k)) > 1.0d0)) then
                         call logger%fatal(fname, "Sorry, scalar retrieval boundary needs to be [0,1]!")
                         call logger%fatal(fname, split_string(i)%chars())
                         stop 1
                      end if

                      if ((MCS%window(i_win)%gas_retrieve_scale_stop(k, gas_retr_count(k)) < 0.0d0) .or. &
                           (MCS%window(i_win)%gas_retrieve_scale_stop(k, gas_retr_count(k)) > 1.0d0)) then
                         call logger%fatal(fname, "Sorry, scalar retrieval boundary needs to be [0,1]!")
                         call logger%fatal(fname, split_string(i)%chars())
                         stop 1
                      end if

                   else
                      call logger%fatal(fname, "Gas state vector needs to be stated as " &
                           // "[gas]|profile or [gas]|scale.")
                      call logger%fatal(fname, " .. instead I got: " // split_string(i)%chars())
                      stop 1
                   end if

                end if
             end do
          end if

          ! Is this a known SV element, and it it related to an aerosol?
          if ((check_sv_name == known_SV(j)) .and. (is_aerosol_SV(j))) then

             do k=1, MCS%window(i_win)%num_aerosols
                aero_index = MCS%window(i_win)%aerosol_index(k)

                ! If this aerosol isn't used, then simply skip
                if (.not. MCS%aerosol(aero_index)%used) cycle

                if (MCS%window(i_win)%aerosols(k) == check_sv_name) then

                   if ((check_sv_retr_type == 'aod') .or. &
                        (check_sv_retr_type == 'aod-log')) then

                      write(tmp_str, '(A, A)') "We are retrieving the optical depth from aerosol: ", &
                           check_sv_name%chars()
                      call logger%debug(fname, trim(tmp_str))

                      MCS%window(i_win)%aerosol_retrieve_aod(k) = .true.
                      num_aerosol_aod_parameters = num_aerosol_aod_parameters + 1
                      
                      if (check_sv_retr_type == 'aod-log') then
                         call logger%debug(fname, ".. in log-space.")
                         MCS%window(i_win)%aerosol_retrieve_aod_log(k) = .true.
                      end if

                      call split_sv_string(3)%split(tokens=split_svval_string, sep=':')

                      if (size(split_svval_string) /= 2) then
                         call logger%fatal(fname, "Error in aerosol string. Must have exactly 1 colon.")
                         stop 1
                      end if

                      write(tmp_str, *) split_svval_string(1)%chars()
                      read(tmp_str, *) MCS%window(i_win)%aerosol_prior_aod(k)
                      write(tmp_str, *) split_svval_string(2)%chars()
                      read(tmp_str, *) MCS%window(i_win)%aerosol_aod_cov(k)

                   else if (check_sv_retr_type == 'height') then

                      write(tmp_str, '(A, A)') "We are retrieving the layer height from aerosol: ", &
                           check_sv_name%chars()
                      call logger%debug(fname, trim(tmp_str))

                      MCS%window(i_win)%aerosol_retrieve_height(k) = .true.
                      num_aerosol_height_parameters = num_aerosol_height_parameters + 1

                      call split_sv_string(3)%split(tokens=split_svval_string, sep=':')

                      if (size(split_svval_string) /= 2) then
                         call logger%fatal(fname, "Error in aerosol string. Must have exactly 1 colon.")
                         stop 1
                      end if

                      write(tmp_str, *) split_svval_string(1)%chars()
                      read(tmp_str, *) MCS%window(i_win)%aerosol_prior_height(k)
                      write(tmp_str, *) split_svval_string(2)%chars()
                      read(tmp_str, *) MCS%window(i_win)%aerosol_height_cov(k)

                   else
                      call logger%fatal(fname, "Aerosol state vector needs to be stated as " &
                           // "[aerosol]|aod")
                      call logger%fatal(fname, " .. instead I got: " // split_string(i)%chars())
                      stop 1
                   end if

                end if

             end do

          end if


          if (allocated(split_sv_string)) deallocate(split_sv_string)
       end do

    end do

    ! Once we have all the data, initialize the state vector with the appropriate
    ! number of elements for each type.

    call initialize_statevector( &
         i_win, &
         num_levels, &
         SV, &
         num_albedo_parameters, &
         num_sif_parameters, &
         num_dispersion_parameters, &
         num_psurf_parameters, &
         num_solar_shift_parameters, &
         num_solar_stretch_parameters, &
         num_zlo_parameters, &
         num_temp_parameters, &
         num_ils_stretch_parameters, &
         num_aerosol_aod_parameters, &
         num_aerosol_height_parameters, &
         gas_retr_count)

  end subroutine parse_and_initialize_SV




  subroutine initialize_statevector(i_win, num_levels, sv, &
       count_albedo, count_sif, count_dispersion, count_psurf, &
       count_solar_shift, count_solar_stretch, count_zlo, &
       count_temp, count_ils_stretch, count_aerosol_aod, &
       count_aerosol_height, &
       gas_retr_count)

    implicit none
    integer, intent(in) :: i_win, num_levels
    type(statevector), intent(inout) :: sv
    integer, intent(in) :: count_albedo, count_sif, &
         count_dispersion, count_psurf, count_solar_shift, &
         count_solar_stretch, count_zlo, count_temp, count_ils_stretch, &
         count_aerosol_aod, count_aerosol_height, gas_retr_count(:)

    integer :: count_gas
    character(len=*), parameter :: fname = "initialize_statevector"
    character(len=999) :: tmp_str
    integer :: sv_count
    integer :: i, j, cnt

    ! Set The Number of Statevector parameters
    sv%num_albedo = count_albedo
    sv%num_sif = count_sif
    sv%num_dispersion = count_dispersion
    sv%num_psurf = count_psurf
    sv%num_solar_shift = count_solar_shift
    sv%num_solar_stretch = count_solar_stretch
    sv%num_zlo = count_zlo
    sv%num_temp = count_temp
    sv%num_ils_stretch = count_ils_stretch
    sv%num_aerosol_aod = count_aerosol_aod
    sv%num_aerosol_height = count_aerosol_height
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

    if (sv%num_zlo > 0) then
       write(tmp_str, '(A, G0.1)') "Number of ZLO SV elements: ", sv%num_zlo
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_zlo(sv%num_zlo))
       do i=1, sv%num_zlo
          sv_count = sv_count + 1
          sv%idx_zlo(i) = sv_count
       end do
    else
       allocate(sv%idx_zlo(1))
       sv%idx_zlo(1) = -1
    end if

    if (sv%num_temp > 0) then
       write(tmp_str, '(A, G0.1)') "Number of temperature SV elements: ", sv%num_temp
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_temp(sv%num_temp))
       do i=1, sv%num_temp
          sv_count = sv_count + 1
          sv%idx_temp(i) = sv_count
       end do
    else
       allocate(sv%idx_temp(1))
       sv%idx_temp(1) = -1
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


    if (sv%num_ils_stretch > 0) then

       write(tmp_str, '(A, G0.1)') "Number of ILS stretch SV elements: ", sv%num_dispersion
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_ils_stretch(sv%num_ils_stretch))
       do i=1, sv%num_ils_stretch
          sv_count = sv_count + 1
          sv%idx_ils_stretch(i) = sv_count
       end do
    else
       allocate(sv%idx_ils_stretch(1))
       sv%idx_ils_stretch(1) = -1
    end if

    ! Solar shift and stretch
    allocate(sv%idx_solar_shift(1))
    if (sv%num_solar_shift == 1) then
       write(tmp_str, '(A, G0.1)') "Number of solar shift elements: ", sv%num_solar_shift
       call logger%info(fname, trim(tmp_str))
       sv_count = sv_count + 1
       sv%idx_solar_shift(1) = sv_count
    else
       sv%idx_solar_shift(1) = -1
    end if

    allocate(sv%idx_solar_stretch(1))
    if (sv%num_solar_stretch == 1) then
       write(tmp_str, '(A, G0.1)') "Number of solar stretch elements: ", sv%num_solar_stretch
       call logger%info(fname, trim(tmp_str))
       sv_count = sv_count + 1
       sv%idx_solar_stretch(1) = sv_count
    else
       sv%idx_solar_stretch(1) = -1
    end if

    ! Surface pressure
    if (sv%num_psurf == 1) then
       write(tmp_str, '(A)') "Number of surface pressure SV elements: 1"
       call logger%info(fname, trim(tmp_str))
       allocate(sv%idx_psurf(1))
       sv_count = sv_count + 1
       sv%idx_psurf = sv_count
    end if


    ! Aerosols
    ! AOD
    if (sv%num_aerosol_aod > 0) then

       allocate(sv%aerosol_aod_idx_lookup(sv%num_aerosol_aod))
       sv%aerosol_aod_idx_lookup(:) = -1

       write(tmp_str, '(A, G0.1)') "Number of aerosol AOD elements: ", sv%num_aerosol_aod
       call logger%info(fname, trim(tmp_str))
       allocate(sv%idx_aerosol_aod(sv%num_aerosol_aod))

       cnt = 1

       do i = 1, MCS%window(i_win)%num_aerosols
          if (MCS%window(i_win)%aerosol_retrieve_aod(i)) then

             sv_count = sv_count + 1
             sv%idx_aerosol_aod(i) = sv_count
             sv%aerosol_aod_idx_lookup(cnt) = i

             cnt = cnt + 1

          end if
       end do
    else
       allocate(sv%idx_aerosol_aod(1))
       sv%idx_aerosol_aod(1) = -1
    end if

    ! Aerosol Height
    if (sv%num_aerosol_height > 0) then

       allocate(sv%aerosol_height_idx_lookup(sv%num_aerosol_height))
       sv%aerosol_height_idx_lookup(:) = -1

       write(tmp_str, '(A, G0.1)') "Number of aerosol height elements: ", sv%num_aerosol_height
       call logger%info(fname, trim(tmp_str))
       allocate(sv%idx_aerosol_height(sv%num_aerosol_height))

       cnt = 1

       do i = 1, MCS%window(i_win)%num_aerosols
          if (MCS%window(i_win)%aerosol_retrieve_height(i)) then

             sv_count = sv_count + 1
             sv%idx_aerosol_height(i) = sv_count
             sv%aerosol_height_idx_lookup(cnt) = i

             cnt = cnt + 1

          end if
       end do
    else
       allocate(sv%idx_aerosol_height(1))
       sv%idx_aerosol_height(1) = -1
    end if


    ! Gases

    ! In order to find out which gases are retrieved, we access the MCS
    ! rather than have it passed into this function. We also have to find out
    ! how many levels the atmosphere currently has, if we run profile
    ! retrievals.

    ! Simply just count the number of gases and add either 1
    ! or the number of levels, depending on retrieval type.

    count_gas = sum(gas_retr_count)

    sv%num_gas = count_gas

    if (count_gas > 0) then
       allocate(sv%idx_gas(count_gas, num_levels))
       allocate(sv%gas_idx_lookup(count_gas))
       allocate(sv%gas_retrieve_scale_start(count_gas))
       allocate(sv%gas_retrieve_scale_stop(count_gas))
       allocate(sv%gas_retrieve_scale_cov(count_gas))
       allocate(sv%s_start(count_gas))
       allocate(sv%s_stop(count_gas))

       sv%idx_gas(:,:) = -1
       sv%gas_idx_lookup(:) = -1
       sv%gas_retrieve_scale_start(:) = -1.0d0
       sv%gas_retrieve_scale_stop(:) = -1.0d0
       sv%gas_retrieve_scale_cov(:) = -1.0d0
       sv%s_start(:) = -1.0d0
       sv%s_stop(:) = -1.0d0

       write(tmp_str, '(A,G0.1)') "Number of gas SV elements: ", count_gas
       call logger%info(fname, trim(tmp_str))

       cnt = 1
       do i = 1, MCS%window(i_win)%num_gases
          if (MCS%window(i_win)%gas_retrieved(i)) then

             if (MCS%window(i_win)%gas_retrieve_profile(i)) then
                do j = 1, num_levels
                   sv_count = sv_count + 1
                   sv%idx_gas(cnt, j) = sv_count
                end do

                sv%gas_idx_lookup(cnt) = i
                cnt = cnt + 1

             end if


             if (MCS%window(i_win)%gas_retrieve_scale(i)) then
                ! Retrieving scale factor(s)? Set idx_gas to -1 first.
                do j = 1, gas_retr_count(i)
                   ! Loop through the number of retrieved elements per
                   ! gas to grab all scale limits for each gas i and each
                   ! retrieved parameter j per gas i
                   sv%idx_gas(cnt, :) = -1

                   sv_count = sv_count + 1
                   sv%idx_gas(cnt, 1) = sv_count
                   sv%gas_idx_lookup(cnt) = i

                   sv%gas_retrieve_scale_start(cnt) = &
                        MCS%window(i_win)%gas_retrieve_scale_start(i, j)
                   sv%gas_retrieve_scale_stop(cnt) = &
                        MCS%window(i_win)%gas_retrieve_scale_stop(i, j)
                   sv%gas_retrieve_scale_cov(cnt) = &
                        MCS%window(i_win)%gas_retrieve_scale_cov(i, j)
                   cnt = cnt + 1
                end do
             end if


          end if
       end do

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




end module statevector_mod
