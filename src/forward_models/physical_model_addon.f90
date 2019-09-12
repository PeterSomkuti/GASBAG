

module physical_model_addon_mod

  use math_utils_mod

  ! Third-party modules
  use stringifor
  use mod_datetime
  use logger_mod, only: logger => master_logger

  ! System modules
  use ISO_FORTRAN_ENV

  implicit none

  !> A simple structure to keep the atmosphere data nice and tidy
  type atmosphere
     integer :: num_levels
     integer :: num_gases
     !> The name(s) of the gas(es), (gas number)
     type(string), allocatable :: gas_names(:)
     !> Gas mixing ratios (level, gas number)
     double precision, allocatable :: gas_vmr(:,:)
     !> To which spectroscopy data does this gas correspond to? (gas number)
     integer, allocatable :: gas_index(:)
     !> Model atmosphere temperature
     double precision, allocatable :: T(:)
     !> Model atmosphere pressure
     double precision, allocatable :: p(:)
     !> Model atmosphere specific humidity
     double precision, allocatable :: sh(:)
     !> Model atmosphere gravity at pressure levels
     double precision, allocatable :: grav(:)
     !> Model atmosphere altitude at pressure levels
     double precision, allocatable :: altitude_levels(:)
     !> Model atmosphere dry air column at layers
     double precision, allocatable :: ndry(:)
  end type atmosphere

  !> Holds optical gas/aerosol properties of the scene
  type optical_properties
     !> Wavelengths
     double precision, allocatable :: wl(:)
     !> Surface albedo
     double precision, allocatable :: albedo(:)

     !> Gas optical depth (spectral, layer, gas number)
     double precision, allocatable :: gas_tau(:,:,:)
     !> dtau / dtemp (spectral, layer, gas number)
     double precision, allocatable :: gas_tau_dtemp(:,:,:)
     !> dtau / dpsurf (spectral, layer, gas_number)
     double precision, allocatable :: gas_tau_dpsurf(:,:,:)
     !> Perturbed gas optical depth (spectral, layer, gas number)
     double precision, allocatable :: gas_tau_pert(:,:,:,:)

     !> Rayleigh extinction optical depth (spectral, layer)
     double precision, allocatable :: ray_tau(:,:)
     !> Rayleigh depolarization factor (spectral)
     double precision, allocatable :: ray_depolf(:)

     !> Total column optical depth (spectral)
     double precision, allocatable :: total_tau(:)
     !> Total single scatter albedo (spectral, layer)
     double precision, allocatable :: omega(:,:)

  end type optical_properties


  ! The "biggest" custom type that contains most per-scene
  ! information you'll ever need to run a physical retrieval.
  type, public :: scene
     !> The atmosphere of the scene
     type(atmosphere) :: atm
     !> Optical properties needed for RT
     type(optical_properties) :: op
     !> Date object
     type(datetime) :: date
     !> Epoch (Year, Month, Day, Hour, Minute, Second, Millisecond)
     integer :: epoch(7)
     !> Longitude
     double precision :: lon
     !> Latitude
     double precision :: lat
     !> Altitude at surface
     double precision :: alt
     !> Solar zenith angle
     double precision :: SZA
     !> Cosine of solar zenith angle
     double precision :: mu0
     !> Viewing zenith angle
     double precision :: VZA
     !> Cosine of viewing zenith angle
     double precision :: mu
     !> Solar azimuth angle
     double precision :: SAA
     !> Viewing azimuth angle
     double precision :: VAA
  end type scene

  public scene_altitude

contains


  subroutine allocate_optical_properties(scn, N_hires, N_gases)

    type(scene), intent(inout) :: scn
    integer, intent(in) :: N_hires
    integer, intent(in) :: N_gases

    allocate(scn%op%wl(N_hires))
    scn%op%wl(:) = 0.0d0
    allocate(scn%op%albedo(N_hires))
    scn%op%albedo(:) = 0.0d0

    ! These arrays are only allocated if we have an atmosphere with at
    ! least one gas in it.
    if (N_gases > 0) then
       allocate(scn%op%gas_tau(N_hires, scn%atm%num_levels-1, N_gases))
       scn%op%gas_tau(:,:,:) = 0.0d0
       allocate(scn%op%gas_tau_dpsurf(N_hires, scn%atm%num_levels-1, N_gases))
       scn%op%gas_tau_dpsurf(:,:,:) = 0.0d0
       allocate(scn%op%gas_tau_dtemp(N_hires, scn%atm%num_levels-1, N_gases))
       scn%op%gas_tau_dtemp(:,:,:) = 0.0d0
       allocate(scn%op%ray_tau(N_hires, scn%atm%num_levels-1))
       scn%op%ray_tau(:,:) = 0.0d0
       allocate(scn%op%ray_depolf(N_hires))
       scn%op%ray_depolf(:) = 0.0d0
       allocate(scn%op%total_tau(N_hires))
       scn%op%total_tau(:) = 0.0d0
       allocate(scn%op%omega(N_hires, scn%atm%num_levels-1))
    end if

  end subroutine allocate_optical_properties


  subroutine destroy_optical_properties(scn)
    type(scene), intent(inout) :: scn

    deallocate(scn%op%wl)
    deallocate(scn%op%albedo)

    if (allocated(scn%op%gas_tau)) deallocate(scn%op%gas_tau)
    if (allocated(scn%op%gas_tau_dpsurf)) deallocate(scn%op%gas_tau_dpsurf)
    if (allocated(scn%op%gas_tau_dtemp)) deallocate(scn%op%gas_tau_dtemp)
    if (allocated(scn%op%ray_tau)) deallocate(scn%op%ray_tau)
    if (allocated(scn%op%ray_depolf)) deallocate(scn%op%ray_depolf)
    if (allocated(scn%op%total_tau)) deallocate(scn%op%total_tau)
    if (allocated(scn%op%omega)) deallocate(scn%op%omega)

  end subroutine destroy_optical_properties

  !> @brief Calculate altitude levels, gravity levels and ndry
  !> @detail This code is borrowed from MS3 (O'Brien, O'Dell et al.)
  !> @param scn Scene object
  subroutine scene_altitude(scn)
    type(scene), intent(inout) :: scn

    ! Local stuff
    double precision :: SH_layer, g_layer, p_layer, T_layer, Tv, dP
    double precision :: logratio, dz, constant
    ! Loop variable
    integer :: i

    ! Set altitudes to zero, and the lowest level to the altitude
    scn%atm%altitude_levels(:) = 0.0d0
    scn%atm%altitude_levels(scn%atm%num_levels) = scn%alt
    scn%atm%ndry(:) = 0.0d0

    ! Loop through layers, starting with the bottom-most (surface) one
    do i = scn%atm%num_levels - 1, 1, -1

       SH_layer = (scn%atm%sh(i) + scn%atm%sh(i+1)) * 0.5d0
       g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1))
       p_layer = (scn%atm%p(i) + scn%atm%p(i+1)) * 0.5d0
       dP = scn%atm%p(i+1) - scn%atm%p(i)
       T_layer = (scn%atm%T(i) + scn%atm%T(i+1)) * 0.5d0
       Tv = T_layer * (1.0d0 + SH_layer * (1.0d0 - EPSILON) / EPSILON)
       logratio = log(scn%atm%p(i+1) / scn%atm%p(i))
       dz = logratio * Tv * Rd / g_layer
       g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1) + 0.5d0 * dz)
       dz = logratio * Tv * Rd / g_layer
       constant = dP / (DRY_AIR_MASS * g_layer)

       ! Write the important stuff back into the
       !----------------------------------------------
       scn%atm%grav(i) = jpl_gravity(scn%lat, scn%atm%altitude_levels(i))
       scn%atm%ndry(i) = constant * (1.0d0 - SH_layer)
       scn%atm%altitude_levels(i) = scn%atm%altitude_levels(i+1) + dz
       !----------------------------------------------

    end do

  end subroutine scene_altitude

  !> @brief Calculate acceleration g(lat, altitude)
  !> @param gdlat Geodetic latitude [deg]
  !> @param altit Geometric altitude [m]
  pure function jpl_gravity(gdlat,altit) result(gravity)
    ! NOTE
    ! This code was taken from MS3 and just modified slightly for data
    ! types.
    !
    ! Computes the effective Earth gravity at a given latitude and altitude.
    ! This is the sum of the gravitational and centripital accelerations.
    ! These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
    ! The Earth is assumed to be an oblate ellipsoid, with a ratio of the
    ! major to minor axes = sqrt(1+con) where con=.006738
    ! This eccentricity makes the Earth's gravititational field smaller at the
    ! poles and larger at the equator than if the Earth were a sphere of the
    ! same mass. It also makes the local mid-latitude gravity field not point
    ! toward the center of mass.
    !
    ! Input Parameters:
    ! gdlat       GeoDetric Latitude (degrees)
    !	altit       Geometric Altitude (meters) ! CHANGED CWO 6/22/2009
    !
    ! Output Parameter:
    !	gravity     Effective Gravitational Acceleration (m/s2)
    !
    ! Interestingly, since the centripital effect of the Earth's rotation
    ! (-ve at equator, 0 at poles) has almost the opposite shape to the
    ! second order gravitational field (+ve at equator, -ve at poles), their
    ! sum is almost constant so that the surface gravity can be approximated
    ! (to .07%) by the simple expression g = 0.99746*GM/radius**2, the latitude
    ! variation coming entirely from the variation of surface r with latitude.

    implicit none

    ! In/Out variables
    double precision, intent(in) :: gdlat   ! geodetic latitude [degrees]
    double precision, intent(in) :: altit   ! geometric altitude [meters]
    double precision :: gravity ! gravitational acceleration [m/s^2]

    ! Local Variables
    double precision :: radius        ! radial distance (metres)
    double precision :: gclat         ! geocentric latitude [radians]
    double precision :: ff, hh, ge       ! scratch variables

    ! Parameters
    ! Gravitational constant times Earth's Mass (m3/s2)
    double precision, parameter  :: gm = 3.9862216d14
    ! Earth's angular rotational velocity (radians/s)
    double precision, parameter  :: omega = 7.292116d-5
    ! (a/b)**2-1 where a & b are equatorial & polar radii
    double precision, parameter  :: con = 0.006738d0
    ! 2nd harmonic coefficient of Earth's gravity field
    double precision, parameter  :: shc = 1.6235d-3

    ! Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
    gclat=atan(tan(DEG2RAD*gdlat)/(1+con))  ! radians
    ! On computers which crash at the poles try the following expression
    ! gclat=d2r*gdlat-con*sin(d2r*gdlat)*cos(d2r*gdlat)/(1+con*cos(d2r*gdlat)**2)
    radius= altit + EARTH_EQUATORIAL_RADIUS/sqrt(1.0d0+con*sin(gclat)**2)
    ff=(radius/ EARTH_EQUATORIAL_RADIUS)**2
    hh=radius*omega**2
    ge=gm/ EARTH_EQUATORIAL_RADIUS**2                       ! = gravity at Re
    gravity=(ge*(1-shc*(3.0d0*sin(gclat)**2-1.0d0)/ff)/ff-hh*cos(gclat)**2) &
         *(1.0d0+0.5d0*(sin(gclat)*cos(gclat)*(hh/ge+2.0d0*shc/ff**2))**2)

  end function jpl_gravity


  !> @brief Reads in the atmosphere file, which contains column-based profiles.
  !>
  !> @param filename Path to the atmosphere file
  !> @param Array of 'strings' containing names of required gases
  !> @param Atmosphere object that will be populated by this function
  !>
  !> @detail We supply here a filename, a list of gas strings and the
  !> atmosphere-object. The function checks whether all required gases
  !> are present in the atmosphere file, and will throw a fatal error if
  !> that is not the case.
  subroutine read_atmosphere_file(filename, gas_strings, atm)

    implicit none
    character(len=*), intent(in) :: filename
    type(string), intent(in) :: gas_strings(:)
    type(atmosphere), intent(inout) :: atm

    ! Function name
    character(len=*), parameter :: fname = "read_atmosphere_file"
    ! File handler and IO stat variable
    integer :: funit, iostat
    ! Whether file exists
    logical :: file_exist
    ! Various counting variables to figure out where the
    ! contents of the atmosphere file start.
    integer :: line_count, nonempty_line_count, file_start, level_count
    ! Indices which tell us which columns are pressure and temperature
    integer :: idx_p, idx_t
    ! Character variables to store lines
    character(len=999) :: dummy, tmp_str
    ! This dummy is for reading in numerical values
    double precision :: dummy_dp
    ! This dummy is for reading in strings
    type(string) :: dummy_string
    ! Lines will be split using this dummy variable
    type(string), allocatable :: split_string(:)
    ! Various loop counters
    integer :: i, j, cnt
    ! Too keep track of whether we have found the required gases, and
    ! at which position/column they are.
    integer, allocatable :: this_gas_index(:)

    integer :: num_gases

    ! Check whether the file exists or not.
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Atmosphere file does not exist: " // filename)
       stop 1
    else
       call logger%debug(fname, "File does exist.")
    end if

    ! First pass: we scan the file to see how many levels our atmosphere has
    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    level_count = 0
    nonempty_line_count = 0
    file_start = -1

    ! Loop through the file until we have reached the end
    do
       read(funit, '(A)', iostat=iostat) dummy

       if (iostat == iostat_end) then
          ! End of file?
          exit
       end if

       ! Keep track of how many lines we have traversed
       line_count = line_count + 1

       if (scan(dummy, "!#;") > 0) then
          ! Skip this line, as it's commented
          cycle
       else if (trim(dummy) == "") then
          ! Skip empty lines
          cycle
       end if

       ! And keep track of now many lines are non-empty
       nonempty_line_count = nonempty_line_count + 1

       ! First non-empty line should be saved here.
       ! file_start is where the real file contents begin
       if (file_start == -1) then
          file_start = line_count
       else
          if (nonempty_line_count > 1) then
             level_count = level_count + 1
          end if
       end if

    end do

    ! Go back to the top of the file.
    rewind(unit=funit, iostat=iostat)

    idx_p = -1
    idx_t = -1

    ! Reset the line counter since we're starting from the top again.
    line_count = 0
    do
       ! Read the line
       read(funit, '(A)', iostat=iostat) dummy
       line_count = line_count + 1

       ! .. and immediately skip until we are at the
       ! beginning of the contents of the file.
       if (line_count < file_start) cycle

       if (line_count == file_start) then
          ! This is the proper atmosphere header that should contain
          ! the information about the gases. So first, let's check if
          ! the numbers match

          ! Read the line into a string object and split it by whitespaces
          dummy_string = dummy
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now that we know both the number of levels and gases, we can allocate the
          ! arrays in the atmosphere structure.
          num_gases = 0
          idx_p = -1
          idx_t = -1
          do j=1, size(split_string)
             ! Skip temp or pressure - this requires the pressure and temperature
             ! columns to be explicitly labeled p/P and t/T
             if (split_string(j)%lower() == "p") then
                idx_p = j
                cycle
             end if
             if (split_string(j)%lower() == "t") then
                idx_t = j
                cycle
             end if
             num_gases = num_gases + 1
          end do

          ! Let the user know how many gases and levels we have found
          write(tmp_str, '(A,G0.1,A,A)') "There seem to be ", num_gases, " gases in ", filename
          call logger%info(fname, trim(tmp_str))
          write(tmp_str, '(A, G0.1)') "The number of atmospheric levels is: ", level_count
          call logger%info(fname, trim(tmp_str))

          ! If the atmosphere structure was already allocated, deallocate it first
          if (allocated(atm%p)) deallocate(atm%p)
          if (allocated(atm%T)) deallocate(atm%T)
          if (allocated(atm%sh)) deallocate(atm%sh)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)
          if (allocated(atm%altitude_levels)) deallocate(atm%altitude_levels)
          if (allocated(atm%grav)) deallocate(atm%grav)
          if (allocated(atm%ndry)) deallocate(atm%ndry)

          ! Allocate according to the file structure
          atm%num_levels = level_count
          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))
          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))
          allocate(atm%altitude_levels(level_count))
          allocate(atm%grav(level_count))
          allocate(atm%ndry(level_count))

          allocate(this_gas_index(size(gas_strings)))

          ! But we also want to know what gas index to use for storage
          do i=1, size(gas_strings)
             this_gas_index(i) = -1
             cnt = 1
             do j=1, size(split_string)
                ! Skip temp or pressure
                if (split_string(j)%lower() == "p") cycle
                if (split_string(j)%lower() == "t") cycle
                ! Check if gas description matches gases we know
                if (split_string(j) == gas_strings(i)) then
                   this_gas_index(i) = j
                   write(tmp_str, '(A,A,A,G0.1)') "Index for atmosphere gas ", &
                        split_string(j)%chars(), ": ", j
                   call logger%debug(fname, trim(tmp_str))
                   exit
                end if
                cnt = cnt + 1
             end do
          end do

          ! And last check - do all required gases have a VMR column in the
          ! atmosphere file.

          do i=1, size(gas_strings)
             if (this_gas_index(i) == -1) then
                ! Uh-oh, gas that was speficied in the "gases" option of the window
                ! could not be found in the atmosphere! Hard exit.
                write(tmp_str, "(A,A)") "The following gas was not found in the atmosphere file: " &
                     // gas_strings(i)%chars()
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if
          end do

       end if


       if (line_count > file_start) then
          ! Right after the header, we should have the data in rows. We
          ! use the string split option to split the row string into substrings,
          ! and then convert each into double precision values and feed them into
          ! the atmosphere structure - at the right position!

          dummy_string = dummy
          ! Need to deallocate the split_string object first
          if (allocated(split_string)) deallocate(split_string)
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now here we need to check again whether a certain line has more than
          ! num_gases+1 columns.
          if (size(split_string) /= (num_gases + 2)) then
             write(tmp_str, '(A, G0.1)') "Too many values in line ", line_count
             call logger%fatal(fname, trim(tmp_str))
             stop 1
          end if

          ! Get the pressure value
          tmp_str = split_string(idx_p)%chars()
          read(tmp_str, *) dummy_dp
          atm%p(line_count - file_start) = dummy_dp

          ! Get the temperature value
          tmp_str = split_string(idx_t)%chars()
          read(tmp_str, *) dummy_dp
          atm%T(line_count - file_start) = dummy_dp

          ! And the gas value(s) from the other column(s)
          do i=1, size(gas_strings)
             tmp_str = split_string(this_gas_index(i))%chars()
             read(tmp_str, *) dummy_dp
             atm%gas_vmr(line_count - file_start, i) = dummy_dp
          end do

       end if

       if (line_count == (file_start + level_count)) exit

    end do

    close(funit)

  end subroutine read_atmosphere_file


end module physical_model_addon_mod
