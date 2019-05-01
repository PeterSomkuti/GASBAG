MODULE doppler_solar_module


IMPLICIT NONE

  double precision :: LOCAL_EARTH_ROTATION_SPEED ! in m/s
  double precision :: local_earth_radius ! local radius of the earth in meters
  double precision :: geocentric_latitude ! local geocentric latitude in radians
  double precision, parameter :: TWO = 2.0d0
  double precision, parameter :: con = 0.006738d0 ! (a/b)**2-1 where a & b are equatorial & polar radii
  double precision, parameter    :: AU = 149597870691.0 ! 1 AU (meters)
  double precision, parameter :: SOLAR_RADIUS = 695.5d6  ! Solar Radius in meters
  double precision, parameter :: PI = 3.14159265358979323846D0
  double precision, parameter :: d2r = PI / 180.0d0
  double precision, parameter :: SPEED_OF_LIGHT = 2.99792458D8 ! m/s
  double precision, parameter :: EARTH_EQUATORIAL_RADIUS = 6378178.0d0 ! Equatorial Radius (m)
  double precision, parameter :: EARTH_ROTATION_PERIOD = 86164.1d0 ! # of seconds for the earth to rotate once
  ! Boltzmann's constant, in J K^{-1}
  double precision, parameter :: BOLTZMANN = 1.3806505D-23

  ! Avogadro's number, molecules/mole
  double precision, parameter :: NA = 6.0221415D23

  ! Planck's constant, in J s
  double precision, parameter :: PLANCK = 6.6260693D-34
  double precision, parameter :: earth_rotation_direction = 90.0d0 ! i.e., east
  integer, parameter :: solar_distance_version = 2


PRIVATE
PUBLIC :: solar_doppler_velocity

CONTAINS

  SUBROUTINE solar_doppler_velocity(solzen, solaz, epoch, &
       latitude, elevation, &
       solar_velocity_wrt_ground, solar_distance)
!     Returns the velocity of the sun with respect to the
!     observation point on the ground.  Positive is motion
!     of the sun TOWARDS the earth, and results in a blue shift.
!     Negative results in a red shift.
!
!     This velocity includes two effects:
!     1) Radial component of the earth's revolution velocity about the sun
!     2) Local rotation of the earth wrt to incoming sunlight direction
!
!     Note : All velocities are in m/s.

      implicit none

      double precision, intent(out) :: solar_velocity_wrt_ground ! m/s
      double precision, intent(out) :: solar_distance ! earth-sun distance in meters
      double precision, intent(in) :: solzen, solaz
      double precision, intent(in) :: latitude, elevation
      double precision :: solar_velocity, earth_rot_wrt_sun

      integer, dimension(7) :: epoch

      !solzen = scene%angles%sol_zenith
      !solaz = scene%angles%sol_azimuth
      !epoch = scene%epoch

      select case(solar_distance_version)
      case(1)
        CALL solar_distance_and_velocity(epoch, solar_distance, solar_velocity)
      case(2)
        CALL solar_distance_and_velocity_v2(epoch, solar_distance, solar_velocity)
      end select

!     Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
      geocentric_latitude = atan(tan(d2r*latitude)/(1+con))  ! radians
      local_earth_radius = elevation + EARTH_EQUATORIAL_RADIUS/sqrt(1.0d0+con*sin(geocentric_latitude)**2)
!     Calculate local rotation velocity in m/s :
      LOCAL_EARTH_ROTATION_SPEED = (TWO*PI*local_earth_radius)/EARTH_ROTATION_PERIOD*abs(cos(geocentric_latitude))
!     Calculate this rate with respect to the solar direction
      earth_rot_wrt_sun = LOCAL_EARTH_ROTATION_SPEED * sin(solzen*d2r) * cos(d2r*(solaz - earth_rotation_direction))

!!$      if (op%gen%verbose) then
!!$          write(op%files%lunit, *)
!!$          write(op%files%lunit, *) 'Calculating solar velocity for doppler shifting solar spectrum.'
!!$          write(op%files%lunit, *) 'Velocity of earth towards sun [m/s] = ', solar_velocity
!!$          write(op%files%lunit, *) 'Component of earth rotational velocity towards sun [m/s] = ', earth_rot_wrt_sun
!!$      endif
      solar_velocity_wrt_ground = solar_velocity + earth_rot_wrt_sun

      RETURN
      END SUBROUTINE solar_doppler_velocity


      SUBROUTINE solar_distance_and_velocity(epoch, solar_distance, solar_velocity)

      integer, dimension(:), intent(in) :: epoch ! 7-element vector: year, month, day, hour, min, sec, millisecond
      double precision, intent(out) :: solar_distance
      double precision, intent(out), OPTIONAL :: solar_velocity

      double precision, dimension(0:2), parameter :: a = (/ 1.000110d0, 0.034221d0, 0.000719d0 /)
      double precision, dimension(2), parameter :: b = (/ 0.001280d0, 0.000077d0 /)
      integer, dimension(12), parameter :: days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      double precision :: doy, ndays, F, F_t, c1, c2, s1, s2, thetay
      logical :: leap_year
      integer :: m

!     determine if number of days in this year
      leap_year = .FALSE.
      if (MOD(epoch(1),4) == 0) leap_year = .TRUE.
      if (MOD(epoch(1),100) == 0) leap_year = .FALSE.
      if (MOD(epoch(1),400) == 0) leap_year = .TRUE.
      ndays = 365.0d0
      if (leap_year) ndays = 366.0d0

!     determine day of year, 0-ndays
      doy = 0.0d0 ! this will be >= 0, <= ndays
      do m = 1, epoch(2)- 1 ! add days in preceding months
          doy = doy + days_in_month(m)
          if (m == 2 .AND. leap_year) doy = doy + 1.0d0
      enddo
      doy = doy + (epoch(3)-1) ! add preceding days
      doy = doy + (epoch(4) + (epoch(5)+epoch(6)/60.0d0)/60.0d0)/24.0d0

!     determine angular position of earth in its revolution about sun, in radians
      thetay = TWO*PI*doy/ndays
      s1 = sin(thetay)
      s2 = sin(TWO*thetay)
      c1 = cos(thetay)
      c2 = cos(TWO*thetay)
      F = a(0) + a(1)*c1 + a(2)*c2 + b(1)*s1 + b(2)*s2
      solar_distance = AU / sqrt(F)  ! solar distance in m

      if (PRESENT(solar_velocity)) then
!         calculate derivative if F wrt time in seconds
          F_t = (b(1)*c1 + TWO*b(2)*c2 - a(1)*s1 - TWO*a(2)*s2) * (TWO*PI)/(ndays*86400.0d0)
          solar_velocity = solar_distance / (TWO*F) * F_t ! speed of sun TOWARDS the earth, in m/s
      endif

      RETURN
      END SUBROUTINE solar_distance_and_velocity

      SUBROUTINE solar_distance_and_velocity_v2(epoch, solar_distance, solar_velocity)

      ! Uses a much improved solar ephemeris from http://eclipse.gsfc.nasa.gov/TYPE/TYPE.html
      ! Hartmut did a 6th order polynomial fit in day-of-year to tabulated emphemeris data.
      ! therefore, the leap-year effect is completely ignored and will result in slight errors.

      integer, dimension(:), intent(in) :: epoch ! 7-element vector: year, month, day, hour, min, sec, millisecond
      double precision, intent(out) :: solar_distance
      double precision, intent(out), OPTIONAL :: solar_velocity

      double precision, dimension(0:6), parameter :: &
                a = (/ 0.98334d0, -1.82823d-5, 2.30179d-6, 6.62402d-9, &
                       -1.33287d-10, 3.98445d-13, -3.54239d-16 /)

      integer, dimension(12), parameter :: days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      double precision :: doy, F, F_t, doy_i(6)
      logical :: leap_year
      integer :: m, i

!     determine if number of days in this year
      leap_year = .FALSE.
      if (MOD(epoch(1),4) == 0) leap_year = .TRUE.
      if (MOD(epoch(1),100) == 0) leap_year = .FALSE.
      if (MOD(epoch(1),400) == 0) leap_year = .TRUE.

!     determine day of year, 1-366
      doy = 0.0d0 ! this will be >= 0, <= ndays
      do m = 1, epoch(2)- 1 ! add days in preceding months
          doy = doy + days_in_month(m)
          if (m == 2 .AND. leap_year) doy = doy + 1.0d0
      enddo
      doy = doy + epoch(3) ! add preceding days
      doy = doy + (epoch(4) + (epoch(5)+epoch(6)/60.0d0)/60.0d0)/24.0d0
      doy = doy + 0.5d0 ! this makes Hartmut's formula agree with truth together.

      ! calculate powers of "doy"
      doy_i(1) = doy
      do i = 2, 6
        doy_i(i) = doy_i(i-1) * doy
      enddo

      F = a(0) + a(1) * doy_i(1) + a(2) * doy_i(2) &
               + a(3) * doy_i(3) + a(4) * doy_i(4) &
               + a(5) * doy_i(5) + a(6) * doy_i(6)

      solar_distance = F * AU  ! solar distance in m

      if (PRESENT(solar_velocity)) then
!         calculate derivative if F wrt time in seconds
          F_t = a(1) + 2.d0 * a(2) * doy_i(1) + 3.d0 * a(3) * doy_i(2) &
               +  4.d0*a(4)*doy_i(3) + 5.d0*a(5)*doy_i(4) + 6.d0*a(6)*doy_i(5)
          solar_velocity = F_t * (AU / 86400.0d0) ! now speed in m/s AWAY FROM EARTH
          solar_velocity = -solar_velocity ! now speed in m/s TOWARDS EARTH
      endif

      RETURN
      END SUBROUTINE solar_distance_and_velocity_v2


END MODULE doppler_solar_module
