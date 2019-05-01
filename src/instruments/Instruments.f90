!> @brief Generic instrument type module
!> @file Instruments.f90
!> @author Peter Somkuti
!>
!! This is an embarrasingly tiny generic instrument class. It is simply
!! required due to the way Fortran does its run-time polymorphism. Based
!! on this generic type, we define various other instrument types, such
!! as the OCO-2-like instrument in OCO2-Instrument.f90.

module instruments_mod
  type, abstract :: generic_instrument
     ! There is really nothing to be put here..
  end type generic_instrument
end module instruments_mod
