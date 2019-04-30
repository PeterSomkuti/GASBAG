module spectroscopy_utils_mod

  use math_utils_mod
  use control_mod, only: MCS, CS_gas
  use logger_mod, only: logger => master_logger

contains

  subroutine regrid_spectroscopy(gas, wl_grid)
    implicit none
    type(CS_gas), intent(inout) :: gas
    double precision, intent(in) :: wl_grid(:)

    integer :: N_wl, d1, d2, d3, d4, wl_idx_left
    double precision :: fac
    double precision, allocatable :: tmp_array(:,:,:,:)
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "regrid_spectroscopy"

    N_wl = size(wl_grid)

    call logger%info(fname, "Re-gridding spectroscopy for gas: " &
         // trim(gas%name%chars()))

    ! Fairly simple stuff here, just loop through every non-wavelength
    ! dimension of the spectroscopy table and stick the interpolated
    ! values into tmp_array.

    allocate(tmp_array(N_wl, &
         size(gas%cross_section, 2), &
         size(gas%cross_section, 3), &
         size(gas%cross_section, 4)))

    tmp_array(:,:,:,:) = 0.0d0

    ! Wavelength loop
    do d1=1, size(tmp_array, 1)

       ! Find left wavelength index in spectroscopy grid
       wl_idx_left = searchsorted_dp(gas%wavelength, wl_grid(d1))

       if ((wl_idx_left < 1) .or. (wl_idx_left > size(gas%wavelength))) then
          ! If the hires WL grid is requested beyond the limits of spectroscopy
          ! we have to fill it up with zeros.
          tmp_array(d1,:,:,:) = 0.0d0

          write(tmp_str, '(A, F10.7, A, F10.7, A, F10.7, A)') "Spectroscopy data is valid between: ", &
               gas%wavelength(1), " and ", gas%wavelength(N_wl), &
               " but you requested: ", wl_grid(d1), ". Filling with 0."
          call logger%debug(fname, trim(tmp_str))
          
          cycle
       endif
       
       fac = (wl_grid(d1) - gas%wavelength(wl_idx_left)) &
            / (gas%wavelength(wl_idx_left+1) - gas%wavelength(wl_idx_left))

       do d2=1, size(tmp_array, 2)
          do d3=1, size(tmp_array, 3)
             do d4=1, size(tmp_array, 4)

                tmp_array(d1,d2,d3,d4) = (1.0d0 - fac) * gas%cross_section(wl_idx_left, d2, d3, d4) &
                     + fac * gas%cross_section(wl_idx_left+1, d2, d3, d4)

             end do
          end do
       end do
    end do
    ! Now we have to replace the old spectroscopy array as well as the
    ! wavelength coordinate.

    deallocate(gas%cross_section)
    allocate(gas%cross_section, mold=tmp_array)
    gas%cross_section(:,:,:,:) = tmp_array(:,:,:,:)

    deallocate(gas%wavelength)
    allocate(gas%wavelength, mold=wl_grid)
    gas%wavelength(:) = wl_grid(:)

  end subroutine regrid_spectroscopy


end module spectroscopy_utils_mod
