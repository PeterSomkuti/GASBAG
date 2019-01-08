module gas_tau_mod

  use control_mod, only: CS_gas
  use math_utils_mod

  public calculate_gas_tau

contains

  subroutine calculate_gas_tau(wl, gas_vmr, psurf, p, T, H2O, gas, N_sublayer, gas_tau)

    implicit none
    double precision, intent(in) :: wl(:) ! Wavelength array
    double precision, intent(in) :: gas_vmr(:) ! Gas volume mixing ratio
    double precision, intent(in) :: psurf ! Surface pressure
    double precision, intent(in) :: p(:), T(:), H2O(:) ! Pressure, temperature and water vapor profiles
    type(CS_gas), intent(in) :: gas ! Gas structure which contains the cross sections
    integer, intent(in) :: N_sublayer ! How many sublayers for more accurate calculation?

    double precision, intent(inout) :: gas_tau(:,:) ! Gas optical depth - wavelength, layer

    double precision :: p_lower, p_higher, T_lower, T_higher, H2O_lower, H2O_higher
    double precision :: VMR_lower, VMR_higher
    double precision :: this_p, p_fac, this_p_fac, this_T, this_H2O, this_VMR, this_M

    double precision :: GK_abscissae(N_sublayer+1), GK_weights(N_sublayer+1), G_weights(N_sublayer+1)

    integer :: N_lay, N_lev, N_wl
    integer :: i,j,k,l

    N_lev = size(gas_vmr)
    N_lay = N_lev - 1
    N_wl = size(wl)

    ! Given the number of sublayers, we precompute the Gauss-Kronrod weights for
    ! the integral approximation. Only if N_sublayer > 0 obviously
    if (N_sublayer > 0) then
       call kronrod(N_sublayer, 1d-6, GK_abscissae, GK_weights, G_weights)
    end if

    ! Just to make sure there's nothing in there already..
    gas_tau(:,:) = 0.0d0

    ! Traverse the atmosphere, starting from the bottom to the top
    do l=N_lev,2,-1

       ! First, grab the lower (altitude-wise) and higher values for
       ! p,T,H2O from the atmosphere profiles.

       if (l == N_lev) then
          p_lower = psurf
          p_fac = (p_lower - p(l-1)) / (p(l) - p(l-1))
          T_lower = (1.0d0 - p_fac) * T(l-1) + p_fac * T(l)
          H2O_lower = (1.0d0 - p_fac) * H2O(l-1) + p_fac * H2O(l)
          VMR_lower = (1.0d0 - p_fac) * gas_vmr(l-1) + p_fac * gas_vmr(l)
       else
          p_lower = p(l)
          T_lower = T(l)
          H2O_lower = H2O(l)
          VMR_lower = gas_vmr(l)
       end if

       p_higher = p(l-1)
       T_higher = T(l-1)
       H2O_higher = H2O(l-1)
       VMR_higher = gas_vmr(l-1)

       do k=1, N_sublayer+1

          ! Re-scale the value of pressure, as we are doing a GK-integration
          ! that takes place from -1 to 1.
          this_p = GK_abscissae(k) * (p_lower - p_higher) / 2.0d0 &
               + (p_higher + p_lower) / 2.0d0

          ! And get corresponding values for T,sh and the VMR
          this_p_fac = ((this_p - p_higher) / (p_lower - p_higher))
          this_T = (1.0d0 - this_p_fac) * T_lower + this_p_fac * T_higher
          this_H2O = (1.0d0 - this_p_fac) * H2O_lower + this_p_fac * H2O_higher
          this_VMR = (1.0d0 - this_p_fac) * VMR_lower + this_p_fac * VMR_higher
          this_M = 1d3 * (((1 - this_H2O) * dry_air_mass) + (H2Om * this_H2O))




       end do



    end do




 end subroutine calculate_gas_tau





end module gas_tau_mod
