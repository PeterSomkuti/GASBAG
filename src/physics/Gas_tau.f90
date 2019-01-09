module gas_tau_mod

  use control_mod, only: CS_gas
  use math_utils_mod

  public calculate_gas_tau, get_CS_value_at

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
    double precision :: CS_value_grid(size(wl))


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

    ! Traverse the atmosphere layers(!), starting from the bottom to the top
    do l=N_lev,2,-1

       ! First, grab the lower (altitude-wise) and higher values for
       ! p,T,H2O from the atmosphere profiles for whatever layer l we're in.

       if (l == N_lev) then
          ! BOA layer - psurf should be between lowermost level and
          ! the level above (ideally). But this will extrapolate linearly
          ! anyway, which might cause some trouble.
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

          write(*,* ) "am in sublayer ", k, "for layer", l

          ! Get the value of pressure according to the Gauss-Kronrod abscissa
          ! value, at which we evaluate the cross section.
          this_p = GK_abscissae(k) * (p_lower - p_higher) / 2.0d0 &
               + (p_higher + p_lower) / 2.0d0

          ! And get corresponding values for T,sh and the VMR at this (scaled)
          ! pressure value, which is obtained via simple linear interpolation.
          this_p_fac = ((this_p - p_higher) / (p_lower - p_higher))
          this_T = (1.0d0 - this_p_fac) * T_lower + this_p_fac * T_higher
          this_H2O = (1.0d0 - this_p_fac) * H2O_lower + this_p_fac * H2O_higher
          this_VMR = (1.0d0 - this_p_fac) * VMR_lower + this_p_fac * VMR_higher
          this_M = 1d3 * (((1 - this_H2O) * dry_air_mass) + (H2Om * this_H2O))

          do j=1, size(wl)
             CS_value_grid(j) = get_CS_value_at(gas, wl(j), this_p, this_T, this_H2O)
             !gas_tau(j,l-1) = gas_tau(j,l-1) + (CS_value_grid(j) + this_vmr * (1.0d0 - this_H2O) &
             !     / (9.81 * this_M) * NA * GK_weights(k) * 0.1d0)
          end do



       end do

       !gas_tau(:,l-1) = gas_tau(:,l-1) * (p_lower - p_higher)

    end do




 end subroutine calculate_gas_tau



 pure function get_CS_value_at(gas, wl, p, T, H2O) result(CS_value)

   type(CS_gas), intent(in) :: gas
   double precision, intent(in) :: wl, p, T, H2O

   double precision :: CS_value

   ! These here are the indices between which the CS array will
   ! be linearly interpolated: e.g. idx_(left, right)_(pressure)
   integer :: idx_closest
   integer :: idx_l_p, idx_r_p ! Pressure
   integer :: idx_lr_T, idx_ll_T, idx_rl_T, idx_rr_T ! Temperature has two dimensions
   integer :: idx_l_H2O, idx_r_H2O
   integer :: idx_l_wl, idx_r_wl
   double precision :: N(16)
   double precision :: C3(0:1,0:1,0:1) ! 0 is 'left', 1 is 'right'
   double precision :: C2(0:1, 0:1), C1(0:1)
   double precision :: wl_d, p_d, T_d_l, T_d_r, H2O_d
   integer :: i,j


   ! Get the pressure indices
   idx_l_p = get_left_index(p, gas%p)
   idx_r_p = idx_l_p + 1

   ! Get the temperature indices, which depend on P - so we have two
   ! indices for the temperature dimension. One set of T indices (idx_l(l,r)_T)
   ! corresponds to the left index of pressure (idx_l_p), the other set
   ! (idx_r(l,r)_T) corresponds to the right index of pressure (idx_r_p).
   idx_ll_T = get_left_index(T, gas%T(:, idx_l_p))
   idx_lr_T = idx_ll_T + 1

   idx_rl_T = get_left_index(T, gas%T(:, idx_r_p))
   idx_rr_T = idx_rl_T + 1

   ! Get the water vapor indices
   if (gas%has_H2O) then
      idx_l_H2O = get_left_index(H2O, gas%H2O)
      idx_r_H2O = idx_l_H2O + 1
   else
      idx_l_H2O = 1
      idx_r_H2O = 1
   end if

   ! Get the wavelength indices
   idx_l_wl = get_left_index(wl, gas%wavelength)

   ! And perform the linear interpolation in now 3 or 4 dimensions!
   ! Remember, we store the CS grid in the following way:
   ! wavelength, h2o, temperature, pressure
   ! .. so we order the vertices of our polytope in the same way.

   wl_d = (wl - gas%wavelength(idx_l_wl)) / &
        (gas%wavelength(idx_r_wl) - gas%wavelength(idx_l_wl))
   p_d = (p - gas%p(idx_l_p)) / (gas%p(idx_r_p) - gas%p(idx_l_p))
   if (gas%has_H2O) then
      H2O_d = (H2O - gas%H2O(idx_l_H2O)) / (gas%H2O(idx_r_H2O) - gas%H2O(idx_l_H2O))
   else
      H2O_d = 0.0d0
   end if
   T_d_l = (T - gas%T(idx_ll_T, idx_l_p)) / &
        (gas%T(idx_lr_T, idx_l_p) - gas%T(idx_ll_T, idx_l_p))
   T_d_r = (T - gas%T(idx_lr_T, idx_r_p)) / &
        (gas%T(idx_rr_T, idx_r_p) - gas%T(idx_lr_T, idx_r_p))

   ! Start interpolating along wavelength, meaning that C3 = C3(h2o, temperature, pressure)
   C3(0,0,0) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_ll_T, idx_l_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_l_H2O, idx_ll_T, idx_l_p) * wl_d

   C3(1,0,0) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_ll_T, idx_l_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_r_H2O, idx_ll_T, idx_l_p) * wl_d

   C3(0,1,0) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_lr_T, idx_l_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_l_H2O, idx_lr_T, idx_l_p) * wl_d

   C3(0,0,1) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_rl_T, idx_l_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_l_H2O, idx_rl_T, idx_r_p) * wl_d

   C3(1,1,0) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_lr_T, idx_l_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_l_H2O, idx_lr_T, idx_l_p) * wl_d

   C3(0,1,1) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_rr_T, idx_r_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_l_H2O, idx_rr_T, idx_r_p) * wl_d

   C3(1,0,1) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_rl_T, idx_r_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_r_H2O, idx_rl_T, idx_r_p) * wl_d

   C3(1,1,1) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_rr_T, idx_r_p) * (1.0d0 - wl_d) &
        + gas%cross_section(idx_r_wl, idx_r_H2O, idx_rr_T, idx_r_p) * wl_d

   ! Next, go across the H2O dimension, such that C2 = C2(temperature, pressure). Since we do not
   ! really need the cross section indices here, we can convientiently write this section as a loop

   do i=0,1
      do j = 0,1
         C2(i,j) = C3(0,i,j) * (1.0d0 - H2O_d) + C3(1,i,j) * H2O_d
      end do
   end do

   ! Next, go across the temperature dimension
   do i=0,1
      C1(i) = C2(0,i) * (1.0d0 - T_d_l) + C2(1,i) * T_d_r
   end do

   ! Finally, we interpolate along pressure, which is the final result that
   ! we pass back.
   CS_value = C1(0) * (1.0d0 - p_d) + C1(1) * p_d

 end function

 pure function get_left_index(x, array) result(left_index)

   double precision, intent(in) :: x, array(:)
   integer :: left_index

   integer :: idx_closest

   idx_closest = minloc(abs(x - array), dim=1)

   if ((x < array(1)) .or. (idx_closest == 1)) then
      left_index = 1
      return
   else if ((x > array(size(array))) .or. (idx_closest == size(array))) then
      left_index = size(array) - 1
      return
   else
      if (abs(x - array(idx_closest)) < abs(x - array(idx_closest - 1))) then
         left_index = idx_closest
      else
         left_index = idx_closest - 1
      end if
   end if

 end function get_left_index



end module gas_tau_mod
