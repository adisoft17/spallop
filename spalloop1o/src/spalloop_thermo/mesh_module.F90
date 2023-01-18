MODULE MESH_MODULE

    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutines for setting dimensions and meshes. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  real, save, public     :: rr_metal_in_initial
  real, save, public     :: total_oxide_thickness_ref

  ! Private Module
  private

  ! Public Procedures
  public :: GET_RADIUS, GET_MESH, GET_TUBE_DIMENSIONS, MESH_POINTS

  CONTAINS

  SUBROUTINE GET_TUBE_DIMENSIONS(id_input)
  
    !=======================================================================
    ! Purpose(s):
    ! set tube_outer_radius and tube_thickness such that minimal changes are required
    !  
    !=======================================================================

  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & oxide_location
  use waterwall_data_module,    only: radius_outer_tube, &
       & thickness_tube, id_length_tube_interval, no_height_oxide_output
  
  integer, intent(IN)  :: id_input

  ! local variables
  integer              :: tube_interval_id, k

  ! id_input is k = 1, no_height_oxide_output

  ! tube dimensions would be those within the tube_interval_id interval for length_tube
  if (id_input < no_height_oxide_output .and. id_input > 0)  then
    tube_interval_id = id_length_tube_interval(id_input)
  else if (id_input == no_height_oxide_output)  then
    tube_interval_id = id_length_tube_interval(id_input-1)
  else if (id_input <= 0)  then
    write(6, *)  'ERROR: interval is 0!! ', id_input, no_height_oxide_output
    STOP
  else
    write(6, *)  'ERROR: interval not valid ', id_input, no_height_oxide_output
    STOP
  endif
  
  tube_outer_radius = radius_outer_tube(tube_interval_id)
  tube_thickness = thickness_tube(tube_interval_id)

  if (tube_outer_radius <= 0.0001)  then ! units are m
    write(6, *)  'ERROR: tube_outer_radius = 0 ', id_input, tube_interval_id
    STOP
  endif
 
  if (tube_thickness <= 0.00001)  then ! units are m
    write(6, *)  'ERROR: tube_thickness = 0 ', id_input, tube_interval_id
    STOP
  endif  

  END   SUBROUTINE GET_TUBE_DIMENSIONS

  SUBROUTINE GET_RADIUS(N, thickness_oxide_layer, &
       & time_now, if_ave_scale, iter)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & oxide_location
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_int, rad_mean
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value, no_metal_layers
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide, update_sola_type
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  integer, intent(IN)  :: N, iter
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  real, intent(IN)      :: time_now
  logical, intent(IN)  :: if_ave_scale
  character(LEN = 80)  :: char1  

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0
  real     :: tube_outer_radius_now, tube_thickness_now, &
         & first_spacing, ratio_metal_layer
  real, dimension(moxide)  :: thickness_oxide_layer_now

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! initialize

  if (N + 1 > moxide)  then
    write(tty_lun, *) 'ERROR: radius cannot be computed '
    stop
  endif

  ! initialize the initial inner radius of the tube; reference radius
  rr_metal_in_initial = tube_outer_radius - tube_thickness ! * &
      ! & (1.0 + tau_th_st(imetal, 1)*(1.0 + poisson_ratio(imetal)))

  ! get position of metal-oxide interface
  call GET_METAL_DIM(N, tube_thickness_now, &
           & thickness_oxide_layer, time_now, iter)

  ! tube_outer_radius - is in fact tube radius at RT
  ! tube_thickness - is in fact tube thickness at RT

      tube_outer_radius_now = tube_outer_radius   ! r_o
      thickness_oxide_layer_now = thickness_oxide_layer

    if (TRIM(oxide_location) == 'inner_metal_surface')  then
    
      rad_int(1) = tube_outer_radius_now

    else if (TRIM(oxide_location) == 'outer_metal_surface')  then

      rad_int(1) = - tube_outer_radius_now

    else

      write(6, *) 'oxide_location must be inner_metal_surface or ', &
        & ' outer_metal_surface'
      write(6, *) 'oxide_location =', TRIM(oxide_location)
      STOP

    endif
  
    ! valid only for no_metal_layers = 1
    ! rad_int(2) = rad_int(1) - tube_thickness_now

    ratio_metal_layer = 0.8

    if (no_metal_layers == 1)  then

      do i = 1, no_metal_layers
        rad_int(i+1) = rad_int(i) - tube_thickness_now / no_metal_layers
      end do

    else if (no_metal_layers > 1)  then

      if (ABS(ratio_metal_layer-1.0) <= 1.0e-4)  then

        first_spacing = tube_thickness_now / no_metal_layers
        do i = 1, no_metal_layers
          rad_int(i+1) = rad_int(i) - tube_thickness_now / no_metal_layers
        end do

      else

        first_spacing = tube_thickness_now *(1.0 /ratio_metal_layer-1.0) / &
           & (1.0/ratio_metal_layer**no_metal_layers-1.0)
        do i = 1, no_metal_layers
          rad_int(i+1) = rad_int(i) -first_spacing/ratio_metal_layer**(i-1)
        end do

      endif

    endif

    do i = 1, no_metal_layers
 
      rad_mean(i) = 0.5 * (rad_int(i+1) + rad_int(i))

    end do

    ! deal with the scale that got average properties
    do i = no_metal_layers + 2, N+1

      if (if_ave_scale)  then

        ! in this case N = 2 (tube and one average oxide)
        rad_int(i) = rad_int(i-1) - thickness_oxide_layer(no_oxide_ave)  ! only for average

        if (ip <= 6)  then
          write(aux_lun, 7) i, thickness_oxide_layer(no_oxide_ave), &
             & rad_int(i-1), rad_int(i), &
             & thickness_oxide_layer(no_oxide_ave)
        endif

      else

        ! this will work in general; thickness_oxide_layer(1) is zero by default (metal)
        rad_int(i) = rad_int(i-1) - thickness_oxide_layer_now(i-1)  ! oxide layer i-2

        if (ip <= 6)  then

          ! fix this no_oxide_ave later on; if_ave_scale=FALSE no_oxide_ave=0
          no_oxide_ave = no_layer + 1  
          write(aux_lun, 7) i, thickness_oxide_layer_now(i-1), &
             & rad_int(i-1), rad_int(i), &
             & thickness_oxide_layer(no_oxide_ave)
        ! write(aux_lun, 7) i,  (thickness_oxide_layer(j), j=1, moxide)
 7        format('mat_interf ', i2, 1x, 20(1pe13.6, 1x))

        endif

      endif

    end do

    do i = no_metal_layers + 1, N
 
      rad_mean(i) = 0.5 * (rad_int(i+1) + rad_int(i))

    end do

    if (iter > 0 .and. time_now > 1.0e-3 .and. ip < 2000)  then

      write(aux_lun, 8) time_now, rad_int(1) - rad_int(N+1), &
         & (rad_int(i), i=1, N+1)
 8    format('rad_interf ', 20(1pe13.6, 1x))

    endif

  return

  END SUBROUTINE GET_RADIUS

  SUBROUTINE GET_METAL_DIM(N, tube_thickness_now, thickness_oxide_layer, &
      & time_now, iter)

    !=======================================================================
    ! Purpose(s):
    !
    !   get metal dimensions
    !  
    !=======================================================================

  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun
  use oxide_data_module,   only: pilling_bedworth_ratio, thickness_fr, &
      & first_oxide_layer, no_metal_layers
  use boiler_data_module, only: tube_outer_radius, tube_thickness
  use parameter_module,   only: moxide
  use solver_data_module, only  : metal_recession_type

  ! Argument List
  integer,    intent(IN)  :: N, iter
  real,    intent(OUT)  :: tube_thickness_now
  real, intent(IN)      :: time_now
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer

  ! Local Variables
  integer  :: i, k, j, mn, kp
  real, dimension(N+1)      :: psi
  real     :: b_term, c_term, r_metal_oxide, u_jump_fr, &
           & total_oxide_thickness_ref_si, r_approx

  if (TRIM(metal_recession_type) == 'constant')  then
    ! unchanged geometry
    tube_thickness_now = tube_thickness
    RETURN
  endif

  ! total oxide thickness in SI units
  total_oxide_thickness_ref_si = &
       & SUM(thickness_oxide_layer(first_oxide_layer:N))

  if (no_metal_layers > 1)  then
    ! consider that metal is one layer since no details is required here
    mn = N - no_metal_layers + 1
  else if (no_metal_layers == 1)  then
    mn = N
  endif

  ! get psi for each layer
  psi(1:mn+1) = 0.0
  psi(2) = thickness_fr(2) / 2.0
  do k = 3, mn   ! N layers including the metal
    kp = k -2 + first_oxide_layer
    psi(k) = psi(k-1) + 0.5 * (thickness_fr(kp-1) + thickness_fr(kp))
  end do

  ! get ri term and free term
  b_term = 0.0
  c_term = 0.0
  do k = 2, mn
    kp = k -2 + first_oxide_layer
    b_term = b_term + thickness_fr(kp) / pilling_bedworth_ratio ! (k)
    c_term = c_term + psi(k) * thickness_fr(kp) / pilling_bedworth_ratio ! (k)
  end do

  b_term = b_term * total_oxide_thickness_ref_si
  c_term = c_term * 2.0 * total_oxide_thickness_ref_si**2 - &
         & rr_metal_in_initial**2

  ! solve the quadratic equation for the inner metal radius (metal_oxide interf)
  ! r^2 - 2 * b * r + c = 0
  r_metal_oxide = b_term + sqrt(b_term**2 - c_term)  ! this should be positive
  ! r_sol(2) = b_term - sqrt(b_term**2 - c_term), this is negative
  
  ! small oxide approximation
  r_approx = 0.5 * (rr_metal_in_initial**2 - c_term) / &
           & (rr_metal_in_initial - b_term)

  if (iter < 5)  then

    write(out_lun, 2) iter, b_term, c_term, &
        & r_approx, b_term - sqrt(b_term**2 - c_term), &
        & b_term + sqrt(b_term**2 - c_term)
2 format('check_met_recess ', i4, 1x, 20(1pe13.6, 1x))

  endif

  tube_thickness_now = tube_outer_radius - r_metal_oxide

  if (iter > 0)  then
    ! displacement of metal_oxide interface
    write(out_lun, 1) iter, time_now, r_metal_oxide - &
      & (tube_outer_radius - tube_thickness), &
      & tube_thickness_now, tube_thickness, &
      & r_metal_oxide, r_approx, tube_outer_radius, &
      & total_oxide_thickness_ref_si, rr_metal_in_initial
1   format('met_ox_displ ', i4, 1x, 20(1pe13.6, 1x))
  endif

  return

  END SUBROUTINE GET_METAL_DIM

  SUBROUTINE GET_MESH(N, thickness_oxide_layer, if_ave_scale)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !    id_high_or_low = 1 for high temp
    !    id_high_or_low = 2 for low temp
    !  
    !=======================================================================

  use parameter_module,   only: moxide, mave
  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_st_old, rad_int, n_mean_even_vm, n_mean_odd_vm
  use oxide_data_module,   only: no_oxide, no_layer, &
      & no_oxide_ave, cond_value
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,            only: tty_lun, out_lun, inp_lun, aux_lun

  integer, intent(IN)  :: N
  real, dimension(moxide), intent(IN)  :: thickness_oxide_layer
  ! real, intent(IN)      :: time_now
  logical, intent(IN)  :: if_ave_scale
  character(LEN = 80)  :: char1  

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  do i = 1, N

    ! get the radius first
    do j = 1, npr_temp(i)
      rad_temp(i, j) = rad_int(i) - (rad_int(i) - rad_int(i+1)) * &
         & (j - 1) / (npr_temp(i) - 1)
    end do

    ! store the old radii
    rad_st_old(i, 1:moxide) = rad_st(i, 1:moxide)

    do j = 1, npr_st(i)
      rad_st(i, j) = rad_int(i) - (rad_int(i) - rad_int(i+1)) * &
         & (j - 1) / (npr_st(i) - 1)
    end do

    if (ip < 6)  then

      do j = 1, npr_temp(i)
          write(aux_lun, 7) i, j, rad_temp(i, j)
 7        format('rad_temp_out ', 2(i2, 1x), 20(1pe13.6, 1x))

      end do

    endif

  end do

  return

  END SUBROUTINE GET_MESH

  SUBROUTINE MESH_POINTS(N)
  
    !=======================================================================
    ! Purpose(s):
    !
    !   Obtain the temperature profile
    !  
    !=======================================================================

  use solution_data_module, only: npr_temp, npr_st, n_mean_even_vm, &
      & n_mean_odd_vm
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,      only: tty_lun, out_lun, inp_lun, aux_lun

  integer, intent(IN)  :: N

  ! Local Variables
  integer  :: i, j, loc1, k
  integer, save  :: ip = 0

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  do i = 1, N

    ! get number of points in the radial direction
    if (i == 1)  then
      npr_temp(i) = no_temp_rad_points_tube
      npr_st(i) = no_st_rad_points_tube
    else 
      npr_temp(i) = no_temp_rad_points_oxide
      npr_st(i) = no_st_rad_points_oxide
    endif
    
    n_mean_even_vm(i) = INT(npr_st(i)/2)
    n_mean_odd_vm(i) = INT((npr_st(i)+1)/2)

    if (npr_st(i) == 2 * n_mean_even_vm(i))  then
      ! npr_st is an even number
      ! n_mean1 = INT(npr_st(i)/2)
      ! integral_th_mean(i) = (integral_th(i, n_mean) + &
      !        & integral_th(i, n_mean+1)) / 2.0
    else
      ! n_mean1 = INT((npr_st(i)+1)/2)
      if (npr_st(i) + 1 == 2 * n_mean_odd_vm(i))  then
        ! npr_st is odd
        ! integral_th_mean(i) = integral_th(i, n_mean1)
      endif
    endif

  end do

  return

  END SUBROUTINE MESH_POINTS

END MODULE MESH_MODULE
