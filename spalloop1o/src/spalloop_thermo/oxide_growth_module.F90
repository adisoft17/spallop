MODULE OXIDE_GROWTH_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   subroutines for oxide growth. 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: OXIDE_GROWTH_LOCATION, OXIDE_GROWTH_INIT

  CONTAINS

  SUBROUTINE OXIDE_GROWTH_INIT()

    !=======================================================================
    ! Purpose(s):
    !
    ! Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, oxide_thick2si, id_high, id_low, &
          & pi
    use output_module,      only: tty_lun, out_lun, aux_lun, enrg_lun
    use boiler_data_module, only: Time_Boiler_p, Temp_gas_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, press_inner_p, &
          & press_inner_slope_p, oxide_thickness, id_temp_op, &
          & thickness_fr_var, press_outer_p, press_outer_slope_p, &
          & tube_thickness, tube_outer_radius
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, const_energy_oxide, &
          & if_temp_steam_or_oxide, thickness_fr, no_oxide_ave, &
          & growth_ox_temperature, thickness_growth_fr, first_oxide_layer, &
          & n_it_max, error_temp_max, error_temp_norm,  &
          & error_oxide_thick_ratio
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use mesh_module, only: GET_RADIUS, GET_TUBE_DIMENSIONS
    use solution_data_module, only: Temp, npr_temp, no_oxide_layers, &
          & oxide_thickness_hf_out, time_hf_out, no_full2low_output, &
          & no_point_full2low, no_pulse_full2low_output, &
          & id_full_load_hf, n_hf_out, rad_int, &
          & time_outage, id_time_outage, no_outages
    use htc_utility_module, only: GET_TEMP_HTC_INNER, GET_ASH_THICKNESS, &
          & LOAD_FUNCTION, GET_TEMP_METAL, GET_STEAM_TEMP_HTC_INNER, &
          & OUT_CONDUCT_RESISTANCE
    use waterwall_data_module, only: no_height_oxide_output, nd1_ash_max, &
          & Temp_growth_layer_height_old, flow_rate_steam_high, &
          & flow_rate_steam_low, Temp_steam_height, Temp_gas_height, &
          & heat_flux_output, Temp_growth_layer_height_old, &
          & delta_oxide_thickness_layer_old, oxide_thickness_layer_height, &
          & no_thick_interval_now, time_now, length_tube_oxide_output, &
          & delta_oxide_thickness_layer_new, Temp_growth_layer_height_new, &
          & htc_solve_steam_or_wall, if_debug_oxide_growth_conv, &
          & no_time_intervals_ash, if_output_height_time_all, &
          & hflux_htc_ash_estimate, heat_flux_output_ash, &
          & htc_tube_outer_ash
    use waterwall_data_module, only: thickness_oxide_height_time
    use blockage_data_module,  only: dthick_ox

    implicit none

    ! Argument List

    ! Local Variables
    integer :: i, j, k, Status, n1, n2, n21, ik, ilay, &
             & ihf, it_dummy = 1, it_dummy10 = 10, k_load_start, k_load_end
    real    :: delta_t, temp1w, temp2w, &
             & one_over_temp, temp_ave, time2
    real    :: time_dummy = 0.0, oxide_thickness_old

    real, dimension(moxide) :: thickness_layer
    real, dimension(3) :: Temp_metal
    integer  :: convergent, j_output, id_ash

    logical :: if_output_event, if_time_output_event
    logical, save :: if_first_ox = .true.

    if (if_first_ox)  then

      time_outage = 0.0
      id_time_outage = 0
      k = 0
      id_time_outage(k) = 1
      time_outage(k) = Time_Boiler_p(1)

      ! identify all the outage times
      do i= 2, TEMP_BOILER_data_no-1

        do j = 1, no_full2low_output
          if (i == no_point_full2low(j) .and. &
            & no_pulse_full2low_output(j) == 0)  then
            k = k + 1
            time_outage(k) = Time_Boiler_p(i)
            id_time_outage(k) = i
            write(6, *) 'outage ', k, i, id_temp_op(i), id_time_outage(k), time_outage(k)
          endif
        end do

      end do
      no_outages = k

      if (no_outages > 0)  then

        ! ALLOCATE(dthick_ox_height_time(0:no_outages, &
        !  & 1:no_height_oxide_output, 1:TEMP_BOILER_data_no), STAT = Status)
        ! dthick_ox_height_time = 0.0

        do n1 = 1, no_outages
          do i= 1 + id_time_outage(n1-1), id_time_outage(n1)
            write(out_lun, 20) n1, Time_Boiler_p(i-1), &
             & Time_Boiler_p(i)
 20          format('check_outage_2_follow ', 1(i2, 1x), 30(1pe13.6, 1x))
          end do
        end do

      else
        write(6, *) 'ERROR: no outages'
        STOP
      endif

      ALLOCATE(thickness_oxide_height_time(first_oxide_layer:no_oxide_layers, &
          & 1:no_height_oxide_output, 1:TEMP_BOILER_data_no), STAT = Status)

      if_first_ox = .false.
   
    endif

    ! units are in microns
    ALLOCATE(oxide_thickness(1:TEMP_BOILER_data_no), STAT = Status)

    ! for metal only do not obtain the oxide thickness
    if (no_oxide_layers == 1)  then

      oxide_thickness = 0.0
      ALLOCATE(thickness_fr_var(1:TEMP_BOILER_data_no, 2:2), &
             & STAT = Status)
      thickness_fr_var(:, 2) = 1.0

      RETURN
    endif

    if (SUM(thickness_growth_fr) < 1.0e-6)  then

      write(out_lun, *) 'ERROR: old run needs thickness_growth_fr', &
         & 'default growth_ox_temperature = oxide_surface'
      STOP

    endif

    ALLOCATE(thickness_fr_var(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)

    error_temp_max = 0.1
    error_temp_norm = 0.01
    error_oxide_thick_ratio = 1.0e-6

    n_it_max = 20

    ! if (ABS(rad_int(2)) < 1.0e-9)  then
    if (ABS(tube_outer_radius - tube_thickness) < 1.0e-9)  then
      write(6, *)  'ERROR: rad_inner_tube = 0 ', tube_outer_radius, tube_thickness
      stop
    endif

    ! if (ALL(
    do i = first_oxide_layer, no_oxide_layers

      if (TRIM(growth_ox_temperature(i)) == 'steam')  then
        ! oxide grows at steam temperature
        n_it_max = 1
      else
        n_it_max = 20
      endif

    end do

   dthick_ox = 0.0

   ! print the header only; formatt 9 goes with format 8
   write(2, 9) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
   write(enrg_lun, 9) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
 9        format('metal_temp_min ', i5, 1x, 30(1pe13.6, 1x))
   write(2, 11) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
   write(enrg_lun, 11) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
 11        format('metal_temp_ave ', i5, 1x, 30(1pe13.6, 1x))
   write(2, 13) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
   write(enrg_lun, 13) i-i, error_temp_max-error_temp_max, &
          & (length_tube_oxide_output(j), j = 1, no_height_oxide_output)
 13        format('metal_temp_max ', i5, 1x, 30(1pe13.6, 1x))
      
    return

  END SUBROUTINE OXIDE_GROWTH_INIT

  SUBROUTINE OXIDE_GROWTH_LOCATION(i_time_start, i_time_end, &
      &   if_outage_ramp_exfol, id_previous_outage, id_next_outage)

    !=======================================================================
    ! Purpose(s):
    !
    ! Compute the oxide thickness increment within this current time interval
    !  based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: moxide, oxide_thick2si, id_high, id_low, &
          & pi
    use output_module,      only: tty_lun, out_lun, aux_lun, enrg_lun
    use boiler_data_module, only: Time_Boiler_p, Temp_gas_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, press_inner_p, &
          & press_inner_slope_p, oxide_thickness, id_temp_op, &
          & thickness_fr_var, press_outer_p, press_outer_slope_p, &
          & tube_thickness, tube_outer_radius
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, const_energy_oxide, &
          & if_temp_steam_or_oxide, thickness_fr, no_oxide_ave, &
          & growth_ox_temperature, thickness_growth_fr, first_oxide_layer, &
          & n_it_max, error_temp_max, error_temp_norm,  &
          & error_oxide_thick_ratio
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use mesh_module, only: GET_RADIUS, GET_TUBE_DIMENSIONS
    use solution_data_module, only: Temp, npr_temp, no_oxide_layers, &
          & oxide_thickness_hf_out, time_hf_out, no_full2low_output, &
          & no_point_full2low, no_pulse_full2low_output, &
          & id_full_load_hf, n_hf_out, rad_int, &
          & time_outage, id_time_outage, no_outages
    use htc_utility_module, only: GET_TEMP_HTC_INNER, GET_ASH_THICKNESS, &
          & LOAD_FUNCTION, GET_TEMP_METAL, GET_STEAM_TEMP_HTC_INNER, &
          & OUT_CONDUCT_RESISTANCE
    use waterwall_data_module, only: no_height_oxide_output, nd1_ash_max, &
          & Temp_growth_layer_height_old, flow_rate_steam_high, &
          & flow_rate_steam_low, Temp_steam_height, Temp_gas_height, &
          & heat_flux_output, Temp_growth_layer_height_old, &
          & delta_oxide_thickness_layer_old, oxide_thickness_layer_height, &
          & no_thick_interval_now, time_now, length_tube_oxide_output, &
          & delta_oxide_thickness_layer_new, Temp_growth_layer_height_new, &
          & htc_solve_steam_or_wall, if_debug_oxide_growth_conv, &
          & no_time_intervals_ash, if_output_height_time_all, &
          & hflux_htc_ash_estimate, heat_flux_output_ash, &
          & htc_tube_outer_ash, temp_steam_height_time, &
          & id_height_now
    use waterwall_data_module, only: thickness_oxide_height_time
    use blockage_data_module,  only: dthick_ox, &
          & thickness_new_ox_growth_init, if_exfoliation_oxide_growth, &
          & thickness_ave_ox_growth_init, start_exfoliation_layers 
    use oxide_stress_data_module,   only: ispinel, imagnetite
    use steam_temp_module,         only: GET_HEAT_FLUX_ST_TEMP_OX, &
          & TEMP_OUT_STEAM_OXIDE_METAL

    implicit none

    ! Argument List
    integer, intent(IN) :: i_time_start, i_time_end, id_previous_outage, &
          & id_next_outage
    logical, intent(INOUT) :: if_outage_ramp_exfol

    ! Local Variables
    integer :: i, j, k, Status, n1, n2, n21, ik, ilay, i_time, j_height, &
             & ihf, it_dummy = 1, it_dummy10 = 10, k_load_start, k_load_end
    integer, dimension(3) :: iex_out
    real    :: delta_t, temp1w, temp2w, &
             & one_over_temp, temp_ave, time2
    real    :: time_dummy = 0.0, oxide_thickness_old
    real, dimension(2:no_oxide_layers)    :: kp, kp_log, temp2, temp1,  &
             & temp_slope, error_temp, dtemp, Temp_growth_now, dz2_old, &
             & Temp_growth_old, dz2, Temp_growth_end, Temp_growth_start, &
             & dz2_skin

    real, dimension(moxide) :: thickness_layer
    real, dimension(3) :: Temp_metal
    real, pointer, dimension(:,:) :: Temp_metal_min, Temp_metal_mean, &
             & Temp_metal_max
    integer  :: convergent, j_output, id_ash

    real :: flow_rate_steam_slope, flow_rate_steam_end, height_now, &
          & flow_rate_steam_start, Temp_steam_start, Temp_steam_end, &
          & Temp_steam_slope, Temp_gas_start, Temp_gas_end, Temp_gas_slope, &
          & heat_flux_start, heat_flux_end, heat_flux_slope, Temp_steam_now, &
          & Temp_gas_now, heat_flux, press_inner_now, press_outer_now, &
          & flow_rate_steam_now, htc_inner_new, Temp_wall_now, &
          & ash_thickness_now, oxide_thickness_now, &
          & delta_oxide_thickness_now, delta_oxide_thickness_old, &
          & temp_growth_error_max, htc_outer_now, dz2_increment
    real :: temp_growth_error, htc_inner_start, &
          & Temp_wall_start, Temp_wall_end, htc_inner_end, time_mean, &
          & flow_rate_steam_now_m2, oxide_thickness_not_exfol, &
          & oxide_thickness_exfol
    real, dimension(4, no_full2low_output, no_height_oxide_output) :: Temp_steam_out, &
          & htc_inner_out, Temp_metal_mean_out, Temp_metal_max_out, &
          & heat_flux_out, Temp_gas_out, flow_rate_out
    real, dimension(no_height_oxide_output) :: oxide_thickness_out
    real, dimension(20)   :: time_j_output

    logical :: if_output_event, if_time_output_event
    logical, save :: if_first_dox = .true.
  character(LEN = 1), dimension(0:9)   :: ch_index = (/'0','1','2','3','4','5', &
                  & '6','7','8','9'/)

    if (if_first_dox)  then

      ALLOCATE(Temp_metal_min(1:no_height_oxide_output, &
        & nd1_ash_max), STAT = Status)
      ALLOCATE(Temp_metal_mean(1:no_height_oxide_output, &
        & nd1_ash_max), STAT = Status)
      ALLOCATE(Temp_metal_max(1:no_height_oxide_output, &
        & nd1_ash_max), STAT = Status)

  ! allocate steam temperature
  ALLOCATE(TEMP_STEAM_HEIGHT_TIME(1:no_height_oxide_output, &
      & 0: id_time_outage(no_outages), 1:2), STAT = Status)

      if_first_dox = .false.

    endif

      if (n_it_max == 1)  then
        convergent = 1
      else
        convergent = 0  ! initialize convergence flag
      endif

  ! allocate steam temperature
  ! ALLOCATE(TEMP_STEAM_HEIGHT_TIME(1:no_height_oxide_output, &
  !    & i_time_start:i_time_end, 1:2), STAT = Status)

  ! get steam temperature at all the heights and all times between outages
  ! using the oxide thickness at previous time
  ! get Temp_steam_height_time
  ! for first outage i_time_start = 2
  ! this is needed to initialize the Tsteam(i_time -1)
  ! not sure that this is needed to recalculate after exfoliation
  call GET_HEAT_FLUX_ST_TEMP_OX(i_time_start - 1)

  TIME_LOOP_HOX1: do i= i_time_start, i_time_end  ! no of data intervals

    ! get the id time instead of the time loop
    i_time = i

    ! write(2, 9) i, j_height, Time_Boiler_p(i-1), Time_Boiler_p(i)
    ! write(6, 9) i, j_height, Time_Boiler_p(i-1), Time_Boiler_p(i)
 9  format('start_time_loop ', i5, 1x, i2, 1x, 10(1pe13.6, 1x))

    ! get the time at the start of the i-th interval
    time2 = Time_Boiler_p(i_time-1)

    if (id_temp_op(i_time) == 1)  then 
      ! interval 4-1, isothermal low load
      k_load_start = id_low
      k_load_end = id_low
    else if (id_temp_op(i_time) == 2)  then 
        ! interval 1-2, transition partial - full load
      k_load_start = id_low
      k_load_end = id_high
    else if (id_temp_op(i_time) == 3)  then 
        ! interval 2-3, isothermal high load
      k_load_start = id_high
      k_load_end = id_high
    else if (id_temp_op(i_time) == 4)  then 
        ! interval 3-4, transition full - partial load
      k_load_start = id_high
      k_load_end = id_low
    endif

    if_output_event = .false.  ! regular event
    if_time_output_event = .false. ! to save disk space height_event_output_ash

        if (i_time == i_time_end)  then  ! i_time_end should be id_time_outage(n1)
          ! baseline used should be id_previous_outage = id_next_outage - 1
          if_outage_ramp_exfol = .true.

          if (k == 1)  then
            write(out_lun, 20) id_next_outage, id_previous_outage, time2
 20          format('outage_2_follow ', 2(i2, 1x), 30(1pe13.6, 1x))
          endif

        endif

    ! check if this is in an output time
    do j = 1, no_full2low_output

      ! print out full three cycles as it looks better on figures
      if (i_time >= no_point_full2low(j) -1 .and. i <= no_point_full2low(j) + 7)  then
        ! this starts at [1:2]
        if (no_pulse_full2low_output(j) > 0)  then  ! not an outage
          if_time_output_event = .true. 
        endif

      endif

      ! 2-3 full load loop
      if (i_time == no_point_full2low(j))  then
 
        if (no_pulse_full2low_output(j) > 0)  then
          ! this is not an outage; 
          ! if id_temp_op(i) = 3 or k_load_end = id_high -> full load
          if_output_event = .true.
          j_output = j
          time_j_output(j_output) = Time_Boiler_p(i_time-1)
          ! write(6, *) 'output event ', id_temp_op(i-1), id_temp_op(i)
          ! confirmed that id_temp_op(i-1) = 2, id_temp_op(i) = 3

          ! initiate the new oxide thickness to the current one
          do k = 1, no_height_oxide_output
            oxide_thickness_out(k) = &
              & SUM(oxide_thickness_layer_height(first_oxide_layer: &
              & no_oxide_layers, k))
          end do

          write(enrg_lun, 17) j_output, Time_Boiler_p(i_time-1), &
              & (oxide_thickness_out(k), &
              & k = 1, no_height_oxide_output)
          write(out_lun, 17) j_output, Time_Boiler_p(i_time-1), &
              & (oxide_thickness_out(k), &
              & k = 1, no_height_oxide_output)
 17       format('output event ', i1, 1x, 50(1pe13.6, 1x))

        endif
      endif

      ! 4-1 low load loop; following the previous full load 2-3
      if (i_time -2 == no_point_full2low(j))  then
        ! if id_temp_op(i_time) = 1 or k_load_end = id_low -> full load
 
        if (no_pulse_full2low_output(j) > 0)  then
          ! this is not an outage 
          if_output_event = .true.
          j_output = j
        endif
      endif

    end do

    ! get steam temperature at all the heights and all times between outages
    ! using the oxide thickness at previous time
    ! get Temp_steam_height_time
    call GET_HEAT_FLUX_ST_TEMP_OX(i_time)

    ! set steam flow rate; does not depend on height  
    ! ACTION 1: flow_rate_steam quantities are defined
    flow_rate_steam_start = LOAD_FUNCTION(k_load_start, &
     & flow_rate_steam_high, flow_rate_steam_low)
    flow_rate_steam_end = LOAD_FUNCTION(k_load_end, &
       & flow_rate_steam_high, flow_rate_steam_low)
    flow_rate_steam_slope = (flow_rate_steam_end - &
       & flow_rate_steam_start) / &
       & (Time_Boiler_p(i_time) - Time_Boiler_p(i_time-1))

    ! initialize the growth temperature at the previous time level
    do j = first_oxide_layer, no_oxide_layers
      Temp_growth_layer_height_old(j, 1:no_height_oxide_output) = &
       & Temp_steam_height_time(1:no_height_oxide_output, i_time, k_load_end)
    end do

    ! initialize the thickness of the oxide that grows during this 
    ! micro cycle

    delta_oxide_thickness_layer_old(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output) = 0.0

    ! at each height determine the oxide thickness
    HEIGHT_LOOP1: do j = 1, no_height_oxide_output

      ! get the height instead of the the loop over all the heights
      j_height = j

      id_height_now = j      

      height_now = length_tube_oxide_output(j_height)

      ! set steam temperature at this height; begin, end, slope
      Temp_steam_start = Temp_steam_height_time(j_height, i_time - 1, k_load_start)
      Temp_steam_end = Temp_steam_height_time(j_height,i_time, k_load_end)
      Temp_steam_slope = (Temp_steam_end - Temp_steam_start) / &
          & (Time_Boiler_p(i_time) - Time_Boiler_p(i_time-1))

      ! set gas temperature at this height; begin, end, slope
      Temp_gas_start = Temp_gas_height(j_height, k_load_start)
      Temp_gas_end = Temp_gas_height(j_height, k_load_end)
      Temp_gas_slope = (Temp_gas_end - Temp_gas_start) / &
          & (Time_Boiler_p(i_time) - Time_Boiler_p(i_time-1))

      ! set heat flux begin, end, slope
      heat_flux_start = heat_flux_output(j_height, k_load_start)
      heat_flux_end = heat_flux_output(j_height, k_load_end)
      heat_flux_slope = (heat_flux_end - heat_flux_start) / &
          & (Time_Boiler_p(i_time) - Time_Boiler_p(i_time-1))

      ! get radii at the interfaces: 
      ! steam-metal, metal-oxide, oxide-oxide, oxide-steam
      ! assumption: 
      ! oxide thickness - for temperature computations does not change
      !     within the micro-time loop
      ! is redundant with the ones from stress_driver_module; 
      ! keep it here, too, such that oxide_growth can be separated, if needed later on
      thickness_layer(1:first_oxide_layer-1) = 0.0 ! metal layers

      if (.not. if_exfoliation_oxide_growth)  then

        thickness_layer(first_oxide_layer:no_oxide_layers) = oxide_thick2si * &
         & oxide_thickness_layer_height(first_oxide_layer:no_oxide_layers, j_height)

      else if (if_exfoliation_oxide_growth)  then

        thickness_layer(first_oxide_layer:start_exfoliation_layers - 1) = &
          & oxide_thick2si * oxide_thickness_layer_height(first_oxide_layer: &
             & start_exfoliation_layers - 1, j_height)

        ! oxide thickness changed after exfoliation
        ! exfoliation layers with the ave_ox baseline

        if (i_time == i_time_start)  then
          ! magnetite was exfoliated so cannot use the regular formula
          ! the oxide increment
          thickness_layer(imagnetite) = oxide_thick2si * &
             & thickness_ave_ox_growth_init(j_height, id_previous_outage, imagnetite)
        else
          thickness_layer(imagnetite) = oxide_thick2si * &
            & oxide_thickness_layer_height(imagnetite, j_height)
        endif

      endif

      ! set tube_outer_radius and tube_thickness at this location
      call GET_TUBE_DIMENSIONS(j_height)

      call GET_RADIUS(no_oxide_layers, thickness_layer, time2, &
        & if_average_oxide, it_dummy10)

      if (.not. if_exfoliation_oxide_growth)  then

        ! initiate the new oxide thickness to the current one
        oxide_thickness_old = SUM(oxide_thickness_layer_height(first_oxide_layer &
          & :no_oxide_layers, j_height))

        ! initialize to previous value of the oxide thickness
        dz2_old(first_oxide_layer:no_oxide_layers) = &
          & oxide_thickness_layer_height(first_oxide_layer &
          & :no_oxide_layers, j_height)**2
 
      else if (if_exfoliation_oxide_growth)  then

        ! initiate the new oxide thickness to the current one
        oxide_thickness_old = SUM(oxide_thickness_layer_height(first_oxide_layer &
          & :start_exfoliation_layers - 1, j_height)) + &
          & thickness_ave_ox_growth_init(j_height, id_previous_outage, imagnetite)

        ! initialize to previous value of the oxide thickness
        dz2_old(first_oxide_layer:start_exfoliation_layers - 1) = &
          & oxide_thickness_layer_height(first_oxide_layer &
          & :start_exfoliation_layers - 1, j_height)**2

        dz2_old(imagnetite) = thickness_ave_ox_growth_init( &
          & j_height, id_previous_outage, imagnetite)**2

      endif

      ! this loop is needed since the growth temperature is known apriori
      ! it may be easier to test the convergence for entire oxide thickness
      ! as the oxide thickness is small; or, use the 
      CONV_LOOP1:  do ik = 1, n_it_max 
 
        ! store the error in the growth temperature over all the layers
        ! and during the entire micro time and
       
        temp_growth_error = 0.0
        temp_growth_error_max = -10.0

        do ilay = first_oxide_layer, no_oxide_layers 

          ! initialize to previous value of the oxide thickness
          ! initialize to previous value of the increment in oxide thickness
          dz2_skin(ilay) =  (thickness_new_ox_growth_init(j_height, &
             & id_previous_outage, ilay) + &
             & dthick_ox(j_height, id_previous_outage, ilay))**2

        end do

          if ((i_time == i_time_start .or. i_time == i_time_end) .and. &
            & ik == 1 .and. &
            & (j_height == 1 .or. j_height == 2))  then
            
            write(out_lun, 19) id_next_outage, id_previous_outage,  &
              & j_height, time_now(i_time, 1), &
              & (thickness_new_ox_growth_init(j_height, &
              &    id_previous_outage, ilay), &
              &  oxide_thickness_layer_height(ilay, j_height), &
              &  dthick_ox(j_height, id_previous_outage, ilay), &
              &   ilay = first_oxide_layer, no_oxide_layers)

 19      format('dox1_init_location ', 3(i2, 1x), 35(1pe13.6, 1x))

          endif

          if (i_time == i_time_start .and. ik == 1 .and. j_height == 2)  then
            
            write(out_lun, 33) id_next_outage, id_previous_outage,  &
              & j_height, time_now(i_time, 1), &
              & (thickness_new_ox_growth_init(j_height, &
              &    id_previous_outage, ilay), &
              &  oxide_thickness_layer_height(ilay, j_height), &
              &  dthick_ox(j_height, id_previous_outage, ilay), &
              &   ilay = first_oxide_layer, no_oxide_layers)

 33      format('dox2_init_location ', 3(i2, 1x), 35(1pe13.6, 1x))

          endif

        ! there will be no_thick_interval_now(i) time subintervals now
        ! loop over the time subincrements of this interval
        TIME_MICRO: do k = 1, no_thick_interval_now(i_time)

          ! write(6, *) 'outx ', i, k, time_now(i_time, k)

          ! get the absolute time at this subinterval
          time2 = time_now(i_time, k)

          time_mean = 0.5 * (time_now(i_time, k) + time_now(i_time, k-1))

          ! get the time increment
          delta_t = time_now(i_time, k) - time_now(i_time, k-1)

          ! get the steam temperature at this height and at this time
          ! get average quantities
          Temp_steam_now = Temp_steam_start + Temp_steam_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! get the gas temperature at this height and at this time
          Temp_gas_now = Temp_gas_start + Temp_gas_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! make sure that ash_thickness was computed or inputed before
          ! ash thickness will be considered only when
          ! id_temp_op(i_time - 1) = id_temp_op(i_time) or when k_load_start = k_load_end
          ! ash thickness will affect only the outside metal temperature
          ash_thickness_now = GET_ASH_THICKNESS(k, k_load_start, k_load_end)

          ! get the heat flux at this height and at this time
          heat_flux = heat_flux_start + heat_flux_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! the heat transfer rate (heat flux * area) is constant
          heat_flux = heat_flux / (1.0 + ash_thickness_now / rad_int(1))

          ! the above heat flux will work for f2l and l2f
          if (k_load_start == k_load_end .and. hflux_htc_ash_estimate > 0)  then

            ! deal with the ash;    heat_flux_start = heat_flux_end
            id_ash = MODULO(k, no_time_intervals_ash)
            if (id_ash == 0)  id_ash = no_time_intervals_ash
            htc_outer_now = htc_tube_outer_ash(j_height, id_ash, k_load_end)  
            heat_flux = heat_flux_output_ash(j_height, id_ash, k_load_end)

          endif   

          ! set steam pressure; does not depend on height; slope is defined diferently
          press_inner_now = press_inner_p(i_time-1) + press_inner_slope_p(i_time-1) * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! Set flue gas pressure; does not depend on height
          press_outer_now = press_outer_p(i_time-1) + press_outer_slope_p(i_time-1) * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! set current flow rate in the steam
          flow_rate_steam_now = flow_rate_steam_start + flow_rate_steam_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! flow_rate_steam_now_m2 = flow_rate_steam_now / (pi * rad_inner_tube**2)
          flow_rate_steam_now_m2 = flow_rate_steam_now / (pi * rad_int(2)**2)

          ! get heat transfer coefficient on the steam-metal interface and Temp_wall
          ! htc_inner=f(diam, flow_rate, heat_flux, Temp_steam, Temp_wall, Press_steam)

          if (htc_solve_steam_or_wall)  then

            call GET_STEAM_TEMP_HTC_INNER(height_now, Temp_steam_now, &
             & press_inner_now, flow_rate_steam_now_m2, heat_flux, Temp_gas_now, &
             & time_now(i_time, k), ash_thickness_now, htc_inner_new, Temp_wall_now)

          else if (.not. htc_solve_steam_or_wall)  then

            call GET_TEMP_HTC_INNER(height_now, Temp_steam_now, press_inner_now, &
             & flow_rate_steam_now_m2, heat_flux, Temp_steam_now, time_now(i_time, k), &
             & ash_thickness_now, htc_inner_new, Temp_wall_now)
          endif

          if (oxide_thickness_old < 5.0)  then  ! microns

            ! oxide too thin for growth temperature to vary a great deal
            Temp_growth_now = Temp_steam_now

          else if (oxide_thickness_old >= 5.0)  then

            call GET_GROWTH_TEMP_SIMPLE(no_oxide_layers, Temp_steam_now, &
               & Temp_wall_now, htc_inner_new, heat_flux, ash_thickness_now, &
               & Temp_growth_now)

          endif

          if (k == 1)  then

            Temp_wall_start = Temp_wall_now
            Temp_growth_start(first_oxide_layer:no_oxide_layers) = &
               & Temp_growth_now(first_oxide_layer:no_oxide_layers) 
            htc_inner_start = htc_inner_new

            if (no_thick_interval_now(i_time) == 1)  then
              Temp_wall_end = Temp_wall_start
              Temp_growth_end(first_oxide_layer:no_oxide_layers) = &
               & Temp_growth_start(first_oxide_layer:no_oxide_layers) 
              htc_inner_end = htc_inner_start
            endif

          endif

          if (k > 1 .and. k == no_thick_interval_now(i_time))  then
            Temp_wall_end = Temp_wall_now
            Temp_growth_end(first_oxide_layer:no_oxide_layers) = &
               & Temp_growth_now(first_oxide_layer:no_oxide_layers) 
            htc_inner_end = htc_inner_new
          endif

          NOG_OX1: if (const_energy_oxide <= 0.0)  then

             dz2_skin(first_oxide_layer: no_oxide_layers) = 0.0

             temp_growth_error = 0.0
             temp_growth_error_max = 0.0
   
          else if (const_energy_oxide > 0.0)  then

          ! write(2, *) 'thickness_growth_frout ', thickness_growth_fr(first_oxide_layer:no_oxide_layers) 

          LAYER_LOOP1: do ilay = first_oxide_layer, no_oxide_layers 

            ! get the oxidation rate constant for this time interval
            ! linear interpolation based on 1/T 

            ! oxide growth based on average growth temperature Temp_grow_now
            call KP_OXIDE_AVE(Temp_growth_now(ilay), kp(ilay))

            ! thickness_growth_fr should be factored in the A - coeff
            ! thickness_growth_fr(mag) = d(mag) / dox same for spinel
            dz2_increment = 2.0 * delta_t * kp(ilay) * &
                 & thickness_growth_fr(ilay)**2
            dz2_skin(ilay) = dz2_skin(ilay) + dz2_increment

            if (ik > 1)  then
              temp_growth_error = temp_growth_error + &
                 & (Temp_growth_now(ilay) - Temp_growth_old(ilay))**2

              if (temp_growth_error_max < ABS(Temp_growth_now(ilay) - &
                & Temp_growth_old(ilay)))  then
                temp_growth_error_max = ABS(Temp_growth_now(ilay) - &
                & Temp_growth_old(ilay))
              endif
            endif

          end do LAYER_LOOP1

          endif NOG_OX1

          if (if_debug_oxide_growth_conv)  then
            write(2, 6) i_time, j_height, ik, n_it_max, k_load_start, k_load_end, &
              & (dz2_old(ilay), dz2_skin(ilay), dz2_skin(ilay)- &
              & dz2_old(ilay), ilay = first_oxide_layer, no_oxide_layers), &
              & delta_oxide_thickness_now, oxide_thickness_now, &
              & temp_growth_error / error_temp_norm, &
              & temp_growth_error_max / error_temp_max
 6          format('conv_iter_micro ', i5, 1x, 5(i2, 1x), 30(1pe13.6, 1x))
          endif

        end do TIME_MICRO

        ! get the total oxide thickness of all the layers
        if (.not. if_exfoliation_oxide_growth)  then

           oxide_thickness_now = SUM(SQRT(dz2_skin(first_oxide_layer: &
                & no_oxide_layers)))

        else if (if_exfoliation_oxide_growth)  then

           oxide_thickness_not_exfol = SUM(SQRT(dz2_skin(first_oxide_layer: &
                & start_exfoliation_layers - 1)))
           ! should also be a sum over the exfoliated layers
           ! oxide_thickness_exfol = SUM(SQRT(dz2_skin(imagnetite))) - &
           oxide_thickness_exfol = SQRT(dz2_skin(imagnetite)) - &
    & thickness_new_ox_growth_init(j_height, id_previous_outage, imagnetite) + &
    & thickness_ave_ox_growth_init(j_height, id_previous_outage, imagnetite)

           oxide_thickness_now = oxide_thickness_not_exfol + &
                & oxide_thickness_exfol

        endif

        delta_oxide_thickness_now = oxide_thickness_now - oxide_thickness_old

        if (if_debug_oxide_growth_conv)  then
          write(2, 5) i_time, j_height, ik, n_it_max, k_load_start, k_load_end, &
           & (dz2_old(ilay), dz2_skin(ilay), dz2_skin(ilay)- &
           & dz2_old(ilay), ilay = first_oxide_layer, no_oxide_layers), &
           & delta_oxide_thickness_now, oxide_thickness_now
 5        format('conv_iter ', i5, 1x, 5(i2, 1x), 20(1pe13.6, 1x))
        endif

        ! test convergence 
        ! if delta_thickness is too small then test the growth temperature
        temp_growth_error = SQRT(temp_growth_error / &
          & (no_oxide_layers - first_oxide_layer + 1) * no_thick_interval_now(i_time))

        if (ik == 1)  then

          ! do not test for the temperatures since 
          ! Temp_growth_old is not available for the first iteration
          if (delta_oxide_thickness_now < oxide_thickness_now * &
              & error_oxide_thick_ratio)  then
              convergent = 4
              exit CONV_LOOP1
          else 
              ! prepare for a new iteration
              oxide_thickness_old = oxide_thickness_now
              Temp_growth_old = Temp_growth_now
          endif

        else if (ik > 1)  then

          if (temp_growth_error_max < error_temp_max)  then
            convergent = 1
            exit CONV_LOOP1
          else if (temp_growth_error_max >= error_temp_max)  then

            if (temp_growth_error < error_temp_norm)  then
              convergent = 2
              exit CONV_LOOP1
            else if (temp_growth_error >= error_temp_norm)  then
          
              if (delta_oxide_thickness_now < oxide_thickness_now * &
                & error_oxide_thick_ratio)  then
                convergent = 3
                exit CONV_LOOP1
              else 
                ! prepare for a new iteration
                oxide_thickness_old = oxide_thickness_now
                Temp_growth_old = Temp_growth_now
              endif

            endif

          endif

        endif

      end do CONV_LOOP1

        TIME_MICRO_METAL: do k = 1, no_thick_interval_now(i_time)

          ! write(6, *) 'outx ', i, k, time_now(i, k)

          ! get the absolute time at this subinterval
          time2 = time_now(i_time, k)

          time_mean = 0.5 * (time_now(i_time, k) + time_now(i_time, k-1))

          ! get the steam temperature at this height and at this time
          ! get average quantities
          Temp_steam_now = Temp_steam_start + Temp_steam_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! get the gas temperature at this height and at this time
          Temp_gas_now = Temp_gas_start + Temp_gas_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! make sure that ash_thickness was computed or inputed before
          ! ash thickness will be considered only when
          ! id_temp_op(i_time - 1) = id_temp_op(i_time) or when k_load_start = k_load_end
          ! ash thickness will affect only the outside metal temperature
          ash_thickness_now = GET_ASH_THICKNESS(k, k_load_start, k_load_end)

          ! get the heat flux at this height and at this time
          heat_flux = heat_flux_start + heat_flux_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))
          
          ! the heat transfer rate (heat flux * area) is constant
          heat_flux = heat_flux / (1.0 + ash_thickness_now / rad_int(1))

          ! set steam pressure; does not depend on height; slope is defined diferently
          press_inner_now = press_inner_p(i_time-1) + press_inner_slope_p(i_time-1) * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! Set flue gas pressure; does not depend on height
          press_outer_now = press_outer_p(i_time-1) + press_outer_slope_p(i_time-1) * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! set current flow rate in the steam
          flow_rate_steam_now = flow_rate_steam_start + flow_rate_steam_slope * &
               & (time_mean - Time_Boiler_p(i_time-1))

          ! flow_rate_steam_now_m2 = flow_rate_steam_now / (pi * rad_inner_tube**2)
          flow_rate_steam_now_m2 = flow_rate_steam_now / (pi * rad_int(2)**2)
 
          if (ABS(flow_rate_steam_now_m2) < 1.0e-12)  then
            write(6, *)  'ERROR:  flow_rate_steam_now_m2= 0 ', pi, &
              & flow_rate_steam_now, rad_int(2)
            stop
          endif

          ! get heat transfer coefficient on the steam-metal interface and Temp_wall
          ! htc_inner=f(diam, flow_rate, heat_flux, Temp_steam, Temp_wall, Press_steam)
          if (htc_solve_steam_or_wall)  then
            call GET_STEAM_TEMP_HTC_INNER(height_now, Temp_steam_now, &
             & press_inner_now, flow_rate_steam_now_m2, heat_flux, Temp_gas_now, &
             & time_now(i_time, k), ash_thickness_now, htc_inner_new, Temp_wall_now)

          else if (.not. htc_solve_steam_or_wall)  then
            call GET_TEMP_HTC_INNER(height_now, Temp_steam_now, press_inner_now, &
             & flow_rate_steam_now_m2, heat_flux, Temp_steam_now, time_now(i_time, k), &
             & ash_thickness_now, htc_inner_new, Temp_wall_now)
          endif

          call GET_TEMP_METAL(no_oxide_layers, first_oxide_layer, Temp_gas_now, &
            & heat_flux, ash_thickness_now, Temp_metal)

          if (k == 1) then
            !write(6,*) ASSOCIATED(Temp_metal_min), no_height_oxide_output, &
! & j_height, k
            if (.not. ASSOCIATED(Temp_metal_min)) then
              ALLOCATE(Temp_metal_min(1:no_height_oxide_output, &
                      & nd1_ash_max), STAT = Status)
            endif
            if (.not. ASSOCIATED(Temp_metal_mean)) then
              ALLOCATE(Temp_metal_mean(1:no_height_oxide_output, &
              & nd1_ash_max), STAT = Status)
            endif
            if (.not. ASSOCIATED(Temp_metal_max)) then
              ALLOCATE(Temp_metal_max(1:no_height_oxide_output, &
              & nd1_ash_max), STAT = Status)
            endif
          endif

          Temp_metal_min(j_height, k) = Temp_metal(1)
          Temp_metal_mean(j_height, k) = Temp_metal(2)
          Temp_metal_max(j_height, k) = Temp_metal(3)
          
          if (if_output_height_time_all)  then
            write(out_lun, 7) j_height, id_temp_op(i_time), time_now(i_time, k), &
              & Time_Boiler_p(i_time-1), Time_Boiler_p(i_time), &
              & height_now, heat_flux, Temp_steam_now, &
              & Temp_gas_now, Temp_wall_now, Temp_metal_min(j_height, k), &
              & Temp_metal_mean(j_height, k), Temp_metal_max(j_height, k), &
              & ash_thickness_now, htc_inner_new, press_inner_now, &
              & flow_rate_steam_now_m2
            write(enrg_lun, 7) j_height, id_temp_op(i_time), time_now(i_time, k), &
              & Time_Boiler_p(i_time-1), Time_Boiler_p(i_time), &
              & height_now, heat_flux, Temp_steam_now, &
              & Temp_gas_now, Temp_wall_now, Temp_metal_min(j_height, k), &
              & Temp_metal_mean(j_height, k), Temp_metal_max(j_height, k), &
              & ash_thickness_now, htc_inner_new, press_inner_now, &
              & flow_rate_steam_now_m2
 7          format('height_output_ash ', i2, 1x, i1, 1x, 30(1pe13.6, 1x))
          endif

          if (if_time_output_event)  then

            write(out_lun, 18) j_height, id_temp_op(i_time), time_now(i_time, k), &
              & Time_Boiler_p(i_time-1), Time_Boiler_p(i_time), &
              & height_now, heat_flux, Temp_steam_now, &
              & Temp_gas_now, Temp_wall_now, Temp_metal_min(j_height, k), &
              & Temp_metal_mean(j_height, k), Temp_metal_max(j_height, k), &
              & ash_thickness_now, htc_inner_new, press_inner_now, &
              & flow_rate_steam_now_m2
            write(enrg_lun, 18) j_height, id_temp_op(i_time), time_now(i_time, k), &
              & Time_Boiler_p(i_time-1), Time_Boiler_p(i_time), &
              & height_now, heat_flux, Temp_steam_now, &
              & Temp_gas_now, Temp_wall_now, Temp_metal_min(j_height, k), &
              & Temp_metal_mean(j_height, k), Temp_metal_max(j_height, k), &
              & ash_thickness_now, htc_inner_new, press_inner_now, &
              & flow_rate_steam_now_m2
 18        format('height_event_output_ash ', i2, 1x, i1, 1x, 30(1pe13.6, 1x))

          endif

          EVENT_STORE: if (if_output_event)  then

            ! if (k_load_end = id_high)  then
              ! full load: store maximum and miniumum temperatures; one value per time period
              if (k == 1)  then
                ! ash thickness should be zero; maximum metal temperatures
                ! (1, j_o, j) for full load, ash = 0
                ! (2, j_o, j) for partial load, ash = 0
                Temp_steam_out(k_load_end, j_output, j_height) = Temp_steam_now
                Temp_gas_out(k_load_end, j_output, j_height) = Temp_gas_now
                heat_flux_out(k_load_end, j_output, j_height) = heat_flux
                htc_inner_out(k_load_end, j_output, j_height) = htc_inner_new
                Temp_metal_max_out(k_load_end, j_output, j_height) = Temp_metal(3)
                Temp_metal_mean_out(k_load_end, j_output, j_height) = Temp_metal(2)
                flow_rate_out(k_load_end, j_output, j_height) = flow_rate_steam_now_m2
              endif

              if (j_height == 1)  then
                ! check ash temperature drop
                
                if (k == 1)  then

                  call OUT_CONDUCT_RESISTANCE(no_oxide_layers, first_oxide_layer, &
                    & heat_flux, ash_thickness_now, Temp_gas_now)

                endif

                if (no_time_intervals_ash > 0 .and. &
                  & k == no_time_intervals_ash)  then

                  call OUT_CONDUCT_RESISTANCE(no_oxide_layers, first_oxide_layer, &
                    & heat_flux, ash_thickness_now, Temp_gas_now)
                endif

              endif

              if (no_time_intervals_ash == 0)  then

                ! just assigned them the previous values such that there no bug
                  Temp_steam_out(k_load_end + 2, j_output, j_height) = &
                    & Temp_steam_out(k_load_end, j_output, j_height)
                  htc_inner_out(k_load_end + 2, j_output, j_height) = &
                    & htc_inner_out(k_load_end, j_output, j_height)
                  Temp_metal_max_out(k_load_end + 2, j_output, j_height) = &
                    & Temp_metal_max_out(k_load_end, j_output, j_height)
                  Temp_metal_mean_out(k_load_end + 2, j_output, j_height) = &
                    & Temp_metal_mean_out(k_load_end, j_output, j_height)

                  ! just for checking
                  Temp_gas_out(k_load_end + 2, j_output, j_height) = &
                    & Temp_gas_out(k_load_end, j_output, j_height)
                  heat_flux_out(k_load_end + 2, j_output, j_height) = &
                    & heat_flux_out(k_load_end, j_output, j_height)
                  flow_rate_out(k_load_end + 2, j_output, j_height) = &
                    & flow_rate_out(k_load_end, j_output, j_height)

              else if (no_time_intervals_ash > 0)  then

                if (k == no_time_intervals_ash)  then
                ! (3, j_o, j) for full load, ash = max
                ! (4, j_o, j) for partial load, ash = max
                  Temp_steam_out(k_load_end + 2, j_output, j_height) = Temp_steam_now
                  htc_inner_out(k_load_end + 2, j_output, j_height) = htc_inner_new
                  Temp_metal_max_out(k_load_end + 2, j_output, j_height) = Temp_metal(3)
                  Temp_metal_mean_out(k_load_end + 2, j_output, j_height) = Temp_metal(2)

                  ! just for checking
                  Temp_gas_out(k_load_end + 2, j_output, j_height) = Temp_gas_now
                  heat_flux_out(k_load_end + 2, j_output, j_height) = heat_flux
                  flow_rate_out(k_load_end + 2, j_output, j_height) = flow_rate_steam_now_m2

                endif

              endif

          endif EVENT_STORE

        end do TIME_MICRO_METAL

      ! store the growth temperature at this height
      Temp_growth_layer_height_new(first_oxide_layer:no_oxide_layers, j_height) = &
         & Temp_growth_now(first_oxide_layer:no_oxide_layers)

      ! store the oxide growth increment (between outages)
      dthick_ox(j_height, id_previous_outage, first_oxide_layer:no_oxide_layers) = &
         & SQRT(dz2_skin(first_oxide_layer:no_oxide_layers)) - &
         & thickness_new_ox_growth_init(j_height, id_previous_outage, &
         & first_oxide_layer:no_oxide_layers)

      if (.not. if_exfoliation_oxide_growth)  then

        oxide_thickness_layer_height(first_oxide_layer:no_oxide_layers, j_height) = &
          & dthick_ox(j_height, id_previous_outage, &
          & first_oxide_layer:no_oxide_layers) + &
          & thickness_new_ox_growth_init(j_height, id_previous_outage, &
          & first_oxide_layer:no_oxide_layers)

      else if (if_exfoliation_oxide_growth)  then

        oxide_thickness_layer_height(first_oxide_layer:&
          & start_exfoliation_layers - 1, j_height) = &
          & dthick_ox(j_height, id_previous_outage, &
          & first_oxide_layer:start_exfoliation_layers - 1) + &
          & thickness_new_ox_growth_init(j_height, id_previous_outage, &
          & first_oxide_layer:start_exfoliation_layers - 1)

        oxide_thickness_layer_height(imagnetite, j_height) = &
          & dthick_ox(j_height, id_previous_outage, imagnetite) + &
          & thickness_ave_ox_growth_init(j_height, id_previous_outage, &
          & imagnetite)

      endif

      if (oxide_thickness_now > 0.0)  then
        thickness_fr_var(first_oxide_layer:no_oxide_layers, j_height) = &
          &  oxide_thickness_layer_height(first_oxide_layer:no_oxide_layers, j_height) / &
          &  oxide_thickness_now
      else
        thickness_fr_var(first_oxide_layer:no_oxide_layers, j_height) = &
        &  1.0
      endif

      if (convergent == 0)  then

        if (if_debug_oxide_growth_conv)  then

          ! convergence was not reached in n_it_max iterations
          write(out_lun, 2) i_time, j_height, time2, height_now, Temp_steam_start, &
          & Temp_steam_end, heat_flux_start, heat_flux_end, &
          & Temp_gas_start, Temp_gas_end, &
          & oxide_thickness_old, oxide_thickness_now, &
          & temp_growth_error / error_temp_norm, &
          & temp_growth_error_max / error_temp_max, &
          & delta_oxide_thickness_now / oxide_thickness_now
 2        format('convergence_not_reached ', i5, 1x, i2, 1x, 20(1pe13.6, 1x))

        endif

      else if (convergent > 0)  then

        if (oxide_thickness_now > 0)  then

         write(out_lun, 3) i_time, j_height, ik, convergent, time2, height_now, &
          & Temp_steam_start, Temp_steam_end, heat_flux_start,  &
          & heat_flux_end, Temp_gas_start, Temp_gas_end, &
          & oxide_thickness_old, oxide_thickness_now, &
          & temp_growth_error / error_temp_norm, &
          & temp_growth_error_max / error_temp_max, &
          & delta_oxide_thickness_now / oxide_thickness_now, &
          & oxide_thickness_layer_height(imagnetite, j_height), &
          & dthick_ox(j_height, id_previous_outage, imagnetite), &
          & thickness_new_ox_growth_init(j_height, id_previous_outage, imagnetite)

        else

         write(out_lun, 3) i_time, j_height, ik, convergent, time2, height_now, &
          & Temp_steam_start, Temp_steam_end, heat_flux_start,  &
          & heat_flux_end, Temp_gas_start, Temp_gas_end, &
          & oxide_thickness_old, oxide_thickness_now, &
          & temp_growth_error / error_temp_norm, &
          & temp_growth_error_max / error_temp_max

        endif

 3      format('convergence_reached ', i5, 1x, 3(i2, 1x), 35(1pe13.6, 1x))

         write(out_lun, 4) i_time, j_height, ik, convergent, time2, height_now, &
          & Temp_steam_start, Temp_steam_end, heat_flux_start,  &
          & heat_flux_end, Temp_gas_start, Temp_gas_end, &
          & oxide_thickness_old, oxide_thickness_now, &
          & Temp_wall_start, Temp_wall_end, htc_inner_start, &
          & htc_inner_end, (Temp_growth_start(ilay), &
          & Temp_growth_end(ilay), ilay = first_oxide_layer, no_oxide_layers)
 4      format('conv_done_ox_thick ', i5, 1x, 3(i2, 1x), 30(1pe13.6, 1x))

      endif

      if (j_height == 2)  then
        ! keep this only for debugging
        write(out_lun,52) i_time, j_height, ik, (thickness_layer(n1), Temp_growth_now(n1), &
           & Temp_growth_old(n1),  n1 = first_oxide_layer, no_oxide_layers)
 52     format(' temp_conv1 ', i5, 1x, 2(i2, 1x), 10(1pe13.6, 1x))

            write(out_lun, 31) i_time, time2, kp(first_oxide_layer), &
               & 2.0*sqrt(kp(first_oxide_layer)), &
               & dz2_skin(first_oxide_layer:no_oxide_layers), &
               & sqrt(dz2_skin(first_oxide_layer:no_oxide_layers)), &
               & oxide_thickness_now 
 31      format('dz_check ', i4, 1x, 30(1pe13.6, 1x))
       endif

       if (oxide_thickness_now > 1.0e+6 .or. &
             & ANY(ABS(Temp_growth_now(:)) > 2.0e+3))  then
          write(tty_lun, *) 'ERROR: oxide_thickness update error; iter= ', ik, i
          write(tty_lun,52) i_time, j_height, ik, oxide_thickness_now, (Temp_growth_now(n1), &
           & Temp_growth_old(n1),  n1 = first_oxide_layer, no_oxide_layers)
            stop
        endif

        ! dthick_ox_height_time(id_previous_outage, j_height, i_time) = &
        !   & dthick_ox(j_height, id_previous_outage, imagnetite)

        call TEMP_OUT_STEAM_OXIDE_METAL(i_time, id_next_outage, j_height, &
             & Temp_metal, heat_flux)

      end do HEIGHT_LOOP1

      do k = 1, no_thick_interval_now(i_time)

        ! write the temperatures at oxide-metal interface for all heights
        write(2, 8) i_time, time_now(i_time, k), (Temp_metal_min(j, k), &
          & j = 1, no_height_oxide_output)
        write(enrg_lun, 8) i_time, time_now(i_time, k), (Temp_metal_min(j, k), &
          & j = 1, no_height_oxide_output)
 8        format('metal_temp_min ', i5, 1x, 30(1pe13.6, 1x))

        ! write the mean metal temperatures  for all heights
        write(2, 10) i_time, time_now(i_time, k), (Temp_metal_mean(j, k), &
          & j = 1, no_height_oxide_output)
        write(enrg_lun, 10) i_time, time_now(i_time, k), (Temp_metal_mean(j, k), &
          & j = 1, no_height_oxide_output)
 10        format('metal_temp_ave ', i5, 1x, 30(1pe13.6, 1x))

        ! write the temperatures at metal-ash interface for all heights
        write(2, 12) i_time, time_now(i_time, k), (Temp_metal_max(j, k), &
          & j = 1, no_height_oxide_output)
        write(enrg_lun, 12) i_time, time_now(i_time, k), (Temp_metal_max(j, k), &
          & j = 1, no_height_oxide_output)
 12        format('metal_temp_max ', i5, 1x, 30(1pe13.6, 1x))

      end do

      ! write out the height dependence; wait till the low load output is stored
      if (if_output_event .and. k_load_end == id_low)  then

        do j = 1, no_height_oxide_output
          write(out_lun, 16) j_output, time_j_output(j_output), &
           & length_tube_oxide_output(j), &
           & Temp_gas_out(1, j_output, j), heat_flux_out(1, j_output, j), &
           & Temp_steam_out(1, j_output, j), Temp_steam_out(3, j_output, j), &
           & Temp_metal_mean_out(1, j_output, j), &
           & Temp_metal_mean_out(3, j_output, j), &
           & Temp_metal_max_out(1, j_output, j), &
           & Temp_metal_max_out(3, j_output, j), &
           & Temp_gas_out(2, j_output, j), heat_flux_out(2, j_output, j), &
           & Temp_steam_out(2, j_output, j), Temp_steam_out(4, j_output, j), &
           & Temp_metal_mean_out(2, j_output, j), &
           & Temp_metal_mean_out(4, j_output, j), &
           & Temp_metal_max_out(2, j_output, j), &
           & Temp_metal_max_out(4, j_output, j), &
           & htc_inner_out(1, j_output, j), htc_inner_out(3, j_output, j), &
           & htc_inner_out(2, j_output, j), htc_inner_out(4, j_output, j), &
           & Temp_gas_out(3, j_output, j), Temp_gas_out(4, j_output, j), &
           & heat_flux_out(3, j_output, j), heat_flux_out(4, j_output, j), &
           & flow_rate_out(1, j_output, j), flow_rate_out(2, j_output, j), &
           & flow_rate_out(3, j_output, j), flow_rate_out(4, j_output, j)
          write(enrg_lun, 16) j_output, time_j_output(j_output), &
           & length_tube_oxide_output(j), &
           & Temp_gas_out(1, j_output, j), heat_flux_out(1, j_output, j), &
           & Temp_steam_out(1, j_output, j), Temp_steam_out(3, j_output, j), &
           & Temp_metal_mean_out(1, j_output, j), &
           & Temp_metal_mean_out(3, j_output, j), &
           & Temp_metal_max_out(1, j_output, j), &
           & Temp_metal_max_out(3, j_output, j), &
           & Temp_gas_out(2, j_output, j), heat_flux_out(2, j_output, j), &
           & Temp_steam_out(2, j_output, j), Temp_steam_out(4, j_output, j), &
           & Temp_metal_mean_out(2, j_output, j), &
           & Temp_metal_mean_out(4, j_output, j), &
           & Temp_metal_max_out(2, j_output, j), &
           & Temp_metal_max_out(4, j_output, j), &
           & htc_inner_out(1, j_output, j), htc_inner_out(3, j_output, j), &
           & htc_inner_out(2, j_output, j), htc_inner_out(4, j_output, j), &
           & Temp_gas_out(3, j_output, j), Temp_gas_out(4, j_output, j), &
           & heat_flux_out(3, j_output, j), heat_flux_out(4, j_output, j), &
           & flow_rate_out(1, j_output, j), flow_rate_out(2, j_output, j), &
           & flow_rate_out(3, j_output, j), flow_rate_out(4, j_output, j)
        end do

 16     format('height_output_ox_gr_event ', i1, 1x, 50(1pe13.6, 1x))

      endif

      thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
        & 1: no_height_oxide_output, i_time) = &
        & oxide_thickness_layer_height(first_oxide_layer: no_oxide_layers, &
        & 1: no_height_oxide_output)

      ! output for height 1
      !  & thickness_fr_var(i_time, first_oxide_layer:no_oxide_layers), &
      if (i_time == 2)  then

        ! store height id at which the exfoliation detail will be given
        iex_out(1) = 1 ! - inlet on first loop or inlet to "inlet loop"

        ! midloop if one loop or inlet to the "outlet loop"
        iex_out(2) = (1 + no_height_oxide_output) / 2
        
        ! outlet
        iex_out(3) = no_height_oxide_output

        ! write header
        write(out_lun, 210) (' d'//ch_index(1+3*j)//'\dox', &
      &   ' d'//ch_index(1+3*j)//'\dox,m', &
      &   ' d'//ch_index(1+3*j)//'\dox,s', &
      &   j = 1, MIN(2, no_height_oxide_output/3 - 1)), &
      &   (' d'//'1'//ch_index(1+3*j-10)//'\dox', &
      &   ' d'//'1'//ch_index(1+3*j-10)//'\dox,m', &
      &   ' d'//'1'//ch_index(1+3*j-10)//'\dox,s', &
      &   j = 1+MIN(2, no_height_oxide_output/3 - 1), &
      &         MIN(6, no_height_oxide_output/3 - 1)), &
      &   (' d'//'2'//ch_index(1+3*j-20)//'\dox', &
      &   ' d'//'2'//ch_index(1+3*j-20)//'\dox,m', &
      &   ' d'//'2'//ch_index(1+3*j-20)//'\dox,s', &
      &   j = 1+MIN(6, no_height_oxide_output/3-1), no_height_oxide_output/3-1)
 210    format('oxide_val1 it id_op Time[h] Tst ', &
          & 'd_i\dox d_i\dox,mag \fD\nd_i\dox,mag d_i\dox,mag\n-kin ', &
          & 'd_i\dox,mag\n-exfol d_i\dox,sp ', &
          & 'd_io\dox d_io\dox,mag \fD\nd_io\dox,mag d_io\dox,mag\n-kin ', &
          & 'd_io\dox,mag\n-exfol d_io\dox,sp ', &
          & 'd_o\dox d_o\dox,mag \fD\nd_o\dox,mag d_o\dox,mag\n-kin ', &
          & 'd_o\dox,mag\n-exfol d_o\dox,sp ', 100(a))

      endif

      write(out_lun,209) i_time, id_temp_op(i_time), Time_Boiler_p(i_time),  &
        & TEMP_STEAM_p(i_time), &
        & (SUM(oxide_thickness_layer_height(first_oxide_layer: &
           & no_oxide_layers, iex_out(j))), &
        & oxide_thickness_layer_height(imagnetite, iex_out(j)), &
        & dthick_ox(iex_out(j), id_previous_outage, imagnetite), &
        & thickness_new_ox_growth_init(iex_out(j), id_previous_outage, &
           & imagnetite), &
        & thickness_ave_ox_growth_init(iex_out(j), id_previous_outage, &
           & imagnetite), &
        & oxide_thickness_layer_height(first_oxide_layer, iex_out(j)), &
           & j = 1, 3), &
        & (SUM(oxide_thickness_layer_height(first_oxide_layer: &
          & no_oxide_layers, 1+3*j)), &
        & oxide_thickness_layer_height(imagnetite, 1+3*j), &
        & oxide_thickness_layer_height(first_oxide_layer, 1+3*j), &
          & j = 1, no_height_oxide_output/3)
        !&   j = 0, MIN(2, no_height_oxide_output/3)), &

 ! 208    format('oxide_val ', 35(1pe13.6, 1x))
 209    format('oxide_val1 ', i4, 1x, i2, 1x, 100(1pe13.6, 1x))

    end do TIME_LOOP_HOX1

    ! allocate steam temperature
    ! to save space it is allocated only for intervals between outages
    ! DEALLOCATE(TEMP_STEAM_HEIGHT_TIME, STAT = Status)

        if (ANY(thickness_layer > 1.0) .or. &
           & ANY(ABS(Temp_growth_now(:)) > 2.0e+3) .or. &
           & ANY(ABS(Temp_growth_old(:)) > 2.0e+3))  then
          write(tty_lun, *) 'ERROR: Temp_growth update error; iter= ', ik, j
          write(tty_lun,52) ik, j, (thickness_layer(i), Temp_growth_now(i), &
           & Temp_growth_old(i),  i = first_oxide_layer, no_oxide_layers)
          STOP
        endif
      
    return

  END SUBROUTINE OXIDE_GROWTH_LOCATION

  SUBROUTINE GET_GROWTH_TEMP_SIMPLE(N, steam_temperature, Temp_wall_now, &
         & htc_inner, heat_flux, ash_thickness_now, growth_temp_this_layer)
    !=======================================================================
    ! Purpose(s):
    !
    !   obtain the growth temperature per each layer 
    !   a simplified heat transfer is used 
    ! 
    !=======================================================================

  use solution_data_module, only: rad_int
  use oxide_data_module,    only: growth_ox_temperature, first_oxide_layer, &
          & cond_value
  use htc_utility_module,   only: GET_TEMP_INTERFACE

  implicit none

  ! Argument List
  integer,               intent(IN)  :: N
  real,                  intent(IN)  :: steam_temperature, htc_inner, &
         & heat_flux, Temp_wall_now, ash_thickness_now
  real, dimension(2:N),  intent(OUT) :: growth_temp_this_layer

  real, dimension(1:N+1) :: Temp_int
  real, dimension(2:N) :: dummy

  integer :: j, i  

  ! get the temperatures at oxide interfaces
  call GET_TEMP_INTERFACE(N, first_oxide_layer, steam_temperature, &
     & Temp_wall_now, htc_inner, heat_flux, ash_thickness_now, Temp_int)
     
     do j = first_oxide_layer, N

      SELECT CASE (TRIM(growth_ox_temperature(j)))

        CASE ('inner_surface')

          growth_temp_this_layer(j) = Temp_int(j) ! Temp(j, 1)

        CASE ('outer_surface')

          growth_temp_this_layer(j) = Temp_int(j + 1) ! Temp(j, npr_temp(j))

        CASE ('average_layer')

          growth_temp_this_layer(j) = 0.5 * (Temp_int(j) + Temp_int(j + 1))

          ! SUM(Temp(j, 1:npr_temp(j))) / npr_temp(j)

        CASE ('average_all_layers')

          do i = first_oxide_layer, N
            ! dummy(i) = SUM(Temp(i, 1:npr_temp(i))) / npr_temp(i)
            dummy(i) = 0.5 * (Temp_int(i) + Temp_int(i + 1))
          end do
          growth_temp_this_layer(j) = SUM(dummy(first_oxide_layer:N)) / &
             & (N-first_oxide_layer + 1)

        CASE ('steam')

          growth_temp_this_layer(j) = steam_temperature ! TEMP_STEAM_p(cycle_no)

        CASE ('oxide_surface')

          ! this would be the wall temperature
          growth_temp_this_layer(j) = Temp_wall_now  ! Temp(N, npr_temp(N))

        CASE DEFAULT

          write(6, *) 'ERROR: growth_ox_temperature not ok'
          STOP

      END SELECT

    end do

    RETURN

  END SUBROUTINE GET_GROWTH_TEMP_SIMPLE

  SUBROUTINE KP_OXIDE_AVE(Temp_ave, KP)
    !=======================================================================
    ! Purpose(s):
    !
    !  Compute the oxide thickness based on kinetics and average temperature
    ! 
    !=======================================================================
    use parameter_module,   only: mfuel
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness
    use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, log_const_energy_oxide, &
          & activ_energy_oxide_r, if_function_or_table_oxide
    use solver_data_module, only: no_thick_interval

    ! Arguments
    real, intent(IN)     :: Temp_ave
    real, intent(OUT)    :: kp

    ! Local Variables
    integer :: i, j, k, Status
    real    :: delta_t, kp_log, one_over_temp
    logical, save :: if_function_kp = .true.
 
  if (if_function_kp)  then

    kp_log = log_const_energy_oxide - activ_energy_oxide_r / &
           & (temp_kelvin_unit + temp_ave)

  else if (.not. if_function_kp)  then

    ! constant temperature
    LOOP1: do j = 1, nrate - 1

      if (Temp_ave >= oxide_rate_temp(j) .and. &
          & Temp_ave <= oxide_rate_temp(j+1))  then

        one_over_temp = 1.0 / (temp_kelvin_unit + temp_ave)
        kp_log = (oxide_rate_log(j) * &
             & (oxide_rate_1temp(j+1) - one_over_temp) + &
             & oxide_rate_log(j+1) * &
             & (one_over_temp - oxide_rate_1temp(j))) / &
             & (oxide_rate_1temp(j+1) - oxide_rate_1temp(j))

      !    write(out_lun, 11) i, k, j, oxide_rate_temp(j), temp1w, temp2w, &
      !   & oxide_rate_temp(j+1), Time_Boiler_p(i), temp_ave, kp_log
      !   11 format('check1 ', 3(i4, 1x), 10(1pe13.6, 1x))

         exit LOOP1

       endif

     end do LOOP1

   endif

   kp = EXP(kp_log)

   return

  END SUBROUTINE KP_OXIDE_AVE

END MODULE OXIDE_GROWTH_MODULE
