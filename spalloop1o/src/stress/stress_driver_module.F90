MODULE STRESS_DRIVER_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   stress analyis during boiler operation
    !
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: STRESS_LOCATION

  CONTAINS

  SUBROUTINE STRESS_LOCATION()

    !=======================================================================
    ! Purpose(s):
    !
    !  obtain strains, stress, at different temperatures 
    !  for the tube, i.e., hollow cylinder case
    !  add 1-2 and 4-1 for creep; used to be CYCLE_STRAIN_CYL_ALL in spallation
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
    use oxide_data_module,  only: temp_kelvin_unit, &
          & if_temp_steam_or_oxide, thickness_fr, no_oxide_ave, &
          & growth_ox_temperature, thickness_growth_fr, first_oxide_layer
    use oxide_data_module,  only: e_th, e_th_mat, id_mat, no_layer, &
          & material_name
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use solver_data_module, only: temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, no_out_f2l, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use mesh_module, only: GET_RADIUS, GET_TUBE_DIMENSIONS
    use solution_data_module, only: Temp, npr_temp, no_oxide_layers, &
          & oxide_thickness_hf_out, time_hf_out, no_full2low_output, &
          & no_point_full2low, no_pulse_full2low_output, &
          & id_full_load_hf, n_hf_out, rad_int
    use solution_data_module, only:  strain_elastic_event_out, &
      & energy_elastic_event_out, &
      & temperature_event_out, no_pulse_full2low_output, &
      & max_energy_elastic_event, max_energy_elastic_event_out, &
      & min_strain_elast_event, min_strain_elast_event_out, &
      & max_strain_elast_event, max_strain_elast_event_out, &
      & strain_elastic_event_out_g, min_strain_elast_event_gen, &
      & max_strain_elast_event_gen
    use solution_data_module, only: rad_int_event, time_event, &
      & oxide_thickness_event, no_f2l_event_out, id_outage, &
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
    use waterwall_data_module, only: thickness_oxide_height_time, &
        & height_id_stress_out, id_height_now
    use transit_out_module,    only: F2L_TRANSITION_OUT
    use solver_data_module, only: if_average_oxide
    use oxide_stress_data_module,   only: total_oxide_thickness_ref, &
      & d_oxide_thickness_4_strain_hoop, ispinel, imagnetite
    use creep_data_module, only: if_creep_ox_scale, if_creep_metal, &
      & operation_id_current
    use oxide_utility_module, only: COOLING_STRAIN_PLANE
    use blockage_data_module, only: max_energy_elastic_event_block, &
      & min_strain_elast_event_block, max_strain_elast_event_block, &
      & id_height_exfoliate2deposit, rad_int_block, oxide_thickness_block
    use creep_module,              only: INIT_CREEP
    use waterwall_data_module,     only: if_new_height, if_outage_ramp_exfol, &
      & jheight_creep, ioutage_creep, ioutage_previous_creep
    use oxide_growth_module,       only: OXIDE_GROWTH_LOCATION
    use blockage_data_module,  only: dthick_ox, &
      & growth_mode_after_exfoliation, thickness_new_ox_growth, &
      & thickness_new_ox_growth_init, if_exfoliation_oxide_growth, &
      & thickness_ave_ox_growth_init, start_exfoliation_layers, &
      & end_exfoliation_layers
    use exfoliation_module,  only:  THICKNESS_INIT_EXFOL, EXFOLIATED_VOLUME
    use creep_exfoliation_module, only: CREEP_RESET_EXFOLIATION, &
      & STRESS_STRAIN_CREEP_OUT_EXFOL

    implicit none

    ! Argument List

    ! Local Variables
    integer :: i, j, k, Status, n_it_max, n1, n2, n21, ik, ilay, &
             & ihf, it_dummy = 1, it_dummy10 = 10, k_load_start, &
             & k_load_end, step_cycle, nstart
    real    :: delta_t, temp1w, temp2w, &
             & one_over_temp, temp_ave, time2
    real    :: time_dummy = 0.0, oxide_thickness_old
    real, dimension(2:no_oxide_layers)    :: kp, kp_log, temp2, temp1,  &
             & temp_slope, error_temp, dtemp, Temp_growth_now, dz2_old, &
             & Temp_growth_old, dz2, Temp_growth_end, Temp_growth_start

    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(3) :: Temp_metal
    integer  :: convergent, j_output, id_ash, id_next_outage, &
          & id_previous_outage, k_height

    real :: flow_rate_steam_slope, flow_rate_steam_end, height_now, &
          & flow_rate_steam_start, Temp_steam_start, Temp_steam_end, &
          & Temp_steam_slope, Temp_gas_start, Temp_gas_end, Temp_gas_slope, &
          & heat_flux_start, heat_flux_end, heat_flux_slope, Temp_steam_now, &
          & Temp_gas_now, heat_flux, press_inner_now, press_outer_now, &
          & flow_rate_steam_now, htc_inner_new, Temp_wall_now, &
          & ash_thickness_now, oxide_thickness_now, &
          & delta_oxide_thickness_now, delta_oxide_thickness_old, &
          & temp_growth_error_max, htc_outer_now
    real :: temp_growth_error, htc_inner_start, &
          & Temp_wall_start, Temp_wall_end, htc_inner_end, time_mean, &
          & flow_rate_steam_now_m2, rad_inner_tube, press_outer_start, &
          & press_inner_start, press_outer_end, press_inner_end
    real, dimension(4, no_full2low_output, no_height_oxide_output) :: Temp_steam_out, &
          & htc_inner_out, Temp_metal_mean_out, Temp_metal_max_out, &
          & heat_flux_out, Temp_gas_out, flow_rate_out
    real, dimension(no_height_oxide_output) :: oxide_thickness_out

    logical :: if_output_event, if_time_output_event
    logical, save:: if_outage_ramp

    if (if_exfoliation_oxide_growth .and. &
      & end_exfoliation_layers /= no_oxide_layers)  then

      write(tty_lun, *) 'ERROR: end_exfoliation_layers= ', &
         & end_exfoliation_layers, ' must be no_oxide_layers= ', no_oxide_layers
      STOP

    endif

    if (if_exfoliation_oxide_growth .and. &
      & start_exfoliation_layers > no_oxide_layers)  then

      write(tty_lun, *) 'ERROR: start_exfoliation_layers > no_oxide_layers'
      STOP

    endif

    ! usually no_thick_interval = 10 -> 6 min output interval
    ! use 2 min output 
    no_out_f2l = 3 * no_thick_interval

    nstart = 2  ! time iteration at which time loop starts
    ! creep should be ran initially; before oxide starts to grow
    ! the stress is made to ran with a certain oxide layer

    ! initialize
    id_outage = 0

  id_next_outage = 1
  id_previous_outage = 0  ! this would be the baseline id

  thickness_new_ox_growth_init(:, id_previous_outage, :) = 0.0
  thickness_ave_ox_growth_init(:, id_previous_outage, :) = 0.0
  thickness_new_ox_growth_init(:, id_next_outage, :) = 0.0
  thickness_ave_ox_growth_init(:, id_next_outage, :) = 0.0

  if_outage_ramp_exfol = .false.
    
  ! reset the oxide thickness between outages
  dthick_ox = 0.0

  OPER2: do n1 = 1, no_outages
    ! this is because oxide_thickness(i+1) is required at SET_OX_STRAINS
    ! oxide_thickness(cycle_no) - oxide_thickness(cycle_no-1)

    id_previous_outage = n1 - 1
    id_next_outage = n1

    ! data for the creep
    ioutage_creep = id_next_outage
    ioutage_previous_creep = ioutage_creep - 1

    ! id_time_outage(k) = 2, so it is consistent with nstart = 2

    ! step 1 - get oxide thickness between outages for all the heights
    call OXIDE_GROWTH_LOCATION(1 + id_time_outage(n1-1), id_time_outage(n1), &
           & if_outage_ramp_exfol, id_previous_outage, id_next_outage)

    ! at each height determine the oxide thickness
    HEIGHT_LOOP1: do k = 1, no_height_oxide_output

      k_height = k

      id_height_now = k

      if_new_height = .true.

      ! data for the creep
      jheight_creep = k_height

      height_now = length_tube_oxide_output(k_height)

      ! set tube_outer_radius and tube_thickness at this location
      call GET_TUBE_DIMENSIONS(k_height)

      write(2, *) 'start_height_loop ', k_height, height_now

      ! initialize creep data;
      ! even for cases when there is not creep some data must be set
      ! creep continuity in time must be considered
      call INIT_CREEP(no_oxide_layers, k_height, n1)

      if_outage_ramp_exfol = .false.

      ! OPER1: do i= nstart, TEMP_BOILER_data_no-1 ! no of data intervals
      OPER1: do i= 1 + id_time_outage(n1-1), id_time_outage(n1)

        ! operation_id_current = stores id_temp_op(i) current
        operation_id_current = id_temp_op(i)

        if (operation_id_current <= 0)  CYCLE OPER1

        ! get the time at the start of the i-th interval
        time2 = Time_Boiler_p(i-1)

        if (i == id_time_outage(n1))  then
          ! baseline used should be id_previous_outage = id_next_outage - 1
          if_outage_ramp_exfol = .true.

          if (k_height == 1)  then
            write(out_lun, 20) id_next_outage, id_previous_outage, time2, &
             & Time_Boiler_p(i)
 20          format('outage_2_follow ', 2(i2, 1x), 30(1pe13.6, 1x))
          endif

        endif

        ! step 2 - 
        ! this does not deal with the outage; but it can be made to in order
        ! to simplify things
        call SETUP_LOAD_OUT(k_height, i, k_load_start, k_load_end, if_output_event, &
            & if_time_output_event, j_output, time2)

        ! set steam flow rate; does not depend on height  
        ! ACTION 1: flow_rate_steam quantities are defined
        flow_rate_steam_start = LOAD_FUNCTION(k_load_start, &
          & flow_rate_steam_high, flow_rate_steam_low)
        flow_rate_steam_end = LOAD_FUNCTION(k_load_end, &
          & flow_rate_steam_high, flow_rate_steam_low)
        flow_rate_steam_slope = (flow_rate_steam_end - &
          & flow_rate_steam_start) / &
          & (Time_Boiler_p(i) - Time_Boiler_p(i-1))

        ! set steam temperature at this height; begin, end, slope
        Temp_steam_start = Temp_steam_height(k_height, k_load_start)
        Temp_steam_end = Temp_steam_height(k_height, k_load_end)
        Temp_steam_slope = (Temp_steam_end - Temp_steam_start) / &
          & (Time_Boiler_p(i) - Time_Boiler_p(i-1))

        ! set gas temperature at this height; begin, end, slope
        Temp_gas_start = Temp_gas_height(k_height, k_load_start)
        Temp_gas_end = Temp_gas_height(k_height, k_load_end)
        Temp_gas_slope = (Temp_gas_end - Temp_gas_start) / &
          & (Time_Boiler_p(i) - Time_Boiler_p(i-1))

        ! set heat flux begin, end, slope
        heat_flux_start = heat_flux_output(k_height, k_load_start)
        heat_flux_end = heat_flux_output(k_height, k_load_end)
        heat_flux_slope = (heat_flux_end - heat_flux_start) / &
          & (Time_Boiler_p(i) - Time_Boiler_p(i-1))

        press_outer_start = press_outer_p(i-1)
        press_inner_start = press_inner_p(i-1)
        press_outer_end = press_outer_p(i)
        press_inner_end= press_inner_p(i)

        ! get radii at the interfaces: 
        ! steam-metal, metal-oxide, oxide-oxide, oxide-steam
        ! assumption: 
        ! oxide thickness - for temperature computations does not change
        !     within the micro-time loop
        thickness_layer(1:first_oxide_layer-1) = 0.0 ! metal layers

        thickness_layer(first_oxide_layer:no_oxide_layers) = oxide_thick2si * &
            & thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
            & k_height, i)

        ! oxide_thickness is in microns while thickness_layer in SI 
        ! oxide_thickness(i) = SUM(thickness_oxide_height_time(  &
        !    & first_oxide_layer:no_oxide_layers, k_height, i))
        oxide_thickness(i) = SUM(thickness_layer( &
            & first_oxide_layer:no_oxide_layers)) / oxide_thick2si

        thickness_fr(first_oxide_layer:no_oxide_layers) = &
            & thickness_layer(first_oxide_layer:no_oxide_layers) / &
            & SUM(thickness_layer(first_oxide_layer:no_oxide_layers))

        ! oxide thickness [micron]
        ! without thermal expasion; just given by the kinetics
        total_oxide_thickness_ref = oxide_thickness(i)

        write(out_lun, 9)  k_height, i, operation_id_current, id_next_outage, &
           & Time_Boiler_p(i-1), Time_Boiler_p(i), oxide_thickness(i), &
           & thickness_fr(imagnetite), &
           & thickness_oxide_height_time(imagnetite, k_height, i), &
           & thickness_ave_ox_growth_init(k_height, id_previous_outage, imagnetite), &
           & Temp_steam_start, Temp_steam_end, Temp_gas_start, Temp_gas_end, &
           & press_outer_start, press_outer_end, &
           & press_inner_start, press_inner_end
 9      format('start_time_loop ', i2, 1x, i5, 1x, 2(i2, 1x), 20(1pe13.6, 1x))

        ! if_outage_ramp_b = if_outage_ramp

        if (operation_id_current == 1 .or. operation_id_current == 5)  then

          ! end of idle and beginning of a new cycle ramp up

          e_th_scale = e_th_scale
          ! sigma_th_scale = sigma_th_scale

          call ISOTHERM_STATE_STRAIN(k_height, i, if_outage_ramp, thickness_layer, &
              & Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Time_Boiler_p(i), Time_Boiler_p(i) - Time_Boiler_p(i-1))

        else if (operation_id_current == 2 .or. operation_id_current == 6)  then
    
          ! end of ramp-up and beginning of a new pulse segment

          e_th_scale = 0.0
          ! sigma_th_scale = 0.0

          call TRANSITION_STATE_STRAIN(k_height, i, if_outage_ramp, thickness_layer, &
              & Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Temp_gas_end, Temp_steam_end, &
              & press_outer_end, press_inner_end, &
              & Time_Boiler_p(i-1), Time_Boiler_p(i))

        else if (operation_id_current == 3 .or. operation_id_current == 7)  then

          ! end of pulse and beginning of a new ramp-down segment

          e_th_scale = 0.0
          ! sigma_th_scale = 0.0

          call ISOTHERM_STATE_STRAIN(k_height, i, if_outage_ramp, thickness_layer, &
              & Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Time_Boiler_p(i), Time_Boiler_p(i) - Time_Boiler_p(i-1))

           ! if (k == -1)  then !  height_id_stress_out)  then
            ! print output only for one height
            ! if_outage_ramp = .true. if outage now
            call F2L_TRANSITION_OUT(k_height, i, id_next_outage, if_outage_ramp, &
              & if_outage_ramp_exfol, thickness_layer, &
              & oxide_thickness(i) - oxide_thickness(i-1))
           ! endif

           if (if_outage_ramp_exfol)  then

             ! store info from outage results for possible tube blockage 
             rad_int_block(k_height, id_next_outage, :) = rad_int_event(id_next_outage, :)
             oxide_thickness_block(k_height, id_next_outage, :) = &
                  & oxide_thickness_event(id_next_outage, :)
             max_energy_elastic_event_block(k_height, id_next_outage, :) = &
                & max_energy_elastic_event_out(id_next_outage, :)

            write(out_lun, 6) id_next_outage, k_height, id_next_outage, &
               & oxide_thickness(i), &
               & Time_Boiler_p(i)
 6          format('must_max_energy_elastic_event_out ', 3(i2, 1x), 20(1pe13.6, 1x))

             ! set this for all layers
             ! for the exfoliated layers this will be overridden in the body of this routine
             if (.not. if_exfoliation_oxide_growth)  then

               thickness_ave_ox_growth_init(k_height, id_next_outage, &
                    & first_oxide_layer:no_oxide_layers) = &
                 & thickness_ave_ox_growth_init(k_height, id_previous_outage, &
                    & first_oxide_layer:no_oxide_layers) + &
                 & dthick_ox(k_height, id_previous_outage, &
                    & first_oxide_layer:no_oxide_layers)

               thickness_new_ox_growth_init(k_height, id_next_outage, &
                 & first_oxide_layer:no_oxide_layers) = &
                 & dthick_ox(k_height, id_previous_outage, &
                 & first_oxide_layer:no_oxide_layers) + &
                 & thickness_new_ox_growth_init(k_height, id_previous_outage, &
                 & first_oxide_layer:no_oxide_layers)

             else if (if_exfoliation_oxide_growth)  then

               thickness_ave_ox_growth_init(k_height, id_next_outage, &
                    & first_oxide_layer:start_exfoliation_layers - 1) = &
                 & thickness_ave_ox_growth_init(k_height, id_previous_outage, &
                    & first_oxide_layer:start_exfoliation_layers - 1) + &
                 & dthick_ox(k_height, id_previous_outage, &
                    & first_oxide_layer:start_exfoliation_layers - 1)
               ! the rest will be taking care later on

               thickness_new_ox_growth_init(k_height, id_next_outage, &
                 & first_oxide_layer:start_exfoliation_layers - 1) = &
                 & dthick_ox(k_height, id_previous_outage, &
                 & first_oxide_layer:start_exfoliation_layers - 1) + &
                 & thickness_new_ox_growth_init(k_height, id_previous_outage, &
                 & first_oxide_layer:start_exfoliation_layers - 1)
               ! the rest will be taking care later on

             endif

             ! do ioutage = 1, no_f2l_event_out  ! no_total_outage
             call EXFOLIATED_VOLUME(no_oxide_layers, k_height, &
                & id_previous_outage, id_next_outage)
             ! end do

             ! set thickness_new_ox_growth_init after new exfoliation event
             call THICKNESS_INIT_EXFOL(k_height, id_previous_outage, id_next_outage)

             ! reset the creep strains in thiner oxide after exfoliation
             ! used also to store von Mises stress
             call CREEP_RESET_EXFOLIATION(no_oxide_layers, time2, &
                & k_height, id_next_outage)

             !  if_outage_ramp_exfol = .false.

             if (k_height == 1)  then

               write(out_lun, 19) id_previous_outage, id_next_outage, &
                 & time2, (dthick_ox(k_height, id_previous_outage, n2), &
                 & thickness_ave_ox_growth_init(k_height, id_previous_outage, n2), &
                 & thickness_ave_ox_growth_init(k_height, id_next_outage, n2), &
                 & thickness_new_ox_growth_init(k_height, id_previous_outage, n2), &
                 & thickness_new_ox_growth_init(k_height, id_next_outage, n2), &
                 & max_energy_elastic_event_block(k_height, id_next_outage, n2), &
                 & n2 = 2, no_oxide_layers)
 19            format('id1_outage ', 2(i2, 1x), 30(1pe10.3, 1x))

             else if (k_height == 2)  then  ! height 2

               write(out_lun, 31) id_previous_outage, id_next_outage, &
                 & time2, (dthick_ox(k_height, id_previous_outage, n2), &
                 & thickness_ave_ox_growth_init(k_height, id_previous_outage, n2), &
                 & thickness_ave_ox_growth_init(k_height, id_next_outage, n2), &
                 & thickness_new_ox_growth_init(k_height, id_previous_outage, n2), &
                 & thickness_new_ox_growth_init(k_height, id_next_outage, n2), &
                 & max_energy_elastic_event_block(k_height, id_next_outage, n2), &
                 & n2 = 2, no_oxide_layers)
 31            format('id2_outage ', 2(i2, 1x), 30(1pe10.3, 1x))

             endif

           endif

        else if (operation_id_current == 4 .or. operation_id_current == 8)  then

          ! end of ramp-down and beginning of a new idle segment

          ! realocate the 
          ! fix no_oxide_ave later on
          no_oxide_ave = no_layer + 1
          j = no_oxide_ave
          e_th_scale(1) = COOLING_STRAIN_PLANE(id_mat(j), &
            & TEMP_STEAM_p(i-1), TEMP_STEAM_p(i))

          do ik= 1, nFe2O3_pct_check_strain
            j = no_oxide_ave + ik
            e_th_scale(ik +1 ) = COOLING_STRAIN_PLANE(id_mat(j), &
            & TEMP_STEAM_p(i-1), TEMP_STEAM_p(i))
          end do

          call TRANSITION_STATE_STRAIN(k_height, i, if_outage_ramp, thickness_layer, &
              & Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Temp_gas_end, Temp_steam_end, &
              & press_outer_end, press_inner_end, &
              & Time_Boiler_p(i-1), Time_Boiler_p(i))

        endif

        write(out_lun, 11) i, Time_Boiler_p(i), TEMP_STEAM_p(i), &
          & oxide_thickness(i), &
          & -e_th_scale(1:nFe2O3_pct_check_strain+1) ! * 1.0e+3
 11     format('eth_cycle ', i4, 1x, 20(1pe13.6, 1x)) 

          EVENT_STORE: if (if_output_event)  then

            ! if (k_load_end = id_high)  then
              ! full load: store maximum and miniumum temperatures; one value per time period
              if (k == 1)  then  ! not the height index
                ! ash thickness should be zero; maximum metal temperatures
                ! (1, j_o, j) for full load, ash = 0
                ! (2, j_o, j) for partial load, ash = 0
                Temp_steam_out(k_load_end, j_output, 1) = Temp_steam_now
                Temp_gas_out(k_load_end, j_output, 1) = Temp_gas_now
                Temp_metal_max_out(k_load_end, j_output, 1) = Temp_metal(3)
                Temp_metal_mean_out(k_load_end, j_output, 1) = Temp_metal(2)
                
              endif

              if (no_time_intervals_ash == 0)  then
                ! just assigned them the previous values such that there no bug
                  Temp_steam_out(k_load_end + 2, j_output, 1) = &
                    & Temp_steam_out(k_load_end, j_output, 1)
                  Temp_metal_max_out(k_load_end + 2, j_output, 1) = &
                    & Temp_metal_max_out(k_load_end, j_output, 1)
                  Temp_metal_mean_out(k_load_end + 2, j_output, 1) = &
                    & Temp_metal_mean_out(k_load_end, j_output, 1)

                  ! just for checking
                  Temp_gas_out(k_load_end + 2, j_output, 1) = &
                    & Temp_gas_out(k_load_end, j_output, 1)

              else if (no_time_intervals_ash > 0)  then

                if (k == no_time_intervals_ash)  then   ! not the height index
                ! (3, j_o, 1) for full load, ash = max
                ! (4, j_o, 1) for partial load, ash = max
                  Temp_steam_out(k_load_end + 2, j_output, 1) = Temp_steam_now
                  Temp_metal_max_out(k_load_end + 2, j_output, 1) = Temp_metal(3)
                  Temp_metal_mean_out(k_load_end + 2, j_output, 1) = Temp_metal(2)

                  ! just for checking
                  Temp_gas_out(k_load_end + 2, j_output, 1) = Temp_gas_now

                endif

              endif

          endif EVENT_STORE

       if (if_output_event .and. k_load_end == id_low)  then

          ! not fixed since ash index must be used when storing these not height k index
          write(enrg_lun, 16) j_output, length_tube_oxide_output(k_height), &
           & Temp_gas_out(1, j_output, 1), &
           & Temp_steam_out(1, j_output, 1), Temp_steam_out(3, j_output, 1), &
           & Temp_metal_mean_out(1, j_output, 1), &
           & Temp_metal_mean_out(3, j_output, 1), &
           & Temp_metal_max_out(1, j_output, 1), &
           & Temp_metal_max_out(3, j_output, 1), &
           & Temp_gas_out(2, j_output, 1),  &
           & Temp_steam_out(2, j_output, 1), Temp_steam_out(4, j_output, 1), &
           & Temp_metal_mean_out(2, j_output, 1), &
           & Temp_metal_mean_out(4, j_output, 1), &
           & Temp_metal_max_out(2, j_output, 1), &
           & Temp_metal_max_out(4, j_output, 1), &
           & Temp_gas_out(3, j_output, 1), Temp_gas_out(4, j_output, 1)

        endif

      end do OPER1

    end do HEIGHT_LOOP1

    END DO OPER2

    call STRESS_STRAIN_CREEP_OUT_EXFOL()

    write(enrg_lun, 27) j, j - j, &
      & (INT(time_event(j)), j = 1, no_f2l_event_out)
  27 format('#oxide_vs_height_outage ', i2, 1x, &
      & 20(i6, 1x))
    do k = 1, no_height_oxide_output
       write(enrg_lun, 26) k, length_tube_oxide_output(k), &
        & (1.0e+6*SUM(oxide_thickness_block(k, j, &
        & first_oxide_layer:no_oxide_layers)), &
        & j = 1, no_f2l_event_out)
    end do
  26 format('oxide_vs_height_outage ', i2, 1x, &
      & 20(1pe13.6, 1x))

      ! write out the height dependence; wait till the low load output is stored
      if (if_output_event .and. k_load_end == id_low)  then

        do k = 1, no_height_oxide_output
          write(out_lun, 16) j_output, length_tube_oxide_output(k), &
           & Temp_gas_out(1, j_output, k), &
           & Temp_steam_out(1, j_output, k), Temp_steam_out(3, j_output, k), &
           & Temp_metal_mean_out(1, j_output, k), &
           & Temp_metal_mean_out(3, j_output, k), &
           & Temp_metal_max_out(1, j_output, k), &
           & Temp_metal_max_out(3, j_output, k), &
           & Temp_gas_out(2, j_output, k),  &
           & Temp_steam_out(2, j_output, k), Temp_steam_out(4, j_output, k), &
           & Temp_metal_mean_out(2, j_output, k), &
           & Temp_metal_mean_out(4, j_output, k), &
           & Temp_metal_max_out(2, j_output, k), &
           & Temp_metal_max_out(4, j_output, k), &
           & Temp_gas_out(3, j_output, k), Temp_gas_out(4, j_output, k)
        end do

 16     format('height_output_event ', i1, 1x, 50(1pe13.6, 1x))

      endif
     
  RETURN

  END SUBROUTINE STRESS_LOCATION

 SUBROUTINE SETUP_LOAD_OUT(k_height_id, i, k_load_start, k_load_end, &
          & if_output_event, if_time_output_event, j_output, time_im1)  
    ! time_im1 = Time_Boiler_p(i-1
    !=======================================================================
    ! Purpose(s):
    !
    ! Compute the oxide thickness based on kinetics and temperature schedule 
    ! 
    !=======================================================================
    use parameter_module,   only: oxide_thick2si, id_high, id_low
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: id_temp_op
    use solver_data_module, only: no_thick_interval, if_average_oxide
    use solution_data_module, only: no_oxide_layers, no_full2low_output, &
          & no_point_full2low, no_pulse_full2low_output
    use waterwall_data_module, only: no_height_oxide_output, &
          & oxide_thickness_layer_height, no_thick_interval_now
    use oxide_data_module,  only: first_oxide_layer

    implicit none

    ! Argument List
    integer, INTENT(IN)  :: i, k_height_id
    integer, INTENT(OUT) :: k_load_start, k_load_end, j_output
    logical, INTENT(OUT) :: if_output_event, if_time_output_event
    real, INTENT(IN)    ::  time_im1

    ! Local Variables
    integer :: j, k
    real, dimension(no_height_oxide_output) :: oxide_thickness_out

    if (id_temp_op(i) == 1)  then 
      ! interval 4-1, isothermal low load
      k_load_start = id_low
      k_load_end = id_low
    else if (id_temp_op(i) == 2)  then 
        ! interval 1-2, transition partial - full load
      k_load_start = id_low
      k_load_end = id_high
    else if (id_temp_op(i) == 3)  then 
        ! interval 2-3, isothermal high load
      k_load_start = id_high
      k_load_end = id_high
    else if (id_temp_op(i) == 4)  then 
        ! interval 3-4, transition full - partial load
      k_load_start = id_high
      k_load_end = id_low
    endif

    if_output_event = .false.  ! regular event
    if_time_output_event = .false. ! to save disk space height_event_output_ash

    ! check if this is in an output time
    do j = 1, no_full2low_output

      ! print out full three cycles as it looks better on figures
      if (i >= no_point_full2low(j) -1 .and. i <= no_point_full2low(j) + 7)  then
        ! this starts at [1:2]
        if (no_pulse_full2low_output(j) > 0)  then  ! not an outage
          if_time_output_event = .true. 
        endif

      endif

      ! 2-3 full load loop
      if (i == no_point_full2low(j))  then
 
        if (no_pulse_full2low_output(j) > 0)  then

          ! this is not an outage; 
          ! if id_temp_op(i) = 3 or k_load_end = id_high -> full load

          if_output_event = .true.
          j_output = j

          ! write(6, *) 'output event ', id_temp_op(i-1), id_temp_op(i)
          ! confirmed that id_temp_op(i-1) = 2, id_temp_op(i) = 3

          ! initiate the new oxide thickness to the current one
          do k = 1, no_height_oxide_output
            oxide_thickness_out(k) = &
              & SUM(oxide_thickness_layer_height(first_oxide_layer: &
              & no_oxide_layers, k))
          end do

          ! write(tty_lun, 17) j_output, time_im1, (oxide_thickness_out(k), &
          !    & k = 1, no_height_oxide_output)
          write(out_lun, 17) j_output, time_im1, (oxide_thickness_out(k), &
              & k = 1, no_height_oxide_output)
 17       format('output_event ', i1, 1x, 50(1pe13.6, 1x))

        endif
      endif

      ! 4-1 low load loop; following the previous full load 2-3
      if (i -2 == no_point_full2low(j))  then
        ! if id_temp_op(i) = 1 or k_load_end = id_low -> full load
 
        if (no_pulse_full2low_output(j) > 0)  then
          ! this is not an outage 
          if_output_event = .true.
          j_output = j
        endif
      endif

    end do

 RETURN

 END SUBROUTINE SETUP_LOAD_OUT

  SUBROUTINE ISOTHERM_STATE_STRAIN(k_height_id, i, if_outage_ramp,  &
        & thickness_layer, Temp_gas, Temp_steam, &
        & press_outer, press_inner, Time_i, Delta_Time_i)
    ! Delta_Time_i = Time_i - Time_(i-1)
    !=======================================================================
    ! Purpose(s):
    !
    !  obtain strains, stress, at 2-3 and 4-1 isotherms
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & gen_lun
    use boiler_data_module, only: oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & thickness_fr_var, no_pulse_per_cycle
    use oxide_data_module,  only: temp_kelvin_unit, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, first_oxide_layer
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, no_out_f2l, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use oxide_data_module,   only: id_fe2o3, id_fe3o4, no_oxide_ave, &
           & thickness_fr
    use boiler_data_module, only: id_temp_op
    use solution_driver_module, only: DRIVER_SOLUTION
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & strain_hoop_gen_ave
    use solution_data_module, only: no_oxide_layers
    use solver_data_module, only: if_average_oxide
    use oxide_utility_module, only: OUTPUT_STRAIN_STRESS_ONE, &
      & OUTPUT_STRAIN_STRESS_INT
    use oxide_stress_data_module,   only: total_oxide_thickness_ref, &
      & d_oxide_thickness_4_strain_hoop
    use creep_data_module, only: if_creep_ox_scale, if_creep_metal, &
      & interval_op_now, no_intervals_now, operation_id_current
    use waterwall_data_module, only:  height_id_stress_out 
    use thermal_expansion_module, only: GET_HEAT_FLUX_ALL

    ! Argument
    integer, Intent(IN)  :: i, k_height_id
    logical, Intent(IN)  :: if_outage_ramp
    real, Intent(IN)     :: Temp_gas, Temp_steam, press_outer, press_inner, &
              & Time_i, Delta_Time_i
    real, dimension(moxide), Intent(IN) :: thickness_layer

    ! Local Variables
    integer  :: j, k, step_cycle, nstart
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    real                    :: heat_flux

      !        2------3
      !       /        \
      !      /          \
      !     /            \
      ! -- 1              4------

      d_oxide_thickness_4_strain_hoop = oxide_thickness(i) - oxide_thickness(i-1)

      if (operation_id_current == 1 .or. operation_id_current == 5)  then

        ! end of idle and beginning of a new cycle ramp up

  ! interval_op_now = current interval now considered for the current segment
        interval_op_now = no_intervals_now(operation_id_current)
                
        ! obtain temperature, stress, strain, dimensions -> high temp
        ! if dt = t(i) - t(i-1) would represent the entire partial load
        call DRIVER_SOLUTION(no_oxide_layers, i, &
              & heat_flux, Temp_gas, Temp_steam, &
              & press_outer, press_inner, thickness_layer, &
              & Time_i, Delta_Time_i, id_low, if_average_oxide) 

        ! no need to print out at the end of the partial load

      else if (operation_id_current == 3 .or. operation_id_current == 7)  then

        ! end of pulse and beginning of a new ramp-down segment

  ! interval_op_now = current interval now considered for the current segment
        interval_op_now = no_intervals_now(operation_id_current)

        ! obtain temperature, stress, strain, dimensions -> high temp
        ! if dt = t(i) - t(i-1) would represent the entire full load
        call DRIVER_SOLUTION(no_oxide_layers, i, &
              & heat_flux, Temp_gas, Temp_steam, &
              & press_outer, press_inner, thickness_layer, &
              & Time_i, Delta_Time_i, id_high, if_average_oxide) 

        ! if (k_height_id == height_id_stress_out)  then
          ! print output only for one height
          ! write min, max, average for all variables
          call OUTPUT_STRAIN_STRESS_ONE(i, thickness_layer, id_high, &
            & no_oxide_layers, if_average_oxide)

          call OUTPUT_STRAIN_STRESS_INT(i, thickness_layer, id_high, &
            & no_oxide_layers, if_average_oxide, if_outage_ramp)
        ! endif

      endif

  RETURN

  END SUBROUTINE ISOTHERM_STATE_STRAIN

  SUBROUTINE TRANSITION_STATE_STRAIN(k_height_id, i, if_outage_ramp, &
              & thickness_layer, Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Temp_gas_end, Temp_steam_end, &
              & press_outer_end, press_inner_end, &
              & Time_start, Time_end)
    !=======================================================================
    ! Purpose(s):
    !
    !  obtain strains, stress, at 3-4 and 1-2 transitions (F2L and L2F)
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & gen_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p, thickness_fr_var, &
          & no_pulse_per_cycle
    use oxide_data_module,  only: temp_kelvin_unit, e_th, e_th_mat, &
          & id_mat, no_layer, material_name, first_oxide_layer
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, no_out_f2l, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use oxide_data_module,   only: id_fe2o3, id_fe3o4, no_oxide_ave, &
           & thickness_fr
    use boiler_data_module, only: id_temp_op
    use solution_driver_module, only: DRIVER_SOLUTION
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & strain_hoop_gen_ave
    use solution_data_module, only: no_oxide_layers
    use solver_data_module, only: if_average_oxide, no_intervals_f2l, &
           & no_intervals_l2f
    use oxide_test_module, only: CHECK_EXACT_STRAIN_CYL
    use oxide_utility_module, only: COOLING_STRAIN_PLANE, &
      & OUTPUT_STRAIN_STRESS_ONE, OUTPUT_STRAIN_STRESS_INT, &
      & OUTPUT_STRAIN_STRESS_F2L
    use oxide_stress_data_module,   only: total_oxide_thickness_ref, &
      & d_oxide_thickness_4_strain_hoop
    use creep_data_module, only: if_creep_ox_scale, if_creep_metal, &
      & interval_op_now, operation_id_current
    use waterwall_data_module, only:  height_id_stress_out 
    use thermal_expansion_module, only: GET_HEAT_FLUX_ALL

    ! Argument
    integer, Intent(IN)     :: i, k_height_id
    logical, Intent(INOUT)  :: if_outage_ramp
    real, Intent(IN)        :: Temp_gas_start, Temp_steam_start, &
              & press_outer_start, press_inner_start, &
              & Temp_gas_end, Temp_steam_end, &
              & press_outer_end, press_inner_end, &
              & Time_start, Time_end
    real, dimension(moxide), Intent(IN) :: thickness_layer

    ! Local Variables
    integer  :: j, k, step_cycle, nstart, no_cycle
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
         & Press_gas_now, Press_steam_now, slope_time, &
         & slope_temp_gas, slope_temp_steam, slope_press_gas, &
         & slope_press_steam, time_now, heat_flux
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int

      !   1, 2, 3, 4 are for daily cycle 1-8 are for weekly cycle
      !
      !        2------3             6------7
      !       /        \           /        \
      !      /          \         /          \
      !     /            \       /            \
      ! -- 1              4-----5              8------
      
     if (operation_id_current == 2 .or. operation_id_current == 6)  then

        ! end of ramp-up and beginning of a new pulse segment

        if (no_intervals_l2f == 0 .or. i <= 8)  then
          interval_op_now = 0
          RETURN
        endif

        no_cycle = TEMP_BOILER_data_no + 10  ! to avoid excessive printout

        slope_time = (Time_end - Time_start) / &
                  & no_intervals_l2f
        slope_temp_gas = (Temp_gas_end - Temp_gas_start) / &
                  &  no_intervals_l2f
        slope_temp_steam = (Temp_steam_end - Temp_steam_start) / &
                  &  no_intervals_l2f
        slope_press_gas = (press_outer_end - press_outer_start) / &
                  & no_intervals_l2f
        slope_press_steam = (press_inner_end - press_inner_start) / &
                  & no_intervals_l2f  
        d_oxide_thickness_4_strain_hoop = (oxide_thickness(i) - oxide_thickness(i-1)) / &
                  & no_intervals_l2f  

        ! print out the stress strain state during ramp down
        ! this processing is done at very few cycles; so the 
        ! creep accumulations will not be that great; but keep in mind
        ! that creep here is counted at least twice during this ramp down
        do k = 1, no_intervals_l2f

  ! interval_op_now = current interval now considered for the current segment
          interval_op_now = k

          ! update temperature profile; keep the same radius
          time_now = Time_start + slope_time * k
          temp_gas_now = Temp_gas_start + slope_temp_gas * k
          temp_steam_now = Temp_steam_start + slope_temp_steam * k
          press_gas_now = press_outer_start + slope_press_gas * k
          press_steam_now = press_inner_start + slope_press_steam * k

          if (k == no_intervals_l2f)  no_cycle = i
          ! consider each state to be the low temperature state
          ! obtain temperature, stress, strain, dimensions -> like-low temp

          call DRIVER_SOLUTION(no_oxide_layers, no_cycle, &
                 & heat_flux, temp_gas_now, temp_steam_now, &
                 & press_gas_now, press_steam_now, thickness_layer, &
                 & time_now, slope_time, id_low, if_average_oxide) 

        end do 

        ! no need to print out at the end of the transition from partial to full load

      else if (operation_id_current == 4 .or. operation_id_current == 8)  then

        ! end of ramp-down and beginning of a new idle segment

        no_cycle = TEMP_BOILER_data_no+10  ! to avoid excessive printout

        slope_time = (Time_end - Time_start) / &
                  & no_intervals_f2l
        slope_temp_gas = (Temp_gas_end - Temp_gas_start) / &
                  & no_intervals_f2l
        slope_temp_steam = (Temp_steam_end - Temp_steam_start) / &
                  & no_intervals_f2l
        slope_press_gas = (press_outer_end - press_outer_start) / &
                  & no_intervals_f2l
        slope_press_steam = (press_inner_end - press_inner_start) / &
                  & no_intervals_f2l  

        ! print out the stress strain state during ramp down
        ! this processing is done at very few cycles; so the 
        ! creep accumulations will not be that great; but keep in mind
        ! that creep here is counted at least twice during this ramp down
        do k = 1, no_intervals_f2l

  ! interval_op_now = current interval now considered for the current segment
          interval_op_now = k

          ! update temperature profile; keep the same radius
          time_now = Time_start + slope_time * k
          temp_gas_now = Temp_gas_start + slope_temp_gas * k
          temp_steam_now = Temp_steam_start + slope_temp_steam * k
          press_gas_now = press_outer_start + slope_press_gas * k
          press_steam_now = press_inner_start + slope_press_steam * k

          if (k == no_intervals_f2l)  no_cycle = i
  
          d_oxide_thickness_4_strain_hoop = (oxide_thickness(i) - oxide_thickness(i-1)) / &
                  & no_intervals_f2l  
          ! consider each state to be the low temperature state
          ! obtain temperature, stress, strain, dimensions -> like-low temp

          call DRIVER_SOLUTION(no_oxide_layers, no_cycle, &
                 & heat_flux, temp_gas_now, temp_steam_now, &
                 & press_gas_now, press_steam_now, thickness_layer, &
                 & time_now, slope_time, id_low, if_average_oxide) 

        end do 

        if (k_height_id == height_id_stress_out)  then
          ! print output only for one height
          ! write min, max, average for all variables
          call OUTPUT_STRAIN_STRESS_ONE(i, thickness_layer, id_low, &
            & no_oxide_layers, if_average_oxide)

          ! if_outage_ramp not changed here
          call OUTPUT_STRAIN_STRESS_INT(i, thickness_layer, id_low, &
            & no_oxide_layers, if_average_oxide, if_outage_ramp)

        endif

        ! reset outage flag after printing out 
        if (if_outage_ramp)  if_outage_ramp = .false.

      endif

  RETURN

  END SUBROUTINE TRANSITION_STATE_STRAIN

END MODULE STRESS_DRIVER_MODULE
