MODULE TRANSIT_OUT_MODULE

    !=======================================================================
    ! Purpose(s):
    !
    !   output or storage for shutdown or full-to-low load transition
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: F2L_TRANSITION_OUT, CHECK_CYCLE_STRAIN

 CONTAINS


  SUBROUTINE F2L_TRANSITION_OUT(k_height, i, id_next_outage,  &
      & if_outage_ramp, if_outage_ramp_exfol, thickness_layer, &
      & delta_oxide_thickness)
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
    use solution_driver_module, only: DRIVER_SOLUTION, &
           & DRIVER_SOLUTION_SIMPLE
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, eps0, &
      & no_pulse_full2low_output, &
      & no_full2low_output, no_point_full2low, time_full2low, &
      & time_full2low, temp_gas_full2low, &
      & temp_steam_full2low, press_gas_full2low, press_steam_full2low, &
      & strain_hoop_gen_ave, temp_oxide_high_debug
    use solution_data_module, only: no_oxide_layers, Temp, &
      & temp_gas_high_load, temp_steam_high_load
    use solver_data_module, only: if_average_oxide, no_intervals_f2l, &
           & no_intervals_l2f
    use oxide_test_module, only: CHECK_EXACT_STRAIN_CYL
    use oxide_utility_module, only: COOLING_STRAIN_PLANE, &
      & OUTPUT_STRAIN_STRESS_F2L
    use oxide_stress_data_module,   only: total_oxide_thickness_ref, &
      & d_oxide_thickness_4_strain_hoop
    use creep_data_module, only: if_creep_ox_scale, if_creep_metal
    use boiler_shutdown_module, only:  Temp_Shut_Down
    use shutdown_data_module, only:  time_outage_period, &
         & duration_outage_period, delta_time_outage_period, &
         & temp_steam_outage_period, number_outage_periods, &
         & no_out_f2l_period, no_out_f2l_total, time_now_f2l, &
         & temp_steam_now_f2l, temp_gas_now_f2l, &
         & press_steam_now_f2l, press_gas_now_f2l, no_out_f2l_total
    use thermal_expansion_module, only: GET_HEAT_FLUX_ALL
    use waterwall_data_module, only: no_height_oxide_output, &
          & Temp_growth_layer_height_old, flow_rate_steam_high, &
          & flow_rate_steam_low, Temp_steam_height, Temp_gas_height

    ! Argument
    integer, Intent(IN)  :: i, k_height, id_next_outage
    logical, Intent(INOUT)  :: if_outage_ramp
    logical, Intent(IN)  :: if_outage_ramp_exfol
    real, dimension(moxide), Intent(IN) :: thickness_layer
    real, intent(IN)  :: delta_oxide_thickness  ! for stress induced oxide growth

    ! Local Variables
    integer  :: j, k, step_cycle, nstart, no_cycle, kper, count1
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now, heat_flux_now, &
              & Temp_start1

    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int

    ! temporary fix for shutdown
    real    :: duration_outage_ramp_down, time_relative_now
    logical, save :: if_long_shutdown_first = .true.
    logical  :: if_const_creep_or_oig_trans = .false. ! .true.
    
    ! in order not to duplicate the ramp down, also during a shut down this would not be ok
    d_oxide_thickness_4_strain_hoop = 0.0

        ! print out the stress strain state during ramp down
        ! this processing is done at very few cycles; so the 
        ! creep accumulations will not be that great; but keep in mind
        ! that creep here is counted at least twice during this ramp down
        do j = 1, no_full2low_output

          if (i == no_point_full2low(j))  then

            ! store pressures at full load at oxide metal interface
            press_full_int(j) = stress_rad(2, 1)

           ! done with the paper output

            if (no_pulse_full2low_output(j) > 0)  then

              write(tty_lun, 3) j, no_pulse_full2low_output(j), &
                & i, oxide_thickness(i)

              if (gen_lun > 0)  then
                write(gen_lun, 3) j, no_pulse_full2low_output(j), &
                  & i, oxide_thickness(i)

              endif

 3            format('f2l_cycle ',i2, ' cycle= ', i5, 1x, &
                & ' no_point ', i5, 1x, ' oxide_thickness= ', 1pe13.6)

            else if (no_pulse_full2low_output(j) == 0)  then

              write(tty_lun, 4) j, no_pulse_full2low_output(j), &
                & i, oxide_thickness(i)
              if (gen_lun > 0)  then
                write(gen_lun, 4) j, no_pulse_full2low_output(j), &
                  & i, oxide_thickness(i)
              endif
 4            format('f2l_outage ',i2, ' cycle= ', i5, 1x, &
                & ' no_point ', i5, 1x, ' oxide_thickness= ', 1pe13.6)

              ! this flag will be used at the end of ramp to print the highest 
              ! elastic energy during the outage transition, instead of the 
              ! value at the end of the outage
              if_outage_ramp = .true.

            endif

            ! initialize the low load quantities
            ! time_full2low(j, no_thick_interval) = time_full2low(j, mramp)
            ! temp_gas_full2low(j, no_thick_interval) = &
            !       & temp_gas_full2low(j, mramp)

              temp_gas_high_load = Temp_gas_height(k_height, id_high)
              temp_steam_high_load = Temp_steam_height(k_height, id_high)

            slope_time = (time_full2low(j, mramp) - &
                    & time_full2low(j, 1)) / no_out_f2l
              slope_temp_gas = (temp_gas_full2low(j, mramp) - &
                  & temp_gas_high_load) / no_out_f2l
              slope_temp_steam = (temp_steam_full2low(j, mramp) - &
                  & temp_steam_high_load) / no_out_f2l
            slope_press_gas = (press_gas_full2low(j, mramp) - &
                  & press_gas_full2low(j, 1)) / no_out_f2l
            slope_press_steam = (press_steam_full2low(j, mramp) - &
                  & press_steam_full2low(j, 1)) / no_out_f2l   
   
            time_now = 0.0

            temp_gas_now = temp_gas_high_load
            temp_steam_now = temp_steam_high_load

            ! store these for the output
            time_now_f2l(1) = time_now
            temp_steam_now_f2l(1) = temp_steam_high_load
            temp_gas_now_f2l(1) = temp_gas_high_load
            press_steam_now_f2l(1) = press_steam_full2low(j, 1)
            press_gas_now_f2l(1) = press_gas_full2low(j, 1)

                ! this is used for debugging purposes
                if (id_next_outage == 1 .and. k_height == 1)  then

                  temp_oxide_high_debug = &
                     & (Temp(no_oxide_layers, npr_temp(no_oxide_layers)) + &
                     & Temp(1, npr_temp(1))) / 2.0 

                endif                

            ! write the output for the transition full to low load
            call OUTPUT_STRAIN_STRESS_F2L(i, j, k_height, id_next_outage, &
               & thickness_layer, time_now, &
               & temp_gas_now, temp_steam_now, no_oxide_layers, &
               & if_average_oxide, if_outage_ramp_exfol, .true.)  ! new output

            duration_outage_ramp_down = time_full2low(j, mramp) - time_full2low(j, 1)

            OUTAGE_PR: if (duration_outage_ramp_down > 90.0)  then

              ! temporary code until input data for new outage is written
              ! new outage
              ! steam pressure drops to 1 atm in 3-4 h (very fast)
              ! gas temperature drops to 300 C in 3-4 h (very fast)
              ! metal pressure drops to 100 C in about 100 h total
          
              ! initialize relative time
              time_now = 0.0

              number_outage_periods = 4             
              ! 3 superheater 3 for reheater
              ! include a dummy fourth step to RT

              ! first period 
              ! 5 superheater; 2 h reheater; keep it to 3 to correlated with Tg
              time_outage_period(1) = 0.0          
              time_outage_period(2) = 3.0           
              ! duration_outage_period(1) = 3.0           ! 5 superheater; 2 h reheater
              ! time response formulation is valid after 5 h but the 
              ! Tg drops to 300 in 3 h and there is no data at 4 h; continuity
              ! solution; prescribe the Tm directly at 4 
              delta_time_outage_period(1) = 1.0/4.0 ! 15 minutes
              temp_steam_outage_period(2) = 414.0   ! 414.0  superheater 485 reheater

              ! second period
              ! here only the tube temperature varies; the pressure is atmospheric
              time_outage_period(3) = 5.0           ! 5 superheater 31 h reheater
              ! duration_outage_period(2) = 2.0           ! 2 superheater 28 h reheater
              delta_time_outage_period(2) = 1.0/4.0 ! 15 min
              temp_steam_outage_period(3) = 353.0   ! 353.0  superheater 472 reheater

              ! third period 
              time_outage_period(4) = 100.0
              ! duration_outage_period(3) = 100.0  ! 100 superheater 70.0 reheater
              ! here only the tube temperature varies; the pressure is atmospheric
              ! duration_outage_period(3) = duration_outage_ramp_down - time_outage_period(1) - &
              !    & time_outage_period(2) 
              ! 100.0 h for superheater
              ! 28 h for reheater

              delta_time_outage_period(3) = 2.0  !2 superheater; 2 reheater

              ! not used yet since exponential decay
              temp_steam_outage_period(4) = 100.0   ! 100.0  superheater 100 reheater

              ! dummy step only for cooling to RT
              time_outage_period(5) = 112.0          
              delta_time_outage_period(4) = 2.0
              temp_steam_outage_period(4) = 25.0

              do kper = 1, number_outage_periods
                duration_outage_period(kper) = time_outage_period(kper+1) - &
                      & time_outage_period(kper)
              end do

              where (delta_time_outage_period > 0.001)
                no_out_f2l_period = INT(duration_outage_period / &
                 & delta_time_outage_period)
              end where

              ! initialize the count of outputs
              ! no_out_f2l_total = SUM(no_out_f2l_period) + 1
              no_out_f2l_total = 1
              do kper = 1, number_outage_periods
                no_out_f2l_total = no_out_f2l_total + no_out_f2l_period(kper)
              end do

              count1 = 1

              if (if_long_shutdown_first)  then

                kper = 1
                k = 1
                write(2, 53) 
  53            format('long_shut_down j kp k time Tst Tox-m Tg pst pg')
                write(2, 51) j, kper, k, time_now, temp_steam_now, &
                  & Temp(1, npr_temp(1)), temp_gas_now, press_steam_now, &
                  & press_gas_now
  51            format('long_shut_down ', i2, 1x, i1, 1x, i3, 1x, 20(1pe13.6, 1x))

              endif

              do kper = 1, number_outage_periods

                write(out_lun, 52) kper, no_out_f2l_period(kper), &
                  & no_out_f2l_total, time_outage_period(kper), &
                  & time_outage_period(kper+1), &
                  & duration_outage_period(kper)
  52            format('long_shut_down1 ', i1, 1x, 2(i3, 1x), 20(1pe13.6, 1x))

                do k = 2, no_out_f2l_period(kper) + 1

                  ! current output count
                  count1 = count1 + 1

                  time_now = time_now + delta_time_outage_period(kper)
                  ! time_now = time_outage_period(kper) + delta_time_outage_period(kper)

                  if (kper == 1)  then

                    ! Tg drops to 300 C in 3 h and the htc_out may not be usable
                    ! actually it is kind of low any way, in this case
                    ! the steam temperature would take over
                    temp_gas_now = temp_gas_high_load + &
                      & (300.0 - temp_gas_high_load) * &
                      & time_now / duration_outage_period(kper)

                    ! sh 414 C; for reheater 485 C
                    temp_steam_now = temp_steam_high_load + &
                      & (temp_steam_outage_period(kper+1) - &
                      & temp_steam_high_load) * &
                      & time_now / duration_outage_period(kper)

                    ! the entire pressure drop occurs in the first period
                    press_gas_now = press_gas_full2low(j, 1) + &
              & (press_gas_full2low(j, mramp) - press_gas_full2low(j, 1)) * &
                      & time_now / duration_outage_period(kper)

                    ! the entire pressure drop occurs in the first period
                    press_steam_now = press_steam_full2low(j, 1) + &
              & (press_steam_full2low(j, mramp) - press_steam_full2low(j, 1)) * &
                      & time_now / duration_outage_period(kper)

                  else if (kper > 1)  then

                    ! assign the same temperature in/out since htc is not available
                    if (kper == 2)  then

                      ! linear variation; 414-353 superheater; 485-472 reheater
                      temp_steam_now = temp_steam_outage_period(kper) + &
                        & (temp_steam_outage_period(kper+1) - &
                        &  temp_steam_outage_period(kper)) * &
                        & (time_now - time_outage_period(kper)) / &
                        &  duration_outage_period(kper)

                    else if (kper == 3)  then

                      ! exponential decay variation
                      time_relative_now = time_now - time_outage_period(kper)
                      temp_steam_now = Temp_Shut_Down(1, &
                         & temp_steam_outage_period(kper), time_relative_now)

                    else if (kper == 4)  then
                      ! linear variation to RT
                      if (k == no_out_f2l_period(kper) + 1)  &
                         & time_now = time_outage_period(kper+1)
                      time_relative_now = (time_now - time_outage_period(kper)) / &
                        & (time_outage_period(kper+1)-time_outage_period(kper))
                      Temp_start1 = Temp_Shut_Down(1, &
                         & temp_steam_outage_period(kper-1), &
                         & time_outage_period(kper)-time_outage_period(kper-1))
                      temp_steam_now = Temp_start1 * (1.0-time_relative_now) + &
                         & temp_steam_outage_period(kper) * time_relative_now

                    endif

                    temp_gas_now = temp_steam_now
                    ! assign the same pressure in/out
                    ! check main stream pressure to be used ? 
                    !    this pressure drops in 5-6 h asymptotically 
                    ! or first stage pressure is used; drops in 2 h linearly
                    press_steam_now = press_steam_full2low(j, mramp) 
                    press_gas_now = press_steam_full2low(j, mramp) 

                  endif

                  ! store time_now for output
                  time_now_f2l(count1) = time_now
                  temp_steam_now_f2l(count1) = temp_steam_now
                  temp_gas_now_f2l(count1) = temp_gas_now
                  press_steam_now_f2l(count1) = press_steam_now
                  press_gas_now_f2l(count1) = press_gas_now

                  ! heat_flux_now = 
                  call GET_HEAT_FLUX_ALL(no_oxide_layers, &
                     & temp_gas_now, temp_steam_now, heat_flux_now)

                ! consider each state to be the low temperature state
                ! obtain temperature, stress, strain, dimensions -> like-low temp
! consider that during an outage the oxide induced growht does not vary since the oxide does not grow as much; for long outage we do not have the data

                if (if_const_creep_or_oig_trans) then

                  call DRIVER_SOLUTION_SIMPLE(no_oxide_layers, &
                  & TEMP_BOILER_data_no+10, &
                  & heat_flux_now, temp_gas_now, temp_steam_now, &
                  & press_gas_now, press_steam_now, thickness_layer, &
                  & time_now + Time_Boiler_p(i), slope_time, &
                  & id_low, if_average_oxide) 

                else

                  call DRIVER_SOLUTION(no_oxide_layers, TEMP_BOILER_data_no+10, &
                  & heat_flux_now, temp_gas_now, temp_steam_now, &
                  & press_gas_now, press_steam_now, thickness_layer, &
                  & time_now + Time_Boiler_p(i), slope_time, &
                  & id_low, if_average_oxide) 

                endif

                  if (if_long_shutdown_first)  then

                    write(2, 51) j, kper, k, time_now, temp_steam_now, &
                     & Temp(1, npr_temp(1)), temp_gas_now, press_steam_now, &
                     & press_gas_now

                  endif

                ! write the output for the transition full to low load
                call OUTPUT_STRAIN_STRESS_F2L(i, j, k_height, id_next_outage, &
                  & thickness_layer, time_now, &
                  & temp_gas_now, temp_steam_now, no_oxide_layers, &
                  & if_average_oxide, if_outage_ramp_exfol, .true.)  ! new output

              end do

              end do

              if (if_long_shutdown_first)  if_long_shutdown_first = .false.

             else if (duration_outage_ramp_down < 90.0)  then

              ! old outage; 
              ! steam temperature, gas temperature, and steam pressure all vary
              ! linear in time and duration for the variation of all them is the same
              ! also the htc's are the same
              no_out_f2l_total = no_out_f2l

              slope_temp_gas = (temp_gas_full2low(j, mramp) - &
                  & temp_gas_high_load) / no_out_f2l
              slope_temp_steam = (temp_steam_full2low(j, mramp) - &
                  & temp_steam_high_load) / no_out_f2l

              do k = 2, no_out_f2l + 1

                ! get temperature profile; keep the same radius

                time_now = slope_time * (k - 1)

                temp_gas_now = temp_gas_high_load + &
                    & slope_temp_gas * (k - 1)
                temp_steam_now = temp_steam_high_load + &
                    & slope_temp_steam * (k - 1)

                press_gas_now = press_gas_full2low(j, 1) + &
                    & slope_press_gas * (k - 1)
                press_steam_now = press_steam_full2low(j, 1) + &
                    & slope_press_steam * (k - 1)

                ! heat_flux_now = 
                call GET_HEAT_FLUX_ALL(no_oxide_layers, &
                     & temp_gas_now, temp_steam_now, heat_flux_now)

                ! consider each state to be the low temperature state
                ! obtain temperature, stress, strain, dimensions -> like-low temp
 ! consider that during an outage the oxide induced growht does not vary since the oxide does not grow as much; for long outage we do not have the data
  
                if (if_const_creep_or_oig_trans) then

                  call DRIVER_SOLUTION_SIMPLE(no_oxide_layers, &
                  & TEMP_BOILER_data_no+10, &
                  & heat_flux_now, temp_gas_now, temp_steam_now, &
                  & press_gas_now, press_steam_now, thickness_layer, &
                  & time_now + Time_Boiler_p(i), slope_time, &
                  & id_low, if_average_oxide) 

                else

                  call DRIVER_SOLUTION(no_oxide_layers, TEMP_BOILER_data_no+10, &
                  & heat_flux_now, temp_gas_now, temp_steam_now, &
                  & press_gas_now, press_steam_now, thickness_layer, &
                  & time_now + Time_Boiler_p(i), slope_time, &
                  & id_low, if_average_oxide) 

                endif

                  time_now_f2l(k) = time_now
                  temp_steam_now_f2l(k) = temp_steam_now
                  temp_gas_now_f2l(k) = temp_gas_now
                  press_steam_now_f2l(k) = press_steam_now
                  press_gas_now_f2l(k) = press_gas_now

                ! write the output for the transition full to low load
                call OUTPUT_STRAIN_STRESS_F2L(i, j, k_height, id_next_outage, &
                  & thickness_layer, time_now, &
                  & temp_gas_now, temp_steam_now, no_oxide_layers, &
                  & if_average_oxide, if_outage_ramp_exfol, .true.)  ! new output

              end do

            endif OUTAGE_PR

            ! write pressures, full load, low load at oxide metal interface
            write(aux_lun, 23) j, oxide_thickness(i), &
              & -press_full_int(j), -stress_rad(first_oxide_layer, 1)
 23     format('full_low_load_press ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe10.3, 1x))

            write(aux_lun, 24) j, oxide_thickness(i), &
              & strain_th_plate_ave(first_oxide_layer, 1)*ten_p3, &
              & strain_th_plate_ave(first_oxide_layer, 2)*ten_p3, &
              & strain_th_plate_ave(first_oxide_layer, 3)*ten_p3,  &
              & stress_hoop_th_ave(first_oxide_layer, 1), &
              & stress_hoop_th_ave(first_oxide_layer, 3), &
              & stress_hoop_th_ave(first_oxide_layer, 2)
 24     format('low_load_hoop_plate ', i4, 1x, 1(1pe13.6, 1x), &
                 & 3(1pe11.4, 1x), 1x, 3(1pe10.3, 1x))

          ! endif NO_CREEP2

          endif

        end do

  RETURN

  END SUBROUTINE F2L_TRANSITION_OUT

  SUBROUTINE CHECK_CYCLE_STRAIN()
    !=======================================================================
    ! Purpose(s):
    !
    !  check the strains at different temperatures and %Fe2O3 
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide
    use output_module,      only: tty_lun, out_lun, aux_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse
    use oxide_data_module,  only: temp_kelvin_unit, e_th, e_th_mat, &
          & id_mat, no_layer, material_name
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use oxide_data_module,   only: id_fe2o3, id_fe3o4, no_oxide_ave
    use boiler_data_module, only: id_temp_op
    use oxide_utility_module, only: COOLING_STRAIN_PLANE

    ! Local Variables
    integer  :: i, j, k
    real     :: Temp_High, Temp_Low
    real, dimension(3)     :: fr, e_check, energy
    real, dimension(mstep)     :: e_ave, energy_ave
    character(LEN = 80)    :: mode_stress = 'biaxial_plane'
    real, dimension(moxide) :: e_th_scale

    if (ntemp_check_strain == 0 .or. nFe2O3_pct_check_strain == 0) return

    write(tty_lun, *) 'Fe2O3',  id_fe2o3, 'Fe3O4',  id_fe3o4

    OPER1: do i= 2, TEMP_BOILER_data_no ! no of data intervals

      if (oxide_thickness(i) < small_oxide_thickness)  cycle OPER1

      if (id_temp_op(i) == 2 .or. id_temp_op(i) == 3)  then

        e_th_scale = 0.0
        ! sigma_th_scale = 0.0

      else if (id_temp_op(i) == 4)  then

        ! realocate the 
        j = no_oxide_ave
        e_th_scale(1) = COOLING_STRAIN_PLANE(id_mat(j), &
         & TEMP_STEAM_p(i-1), TEMP_STEAM_p(i))

        do k= 1, nFe2O3_pct_check_strain
          j = no_oxide_ave + k
          e_th_scale(k+1) = COOLING_STRAIN_PLANE(id_mat(j), &
            & TEMP_STEAM_p(i-1), TEMP_STEAM_p(i))
        end do

      else if (id_temp_op(i) == 1)  then

        e_th_scale = e_th_scale
        ! sigma_th_scale = sigma_th_scale

      endif

      write(out_lun, 11) i, Time_Boiler_p(i), TEMP_STEAM_p(i), &
       & oxide_thickness(i), &
       & -e_th_scale(1:nFe2O3_pct_check_strain+1) ! * 1.0e+3
 11   format('ecth_cycle ', i4, 1x, 20(1pe13.6, 1x)) 

   end do OPER1

   return

  END SUBROUTINE CHECK_CYCLE_STRAIN
    
END MODULE TRANSIT_OUT_MODULE
