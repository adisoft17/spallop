MODULE STEAM_TEMP_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   steam temperature as a function of tube height, 
    !    heat flux is intermediary
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: GET_HEAT_FLUX_ST_TEMP, GET_HEAT_FLUX_ST_TEMP_OX, &
    & TEMP_OUT_STEAM_OXIDE_METAL

 CONTAINS


      !   1, 2, 3, 4 are for daily cycle 1-8 are for weekly cycle
      !
      !        2------3 
      !       /        \ 
      !      /          \
      !     /            \
      ! -- 1              4
    
  SUBROUTINE GET_HEAT_FLUX_ST_TEMP()
    !=======================================================================
    ! Purpose(s):
    !
    !  get the flue gas temperature and the htc inner, these are coupled 
    !  together as a function of height, 
    !  also the steam temperature is obtained
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & gen_lun, enrg_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p, thickness_fr_var, &
          & no_pulse_per_cycle
   use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer, temp_gas_pulse, &
      & temp_steam_pulse, temp_gas_idle, temp_steam_idle, &
      & heat_flux_fraction_idle
   use oxide_data_module,  only: cond_value

   use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, no_out_f2l, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use oxide_data_module,   only: cond_value, &
           & thickness_fr, first_oxide_layer
    use boiler_data_module, only: id_temp_op
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
    use waterwall_data_module,  only: duration_ash_deslaging, &
       & no_time_intervals_ash, if_constant_htc_inner, &
       & ratio_time_intervals_ash, heat_flux_height, length_tube_points, &
       & no_heat_flux_height, if_ash, no_heat_flux_height, &
       & height_boiler, if_ash, Temp_steam_height, &
       & Temp_steam_inlet, flow_rate_steam_high, flow_rate_steam_low, &
       & Temp_gas_height, htc_inner_height, no_height_oxide_output, &
       & length_tube_oxide_output, heat_flux_output, &
       & Delta_Temp_metal_mean_ash_max, ash_thickness_max, &
       & hflux_htc_ash_estimate, htc_tube_outer_ash, k_ash, &
       & heat_flux_output_ash, no_time_intervals_ash
    use waterwall_data_module,  only: Temp_gas_height_input, &
       & if_Tg_or_heat_flux, if_constant_steam_temp, &
       & metal_temperature_model_type, htc_inner_height_input, &
       & htc_outer_height_input, Temp_steam_height_input
    use mesh_module,   only: total_oxide_thickness_ref
    use htc_utility_module, only: GET_ASH_THICKNESS
    use steam_module, only:  Calc_dTSteam_dHeight, STEAM_PIPE_EXPON
    use mesh_module, only: GET_TUBE_DIMENSIONS

    implicit none

    ! Local Variables
    integer  :: i, j, k, step_cycle, nstart, convergent, max_iter_htc_temp
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now
    real, dimension(3)     :: fr
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    real, dimension(2)  :: din, dout, htc_tube_inner_now, &
       & htc_tube_outer_now, htc_tube_inner_old
    real :: rad_out, rad_in, temp_max_full, temp_max_idle, a1, &
       & temp_ave_full, temp_ave_idle, temp_min_full, temp_min_idle, &
       & delta_temp_full, height_now, Temp_wall_old, Temp_wall_new, &
       & dhtc_tube_inner_dT_wall, delta_Temp_wall, dhtc_tube_inner_dT_wall_old, &
       & fval, df_dT_wall, dTg_now
    real :: HTC_now, dTSteam_dHeight, dHeight, dTSteam, &
       & c2, d2, x_ash, dT_metal_mean, c1, htc_total, Tg_now, beta_exponent
    real, dimension(20, 2)     :: dT_MW_estimate  ! mean wall temperature variation with ash

    if (.not. if_Tg_or_heat_flux)  then
      write(6, *) 'ERROR: Tg must be used at input not heat flux'
      STOP
    endif

    if (.not. if_constant_htc_inner)  then

      write(6, *) 'ERROR: exfoliation_loop only with constant_htc_inner'
      stop

     endif

     ! maximum iterations to iterate for flue gas temperature and htc_in 
     max_iter_htc_temp = 10

     ! find the distribution of the steam temperature when heat flux is known
     ! neglect the oxide inside the tube and neglect the ash deposit outside the tube
     ! id_low = 2, id_high = 1 to keep track of low load and full load
      if (TRIM(metal_temperature_model_type) == 'design')  then

        Temp_steam_height(:, id_high) = Temp_steam_height_input
        Temp_steam_height(:, id_low) = Temp_steam_height_input

      else

        Temp_steam_height(1, id_high) = Temp_steam_inlet
        Temp_steam_height(1, id_low) = Temp_steam_inlet

      endif

      ! need to reduce Tg at the low load
      Temp_gas_height(:, id_high) = &
                     & Temp_gas_height_input(1:no_height_oxide_output)
      Temp_gas_height(:, id_low) = &
                     & Temp_gas_height_input(1:no_height_oxide_output)

      j = 1

      ! at each height determine the steam temperature
      HEIGHT_LOOP1: do j = 2, no_height_oxide_output

        ! set tube_outer_radius and tube_thickness at this location
        call GET_TUBE_DIMENSIONS(j-1)

        ! parameters from the log profile of temperature in the tube without oxide/ash
        rad_out = tube_outer_radius
        rad_in = rad_out - tube_thickness
        a1 = LOG(rad_out / rad_in)

        height_now = length_tube_oxide_output(j)
        dHeight = length_tube_oxide_output(j) - length_tube_oxide_output(j-1)

        if (TRIM(metal_temperature_model_type) == 'design')  then

          htc_tube_outer_now(id_high) =  htc_outer_height_input(j)
          htc_tube_outer_now(id_low) =  htc_outer_height_input(j)

        else

          htc_tube_outer_now(id_high) =  htc_tube_outer
          htc_tube_outer_now(id_low) =  htc_tube_outer

        endif

        ! for each load type: full (k=1) and partial (k=2)
        LOAD_LOOP1: do k = 1, 2

          ! steam temperature is UNKNOWN
          ! gas temperature and htc_inner are KNOWN
          if (TRIM(metal_temperature_model_type) == 'design')  then

            htc_tube_inner_now(k) =  htc_inner_height_input(j)

          else

            htc_tube_inner_now(k) =  htc_tube_inner

          endif

          htc_total = 1.0 / htc_tube_outer_now(k) + &
                    & rad_out / (htc_tube_inner_now(k) * rad_in) + &
                    & (rad_out / cond_value(1, 1)) * LOG(rad_out / rad_in)

          htc_total = 1.0 / htc_total

          beta_exponent = STEAM_PIPE_EXPON(k, Temp_steam_height(j-1,k), &
                & htc_total)

          ! constant gas temprature
          Tg_now = (Temp_gas_height(j-1, k) + Temp_gas_height(j, k)) / 2.0

          dTg_now = Temp_gas_height(j, k) - Temp_gas_height(j - 1, k)

          if (if_constant_steam_temp)  then

            ! constant in time; to be also constant spatially, make Temp_steam_height_input=const.
            dTSteam = 0.0
            if (TRIM(metal_temperature_model_type) == 'design')  then
              ! Temp_steam_height(j,k) was set at input
            else
              Temp_steam_height(j,k) = Temp_steam_height(j-1,k)
            endif

          else
            ! old; before 08/27/2012 
            ! dTSteam = (Tg_now - Temp_steam_height(j-1,k)) * &
            !      & (1.0 - EXP(-beta_exponent * dHeight))

            dTSteam = dTg_now + (Temp_gas_height(j, k) - &
               & Temp_steam_height(j-1, k) - &
               & dTg_now / (dHeight * beta_exponent)) * &
               & (1.0 - EXP(-beta_exponent * dHeight))

            Temp_steam_height(j,k) = Temp_steam_height(j-1,k) + dTSteam

          endif

          ! STEP H1
          ! get the heat flux at this height
          ! use interpolation between tabulated heat_flux_height(length_tube_points)
          !  if length_tube_points is different than length_tube_oxide_output
          ! heat_flux_output(j, k) = heat_flux_height(j) ! ACTION 1
          heat_flux_output(j, k) = htc_total * (Temp_gas_height(j, k) - &
                 & Temp_steam_height(j, k))

          if (j == 2)  then
            ! get the heat flux for j = 1
            heat_flux_output(j-1, k) = htc_total * (Temp_gas_height(j-1, k) - &
                 & Temp_steam_height(j-1, k))
          endif

          ! different heat fluxes are used at low/high load
          ! heat_flux_output(j,id_low)= heat_flux_output(j,id_high) * &
          !    & heat_flux_fraction_idle(1)

          write(out_lun, 1) j, k, height_now, Temp_steam_height(j-1, k), &
              & Temp_steam_height(j, k), Temp_gas_height(j, k), &
              & dTSteam, htc_total, htc_total * (Temp_gas_height(j-1, k) - &
              & Temp_steam_height(j-1, k)), htc_total * (Temp_gas_height(j, k) - &
              & Temp_steam_height(j, k)), Temp_gas_height(j-1, k), tg_now
          write(enrg_lun, 1) j, k, height_now, Temp_steam_height(j-1, k), &
              & Temp_steam_height(j, k), Temp_gas_height(j, k), &
              & dTSteam, htc_total, htc_total * (Temp_gas_height(j-1, k) - &
              & Temp_steam_height(j-1, k)), htc_total * (Temp_gas_height(j, k) - &
              & Temp_steam_height(j, k)), Temp_gas_height(j-1, k), tg_now
 1        format('Tsteam_calc ', i2, 1x, i1, 1x, &
                  &  14(1pe13.6, 1x))

        end do LOAD_LOOP1

      end do HEIGHT_LOOP1

   return

  END SUBROUTINE GET_HEAT_FLUX_ST_TEMP

  SUBROUTINE GET_HEAT_FLUX_ST_TEMP_OX(i_time)
    !=======================================================================
    ! Purpose(s):
    !
    !  get the flue gas temperature and the htc inner, these are coupled 
    !  oxide thickness is that at previous time step
    !  together as a function of height, 
    !  also the steam temperature is obtained
    ! 
    !=======================================================================
    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high, mave, ten_p3, mramp, zero, mramp_out
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & gen_lun, enrg_lun
    use boiler_data_module, only: Time_Boiler_p, TEMP_STEAM_p, &
          & TEMP_STEAM_Slope_p, TEMP_BOILER_data_no, oxide_thickness, &
          & temp_steam_pulse, temp_steam_idle, no_total_pulse, &
          & Temp_gas_p, press_inner_slope_p, press_outer_slope_p, &
          & press_inner_p, press_outer_p, thickness_fr_var, &
          & no_pulse_per_cycle
   use boiler_data_module, only: tube_thickness, tube_outer_radius, &
      & htc_tube_inner, htc_tube_outer, temp_gas_pulse, &
      & temp_steam_pulse, temp_gas_idle, temp_steam_idle, &
      & heat_flux_fraction_idle
   use oxide_data_module,  only: cond_value

   use oxide_data_module,  only: temp_kelvin_unit, oxide_rate_log, &
          & oxide_rate_1temp, nrate, oxide_rate_temp, e_th, e_th_mat, &
          & id_mat, no_layer, material_name
    use solver_data_module, only: no_thick_interval, temp_check_strain, &
           & Fe2O3_pct_check_strain, ntemp_check_strain, no_out_f2l, &
           & nFe2O3_pct_check_strain, small_oxide_thickness
    use oxide_data_module,   only: cond_value, &
           & thickness_fr, first_oxide_layer
    use boiler_data_module, only: id_temp_op
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
    use waterwall_data_module,  only: duration_ash_deslaging, &
       & no_time_intervals_ash, if_constant_htc_inner, &
       & ratio_time_intervals_ash, heat_flux_height, length_tube_points, &
       & no_heat_flux_height, if_ash, no_heat_flux_height, &
       & height_boiler, if_ash, Temp_steam_height, &
       & Temp_steam_inlet, flow_rate_steam_high, flow_rate_steam_low, &
       & Temp_gas_height, htc_inner_height, no_height_oxide_output, &
       & length_tube_oxide_output, heat_flux_output, &
       & Delta_Temp_metal_mean_ash_max, ash_thickness_max, &
       & hflux_htc_ash_estimate, htc_tube_outer_ash, k_ash, &
       & heat_flux_output_ash, no_time_intervals_ash, temp_steam_height_time, &
       & thickness_oxide_height_time, if_temp_steam_depend_oxide
    use waterwall_data_module,  only: Temp_gas_height_input, &
       & if_Tg_or_heat_flux, if_constant_steam_temp, &
       & metal_temperature_model_type, htc_inner_height_input, &
       & htc_outer_height_input, Temp_steam_height_input
    use mesh_module,   only: total_oxide_thickness_ref
    use htc_utility_module, only: GET_ASH_THICKNESS
    use steam_module, only:  Calc_dTSteam_dHeight, STEAM_PIPE_EXPON
    use mesh_module, only: GET_TUBE_DIMENSIONS

    implicit none

    ! Argument List
    integer, intent(IN) :: i_time

    ! Local Variables
    integer  :: i, j, k, step_cycle, nstart, convergent, max_iter_htc_temp
    real     :: Temp_High, Temp_Low, Temp_gas_now, Temp_steam_now, &
              & Press_gas_now, Press_steam_now, slope_time, &
              & slope_temp_gas, slope_temp_steam, slope_press_gas, &
              & slope_press_steam, time_now
    real, dimension(3)     :: fr
    real, dimension(moxide) :: e_th_scale, thickness_layer
    real, dimension(moxide) :: dtemp, dstrain_th, &
          & dstrain_hoop, dstress_hoop
    real, dimension(mramp_out)  :: press_full_int
    real, dimension(2)  :: din, dout, htc_tube_inner_now, &
       & htc_tube_outer_now, htc_tube_inner_old
    real :: rad_out, rad_in, temp_max_full, temp_max_idle, a1, &
       & temp_ave_full, temp_ave_idle, temp_min_full, temp_min_idle, &
       & delta_temp_full, height_now, Temp_wall_old, Temp_wall_new, &
       & dhtc_tube_inner_dT_wall, delta_Temp_wall, dhtc_tube_inner_dT_wall_old, &
       & fval, df_dT_wall, rad_in_ox, dTg_now
    real :: HTC_now, dTSteam_dHeight, dHeight, dTSteam, &
       & c2, d2, x_ash, dT_metal_mean, c1, htc_total, Tg_now, beta_exponent
    real, dimension(20, 2)     :: dT_MW_estimate  ! mean wall temperature variation with ash

    if (.not. if_Tg_or_heat_flux)  then
      write(6, *) 'ERROR: Tg must be used at input not heat flux'
      STOP
    endif

    if (.not. if_constant_htc_inner)  then

      write(6, *) 'ERROR: exfoliation_loop only with constant_htc_inner'
      stop

     endif

     ! maximum iterations to iterate for flue gas temperature and htc_in 
     max_iter_htc_temp = 10

     ! find the distribution of the steam temperature when heat flux is known
     ! neglect the oxide inside the tube and neglect the ash deposit outside the tube
     ! id_low = 2, id_high = 1 to keep track of low load and full load
      if (TRIM(metal_temperature_model_type) == 'design')  then

        if (i_time <= 2)  then
          Temp_steam_height_time(:, i_time - 1, id_high) = Temp_steam_height_input
          Temp_steam_height_time(:, i_time - 1, id_low) = Temp_steam_height_input
        endif

      else

        Temp_steam_height_time(1, i_time, id_high) = Temp_steam_inlet
        Temp_steam_height_time(1, i_time, id_low) = Temp_steam_inlet

      endif

      ! need to reduce Tg at the low load
      Temp_gas_height(:, id_high) = &
                     & Temp_gas_height_input(1:no_height_oxide_output)
      Temp_gas_height(:, id_low) = &
                     & Temp_gas_height_input(1:no_height_oxide_output)

      j = 1
      ! at inlet steam temperature stays constant irrespective of the model
      if (TRIM(metal_temperature_model_type) == 'design')  then
        Temp_steam_height_time(j, i_time, 1:2) =  &
            & Temp_steam_height_input(1)
      else
        Temp_steam_height_time(j, i_time, 1:2) =  &
            & Temp_steam_inlet  !  Temp_steam_height_input(1)
      endif

      ! at each height determine the steam temperature
      HEIGHT_LOOP1: do j = 2, no_height_oxide_output

        ! set tube_outer_radius and tube_thickness at this location
        call GET_TUBE_DIMENSIONS(j-1)

        ! parameters from the log profile of temperature in the tube without oxide/ash
        rad_out = tube_outer_radius
        rad_in = rad_out - tube_thickness
        a1 = LOG(rad_out / rad_in)

        height_now = length_tube_oxide_output(j)
        dHeight = length_tube_oxide_output(j) - length_tube_oxide_output(j-1)

        if (TRIM(metal_temperature_model_type) == 'design')  then

          htc_tube_outer_now(id_high) =  htc_outer_height_input(j)
          htc_tube_outer_now(id_low) =  htc_outer_height_input(j)

        else

          htc_tube_outer_now(id_high) =  htc_tube_outer
          htc_tube_outer_now(id_low) =  htc_tube_outer

        endif

        ! for each load type: full (k=1) and partial (k=2)
        LOAD_LOOP1: do k = 1, 2

          ! steam temperature is UNKNOWN
          ! gas temperature and htc_inner are KNOWN
          if (TRIM(metal_temperature_model_type) == 'design')  then

            htc_tube_inner_now(k) =  htc_inner_height_input(j)

          else

            htc_tube_inner_now(k) =  htc_tube_inner

          endif

          if (i_time <= 2)  then
            ! no oxide grown
            htc_total = 1.0 / htc_tube_outer_now(k) + &
                    & rad_out / (htc_tube_inner_now(k) * rad_in) + &
                    & (rad_out / cond_value(1, 1)) * LOG(rad_out / rad_in)

          else if (i_time > 2)  then

            htc_total = 1.0 / htc_tube_outer_now(k) + &
                    & rad_out / (htc_tube_inner_now(k) * rad_in) + &
                    & (rad_out / cond_value(1, 1)) * LOG(rad_out / rad_in)

            if (if_temp_steam_depend_oxide) then
              ! steam depends on the oxide growth
              thickness_layer(first_oxide_layer:no_oxide_layers) = oxide_thick2si * &
               & thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
               & j, i_time - 1)  ! previous time
              rad_in_ox = rad_in

              do i = first_oxide_layer, no_oxide_layers

                rad_in_ox = rad_in_ox - thickness_layer(i)
                htc_total = htc_total + rad_out * &
                    & LOG((rad_in_ox + thickness_layer(i))/ rad_in_ox) / &
                    & cond_value(i, 1)

              end do

            endif

          endif

          htc_total = 1.0 / htc_total

          beta_exponent = STEAM_PIPE_EXPON(k,  &
                & Temp_steam_height_time(j-1, i_time, k), htc_total)

          ! constant gas temprature
          Tg_now = (Temp_gas_height(j-1, k) + Temp_gas_height(j, k)) / 2.0

          dTg_now = Temp_gas_height(j, k) - Temp_gas_height(j - 1, k)

          if (if_constant_steam_temp)  then

            ! constant in time; to be also constant spatially, make Temp_steam_height_input=const.
            dTSteam = 0.0

            if (TRIM(metal_temperature_model_type) == 'design')  then
              Temp_steam_height_time(j, i_time, k) =  dTSteam + &
                & Temp_steam_height_time(j, i_time - 1, k)
            else
              Temp_steam_height_time(j, i_time, k) =  dTSteam + &
                & Temp_steam_height_time(j - 1, i_time, k)
            endif

          else
            ! old; before 08/27/2012 
            ! dTSteam = (Tg_now - Temp_steam_height_time(j-1, i_time, k)) * &
            !      & (1.0 - EXP(-beta_exponent * dHeight))

            dTSteam = dTg_now + (Temp_gas_height(j, k) - &
               & Temp_steam_height_time(j-1, i_time, k) - &
               & dTg_now / (dHeight * beta_exponent)) * &
               & (1.0 - EXP(-beta_exponent * dHeight))

            Temp_steam_height_time(j, i_time, k) =  dTSteam + &
               & Temp_steam_height_time(j - 1, i_time, k)

          endif


          ! STEP H1
          ! get the heat flux at this height
          ! use interpolation between tabulated heat_flux_height(length_tube_points)
          !  if length_tube_points is different than length_tube_oxide_output
          ! heat_flux_output(j, k) = heat_flux_height(j) ! ACTION 1
          heat_flux_output(j, k) = htc_total * (Temp_gas_height(j, k) - &
                 & Temp_steam_height_time(j, i_time, k))

          if (j == 2)  then
            ! get the heat flux for j = 1
            heat_flux_output(j-1, k) = htc_total * (Temp_gas_height(j-1, k) - &
                 & Temp_steam_height_time(j-1, i_time, k))
          endif

          ! different heat fluxes are used at low/high load
          ! heat_flux_output(j,id_low)= heat_flux_output(j,id_high) * &
          !    & heat_flux_fraction_idle(1)

          if (i_time == 1)  then

            write(out_lun, 1) j, k, height_now, &
              & Temp_steam_height_time(j-1, i_time, k), &
              & Temp_steam_height_time(j, i_time, k), Temp_gas_height(j, k), &
              & dTSteam, htc_total, htc_total * (Temp_gas_height(j-1, k) - &
              & Temp_steam_height_time(j-1, i_time, k)), htc_total * &
              & (Temp_gas_height(j, k) - Temp_steam_height_time(j, i_time, k)), &
              & Temp_gas_height(j-1, k), tg_now
            write(enrg_lun, 1)  j, k, height_now, &
              & Temp_steam_height_time(j-1, i_time, k), &
              & Temp_steam_height_time(j, i_time, k), Temp_gas_height(j, k), &
              & dTSteam, htc_total, htc_total * (Temp_gas_height(j-1, k) - &
              & Temp_steam_height_time(j-1, i_time, k)), htc_total * &
              & (Temp_gas_height(j, k) - Temp_steam_height_time(j, i_time, k)), &
              & Temp_gas_height(j-1, k), tg_now

 1          format('Tsteam_1_calc ', i2, 1x, i1, 1x, &
                  &  14(1pe13.6, 1x))

          endif

        end do LOAD_LOOP1

      end do HEIGHT_LOOP1

   return

  END SUBROUTINE GET_HEAT_FLUX_ST_TEMP_OX

  SUBROUTINE TEMP_OUT_STEAM_OXIDE_METAL(i_time, id_next_outage, j_height, &
       & Temp_metal, heat_flux)

  ! oxide thickness calculated at 
  ! 1 + id_time_outage(n1-1), id_time_outage(n1)
  ! steam temperature uses dox at id_time_outage(n1-1), id_time_outage(n1) - 1
  ! dox[id_time_outage(n1-1)] is without exfoliation
  ! dox[1 + id_time_outage(n1-1)] is with exfoliation; thus write T_steam at 1,
  ! 2 + id_time_outage(n1=1, no_outages-1) and at the end id_time_outage(no_outages)

    use parameter_module,   only: mstep, moxide, oxide_thick2si, &
          & id_low, id_high
    use output_module,      only: tty_lun, out_lun, aux_lun, &
          & gen_lun, enrg_lun
    use boiler_data_module, only: Time_Boiler_p
    use boiler_data_module, only: tube_thickness, tube_outer_radius
    use oxide_data_module,   only: cond_value, &
           & thickness_fr, first_oxide_layer
    use boiler_data_module, only: id_temp_op
    use solution_data_module, only: Temp, Temp_ave, &
      & strain_rad_ave, strain_hoop_ave, displ_rad_ave, &
      & stress_axial_ave, stress_hoop_ave, stress_rad_ave, &
      & strain_th_ave, stress_hoop_th_ave, npr_temp, stress_rad, &
      & strain_th_plate_ave, strain_hoop_gen_ave
    use solution_data_module, only: no_oxide_layers, id_time_outage, no_outages
    use solver_data_module, only: if_average_oxide
    use waterwall_data_module,  only: Temp_gas_height, no_height_oxide_output, &
       & temp_steam_height_time, length_tube_oxide_output, &
       & thickness_oxide_height_time

    implicit none

    ! Argument List
    integer, intent(IN) :: id_next_outage, i_time, j_height
    real, intent(IN) :: heat_flux
    real, dimension(3), intent(IN) :: Temp_metal

    ! Local Variables
    integer  :: i, j, k, n1
    integer, dimension(30), save :: id_time
    real :: height_now
    real, dimension(30, 100), save :: Temp_metal_max_height_out, &
        & Temp_metal_mean_height_out, heat_flux_height_out
    character(LEN = 80) :: ch_repeat
    character(LEN = 1), dimension(0:9)   :: ch_index = (/'0','1','2','3','4','5', &
                  & '6','7','8','9'/)

    k = 1 ! high load and low load are the same

  if (j_height == 1 .and. id_next_outage == 1 .and. i_time == 2)  then

    do n1 = 1, no_outages + 1

      if (n1 == no_outages + 1)  then

        id_time(n1) = id_time_outage(n1 - 1)  ! (no_outages)
        id_time(n1) = id_time_outage(n1 - 1) - 1 ! (no_outages)

      else

        id_time(n1) = 2 + id_time_outage(n1 - 1)

      endif

    end do

  end if

  do n1 = 1, no_outages + 1
 
    if (id_time(n1) == i_time)  then

      Temp_metal_max_height_out(n1, j_height) = Temp_metal(3)
      Temp_metal_mean_height_out(n1, j_height) = Temp_metal(2)
      heat_flux_height_out(n1, j_height) = heat_flux

      ! write the output after data at all the heights was obtained
      if (j_height == no_height_oxide_output)  then

        if (n1 < 10)  then
          write(out_lun, 4) n1
          write(enrg_lun, 4) n1
 4        format('Tsteam_dox_height_event ', i2, ' t h dox d_sp d_mag ', &
            &  'Tst Tmw Tmax Tg hf')
        else 
          write(out_lun, 5) n1
          write(enrg_lun, 5) n1
 5        format('Tsteam_dox_height_event ', i2, ' t h dox d_sp d_mag ', &
            &  'Tst Tmw Tmax Tg hf')
        endif

        do j = 1, no_height_oxide_output

          height_now = length_tube_oxide_output(j)

          write(out_lun, 2) n1, Time_Boiler_p(i_time), height_now, &
            & SUM(thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
            & j, id_time(n1) - 1)), &
            & thickness_oxide_height_time(first_oxide_layer, j, id_time(n1) - 1), &
            & thickness_oxide_height_time(no_oxide_layers, j, id_time(n1) - 1), &
            & Temp_steam_height_time(j, id_time(n1), k), &
            & Temp_metal_mean_height_out(n1, j), &
            & Temp_metal_max_height_out(n1, j), Temp_gas_height(j, k), &
            & heat_flux_height_out(n1, j)

          write(enrg_lun, 2)  n1, Time_Boiler_p(i_time), height_now, &
            & SUM(thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
            & j, id_time(n1) - 1)), &
            & thickness_oxide_height_time(first_oxide_layer, j, id_time(n1) - 1), &
            & thickness_oxide_height_time(no_oxide_layers, j, id_time(n1) - 1), &
            & Temp_steam_height_time(j, id_time(n1), k), &
            & Temp_metal_mean_height_out(n1, j), &
            & Temp_metal_max_height_out(n1, j), Temp_gas_height(j, k), &
            & heat_flux_height_out(n1, j)

 2        format('Tsteam_dox_height_event ', i2, 20(1pe13.6, 1x))

        end do

      endif

      EXIT

    endif
  
  end do

  ! write the output after data at all the heights and all outages was obtained
  if (id_next_outage == no_outages .and. id_time(no_outages+1) == i_time .and. &
     & j_height == no_height_oxide_output)  then

    ch_repeat = ' d_ox d_sp d_mag Tst'
    write(out_lun, 3)  (ch_repeat(1:20)//ch_index(j), &
         &     j = 1, MIN(9, no_outages)), &
         & (ch_repeat(1:20)//'1'//ch_index(j-10), &
         &     j = MIN(9, no_outages) + 1, no_outages)
    write(enrg_lun, 3) (ch_repeat(1:20)//ch_index(j), &
         &     j = 1, MIN(9, no_outages)), &
         & (ch_repeat(1:20)//'1'//ch_index(j-10), &
         &     j = MIN(9, no_outages) + 1, no_outages)
 3  format('Tsteam_dox_height_final j h Tg ', 100(a))
 ! 3  format('Tsteam_dox_height_final j h Tg ', 12('d_ox Tst '))

    write(out_lun, 6) (Time_Boiler_p(id_time(n1)), n1 = 1, no_outages)
 6  format('Tsteam_time_dox_height_final ', 100(1pe13.6, 1x))

    do j = 1, no_height_oxide_output

      height_now = length_tube_oxide_output(j)

      write(out_lun, 1) j, height_now, Temp_gas_height(j, k), &
        & (SUM(thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
        & j, id_time(n1) - 1)), thickness_oxide_height_time(first_oxide_layer, &
        & j, id_time(n1) - 1), thickness_oxide_height_time(first_oxide_layer+1, &
        & j, id_time(n1) - 1), Temp_steam_height_time(j, id_time(n1), k), &
        & n1 = 1, no_outages + 1)
      write(enrg_lun, 1)  j, height_now, Temp_gas_height(j, k), &
        & (SUM(thickness_oxide_height_time(first_oxide_layer: no_oxide_layers, &
        & j, id_time(n1) - 1)), thickness_oxide_height_time(first_oxide_layer, &
        & j, id_time(n1) - 1), thickness_oxide_height_time(first_oxide_layer+1, &
        & j, id_time(n1) - 1), Temp_steam_height_time(j, id_time(n1), k), &
        & n1 = 1, no_outages + 1)

 1    format('Tsteam_dox_height_final ', i2, 1x, &
                  &  100(1pe13.6, 1x))

    end do

  endif

  return

  END SUBROUTINE TEMP_OUT_STEAM_OXIDE_METAL

END MODULE STEAM_TEMP_MODULE
