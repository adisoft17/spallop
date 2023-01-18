 MODULE WATERWALL_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Input all the data for the waterwall - height dependence, ash, etc
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov waterwall_input_module.F90
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: WATERWALL_INPUT, WATERWALL_INIT

 CONTAINS

 SUBROUTINE WATERWALL_INPUT(fatal)

  use parameter_module,     only : mstep, zero
  use input_utilities_module,   only: SEARCH_NML
  use output_module,            only: tty_lun, out_lun, inp_lun
  use waterwall_data_module,    only: ash_thickness_max, &
       & duration_ash_deslaging, no_ash_deslaging_per_load, &
       & no_ash_deslaging_per_load, no_time_intervals_ash, &
       & ratio_time_intervals_ash, heat_flux_height, length_tube_points, &
       & no_heat_flux_height, if_ash, no_heat_flux_height, no_delta_height, &
       & height_boiler, length_tube_oxide_output, k_ash, &
       & Temp_steam_inlet, flow_rate_steam_high, flow_rate_steam_low, &
       & no_heat_flux_height, if_height_oxide_output, &
       & no_height_oxide_output, if_constant_htc_inner, &
       & htc_solve_steam_or_wall, if_debug_oxide_growth_conv, &
       & if_debug_steam_htc_inner, if_output_height_time_all, &
       & hflux_htc_ash_estimate, Delta_Temp_metal_mean_ash_max
 use waterwall_data_module,    only: Temp_gas_height_input, if_Tg_or_heat_flux
 use waterwall_data_module,    only: length_tube, radius_outer_tube, &
       & thickness_tube, orientation_tube, number_intervals_per_length, &
       & no_length_tube, id_length_oxide_output, id_length_tube_interval, &
       & height_id_stress_out, if_constant_steam_temp, &
       & if_temp_steam_depend_oxide, loop_transition_length, &
       & metal_temperature_model_type, temp_design_file, no_loops, &
       & Temp_steam_height_input,  htc_inner_height_input, htc_outer_height_input

     implicit none

   ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i, k, j
    logical :: no_solv_namelist, solv_namelist
    real    :: factor1, factor2, min1, max1, delta_length1, small_length

  ! default

    namelist /waterwall/  heat_flux_height, Temp_gas_height_input, &
       & Temp_steam_height_input, htc_inner_height_input, htc_outer_height_input, &
       & length_tube, number_intervals_per_length, radius_outer_tube, &
       & thickness_tube, orientation_tube, length_tube_oxide_output, &
       & if_Tg_or_heat_flux, no_delta_height, height_boiler, &
       & Temp_steam_inlet, flow_rate_steam_high, flow_rate_steam_low, &
       & if_height_oxide_output, if_constant_htc_inner, &
       & if_debug_oxide_growth_conv, &
       & if_debug_steam_htc_inner, if_output_height_time_all, &
       & height_id_stress_out, if_constant_steam_temp, &
       & if_temp_steam_depend_oxide, loop_transition_length, &
       & metal_temperature_model_type, temp_design_file

  if_constant_steam_temp = .false.

  duration_ash_deslaging = 0.0
  no_time_intervals_ash = 0
  ash_thickness_max = 0.0

  if_ash = .false.
  if_constant_htc_inner = .true.
  htc_solve_steam_or_wall = .true.

  no_heat_flux_height = 0
  ratio_time_intervals_ash = 1.0 
 
  length_tube_oxide_output = 0.0
  if_height_oxide_output = .false.

  ! no debug info
  if_debug_oxide_growth_conv = .false.
  if_debug_steam_htc_inner = .false.

  if_output_height_time_all = .false.
  hflux_htc_ash_estimate = 0  ! estimate Tst, same htc_outer, q, Tg
  hflux_htc_ash_estimate = 1  ! estimate heat_flux(ash), same q, Tg, Tst, htc_outer
  hflux_htc_ash_estimate = 2  ! estimate heat_flux(ash) & htc_outer, same q, Tg, Tst
  hflux_htc_ash_estimate = 0  ! estimate Tst, same htc_outer, q, Tg

  heat_flux_height = 0.0
  Temp_gas_height_input = 0.0
  Temp_steam_height_input = 0.0
  if_Tg_or_heat_flux = .true. ! use Tg at the input

  ! height_id_stress_out is the interval id for the height at which the 
  ! stress output will be written
  height_id_stress_out = 1

  orientation_tube = 'none'

  Delta_Temp_metal_mean_ash_max = 100.0  ! default

  if_temp_steam_depend_oxide = .false.

  loop_transition_length = 0.0
  no_loops = 1

  metal_temperature_model_type = 'model_steam_inlet_gas_temp'

  temp_design_file = 'none'

       ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'waterwall', 'WATERWALL')
       solv_namelist = .NOT. no_solv_namelist

       if (solv_namelist) then
          read (inp_lun, NML= waterwall, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (.not. fatal)  then
            write (tty_lun, 15)
            write (out_lun, 15)
15          format (/' DONE Reading WATERWALL Namelist ...')
          else
            write(tty_lun, *) 'STOP: WATERWALL Input'
            stop
          endif
       end if

      no_length_tube = 0
      do i = 40, 1, -1
        if (no_length_tube == 0 .and. length_tube(i) /= zero) &
           & no_length_tube = i
      end do

      if (no_length_tube <=0)  then
        write(tty_lun, *) 'ERROR: no_length_tube must be > 0'
        stop
      else
        no_heat_flux_height = 1
        write(out_lun, *) 'no_length_tube = ', no_length_tube
      endif

      do i = 1, no_length_tube - 1

        write(out_lun, *) 'orientation_tube ', i, TRIM(orientation_tube(i))

        if (TRIM(orientation_tube(i))=='vertical' .or. &
          & TRIM(orientation_tube(i))=='VERTICAL' .or. &
          & TRIM(orientation_tube(i))=='v' .or. &
          & TRIM(orientation_tube(i))=='V' .or. &
          & TRIM(orientation_tube(i))=='horizontal' .or. &
          & TRIM(orientation_tube(i))=='HORIZONTAL' .or. &
          & TRIM(orientation_tube(i))=='h' .or. &
          & TRIM(orientation_tube(i))=='H') then

          no_heat_flux_height = no_heat_flux_height + number_intervals_per_length(i)

          CYCLE 
        else
          write(tty_lun, *) 'ERROR:orientation_tube must be vertical, v, horizontal, or h'
          write(tty_lun, *) 'ERROR:orientation_tube must be VERTICAL, V, HORIZONTAL, or H'
          STOP
        endif

      end do

      if (no_heat_flux_height <= 1)  then
        write(tty_lun, *) 'ERROR: no_heat_flux_height must be > 1'
        stop
      else
        write(out_lun, *) 'no_heat_flux_height = ', no_heat_flux_height
      endif

      length_tube_points(1) = 0.0
      j = 1
      write(out_lun, *) 'length_tube_points set ', j, length_tube_points(j) 

      do i = 1, no_length_tube - 1

        delta_length1 = (length_tube(i+1) - length_tube(i)) / &
           & number_intervals_per_length(i)
        do k = 1, number_intervals_per_length(i)
          j = j + 1
          length_tube_points(j) = length_tube(i) + delta_length1 * k
          write(out_lun, *) 'length_tube_points set ', j, length_tube_points(j) 
        end do
      end do

     if (no_ash_deslaging_per_load(1) == 0 .or. &
       & no_time_intervals_ash == 0)  then
       if_ash = .false.
       write(out_lun, *) 'if_ash = F since no_time_intervals_ash = 0', &
         & ' or no_ash_deslaging_per_load(1) = 0 '

     else

       if_ash = .true.

       if (ABS(ash_thickness_max) < 1.0e-4)  then
         if_ash = .false.
         write(out_lun, *) 'if_ash = F since ash_thickness_max < 1.0e-4 '
       else 
         write(out_lun, *) 'if_ash = T since ash_thickness_max and ', &
           & ' no_time_intervals_ash are given'
       endif

     endif

      no_height_oxide_output = 0

       do i = 100, 1, -1
         if (no_height_oxide_output == 0 .and. &
           & length_tube_oxide_output(i) /= zero) &
           & no_height_oxide_output = i
       end do

       if (no_height_oxide_output <=0)  then
         write(tty_lun, *) 'WARNING: no_height_oxide_output must be > 0'
         if_height_oxide_output = .false.

       else
         if_height_oxide_output = .true.
         ! in this case length_tube_oxide_output could be different than length_tube_points
         write(out_lun, *) 'no_height_oxide_output = ', no_height_oxide_output

       endif

     if (no_height_oxide_output == 0)  then

       if (no_heat_flux_height > 10)  then

         write(tty_lun, *) 'ERROR: no_height_oxide_output would be too large'
         write(tty_lun, *) 'use input no_heat_flux_height - less than 30 values'
         STOP

       else
         
         length_tube_oxide_output(1: no_heat_flux_height) = &
            & length_tube_points(1: no_heat_flux_height)
         no_height_oxide_output = no_heat_flux_height
         if_height_oxide_output = .false.

       endif

     else if (no_height_oxide_output > 0)  then

       ! identify the closest length_tube_points value
       do k = 1, no_height_oxide_output

         id_length_oxide_output(k) = 0
         small_length = length_tube_points(no_heat_flux_height)
         do i = 1, no_heat_flux_height
           if (ABS(length_tube_oxide_output(k) - length_tube_points(i)) < &
              & small_length) then

             ! change the ID 
             id_length_oxide_output(k) = i
             small_length = ABS(length_tube_oxide_output(k) - length_tube_points(i))
           endif
         end do

         if (id_length_oxide_output(k) == 0)  then
           write(tty_lun, *) 'ERROR:id_length_oxide_output = 0 for ', k
           STOP
         endif

         length_tube_oxide_output(k) = length_tube_points(id_length_oxide_output(k))

         ! identify the kind of tube, i.e., the tube dimensions might change
         do i = 1, no_length_tube - 1

           if (length_tube_oxide_output(k) >= length_tube(i) .and. &
             & length_tube_oxide_output(k) < length_tube(i+1)) then

             ! change the ID for tube_points [ltp(i):ltp(i+1)]
             id_length_tube_interval(k) = i

           endif

         end do

         i = no_length_tube
         if (ABS(length_tube_oxide_output(k) - length_tube(i)) <= 0.01)  then
           id_length_tube_interval(k) = i
         endif

       end do

       write(out_lun, 2) (length_tube_oxide_output(k), &
          & k = 1, no_height_oxide_output)
  2    format('length_tube_oxide_output = ', 50(1pe13.6, 1x))
       write(out_lun, 3) (id_length_oxide_output(k), &
          & k = 1, no_height_oxide_output)
  3    format('id_length_oxide_output = ', 50(i2, 1x))

       write(out_lun, 5) (id_length_tube_interval(k), &
          & k = 1, no_height_oxide_output)
  5    format('id_length_tube_interval = ', 50(i2, 1x))


     endif

    if (if_constant_htc_inner)  then

    endif

    if (if_constant_steam_temp) then
      write(out_lun, *) 'constant_steam_temperature'
    else
      write(out_lun, *) 'constant_steam_temperature'
    endif

    if (loop_transition_length > length_tube(no_length_tube) / 20.0)  then
      no_loops = 2
    endif

    if  (.not. (TRIM(metal_temperature_model_type) == 'model_steam_inlet_gas_temp' .or. &
      & TRIM(metal_temperature_model_type) == 'design'))  then
      write(6, *) 'ERROR: metal_temperature_model_type must be ', &
        & 'model_steam_inlet_gas_temp or design'
      STOP

    endif


  return

  END SUBROUTINE WATERWALL_INPUT
  
  SUBROUTINE WATERWALL_INIT()
  !=======================================================================
  ! Purpose(s):
  !
  !   Initialize variables for the waterwall
  !
  !=======================================================================

  use parameter_module,     only : mstep, zero, id_low, id_high
  use input_utilities_module,   only: SEARCH_NML
  use output_module,            only: tty_lun, out_lun, inp_lun
  use waterwall_data_module,    only: duration_ash_deslaging, &
       & no_ash_deslaging_per_load, no_time_intervals_ash, &
       & ratio_time_intervals_ash, heat_flux_height, length_tube_points, &
       & no_heat_flux_height, if_ash, no_heat_flux_height, no_delta_height, &
       & height_boiler, Temp_steam_inlet, &
       & flow_rate_steam_high, flow_rate_steam_low, TEMP_STEAM_HEIGHT, &
       & TEMP_GAS_HEIGHT, oxide_thickness_layer_height, ash_thickness, &
       & dt_ash, time_relative_ash, ash_thickness_max, &
       & htc_inner_height, no_height_oxide_output, &
       & no_thick_interval_now, time_now
  use waterwall_data_module,    only: delta_oxide_thickness_layer_old, &
       & delta_oxide_thickness_layer_new, Temp_growth_layer_height_new, &
       & Temp_growth_layer_height_old, nd1_ash_max, if_constant_steam_temp
  use boiler_data_module, only: TEMP_BOILER_data_no, time_boiler_idle, &
       & time_boiler_pulse, Time_Boiler_p, id_temp_op
  use oxide_data_module,        only: no_oxide, no_layer, &
        & no_metal_layers, first_oxide_layer, no_layers_oxide
  use solution_data_module,      only: no_oxide_layers
  use solver_data_module,       only: no_thick_interval
  use htc_utility_module, only: GET_ASH_THICKNESS

  implicit none

  ! Local Variables
  integer :: ioerror, i, k, j, Status, nd1_max
  integer :: i1, i1c, i2, k_load_start, k_load_end
  integer, dimension(3)  :: nd1
  real, dimension(2) :: time_load
  real    :: ash_thickness_now

  write(out_lun,*) 'allocate waterwall data ', no_height_oxide_output
  ALLOCATE(TEMP_STEAM_HEIGHT(1:no_height_oxide_output, 1:2), STAT = Status)
  ALLOCATE(TEMP_GAS_HEIGHT(1:no_height_oxide_output, 1:2), STAT = Status)

  ! do not store oxide thickness in time; too much data
  ALLOCATE(oxide_thickness_layer_height(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)
  oxide_thickness_layer_height = 0.0

  ALLOCATE(delta_oxide_thickness_layer_old(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)
  ALLOCATE(delta_oxide_thickness_layer_new(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)

  ! allocate the growth temperature
  ALLOCATE(Temp_growth_layer_height_new(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)
  ALLOCATE(Temp_growth_layer_height_old(first_oxide_layer:no_oxide_layers, &
       & 1:no_height_oxide_output), STAT = Status)

  ALLOCATE(htc_inner_height(1: no_height_oxide_output, 1:2), STAT = Status)

  ! allocate data for the strain-stress analysis, the creep data needs to be stored

  ! allocate data for the blockage

  ! ash_thickness does not depend on height since we do not have an equation
  ! for its growth based on metal temperature

  write(out_lun,*) 'setup ash data '
  ! duration_ash_deslaging = time between deslaging
  
  ! if (if_ash)  then

  ASH_LOOP1: do k = 1, 2
  
    if (k == 1)  then   ! full load
      time_load(k) = time_boiler_pulse(1)
    else if (k == 2)  then   ! low load
      time_load(k) = time_boiler_idle(1)
    endif

    if (no_time_intervals_ash == 0)  then 

      if (k == 1)  then   ! full load
        if (no_ash_deslaging_per_load(k) <= 0)  & 
          & no_ash_deslaging_per_load(k) = no_thick_interval
      else if (k == 2)  then ! low load
        if (no_ash_deslaging_per_load(k) <= 0)  & 
          & no_ash_deslaging_per_load(k) = 2
      endif

    endif

    if (no_ash_deslaging_per_load(k) > 0)  then
      
       duration_ash_deslaging(k) = time_load(k) / no_ash_deslaging_per_load(k)
  
    else
  
       write(2, *)  ' use duration_ash_deslaging >= 4 but <= 8 '
       write(2, *)  ' better use no_ash_deslaging_per_load instead of duration_ash_deslaging such that ', & 
          & 'duration_ash_deslaging = time_load / no_ash_deslaging_per_load'
       
       if (k == 1)  then
         write(2, *) 'choose no_ash_deslaging_per_load_full = 36 or 18 for WEEKLY schedule such that ', &
          &  'duration_ash_deslaging_full = time_load / 36 or 18'
         write(2, *) 'choose no_ash_deslaging_per_load_full = 2 or 3 for DAILY schedule such that ', &
          &  'duration_ash_deslaging_full = time_load / 2 or 3'
       else if (k == 2)  then
          write(2, *) 'choose no_ash_deslaging_per_load_low = 3 or 5 for WEEKLY schedule such that ', &
          &  'duration_ash_deslaging_low = time_load / 3 or 5'
         write(2, *) 'choose no_ash_deslaging_per_load_low = 1 or 2 for daily schedule such that ', &
          &  'duration_ash_deslaging_low = time_load / 1 or 2'   
       endif
    
       if (no_time_intervals_ash > 0 .and. &
           & (duration_ash_deslaging(k) < 2.0 .or. &
           & time_load(k) / duration_ash_deslaging(k) < 1))  then
    
         write(2, *) 'ERROR: too large time between deslaging: use the following'
         write(6, *) 'ERROR: too large time between deslaging: use the following'
      
         if (k == 1)  then
           write(6, *) 'choose no_ash_deslaging_per_load_full = 36 or 18 for WEEKLY schedule such that ', &
          &  'duration_ash_deslaging_full_full = time_load / 36 or 18'
           write(6, *) 'choose no_ash_deslaging_per_load = 2 or 3 for DAILY schedule such that ', &
          &  'duration_ash_deslaging_full = time_load / 2 or 3'
         else if (k == 2)  then
           write(6, *) 'choose no_ash_deslaging_per_load_low = 3 or 5 for WEEKLY schedule such that ', &
          &  'duration_ash_deslaging_low = time_load / 3 or 5'
           write(6, *) 'choose no_ash_deslaging_per_load_low = 1 or 2 for DAILY schedule such that ', &
          &  'duration_ash_deslaging_low = time_load / 1 or 2'   
         endif
      
         stop
      
       else if (no_time_intervals_ash > 0 .and. &
           & duration_ash_deslaging(k) >= 2.0 .and. &
           & time_load(k) / duration_ash_deslaging(k) >= 1)  then
    
         no_ash_deslaging_per_load(k) = INT(time_load(k) / &
             & duration_ash_deslaging(k))
      
       endif
    
     endif 
  
     write(2, *) 'ash_duration_number_deslags ', k, duration_ash_deslaging(k),  &
          & no_ash_deslaging_per_load(k)
  
     time_relative_ash(k, 1) = 0.0
  
     if (ABS(ratio_time_intervals_ash - 1.0)  > 1.0e-5)  then
  
       dt_ash(k, 1) = duration_ash_deslaging(k) * &
         &  (ratio_time_intervals_ash - 1.0) / &
         &  (ratio_time_intervals_ash**(no_time_intervals_ash+1) - 1.0)
        
     else
  
       dt_ash(k, 1) = duration_ash_deslaging(k) / (no_time_intervals_ash+1)
    
     endif

     write(2, *) 'dt_ash_first ', k, dt_ash(k, 1), duration_ash_deslaging(k), &
          & no_time_intervals_ash
  
     do i = 2, no_time_intervals_ash+1
       ! no_time_intervals_ash intervals but points are no_time_intervals_ash + 1
       dt_ash(k, i) =  dt_ash(k, i-1) * ratio_time_intervals_ash
       time_relative_ash(k, i) = time_relative_ash(k, i-1) + dt_ash(k, i) 
    
     end do
  
  end do ASH_LOOP1

  ! else if (.not. if_ash)  then

  ! endif
    
  ! set up the ash thickness

  ! if (no_time_intervals_ash > 0)  then
    
    ! linear profile
    do i = 1, no_time_intervals_ash
      ash_thickness(i) = ash_thickness_max * (i-1) / &
         & (no_time_intervals_ash-1)
    end do
    ash_thickness(no_time_intervals_ash+1) = 0.0 !  1 and N coincide

    do k = 1, 2
      do i = 1, no_time_intervals_ash+1
        if (k == id_high)  then
          write(out_lun, 2) time_relative_ash(k, i), ash_thickness(i)
 2        format('ash_thikness(t)_full_load', 3(1pe13.6, 1x))
        else if (k == id_low)  then
          write(out_lun, 4) time_relative_ash(k, i), ash_thickness(i)
 4        format('ash_thikness(t)_partial_load', 3(1pe13.6, 1x))
        endif
      end do
    end do

  ! endif

  ALLOCATE(no_thick_interval_now(1:TEMP_BOILER_data_no), STAT = Status)

  ! get the maximum number of subintervals; use nd1 as dummy
  nd1(1) = no_thick_interval
  ! works even when no_time_intervals_ash = 1 since 
  !    no_ash_deslaging_per_load was set non zero
  nd1(2) = (no_time_intervals_ash + 1) * &
            & no_ash_deslaging_per_load(id_low)
  nd1(3) = (no_time_intervals_ash + 1) * &
            & no_ash_deslaging_per_load(id_high)

  nd1_ash_max = MAXVAL(nd1)
  ALLOCATE(time_now(1:TEMP_BOILER_data_no, 0:nd1_ash_max), STAT = Status)

  ! initialize for the sake of completion; needed for i = 2
  ! no_thick_interval_now(1) = 1 not needed
  ! time_now(1, 1) = 0.0

  ! setting up the time intervals that considers ash growth and deslaging
  ! this is time interval [i-1, i]
      TIME_LOOP_HOX1: do i= 2, TEMP_BOILER_data_no ! no of data intervals

        ! get the time at the start of the i-th interval; these should be equal
        ! time_now(i, 0) = time_now(i-1, no_thick_interval_now(i-1))
        time_now(i, 0) = Time_Boiler_p(i-1)

        ! find the appropriate numbers of time intervals
        if (id_temp_op(i) == 1)  then 
          ! interval 4-1, isothermal low load
          ! if_ash = F then no_thick_interval_now = 1

          if (i > 2)  then
            no_thick_interval_now(i) = (no_time_intervals_ash +1) * &
              & no_ash_deslaging_per_load(id_low)
            
            if (no_thick_interval_now(i) > 0)  then

              do i1 = 1, no_ash_deslaging_per_load(id_low)
                i1c = (i1-1) * (no_time_intervals_ash + 1)
                do i2 = 1, no_time_intervals_ash + 1
                  time_now(i, i1c+i2) = time_now(i, i1c+i2-1) + dt_ash(id_low, i2)
                end do
              end do
            
            else 
              no_thick_interval_now(i) = 1
              time_now(i, 1) = Time_Boiler_p(i)
            endif

          else if (i == 2)  then
            ! first setup is different
            no_thick_interval_now(i) = 1
            time_now(i, 1) = Time_Boiler_p(i)
          endif

        else if (id_temp_op(i) == 2)  then 
          ! interval 1-2, transition partial - full load
          ! no deslaging but temperature variation
          no_thick_interval_now(i) = no_thick_interval
          do i1 = 1, no_thick_interval_now(i)
            time_now(i, i1) = time_now(i, i1-1) + &
                & (Time_Boiler_p(i) - Time_Boiler_p(i-1)) / &
                & no_thick_interval_now(i)
          end do
          
        else if (id_temp_op(i) == 3)  then 
          ! interval 2-3, isothermal high load
          ! if_ash = F then no_thick_interval_now = 1
          no_thick_interval_now(i) = (no_time_intervals_ash + 1) * &
            & no_ash_deslaging_per_load(id_high)
            
          if (no_thick_interval_now(i) > 0)  then
          
            do i1 = 1, no_ash_deslaging_per_load(id_high)
              i1c = (i1-1) * (no_time_intervals_ash + 1)
              do i2 = 1, no_time_intervals_ash + 1
                time_now(i, i1c+i2) = time_now(i, i1c+i2-1) + dt_ash(id_high, i2)
              end do
            end do
            
          else 
            no_thick_interval_now(i) = 1
            time_now(i, 1) = Time_Boiler_p(i)
          endif

        else if (id_temp_op(i) == 4)  then 
          ! interval 3-4, transition full - partial load
          ! no deslaging but temperature variation
          no_thick_interval_now(i) = no_thick_interval
          do i1 = 1, no_thick_interval_now(i)
            time_now(i, i1) = time_now(i, i1-1) + &
                & (Time_Boiler_p(i) - Time_Boiler_p(i-1)) / &
                & no_thick_interval_now(i)
          end do

        endif

        if (i <= 10)  then

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
          do i1 = 1, no_thick_interval_now(i)
            ash_thickness_now = GET_ASH_THICKNESS(i1, k_load_start, k_load_end)
            write(out_lun, 3) i, i1, id_temp_op(i), k_load_start, k_load_end, &
              & time_now(i, i1) - time_now(i, i1-1), time_now(i, i1), &
              & Time_Boiler_p(i-1), Time_Boiler_p(i), ash_thickness_now
 3          format('time_now_out ', i3, 1x, i2, 1x, 3(i1, 1x), 30(1pe13.6, 1x))
          end do
        endif

    END DO TIME_LOOP_HOX1

  RETURN

 END SUBROUTINE WATERWALL_INIT
 

END MODULE WATERWALL_INPUT_MODULE
