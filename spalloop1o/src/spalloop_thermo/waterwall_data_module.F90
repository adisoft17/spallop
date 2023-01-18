MODULE WATERWALL_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define all data structures integer and scalar parameters
  !   for boiler tube height variation of temperatures
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov waterwall_data_module.F90
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================
  use parameter_module, only: moxide_layer, moxide, mtime, mrate, mgrid

  implicit none
  save
  ! Public Module
  public

  ! if_new_height new height flag
  ! if_outage_ramp_exfol indicate that at the end of this time step 
  !   exfoliation is going on
  logical :: if_new_height, if_outage_ramp_exfol

  ! if_constant_steam_temp = T, then T_steam = T_steam_inlet
  logical   ::  if_constant_steam_temp 

  ! if_temp_steam_depend_oxide = T then T_steam(dox)
  ! if_temp_steam_depend_oxide = F then T_steam(dox = 0) no oxide
  logical   :: if_temp_steam_depend_oxide

  ! Tsteam(h) = TEMP_STEAM_HEIGHT
  ! T_steam = T_steam(height, 1) - for full load
  ! T_steam = T_steam(height, 2) - for partial load
  ! T_steam = T_steam(height, time, load)
  real, dimension(:, :), pointer   :: TEMP_STEAM_HEIGHT, &
       & TEMP_GAS_HEIGHT, htc_inner_height, htc_outer_height

  ! Tsteam(h) = TEMP_STEAM_HEIGHT
  ! T_steam = T_steam(height, time, load)
  real, dimension(:, :, :), pointer   :: TEMP_STEAM_HEIGHT_TIME

  integer :: height_id_stress_out = 0

  ! use htc_tube_inner - a constant value for the htc not 
  ! a variable one with temperature and pressure that considers the
  ! subcritical to supercritical transition
  logical   :: if_constant_htc_inner

  ! htc_solve_steam_or_wall = T, then htc_inner and Tsteam are unknowns
  ! htc_solve_steam_or_wall = F, then htc_inner and Twall are unknowns
  logical   :: htc_solve_steam_or_wall

  ! oxide_thickness_height(layer, height)
  real, dimension(:, :), pointer :: oxide_thickness_layer_height

  ! delta_oxide_thickness_layer_new - store delta_h during a macro cycle
  real, dimension(:, :), pointer :: delta_oxide_thickness_layer_old
  ! delta_oxide_thickness_layer_old - store delta_h during a macro cycle
  ! used for convergence
  real, dimension(:, :), pointer :: delta_oxide_thickness_layer_new
  
  ! thickness_oxide_height_time(layer, height, time)
  real, dimension(:, :, :), pointer :: thickness_oxide_height_time

  ! growth temperature for this layer at this height (layer, height)
  real, dimension(:, :), pointer :: Temp_growth_layer_height_new
  ! Temp_growth_layer_height_old is used to check the convergence
  real, dimension(:, :), pointer :: Temp_growth_layer_height_old

  ! ash thickness is dependent on height
  real, dimension(200)    :: ash_thickness

  ! ash_thickness_max = maximum ash thickness [in SI units; meters]
  real                    :: ash_thickness_max
 
  ! k_ash = thermal conductivity of ash [SI units; W/mK]
  real :: k_ash

  ! duration_ash_deslaging = time between deslaging
  ! duration_ash_deslaging(2) = time_low_load / no_ash_deslaging_per_load(2)
  ! duration_ash_deslaging(1) = time_full_load / no_ash_deslaging_per_load(1)
  real, dimension(2) :: duration_ash_deslaging
  ! number of deslaging during the low load and full load
  integer, dimension(2) :: no_ash_deslaging_per_load
 
  ! nd1_ash_max = maxiumum number of time increments within each load
  integer               :: nd1_ash_max

  ! start input for the tube U loop
  real, dimension(81) :: length_tube

  ! radius and tube thickness is constant for [length_tube(i):length_tube(i+1)]
  real, dimension(81) :: radius_outer_tube, thickness_tube

  ! orientation_tube is 'vertical', 'v', 'horizontal', or 'h' for [length_tube(i):length_tube(i+1)]
  character(LEN = 80), dimension(1:81) ::  orientation_tube

  ! length_tube_points will be set such that there will be 
  ! number_intervals_per_length(i) intervals for [length_tube(i):length_tube(i+1)]
  integer, dimension(81)               :: number_intervals_per_length
  
  ! no_length_tube = number of length points, there will be no_length_tube segments
  integer               :: no_length_tube

  ! id_length_oxide_output(k) indicate the corresponding length_tube_points
  integer, dimension(81)   ::  id_length_oxide_output

  ! id_length_tube_interval(k) indicate the length_tube interval to get dimensions
  integer, dimension(81)   ::  id_length_tube_interval

  ! end input for the tube U loop

  ! input heat flux as a function of height
  ! heat_flux_height(length_tube_points) as table input
  real, dimension(100) :: Temp_gas_height_input, Temp_steam_height_input, &
     & htc_inner_height_input, htc_outer_height_input, heat_flux_height, &
     & length_tube_points

  logical :: if_Tg_or_heat_flux

  ! heat_flux_output = heat flux at the height = length_tube_oxide_output 
  real, dimension(100, 2) :: heat_flux_output

  ! heat_flux_output_ash = heat flux at different ash thicknesses
  ! heat_flux_output_ash(height, ash, F/P)
  real, dimension(100, 40, 2) :: heat_flux_output_ash

  ! htc_tube_outer_ash = htc_tube_outer at different ash thicknesses
  ! htc_tube_outer_ash(height, ash, F/P)
  real, dimension(100, 40, 2) :: htc_tube_outer_ash

  ! height at which the oxide thickness is computed; 
  ! used if length_tube_points has too many data points, i.e., if no_heat_flux_height > 20
  real, dimension(100) :: length_tube_oxide_output

  ! number of data points in the length_tube_oxide_output
  integer :: no_height_oxide_output

  ! if_height_oxide_output = T, then use length_tube_oxide_output
  logical :: if_height_oxide_output

  ! number of data points in the heat flux data as a function of height
  integer :: no_heat_flux_height

  ! number of height intervals in the boiler
  integer :: no_delta_height

  ! height_boiler - the height of the boiler
  ! delta_height = height_boiler / no_delta_height the height increment
  real    :: height_boiler, delta_height

  ! duration_ash_deslaging will be divided in no_time_intervals_ash
  integer :: no_time_intervals_ash

  ! ratio between adjacent ash intervals: 
  ! dt_ash(i) = dt_ash(i-1) * ratio_time_intervals_ash
  ! time_relative_ash(i) = time_relative_ash(i-1) + dt_ash(i)
  real    :: ratio_time_intervals_ash
  real, dimension(2, 10)    :: dt_ash, time_relative_ash

  ! if_ash = F do not consider the ash effects
  logical  :: if_ash

  ! temperature of the steam at the inlet
  real     :: Temp_steam_inlet

  ! flow rate at full load
  real     :: flow_rate_steam_high, flow_rate_steam_low

  ! time_now = time at which the h_ox, htc_inner, T_steam, etc are computed
  ! by accounting for ash/deslaging; this includes subintervals to boiler load regimes
  real, dimension(:, :), pointer  :: time_now

  ! how many time intervals are within one segment (full, full 2 low, low, low 2 full)
  integer, dimension(:), pointer  :: no_thick_interval_now

!!!5wz - geometry factor for the heat flux
  ! whether or not to use the geometry factor
  logical :: if_use_geom_factor = .TRUE.
  ! spacing between two neighbor tubes (m)
  real :: tube_spacing_m  = 0.045
  ! thickness of rib (m)
  real :: rib_thickness_m = 0.008
  ! geometry factor (unitless)
  real :: geom_factor_heat_flux = 1.0

  ! if_debug_oxide_growth_conv = T for printing debug info
  logical  :: if_debug_oxide_growth_conv
  logical  :: if_debug_steam_htc_inner

  ! if_output_height_time_all = T then height_output_ash will be written at output
  logical  :: if_output_height_time_all

  ! hflux_htc_ash_estimate = 1, then heat_flux(ash) will be estimated
  !   htc_gas(ash) = htc_gas(no ash)
  !_hflux_htc_ash_estimate = 2, htc_gas(ash) will be estimated by using Delta_Tmean_metal 
  integer  :: hflux_htc_ash_estimate

  ! Delta_Temp_metal_mean_ash_max = maximum temperature drop in the mean wall temperature due to ash= Temp_metal_mean_ash_max(no ash) - Temp_metal_mean_ash_max(ash)
  real     :: Delta_Temp_metal_mean_ash_max

  ! if loop_transition_length <= 0 or max(length_tube) < 10 then only one loop
  real     :: loop_transition_length

  integer  :: no_loops

  character(LEN = 80), dimension(1:5) :: temp_design_data

  ! metal_temperature_model_type = 'design' or 'model_steam_inlet_gas_temp'
  character(LEN = 80)  :: metal_temperature_model_type

  ! if metal_temperature_model_type = 'design' metal temperature must be read
  ! from temp_design_file
  character(LEN = 80)  :: temp_design_file

  integer  :: id_height_now
  
  ! indices to be used for storring creep related quantities
  integer  :: jheight_creep, ioutage_creep, ioutage_previous_creep

END MODULE WATERWALL_DATA_MODULE
