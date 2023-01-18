 MODULE BLOCKAGE_DATA_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Input all the data for the blockage - height dependence
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================

  use parameter_module, only: moxide_layer, moxide, mramp_out

  implicit none

  save
  ! Public Module
  public

  ! if_blockage_area = T then all the data need to be input for the volume exfoliation and blockage area
  logical            :: if_blockage_area

  ! if_exfoliation_volume = T, then all the data need to be input for the volume exfoliation 
  logical             :: if_exfoliation_volume

  ! if_exfoliation_oxide_growth = T, then oxide growth is affected by exfoliation
  logical             :: if_exfoliation_oxide_growth

  integer, parameter  :: mblock = 30

  integer, parameter  :: mskin = 100  ! maximum numbers of skins to exfoliate

  ! area_fraction_spalled(energy_elastic_oxide) as given by Armitt report
  real, dimension(mblock) :: area_fraction_spalled, energy_elastic_oxide

  ! numbers of data segments in the area_fraction_spalled(energy_elastic_oxide)
  integer             :: no_area_fr_spalled

   ! energy_armitt_factor_scale - 
   ! energy_elastic_oxide will be multiplied with energy_armitt_factor_scale
   ! in order to use only the energy in the magnetite for magnetite exfoliation
   ! or use energy in the magnetite and haematite layers (all but spinel) 
   !                      since they exfoliate only
   ! energy_armitt_factor_scale = 0.5 since energy_elastic_oxide(total_oxide_thickness)
   ! and magnetite thickness = 0.5 total_oxide_thickness
   ! energy_armitt_factor_scale = x when
   ! x = (magnetite+haematite)_thickness / total_oxide_thickness
   real      :: energy_armitt_factor_scale

   ! translate the armitt's energy levels to lower ones
   ! thus, the shape of the fraction_area(energy) is the same
   ! this is needed for 347HFG as the oxide layer is too thin
   real      :: energy_at_zero_area_fraction_spalled

  ! accumulation_deposit_area - an equivalent area over which the deposits will be accumulated
  ! deposit_porosity_fraction = porosity in the deposits since the oxide
  !     flakes cannot be compressed fully dense
  ! non_blocking_vol_fraction = volumetric fraction of the exofoliated deposits
  !    that are dispersed horizontally and do not account for blockage area
  real, dimension(mblock) :: accumulation_deposit_area,  &
      & non_blocking_vol_fraction, &
      & radius_deposit_cylinder, volume_deposit_2cyl, volume_deposit_cube, &
      & volume_deposit_bend, volume_deposit_1cyl
  real, dimension(mblock) :: deposit_fluid_density, deposit_oxide_density, &
      & deposit_porosity_fraction
 
  ! start data for deposits not in the bend but along the cylinder
  ! 
  real       :: radius_deposit_tube
  real, dimension(10) :: length_deposit_tube

  ! no_length_deposit_tube - number of length_deposit_cyl data
  integer :: no_length_deposit_tube

  ! number_deposits_along_tube - number of simultaneous blockage sites,
  ! i.e., the exfoliated volume will be uniformly distributed among several sits
  !  this is used as a parameter by Armitt
  ! number_deposits_along_tube = 1, 2, 3, 5, 10, 20
  integer, dimension(7)  :: number_deposits_along_tube

  ! no_deposits_along_tube = number of data in number_deposits_along_tube
  integer     :: no_deposits_along_tube

  ! volume_deposit_tube(outage, layer, length_deposit, # of identical deposits)
  real, dimension(mramp_out, 10, 10, 7) :: volume_deposit_tube, &
      & area_deposit_tube, &
      & area_fraction_blocked_tube, area_fraction_blocked_tube0

  ! area_exfoliated_fr_block(height, outage, layer) for output
  real, dimension(100, mramp_out, moxide_layer)  :: area_exfoliated_fr_block

  ! END data for deposits not in the bend but along the cylinder

  ! area_fraction_blocked = fractional area covered with exfoliated oxide deposits
  real, dimension(10, mramp_out, moxide_layer) :: area_fraction_blocked_bend,  &
      & volume_fraction_blocked_bend, height_deposit_bend, &
      & area_fraction_blocked_cube, volume_fraction_blocked_cube, &
      & height_deposit_cube, area_fraction_blocked_2cyl,  &
      & volume_fraction_blocked_2cyl, height_deposit_2cyl, &
      & area_fraction_blocked_1cyl,  &
      & volume_fraction_blocked_1cyl, height_deposit_1cyl
  !  real, dimension(mblock, moxide_layer) :: 

  ! no_deposit_sites = number of accumulation_deposit_area or sites that can be blocked by exofoliated oxides
  ! no_deposit_sites = 2 for superheater - 2 corners
  integer    :: no_deposit_sites  

  ! deposit_origin_id - identifies the id of the exfoliated mass
  !  deposit_origin_id points to the [1:no_height_oxide_output] or 
  !  deposit_origin_id points(bend, height) = t or falso
  ! integer, dimension(mblock, 100) :: deposit_origin_id
  ! id_height_exfoliate2deposit(height) = deposit_id - or bend id
  integer, dimension(100) :: id_height_exfoliate2deposit

  integer    :: no_deposit_origin ! numbers of height segments where exfoliation occurs

  ! area_exfoliated_origin = area that is exfoliated at the origin, i.e., at a certain location inside the tube
  ! thickness_exfoliated_origin = thickness of the exfoliated mass at each location
  ! volume_exfoliated_origin = volume that is exfoliated at each location
  real, dimension(mblock, moxide_layer) :: area_exfoliated_origin
  real, dimension(mblock, moxide_layer) :: thickness_exfoliated_origin
  real, dimension(mblock, moxide_layer) :: volume_exfoliated_origin
  ! volume_exfoliated_bend is fully solid, i.e., 0 porosity
  ! volume_exfoliated_bend(bend, outage, layer)
  ! volume_exfoliated_bend_porous(bend, outage, layer) - accoounts for porosity in the deposit
  real, dimension(10, mramp_out, moxide_layer) :: volume_exfoliated_bend, &
    & volume_exfoliated_bend_porous

   ! (1) will store the total elastic energy, then 2,..N will store for each oxide layer
  ! 100 consistent with length_tube_oxide_output
  ! max_energy_elastic_event_block(height, outage, layer) = max_energy_elastic_event_out(outage, layer)
   real, dimension(100, mramp_out, moxide_layer)   :: max_energy_elastic_event_block

   ! max_strain_elast_event_block(height, outage, layer) = max_strain_elast_event_out(outage, layer) 
   real, dimension(100, mramp_out, moxide_layer)   :: min_strain_elast_event_block, &
       & max_strain_elast_event_block

   ! rad_int_block(height, outage, layer)
   real, dimension(100, mramp_out, moxide_layer)   :: rad_int_block, oxide_thickness_block

   ! skin_no(height, layer) - stores the current skin number that's being exfoliated
   ! fraction_exfoliated_skin(i < skin_no) = 1.0
   ! fraction_exfoliated_skin(skin_no) <= 1.0
   integer, dimension(100, moxide_layer)    :: skin_no   

   !  skin_no_start(height), skin_no_end(height)
   integer, dimension(100)    :: skin_no_start, skin_no_end

   ! fraction_exfoliated_skin(height, skin_id, layer) - stores the fraction of normal area, not in cross-section, 
   ! fraction_exfoliated_skin is an accumulation of the area fractions given by Armitt
   real, dimension(100, mskin, moxide_layer)   :: fraction_exfoliated_skin

   ! dthick_ox(height, outage, layer) - oxide grown between outage events
   real, dimension(100, 0:mramp_out, moxide_layer)   :: dthick_ox

   ! dthick_ox_height_time(outage, height, time) for magnetite
   real, dimension(:, :, :), pointer :: dthick_ox_height_time

   ! growth_mode_after_exfoliation - indicate how does the oxide grows after exfoliation
   ! growth_mode_after_exfoliation - continuous, cont 
   ! growth_mode_after_exfoliation - zero
   ! growth_mode_after_exfoliation - start_first_exfolation
   ! growth_mode_after_exfoliation - from_given_thickness
   character(LEN = 80)   ::  growth_mode_after_exfoliation

   ! growth_mode_after_exfoliation = 'from_given_thickness' then thickness_new_ox_growth in microns will be the initial oxide thickness
   real   :: thickness_new_ox_growth

   ! thickness_new_ox_growth_init(height, outage, layer) - to initialize dz2 
   real, dimension(100, 0:mramp_out, moxide_layer)   :: thickness_new_ox_growth_init

   ! thickness_ave_ox_growth_init(height, outage, layer) - to initialize average magnetite thickness 
   ! 
   real, dimension(100, 0:mramp_out, moxide_layer)   :: thickness_ave_ox_growth_init

   ! thick_ox_skin(skin, height, layer)
   real, dimension(mskin, 100, moxide_layer)  :: thick_ox_skin, rad_ox_skin_large, &
      & rad_ox_skin_small, fraction_skin_on

   ! character(LEN = 80)   ::  exfoliation_layers
   ! exfoliation_layers(1) = 3 - magnetite, for two layers
   ! exfoliation_layers(1) = 4 - haematite, exfoliation_layers(2) = 3 - magnetite
   ! exfoliation_layers(1) = 4 - haematite, exfoliation_layers(2) = 3 - magnetite and exfoliation_layers(3) = 2 - spinel
   integer, dimension(10)  :: exfoliation_layers

   ! no_exfoliation_layers - data points in exfoliation_layers
   integer   ::  no_exfoliation_layers, start_exfoliation_layers, &
       & end_exfoliation_layers

   ! energy used in Armitt's formulation would be based on energy from layers 
   ! energy_exfoliation_start_layer to energy_exfoliation_end_layer
   ! for Armitt energy_exfoliation_start_layer = 2 energy_exfoliation_end_layer = 3
   ! spinel and magnetite
   ! start_energy_exfoliation_layer - for each exfoliating layer
   ! [start_exfoliation_layers: end_exfoliation_layers] the 
   ! energy would be evaluated using [start_energy_exfoliation_layer: end_ ..]
   integer, dimension(10)   ::  start_energy_exfoliation_layer, &
       & end_energy_exfoliation_layer

   ! exfoliate_fr_skin_outage(height, outage) - how much fraction of original oxide exfoliates
   ! fraction_skin_remain_outage(height, outage) - how much fraction of original oxide remains on the tube
   real, dimension(100, 0:mramp_out) :: exfoliate_fr_skin_outage, &
       & fraction_skin_remain_outage
   ! dox_before_exfol(layer, height, outage) layer = 1 - all scale, 2 - magnetite
   real, dimension(2, 100, 0:mramp_out) :: dox_before_exfol, dox_after_exfol

   character(LEN = 80), dimension(15)   :: area_fr_vs_volume_deposit_name
   ! area_fraction_blocked_input at actual deposit volumes of 
   ! volume_reference_input * volume_dimless_deposit_input
   ! area_fraction_blocked_input vs volume_dimless_deposit_input
   ! was obtained for a bend_radius = 2.5 * OD and ID = 0.75 * OD
   real, dimension(15, 20) :: area_fraction_blocked_input, &
      & volume_dimless_deposit_input
   ! better to input volume_reference_input as 1.25*pi*ID^3 as the bend 
   ! cannot be easily identified to get the ID at the bend
   ! bend location could be at length = loop_transition_length*0.5 for loop 1 and
   !  loop_transition_length*1.5 for loop 2
   ! the volume_dimless_deposit_input and area_fr_vs_volume_deposit_name should depend on the loop; 
   ! use the loop name when defining these quantities
   ! volume_reference_input(deposit, function); per each loop should be no_length_deposit_tube simulated
   real, dimension(4, 15) :: volume_reference_input 
   integer, dimension(15) ::  nmax_area_fraction_blocked_input

   ! no_area_fraction_blocked_input = numbers of datasets for 
   ! area_fraction_blocked_input vs volume_dimless_deposit_input
   integer             :: no_area_fraction_blocked_input

   ! deposit_bend_loop_type must be one of area_fr_vs_volume_deposit_name
   character(LEN = 80), dimension(4) :: deposit_bend_loop_type

   ! id_deposit_bend_loop_type(length_deposit id) points to "i" which is the dataset 
   ! for (i, :) area_fraction_blocked_input
   integer, dimension(4) :: id_deposit_bend_loop_type

   ! radius_deposit_tube_loop - radius of the deposit in each loop
   ! radius_deposit_tube - radius at each 90 degree turn in the bend
   real, dimension(2) :: radius_deposit_tube_loop

   ! length_deposit_over_id_tube = L_deposit/ ID_deposit, the B factor
   real, dimension(10)  :: length_deposit_over_id_tube, volume_reference_deposit
 
 END MODULE BLOCKAGE_DATA_MODULE
