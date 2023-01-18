 MODULE BLOCKAGE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Input all the data for the blockage - height dependence
  !
  ! Author: Adrian S. Sabau, sabaua@ornl.gov
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: BLOCKAGE_INPUT, BLOCKAGE_AREA

 CONTAINS

 SUBROUTINE BLOCKAGE_INPUT(fatal)

  use parameter_module,        only : mstep, zero, pi
  use input_utilities_module,   only: SEARCH_NML
  use output_module,            only: tty_lun, out_lun, inp_lun
  use blockage_data_module,     only: area_fraction_spalled, &
      & energy_elastic_oxide, deposit_fluid_density, &
      & no_area_fr_spalled, deposit_porosity_fraction, &
      & accumulation_deposit_area, non_blocking_vol_fraction, &
      & no_deposit_sites, mblock, radius_deposit_cylinder, &
      & if_exfoliation_volume, if_blockage_area, deposit_oxide_density
  use blockage_data_module,     only: length_deposit_tube, &
      & number_deposits_along_tube, no_deposits_along_tube, &
      & radius_deposit_tube, no_length_deposit_tube, &
      & growth_mode_after_exfoliation, thickness_new_ox_growth, &
      & if_exfoliation_oxide_growth, start_exfoliation_layers, &
      & end_exfoliation_layers, energy_armitt_factor_scale, &
      & start_energy_exfoliation_layer, end_energy_exfoliation_layer, &
      & energy_at_zero_area_fraction_spalled, &
      & area_fr_vs_volume_deposit_name, area_fraction_blocked_input, &
      & volume_dimless_deposit_input, volume_reference_input, &
      & no_area_fraction_blocked_input, deposit_bend_loop_type, &
      & id_deposit_bend_loop_type, nmax_area_fraction_blocked_input, &
      & length_deposit_over_id_tube
     
 use waterwall_data_module,    only: length_tube, radius_outer_tube, &
      & thickness_tube, orientation_tube, &
      & no_length_tube
  use property_module, only:  LINEAR_PROPERTY_SIMPLE

  implicit none

  ! Argument List
  logical, intent(INOUT) :: fatal

    ! Local Variables
    integer :: ioerror, l, i, k, j, no_block1
    integer, dimension(mblock) :: id_block_area
    logical :: no_solv_namelist, solv_namelist
    real    :: factor1, energy_check, area_exfoliated_check
  real    :: rad_min
  integer, dimension(6)  :: num_deposits_along_tube_in =(/1, 2, 3, 5, 10, 20/)
  real, dimension(20) :: reorder_array

    ! default armitt data p. 4-21 in the report; the last point was extrapolated
    !  real, dimension(mblock)  gives an error
    ! real, dimension(30)   :: energy_elastic_oxide_armitt = (/5.9060, &
    !    & 7.7050, 8.7210, 9.3470, 9.5820, 9.8160, 10.129, 10.833, &
    !    & 12.788, 27.256, 36.094, 110.0/) ! , 18*0.0/)
    !real, dimension(30)   :: area_fraction_spalled_armitt =(/ 0.0, &
    !    & 0.0430, 0.114, 0.270, 0.440, 0.466, 0.485, 0.506, 0.519, &
    !    & 0.598, 0.633, 1.0/) ! , 18*0.0/)

   namelist /blockage/  area_fraction_spalled, energy_elastic_oxide, &
      & deposit_porosity_fraction, accumulation_deposit_area, &
      & non_blocking_vol_fraction, radius_deposit_cylinder, &
      & if_exfoliation_volume, if_blockage_area, &
      & growth_mode_after_exfoliation, thickness_new_ox_growth, &
      & if_exfoliation_oxide_growth, start_exfoliation_layers, &
      & end_exfoliation_layers, energy_armitt_factor_scale, &
      & start_energy_exfoliation_layer, end_energy_exfoliation_layer, &
      & energy_at_zero_area_fraction_spalled, &
      & deposit_fluid_density, deposit_oxide_density, &
      & energy_at_zero_area_fraction_spalled, &
      & area_fr_vs_volume_deposit_name, area_fraction_blocked_input, &
      & volume_dimless_deposit_input, volume_reference_input, &
      & deposit_bend_loop_type, length_deposit_over_id_tube

     if_exfoliation_volume = .false.
     if_blockage_area = .false.

     non_blocking_vol_fraction = 0.0

     deposit_porosity_fraction = 0.75  ! Armitt report, p. 6-3

     ! 5180.0e+3 - average density of the deposits; use averge density of the '(Fe,Cr)3O4', and 'Fe3O4'
     deposit_oxide_density = 0.0

     deposit_fluid_density = 1000.0 ! density of water

     accumulation_deposit_area = 0.0 ! chosen based on bends

     no_deposits_along_tube = 1
     number_deposits_along_tube(1:no_deposits_along_tube) = &
      & num_deposits_along_tube_in(1:no_deposits_along_tube)

     no_length_deposit_tube = 10   ! as in Armitt
     no_length_deposit_tube = 4
     no_length_deposit_tube = 0  ! when length_deposit_over_id_tube is used

     growth_mode_after_exfoliation = 'continuous'
     thickness_new_ox_growth = 0.0 ! equivalent  growth_mode_after_exfoliation = 'zero'

     if_exfoliation_oxide_growth = .false.

     start_exfoliation_layers = 3
     end_exfoliation_layers = 3

     ! this is default for each oxide layer
     do i = 2, 3
       start_energy_exfoliation_layer(i) = i
       end_energy_exfoliation_layer(i) = i
     end do

     ! this is default for Armitt's energy over the entire oxide scale
     start_energy_exfoliation_layer = 2
     end_energy_exfoliation_layer = 3

     energy_armitt_factor_scale = 1.0

     energy_at_zero_area_fraction_spalled = 0.0

     do i = 1, 15
       area_fr_vs_volume_deposit_name(i) = 'none'
     end do

     area_fraction_blocked_input = 0.0
     volume_dimless_deposit_input = 0.0

     do i = 1, 4
       deposit_bend_loop_type(i) = 'none'
     end do

       ! Find namelist
       no_solv_namelist = .false.
       call SEARCH_NML (inp_lun, no_solv_namelist, 'blockage', 'BLOCKAGE')
       solv_namelist = .NOT. no_solv_namelist

       if (solv_namelist) then
          read (inp_lun, NML= blockage, IOSTAT=ioerror)
          fatal = .not. (ioerror == 0) ! If read error, then didn't read namelist
          if (.not. fatal)  then
            write (tty_lun, 15)
            write (out_lun, 15)
15          format (/' DONE Reading BLOCKAGE Namelist ...')
          else
            write(tty_lun, *) 'STOP: BLOCKAGE Input'
            stop
          endif
       else
         fatal = .false.
         write(tty_lun, *) 'No BLOCKAGE Namelist: if_blockage_area = F, if_exfoliation_volume = F'
         write(out_lun, *) 'No BLOCKAGE Namelist: if_blockage_area = F, if_exfoliation_volume = F'
         RETURN
       end if

   if (if_blockage_area .and. (.not. if_exfoliation_volume))  then
     write(tty_lun, *) 'if_blockage_area = T, set if_exfoliation_volume to T'
     STOP
   endif

   if (.not. if_exfoliation_volume)  then
     write(tty_lun, *) 'if_blockage_area = F, if_exfoliation_volume = F'
     RETURN
   endif

   no_area_fr_spalled = 0
   do i = mblock, 1, -1
     if (no_area_fr_spalled == 0 .and. energy_elastic_oxide(i) /= 0.0) &
       & no_area_fr_spalled = i
   end do

   if (no_area_fr_spalled <=0)  then
     write(tty_lun, *) 'WARNING: no_area_fr_spalled  must be > 0; use DEFAULT'
     STOP
     no_area_fr_spalled = 12
     ! energy_elastic_oxide(1: no_area_fr_spalled)= energy_elastic_oxide_armitt(1: no_area_fr_spalled)
     ! area_fraction_spalled(1: no_area_fr_spalled) = area_fraction_spalled_armitt(1: no_area_fr_spalled)
   else

     if (ABS(energy_armitt_factor_scale - 1.0) <= 1.0e-2)  then

       write(out_lun, *) 'set default for Armitts energy over the entire oxide scale'
       write(out_lun, *) 'set default for Armitts start_energy_exfoliation_layer = 2'
       write(out_lun, *) 'set default for Armitts end_energy_exfoliation_layer = 3'
       write(out_lun, *) 'if other values are needed then code it'

       do i = 2, 3

         if (start_energy_exfoliation_layer(i) /= 2)  then

           write(out_lun, *) 'start_energy_exfoliation_layer= 2 from ', &
             & start_energy_exfoliation_layer(i) 
           start_energy_exfoliation_layer(i) = 2

         endif
         
         if (end_energy_exfoliation_layer(i) /= 3)  then

           write(out_lun, *) 'end_energy_exfoliation_layer= 3 from ', &
             & end_energy_exfoliation_layer(i) 
           end_energy_exfoliation_layer(i) = 3

         endif

       end do
       
     endif

     if (ABS(energy_armitt_factor_scale - 1.0) > 1.0e-2)  then
       write(out_lun, *) 'armit_energy_elastic_oxide changed by a factor of ', &
        & energy_armitt_factor_scale
     endif

     do i = 1, no_area_fr_spalled
       write(out_lun, 5) energy_elastic_oxide(i), &
        & energy_armitt_factor_scale *  &
        & (energy_elastic_oxide - energy_at_zero_area_fraction_spalled), &
        & area_fraction_spalled(i), &
        & energy_at_zero_area_fraction_spalled, &
        & energy_armitt_factor_scale * energy_at_zero_area_fraction_spalled, &
        & energy_armitt_factor_scale * energy_elastic_oxide(i)
 5     format('armitt_area_fr_vs_energy_input ', 5(1pe13.6, 1x))
     end do

     ! energy_elastic_oxide = energy_at_zero_area_fraction_spalled + &
     !        & energy_armitt_factor_scale * energy_elastic_oxide
     energy_elastic_oxide = energy_armitt_factor_scale *  &
             & (energy_elastic_oxide - energy_at_zero_area_fraction_spalled)

     ! testing the energy 
     do i = 1, no_area_fr_spalled -1
       energy_check = 0.5 * (energy_elastic_oxide(i) + energy_elastic_oxide(i+1))
       call LINEAR_PROPERTY_SIMPLE(energy_check, &
         & no_area_fr_spalled-1, 1, mblock, &
         & area_fraction_spalled, energy_elastic_oxide, area_exfoliated_check)
       write(out_lun, *) 'armitt_area_fr_vs_energy_check ', &
         & energy_check, area_exfoliated_check
     end do

   endif

   if (if_blockage_area)  then

     no_deposit_sites = 0
     do i = mblock, 1, -1
       if (no_deposit_sites == 0 .and. accumulation_deposit_area(i) /= 0.0) &
         & no_deposit_sites = i
     end do

     if (no_deposit_sites <=0)  then
       write(tty_lun, *) 'WARNING: no_deposit_sites  must be > 0; use DEFAULT'
       ! STOP

       ! should get the radius in the bend; when the orientation changes from v to h
       no_block1 = 0
  do i = 2, no_length_tube -1
    if (TRIM(orientation_tube(i))/= TRIM(orientation_tube(i-1)))  then
      no_block1 = no_block1 + 1
      id_block_area(no_block1) = i
    endif
  end do

     if (no_block1 == 0)  then
       write(tty_lun, *) 'ERROR: no_deposit_sites is 0; accumulation_deposit_area must be > 0'
       STOP
     endif

  do i = 1, no_block1

    radius_deposit_cylinder(i) = radius_outer_tube(id_block_area(i)) - &
       & thickness_tube(id_block_area(i))
    accumulation_deposit_area(i) = pi * &
       & (radius_outer_tube(id_block_area(i)) - &
       & thickness_tube(id_block_area(i)))**2
    write(out_lun, *) 'acumul_area_calc ', i, &
       & radius_outer_tube(id_block_area(i)) - &
       & thickness_tube(id_block_area(i)), accumulation_deposit_area(i)

  end do

  no_deposit_sites = no_block1

   else if (no_deposit_sites > 0)  then

  do i = 1, no_deposit_sites
    write(out_lun, *) 'acumul_area_read ', i, accumulation_deposit_area(i)
  end do

   endif

     rad_min = 10.0

  do i = 1, no_deposit_sites

    if (rad_min > radius_deposit_cylinder(i))  then
      rad_min = radius_deposit_cylinder(i)
      ! here the transversal area does not work well as filling of the gap
      ! occurs by piling up of scales inside a cylinder
      ! here the area would be A=L * segment(filled_transversal_distance)
      ! Armitt uses L as a parameter L=1 to 100 cm for r=6.5 cm
      ! us L/r = 0.2 to 10
    endif

  end do

     ! the next are if the entire exfoliation for one loop was 
     ! accumulated not in the corner but one location; which radius
     ! use the minimum radius for the next ones; 
     ! this assignment would not work for two bends; each bend would have a different radius
       radius_deposit_tube = rad_min

     ! length_deposit_tube will be set later on when the tube loop dimensions will be available
     length_deposit_tube = 0.0

   do i = 10, 1, -1
     if (length_deposit_over_id_tube(i) /= 0 .and. no_length_deposit_tube ==0) &
     no_length_deposit_tube = i
   end do

   if (no_length_deposit_tube == 0)  then
     no_length_deposit_tube = 2  ! this is per each loop; check if 2 or 4 over 2 loops
     length_deposit_over_id_tube(1)= 2.5  ! to illustrate the radius rather than ID
     length_deposit_over_id_tube(2) = 5.0  ! based on EPRI pictures
   endif

     ! calculated with minimum radiius in the bend; 
     ! (1) use ID instead of the radius
     ! (2) later on correct to use the radius in the bend not the minimum one
     ! later on deal with the second loop as Blockage namelist the no_loops is not available
     write(out_lun, 3) (length_deposit_over_id_tube(i), i = 1, no_length_deposit_tube)
 3   format('dless_length_deposit_tube_input ', 20(1pe13.6, 1x))

    do i = 1, no_length_deposit_tube
       factor1 = EXP((LOG(0.2) * (no_length_deposit_tube-i) + &
              & LOG(10.0)*(i-1)) / &
              & (no_length_deposit_tube - 1))
       length_deposit_tube(i) = radius_deposit_tube * factor1
     end do

     length_deposit_tube(5) = 2.0 * radius_deposit_tube
     ! the above was for armitt
     ! old onelength_deposit_tube
     ! length_deposit_tube(1:4) = (/1.0, 2.0, 5.0, 10.0/)
     ! length_deposit_tube(1:4) = length_deposit_tube(1:4) * 2.0 * radius_deposit_tube
     ! this is overriden as radius_deposit_tube could be different for each loop
     length_deposit_tube = 0.0

   endif 

   if (.not. (TRIM(growth_mode_after_exfoliation) == 'continuous' .or. &
     & TRIM(growth_mode_after_exfoliation) == 'cont' .or. &
     & TRIM(growth_mode_after_exfoliation) == 'zero'  .or. &
     & TRIM(growth_mode_after_exfoliation) == 'start_first_exfolation'  .or. &
     & TRIM(growth_mode_after_exfoliation) == 'from_given_thickness'  .or. &
     & TRIM(growth_mode_after_exfoliation) == 'CONTINUOUS'))  then
     write(6, *) 'ERROR: growth_mode_after_exfoliation must be: ', &
     & 'continuous, cont, zero, start_first_exfolation, OR ', &
     & 'from_given_thickness'
     STOP

   endif 

   if (.not. if_exfoliation_volume)  then 

     if (if_exfoliation_oxide_growth)  then

       if_exfoliation_oxide_growth = .false.

     endif

   endif

   if (.not. if_exfoliation_oxide_growth)  then

     start_exfoliation_layers = 0
     end_exfoliation_layers = 0

   endif

   no_area_fraction_blocked_input = 0
   do i = 15, 1, -1
     if (no_area_fraction_blocked_input == 0 .and. &
       &  TRIM(area_fr_vs_volume_deposit_name(i)) /= 'none') &
       & no_area_fraction_blocked_input = i
   end do

   if (no_area_fraction_blocked_input == 0)  then

     no_area_fraction_blocked_input = 1
     nmax_area_fraction_blocked_input(1)= 2
     area_fr_vs_volume_deposit_name(1) = 'horizontal_tube'
     area_fraction_blocked_input(1, 1) = 0.0
     volume_dimless_deposit_input(1, 1) = 0.0
     area_fraction_blocked_input(1, 2) = 1.0
     volume_dimless_deposit_input(1, 2) = 1.0

     ! later set up volume_reference_input to 5*ID for each loop

     ! later on set these per each loop by default; right now the no_loops is unknown
     id_deposit_bend_loop_type(1:no_length_deposit_tube) = 1
     do i = 1, no_length_deposit_tube
       deposit_bend_loop_type(i) = 'horizontal_tube'
     end do

   else

     do i = 1, no_area_fraction_blocked_input

       nmax_area_fraction_blocked_input(i) = 0
       do j = 20, 1, -1
         if (nmax_area_fraction_blocked_input(i) == 0 .and. &
           & volume_dimless_deposit_input(i, j) > 1.0e-12)  then 
           nmax_area_fraction_blocked_input(i) = j
         endif
       end do

       if (nmax_area_fraction_blocked_input(i) > 0)  then

         ! make sure that the x array is monotonic increasing
         if (volume_dimless_deposit_input(i, 1) > &
           & volume_dimless_deposit_input(i, 2))  then
           ! reorder
           do j = nmax_area_fraction_blocked_input(i), 1, -1
             reorder_array(j) = volume_dimless_deposit_input(i, &
                & nmax_area_fraction_blocked_input(i) - j + 1)
           enddo
           volume_dimless_deposit_input(i, :) = reorder_array(:)

           do j = nmax_area_fraction_blocked_input(i), 1, -1
             reorder_array(j) = area_fraction_blocked_input(i, &
                & nmax_area_fraction_blocked_input(i) - j + 1)
           enddo
           area_fraction_blocked_input(i, :) = reorder_array(:)

         endif

       endif

       write(out_lun, *) 'af_vol_input ', i, nmax_area_fraction_blocked_input(i), &
         & volume_reference_input(1, i), volume_reference_input(2, i)

       do j = 1, nmax_area_fraction_blocked_input(i)
         write(out_lun, 12) TRIM(area_fr_vs_volume_deposit_name(i)), &
           & area_fraction_blocked_input(i, j), volume_dimless_deposit_input(i, j)
 12      format('af_vol_input ', (a), 1x, 20(1pe13.6, 1x))
       end do

     end do

   endif

   ! check later on deposit_bend_loop_type == area_fr_vs_volume_deposit_name
   ! as at input the number of loops may not be known in this subroutine

  RETURN

  END SUBROUTINE BLOCKAGE_INPUT

 SUBROUTINE BLOCKAGE_AREA(N)

  use parameter_module,         only: mstep, zero
  use output_module,            only: tty_lun, out_lun, inp_lun, enrg_lun
  use blockage_data_module,     only: deposit_porosity_fraction, &
      & accumulation_deposit_area, non_blocking_vol_fraction, &
      & no_deposit_sites, mblock, volume_exfoliated_bend, &
      & area_fraction_blocked_bend, height_deposit_bend, &
      & radius_deposit_cylinder, volume_deposit_2cyl, volume_deposit_cube, &
      & volume_deposit_bend, volume_deposit_1cyl, &
      & area_fraction_blocked_cube, volume_fraction_blocked_cube, &
      & height_deposit_cube, area_fraction_blocked_2cyl,  &
      & volume_fraction_blocked_2cyl, height_deposit_2cyl, &
      & area_fraction_blocked_1cyl,  &
      & volume_fraction_blocked_1cyl, height_deposit_1cyl, &
      & volume_fraction_blocked_bend
  use blockage_data_module,     only: volume_exfoliated_bend_porous, &
      & if_blockage_area, volume_deposit_tube, area_deposit_tube, &
      & area_fraction_blocked_tube, area_fraction_blocked_tube0, &
      & no_deposits_along_tube, length_deposit_tube, radius_deposit_tube, &
      & number_deposits_along_tube, no_length_deposit_tube, &
      & exfoliate_fr_skin_outage, fraction_skin_remain_outage, &
      & dox_before_exfol, dox_after_exfol, id_deposit_bend_loop_type, &
      & area_fr_vs_volume_deposit_name, no_area_fraction_blocked_input, &
      & deposit_bend_loop_type, id_deposit_bend_loop_type, &
      & nmax_area_fraction_blocked_input, &
      & volume_dimless_deposit_input, area_fraction_blocked_input, &
      & radius_deposit_tube_loop, volume_reference_deposit, &
      & length_deposit_over_id_tube
  use solution_data_module,     only: time_event, no_total_outage, &
      & no_f2l_event_out
  use waterwall_data_module, only:  no_loops, no_height_oxide_output, &
      & length_tube_oxide_output, loop_transition_length

  implicit none

  ! Argument List
  integer, intent(IN) :: N

  ! Local variables
  integer :: ibend, jlayer, ioutage, j, nodep, ileng, i, k, iloop, idep, j1
  real    :: delta_height, area_total_deposit_tube, area_fr_check1, &
          & volume_exfoliated_now, inner_diam_bend
  real    :: pi
  character(LEN = 80), dimension(20) :: ch_repeat, ch_repeat0
  character(LEN = 1), dimension(0:9)   :: ch_index = (/'0','1','2','3','4','5', &
                  & '6','7','8','9'/)
  logical, save :: if_first_out1 = .true.

  if (.not. if_blockage_area)  RETURN

  pi = 4.0 * ATAN(1.0)

  ! bend = short for deposit site; accumulation_deposit_area is given
  ! no need for each deposit site since at each location the exoliated mass
  ! will accumulate at only one bend

  ! cube = bend geometry is equivalent to that of a cube

  ! 2cyl = bend geometry is exactly that given by the intersection of 2 cylinders

  ! 1cyl = bend geometry is that given by a cylinder

do ibend = 1, no_deposit_sites

  ! porous deposits, volume_exfoliated_bend is fully solid, i.e., 0 porosity
  volume_exfoliated_bend_porous(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) * &
     & deposit_porosity_fraction(ibend)

  volume_deposit_bend(ibend) = 2.0 * &
     & radius_deposit_cylinder(ibend) * &
     & accumulation_deposit_area(ibend)

  volume_deposit_1cyl(ibend) = 2.0 * pi * &
     & radius_deposit_cylinder(ibend)**3

  volume_deposit_2cyl(ibend) = (16.0/3.0) * &
     & radius_deposit_cylinder(ibend)**3

  volume_deposit_cube = volume_deposit_2cyl

  volume_fraction_blocked_bend(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & volume_deposit_bend(ibend) 

  volume_fraction_blocked_cube(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & volume_deposit_cube(ibend) 

  volume_fraction_blocked_1cyl(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & volume_deposit_1cyl(ibend) 

  volume_fraction_blocked_2cyl(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & volume_deposit_2cyl(ibend) 

  ! get the heights

  height_deposit_bend(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & accumulation_deposit_area(ibend) 

  height_deposit_cube(ibend, :, :) = &
     & volume_exfoliated_bend(ibend, :, :) / &
     & ((1.7472 * radius_deposit_cylinder(ibend))**2)

  ! height_deposit_2cyl(ibend, :, :) = &

end do

  write(out_lun, 2) (volume_deposit_cube(j), volume_deposit_bend(j), &
     & volume_deposit_1cyl(j), volume_deposit_2cyl(j), j=1, no_deposit_sites)
  write(enrg_lun, 2) (volume_deposit_cube(j), volume_deposit_bend(j), &
     & volume_deposit_1cyl(j), volume_deposit_2cyl(j), j=1, no_deposit_sites)
 2 format('volume_in_bend ', 10(1pe13.6, 1x))

  do ioutage = 1, no_f2l_event_out  !  no_total_outage

  do ibend = 1, no_deposit_sites

    do jlayer = 1, N

      if (volume_fraction_blocked_cube(ibend, ioutage, jlayer) < &
        & 1.0)  then
        ! assuming an equivalent cube such that a = 1.7472 r to give V2=16/3 * r^3=a^3
        area_fraction_blocked_cube(ibend, ioutage, jlayer) = &
          & volume_fraction_blocked_cube(ibend, ioutage, jlayer) 
      else 
        area_fraction_blocked_cube(ibend, ioutage, jlayer) = 1.0
      endif

    end do

  end do

  ! the area fraction for the cube is equivalent to the volumetric fraction
  write(out_lun, 1) ioutage, time_event(ioutage), &
    & (volume_fraction_blocked_cube(ibend, ioutage, 1), &
    & volume_fraction_blocked_cube(ibend, ioutage, 3), &
    & volume_exfoliated_bend(ibend, ioutage, 1), &
    & volume_exfoliated_bend(ibend, ioutage, 3), &
    & ibend = 1, no_deposit_sites)
  write(enrg_lun, 1) ioutage, time_event(ioutage), &
    & (volume_fraction_blocked_cube(ibend, ioutage, 1), &
    & volume_fraction_blocked_cube(ibend, ioutage, 3), &
    & volume_exfoliated_bend(ibend, ioutage, 1), &
    & volume_exfoliated_bend(ibend, ioutage, 3), &
    & ibend = 1, no_deposit_sites)
 1 format('blockage_out1 ', i2, 1x, 10(1pe13.6, 1x))

  end do

  FIRST_IF: if (if_first_out1)  then

     ! old one for one loop
     ! length_deposit_tube(1:4) = (/1.0, 2.0, 5.0, 10.0/)

     ! calculated with minimum radius in the bend; 
     ! (1) use ID instead of the radius
     ! (2) later on correct to use the radius in the bend not the minimum one
    ! length_deposit_tube(1:2) = (/2.0, 5.0/)
    ! length_deposit_tube(3:4) = length_deposit_tube(1:2)


    if (no_loops == 1)  then

      ! no_length_deposit_tube = 2; desired to be 2
      ! radius_deposit_cylinder(1:2) refers to the first loop; start and end of bend
      radius_deposit_tube_loop(1) = radius_deposit_cylinder(1)
      length_deposit_tube(1:no_length_deposit_tube) =  2.0 * &
        & radius_deposit_tube_loop(1) * &
        & length_deposit_over_id_tube(1:no_length_deposit_tube)

    else if (no_loops == 2)  then

      ! no need for this; no_length_deposit_tube = 4

      ! no_loops is now known and additional setup is made
      do i = 1+no_length_deposit_tube, no_loops*no_length_deposit_tube
        if (TRIM(deposit_bend_loop_type(i)) == 'none')  then
          id_deposit_bend_loop_type(i) = 1
          deposit_bend_loop_type(i) = 'horizontal_tube'   
        endif
      end do

      ! radius_deposit_cylinder(1:2) refers to the first loop
      radius_deposit_tube_loop(1) = radius_deposit_cylinder(1)
      length_deposit_tube(1:no_length_deposit_tube) =  2.0 * &
        & radius_deposit_tube_loop(1) * &
        & length_deposit_over_id_tube(1:no_length_deposit_tube)

      ! radius_deposit_cylinder(3:4) refers to the second loop; start and end of bend
      radius_deposit_tube_loop(2) = radius_deposit_cylinder(3)
      ! length_deposit_tube(3:4) = length_deposit_tube(3:4) * 2.0 * radius_deposit_tube_loop(2)

      length_deposit_over_id_tube(1+no_length_deposit_tube: &
        & no_loops*no_length_deposit_tube) = &
        & length_deposit_over_id_tube(1:no_length_deposit_tube)

      length_deposit_tube(1+no_length_deposit_tube: &
        & no_loops*no_length_deposit_tube) =  &
        & 2.0 * radius_deposit_tube_loop(2) * &
        & length_deposit_over_id_tube(1+no_length_deposit_tube: &
        & no_loops*no_length_deposit_tube)

    else
      write(6, *) 'ERROR: no_loops= ', no_loops, ' < 1 or > 2 '
      STOP
    endif

    ! select the Af(Volume)
    id_deposit_bend_loop_type = 0

    if (TRIM(area_fr_vs_volume_deposit_name(1)) /= 'horizontal_tube')  then
      write(tty_lun, *) 'STOP: area_fr_vs_volume_deposit_name(1) must be horizontal_tube'
      STOP
    endif

    do j = 1, no_loops  ! ideally 
      ! loop 1 is represented by length_deposit_tube(1:2) 
      ! loop 2 is represented by length_deposit_tube(3:4)
      ! do j = 1, no_length_deposit_tube

      do k = 1, no_length_deposit_tube

        idep = no_length_deposit_tube*(j-1) + k

        do i = 1, no_area_fraction_blocked_input

          if (TRIM(deposit_bend_loop_type(idep)) == &
            & TRIM(area_fr_vs_volume_deposit_name(i)))  then

            id_deposit_bend_loop_type(idep) = i

            if (TRIM(deposit_bend_loop_type(idep)) == 'horizontal_tube')  then

              volume_reference_deposit(idep) = length_deposit_tube(idep) * pi * &
                & radius_deposit_tube_loop(j)**2

            else

              volume_reference_deposit(idep) = length_deposit_tube(idep) * pi * &
                & radius_deposit_tube_loop(j)**2
              ! volume_reference_deposit(idep) = 5 * 2 * radius_deposit_tube_loop(j) * pi * &
              !  & radius_deposit_tube_loop(j)**2

            endif 
                         
            write(out_lun, 12) TRIM(area_fr_vs_volume_deposit_name(i)),  &
              & radius_deposit_tube_loop(j), volume_reference_deposit(idep), &
              & volume_reference_deposit(idep)/ &
              & (pi * radius_deposit_tube_loop(j)**2), &
              & volume_dimless_deposit_input(i, j), area_fraction_blocked_input(i, j)
 12          format('length_radius_volume_length_selected ', (a), 1x, 20(1pe13.6, 1x))

            do j1 = 1, nmax_area_fraction_blocked_input(i)
              write(out_lun, 14) TRIM(area_fr_vs_volume_deposit_name(i)),  &
              & volume_dimless_deposit_input(i, j1), area_fraction_blocked_input(i, j1)
 14           format('af_vol_selected ', (a), 1x, 20(1pe13.6, 1x))
            end do

            EXIT
          endif

        end do

        ! if (id_deposit_bend_loop_type(2*(j-1) + k) == 0)  then
        if (id_deposit_bend_loop_type(idep) == 0)  then
          write(tty_lun, *) 'ERROR: deposit_bend_loop_type /= area_fr_vs_volume_deposit_name'
          write(tty_lun, *) 'ERROR: loop', idep, TRIM(deposit_bend_loop_type(idep)), &
             &  (TRIM(area_fr_vs_volume_deposit_name(i)), &
             & i = 1, no_area_fraction_blocked_input)
          STOP
        endif

      end do

    end do

    write(out_lun, 3) (length_deposit_tube(i), i = 1, &
      & no_loops * no_length_deposit_tube)
 3  format('length_deposit_tube_loop ', 20(1pe13.6, 1x))

    write(out_lun, 11) (length_deposit_tube(i)/(2.0*radius_deposit_tube_loop(1)), &
         & i = 1, no_length_deposit_tube), &
         & (length_deposit_tube(i)/(2.0*radius_deposit_tube_loop(2)), &
         & i = 1 + no_length_deposit_tube, no_loops*no_length_deposit_tube)
 11 format('length_fraction_deposit_ID_tube_loop ', 20(1pe13.6, 1x))

    do idep = 1, no_loops*no_length_deposit_tube

      k = id_deposit_bend_loop_type(idep)  !  id_deposit_bend_loop_type(i)
      do j = 1, nmax_area_fraction_blocked_input(k)
        write(out_lun, 13) idep, TRIM(area_fr_vs_volume_deposit_name(k)), &
           & volume_reference_deposit(idep), volume_dimless_deposit_input(k, j), &
           & area_fraction_blocked_input(k, j)
 13     format('af_vol_input_loop_final ', i1, 1x, (a), 20(1pe13.6, 1x))
      end do

    end do

    if_first_out1 = .false.

  endif FIRST_IF

  do ioutage = 1, no_f2l_event_out ! no_total_outage

    do jlayer = 1, N

      do nodep = 1, no_deposits_along_tube
                
        ! each deposit would have vol_exfol * fr_per_each deposit
        ! fr_per_each_deposit = 1.0 / number_deposits_along_tube(nodep)

        do ileng = 1, no_loops * no_length_deposit_tube

          if (ileng <= no_length_deposit_tube) iloop = 1   ! ileng 1 and 2 refers to the first loop
          if (ileng >= no_length_deposit_tube + 1) iloop = 2   ! ileng 3 and 4 to the second loop

          area_total_deposit_tube = pi * radius_deposit_tube_loop(iloop)**2

          ! to be consistent with the volume_all1 and volume_mag1 from exfoliation_modu
          if (no_loops == 1)  then

            ! ileng 1 and 2 refers to the first loop
            volume_deposit_tube(ioutage, jlayer, ileng, nodep) = &
              & volume_deposit_tube(ioutage, jlayer, 1, 1) / &
              & number_deposits_along_tube(nodep)

          else if (no_loops == 2)  then

            if (ileng <= no_length_deposit_tube)  then

              ! ileng 1 and 2 refers to the first loop
              volume_deposit_tube(ioutage, jlayer, ileng, nodep) = &
                & volume_deposit_tube(ioutage, jlayer, 1, 6) / &
                & number_deposits_along_tube(nodep)

            else if (ileng >= 1 + no_length_deposit_tube)  then

              ! ileng 3 and 4 to the second loop
              volume_deposit_tube(ioutage, jlayer, ileng, nodep) = &
                & volume_deposit_tube(ioutage, jlayer, 1, 7) / &
                & number_deposits_along_tube(nodep)

            endif

          endif

          ! without porosity
          area_deposit_tube(ioutage, jlayer, ileng, nodep) = &
            & volume_deposit_tube(ioutage, jlayer, ileng, nodep) / &
            & length_deposit_tube(ileng)

          area_fraction_blocked_tube0(ioutage, jlayer, ileng, nodep) = &
            & area_deposit_tube(ioutage, jlayer, ileng, nodep) / &
            & area_total_deposit_tube

          ! new; keep the above only for reference
          call AREA_BLOCKED_VS_VOLUME(iloop, id_deposit_bend_loop_type(ileng), &
            & volume_deposit_tube(ioutage, jlayer, ileng, nodep), &
            & volume_reference_deposit(ileng), area_fr_check1)

          if (jlayer == 3)  then
            write(out_lun, 17) ioutage, length_deposit_tube(ileng), &
              & area_total_deposit_tube * length_deposit_tube(ileng) / &
              & volume_reference_deposit(ileng), &
              & area_fraction_blocked_tube0(ioutage, jlayer, ileng, nodep), &
              & area_fr_check1, &
              & volume_deposit_tube(ioutage, jlayer, ileng, nodep) / &
              & volume_reference_deposit(ileng), &
              & volume_deposit_tube(ioutage, jlayer, ileng, nodep), &
              & volume_reference_deposit(ileng)
      17    format('check_area ', i2, 10(1pe13.6))
          endif

          ! after the printout
          area_fraction_blocked_tube0(ioutage, jlayer, ileng, nodep) = &
            & area_fr_check1

          ! with porosity
          area_deposit_tube(ioutage, jlayer, ileng, nodep) = &
            & volume_deposit_tube(ioutage, jlayer, ileng, nodep) / &
            & (length_deposit_tube(ileng) * (1.0 - deposit_porosity_fraction(1)))

          area_fraction_blocked_tube(ioutage, jlayer, ileng, nodep) = &
            & area_deposit_tube(ioutage, jlayer, ileng, nodep) / &
            & area_total_deposit_tube

          ! new; keep the above only for reference
          call AREA_BLOCKED_VS_VOLUME(iloop, id_deposit_bend_loop_type(ileng), &
            & volume_deposit_tube(ioutage, jlayer, ileng, nodep) / &
            & (1.0 - deposit_porosity_fraction(1)), &
            & volume_reference_deposit(ileng), &
            & area_fraction_blocked_tube(ioutage, jlayer, ileng, nodep))

          if (area_fraction_blocked_tube(ioutage, jlayer, ileng, nodep) >= &
            &  1.0)  then

            area_fraction_blocked_tube(ioutage, jlayer, ileng, nodep) = 1.0

          endif

        end do

      end do

    end do

  end do

  ! output now
    if (no_loops == 1)  then

      ch_repeat(1) = ' afb_tot afb_magL1'
      ch_repeat(2:no_length_deposit_tube) = ch_repeat(1)
      ch_repeat0(1) = ' afb0_tot afb0_magL1'
      ch_repeat0(2:no_length_deposit_tube) = ch_repeat0(1)
      write(out_lun, 4) (ch_repeat(j)(1:18)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j), j = 1, &
        & no_loops*no_length_deposit_tube)
      write(enrg_lun, 4) (ch_repeat(j)(1:18)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j), j = 1, &
        & no_loops*no_length_deposit_tube)

    else if (no_loops == 2)  then

      ch_repeat(1:no_length_deposit_tube) = ' afb_tot afb_magL1'
      ch_repeat0(1:no_length_deposit_tube) = ' afb0_tot afb0_magL1'
      ch_repeat(1+no_length_deposit_tube:2*no_length_deposit_tube) = ' afb_tot afb_magL2'
      ch_repeat0(1+no_length_deposit_tube:2*no_length_deposit_tube) = ' afb0_tot afb0_magL2'
      write(out_lun, 4) (ch_repeat(j)(1:18)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat(j)(1:18)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube)
      write(enrg_lun, 4) (ch_repeat(j)(1:18)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat(j)(1:18)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat0(j)(1:20)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube)

    endif

 4  format('#blockage_shd_lengthdep outage time_outage', 20(a))

  do ioutage = 1, no_f2l_event_out ! no_total_outage

    write(out_lun, 5) ioutage, time_event(ioutage), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &  ! all layers
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &   ! only magnetite layer
      & ileng = 1, no_loops*no_length_deposit_tube), &
      & (area_fraction_blocked_tube0(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube0(ioutage, 3, ileng, 1), &
      & ileng = 1, no_loops*no_length_deposit_tube)
    write(enrg_lun, 5) ioutage, time_event(ioutage), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &
      & ileng = 1, no_loops*no_length_deposit_tube), &
      & (area_fraction_blocked_tube0(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube0(ioutage, 3, ileng, 1), &
      & ileng = 1, no_loops*no_length_deposit_tube)
 5 format('blockage_shd_lengthdep ', i2, 1x, 50(1pe13.6, 1x))

  end do

    if (no_loops == 1)  then

      ch_repeat(1) = ' a_tot a_magL1'
      ch_repeat(2:no_length_deposit_tube) = ch_repeat(1)
      write(out_lun, 7) (ch_repeat(j)(1:14), j = 1, no_length_deposit_tube)
      write(enrg_lun, 7) (ch_repeat(j)(1:14), j = 1, no_length_deposit_tube)

    else if (no_loops == 2)  then

      ch_repeat(1:no_length_deposit_tube/2) = ' a_tot a_magL1'
      ch_repeat(1+no_length_deposit_tube/2:no_length_deposit_tube) = ' a_tot a_magL2'
      write(out_lun, 7) (ch_repeat(j)(1:14)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat(j)(1:14)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube)
      write(enrg_lun, 7) (ch_repeat(j)(1:14)//'_d'//ch_index(j), &
        &  j = 1, no_length_deposit_tube), &
        & (ch_repeat(j)(1:14)//'_d'//ch_index(j-no_length_deposit_tube), &
        &  j = 1+no_length_deposit_tube, 2*no_length_deposit_tube)

    endif

 7  format('area_deposit_length1 outage time_outage', 10(a))

  do ioutage = 1, no_f2l_event_out ! no_total_outage

    write(out_lun, 8) ioutage, time_event(ioutage), &
      & (area_deposit_tube(ioutage, 1, ileng, 1), &
      & area_deposit_tube(ioutage, 3, ileng, 1), &
      & ileng = 1, no_loops*no_length_deposit_tube)
    write(enrg_lun, 8) ioutage, time_event(ioutage), &
      & (area_deposit_tube(ioutage, 1, ileng, 1), &
      & area_deposit_tube(ioutage, 3, ileng, 1), &
      & ileng = 1, no_loops*no_length_deposit_tube)
 8 format('area_deposit_length1 ', i2, 1x, 10(1pe13.6, 1x))

  end do

  ! this is not fixed; but is not much of use as it will have only two lines and lots of columns
    ch_repeat(1) = ' afb_tot afb_mag'
    ch_repeat(2:10) = ch_repeat(1)
    write(out_lun, 9) (ch_repeat(j)(1:16), j = 1, no_f2l_event_out)
    write(enrg_lun, 9) (ch_repeat(j)(1:16), j = 1, no_f2l_event_out)
 9  format('#blockage_lengthdep_shd outage time_outage', 10(a))

  do ileng = 1, no_length_deposit_tube
    ! first loop
    write(out_lun, 10) ileng, length_deposit_tube(ileng) / &
      & (2.0*radius_deposit_tube_loop(1)), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &
      & ioutage = 1, no_f2l_event_out)
    write(enrg_lun, 10) ileng, length_deposit_tube(ileng) / &
      & (2.0*radius_deposit_tube_loop(1)), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &
      & ioutage = 1, no_f2l_event_out)
 10 format('blockage_lengthdep_shd ', i2, 1x, 50(1pe13.6, 1x))
  end do

  do ileng = 1 + no_length_deposit_tube, no_loops*no_length_deposit_tube
    ! second loop
    write(out_lun, 10) ileng, length_deposit_tube(ileng) / &
      & (2.0*radius_deposit_tube_loop(2)), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &
      & ioutage = 1, no_f2l_event_out)
    write(enrg_lun, 10) ileng, length_deposit_tube(ileng) / &
      & (2.0*radius_deposit_tube_loop(2)), &
      & (area_fraction_blocked_tube(ioutage, 1, ileng, 1), &
      & area_fraction_blocked_tube(ioutage, 3, ileng, 1), &
      & ioutage = 1, no_f2l_event_out)
  end do

  write(enrg_lun, 91) 'exfol_fr_skinoutage height ', &
    & (' efr_skin'//ch_index(j), &
    & ' skin1_on'//ch_index(j), j = 1, MIN(9, no_f2l_event_out)), &
    & (' efr_skin'//'1'//ch_index(j-10), ' skin1_on'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)
  write(out_lun, 91) 'exfol_fr_skinoutage height ', &
    & (' efr_skin'//ch_index(j), &
    & ' skin1_on'//ch_index(j), j = 1, MIN(9, no_f2l_event_out)), &
    & (' efr_skin'//'1'//ch_index(j-10), ' skin1_on'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)
 91 format(150(a))

   do j = 1, no_height_oxide_output
     write(enrg_lun, 92) length_tube_oxide_output(j), &
       & (exfoliate_fr_skin_outage(j, i), &
       & fraction_skin_remain_outage(j, i), i = 1, no_f2l_event_out)
     write(out_lun, 92) length_tube_oxide_output(j), &
       & (exfoliate_fr_skin_outage(j, i), &
       & fraction_skin_remain_outage(j, i), i = 1, no_f2l_event_out)
   end do
 92    format('exfol_fr_skinoutage ', 100(1pe13.6, 1x))

  write(enrg_lun, 93) 'dox_outage height ', (' sp'//ch_index(j), &
    & ' mag_b'//ch_index(j), ' mag_a'//ch_index(j), &
    & ' dox_b'//ch_index(j), ' dox_a'//ch_index(j), j = 1, MIN(9, no_f2l_event_out)), &
    & (' sp'//'1'//ch_index(j-10), ' mag_b'//'1'//ch_index(j-10), &
    & ' mag_a'//'1'//ch_index(j-10), ' dox_b'//'1'//ch_index(j-10), &
    & ' dox_a'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)
  write(out_lun, 93) 'dox_outage height ', (' sp'//ch_index(j), &
    & ' mag_b'//ch_index(j), ' mag_a'//ch_index(j), &
    & ' dox_b'//ch_index(j), ' dox_a'//ch_index(j), j = 1, MIN(9, no_f2l_event_out)), &
    & (' sp'//'1'//ch_index(j-10), ' mag_b'//'1'//ch_index(j-10), &
    & ' mag_a'//'1'//ch_index(j-10), ' dox_b'//'1'//ch_index(j-10), &
    & ' dox_a'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)
 93 format(150(a))

   do j = 1, no_height_oxide_output
     write(enrg_lun, 94) length_tube_oxide_output(j), &
       & (dox_before_exfol(1, j, i) - dox_before_exfol(2, j, i), &
       & dox_before_exfol(2, j, i), dox_after_exfol(2, j, i), &
       & dox_before_exfol(1, j, i), dox_after_exfol(1, j, i), &
       & i = 1, no_f2l_event_out)
     write(out_lun, 94) length_tube_oxide_output(j), &
       & (dox_before_exfol(1, j, i) - dox_before_exfol(2, j, i), &
       & dox_before_exfol(2, j, i), dox_after_exfol(2, j, i), &
       & dox_before_exfol(1, j, i), dox_after_exfol(1, j, i), &
       & i = 1, no_f2l_event_out)
   end do
 94    format('dox_outage ', 100(1pe13.6, 1x))

  RETURN

 END SUBROUTINE BLOCKAGE_AREA

 SUBROUTINE AREA_BLOCKED_VS_VOLUME(iloop, id_function, volume_deposit, &
     & volume_reference, area_fr_blocked_func)
     ! & volume_reference, area_fr_blocked_func)

  use blockage_data_module,     only: area_fraction_blocked_input, &
      & volume_dimless_deposit_input, nmax_area_fraction_blocked_input
  use property_module,    only: LINEAR_PROPERTY_SIMPLE

  implicit none

  integer, intent(IN) :: iloop, id_function
  real, intent(IN):: volume_deposit, volume_reference
  real, intent(OUT):: area_fr_blocked_func

  real :: volume_deposit_dimless, slope
  integer :: i, data_no

  volume_deposit_dimless = volume_deposit / volume_reference

  call LINEAR_PROPERTY_SIMPLE(volume_deposit_dimless, &
       & nmax_area_fraction_blocked_input(id_function), 1, 20, &
       & area_fraction_blocked_input(id_function, :), &
       & volume_dimless_deposit_input(id_function, :), &
       & area_fr_blocked_func)

 data_no = nmax_area_fraction_blocked_input(id_function)

 if (volume_deposit_dimless <= volume_dimless_deposit_input(id_function, 1))  then
   area_fr_blocked_func = area_fraction_blocked_input(id_function, 1)
   RETURN
 endif
 if (volume_deposit_dimless >= volume_dimless_deposit_input(id_function, data_no))  then
   area_fr_blocked_func = area_fraction_blocked_input(id_function, data_no)
   RETURN
 endif

 do i = 1, data_no - 1
   if (volume_deposit_dimless >= volume_dimless_deposit_input(id_function, i) .and. &
     & volume_deposit_dimless < volume_dimless_deposit_input(id_function, i+1))  then
     slope = & (area_fraction_blocked_input(id_function, i+1) - &
           &    area_fraction_blocked_input(id_function, i)) / &
           &  (volume_dimless_deposit_input(id_function, i+1) - &
           &   volume_dimless_deposit_input(id_function, i))
           
     area_fr_blocked_func = area_fraction_blocked_input(id_function, i) + &
           &  slope * (volume_deposit_dimless - volume_dimless_deposit_input(id_function, i))
     RETURN
   endif

 end do

 RETURN

 END SUBROUTINE AREA_BLOCKED_VS_VOLUME

END MODULE BLOCKAGE_MODULE


