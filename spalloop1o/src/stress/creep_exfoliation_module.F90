MODULE CREEP_EXFOLIATION_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   Reset the creep quantities after exfoliation 
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: CREEP_RESET_EXFOLIATION, STRESS_STRAIN_CREEP_OUT_EXFOL

 CONTAINS

  SUBROUTINE CREEP_RESET_EXFOLIATION(N, time_now, iheight, ioutage)

    !=======================================================================
    ! Purpose(s):
    !
    ! Reinitialize the creep strain based on the reduced oxide thickness
    ! Rule: old oxide that remained after exfololiation has the same creep strain
    !=======================================================================

  use parameter_module,   only: moxide, mave, mgrid
  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only: rad_temp, npr_temp, &
      & npr_st, rad_st, rad_st_old, rad_int
  use oxide_data_module,   only: micron2m, no_oxide, no_layer, &
      & no_oxide_ave, cond_value, first_oxide_layer
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,      only: tty_lun, out_lun, inp_lun, aux_lun
  use creep_data_module, only: eps_oxide_grid, creep_strain_eq, &
      & creep_strain_old, creep_strain_new, stress_eq_conv, &
      & i_cr_start, i_cr_end, creep_strain_ave, stress_eq_ave, &
      & if_creep_ox_scale, creep_strain_old_height_shutd, &
      & creep_strain_eq_height_shutd, stress_eq_conv_height_shutd, &
      & creep_strain_eq_before_exfol, creep_strain_eq_ave, &
      & creep_strain_eq_after_exfol, stress_eq_ave_before_exfol, &
      & stress_eq_ave_after_exfol, stress_vm_mean_before_exfol, &
      & stress_vm_mean_after_exfol, creep_strain_before_exfol, &
      & creep_strain_after_exfol, stress_eq, &
      & if_creep_ox_scale, if_creep_metal
  use property_module,    only: LINEAR_PROPERTY_SIMPLE, &
      & OPERAND_ARRAY
  use blockage_data_module,     only: dthick_ox, &
     & thick_ox_skin, rad_ox_skin_large, rad_ox_skin_small, &
     & rad_int_block, thickness_ave_ox_growth_init, &
     & if_exfoliation_oxide_growth, &
     & start_exfoliation_layers, end_exfoliation_layers, &
     & thickness_new_ox_growth_init
  use oxide_stress_data_module, only: ispinel, imagnetite

  implicit none

  ! Argument List
  integer, intent(IN)  :: N, iheight, ioutage
  real, intent(IN)     :: time_now

  ! Local Variables
  integer  :: i, j, io, jo, loc1
  integer, save  :: ip = 0
  real     :: rad_diff, displ_max, slope, summ1, rad_now
  logical  :: reset_creep_strain_old, metal_oxide_interface, &
      & oxide_steam_interface
  real, dimension(moxide)   :: rad_int_after_e  ! radius at the interface
  real, dimension(moxide, mgrid) :: rad_st_after_e, var
  real, dimension(N)  :: stress_vm_mean, temp_cr_mean
  character(LEN = 80), dimension(mave), save    :: operation_type

  operation_type(3) = 'ave_sum'

  ! N = total number of layers; minimum 2 (substrate + average scale)

  ! rad = radius at which temp, temp_change, strain, displ, stress
  ! will be computed

  ! hopefully ioutage is defined when exfoliation is not going on
  stress_eq_ave_before_exfol(1, iheight, ioutage) = &
        & stress_eq_ave(1, 3)
  stress_eq_ave_before_exfol(2, iheight, ioutage) = &
        & stress_eq_ave(imagnetite, 3)

  ! state of stress in the metal does not change after exfoiation
  stress_eq_ave_after_exfol(1, iheight, ioutage) = &
        & stress_eq_ave(1, 3)

  ! a partial re-meshing of the outer layer that exfoliates
  rad_int_after_e = rad_int
  rad_int_after_e(N+1) = rad_int(N) - micron2m * &
     & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite)

  rad_st_after_e(N, 1:moxide) = rad_st(N, 1:moxide)
  do j = 1, npr_st(N)
    rad_st_after_e(N, j) = rad_int_after_e(N) - &
       & (rad_int_after_e(N) - rad_int_after_e(N+1)) * &
       & (j - 1) / (npr_st(N) - 1)
  end do

  ! stresses are set irrespective of creep or not
  var = stress_eq
  do j = 1, npr_st(N)

    call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
    & npr_st(N), 1, mgrid, var(N, :), &
    & rad_st(N, :), stress_eq(N, j))

  end do
  call OPERAND_ARRAY(1, mgrid, stress_eq(N, :),  &
      & 1, npr_st(N), operation_type(3), &
      & stress_eq_ave(N, 3), loc1)
  stress_eq_ave_after_exfol(2, iheight, ioutage) = &
      & stress_eq_ave(imagnetite, 3)

   if (i_cr_start > i_cr_end)  then
    ! there is no creep
    RETURN !  to store info on stresses
  endif

  if (.not. if_creep_ox_scale)  then
    ! there is no creep in the oxide scale
    creep_strain_before_exfol(2, iheight, ioutage)%rad = 0.0
    creep_strain_before_exfol(2, iheight, ioutage)%hoop = 0.0
    creep_strain_before_exfol(2, iheight, ioutage)%ax = 0.0
    creep_strain_eq_before_exfol(2, iheight, ioutage) = 0.0

    creep_strain_after_exfol(2, iheight, ioutage)%rad = 0.0
    creep_strain_after_exfol(2, iheight, ioutage)%hoop = 0.0
    creep_strain_after_exfol(2, iheight, ioutage)%ax = 0.0
    creep_strain_eq_after_exfol(2, iheight, ioutage) = 0.0

    RETURN

  endif

  if (.not. if_exfoliation_oxide_growth)  then

    ! RETURN to store info even when there is no exfoliation

  else if (if_exfoliation_oxide_growth)  then

    if (start_exfoliation_layers /= imagnetite)  then

      write(tty_lun, *) 'ERROR: creep start_exfoliation_layers /= imagnetite', &
        &  ' exfoliation of only one layer implemented '
      stop

    endif

  endif

  ! no need for this one
  ! if (i_cr_end /= imagnetite)  RETURN  ! assume i_cr_end == imagnetite == N

  ! extend this later on to start_exfoliation_layers, end_exfoliation_layers
  ! not only to imagnetite

  ! ip - print id; print only for the first 6 calls of this routine
  ip = ip + 1

  ! store average quantities before exfoliation
  if (if_creep_metal)  then

    creep_strain_before_exfol(1, iheight, ioutage) = creep_strain_ave(1, 3)
    creep_strain_eq_before_exfol(1, iheight, ioutage) = &
        & creep_strain_eq_ave(1, 3)

  else

    creep_strain_before_exfol(1, iheight, ioutage)%rad = 0.0
    creep_strain_before_exfol(1, iheight, ioutage)%hoop = 0.0
    creep_strain_before_exfol(1, iheight, ioutage)%ax = 0.0
    creep_strain_eq_before_exfol(1, iheight, ioutage) = 0.0
    
  endif

  if (if_creep_ox_scale)  then

    creep_strain_before_exfol(2, iheight, ioutage) = &
        & creep_strain_ave(imagnetite, 3)
    creep_strain_eq_before_exfol(2, iheight, ioutage) = &
        & creep_strain_eq_ave(imagnetite, 3)

  endif

  ! check only the displacement between the previous and current position
  ! for the radii only in the exfoliating layers

  displ_max = ABS(rad_st(first_oxide_layer, 1) - &
        & rad_st(first_oxide_layer, npr_st(first_oxide_layer)))
  rad_diff = SUM((rad_st(N, 1:npr_st(N)) - rad_st_after_e(N, 1:npr_st(N)))**2) 

  if (first_oxide_layer < N)  then
    rad_diff = SQRT(rad_diff / npr_st(N))
  else
    rad_diff = 0.0
  endif

  if (rad_diff < eps_oxide_grid * displ_max)  then

    reset_creep_strain_old = .false.

  else

    reset_creep_strain_old = .true.

  endif

  if (.not. reset_creep_strain_old)  RETURN

  i = N

  do j = 1, npr_st(N)

    call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
      & npr_st(N), 1, mgrid, creep_strain_old(N, :)%rad, &
      & rad_st(N, :), creep_strain_new(N, j)%rad)

    call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
      & npr_st(N), 1, mgrid, creep_strain_old(N, :)%hoop, &
      & rad_st(N, :), creep_strain_new(N, j)%hoop)

    call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
      & npr_st(N), 1, mgrid, creep_strain_old(N, :)%ax, &
      & rad_st(N, :), creep_strain_new(N, j)%ax)

      ! write(aux_lun, 4) i, time_now, rad_st(N, j), &
      !    & rad_st_after_e(No, npr_st(No)), &
      !    & rad_now - rad_st_after_e(No, npr_st(No))
 4    format('steam_check_reset_cr ', i2, 1x, 20(1pe13.6, 1x))

  end do

  var = creep_strain_eq

  do j = 1, npr_st(N)

     call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
      & npr_st(N), 1, mgrid, var(N, :), &
      & rad_st(N, :), creep_strain_eq(N, j))

  end do

  ! creep_strain_new was used as dummy here
  creep_strain_old(N, :) = creep_strain_new(N, :)

  ! reset variables that must be stored in magnetite
  creep_strain_old_height_shutd(iheight, ioutage, N, :)%rad = &
        & creep_strain_old(N, :)%rad
  creep_strain_old_height_shutd(iheight, ioutage, N, :)%hoop = &
        & creep_strain_old(N, :)%hoop
  creep_strain_old_height_shutd(iheight, ioutage, N, :)%ax = &
        & creep_strain_old(N, :)%ax
  creep_strain_eq_height_shutd(iheight, ioutage, N, :) = &
        & creep_strain_eq(N, :)

  ! this is only to store "the initial stress" for convergence purposes
  ! when the new interval outage starts
  do i = 1, 10

    var(:, :) = stress_eq_conv(i, :, :)

    do j = 1, npr_st(N)

      call LINEAR_PROPERTY_SIMPLE (rad_st_after_e(N, j), &
      & npr_st(N), 1, mgrid, var(N, :), &
      & rad_st(N, :), stress_eq_conv(i, N, j))

      stress_eq_conv_height_shutd(iheight, ioutage, i, N, j) = &
        & stress_eq_conv(i, N, j)

    end do

  end do

  ! stresses need to be reset just for the output

  ! store average quantities after exfoliation
  if (if_creep_metal)  then
    ! metal does not change after exfoliation
    creep_strain_after_exfol(1, iheight, ioutage) = creep_strain_ave(1, 3)
    creep_strain_eq_after_exfol(1, iheight, ioutage) = &
        & creep_strain_eq_ave(1, 3)
  else

    creep_strain_after_exfol(1, iheight, ioutage)%rad = 0.0
    creep_strain_after_exfol(1, iheight, ioutage)%hoop = 0.0
    creep_strain_after_exfol(1, iheight, ioutage)%ax = 0.0
    creep_strain_eq_after_exfol(1, iheight, ioutage) = 0.0
    
  endif

  if (if_creep_ox_scale)  then

    ! get average creep in each layer
    call OPERAND_ARRAY(1, mgrid, creep_strain_old(N, :)%rad,  &
        & 1, npr_st(N),  operation_type(3), &
        & creep_strain_ave(N, 3)%rad, loc1)
    call OPERAND_ARRAY(1, mgrid, creep_strain_old(N, :)%hoop,  &
        & 1, npr_st(N),  operation_type(3), &
        & creep_strain_ave(N, 3)%hoop, loc1)
    call OPERAND_ARRAY(1, mgrid, creep_strain_old(N, :)%ax,  &
        & 1, npr_st(N),  operation_type(3), &
        & creep_strain_ave(N, 3)%ax, loc1)
    call OPERAND_ARRAY(1, mgrid, creep_strain_eq(N, :),  &
        & 1, npr_st(N),  operation_type(3), &
        & creep_strain_eq_ave(N, 3), loc1)

    creep_strain_after_exfol(2, iheight, ioutage) = &
        & creep_strain_ave(imagnetite, 3)
    creep_strain_eq_after_exfol(2, iheight, ioutage) = &
        & creep_strain_eq_ave(imagnetite, 3)

  endif

  write(out_lun, 1) ioutage, iheight, &
     & 1.0e+6*(rad_st_after_e(N, npr_st(N)) - rad_st(N, npr_st(N))), &
     & creep_strain_before_exfol(2, iheight, ioutage)%hoop, &
     & creep_strain_after_exfol(2, iheight, ioutage)%hoop, &
     & stress_eq_ave_before_exfol(2, iheight, ioutage), &
     & stress_eq_ave_after_exfol(2, iheight, ioutage)
1 format('creep_reset ', i2, 1x, i2, 1x, 100(1pe13.6, 1x))

  ! stop

  RETURN

  END SUBROUTINE CREEP_RESET_EXFOLIATION

  SUBROUTINE STRESS_STRAIN_CREEP_OUT_EXFOL()

    !=======================================================================
    ! Purpose(s):
    !
    ! output the creep rate and von Mises stress before/after exfoliation
    !=======================================================================

  use parameter_module,   only: moxide, mave, mgrid
  use boiler_data_module, only: tube_thickness, tube_outer_radius
  use solution_data_module, only: time_event, no_total_outage, &
      & no_f2l_event_out
  use oxide_data_module,   only: micron2m, no_oxide, no_layer, &
      & no_oxide_ave, cond_value, first_oxide_layer
  use solver_data_module, only: no_temp_rad_points_tube, &
      & no_temp_rad_points_oxide, no_st_rad_points_tube, &
      & no_st_rad_points_oxide 
  use output_module,      only: tty_lun, out_lun, aux_lun, enrg_lun
  use creep_data_module, only: i_cr_start, i_cr_end, &
      & if_creep_ox_scale, creep_strain_old_height_shutd, &
      & creep_strain_eq_height_shutd, stress_eq_conv_height_shutd, &
      & creep_strain_eq_before_exfol, &
      & creep_strain_eq_after_exfol, stress_eq_ave_before_exfol, &
      & stress_eq_ave_after_exfol, stress_vm_mean_before_exfol, &
      & stress_vm_mean_after_exfol, creep_strain_before_exfol, &
      & creep_strain_after_exfol, if_creep_metal
  use blockage_data_module,     only: dthick_ox, &
     & thickness_ave_ox_growth_init, &
     & if_exfoliation_oxide_growth, &
     & start_exfoliation_layers, end_exfoliation_layers, &
     & thickness_new_ox_growth_init
  use oxide_stress_data_module, only: ispinel, imagnetite
  use waterwall_data_module, only:  no_loops, no_height_oxide_output, &
      & length_tube_oxide_output

  implicit none

  ! Argument List

  ! Local Variables
  integer  :: i, j, iheight, ioutage
  integer, save  :: ip = 0
  character(LEN = 1), dimension(0:9) :: ch_index = (/'0','1','2','3','4','5', &
                  & '6','7','8','9'/)

  write(enrg_lun, 93) 'vm_creep_outage height ', (' vm_met'//ch_index(j), &
    & ' vm_mag_b'//ch_index(j), ' vm_mag_a'//ch_index(j), &
    & ' scr_hoop_met'//ch_index(j), &
    & ' scr_hoop_mag_b'//ch_index(j), ' scr_hoop_mag_a'//ch_index(j), &
    & j = 1, MIN(9, no_f2l_event_out)), &
    & (' vm_met'//'1'//ch_index(j-10), &
    & ' vm_mag_b'//'1'//ch_index(j-10), &
    & ' vm_mag_a'//'1'//ch_index(j-10), &
    & ' scr_hoop_met'//'1'//ch_index(j-10), &
    & ' scr_hoop_mag_b'//'1'//ch_index(j-10), &
    & ' scr_hoop_mag_a'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)

  write(out_lun, 93)  'vm_creep_outage height ', (' vm_met'//ch_index(j), &
    & ' vm_mag_b'//ch_index(j), ' vm_mag_a'//ch_index(j), &
    & ' scr_hoop_met'//ch_index(j), &
    & ' scr_hoop_mag_b'//ch_index(j), ' scr_hoop_mag_a'//ch_index(j), &
    & j = 1, MIN(9, no_f2l_event_out)), &
    & (' vm_met'//'1'//ch_index(j-10), &
    & ' vm_mag_b'//'1'//ch_index(j-10), &
    & ' vm_mag_a'//'1'//ch_index(j-10), &
    & ' scr_hoop_met'//'1'//ch_index(j-10), &
    & ' scr_hoop_mag_b'//'1'//ch_index(j-10), &
    & ' scr_hoop_mag_a'//'1'//ch_index(j-10), &
    & j = 1 + MIN(9, no_f2l_event_out), no_f2l_event_out)
 93 format(150(a))

   do j = 1, no_height_oxide_output
     write(enrg_lun, 94) length_tube_oxide_output(j), &
       & (stress_eq_ave_before_exfol(1, j, i), &
       & stress_eq_ave_before_exfol(2, j, i), &
       & stress_eq_ave_after_exfol(1, j, i), &
       & creep_strain_before_exfol(1, j, i)%hoop, &
       & creep_strain_before_exfol(2, j, i)%hoop, &
       & creep_strain_after_exfol(2, j, i)%hoop, &
       & i = 1, no_f2l_event_out)
     write(out_lun, 94) length_tube_oxide_output(j), &
       & (stress_eq_ave_before_exfol(1, j, i), &
       & stress_eq_ave_before_exfol(2, j, i), &
       & stress_eq_ave_after_exfol(1, j, i), &
       & creep_strain_before_exfol(1, j, i)%hoop, &
       & creep_strain_before_exfol(2, j, i)%hoop, &
       & creep_strain_after_exfol(2, j, i)%hoop, &
       & i = 1, no_f2l_event_out)
   end do
 94    format('vm_creep_outage ', 100(1pe13.6, 1x))

  RETURN

  END SUBROUTINE STRESS_STRAIN_CREEP_OUT_EXFOL

END MODULE CREEP_EXFOLIATION_MODULE
