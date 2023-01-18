 MODULE EXFOLIATION_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Input all the data for the blockage - height dependence
  !
  !! Author: Adrian S. Sabau, sabaua@ornl.gov
  !
  !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: EXFOLIATED_VOLUME, EXFOLIATED_SKIN, THICKNESS_INIT_EXFOL, &
      & EXFOLIATED_AREA_ARMITT

 CONTAINS

 SUBROUTINE EXFOLIATED_VOLUME(N, iheight, ioutage_old, ioutage)

  use parameter_module,         only: mstep, zero, pi
  use output_module,            only: tty_lun, out_lun, inp_lun, enrg_lun
  use blockage_data_module,     only: area_fraction_spalled, &
      & energy_elastic_oxide, no_area_fr_spalled, no_deposit_sites, mblock
  use blockage_data_module,     only: max_energy_elastic_event_block, &
      & min_strain_elast_event_block, max_strain_elast_event_block, &
      & id_height_exfoliate2deposit, rad_int_block, oxide_thickness_block, &
      & volume_exfoliated_bend, if_exfoliation_volume, &
      & area_exfoliated_fr_block, deposit_fluid_density, &
      & deposit_oxide_density, deposit_porosity_fraction
  use blockage_data_module,     only: volume_deposit_tube, &
      & radius_deposit_tube, if_exfoliation_oxide_growth
  use solution_data_module,     only: rad_int_event, time_event, &
      & no_f2l_event_out
  use waterwall_data_module, only: no_height_oxide_output, &
      & length_tube_oxide_output
  use property_module, only:  LINEAR_PROPERTY_SIMPLE
  use oxide_data_module,  only: first_oxide_layer, poisson_ratio, &
      & no_oxide_ave, rho_mat, no_oxide, material_name
  use oxide_stress_data_module,   only: ispinel, imagnetite
  use waterwall_data_module, only:  no_loops, loop_transition_length

  implicit none

  ! Argument List
  integer, intent(IN) :: N, iheight, ioutage_old, ioutage

  ! Local variables
  integer :: i, j
  real    :: delta_height, oxide_thickness1, segment_exfoliated_fr, &
          & volume_exfoliated_now, area_exfoliated_fr, &
          & energy_check, area_exfoliated_check, oxide_thickness2, &
          & volume_ref_one_cyl
  real, dimension(2) :: volume_all1, volume_mag1, mass_dry_all, &
          & mass_wet_all, mass_dry_mag, mass_wet_mag
  real, dimension(N) :: rad_large, rad_small, area_exfoliated_fr_layer
  logical, save :: if_first_exfol_volume = .true.
  character(LEN = 80), dimension(20) :: ch_repeat, ch_repeat1
  character(LEN = 1), dimension(0:9)   :: ch_index = (/'0','1','2','3','4','5', &
                  & '6','7','8','9'/)

  if (.not. if_exfoliation_volume)  RETURN

  if (if_first_exfol_volume)  then

    volume_exfoliated_bend = 0.0
    volume_deposit_tube = 0.0

    if_first_exfol_volume = .false.

    if (ABS(deposit_oxide_density(1)) < 0.01)  then

      do i = 1, no_oxide

        if (TRIM(material_name(i)) == 'Fe3O4')  then
          deposit_oxide_density(1) = rho_mat(i)
          write(out_lun, *) 'deposit_oxide_density = ',  rho_mat(i), &
            & 'density_fluid_add ', deposit_fluid_density(1) * &
            & deposit_porosity_fraction(1) / &
            & (deposit_oxide_density(1) * (1.0 - deposit_porosity_fraction(1)))
          EXIT
        endif

      end do

      if (ABS(deposit_oxide_density(1)) < 0.01)  then
        write(out_lun, *) 'ERROR: deposit_oxide_density = ', deposit_oxide_density(1)
        STOP
      endif

    endif

  endif

  ! vertical and horizontal orientation do no work well for identifying 
  ! where the exfoliated mass will go since the tubes can also be h not only v
  ! very simple now, use also 2 coordinates to specify the tube not only height
  ! simple includes also horizontal positions
  if (iheight <= no_height_oxide_output / 2)  then
    id_height_exfoliate2deposit(iheight) = 1 ! or 2
  else 
    id_height_exfoliate2deposit(iheight) = 2
  endif

  ! get the height
  if (iheight == 1)  then
    delta_height = (length_tube_oxide_output(iheight+1) - &
       & length_tube_oxide_output(iheight)) / 2.0
  else if (iheight == no_height_oxide_output)  then
    delta_height = (length_tube_oxide_output(iheight) - &
       & length_tube_oxide_output(iheight-1)) / 2.0
  else
    delta_height = (length_tube_oxide_output(iheight+1) - &
       & length_tube_oxide_output(iheight)) / 2.0 + &
       & (length_tube_oxide_output(iheight) - &
       & length_tube_oxide_output(iheight-1)) / 2.0
  endif

  ! get radii at steam-oxide interface, use rad_int_event(height, outage)
  ! exfoliated thickness = rad_outer - rad_inner
  ! entire oxide scale
  rad_large(1) = rad_int_block(iheight, ioutage, 2)
  rad_small(1) = rad_int_block(iheight, ioutage, N+1)
  do j = 2, N
    rad_large(j) = rad_int_block(iheight, ioutage, j)
    rad_small(j) = rad_int_block(iheight, ioutage, j+1)
  end do  

  if (iheight == 1 .or. iheight == 2)  then
    do j = 1, N
      write(out_lun, 2) iheight, ioutage, &
       & rad_large(j) - rad_small(j), rad_large(j)**2 - rad_small(j)**2
  2   format('rad_large_small ', 2(i2, 1x), 30(1pe13.6, 1x))
    end do
  endif

  ! get the oxide thickness at this outage, 
  ! use oxide_thickness_event(height, outage)
  oxide_thickness1 = rad_int_block(iheight, ioutage, 2) - &
       & rad_int_block(iheight, ioutage, N+1)
  oxide_thickness2 = SUM(oxide_thickness_block(iheight, ioutage, &
         & first_oxide_layer:N))

  ! j = 1 indicates the total energy; as j=1 would indicate metal 
  ! j = 1 the entire scale can exfoliate (all layers) 
  ! j = 3 only magnetite exfoliation (the top first layer)
  ! (j=2 is spinel; but spinel cannot exfoliate without magnetite)
  LAYER1: do j = 1, N

    ! obtain area_exfoliated_fr_block(iheight, ioutage, j_layer)
    call EXFOLIATED_AREA_ARMITT(N, iheight, ioutage, j)

    area_exfoliated_fr = area_exfoliated_fr_block(iheight, ioutage, j)

    area_exfoliated_fr_layer(j) = area_exfoliated_fr  ! for output only

    ! equivalent circumferential length
    segment_exfoliated_fr = 2.0 * pi * area_exfoliated_fr

    if (if_exfoliation_oxide_growth)  then

      if (j == imagnetite)  then

        call EXFOLIATED_SKIN(N, iheight, ioutage_old, ioutage, &
          & area_exfoliated_fr, delta_height, &
          & volume_exfoliated_now)

      else

        ! (rad_large(j)**2 - rad_small(j)**2) - gives a bug
        ! oxide_thickness_block
        ! volumetric fraction of exofilated oxide at this height
        volume_exfoliated_now = 0.5 * delta_height * segment_exfoliated_fr * &
         & (rad_large(j)**2 - rad_small(j)**2)

      endif

    else if (.not. if_exfoliation_oxide_growth)  then

      ! (rad_large(j)**2 - rad_small(j)**2) - gives a bug
      ! oxide_thickness_block
      ! volumetric fraction of exofilated oxide at this height
      volume_exfoliated_now = 0.5 * delta_height * segment_exfoliated_fr * &
         & (rad_large(j)**2 - rad_small(j)**2)

    endif

    ! bend = short for deposit site
    ! no need for each deposit site since at each location the exoliated mass will       
    ! accumulate at only one bend; bends are not corrected for no_layers = 2
    volume_exfoliated_bend(id_height_exfoliate2deposit(iheight), ioutage, j) = &
      & volume_exfoliated_bend(id_height_exfoliate2deposit(iheight), ioutage, j) + &
      & volume_exfoliated_now 

    ! this is actually solid volume not accounting for porosity in the deposit
 
    if (no_loops == 1)  then

      volume_deposit_tube(ioutage, j, 1, 1) = volume_exfoliated_now + &
        & volume_deposit_tube(ioutage, j, 1, 1)

    else if (no_loops == 2)  then

      if (length_tube_oxide_output(iheight) <= loop_transition_length)  then

        ! first loop volume_deposit_tube(ioutage, j, 1, 6)
        volume_deposit_tube(ioutage, j, 1, 6) = volume_exfoliated_now + &
          & volume_deposit_tube(ioutage, j, 1, 6)

      else if (length_tube_oxide_output(iheight) > loop_transition_length)  then

        ! second loop volume_deposit_tube(ioutage, j, 1, 7)
        volume_deposit_tube(ioutage, j, 1, 7) = volume_exfoliated_now + &
          & volume_deposit_tube(ioutage, j, 1, 7)

      endif

    endif

    write(out_lun, 5) iheight, ioutage, max_energy_elastic_event_block(iheight, &
      & ioutage, j), area_exfoliated_fr, segment_exfoliated_fr, &
      & delta_height, rad_large(j)**2 - rad_small(j)**2, volume_exfoliated_now, &
      & rad_int_block(iheight, ioutage, N+1), &
      & rad_int_block(iheight, ioutage, 2)-rad_int_block(iheight, ioutage, N+1), &
      & oxide_thickness1, oxide_thickness2, volume_deposit_tube(ioutage, j, 1, 6), &
      & volume_deposit_tube(ioutage, j, 1, 7), volume_deposit_tube(ioutage, j, 1, 1)
 5  format('energy_area_dh_r2m_vol_exfol ', i3, 1x, i2, 1x, 20(1pe13.6, 1x))

  end do LAYER1

  ! reference volume of a cylinder of length = radius
  volume_ref_one_cyl = pi * radius_deposit_tube**3

  if (no_loops == 1)  then

    volume_all1 = volume_deposit_tube(ioutage, 1, 1, 1)
    volume_mag1 = volume_deposit_tube(ioutage, 3, 1, 1)

  else if (no_loops == 2)  then

    if (length_tube_oxide_output(iheight) <= loop_transition_length)  then

      ! first loop volume_deposit_tube(ioutage, j, 1, 6)
      volume_all1(1) = volume_deposit_tube(ioutage, 1, 1, 6)
      volume_mag1(1) = volume_deposit_tube(ioutage, 3, 1, 6)

    else if (length_tube_oxide_output(iheight) > loop_transition_length)  then

      ! second loop volume_deposit_tube(ioutage, j, 1, 7)
      volume_all1(2) = volume_deposit_tube(ioutage, 1, 1, 7)
      volume_mag1(2) = volume_deposit_tube(ioutage, 3, 1, 7)

    endif

  end if

  mass_dry_all =  volume_all1 * deposit_oxide_density(1) * 1.0e+3
  mass_dry_mag =  volume_mag1 * deposit_oxide_density(1) * 1.0e+3

  write(out_lun, 4) iheight, ioutage, length_tube_oxide_output(iheight), &
      & mass_dry_all(1), mass_dry_mag(1), mass_dry_all(2), mass_dry_mag(2), &
      & (volume_all1(i),  volume_mag1(i), i= 1, 2), &
      & volume_deposit_tube(ioutage, 1, 1, 1) / volume_ref_one_cyl, &
      & volume_deposit_tube(ioutage, 3, 1, 1) / volume_ref_one_cyl
  write(enrg_lun, 4) iheight, ioutage, length_tube_oxide_output(iheight), &
      & mass_dry_all(1), mass_dry_mag(1), mass_dry_all(2), mass_dry_mag(2), &
      & (volume_all1(i),  volume_mag1(i), i= 1, 2), &
      & volume_deposit_tube(ioutage, 1, 1, 1) / volume_ref_one_cyl, &
      & volume_deposit_tube(ioutage, 3, 1, 1) / volume_ref_one_cyl
  4 format('mass_vol_deposit_partial_h_sd ', i3, 1x, i1, 1x, &
      & 20(1pe13.6, 1x))

  write(out_lun, 3) iheight, ioutage, max_energy_elastic_event_block(iheight, &
      & ioutage, 1), area_exfoliated_fr_layer(1), &
      & max_energy_elastic_event_block(iheight, &
      & ioutage, 3), area_exfoliated_fr_layer(3)
  write(enrg_lun, 3) iheight, ioutage, max_energy_elastic_event_block(iheight, &
      & ioutage, 1), area_exfoliated_fr_layer(1), &
      & max_energy_elastic_event_block(iheight, &
      & ioutage, 3), area_exfoliated_fr_layer(3)
  3 format('max_energy_el_armitt_area_exfoliated_fr ', i3, 1x, i1, 1x, &
      & 20(1pe13.6, 1x))

   if (ioutage == no_f2l_event_out .and. &
      & iheight == no_height_oxide_output)  then

    if (no_loops == 1)  then

      write(enrg_lun, 7) 'v_all1 v_mag1'

    else if (no_loops == 2)  then

      write(enrg_lun, 7)  'm_all_dry2 m_mag_dry2 m_all_wet2 m_mag_wet2 ', &
      &  'v_all1 v_mag1 v_all2 v_mag2'
    endif

  7 format('mass_vol_deposit_total_shd outage t[h] t[d] t[w] ', &
      &  'm_all_dry1 m_mag_dry1 m_all_wet1 m_mag_wet1 ', 2(a)) !  ' vfr_all vfr_mag')

     ! final exfoliated volume, mass, 
     do j = 1, no_f2l_event_out

       if (no_loops == 1)  then

         volume_all1 = volume_deposit_tube(j, 1, 1, 1)
         volume_mag1 = volume_deposit_tube(j, 3, 1, 1)

       else if (no_loops == 2)  then

         ! first loop volume_deposit_tube(j, j, 1, 6)
         volume_all1(1) = volume_deposit_tube(j, 1, 1, 6)
         volume_mag1(1) = volume_deposit_tube(j, 3, 1, 6)

         ! second loop volume_deposit_tube(j, j, 1, 7)
         volume_all1(2) = volume_deposit_tube(j, 1, 1, 7)
         volume_mag1(2) = volume_deposit_tube(j, 3, 1, 7)

       end if
       
       ! mass = vol * rho * 1000 - conversion to grams
       mass_dry_all =  volume_all1 * deposit_oxide_density(1) * 1.0e+3
       mass_wet_all =  mass_dry_all * (1.0 + &
         & deposit_fluid_density(1) * deposit_porosity_fraction(1) / &
         & (deposit_oxide_density(1) * (1.0 - deposit_porosity_fraction(1))))
       mass_dry_mag =  volume_mag1 * deposit_oxide_density(1) * 1.0e+3
       mass_wet_mag =  mass_dry_mag * (1.0 + &
         & deposit_fluid_density(1) * deposit_porosity_fraction(1) / &
         & (deposit_oxide_density(1) * (1.0 - deposit_porosity_fraction(1))))
       write(enrg_lun, 6) j, time_event(j), time_event(j) / 24.0, &
         & time_event(j) / (24.0 * 7.0), (mass_dry_all(i), mass_dry_mag(i), &
         & mass_wet_all(i), mass_wet_mag(i), i= 1, no_loops), &
         & (volume_all1(i),  volume_mag1(i), i= 1, no_loops)
         ! & volume_deposit_tube(j, 1, 1, 1) / volume_ref_one_cyl, &
         ! & volume_deposit_tube(j, 3, 1, 1) / volume_ref_one_cyl

     end do

  6 format('mass_vol_deposit_total_shd ', i2, 1x, &
      & 20(1pe13.6, 1x))

    ! ch_repeat(1) = ' en_mag aX_mag'
    ! ch_repeat1(1) = ' en_tot aX_tot'
    ch_repeat(1) = ' en_m aX_m'
    ch_repeat1(1) = ' en_t aX_t'
    ch_repeat(2:10) = ch_repeat(1)
    ch_repeat1(2:10) = ch_repeat1(1)
    write(out_lun, 21) (ch_repeat(j)(1:10)//ch_index(j), &
      &   j = 1, MIN(9, no_f2l_event_out)), &
      & (ch_repeat(j-9)(1:10)//'1'//ch_index(j-10), &
      &   j = MIN(9, no_f2l_event_out) + 1, no_f2l_event_out), &
      & (ch_repeat1(j)(1:10)//ch_index(j), j = 1, MIN(9, no_f2l_event_out)), &
      & (ch_repeat1(j-9)(1:10)//'1'//ch_index(j-10), &
      &   j = MIN(9, no_f2l_event_out) + 1, no_f2l_event_out)
    write(enrg_lun, 21) (ch_repeat(j)(1:10)//ch_index(j), &
      &   j = 1, MIN(9, no_f2l_event_out)), &
      & (ch_repeat(j-9)(1:10)//'1'//ch_index(j-10), &
      &   j = MIN(9, no_f2l_event_out) + 1, no_f2l_event_out), &
      & (ch_repeat1(j)(1:10)//ch_index(j), &
      &   j = 1, MIN(9, no_f2l_event_out)), &
      & (ch_repeat1(j-9)(1:10)//'1'//ch_index(j-10), &
      &   j = MIN(9, no_f2l_event_out) + 1, no_f2l_event_out)

 21  format('#energy_armitt_areaex ih length', 100(a))
    do i = 1, no_height_oxide_output
     write(enrg_lun, 22) i, length_tube_oxide_output(i), &
      & (max_energy_elastic_event_block(i, &
      & j, 3), area_exfoliated_fr_block(i, j, 3), &
      & j = 1, no_f2l_event_out), &
      & (max_energy_elastic_event_block(i, &
      & j, 1), area_exfoliated_fr_block(i, j, 1), &
      & j = 1, no_f2l_event_out)
    end do
 22  format('energy_armitt_areaex ', i2, 1x, 100(1pe13.6, 1x))

   endif 

  RETURN

 END SUBROUTINE EXFOLIATED_VOLUME

 SUBROUTINE EXFOLIATED_AREA_ARMITT(N, iheight, ioutage, j_layer)

  use output_module,            only: tty_lun, out_lun, inp_lun, enrg_lun
  use blockage_data_module,     only: area_fraction_spalled, &
      & energy_elastic_oxide, &
      & no_area_fr_spalled, no_deposit_sites, mblock
  use blockage_data_module,     only: max_energy_elastic_event_block, &
      & min_strain_elast_event_block, max_strain_elast_event_block, &
      & id_height_exfoliate2deposit, rad_int_block, oxide_thickness_block, &
      & if_exfoliation_volume, area_exfoliated_fr_block, &
      & start_energy_exfoliation_layer, end_energy_exfoliation_layer
  use property_module, only:  LINEAR_PROPERTY_SIMPLE
  use oxide_data_module,  only: first_oxide_layer, poisson_ratio, &
      & no_oxide_ave
  use oxide_stress_data_module,   only: ispinel, imagnetite
  use solution_data_module, only      : no_oxide_layers

  implicit none

  ! Argument List
  integer, intent(IN) :: N, iheight, ioutage, j_layer

  ! Local variables
  integer :: i
  real    :: delta_height, oxide_thickness1, segment_exfoliated_fr, &
          & area_exfoliated_fr, elastic_energy4_armitt, &
          & energy_check, area_exfoliated_check, oxide_thickness2, &
          & volume_ref_one_cyl, poisson
  real, dimension(N) :: rad_large, rad_small, area_exfoliated_fr_layer
  logical, save :: if_first_exfol_volume = .true.
  character(LEN = 80), dimension(20) :: ch_repeat, ch_repeat1

  if (.not. if_exfoliation_volume)  then

    area_exfoliated_fr_block(iheight, ioutage, j_layer) = 0.0
    RETURN

  endif

  if (if_first_exfol_volume)  then

    ! testing the energy 
    do i = 1, no_area_fr_spalled -1

       energy_check = 0.5 * (energy_elastic_oxide(i) + energy_elastic_oxide(i+1))
       call LINEAR_PROPERTY_SIMPLE(energy_check, &
         & no_area_fr_spalled-1, 1, mblock, &
         & area_fraction_spalled, energy_elastic_oxide, area_exfoliated_check)
        write(out_lun, *) 'armitt_area_fr_vs_energy_check2 ', &
         & energy_check, area_exfoliated_check
    end do

    if_first_exfol_volume = .false.

  endif

  ! j = 1 indicates the total energy; as j=1 would indicate metal 
  ! j = 1 the entire scale can exfoliate (all layers) 
  ! j = 3 only magnetite exfoliation (the top first layer)
  ! (j=2 is spinel; but spinel cannot exfoliate without magnetite)
  ! do j = 1, N

  ! ARMITT uses the strain energy for the entire oxide scale for 316 steel 
  ! (1) use the strain energy for the entire oxide scale to get exfol_fraction
  ! (2) use the strain energy for layers that exfoliate but change the 
  !     Armitt's exfol fraction (the x axis since one 1/2 of oxide will exfoliate)

    if (j_layer == 1)  then
      poisson = poisson_ratio(no_oxide_ave)
      poisson = SUM(poisson_ratio(first_oxide_layer:no_oxide_layers) * &
           & oxide_thickness_block(iheight, ioutage, &
           & first_oxide_layer:no_oxide_layers)) / &
           & SUM(oxide_thickness_block(iheight, ioutage, &
           & first_oxide_layer:no_oxide_layers))
      elastic_energy4_armitt = max_energy_elastic_event_block(iheight, &
           & ioutage, j_layer)  ! it was already averaged over all layers
    else

      if (start_energy_exfoliation_layer(j_layer) < &
         & end_energy_exfoliation_layer(j_layer))  then

        ! energy would be calculated for all layers; traditional Armitt's fraction
        poisson = poisson_ratio(no_oxide_ave)
        poisson = SUM(poisson_ratio(start_energy_exfoliation_layer(j_layer) : &
           & end_energy_exfoliation_layer(j_layer))* &
           & oxide_thickness_block(iheight, ioutage, &
           & start_energy_exfoliation_layer(j_layer) : &
           & end_energy_exfoliation_layer(j_layer))) / &
           & SUM(oxide_thickness_block(iheight, ioutage, &
           & start_energy_exfoliation_layer(j_layer) : &
           & end_energy_exfoliation_layer(j_layer)))
        elastic_energy4_armitt = SUM(max_energy_elastic_event_block(iheight, ioutage, &
           & start_energy_exfoliation_layer(j_layer) : &
           & end_energy_exfoliation_layer(j_layer)))

      else if (start_energy_exfoliation_layer(j_layer) == &
           & end_energy_exfoliation_layer(j_layer))  then

        ! energy would be calculated for each layer in part
        poisson = poisson_ratio(j_layer)
        elastic_energy4_armitt = max_energy_elastic_event_block(iheight, &
           & ioutage, j_layer)

      endif

    endif

    call LINEAR_PROPERTY_SIMPLE(elastic_energy4_armitt * (1.0-poisson), &
      & no_area_fr_spalled - 1, 1, mblock, &
      & area_fraction_spalled, energy_elastic_oxide, &
      & area_exfoliated_fr_block(iheight, ioutage, j_layer))

    if (iheight == 2 .and. j_layer == imagnetite)  then

       write(out_lun, 1)  ioutage, max_energy_elastic_event_block(iheight, &
      & ioutage, j_layer), elastic_energy4_armitt, &
      & area_exfoliated_fr_block(iheight, ioutage, j_layer), &
      & poisson
 1    format('armitt_check1 ', i2, 1x, 10(1pe13.6, 1x))

    endif

  RETURN

 END SUBROUTINE EXFOLIATED_AREA_ARMITT

 SUBROUTINE THICKNESS_INIT_EXFOL(iheight, id_previous_outage, id_now_outage)

  use parameter_module,           only: zero, pi
  use output_module,              only: tty_lun, out_lun
  use oxide_stress_data_module,   only: ispinel, imagnetite
  use oxide_data_module,  only: first_oxide_layer
  use solution_data_module,       only: no_oxide_layers
  use blockage_data_module,       only: dthick_ox, &
     & thickness_new_ox_growth_init, thickness_new_ox_growth, &
     & if_exfoliation_oxide_growth, growth_mode_after_exfoliation, &
     & start_exfoliation_layers, end_exfoliation_layers

  implicit none

  ! Argument List
  integer, intent(IN) :: iheight, id_now_outage
  integer, intent(INOUT) :: id_previous_outage

  ! Local variables
  integer :: i, j, k

  ! default for layers that do not exfoliate; using ('continuous') option

  if (if_exfoliation_oxide_growth)  then

  ! below only for layers that exfoliate

        ! no need to reset the oxide thickness between outages since 
        ! the one variable is used per outage
        ! dthick_ox(k, id_now_outage, imagnetite) = 0.0

        ! this has to be used if the check is done at all times; 
        ! before calling this routine the occurance of an outage at this time
        ! iteration was identifyed; thus initiate quantities at the new outage
        ! if (id_previous_outage < id_now_outage)  then
          ! set thickness_new_ox_growth_init

          SELECT CASE (TRIM(growth_mode_after_exfoliation))

          CASE('continuous')

            thickness_new_ox_growth_init(iheight, id_now_outage, &
              & start_exfoliation_layers:no_oxide_layers) = &
              & dthick_ox(iheight, id_previous_outage, &
              & start_exfoliation_layers:no_oxide_layers) + &
              & thickness_new_ox_growth_init(iheight, id_previous_outage, &
              & start_exfoliation_layers:no_oxide_layers)

          CASE('zero')

            thickness_new_ox_growth_init(iheight, id_now_outage, &
              & start_exfoliation_layers:no_oxide_layers) = 0.0

          CASE('start_first_exfolation')

            if (id_previous_outage == 0)  then

              thickness_new_ox_growth_init(iheight, id_now_outage, &
                & start_exfoliation_layers:no_oxide_layers) = &
                & dthick_ox(iheight, id_previous_outage, &
                & start_exfoliation_layers:no_oxide_layers)

            else

              thickness_new_ox_growth_init(iheight, id_now_outage, &
                & start_exfoliation_layers:no_oxide_layers) = &
                & thickness_new_ox_growth_init(iheight, 1, &
                & start_exfoliation_layers:no_oxide_layers)

            endif

          CASE('from_given_thickness')

            thickness_new_ox_growth_init(iheight, id_now_outage, &
              & start_exfoliation_layers:no_oxide_layers) = &
              & thickness_new_ox_growth
          
          CASE DEFAULT
            write(6, *) 'ERROR: growth_mode_after_exfoliation not ok'
            STOP
          END SELECT

    endif

    if (iheight == 1)  then
      write(tty_lun, 1) TRIM(growth_mode_after_exfoliation), &
        & id_previous_outage, id_now_outage, &
        & thickness_new_ox_growth_init(iheight, id_now_outage, imagnetite), &
        & thickness_new_ox_growth_init(iheight, id_previous_outage, imagnetite), &
        & dthick_ox(iheight, id_previous_outage, imagnetite)
      write(out_lun, 1) TRIM(growth_mode_after_exfoliation), &
        & id_previous_outage, id_now_outage, &
        & thickness_new_ox_growth_init(iheight, id_now_outage, imagnetite), &
        & thickness_new_ox_growth_init(iheight, id_previous_outage, imagnetite), &
        & dthick_ox(iheight, id_previous_outage, imagnetite)
    endif

 1  format('growth_mode_after_exfoliation_run ', a10, 1x, 2(i2, 1x), &
              & 10(1pe13.6, 1x))

  RETURN

 END SUBROUTINE THICKNESS_INIT_EXFOL

 SUBROUTINE EXFOLIATED_SKIN(N, iheight, ioutage_old, ioutage, &
    & area_exfoliated_fr, delta_height, volume_exfoliated_now)

  use parameter_module,         only: zero, pi
  use oxide_stress_data_module, only: ispinel, imagnetite
  use blockage_data_module,     only: dthick_ox, &
     & thick_ox_skin, rad_ox_skin_large, rad_ox_skin_small, &
     & fraction_skin_on, rad_int_block, fraction_skin_on, mskin, &
     & thickness_ave_ox_growth_init, if_exfoliation_oxide_growth, &
     & start_exfoliation_layers, end_exfoliation_layers, &
     & skin_no_start, skin_no_end, exfoliate_fr_skin_outage, &
     & fraction_skin_remain_outage, dox_before_exfol, &
     & dox_after_exfol, thickness_new_ox_growth_init
  use output_module,            only: tty_lun, out_lun
  use oxide_data_module,        only: micron2m, first_oxide_layer

  implicit none

  ! Argument List
  integer, intent(IN) :: N, iheight, ioutage, ioutage_old
  real, intent(IN)    :: area_exfoliated_fr, delta_height
  real, intent(OUT)   :: volume_exfoliated_now

  ! Local variables
  integer :: i, j, k, skin_no_start_now, skin_no_end_now

  real :: fraction_to_exfoliate, fraction_skin_total, delta_skin_thickness, &
       & dox_spinel
  real, dimension(mskin) :: exfoliate_fr_skin

  if (.not. if_exfoliation_oxide_growth)  then

    volume_exfoliated_now = 0.0
    exfoliate_fr_skin_outage(iheight, ioutage) = 0.0
    fraction_skin_remain_outage(iheight, ioutage) = 1.0

    ! magnetite before exfoliation
    dox_before_exfol(2, iheight, ioutage) = &
       & dthick_ox(iheight, ioutage_old, imagnetite) + &
       & thickness_new_ox_growth_init(iheight, ioutage_old, imagnetite)
    ! spinel before exfoliation
    dox_spinel = dthick_ox(iheight, ioutage_old, first_oxide_layer) + &
       & thickness_new_ox_growth_init(iheight, ioutage_old, first_oxide_layer)
    ! entire oxide thickness before exfoliation
    dox_before_exfol(1, iheight, ioutage) = dox_spinel + &
       & dox_before_exfol(2, iheight, ioutage)

    ! magnetite after exfoliation
    dox_after_exfol(2, iheight, ioutage) = dox_before_exfol(2, iheight, ioutage) 
    ! entire oxide thickness after exfoliation
    dox_after_exfol(1, iheight, ioutage) = dox_before_exfol(1, iheight, ioutage)

    RETURN

  else if (if_exfoliation_oxide_growth)  then

    if (start_exfoliation_layers /= imagnetite)  then

      write(tty_lun, *) 'ERROR: start_exfoliation_layers /= imagnetite'
      stop

    endif

  endif

  ! extend this later on to start_exfoliation_layers, end_exfoliation_layers
  ! not only to imagnetite

   ! initialize before the first exfoliation event
   if (ioutage == 1)  then

     skin_no_start(iheight) = 1  ! the first active skin
     skin_no_end(iheight) = skin_no_start(iheight)

     ! refer only to magnetite

     ! intialize skin thickness in [microns]
     thick_ox_skin(1, iheight, imagnetite) = &
        & (rad_int_block(iheight, ioutage, imagnetite) - &
        & rad_int_block(iheight, ioutage, 4)) / micron2m

     ! initialize the radius at the spinel-magnetite interface
     ! this would be valid all the time as spinel keeps growing
     rad_ox_skin_large(1, iheight, imagnetite) = &
        & rad_int_block(iheight, ioutage, imagnetite)

     ! initialize the radius at the magnetite surface (steam contact)
     rad_ox_skin_small(1, iheight, imagnetite) = &
        & rad_ox_skin_large(1, iheight, imagnetite) - &
        & thick_ox_skin(1, iheight, imagnetite) * micron2m

     ! initialize the skin fraction
     fraction_skin_on(1, iheight, imagnetite) = 1.0

   else if (ioutage > 1)  then

     ! initialize the radius at the spinel-magnetite interface
     ! this would be valid all the time as spinel keeps growing without exfoliatin
     do k = skin_no_start(iheight), skin_no_end(iheight)
       rad_ox_skin_large(k, iheight, imagnetite) = &
          & rad_int_block(iheight, ioutage, imagnetite)
     end do

     ! skins were thickened between the two exfoliation events by the same amount
     do k = skin_no_start(iheight), skin_no_end(iheight)

       thick_ox_skin(k, iheight, imagnetite) = &
          & thick_ox_skin(k, iheight, imagnetite) + &
          & dthick_ox(iheight, ioutage_old, imagnetite)

       rad_ox_skin_small(k, iheight, imagnetite) = &
          & rad_ox_skin_large(k, iheight, imagnetite) - &
          & thick_ox_skin(k, iheight, imagnetite) * micron2m

     end do

   endif

   ! fr to exfoliate right now
   fraction_to_exfoliate = area_exfoliated_fr

   if (iheight == 2)  then
     write(out_lun, 3) ioutage, skin_no_start(iheight), skin_no_end(iheight), &
       & area_exfoliated_fr, &
       & thickness_ave_ox_growth_init(iheight, ioutage_old, imagnetite), &
       & dthick_ox(iheight, ioutage_old, imagnetite), &
       & (thick_ox_skin(k, iheight, imagnetite), k = skin_no_start(iheight), &
       & skin_no_end(iheight)), &
       & (rad_ox_skin_large(k, iheight, imagnetite), k = skin_no_start(iheight), &
       & skin_no_end(iheight))
 3   format('skin_init1 ', 3(i2, 1x), 20(1pe13.6, 1x))
   endif

   if (fraction_to_exfoliate <= 1.0e-5)  then

     volume_exfoliated_now = 0.0

     ! average initial deposit thickness increases if no exfoliation
     thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) = &
         & thickness_ave_ox_growth_init(iheight, ioutage_old, imagnetite) + &
         & dthick_ox(iheight, ioutage_old, imagnetite)

     if (iheight == 2)  then
       write(out_lun, 4) ioutage, skin_no_start(iheight), skin_no_end(iheight), &
         & area_exfoliated_fr, fraction_to_exfoliate, &
         & dthick_ox(iheight, ioutage_old, imagnetite), &
         & thickness_ave_ox_growth_init(iheight, ioutage_old, imagnetite), &
         & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite)
 4     format('none_exfol_skin_init2 ', 3(i2, 1x), 20(1pe13.6, 1x))
     endif

     exfoliate_fr_skin_outage(iheight, ioutage) = 0.0
     fraction_skin_remain_outage(iheight, ioutage) = &
         & fraction_skin_on(1, iheight, imagnetite)

     ! magnetite before exfoliation
     dox_before_exfol(2, iheight, ioutage) = &
       & dthick_ox(iheight, ioutage_old, imagnetite) + &
       & thickness_ave_ox_growth_init(iheight, ioutage_old, imagnetite)
     ! spinel before exfoliation
     dox_spinel = dthick_ox(iheight, ioutage_old, first_oxide_layer) + &
       & thickness_new_ox_growth_init(iheight, ioutage_old, first_oxide_layer)
     ! entire oxide thickness before exfoliation
     dox_before_exfol(1, iheight, ioutage) = dox_spinel + &
       & dox_before_exfol(2, iheight, ioutage)

     ! magnetite after exfoliation
     dox_after_exfol(2, iheight, ioutage) = dox_before_exfol(2, iheight, ioutage) 
     ! entire oxide thickness after exfoliation
     dox_after_exfol(1, iheight, ioutage) = dox_before_exfol(1, iheight, ioutage)

     RETURN

   endif

   ! initialize volume to be exfoliated now
   volume_exfoliated_now = 0.0

   ! this is at the next outage 
   thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) = 0.0

   ! initialize new skins
   skin_no_start_now = skin_no_start(iheight)
   skin_no_end_now = skin_no_end(iheight)

   SKIN_LOOP1: do k = skin_no_start(iheight), skin_no_end(iheight) ! 1, mskin ! skin_no_max

     ! skin was already exfoliated; for older skins
     ! this should be used with do k = 1, mskin when all skins are included
     if (ABS(fraction_skin_on(k, iheight, imagnetite)) <= 1.0e-4)  cycle SKIN_LOOP1

     ! no more matter to exfoliate
     if (fraction_to_exfoliate <= 1.0e-5)  then

       ! the average skin still needs to be incremented
       thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) = &
           & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) + &
           & fraction_skin_on(k, iheight, imagnetite) * &
           & thick_ox_skin(k, iheight, imagnetite)
       delta_skin_thickness = fraction_skin_on(k, iheight, imagnetite) * &
           & thick_ox_skin(k, iheight, imagnetite)
 
     else if (fraction_to_exfoliate > 1.0e-5)  then

       ! if (ABS(fraction_skin_on(k, iheight, imagnetite)) > 1.0e-4)  then

       if (fraction_skin_on(k, iheight, imagnetite) >= fraction_to_exfoliate)  then

         ! all exfoliation matter will be exfoliated in this skin
         exfoliate_fr_skin(k) = fraction_to_exfoliate

       else if (fraction_skin_on(k, iheight, imagnetite) < fraction_to_exfoliate)  then

         ! a fraction will be exfoliated from this skin
         exfoliate_fr_skin(k) = fraction_skin_on(k, iheight, imagnetite)

       end if

       if (exfoliate_fr_skin(k) > 0.0)  then

         volume_exfoliated_now = volume_exfoliated_now + &
           & VOL_EXFOL(exfoliate_fr_skin(k), &
           &           thick_ox_skin(k, iheight, imagnetite) * micron2m, &
           &           rad_ox_skin_large(k, iheight, imagnetite), delta_height)

         ! fraction of skin that stays on
         fraction_skin_on(k, iheight, imagnetite) = &
           & fraction_skin_on(k, iheight, imagnetite) - exfoliate_fr_skin(k)

         ! fraction of skin to be exfoliated from next skins
         fraction_to_exfoliate = fraction_to_exfoliate - exfoliate_fr_skin(k)

         ! get the average initial deposit thickness after exfoliation
         thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) = &
           & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite) + &
           & fraction_skin_on(k, iheight, imagnetite) * &
           & thick_ox_skin(k, iheight, imagnetite)
         delta_skin_thickness = fraction_skin_on(k, iheight, imagnetite) * &
           & thick_ox_skin(k, iheight, imagnetite)

         ! add a new skin
         skin_no_end_now = skin_no_end_now + 1

         ! initialize the fraction of the new skin
         fraction_skin_on(skin_no_end_now, iheight, imagnetite) = &
            & exfoliate_fr_skin(k)

         ! thickness for the new skin = zero
         thick_ox_skin(skin_no_end_now, iheight, imagnetite) = 0.0

         ! magnetite thickness continues to grow for this skin that partially exfoliated

         ! remove an older skin when it vanishes 
         if (fraction_skin_on(k, iheight, imagnetite) <= 1.0e-5)  then
           skin_no_start_now = skin_no_start_now + 1
         endif

       else if (exfoliate_fr_skin(k) <= 0.0)  then

         delta_skin_thickness = 0.0

       endif

     endif

         if (iheight == 2)  then
           write(out_lun, 1) ioutage, skin_no_start(iheight), k,  &
             & skin_no_end(iheight), skin_no_start_now, skin_no_end_now, &
             & exfoliate_fr_skin(k), fraction_to_exfoliate, &
             & fraction_skin_on(k, iheight, imagnetite), &
             & thick_ox_skin(skin_no_end_now, iheight, imagnetite), &
             & fraction_skin_on(skin_no_end_now, iheight, imagnetite), &
             & delta_skin_thickness
 1         format('skin_data1 ', 6(i2, 1x), 20(1pe13.6, 1x))
         endif

   end do SKIN_LOOP1

   if (skin_no_start(iheight) > 1)  then
     exfoliate_fr_skin_outage(iheight, ioutage) = 0.0
   else
     exfoliate_fr_skin_outage(iheight, ioutage) = exfoliate_fr_skin(1)
   endif

   fraction_skin_remain_outage(iheight, ioutage) = &
       & fraction_skin_on(1, iheight, imagnetite)

   ! magnetite before exfoliation
   dox_before_exfol(2, iheight, ioutage) = &
       & dthick_ox(iheight, ioutage_old, imagnetite) + &
       & thickness_ave_ox_growth_init(iheight, ioutage_old, imagnetite)
   ! spinel before exfoliation
   dox_spinel = dthick_ox(iheight, ioutage_old, first_oxide_layer) + &
       & thickness_new_ox_growth_init(iheight, ioutage_old, first_oxide_layer)
   ! entire oxide thickness before exfoliation
   dox_before_exfol(1, iheight, ioutage) = dox_spinel + &
       & dox_before_exfol(2, iheight, ioutage)

   ! magnetite after exfoliation
   dox_after_exfol(2, iheight, ioutage) = &
       & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite)
   ! entire oxide thickness after exfoliation
   dox_after_exfol(1, iheight, ioutage) = dox_spinel + &
       & dox_after_exfol(2, iheight, ioutage)

   ! update the number of skins
   skin_no_start(iheight) = skin_no_start_now
   skin_no_end(iheight) = skin_no_end_now

   fraction_skin_total = 0.0
   do k = skin_no_start(iheight), skin_no_end(iheight)

     ! the skin fractions should be 1 at all times
     fraction_skin_total = fraction_skin_total + &
        & fraction_skin_on(k, iheight, imagnetite)

   end do

   if (iheight == 2)  then
     write(out_lun, 2) ioutage, fraction_to_exfoliate, &
      & fraction_skin_total, volume_exfoliated_now, &
      & thickness_ave_ox_growth_init(iheight, ioutage, imagnetite)
 2  format('skin_out1 ', 1(i2, 1x), 20(1pe13.6, 1x))
   endif

   RETURN

  END SUBROUTINE EXFOLIATED_SKIN

  REAL FUNCTION VOL_EXFOL(exfoliate_fr, thickness, &
         & rad_large, delta_height)

    ! evaluate the volume exfoliated
  use parameter_module,         only:  pi

  real, Intent(IN):: exfoliate_fr, thickness, &
         & rad_large, delta_height

  real  :: segment_exfoliated_fr, rad_small

    ! equivalent circumferential length
    segment_exfoliated_fr = 2.0 * pi * exfoliate_fr

    rad_small = rad_large - thickness

    ! volumetric fraction of exofilated oxide at this height
    VOL_EXFOL = 0.5 * delta_height * segment_exfoliated_fr * &
       & (rad_large**2 - rad_small**2)

  return

  END FUNCTION VOL_EXFOL

END MODULE EXFOLIATION_MODULE
