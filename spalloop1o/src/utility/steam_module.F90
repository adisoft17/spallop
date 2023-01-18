MODULE STEAM_MODULE
    !==========================================================================================
    !  evaluating steam properties
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !==========================================================================================

IMPLICIT NONE
SAVE

PRIVATE

DOUBLE PRECISION :: Mass_Flow_Rate_kg_m2_s_full         ! mass flow rate (kg/m2-sec) at full load
DOUBLE PRECISION :: Mass_Flow_Rate_kg_m2_s_low          ! mass flow rate (kg/m2-sec) at low  load
DOUBLE PRECISION :: Heat_Flux_kW_m2_full                ! heat flux (kW/m2) at full load (evaluated at pipe OD)
DOUBLE PRECISION :: Heat_Flux_kW_m2_low                 ! heat flux (kW/m2) at low  load
DOUBLE PRECISION :: Pipe_Inner_Diameter_m               ! pipe inner diameter (m)

logical          :: if_deb_Calc_HTC_with_Def = .false.  ! T if want to print out
logical          :: if_deb_Calc_HTC_dHTC = .false.      ! T if want to print out
logical          :: if_deb_Calc_HTC_Now = .false.       ! T if want to print out

!!! Other data
DOUBLE PRECISION :: math_PI = 3.141593D0

PUBLIC :: Pipe_Inner_Diameter_m, math_PI, &
    & Mass_Flow_Rate_kg_m2_s_full, Mass_Flow_Rate_kg_m2_s_low, &
    & Heat_Flux_kW_m2_full, Heat_Flux_kW_m2_low

!! Global variables
!!! Data used by NIST REFPROP v8
INTEGER, PARAMETER :: ncmax = 20                ! number of components (1 for pure fluid)
CHARACTER(LEN=255) :: h_fluid_name(ncmax)       ! file names specifying fluid/mixture components
CHARACTER(LEN=255) :: h_mixture_name            ! file name containing coefficients for mixture model
CHARACTER(LEN=3) :: h_reference_state           ! reference state for thermodynamic calculations
CHARACTER(LEN=255) :: h_nist_err                ! error string
INTEGER :: i_nist_err                           ! error flag:  0 = successful
DOUBLE PRECISION :: xbulk_composition(ncmax)    ! overall (bulk) composition  [array of mol frac]
DOUBLE PRECISION :: xliquid_composition(ncmax)  ! composition of liquid phase [array of mol frac]
DOUBLE PRECISION :: xvapor_composition(ncmax)   ! composition of vapor phase [array of mol frac]
DOUBLE PRECISION :: molecular_weight_kg_mol     ! molecular weight / molar mass  of water [kg/mol]


!!! Other data

! Declare all the public functions
PUBLIC :: Calc_dTSteam_dHeight, STEAM_PIPE_EXPON

CONTAINS

!!!=========================================================================!!!
!! Calculate properties at given temperature and pressure using NIST REFPROP
!
DOUBLE PRECISION FUNCTION Calc_dTSteam_dHeight(id_high_or_low, Temperature_Steam_degC, &
    & Heat_Flux_Fire_W_m2)
USE parameter_module,       ONLY: id_low, id_high
USE boiler_data_module,     ONLY: press_inner_pulse, press_inner_idle, tube_outer_radius
USE waterwall_data_module,  ONLY: flow_rate_steam_high, flow_rate_steam_low
USE output_module,          ONLY: out_lun
IMPLICIT NONE
    INTEGER, INTENT(IN)  :: id_high_or_low
    DOUBLE PRECISION, INTENT(IN) :: Temperature_Steam_degC      ! steam temperature [degC]
    DOUBLE PRECISION, INTENT(IN) :: Heat_Flux_Fire_W_m2         ! Heat flux from the fire at the current height [W m-2]

!   Local variables
    DOUBLE PRECISION :: Pressure_MPa            ! pressure [mega-pascal]
    DOUBLE PRECISION :: Mass_Flow_Rate_kg_s     ! mass flow rate [kg s-1]
!   Dummy variables
    DOUBLE PRECISION :: Density_kg_m3              ! bulk density [kg/m^3]
    DOUBLE PRECISION :: Cp_heat_capacity_J_kg_degC ! isobaric (constant p) heat capacity [J/kg-degC]
    DOUBLE PRECISION :: Viscosity_Pa_s             ! viscosity (Pa s)
    DOUBLE PRECISION :: Thermal_conductivity_W_m_degC  ! thermal conductivity (W/m-degC)
    DOUBLE PRECISION :: Enthalpy_J_kg              ! overall (bulk) enthalpy [J/kg]
    DOUBLE PRECISION :: Prandtl_Number             ! Prandtl number

    IF( id_high_or_low == id_high )THEN
        ! full load
        Pressure_MPa = press_inner_pulse(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_high
    ELSE IF( id_high_or_low == id_low)THEN
        ! low load
        Pressure_MPa = press_inner_idle(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_low
    ELSE
        ! error
        Pressure_MPa = press_inner_pulse(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_high
        WRITE(out_lun,'(A)') "wzHTC_Warning: Neither high nor low load"
    ENDIF

    ! obtain Cp from NIST database
    Cp_heat_capacity_J_kg_degC = 3.0e+3  ! used to be 3.0 until 12/12/13 but routine not used

    Calc_dTSteam_dHeight = 2.0D0 * math_PI * tube_outer_radius * Heat_Flux_Fire_W_m2 / &
        & (Mass_Flow_Rate_kg_s * Cp_heat_capacity_J_kg_degC)

RETURN
END FUNCTION Calc_dTSteam_dHeight

DOUBLE PRECISION FUNCTION STEAM_PIPE_EXPON(id_high_or_low, Temperature_Steam_degC, &
    & htc_total)
USE parameter_module,       ONLY: id_low, id_high
USE boiler_data_module,     ONLY: press_inner_pulse, press_inner_idle, tube_outer_radius
USE waterwall_data_module,  ONLY: flow_rate_steam_high, flow_rate_steam_low
USE output_module,          ONLY: out_lun

IMPLICIT NONE

    INTEGER, INTENT(IN)  :: id_high_or_low
    DOUBLE PRECISION, INTENT(IN) :: Temperature_Steam_degC      ! steam temperature [degC]
    DOUBLE PRECISION, INTENT(IN) :: htc_total         ! Total heat transfer coefficient at the current height [W m-2]

   ! Cp_water_24.13MPa[400:450] = 5.0664535e+6, -3.4463366e+4, 7.8338918e+1,-5.9443173e-2
   ! Cp_water_24.13MPa[450:550] = 1.8631205e+5, -9.9879955e+2, 1.83136, -1.1286825e-3     
   ! Cp_water_24.13MPa[550:700] = 2.7739647e+4, -1.0361113e+2, 1.4484534e-1, -6.8453212e-5
   ! exponents of 0, 1, 2, and 3
   

!   Local variables
    DOUBLE PRECISION :: Pressure_MPa            ! pressure [mega-pascal]
    DOUBLE PRECISION :: Mass_Flow_Rate_kg_s     ! mass flow rate [kg s-1]
!   Dummy variables
    DOUBLE PRECISION :: Density_kg_m3              ! bulk density [kg/m^3]
    DOUBLE PRECISION :: Cp_heat_capacity_J_kg_degC ! isobaric (constant p) heat capacity [J/kg-degC]
    DOUBLE PRECISION :: Viscosity_Pa_s             ! viscosity (Pa s)
    DOUBLE PRECISION :: Thermal_conductivity_W_m_degC  ! thermal conductivity (W/m-degC)
    DOUBLE PRECISION :: Enthalpy_J_kg              ! overall (bulk) enthalpy [J/kg]
    DOUBLE PRECISION :: Prandtl_Number             ! Prandtl number
    DOUBLE PRECISION, dimension(4) :: Cp_water1_24MPa = (/5.0664535e+6, &
        & -3.4463366e+4, 7.8338918e+1,-5.9443173e-2/), &
        & Cp_water2_24MPa = (/1.8631205e+5, -9.9879955e+2, 1.83136, -1.1286825e-3/), &
        & Cp_water3_24MPa = (/2.7739647e+4, -1.0361113e+2, 1.4484534e-1, -6.8453212e-5/)
   ! Cp_water1_19MPa for Czeck boiler [480:690C]
   DOUBLE PRECISION, dimension(4) :: Cp_water1_19MPa = (/2.751125e+4, &
        & -1.106725e+2, 1.659983e-1,-8.393314e-5/)
    integer   :: i

    IF( id_high_or_low == id_high )THEN
        ! full load
        Pressure_MPa = press_inner_pulse(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_high
    ELSE IF( id_high_or_low == id_low)THEN
        ! low load
        Pressure_MPa = press_inner_idle(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_low
    ELSE
        ! error
        Pressure_MPa = press_inner_pulse(1)
        Mass_Flow_Rate_kg_s = flow_rate_steam_high
        WRITE(out_lun,'(A)') "wzHTC_Warning: Neither high nor low load"
    ENDIF

    ! obtain Cp from NIST database
    if (Temperature_Steam_degC < 400.0)  then
      write(6, *) 'NIST Cp(temp) profile only for T > 400 C'
      STOP
    else if (Temperature_Steam_degC > 700.0)  then
      write(6, *) 'NIST Cp(temp) profile only for T < 700 C'
      STOP
    else if (Temperature_Steam_degC >= 400.0 .and. &
      & Temperature_Steam_degC < 450.0)  then

      if (Pressure_MPa > 24.0)  then  ! supercritical, Weston, Genon

      Cp_heat_capacity_J_kg_degC = Cp_water1_24MPa(4)
      do i = 3, 1, -1
        Cp_heat_capacity_J_kg_degC = Cp_water1_24MPa(i) + &
        & Cp_heat_capacity_J_kg_degC * Temperature_Steam_degC
      end do

      else if (Pressure_MPa < 19.2)  then  ! supercritical, Czeck boiler at 19.1 MPa
        write(6, *) 'STOP: subcritical data at P=19.1 MPa at T=[480:690] C'
        STOP
      endif

    else if (Temperature_Steam_degC >= 450.0 .and. &
      & Temperature_Steam_degC < 550.0)  then

      if (Pressure_MPa > 24.0)  then  ! supercritical, Weston, Genon
        Cp_heat_capacity_J_kg_degC = Cp_water2_24MPa(4)
        do i = 3, 1, -1
          Cp_heat_capacity_J_kg_degC = Cp_water2_24MPa(i) + &
          & Cp_heat_capacity_J_kg_degC * Temperature_Steam_degC
        end do

      else if (Pressure_MPa < 19.2)  then  ! subcritical, Czeck boiler at 19.1 MPa
        if (Temperature_Steam_degC >= 450.0 .and. &
      &   Temperature_Steam_degC < 475.0)  then

          write(6, *) 'STOP: subcritical data at P=19.1 MPa at T=[480:690] C'
          STOP

        else if (Temperature_Steam_degC >= 475.0 .and. &
      &   Temperature_Steam_degC < 550.0)  then

          Cp_heat_capacity_J_kg_degC = Cp_water1_19MPa(4)
          do i = 3, 1, -1
            Cp_heat_capacity_J_kg_degC = Cp_water1_19MPa(i) + &
          & Cp_heat_capacity_J_kg_degC * Temperature_Steam_degC
          end do

        endif

      endif

    else if (Temperature_Steam_degC >= 550.0 .and. &
      & Temperature_Steam_degC < 700.0)  then

      if (Pressure_MPa > 24.0)  then  ! supercritical, Weston, Genon
      Cp_heat_capacity_J_kg_degC = Cp_water3_24MPa(4)
      do i = 3, 1, -1
        Cp_heat_capacity_J_kg_degC = Cp_water3_24MPa(i) + &
        & Cp_heat_capacity_J_kg_degC * Temperature_Steam_degC
      end do
      else if (Pressure_MPa < 19.2)  then  ! subcritical, Czeck boiler at 19.1 MPa

          Cp_heat_capacity_J_kg_degC = Cp_water1_19MPa(4)
          do i = 3, 1, -1
            Cp_heat_capacity_J_kg_degC = Cp_water1_19MPa(i) + &
          & Cp_heat_capacity_J_kg_degC * Temperature_Steam_degC
          end do

      endif

    endif

     !  Cp_heat_capacity_J_kg_degC = 3.0e+3 until 12/12/13

    STEAM_PIPE_EXPON = 2.0D0 * math_PI * tube_outer_radius * htc_total / &
        & (Mass_Flow_Rate_kg_s * Cp_heat_capacity_J_kg_degC)

RETURN
END FUNCTION STEAM_PIPE_EXPON

END MODULE STEAM_MODULE
