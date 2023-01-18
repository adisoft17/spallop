MODULE CREEP_UTILITY_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   functions and subroutines for creep modeling
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================


  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: CREEP_NEW_OXIDE, VON_MISES

 CONTAINS

  REAL FUNCTION VON_MISES(stress_r, stress_theta, stress_z)
  
  ! Von Mises sequivalent stress
  real,   Intent(IN)  :: stress_r, stress_theta, stress_z

  real :: dummy

  dummy = (stress_r - stress_theta)**2 + &
     & (stress_r - stress_z)**2 + (stress_z - stress_theta)**2

  if (dummy > 0.0)  then

    VON_MISES = SQRT(0.5 * dummy) 

  else

    VON_MISES = 0.0

  endif

  return

  END FUNCTION VON_MISES

  REAL FUNCTION CREEP_NEW_OXIDE(option, creep_old1, creep_old2, &
        & r1_old, r2_old, r1_new, r2_new)

  ! 2 is at the interface, 1 is inside

  real, Intent(IN)  :: creep_old1, creep_old2, &
        & r1_old, r2_old, r1_new, r2_new
  character(LEN = 80), intent(IN) :: option

  !local vars
  real  :: fnew, new_ox_term, inside_term

  SELECT CASE (TRIM(option))

    CASE ('old')

      CREEP_NEW_OXIDE = creep_old2

    CASE ('zero')

      CREEP_NEW_OXIDE = 0.0

    CASE ('half_old')

      CREEP_NEW_OXIDE = creep_old2 / 2.0

    CASE ('interpolation_zero')

      fnew = 0.0
      new_ox_term = 2.0 * fnew * (r2_new - r2_old)
      ! p1-o2 term; creep_old1 - should be replaced with p1 interpolated value
      inside_term = (creep_old1 + creep_old2) * (r2_old - r1_new)
      CREEP_NEW_OXIDE = -creep_old1 + &
        & (new_ox_term + inside_term) / (r2_new - r1_new)

    CASE ('interpolation_half')

      fnew = creep_old2 / 2.0
      new_ox_term = 2.0 * fnew * (r2_new - r2_old)
      ! p1-o2 term; creep_old1 - should be replaced with p1 interpolated value
      inside_term = (creep_old1 + creep_old2) * (r2_old - r1_new)
      CREEP_NEW_OXIDE = -creep_old1 + &
        & (new_ox_term + inside_term) / (r2_new - r1_new)

    CASE ('exfoliation')

      CREEP_NEW_OXIDE = 0.0

    CASE DEFAULT

      CREEP_NEW_OXIDE = creep_old2
      write(6, *) 'new_oxide_creep_option must be: ', &
       & 'old, zero, half_old, interpolation_zero, or ', &
       & 'interpolation_half '
      stop

  END SELECT
    
  RETURN

  END FUNCTION CREEP_NEW_OXIDE

END MODULE CREEP_UTILITY_MODULE
