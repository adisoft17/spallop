MODULE UTILITIES_MODULE
  !======================================================================
  ! Purpose(s):
  !
  !   Define general utility routines which are used throughout Telluride.
  !
  !   Public Interface(s):
  !
  !     * call TIMESTAMP (date_time)
  !
  !       Returns the date and time in string date_time.
  !
  ! Contains: 
  !           TIMESTAMP
  !
  ! Author(s): Douglas B. Kothe, LANL T-3 (dbk@lanl.gov)
  !            Anand V. Reddy, Cateripillar (reddy_anand_v@cat.com)
  ! Subroutine obtained from Telluride
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: TIMESTAMP

  ! File Version

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE TIMESTAMP (date_time)
    !=======================================================================
    ! Purpose(s):
    !
    !   Subroutine to build a 26-character string containing current
    !   date and time using DATE_AND_TIME intrinsic.
    !
    !                                    12345678901234567890123456
    !   String returned is of the form:  Fri 20 Aug 93 09:33:35.500
    !
    !   Initial routine obtained from ftp.ora.com:/pub/nutshell/fortran90
    !   as part of the examples for the O'Reilly book "Migrating to
    !   Fortran 90" by James F. Kerrigan.
    !
    !=======================================================================

    implicit none

    ! Argument List
    character(LEN = 26), intent(OUT) :: date_time

    ! Local Variables

    character(LEN = 3), dimension(0:6) :: days
    character(LEN = 3), dimension(12)  :: months

    integer, dimension(8) :: elements
    integer               :: m, y, w

    ! Define Days and Months
    data days   / 'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /  
    data months / 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'  /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! System date and time
    call DATE_AND_TIME (VALUES = elements)

    ! Format date and time.
    INVALID: if (elements(1) /= -HUGE(0)) then
       y = elements(1)
       m = elements(2)

       if (m < 3) then
          m = m + 12
          y = y -  1
       end if

       w = MOD((elements(3) + (13*m-27)/5 + y + y/4 - y/100 + y/400 ), 7)  
       CENTURY: if (elements(1) < 2000 ) then
          elements(1) = elements(1) - 1900
       else
          elements(1) = elements(1) - 2000
       end if CENTURY

       write (date_time,10) days(w), elements(3), months (elements(2)), &
                            elements(1), elements(5), elements(6), &
                            elements(7), elements(8)
10     format (a3, 1x, i2.2, 1x, a3, 1x, i2.2, 1x,   &
            i2.2, ':', i2.2, ':', i2.2, '.', i3.3 ) 
    else INVALID
       date_time = ' '
    end if INVALID

    return

  END SUBROUTINE TIMESTAMP

END MODULE UTILITIES_MODULE
