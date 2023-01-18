MODULE SOLVER_MODULE
    !=======================================================================
    ! Purpose(s):
    !
    !   solver interface  
    !
    !  Author(s): Adrian S. Sabau, sabaua@ornl.gov 
    !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
    !=======================================================================

  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: SOLVER_SIMPLE

 CONTAINS

 SUBROUTINE SOLVER_SIMPLE (neq_solver, matrix, solution)

    use solver_data_module, only : nrhs

    implicit none

    external                           DGETRS
    external                           DGETRF

    ! arguments
    integer, intent(IN)  :: neq_solver
    real(kind = 8), dimension(neq_solver, neq_solver), intent(INOUT) :: matrix
    real(kind = 8), dimension(neq_solver), intent(INOUT)      :: solution

    ! local variables
    integer, dimension(neq_solver)                 :: ipiv
    integer                                        :: solver_error

    integer            :: i, j, prepare_in = 1, prepare_out = 2

    do i=1, neq_solver
      ! write(2, 3) (matrix(i, j), j=1, neq_solver)
 3    format(50(1pe9.2, 1x))
    end do
!
!  Factor the matrix.
!
    call DGETRF (neq_solver, neq_solver, matrix, neq_solver, IPIV, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  DGETRF returned INFO = ', solver_error
    write ( *, '(a)' ) '  The matrix is numerically singular.'
    stop
  end if
!
!  Solve the linear system; 
! 
    call DGETRS('N', neq_solver, NRHS, matrix, neq_solver, IPIV, solution, &
      & neq_solver, solver_error)

  if ( solver_error /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DGETRS; Solution procedure failed!'
    write ( *, '(a,i6)' ) '  INFO = ', solver_error
    stop
  end if

    return

 END SUBROUTINE SOLVER_SIMPLE

 SUBROUTINE SET_MATRIX(prepare, nsize, matrix, solution)

    use solver_data_module, only       : neq_solver, fval, jacf, delta_var
    ! arguments
    integer,    intent(IN)            :: prepare
    integer,    intent(IN)            :: nsize
    real(kind = 8), dimension(nsize, nsize),  intent(INOUT)  :: matrix
    real(kind = 8), dimension(nsize),  intent(INOUT)  :: solution

    ! local variables
    integer   :: i, j
    ! allocate matrix and vector for the right size

    if (prepare == 1)  then
      do i = 1, neq_solver
        solution(i) = -fval(i)
        do j = 1, neq_solver
          matrix(i, j) = jacf(i, j)
        end do
      end do
    else if (prepare == 2)  then
      do i = 1, neq_solver
        delta_var(i) = solution(i)
      end do
    else 
      write(6, *) 'SET_MATRIX prepare has wrong option'
      stop
    endif

    return

 END SUBROUTINE SET_MATRIX

 END MODULE SOLVER_MODULE
