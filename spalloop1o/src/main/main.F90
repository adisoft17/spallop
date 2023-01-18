   Program SPALLOOP
  !===============================================================================
  !                                  SPALLOOP
  !                    Spallation Loop Computer Program
  !
  !  SpalLoop calculates the stress/strain evolution not only through thickness, 
  !      but along an entire boiler loop tube, taking into account the
  !      temperature 
  !      increase in the steam and metal temperature along the tubes.  
  !      SpalLoop code includes the modules for handling: 
  !      (a) multiple layered oxide scales, such as magnetite, spinel and/ore
  !      haematite,
  !      (b) temperature gradient through the tube circumference,
  !      (c) flue temperature and/or heat flux distribution due to flue gases
  !      along 
  !          and entire length of a boiler tube,
  !      (d) different growth temperatures for the spinel and magnetite,
  !      (e) load cycling,
  !      (f) creep in the metal and oxide,
  !      (g) blockage model providing estimates for the amount of 
  !          tube cross-sectional area blocked by the exfoliated oxide,
  !      (h) one or two tube bends, differentiating between blockages  
  !          in the inlet loop bend and outlet loop bend
  !        
  !       UT-BATTELLE, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND 
  !       DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.  THERE ARE NO 
  !       EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A 
  !       PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE 
  !       ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT 
  !       THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE 
  !       OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.  THE USER ASSUMES 
  !       RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, 
  !       CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR 
  !       ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF
  !       THE 
  !       SOFTWARE.
  !
  !  Developed by: Adrian S. Sabau, sabaua@ornl.gov and Ian G. Wright 
  !  Copyright 2018 UT-Battelle LLC.  All rights reserved.
  !
  !================================================================================
  use input_module,              only: READ_INPUT
  use output_module,             only: INITIALIZE_IO
  use input_utilities_module,    only: GET_ARGUMENTS
  use solver_data_module,        only: if_average_oxide 
  use solution_data_module,      only: no_oxide_layers
  use oxide_data_module,         only: no_oxide, no_layer
  use oxide_growth_module,       only: OXIDE_GROWTH_INIT
  use steam_temp_module,         only: GET_HEAT_FLUX_ST_TEMP
  use waterwall_input_module,    only: WATERWALL_INIT
  use stress_driver_module,      only: STRESS_LOCATION
  use mesh_module,               only: MESH_POINTS
  use blockage_module,           only: BLOCKAGE_AREA

  implicit none
  ! Local variables
  integer   :: iter, nproblem
  logical   :: convergent

  ! read and parse the command line
  call GET_ARGUMENTS ()

  ! initalized the I/O 
  call INITIALIZE_IO ()

  ! get the credit
  call CODE_INIT ()
  
  ! read input
  call READ_INPUT ()
  write(6, *)  'DONE input'

  ! store the number of layers
  if (if_average_oxide)  then
    no_oxide_layers = 2 
  else
    no_oxide_layers = no_layer  ! check this
  endif

  ! get the mesh points
  call MESH_POINTS(no_oxide_layers)

  ! initialize the waterwall
  call WATERWALL_INIT()

  ! get hf and T_st at various heights
  call GET_HEAT_FLUX_ST_TEMP()

  ! initialize oxide thickness data
  call OXIDE_GROWTH_INIT ()

  ! get the strain-stress
  call STRESS_LOCATION()

  ! get the blockage area based on Armitt stored elastic energy
  call BLOCKAGE_AREA(no_oxide_layers)

  write(6,*) 'end of program'
  stop

  End Program SPALLOOP

  SUBROUTINE CODE_INIT ()
   !=======================================================================
  ! Purpose:
  !
  !   Determine and print the problem information
  !
  !=======================================================================
  use credit_module,          only: architecture, code_date, code_name, &
                                  code_version, host_name, libraries, &
                                  run_date, Copyright
  use output_module,        only: aux_lun, blank_line, inp_lun, &
                                  input_file, out_lun, Output_String, &
                                  prefix, title, tty_lun
  use utilities_module,     only: TIMESTAMP

  implicit none
  ! Local Variables
  integer  ::  i, k

  code_name = 'waterwall'

  ! Get the code version.
  code_version = '3.0.0'

  ! Get the build date.
  code_date = 'Thu Mar 19 11:31:23 EDT 2009'

  ! Get the architecture.
  architecture = 'i686'

  ! Get the host name.
  host_name = 'cast'

  ! Get the libraries.
  libraries = '-'
  libraries = TRIM(ADJUSTL(libraries)) // ', ' // '-'
  ! Write out code/problem info (including Copyright notice).
  write (tty_lun, 2) TRIM(code_name)
  write (out_lun, 2) TRIM(code_name)
  write (aux_lun, 2) TRIM(code_name)
2 format(36x,a)
  ! 2 format('(/,36x,a)')

  do i = 1,SIZE(Copyright)
     write (tty_lun, 6) TRIM(Copyright(i))
     write (out_lun, 6) TRIM(Copyright(i))
     write (aux_lun, 6) TRIM(Copyright(i))
  end do
6 format( 9x, a) ! '(9x,a)')

  write (tty_lun, 3) TRIM(code_version), TRIM(architecture), TRIM(libraries)
  write (out_lun, 3) TRIM(code_version), TRIM(architecture), TRIM(libraries)
  write (aux_lun, 3) TRIM(code_version), TRIM(architecture), TRIM(libraries)
3 format(/,1x,80('='),//,32x,'Problem/Code Specs',//, &
         ' Version: ',a,//,' Architecture: ',a,//,' Libraries: ',a)

  ! Print compilation date, run date, and problem title.

  read (inp_lun, '(a)') title

  call TIMESTAMP (run_date)
  k = LEN_TRIM(code_date) - 9
  Output_String = blank_line
  write (tty_lun, 4) code_date(5:k), run_date(5:22), TRIM(title)
  write (out_lun, 4) code_date(5:k), run_date(5:22), TRIM(title)
  write (aux_lun, 4) code_date(5:k), run_date(5:22), TRIM(title)
4 format(/,' Build Date/Time: ',a,//,' Execution Date/Time: ',a, &
       //,' Problem: ',a,//,1x,80('='))

  ! Print problem initialization header.
  input_file = TRIM(prefix) // '.inp'
  write (tty_lun, 5) TRIM(input_file)
  write (out_lun, 5) TRIM(input_file)
5 format(/,35x,'Input Phase',//,' Parsing input file ',a)

  return
END SUBROUTINE CODE_INIT
