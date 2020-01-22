! =====================================================================
! LICORICE: Ly_alpha transfer part, on a fixed grid 
! =====================================================================

PROGRAM TREECD

USE COSMOLOGY
USE FUNCTIONS
USE IO
USE TR_LY
USE TR_LY_UTILS
USE VARS
USE OMP_LIB
USE MPI

implicit none

integer     :: ix,iy,iz
real(KIND=8):: a,x_in,x_out,Ho,nb_dens,cell_path,Ly_alpha_cross,cell_opt_depth,target_opt_depth,target_x,t_v,aa

!       Open data files, read input parameters and initial system state,
!       initialize system parameters, and start output.
!       ----------------------------------------------------------------

 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,PROVIDED,ierror)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierror)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
 IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
   print*,'MPI initialized'
 END IF

 call INPARM

 !$OMP PARALLEL
 call init_rando()
 !$OMP END PARALLEL

 IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
   open(50,FILE='dtb_moy.dat',STATUS='REPLACE')
   open(44,FILE='sources.dat',STATUS='REPLACE')
 print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
 print*,' LICORICE: LOOP OVER THE SNAPSHOT FILES'
 print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
 END IF

 do num_file=first_file,last_file

   IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
     print*,'processing file:',num_file
   END IF

   call LICORICE_NEW_INPUT(num_file)
   call expansion_factor(tnow,exp_fact,exp_fact_dot)
   ! Photons Lyman-alpha X; factor 0.5 because photons at the line center
   ! ENLEVER CA S'IL N'Y A PAS DE RAYONS X!!!!!!!!
   !do ix=1,array_size
   !  do iy=1,array_size
   !    do iz=1,array_size
   !      diff_num(ix,iy,iz)=diff_num(ix,iy,iz)*1.0d65*0.5*1.3475d-7*cell(ix,iy,iz)%HI_number_density* &
   !                         exp_fact**(-3.0)*exp_fact/exp_fact_dot
   !    enddo
   !  enddo
   !enddo
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

   IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
     print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
     print*,'redshift initial: ',1./exp_fact-1
     print*,'exp_fact: ',exp_fact
     print*,'tnow: ',tnow
     print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
   END IF

   call CAST_PHOTONS_LY

   IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
     print*,''
     print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
     print*,'"""""""""""""""""""""""""""""""""""""""""""""""""""""'
     print*,''
   END IF

 enddo  

 IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
   close(50)
   close(44)
 END IF

 CALL MPI_FINALIZE(ierror)

END PROGRAM
