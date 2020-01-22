MODULE TR_LY
USE VARS

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Casts photons 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CAST_PHOTONS_LY
USE TR_LY_UTILS
USE TSPIN

implicit none

integer :: i,j,k,ind,l,np,ix,iy,iz
character(3) :: num1,num20
real :: t1,t2,t3
REAL(KIND=8), dimension(array_size,array_size) :: diff_num_MPI
!-----------------------------------------------------------------------

 call expansion_factor(tnow,exp_fact,exp_fact_dot)
 T_INIT_BATCH=tnow
 call tabulate_exp_fact(T_INIT_BATCH,exp_fact,T_INIT_BATCH+LY_ALPHA_DT)

 ! A initialiser dans io.f90 à cause des photons Ly-alpha produits par rayons X
! diff_num=0
 ! Vérification que diff_num n'est pas nul
! print*,'Max diff_num avant la propagation des photons:',maxval(diff_num)

 call cpu_time(t1)
 call TRACE_PHOT_LY(1)
 call cpu_time(t3)
 print*,'new phot cpu cost: ',t3-t1,'Process',myrank
 CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

 call cpu_time(t3)
 call TRACE_PHOT_LY(-1)
 call cpu_time(t2)
 print*,'old phot cpu cost: ',t2-t3,'Process',myrank
 ! Vérification que diff_num est différent avant et après la propagation des photons
! print*,'Max diff_num après la propagation des photons:',maxval(diff_num)

 IF((ierror==MPI_SUCCESS).and.(myrank.ne.MASTER))then
!   print*,'Barrier: process',myrank,'is waiting at the barrier while other processes are finishing TRACE_PHOT_LY'
 END IF
 CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
! print*,'MPI_REDUCING diff_num_MPI'

do i=1,array_size
 diff_num_MPI=0.d0
 CALL MPI_REDUCE(diff_num(:,:,i),diff_num_MPI,array_size*array_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
 diff_num(:,:,i)=diff_num_MPI
enddo
! CALL MPI_REDUCE(MPI_IN_PLACE,diff_num,array_size*array_size*array_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
! IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
!   print*,'TEST:',' Process no',myrank,'Diff_num_MPI:',diff_num_MPI(458,137,180),diff_num_MPI(37,160,169)
! END IF

 IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
   call TS                                     ! Computes the spin temperature and the differential brightness temperature
 END IF

 IF((ierror==MPI_SUCCESS).and.(myrank.ne.MASTER))then
!   print*,'Barrier: process',myrank,'is waiting at the barrier while Master is calling subroutine tspin.'
 END IF
 CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

 do i=0,nb_thread-1
  write(num1,'(I3.3)') i
  write(num20,'(I3.3)') myrank
  call execute_command_line('mv inflight/new_inflight_photons_thread='//num1//'_process='//num20//'.dat inflight/old_inflight_photons_thread='//num1//'_process='//num20//'.dat')
 enddo

END SUBROUTINE CAST_PHOTONS_LY

END MODULE TR_LY
