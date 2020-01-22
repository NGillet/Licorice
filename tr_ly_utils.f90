MODULE TR_LY_UTILS
USE FUNCTIONS
USE COSMOLOGY
USE VARS
USE OMP_LIB



logical,parameter :: USE_ACCEL=.FALSE., USE_X_ALPHA_ACCEL=.TRUE.

CONTAINS 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TRACE_PHOT_LY(p)
implicit none

logical, parameter :: TR_PERIODIC=.TRUE.
logical            :: check
integer, INTENT(IN) :: p
integer :: j,ir,i,jj,indice,nb_diff,l,ix,iy,iz,kk
integer :: in_face,out_face,err,ct,err1,n_iter,n_loop,n_loop2
integer, dimension(3) :: cur_cell,next_cell
real    :: r
real(KIND=8) :: total_box_luminosity
real(KIND=8), dimension(3) :: cur_pos,cur_dir,next_pos,tmp_pos,next_dir,vel_cel_macro_over_c,out_pos,inv_cur_dir
real(KIND=8) :: cur_opt_depth,target_opt_depth,loc_opt_depth,cell_opt_depth_at_line_center
real(KIND=8) :: local_Ly_alpha_cross_section,a,x,T4,inv_sqrtT4,density_conv,dist_conv,Ho,v_therm_over_c
real(KIND=8) :: local_Ly_alpha_cross_section_at_line_center,vel_conv,x_in,x_out
real(KIND=8), dimension(2) :: a12,x_in12,x_out12,Ho12,local_Ly_alpha_cross_section_at_line_center12,local_Ly_alpha_cross_section12
real(KIND=8) :: path_length,next_freq,tot_path_length,x_core,atau,local_cell_freq
real(KIND=8) :: t_phot,t_phot_old,t_out,local_out_freq,exp_fact_old,exp_fact_out,exp_fact_dot_out
real(KIND=8) :: Ly_alpha_const=5.9d-14
real(KIND=8) :: Ly_alpha_freq=2.466d+15
real(KIND=8) :: inv_Ly_alpha_freq=4.05515d-16
real(KIND=8) :: Ly_alpha_line_width=9.936d+7
real(KIND=8) :: target_val,target_x,diffusion_nb,exp_fact_tnow
real(KIND=8) :: np,nb_physical_phot_per_phot,t_moy
real(KIND=8) :: real_tab_index,fac_tab_index,m_part
real(KIND=8) :: half_rsize,inv_ly_alpha_dt,vel_cel_macro_over_c_dot_cur_dir,cell_path_length,cpt,tempor
integer,dimension(NB_PHOT_PER_SNAPSHOT) :: phot_source_index
INTEGER :: int_tab_index,thread_num,file_unit,old_phot_nb,err_scatter_freq
character(3) :: num1,num20

! Higher-order Lyman-series variables
integer                     :: jk,kn,indice_freq,indice12
!-----------------------------------------------------------------------

half_rsize=box_size/2.d0
inv_ly_alpha_dt=1.d0/LY_ALPHA_DT

if(mod(NB_PHOT_PER_SNAPSHOT,nb_thread)/=0) then
  write(*,*)'NB_PHOT_PER_SNAPSHOT should be a multiple of the number of threads'
  STOP
endif
if(mod(NB_PHOT_PER_SNAPSHOT/nb_thread,NB_PHOT_PER_BATCH_PER_THREAD)/=0) then
  write(*,*)'NB_PHOT_PER_SNAPSHOT/nb_thread should be a multiple of NB_PHOT_PER_BATCH_PER_THREAD'
  STOP
endif


exp_fact_tnow=exp_fact

 cpt=0.0d0
if (p == 1) then
  total_box_luminosity=sum(dble(cell(:,:,:)%luminosity))
  IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
    print*,'total_box_luminosity',total_box_luminosity*luminosity_unit
  END IF
  do ix=1,array_size
    do iy=1,array_size
      do iz=1,array_size
        do jj=int(cpt)+1,nint(cpt+NB_PHOT_PER_SNAPSHOT*dble(cell(ix,iy,iz)%luminosity)/total_box_luminosity)
          phot_source_index(jj)=(ix-1)*array_size**2+(iy-1)*array_size+iz
!          if(mod(jj,10000000)==0) print*,jj,myrank
        enddo
      cpt=cpt+NB_PHOT_PER_SNAPSHOT*dble(cell(ix,iy,iz)%luminosity)/total_box_luminosity
      enddo
    enddo
  enddo
  ! Juin 2012
  if(count(phot_source_index == 0) .ne. 0) then
    print*,'Error: ProblËme assignation indices des photons!!!',count(phot_source_index == 0)
    do jj=1,NB_PHOT_PER_SNAPSHOT
       if(phot_source_index(jj) == 0) then
         print*,'Error: ProblËme assignation indice photon',jj,'/',NB_PHOT_PER_SNAPSHOT,'domaine',myrank
         phot_source_index(jj) = phot_source_index(jj-1)
         print*,'Error: Indice attribuÈ:',phot_source_index(jj-1),'de photon jj-1=',jj-1
         print*,'Error: correspondant ‡ x=',(phot_source_index(jj-1)-1)/array_size**2+1
         print*,'Error: y=',(phot_source_index(jj-1)-1-(((phot_source_index(jj-1)-1)/array_size**2+1)-1)*array_size**2)/array_size+1
         print*,'Error: z=',phot_source_index(jj-1)-(((phot_source_index(jj-1)-1)/array_size**2+1)-1)*array_size**2-(((phot_source_index(jj-1)-1-(((phot_source_index(jj-1)-1)/array_size**2+1)-1)*array_size**2)/array_size+1)-1)*array_size
       end if
    end do
    if(count(phot_source_index == 0) == 0) then
      print*,'Error: ProblËme rÈsolu'
    else
      print*,'Error: ProblËme non rÈsolu!!!'
      call sleep(20)
      STOP
    end if
  end if
endif

!$OMP PARALLEL PRIVATE(cur_cell,next_cell,in_face,out_face,err,err1,ct,cur_pos,cur_dir,inv_cur_dir,next_pos,tmp_pos,next_dir,vel_cel_macro_over_c,out_pos,cur_opt_depth,target_opt_depth,loc_opt_depth,cell_opt_depth_at_line_center,local_Ly_alpha_cross_section,a,x,T4,inv_sqrtT4,Ho,v_therm_over_c,local_Ly_alpha_cross_section_at_line_center,x_in,x_out,path_length,next_freq,tot_path_length,x_core,atau,local_cell_freq,t_phot,t_phot_old,t_out,exp_fact_old,exp_fact_out,target_val,target_x,ir,exp_fact_dot_out,j,jj,indice,diff_num_lt,nb_diff,check,l,nb_physical_phot_per_phot,real_tab_index,int_tab_index,fac_tab_index,vel_cel_macro_over_c_dot_cur_dir,cell_path_length,thread_num,num1,file_unit,old_phot_nb,n_loop,n_loop2,kk,diff_pos_lt,tempor,jk,kn,indice_freq,a12,x_in12,x_out12,Ho12,local_Ly_alpha_cross_section_at_line_center12,local_Ly_alpha_cross_section12,indice12,err_scatter_freq)


if(omp_get_num_threads() /= nb_thread) then
  print*,'Error: inconsistent threadnumber'
  STOP
endif
if(p==1) nb_inflight_phot_lt=0

if(p == -1 .and. (RESTART.or. NUM_FILE /=FIRST_FILE) ) then
  thread_num=omp_get_thread_num()
  file_unit=250+(thread_num+1)
  write(num1,'(I3.3)') thread_num
  write(num20,'(I3.3)') myrank
  open(file_unit,FILE='inflight/nb_inflight_phot_thread='//trim(adjustl(num1))//'_process='//trim(adjustl(num20))//'.dat',STATUS='OLD')
  read(file_unit,*) old_phot_nb
  close(file_unit)
  n_loop=old_phot_nb/NB_PHOT_PER_BATCH_PER_THREAD+1
  IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
    print*,'Process, thread',myrank,thread_num,'casting ',old_phot_nb,'old photons in ',n_loop,'batches'
  END IF
else if (p ==1) then
  n_loop=NB_PHOT_PER_SNAPSHOT/nb_thread/NB_PHOT_PER_BATCH_PER_THREAD
else
  n_loop=0
endif

! Loop on the batches of photons
do jj=1,n_loop

  IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
!    print*,'Process, Thread',myrank,omp_get_thread_num(),'begin with batch',jj,'of',n_loop
  END IF

  OLD_PHOT_ARRAY_IN_USE_LT=.FALSE.
  n_loop2=NB_PHOT_PER_BATCH_PER_THREAD
  if(p==-1 .and. jj==n_loop) then
    n_loop2=old_phot_nb-(old_phot_nb/NB_PHOT_PER_BATCH_PER_THREAD)*NB_PHOT_PER_BATCH_PER_THREAD
  endif
  if(p==-1) then
    thread_num=omp_get_thread_num()
    file_unit=(thread_num+1)+150
    if(jj==1) then
     write(num1,'(I3.3)') thread_num
     write(num20,'(I3.3)') myrank
     open(file_unit,FILE='inflight/old_inflight_photons_thread='//trim(adjustl(num1))//'_process='//trim(adjustl(num20))//'.dat',FORM='UNFORMATTED',STATUS='OLD')
    endif
    do kk=1,n_loop2
      read(file_unit) old_phot_pos_lt(kk,1:3),old_phot_dir_lt(kk,1:3),old_phot_target_opt_depth_lt(kk), &
                      old_phot_time_lt(kk),old_phot_freq_lt(kk),nb_physical_phot_per_old_phot_lt(kk)
    enddo
  endif
  diff_num_lt=0
  diff_pos_lt=1
!  print*,'Entree boucle photons de ce thread'
! Processing the photons in this batch
  do kk=1,n_loop2
    exp_fact_thread=exp_fact
    if( p == -1) then
      call initialize_old_photon(kk,t_phot,cur_opt_depth,target_opt_depth, &
                             cur_pos,next_pos,cur_dir,in_face,out_face,err,err1)
      inv_cur_dir=1./cur_dir
      nb_physical_phot_per_phot=nb_physical_phot_per_old_phot_lt(kk)
    else
      j=kk+(jj-1)*NB_PHOT_PER_BATCH_PER_THREAD+omp_get_thread_num()*(NB_PHOT_PER_SNAPSHOT/nb_thread)

      call initialize_photon(phot_source_index(j),t_phot,cur_opt_depth,target_opt_depth, &
                           cur_pos,next_pos,cur_dir,in_face,out_face,err,err1)
      inv_cur_dir=1./cur_dir
      nb_physical_phot_per_phot=LY_ALPHA_DT*total_box_luminosity*luminosity_unit &
                               /(planck_cgs*113.d0*ly_alpha_freq/98.d0)       &              ! Entre Lya et Ly_zeta
                               /NB_PHOT_PER_SNAPSHOT &
                               /size                                                     ! MPI
    endif
    real_tab_index=(t_phot-T_INIT_BATCH)/LY_ALPHA_DT*exp_fact_ntab
    int_tab_index=int(real_tab_index)
    fac_tab_index=real_tab_index-int(real_tab_index)
    exp_fact_thread=(1.d0-fac_tab_index)*a_tab(int_tab_index)+fac_tab_index*a_tab(int_tab_index+1)
    inv_exp_fact_thread=1./exp_fact_thread
    exp_fact_dot_thread=(1.d0-fac_tab_index)*adot_tab(int_tab_index)+fac_tab_index*adot_tab(int_tab_index+1)

    if(TR_PERIODIC) then
      do l=1,3
        if(cur_pos(l) >= box_size .and. cur_dir(l) > 0d0) then
          cur_pos(l)=cur_pos(l)-box_size
        end if
        if(cur_pos(l) <= 0.d0 .and. cur_dir(l) < 0d0) then
          cur_pos(l)=cur_pos(l)+box_size
        end if
      enddo
    endif
    next_pos=cur_pos

    cur_cell(1)=int(cur_pos(1)/box_size*array_size)+1
    cur_cell(2)=int(cur_pos(2)/box_size*array_size)+1
    cur_cell(3)=int(cur_pos(3)/box_size*array_size)+1

    next_cell=cur_cell
    ct=0
    nb_diff=0
    check=.FALSE.
!    cur_opt_depth=0.d0  Commented 05/02/2018. Conflicts with initialize_old_phot


! Propagate this photon until next snapshot, or it reaches the line core
    phot_prop_loop: do while( err == 0 .and. err1==0 .and. t_phot <= tnow+LY_ALPHA_DT ) 
      ct=ct+1
      cur_pos=next_pos
      cur_cell=next_cell
      if(TR_PERIODIC) then
        do l=1,3
          if((cur_pos(l) >= box_size*(1.d0-1.d-9) .and. cur_dir(l) > 0d0) .or.cur_cell(l) > array_size) then
            cur_pos(l)=cur_pos(l)-box_size
            cur_cell(l)=1
          end if
          if((cur_pos(l) <= box_size*1.d-9 .and. cur_dir(l) < 0d0) .or.cur_cell(l) < 1) then
            cur_pos(l)=cur_pos(l)+box_size
            cur_cell(l)=array_size
          end if
        enddo
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If cell is empty, go to next cell (zero density and T => bug )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density == 0.) then
        call exit_point_fast(cur_pos,cur_dir,inv_cur_dir,cur_cell,next_pos,in_face,out_face)
        tmp_pos=next_pos
        t_phot_old=t_phot
        t_phot=t_phot+sqrt(sum((next_pos-cur_pos)**2))*exp_fact_thread/3d+10
        exp_fact_old=exp_fact_thread
        real_tab_index=(t_phot-T_INIT_BATCH)/LY_ALPHA_DT*exp_fact_ntab
        int_tab_index=int(real_tab_index)
        print*,'Error: empty',int_tab_index,t_out,T_INIT_BATCH
        fac_tab_index=real_tab_index-int(real_tab_index)
        exp_fact_thread=(1.d0-fac_tab_index)*a_tab(int_tab_index)+fac_tab_index*a_tab(int_tab_index+1)
        inv_exp_fact_thread=1./exp_fact_thread
        exp_fact_dot_thread=(1.d0-fac_tab_index)*adot_tab(int_tab_index)+fac_tab_index*adot_tab(int_tab_index+1)

        PHOT_freq=PHOT_freq*exp_fact_old*inv_exp_fact_thread
        v_therm_over_c=1d-6  ! (for err test later on)

        next_cell=cur_cell
        if(out_face==1) next_cell(1)=next_cell(1)+1
        if(out_face==2) next_cell(1)=next_cell(1)-1
        if(out_face==3) next_cell(2)=next_cell(2)+1
        if(out_face==4) next_cell(2)=next_cell(2)-1
        if(out_face==5) next_cell(3)=next_cell(3)+1
        if(out_face==6) next_cell(3)=next_cell(3)-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If cell is not empty... compute transfer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else

        if(abs(dble(cur_cell(1))-0.5d0 - (cur_pos(1)/box_size*array_size) ) >0.51 .OR. &
           abs(dble(cur_cell(2))-0.5d0 - (cur_pos(2)/box_size*array_size) ) >0.51 .OR. &
           abs(dble(cur_cell(3))-0.5d0 - (cur_pos(3)/box_size*array_size) ) >0.51 ) then
          print*,'Error2: photon not in cell!!!'
          print*,'Error2: PHOT_freq',PHOT_freq
          print*,'Error2: cur_pos',cur_pos
          print*,'Error2: next_pos',next_pos
          print*,'Error2: cur_dir',cur_dir
          print*,'Error2: cur cell',cur_cell
          print*,'Error2: should be',cur_pos(:)/box_size*array_size+1
          print*,'Error2: in_face,out_face,nb_diff',in_face,out_face,nb_diff
          print*,'Error2: previous diffusion:',path_length,clsize
          print*,'Error2: fixing.......'
          if(TR_PERIODIC) then
            do l=1,3
              if (cur_pos(l) >= box_size .and. cur_dir(l) > 0) then
              cur_pos(l)=cur_pos(l)-box_size*int(cur_pos(l)/box_size)
              end if
              if(cur_pos(l) <= 0. .and. cur_dir(l) < 0) then
              cur_pos(l)=cur_pos(l)+box_size*(int(-cur_pos(l)/box_size)+1)
              end if
            enddo
          endif
          cur_cell(1)=int(cur_pos(1)/box_size*array_size)+1
          cur_cell(2)=int(cur_pos(2)/box_size*array_size)+1
          cur_cell(3)=int(cur_pos(3)/box_size*array_size)+1
        endif

        T4=dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%T4)
        inv_sqrtT4=1.d0/sqrt(T4)
        vel_cel_macro_over_c=dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%vel(:))*exp_fact_thread/3.d+10
        vel_cel_macro_over_c_dot_cur_dir=sum(vel_cel_macro_over_c*cur_dir)
        local_cell_freq=PHOT_freq*(1-vel_cel_macro_over_c_dot_cur_dir)
        indice_freq = floor(1.d0/sqrt(1.0d0-local_cell_freq/3.288d15))            ! Indicates toward which transition the photon is redshifting
        if(indice_freq.gt.6) indice_freq=6
        if(local_cell_freq.lt.2.466d15) indice_freq=2
        if (indice_freq.eq.1) then
           print*,'Error: indice_freq = 1',indice_freq,phot_freq
           STOP
        endif

        a12(1)=3.548d-12*inv_sqrtT4*natural(indice_freq)/(1.0d0-1.0d0/dble(indice_freq)**(2.d0))
        a12(2)=3.548d-12*inv_sqrtT4*natural(indice_freq+1)/(1.0d0-1.0d0/dble(indice_freq+1)**(2.d0))

        if(a12(1) < 1.d-6) print*,'Error: Warning 1(1)',a12(1),T4,cell(cur_cell(1),cur_cell(2),cur_cell(3))%T4,cur_cell, &
                                  'indice_freq = ',indice_freq
        if(a12(2) < 1.d-6.and.indice_freq.ne.6) print*,'Error: Warning 1(2)',a12(2),T4*10000.,cur_cell, &
                                  'indice_freq = ',indice_freq+1
        v_therm_over_c=1.d0/(2.33d+4*inv_sqrtT4)

  
        call exit_point_fast(cur_pos,cur_dir,inv_cur_dir,cur_cell,next_pos,in_face,out_face)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute loc_opt_depth, optical depth in curent cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        cell_path_length=sqrt(sum((next_pos-cur_pos)**2)) 
! Compute x at cur_pos
        x_in12(1)=2.33d+4*inv_sqrtT4*(local_cell_freq/(3.288d15*(1.0d0-1.0d0/dble(indice_freq)**(2.d0))) - 1.d0)
        x_in12(2)=2.33d+4*inv_sqrtT4*(local_cell_freq/(3.288d15*(1.0d0-1.0d0/dble(indice_freq+1)**(2.d0))) - 1.d0)
! Compute x at next_pos
        t_out=t_phot+cell_path_length*exp_fact_thread/3d+10

        real_tab_index=(t_out-T_INIT_BATCH)*inv_LY_ALPHA_DT*exp_fact_ntab
        int_tab_index=int(real_tab_index)
        fac_tab_index=real_tab_index-int(real_tab_index)
         
        if(int_tab_index > 2*exp_fact_ntab) then
          print*,'Error2: T_INIT_BATCH,t_out,LY_ALPHA_DT',T_INIT_BATCH,t_out,LY_ALPHA_DT
          print*,'Error2: t_phot,cell_path_length',t_phot,cell_path_length
          print*,'Error2: next_pos',next_pos
          print*,'Error2: cur_pos',cur_pos
          print*,'Error2: cur_dir',cur_dir
          print*,'Error2: in_face,out_face',in_face,out_face
          STOP
        endif

        exp_fact_out=(1.d0-fac_tab_index)*a_tab(int_tab_index)+fac_tab_index*a_tab(int_tab_index+1)

        local_cell_freq=PHOT_freq*exp_fact_thread/exp_fact_out*(1-vel_cel_macro_over_c_dot_cur_dir) 
        x_out12(1)=2.33d+4*inv_sqrtT4*(local_cell_freq/(3.288d15*(1.0d0-1.0d0/dble(indice_freq)**(2.d0))) - 1.d0)
        x_out12(2)=2.33d+4*inv_sqrtT4*(local_cell_freq/(3.288d15*(1.0d0-1.0d0/dble(indice_freq+1)**(2.d0))) - 1.d0)

        if(USE_X_ALPHA_ACCEL) then
          if(x_out12(1) > 4.d0 .or. x_in12(1)< -4.d0) then
            Ho12(1)=a12(1)*inv_sqrt_pi/(x_in12(1)*x_out12(1))
          else
            if(x_in12(1)==x_out12(1)) then
              call Voigt(a12(1),(x_in12(1)+x_out12(1))/2.d0,Ho12(1))
            else
              call average_Voigt(a12(1),x_in12(1),x_out12(1),Ho12(1))
            endif
          endif
        else
          if(abs(x_out12(1)-x_in12(1)) > 0.1d0) then
            call average_Voigt(a12(1),x_in12(1),x_out12(1),Ho12(1))
          else
            call Voigt(a12(1),(x_in12(1)+x_out12(1))/2.d0,Ho12(1))
          endif
        endif

        if(USE_X_ALPHA_ACCEL) then
          if(x_out12(2) > 4.d0 .or. x_in12(2)< -4.d0) then
            Ho12(2)=a12(2)*inv_sqrt_pi/(x_in12(2)*x_out12(2))
          else
            if(x_in12(2)==x_out12(2)) then
              call Voigt(a12(2),(x_in12(2)+x_out12(2))/2.d0,Ho12(2))
            else
              call average_Voigt(a12(2),x_in12(2),x_out12(2),Ho12(2))
            endif
          endif
        else
          if(abs(x_out12(2)-x_in12(2)) > 0.1d0) then
            call average_Voigt(a12(2),x_in12(2),x_out12(2),Ho12(2))
          else
            call Voigt(a12(2),(x_in12(2)+x_out12(2))/2.d0,Ho12(2))
          endif
        endif

        local_Ly_alpha_cross_section_at_line_center12(1)=3.487d2*inv_sqrtT4*osc_strength(indice_freq)/(3.288d15*(1.0d0-1.0d0/dble(indice_freq)**(2.d0)))
        local_Ly_alpha_cross_section_at_line_center12(2)=3.487d2*inv_sqrtT4*osc_strength(indice_freq+1)/(3.288d15*(1.0d0-1.0d0/dble(indice_freq+1)**(2.d0)))
        local_Ly_alpha_cross_section12(1)=local_Ly_alpha_cross_section_at_line_center12(1)*Ho12(1)
        local_Ly_alpha_cross_section12(2)=local_Ly_alpha_cross_section_at_line_center12(2)*Ho12(2)
        local_Ly_alpha_cross_section=max(local_Ly_alpha_cross_section12(1),local_Ly_alpha_cross_section12(2))

        if(local_Ly_alpha_cross_section==local_Ly_alpha_cross_section12(1))then
         indice12 = 1
        else
         indice12 = 2
        endif
        local_Ly_alpha_cross_section_at_line_center=local_Ly_alpha_cross_section_at_line_center12(indice12)
        a=a12(indice12)
        x_in=x_in12(indice12)
        x_out=x_out12(indice12)
        if (indice_freq.eq.6) then
            local_Ly_alpha_cross_section_at_line_center=local_Ly_alpha_cross_section_at_line_center12(1)
            a=a12(1)
            x_in=x_in12(1)
            x_out=x_out12(1)
            local_Ly_alpha_cross_section=local_Ly_alpha_cross_section12(1)
        end if

        loc_opt_depth=cell_path_length*dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density) &
                      *local_Ly_alpha_cross_section*inv_exp_fact_thread**2
        cell_opt_depth_at_line_center=clsize*dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density) &
                                     *local_Ly_alpha_cross_section_at_line_center*inv_exp_fact_thread**2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If no absorption, send on to next cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(cur_opt_depth+loc_opt_depth < target_opt_depth) then

          cur_opt_depth=cur_opt_depth+loc_opt_depth
          t_phot_old=t_phot
          t_phot=t_out
          exp_fact_old=exp_fact_thread

          real_tab_index=(t_phot-T_INIT_BATCH)*inv_LY_ALPHA_DT*exp_fact_ntab
          int_tab_index=int(real_tab_index)
          fac_tab_index=real_tab_index-int(real_tab_index)
          exp_fact_thread=(1.d0-fac_tab_index)*a_tab(int_tab_index)+fac_tab_index*a_tab(int_tab_index+1)
          inv_exp_fact_thread=1.d0/exp_fact_thread
          exp_fact_dot_thread=(1.d0-fac_tab_index)*adot_tab(int_tab_index)+fac_tab_index*adot_tab(int_tab_index+1)

          PHOT_freq=PHOT_freq*exp_fact_old*inv_exp_fact_thread
          tmp_pos=next_pos

          next_cell=cur_cell
          if(out_face==1) next_cell(1)=next_cell(1)+1
          if(out_face==2) next_cell(1)=next_cell(1)-1
          if(out_face==3) next_cell(2)=next_cell(2)+1
          if(out_face==4) next_cell(2)=next_cell(2)-1
          if(out_face==5) next_cell(3)=next_cell(3)+1
          if(out_face==6) next_cell(3)=next_cell(3)-1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If absorption, compute new photon parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else
          if(USE_X_ALPHA_ACCEL) then
            if(x_out /= x_in) then
              target_val=(target_opt_depth-cur_opt_depth)/ &                                     ! target_val = integrale de H
                       ( dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)* &
                       local_Ly_alpha_cross_section_at_line_center*cell_path_length )* &
                       exp_fact_thread**2*(x_out-x_in)
              target_x=1.d0/(1.d0/x_in-target_val/a*sqrt(pi))                                        ! par d√©faut: int√©grale de 1/x2
              if(abs(target_x) < 4d0 .or. abs(x_in) < 4.0d0) then
                target_val=(target_opt_depth-cur_opt_depth)/ &
                       ( dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)* &
                       local_Ly_alpha_cross_section_at_line_center*cell_path_length )* &
                       exp_fact_thread**2*(x_out-x_in)
                call compute_integral_path(target_x,x_in,a,target_val)
              endif
              path_length=(target_x-x_in)/(x_out-x_in)*cell_path_length
            else
              path_length=0.d0
            endif
          else
            if(abs(x_out-x_in) < 0.1d0) then
              path_length=(target_opt_depth-cur_opt_depth)/ &
                        ( dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density) &
                        *local_Ly_alpha_cross_section )*exp_fact_thread**2
            else
              target_val=(target_opt_depth-cur_opt_depth)/ &
                       ( dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)* &
                       local_Ly_alpha_cross_section_at_line_center*cell_path_length )* &
                       exp_fact_thread**2*(x_out-x_in)
              call compute_integral_path(target_x,x_in,a,target_val)
              path_length=(target_x-x_in)/(x_out-x_in)*cell_path_length
            endif
          endif
          if (path_length.lt.0d0) then
             print*,'Error:path_length nÈgatif',PHOT_freq,target_x,x_in,x_out,cell_path_length,path_length,target_val,a,10.*v_therm_over_c
             if(abs(path_length) < 1.d-9*cell_path_length)then
               print*,'Error:fixing to path_length = 0'
               path_length=0d0
             endif
          end if
          if (target_x > x_in) then
             print*,'Error: Problem: target_x > x_in: ',x_in, target_x
             print*,'Error: indice_freq',indice_freq
          end if

          if(path_length > cell_path_length*1.000001d+0) then
            print*,'Error: !!!!!!!!!Warning!!!!!!!!!!!!'
            print*,'Error: PHOT_freq',PHOT_freq
            print*,'Error: path_length,cell_path_length',path_length,cell_path_length,x_in,x_out,target_x,cur_opt_depth,loc_opt_depth,target_opt_depth
            call find_core(a,tempor)
            print*,'Error:',a,tempor,cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density,ct,nb_diff
          else if(path_length > cell_path_length*0.999999d+0) then
            path_length=cell_path_length*0.999999d+0
          endif

          x=target_x
          next_pos=cur_pos+cur_dir*path_length
          next_cell=cur_cell


          t_phot_old=t_phot

          t_phot=t_phot+path_length*exp_fact_thread/3d+10
          exp_fact_old=exp_fact_thread
          real_tab_index=(t_phot-T_INIT_BATCH)*inv_LY_ALPHA_DT*exp_fact_ntab
          int_tab_index=int(real_tab_index)
          fac_tab_index=real_tab_index-int(real_tab_index)
          exp_fact_thread=(1.d0-fac_tab_index)*a_tab(int_tab_index)+fac_tab_index*a_tab(int_tab_index+1)
          inv_exp_fact_thread=1.d0/exp_fact_thread
          exp_fact_dot_thread=(1.d0-fac_tab_index)*adot_tab(int_tab_index)+fac_tab_index*adot_tab(int_tab_index+1)

          PHOT_freq=PHOT_freq*exp_fact_old*inv_exp_fact_thread

          atau=cell_opt_depth_at_line_center*a 
          if(atau < 10.d0) then
            x_core=0.02d0
          elseif (atau < 100.d0) then
            x_core=0.1d0
          elseif (atau < 1000.d0) then
            x_core=0.8d0
          elseif (atau < 10000.d0) then
            x_core=2.d0
          elseif (atau > 10000.d0) then
            x_core=3.d0
          endif
          call scatter_dir(cur_dir,next_dir)
          if(atau < 1d0 .OR. abs(x) > abs(x_core) .OR. .not.(USE_ACCEL)) then
! NORMAL SCATTERING
            call scatter_freq(a,x,v_therm_over_c,cur_cell,cur_dir,next_dir,PHOT_freq,next_freq,err_scatter_freq)
            if(err_scatter_freq /=0) then
              print*,'Error: err_scatter_freq /= 0, killing photon'
              OLD_PHOT_ARRAY_IN_USE_LT(kk)=.FALSE.
              exit phot_prop_loop
            endif
          else
! ACCELERATED SCHEME (skip core scaterring, e.g. TASITSIOMI 2006)
            call scatter_freq_accel(x_core,v_therm_over_c,cur_cell,next_dir,next_freq)
            print*, 'Tasitsiomi!'
            STOP
          endif
          nb_diff=nb_diff+1
          cur_dir=next_dir
          inv_cur_dir=1./cur_dir
          PHOT_freq=next_freq
          cur_opt_depth=0.d0
          r=rando(iseed)
          target_opt_depth=-alog(r)
          in_face=0
          out_face=0

          ! CASCADES
          if(indice_freq.gt.2.and.phot_freq > 3.288d15*(1.0d0-1.0d0/dble(indice_freq)**(2.d0))*(1d0+2.d0*v_therm_over_c).and.&
            phot_freq < 3.288d15*(1.0-1.0/dble(indice_freq+1)**(2.d0))*(1d0-2.d0*v_therm_over_c)) then
            r=rando(iseed)
            if(r.le.prob_casc(indice_freq)) then                          ! cascade to the center of Ly-alpha
              if(indice12==1.or.indice_freq==6) then
                diff_num_lt(kk)=nb_physical_phot_per_phot*1.3475d-7*dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)*exp_fact_thread**(-3.0)* &
                                exp_fact_thread/exp_fact_dot_thread*frec(indice_freq)
              else
                diff_num_lt(kk)=nb_physical_phot_per_phot*1.3475d-7*dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)*exp_fact_thread**(-3.0)* &
                                exp_fact_thread/exp_fact_dot_thread*frec(indice_freq+1)
              end if
              diff_pos_lt(kk)=(cur_cell(1)-1)*array_size**2+(cur_cell(2)-1)*array_size+cur_cell(3)
              err1=1
              OLD_PHOT_ARRAY_IN_USE_LT(kk)=.FALSE.
              check=.TRUE.
            endif
          endif

        ! Fin boucle absorption ou non
        endif
      ! Fin boucle "if cell empty"
      endif

      if(t_phot > T_INIT_BATCH + LY_ALPHA_DT) then
        err1=1
        indice=kk
        OLD_PHOT_time_lt(indice)=t_phot
        OLD_PHOT_ARRAY_IN_USE_LT(indice)=.TRUE.
        OLD_PHOT_FREQ_LT(indice)=phot_freq
        OLD_PHOT_POS_LT(indice,:)=next_pos
        OLD_PHOT_DIR_LT(indice,:)=cur_dir
        old_phot_target_opt_depth_lt(indice)=target_opt_depth-cur_opt_depth
        nb_physical_phot_per_old_phot_lt(indice)=nb_physical_phot_per_phot
      endif

      if(phot_freq < 3.288d15*(1.0d0-1.0d0/dble(indice_freq)**(2.d0))*(1d0+2.d0*v_therm_over_c).or.&
         phot_freq > 3.288d15*(1.0d0-1.0d0/dble(indice_freq+1)**(2.d0))*(1d0-2.d0*v_therm_over_c)) then

        diff_num_lt(kk)=nb_physical_phot_per_phot*1.3475d-7*dble(cell(cur_cell(1),cur_cell(2),cur_cell(3))%HI_number_density)*exp_fact_thread**(-3.0d0)* &
                        exp_fact_thread/exp_fact_dot_thread*frec(indice_freq)
        diff_pos_lt(kk)=(cur_cell(1)-1)*array_size**2+(cur_cell(2)-1)*array_size+cur_cell(3)
        err1=1
        OLD_PHOT_ARRAY_IN_USE_LT(kk)=.FALSE.
        check=.TRUE.
      endif

    ! Fin propagation du photon
    end do phot_prop_loop
  ! Fin de la boucle sur les photons de ce batch
  enddo

  ! Ouverture des fichiers pour photons en vol
  thread_num=omp_get_thread_num()
  file_unit=(thread_num+1)+50
  if(p==1 .and. jj==1) then
   write(num1,'(I3.3)') thread_num
   write(num20,'(I3.3)') myrank
   open(file_unit,FILE='inflight/new_inflight_photons_thread='//trim(adjustl(num1))//'_process='//trim(adjustl(num20))//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
  endif
  ! Ecriture des photons
  do kk=1,NB_PHOT_PER_BATCH_PER_THREAD
    if(OLD_PHOT_ARRAY_IN_USE_LT(kk)) then
      nb_inflight_phot_lt=nb_inflight_phot_lt+1
      write(file_unit) old_phot_pos_lt(kk,1:3),old_phot_dir_lt(kk,1:3),old_phot_target_opt_depth_lt(kk), &
                       old_phot_time_lt(kk),old_phot_freq_lt(kk),nb_physical_phot_per_old_phot_lt(kk)
    endif
  enddo
 ! Ajout des diffusions de ce batch
 !$OMP CRITICAL
  do kk=1,NB_PHOT_PER_BATCH_PER_THREAD
    if(diff_num_lt(kk) /= 0.) then
      ix=(diff_pos_lt(kk)-1)/array_size**2+1
      iy=(diff_pos_lt(kk)-1-(ix-1)*array_size**2)/array_size+1
      iz=diff_pos_lt(kk)-(ix-1)*array_size**2-(iy-1)*array_size
      diff_num(ix,iy,iz)=diff_num(ix,iy,iz)+diff_num_lt(kk)
    endif
  enddo
  !$OMP END CRITICAL

! Fin boucle sur les batches
enddo

if(p==-1) then
  thread_num=omp_get_thread_num()
  if(RESTART .or. num_file/=first_file) then
    file_unit=(thread_num+1)+150
    close(file_unit)
  endif
  file_unit=(thread_num+1)+50
  close(file_unit)
  file_unit=(thread_num+1)+250
  write(num1,'(I3.3)') thread_num
  write(num20,'(I3.3)') myrank
  open(file_unit,FILE='inflight/nb_inflight_phot_thread='//trim(adjustl(num1))//'_process='//trim(adjustl(num20))//'.dat',STATUS='REPLACE')
  write(file_unit,*) nb_inflight_phot_lt
  close(file_unit)
endif

!$OMP END PARALLEL


END SUBROUTINE TRACE_PHOT_LY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCATTER_DIR(cur_dir,next_dir)

implicit none
real(KIND=8), dimension(3), intent(in) :: cur_dir
real(KIND=8), dimension(3), intent(out) :: next_dir
real(KIND=8) :: phi1,theta1,r

! isotropic scaterring

!call random_number(r)
r=rando(iseed)
phi1=dble(r)*2.d0*3.141592653589d0
!call random_number(r)
r=rando(iseed)
theta1=acos(r*2.d0-1.d0)
next_dir(1)=cos(phi1)*sin(theta1)
next_dir(2)=sin(phi1)*sin(theta1)
next_dir(3)=cos(theta1)


END SUBROUTINE SCATTER_DIR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCATTER_FREQ(a,x,v_therm_over_c,cur_cell_loc,cur_dir,next_dir,cur_freq,next_freq,err)
implicit none
integer, dimension(3), intent(in) :: cur_cell_loc
integer, intent(out) :: err
real(KIND=8), intent(in) :: a,x,cur_freq,v_therm_over_c
real(KIND=8), dimension(3), intent(in) :: cur_dir,next_dir
real(KIND=8), intent(out) :: next_freq
real(KIND=8), dimension(3) :: atom_vel_over_c, vel_macro_over_c,u1,u2,u3,atom_vel_over_c_glob_axis
real(KIND=8) :: u,theta1,freq_atom_rest_frame,vel_conv,dist_conv,freq_cell_rest_frame
integer :: err_bis

err=0
! atom_vel_over_c(1) in the direction cur_dir
call prob_func(x,a,u,err_bis)
if(err_bis/=0) then
  err=1
else
  atom_vel_over_c(1)=v_therm_over_c*u
! atom_vel perpendicular to cur_dir
  call gauss_prob_func(u)
  atom_vel_over_c(2)=v_therm_over_c/sqrt(2.0d0)*u
  call gauss_prob_func(u)
  atom_vel_over_c(3)=v_therm_over_c/sqrt(2.0d0)*u

  vel_macro_over_c=cell(cur_cell_loc(1),cur_cell_loc(2),cur_cell_loc(3))%vel*exp_fact_thread/3.d+10 

! frequency shift : double Doppler effect, no recoil effect included

  freq_cell_rest_frame = cur_freq*(1d0-sum(vel_macro_over_c*cur_dir))
  freq_atom_rest_frame = freq_cell_rest_frame*(1-atom_vel_over_c(1))

  u1=cur_dir
  if(sqrt(u1(1)*u1(1)+u1(2)*u1(2)) /= 0) then
    u2(1)=u1(2)/sqrt(u1(1)*u1(1)+u1(2)*u1(2))
    u2(2)=-u1(1)/sqrt(u1(1)*u1(1)+u1(2)*u1(2))
    u2(3)=0.d0
  else
    u2(1)=1.d0
    u2(2)=0.d0
    u2(3)=0.d0
  endif
  u3(1)=u1(2)*u2(3)-u1(3)*u2(2)
  u3(2)=u1(3)*u2(1)-u1(1)*u2(3)
  u3(3)=u1(1)*u2(2)-u1(2)*u2(1)


  atom_vel_over_c_glob_axis(1)=atom_vel_over_c(1)*u1(1)+ &
                             atom_vel_over_c(2)*u2(1)+ &
                             atom_vel_over_c(3)*u3(1)
  atom_vel_over_c_glob_axis(2)=atom_vel_over_c(1)*u1(2)+ &
                             atom_vel_over_c(2)*u2(2)+ &
                             atom_vel_over_c(3)*u3(2)
  atom_vel_over_c_glob_axis(3)=atom_vel_over_c(1)*u1(3)+ &
                             atom_vel_over_c(2)*u2(3)+ &
                             atom_vel_over_c(3)*u3(3)


  freq_cell_rest_frame = freq_atom_rest_frame/(1.d0-sum(atom_vel_over_c_glob_axis*next_dir))

  next_freq = freq_cell_rest_frame/(1.d0-sum(vel_macro_over_c*next_dir))
  !print*,'cur_freq',cur_freq,'next_freq',next_freq
endif


END SUBROUTINE SCATTER_FREQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCATTER_FREQ_ACCEL(x_core,v_therm_over_c,cur_cell_loc,next_dir,next_freq)
implicit none
real(KIND=8), intent(in) :: v_therm_over_c,x_core
real(KIND=8), intent(out) :: next_freq
real(KIND=8), dimension(3), intent(in) :: next_dir
integer, dimension(3), intent(in) :: cur_cell_loc
real(KIND=8) :: atom_vel_over_c
real(KIND=8) :: u,theta1,phi1,freq_atom_rest_frame
real(KIND=8) :: Ly_alpha_freq=2.466d+15
real(KIND=8), dimension(3) :: vel_macro_over_c

call vulcano_prob_func_bis(u,x_core*sqrt(2.d0))
atom_vel_over_c=v_therm_over_c/sqrt(2.d0)*u

vel_macro_over_c=cell(cur_cell_loc(1),cur_cell_loc(2),cur_cell_loc(3))%vel*exp_fact_thread/3.d+10 


next_freq = Ly_alpha_freq/(1d0-atom_vel_over_c)
next_freq = next_freq/(1d0-sum(vel_macro_over_c*next_dir))

END SUBROUTINE SCATTER_FREQ_ACCEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FREQUENCE_LY(f)
implicit none

!integer :: i,iseed
integer :: i
real(KIND=4) :: r
real(KIND=8) :: f,a,b


  r=rando(iseed)

!  f=2.466d+15*(1.d+0+r*5./27.)        ! Entre Lya et Lyb
!  f=3.288d+15*(7.5d-1+r*3./16.)       ! Entre Lya et Lyg
  f=3.288d+15*(3./4.+r*45./196.)       ! Entre Lya et Ly_zeta
!  f=3.288d+15*(7.5d-1+r*2.5d-1)       ! Entre Lya et la limite Ly

!  a = -2.116416d-31                                  ! T=50000 K
!  b = 1.48429d-15                                    ! T=50000 K
!  a = 4.013822959d-32                                ! T=100000 K
!  b = 1.294198948d-15                                ! T=100000 K
!  a = 1.1237d-31                                     ! T=150000 K
!  b = 1.23966d-15                                    ! T=150000 K
!  f = 2.466d15-b/2./a + sqrt(b*b+4.*a*r)/2./a        ! Spectre non plat


END SUBROUTINE FREQUENCE_LY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Draw a realisation from the probability function p(u)=exp(-u*u)/((x-u)**2+a**2)`
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PROB_FUNC(y,a,u,err)
implicit none
real(KIND=8), intent(IN) :: y,a
real(KIND=8), intent(OUT) :: u
integer, intent(OUT) ::err
real(KIND=8) :: coeff,u_0,signe,x,norm
real(KIND=8), parameter :: cutoff=3.
real(KIND=8) :: r,p,theta1
LOGICAL :: success
real(KIND=8), parameter :: pi2=3.141592653589/2.
integer :: ctloc
!integer :: iseed

err=0
if(a < 1.d-6) then
  print*,'Error using PROB_FUNC: a too small'
  err=1
endif
ctloc=0
success=.FALSE.
!
! u_0: parameter used to minize the rejection fraction (Zheng et Miralda-escude (2002))
! Empirical fit okay for 0.0001 < a < 0.1 and x > 4
!
signe=sign(dble(1),y)
x=y*signe

if( x < 8d0 .and. err==0) then
! method 1 (Zheng et Miralda-escude (2002)
  if( x > 3d0) then
    u_0=1.85d0-log(a)/6.73d0+log(log(x))
  else
    u_0=0.
  endif
! print*,'u0',u_0
  coeff=exp(-u_0*u_0)
  theta1=atan((u_0-x)/a)
  p=(theta1+pi2)/((1.d0-coeff)*theta1+(1.d0+coeff)*pi2)

  do while(.not.success .and. ctloc < 10000000)
!  call random_number(r)
    r=rando(iseed)
    if(r<p) then
!    call random_number(r)
      r=rando(iseed)
      r=r*(theta1+pi2)-pi2
      u=a*tan(r)+x
!    print*,'u1',u,a,r,x
!    call random_number(r)
      r=rando(iseed)
      if(r < exp(real(-u*u)) ) then
         success=.TRUE.
      endif
    else
!    call random_number(r)
      r=rando(iseed)
      r=r*(-theta1+pi2)+theta1
      u=a*tan(r)+x
!    print*,'u2',u,a,r,x
!    call random_number(r)
      r=rando(iseed)
      if(r < exp(real(-u*u+u_0*u_0)) ) then
         success=.TRUE.
      endif
    endif
    ctloc=ctloc+1
    if(mod(ctloc,10000000)==0) print*,'Error: rrr',omp_get_thread_num(),ctloc,r,u,u_0,a,x
  enddo
else if(err==0) then
  !methode 2
  norm= atan( (cutoff-x)/a)  - atan( (-cutoff-x)/a )
  coeff=a*a+x*(cutoff+x)

  do while(.not.success .and. ctloc < 10000000)

!  call random_number(p)
    p=rando(iseed)

    r=tan(norm*p)

    u= ( r* coeff -cutoff*a ) / (r* (cutoff+x) +a )

!  call random_number(p)
    p=rando(iseed)

    success=(p < exp(-u*u) )
    ctloc=ctloc+1
    if(mod(ctloc,1000000)==0) print*,'Error: lll',omp_get_thread_num(),ctloc,y,a,u
  enddo

endif

u=u*signe
if (ctloc >= 10000000) err=1

END SUBROUTINE PROB_FUNC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Draw from gaussian prob function, Box-Muller method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GAUSS_PROB_FUNC(u)

IMPLICIT NONE
REAL(kIND=8), INTENT(OUT) :: u
REAL(KIND=8) :: rsq
REAL(KIND=8) :: v1,v2
REAL(KIND=8), SAVE :: g
LOGICAL, SAVE :: gaus_stored=.false.
INTEGER :: ctloc
ctloc=0
!integer :: iseed
if (gaus_stored) then
  u=g
  gaus_stored=.false.
else
  rsq=-1.
  do while(rsq <=0.0d0 .OR. rsq > 1.0d0)
!    call random_number(v1)
!    call random_number(v2)
    v1=rando(iseed)
    v2=rando(iseed)
    v1=2.0d0*v1-1.0d0
    v2=2.0d0*v2-1.0d0
    rsq=v1**2+v2**2
    ctloc=ctloc+1  
    if(mod(ctloc,100000)==0) print*,'GAUSS_PROB_FUNC',omp_get_thread_num(),ctloc
  end do
  rsq=sqrt(-2.0d0*log(rsq)/rsq)
  u=v1*rsq
  g=v2*rsq
  gaus_stored=.true.
endif

END SUBROUTINE GAUSS_PROB_FUNC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VULCANO_PROB_FUNC(u,umin)

IMPLICIT NONE
REAL(kIND=8), INTENT(OUT) :: u
REAL(KIND=8), INTENT(IN) :: umin

REAL(KIND=8) :: v1,v2
integer :: ct
!integer :: iseed

u=0.d0
ct=0
do while(abs(u)<abs(umin))
!call random_number(v1)
!call random_number(v2)
v1=rando(iseed)
v2=rando(iseed)

v1=v1*exp(-umin)
v2=v2*2.d0*3.141592653589d0
u=sqrt(-2.0d0*log(v1))*cos(v2)
ct=ct+1
if(mod(ct,1000)==0) print*,'VULCANO',omp_get_thread_num(),ct
enddo


END SUBROUTINE VULCANO_PROB_FUNC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VULCANO_PROB_FUNC_BIS(u,umin)

IMPLICIT NONE
REAL(kIND=8), INTENT(OUT) :: u
REAL(KIND=8), INTENT(IN) :: umin
REAL(KIND=8) :: pmin,p
INTEGER :: ctloc
!integer :: iseed

ctloc=0
if(umin > 4.5d0) then
  print*,'Pb de pr√©cision num√©rique de erf'
  STOP
endif
pmin=my_erf(umin/sqrt(2.d0))
!call random_number(p)
p=1.
do while(p == 1.)
  p=rando(iseed)
  p=p*(1-pmin)+pmin
  ctloc=ctloc+1
  if(mod(ctloc,10000)==0) print*,'Error: xxx',omp_get_thread_num(),ctloc
enddo
u=inv_erf(p)*sqrt(2.)
!call random_number(p)
p=rando(iseed)
if(p>0.5d0) u=-u

END SUBROUTINE VULCANO_PROB_FUNC_BIS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analytical fit of the Voigt function (Tasitsiomi 2006)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VOIGT(alpha,x,H)

implicit none
real(KIND=8), intent(in) :: alpha,x
real(KIND=8), intent(out) :: H
real(KIND=8) :: x2,z,q
real(KIND=8), parameter :: pi=3.141592653589

!x2=x*x
!z=(x2-0.855)/(x2+3.42)
!if( z < 0 ) then
!  q = 0
!else
!  q = z*(1. + 21./x2)*alpha/pi/(x2+1)*(0.1117+z*(4.421+z*(-9.207+5.674*z)))
!endif
!H=sqrt(pi)*(q+exp(-x2)/1.77245385)

real(KIND=8) :: a,b,fa,fb,c,fc,x_core

call find_core(alpha,x_core)

if(abs(x)> x_core) then
  H=alpha/sqrt(pi)/x**2
else
  H=exp(-x*x)
endif


END SUBROUTINE VOIGT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FIND_CORE(alpha1,x_core)

real(KIND=8), intent(in) :: alpha1
real(KIND=8), intent(out) :: x_core
real(KIND=8) :: a,b,c,fa,fb,fc,alpha
real(KIND=8), parameter :: sqrtpi=1.77245385090551
integer :: ct

alpha=alpha1
if( alpha > sqrtpi*exp(-1.d0) ) then
  print*,'Error: Unable to find core',alpha,sqrtpi*exp(-1.)
  print*,'Error: stop'
  STOP 
  alpha=sqrtpi*exp(-1.)-0.00001
endif
a=1.d0 ! maximum of x*exp(-x)
b=30.d0
fb=b*exp(-b)-alpha/sqrtpi
fa=a*exp(-a)-alpha/sqrtpi
ct=0
do while(b-a>0.01d0)
  c=(a+b)/2.
  fc=c*exp(-c)-alpha/sqrtpi
  if(fc*fb > 0.) then
    b=c
    fb=fc
  else
    a=c
    fa=fc
  endif
  ct=ct+1
  if(ct > 100) then
    print*,'Error: infinite loop FIND_CORE',ct
    STOP
  endif
enddo
x_core=sqrt((a+b)/2.)

if(x_core-0.5 < 0.01 .or. 30.-x_core < 0.01) print*,'Error: ALARM VOIGT',x_core


END SUBROUTINE FIND_CORE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AVERAGE_VOIGT(alpha,x_in,x_out,H)

implicit none
real(KIND=8), intent(in) :: alpha,x_in,x_out
real(KIND=8), intent(out) :: H

!real(KIND=8), parameter :: pi=3.141592653589
real(KIND=8) :: x_min,x_max,x_core,a,b,c,fa,fb,fc,v1,v2,v3,v4

x_min=min(x_in,x_out)
x_max=max(x_in,x_out)
!---------------------------Modif f√©vrier 2010 transitions sup√©rieures-------
if(x_out.gt.x_in) then
print*,'Error: x_out > x_in! ',x_out,x_in
!STOP
end if
!---------------------------Modif f√©vrier 2010 transitions sup√©rieures-------

call find_core(alpha,x_core)

v1=alpha/sqrt(pi)

if(x_min < -x_core) then
  if(x_max < -x_core) then
    H=v1*(1./x_min-1./x_max)
  else if (x_max < x_core) then
    H=v1*(1.d0/x_min+1.d0/x_core)+sqrt(pi)/2.d0*(-my_erf(-x_core)+my_erf(x_max))
  else
    H=v1*(1.d0/x_min+2.d0/x_core-1.d0/x_max)+sqrt(pi)*my_erf(x_core)
  endif
else if (x_min < x_core) then
  if (x_max < x_core) then
    H=sqrt(pi)/2.d0*(-my_erf(x_min)+my_erf(x_max))
  else
    H=v1*(1.d0/x_core-1.d0/x_max)+sqrt(pi)/2.d0*(-my_erf(x_min)+my_erf(x_core))
  endif
else
  H=v1*(1.d0/x_min-1.d0/x_max)
endif

!H=H*alpha/sqrt(pi)/(x_max-x_min)
H=H/(x_max-x_min)

END SUBROUTINE AVERAGE_VOIGT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE compute_integral_path(target_x,x_i,alpha,t_v)

implicit none
real(KIND=8), intent(in) :: x_i,alpha,t_v
real(KIND=8), intent(out) :: target_x
real(KIND=8) :: target_val,x_in,v1,v2,v3,v4
real(KIND=8) :: a,b,fa,fb,c,fc,x_core

if(t_v>0d0) then
  print*,'Error: t_v>0 in compute_integral_path'
endif
target_val=-t_v
x_in=-x_i

  call find_core(alpha,x_core)

  v1=alpha/sqrt(pi)
  v2=2d0/sqrt(pi)
  v3=1.d0/x_in
  v4=1.d0/x_core

  if(x_in < -x_core) then
    if(v1*(v3+v4) > target_val ) then
      target_x=1.d0/(v3-target_val/v1)
    else
      target_val=target_val-v1*(v3+v4)
      if((my_erf(x_core)-my_erf(-x_core))/v2 > target_val ) then
        target_x=inv_erf(v2*target_val+my_erf(-x_core))
      else
        target_val=target_val-(my_erf(x_core)-my_erf(-x_core))/v2
        target_x=1.d0/(v4-target_val/v1)
      endif
    endif
  else if (x_in < x_core) then
    if((my_erf(x_core)-my_erf(x_in))/v2 > target_val ) then
      target_x=inv_erf(v2*target_val+my_erf(x_in))
    else
      target_val=target_val-(my_erf(x_core)-my_erf(x_in))/v2
      target_x=1.d0/(v4-target_val/v1)
    endif
  else
    target_x=1.d0/(v3-target_val/v1)
  endif
target_x=-target_x  

if(target_x > x_i) then
  print*,'Error: target_x > x_i in compute_integral_path',target_x,x_i
endif


END SUBROUTINE compute_integral_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_photon(p,t_phot,cur_opt_depth,target_opt_depth, &
                        cur_pos,next_pos,cur_dir,in_face,out_face,err,err1)
implicit none

integer, intent(in) :: p
integer, intent(out) :: in_face,out_face,err,err1
real(KIND=8), dimension(3), intent(out) :: cur_pos,next_pos,cur_dir
real(KIND=8), intent(out) :: cur_opt_depth,t_phot,target_opt_depth
real(KIND=8) :: phi1, theta1,r
integer :: iix,iiy,iiz

iix=(p-1)/array_size**2+1
iiy=(p-1-(iix-1)*array_size**2)/array_size+1
iiz=p-(iix-1)*array_size**2-(iiy-1)*array_size


call frequence_ly(phot_freq)

r=rando(iseed)
t_phot=T_INIT_BATCH+r*LY_ALPHA_DT
cur_opt_depth=0.d0
r=rando(iseed)
target_opt_depth=-dlog(r)
if(target_opt_depth < 0.d0) then
  print*,'Error: Negative target_opt_depth in initialize photon !!!!!!!!!!!!!!!'
endif
cur_pos(1)=(real(iix)-0.5)*box_size/real(array_size)
cur_pos(2)=(real(iiy)-0.5)*box_size/real(array_size)
cur_pos(3)=(real(iiz)-0.5)*box_size/real(array_size)
next_pos=cur_pos
r=rando(iseed)
phi1=dble(r)*2.d0*3.141592653589d0
r=rando(iseed)
theta1=acos(r*2.d0-1.d0)
cur_dir(1)=cos(phi1)*sin(theta1)
cur_dir(2)=sin(phi1)*sin(theta1)
cur_dir(3)=cos(theta1)
in_face=0
out_face=0
err=0
err1=0

! Convert source freq to global comoving frame freq through dopler shift
PHOT_freq=PHOT_freq*(1+sum(cell(iix,iiy,iiz)%vel(:)*exp_fact_thread/3.d+10*cur_dir))


END SUBROUTINE initialize_photon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_old_photon(phot_num,t_phot,cur_opt_depth,target_opt_depth, &
                        cur_pos,next_pos,cur_dir,in_face,out_face,err,err1)
implicit none

integer, intent(in) :: phot_num
integer, intent(out) :: in_face,out_face,err,err1
real(KIND=8), dimension(3), intent(out) :: cur_pos,next_pos,cur_dir
real(KIND=8), intent(out) :: cur_opt_depth,t_phot,target_opt_depth
integer :: i

phot_freq=old_phot_freq_lt(phot_num)

t_phot=old_phot_time_lt(phot_num)
cur_opt_depth=0.
target_opt_depth=old_phot_target_opt_depth_lt(phot_num)
cur_pos=old_phot_pos_lt(phot_num,:)
next_pos=cur_pos
cur_dir=old_phot_dir_lt(phot_num,:)
in_face=0
out_face=0
err=0
err1=0
OLD_PHOT_ARRAY_IN_USE_LT(phot_num)=.FALSE.


END SUBROUTINE initialize_old_photon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE exit_point_fast(c_pos,c_dir,inv_c_dir,cell_indexes,n_pos,entry_face,exit_face)

implicit none

real(KIND=8), DIMENSION(3), INTENT(IN) :: c_pos,inv_c_dir,c_dir
integer, DIMENSION(3), INTENT(IN) :: cell_indexes
real(KIND=8), DIMENSION(3), INTENT(OUT) :: n_pos
INTEGER, INTENT(IN) :: entry_face
INTEGER, INTENT(OUT) :: exit_face
INTEGER, DIMENSION(1) :: ef
real(KIND=8), DIMENSION(3) :: cell_pos
real(KIND=8), DIMENSION(3) :: mu
real(KIND=8) :: coeff,cs2
INTEGER :: AXIS

mu=-1.
cs2=clsize/2.0d0

cell_pos=(cell_indexes-0.5d0)*clsize

mu=(cell_pos+sign(1.d+0,c_dir)*cs2-c_pos)*inv_c_dir

!ef=MINLOC(mu,mu>=0.d0)
!coeff=MINVAL(mu,mu>=0.d0)
!14/03/2018 Modification attempting to deal with rare execptions, appearing
!with of the order of 1d12 photons. Connected to equal values of two compoennts in pos and dir +
!limited accuracy of double precision.
ef=MINLOC(mu,mu>=-1.d-12*cs2)
coeff=MINVAL(mu,mu>=-1.d-12*cs2)

exit_face=ef(1)*2
if(c_dir(ef(1)) > 0d0 ) exit_face=exit_face-1

if(exit_face == entry_face .and. entry_face /= 0) then
  ef=MINLOC(mu,mu>coeff)
  coeff=MINVAL(mu,mu>coeff)
  exit_face=ef(1)*2
  if(c_dir(ef(1)) > 0d0 ) exit_face=exit_face-1
endif

n_pos=c_pos+coeff*c_dir

END SUBROUTINE exit_point_fast

      

END MODULE TR_LY_UTILS
