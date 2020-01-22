MODULE IO
USE VARS
USE FUNCTIONS
USE OMP_LIB
USE MPI
CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE C2RAY_INPUT(num_input)
USE COSMOLOGY
USE VARS
IMPLICIT NONE

integer, intent(in) :: num_input
CHARACTER(LEN=6) red_stg
real(KIND=4),dimension(200) :: redshift_table
integer :: ix,iy,iz,nb,i,nx,ny,nz,record
real(KIND=8) :: val,tf,a_next,t_next,a,adot,dens,ionization_fraction,temp,ionizing_luminosity,rho_crit
real(KIND=8) :: lscale,tscale,lum_conv
real(KIND=8), dimension(array_size,array_size) :: dum

!open(40,FILE='inputs/red.dat',STATUS='OLD')
!read(40,*) nb
!do i=1,min(num_input+1,nb)
!  read(40,*) redshift_table(i)
!enddo

exp_fact=1./(1.+redshift_table(num_input))

write(red_stg,'(F6.3)') redshift_table(num_input)
print*,'current redshift :',red_stg

!OPEN(34,FILE='inputs/'//trim(adjustl(red_stg))//'n_all.dat', FORM='BINARY',STATUS='OLD')
!read(34) nx,ny,nz
if(nx/=array_size .or. ny/=array_size .or. nz/=array_size) then
  write(*,*) 'Incorrect array_size in GETBDSBIN_C2RAY_OUTPUT ',nx,ny,nz,array_size
  STOP
endif
!read(34) cell(:,:,:)%HI_number_density
!close(34)


!OPEN(34,FILE='inputs/'//trim(adjustl(red_stg))//'v_all.dat', FORM='BINARY',STATUS='OLD')
!read(34) nx,ny,nz
if(nx/=array_size .or. ny/=array_size .or. nz/=array_size) then
  write(*,*) 'Incorrect array_size in GETBDSBIN_C2RAY_OUTPUT ',nx,ny,nz,array_size
  STOP
endif
do iz=1,nz
  do iy=1,ny
    do ix=1,nx
!      read(34) cell(ix,iy,iz)%vel
    enddo
  enddo
enddo
!close(34)

! Unit conversion
rho_crit=3.*(hubble*100.*1d+5/3.085d+24)**2/8/3.141592/6.67259d-8
do i=1,3
  cell(:,:,:)%vel(i)=cell(:,:,:)%vel(i)*8./cell(:,:,:)%HI_number_density
enddo
cell(:,:,:)%HI_number_density=cell(:,:,:)%HI_number_density*(rho_crit* &
                      omega_b*(real(array_size)/cubeP3M_array_size)**3/proton_mass*(1-Y_p))
lscale=box_size/cubeP3M_array_size*exp_fact
tscale=2.d0/(3.d0*sqrt(omega_0)*hubble*100.*1d+5/3.085d+24)*exp_fact**2
do i=1,3
  cell(:,:,:)%vel(i)=cell(:,:,:)%vel(i)*lscale/tscale/exp_fact
enddo


!OPEN(34,FILE='inputs/xfrac3d_'//trim(adjustl(red_stg))//'.bin', FORM='BINARY',STATUS='OLD')
!read(34) record,nx,ny,nz,record
if(nx/=array_size .or. ny/=array_size .or. nz/=array_size) then
  write(*,*) 'Incorrect array_size in GETBDSBIN_C2RAY_OUTPUT',nx,ny,nz,array_size
  STOP
endif
!read(34) record
do iz=1,nz
!  read(34) dum
  cell(:,:,iz)%HI_number_density=cell(:,:,iz)%HI_number_density*(1.-dum)
enddo
!close(34)

! No data from temperature
cell(:,:,:)%T4=20./10000.

!OPEN(34,FILE='inputs/'//trim(adjustl(red_stg))//'-coarsest_sources_used_wfgamma.dat', STATUS='OLD')
!read(34,*) nx
do i=1,nx
!  read(34,*) ix,iy,iz,ionizing_luminosity
  cell(ix,iy,iz)%luminosity=ionizing_luminosity
enddo
!close(34)

val=0.00000
print*,'begin compute UNIV_AGE'
call compute_t(tf,val)
UNIV_AGE=tf
print*,'UNIV_AGE =', tf/86400./365./1d+9,' Gyr'
exp_fact=1./(redshift_table(num_file)+1)
call compute_t(tnow,exp_fact)
tnow=tf-tnow
a_next=1./(redshift_table(num_file+1)+1.)
call compute_t(t_next,a_next)
t_next=tf-t_next
LY_ALPHA_DT=t_next-tnow


call compute_t(t_next,dble(0.1))
t_next=tf-t_next
call expansion_factor(t_next,a,adot)
Hz10=adot/a

print*,'Current time (in Myr):',tnow/86400./365./1d+6
print*,'INtegration time until next output :',LY_ALPHA_DT

! First convert from cubep3m units to mass in g
lum_conv=rho_crit*omega_0*(box_size/array_size)**3
! Then convert to number of photons per sec
lum_conv=lum_conv*omega_b/(omega_0*proton_mass)/LY_ALPHA_DT
!Then convert to erg per sec
lum_conv=lum_conv*planck_cgs*3d+8/2.8977d-3*C2ray_source_blackbody_temp
!To enable use of single precision
lum_conv=lum_conv/luminosity_unit

cell(:,:,:)%luminosity=cell(:,:,:)%luminosity*lum_conv


END SUBROUTINE C2RAY_INPUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LICORICE_NEW_INPUT(num_input)
USE COSMOLOGY
USE VARS
IMPLICIT NONE

integer, intent(in) :: num_input
CHARACTER(LEN=3) stg
integer :: mxbody_data_in,i,ix,iy,iz,ix_min,ix_max,iy_min,iy_max,iz_min,iz_max,number_of_stars,number_of_gas,compte,particle_index,number_of_sources
integer(KIND=8) :: compte_8
integer :: ix_per,iy_per,iz_per,ct,domain,j
logical :: in_cell
real(KIND=8) :: val,tf,a_next,t_next,a,adot,part_mass
real(KIND=8) :: rho_crit
real(KIND=8) :: dist,ker,dx,dy,dz,p_weight_sum
real(KIND=8) :: vol_unit=1d+70,t1,t2
real(KIND=8),dimension(3) :: pos,vel
real(KIND=8) :: h,frac_mass,dens,frac_ion,lum,temp
real(KIND=4) :: star_age,stellar_life_time,mass,nb_lya
integer      :: particle_type
integer(KIND=8) ::  particle_id
integer      :: p_type,ID,numfich,ct_type_1,ct_type_2,pos_in_file,funit1,funit2,funit3,funit4,funit,funit5
!real(KIND=4),dimension(array_size,array_size,array_size) :: weight
real(KIND=4), allocatable :: weight(:,:,:)
real(KIND=4), dimension(12) ::datain
CHARACTER(len=3) :: NUMBER1
CHARACTER(len=5) ::  NUMBER2
real(KIND=4), allocatable :: work_array1(:,:),work_array2(:,:)

integer, dimension(:), allocatable :: indices_sources,indices_stars
integer                            :: indice_courant_sources,indice_courant_stars

allocate(work_array1(array_size,array_size))
allocate(work_array2(array_size,array_size))
allocate(weight(array_size,array_size,array_size))


! Initialisation du tableau de diffusions ici plutôt que dans tr_ly.f90
! pour pouvoir inclure les diffusions dues aux photons Ly-alpha créés par
! les rayons X.
diff_num=0.0

IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'=============================================================='
  print*,'                 LECTURE FICHIER',num_input
  print*,'=============================================================='
END IF

 weight=0.
 cell(:,:,:)%vel(1)=0
 cell(:,:,:)%vel(2)=0
 cell(:,:,:)%vel(3)=0
 cell(:,:,:)%T4=0.
 cell(:,:,:)%HI_number_density=0.
 cell(:,:,:)%luminosity=0.
 cell(:,:,:)%HII_number_density=0.
 ct_type_1=0
 ct_type_2=0
 mxbody_data_in=0
 number_of_sources=0
 number_of_stars=0
 number_of_gas=0

compte=0
do domain=myrank,domain_nb_in_input_data-1,size
  write(NUMBER1,'(I3.3)') num_input 
  write(NUMBER2,'(I5.5)') domain
!  print*,'Treating file','rank'//trim(adjustl(NUMBER2))//'_out'//trim(adjustl(NUMBER1))
  OPEN(34,FILE='../simu_GD/snapshots/rank'//trim(adjustl(NUMBER2))//'_out'//trim(adjustl(NUMBER1)),FORM='unformatted',STATUS = 'OLD')
!  print*,'task',myrank,'openning file','../simu_GD/snapshots/rank'//trim(adjustl(NUMBER2))//'_out'//trim(adjustl(NUMBER1))

  read(34) mxbody_data_in
  read(34) exp_fact,numfich

 call cpu_time(t1)

  do i=1,mxbody_data_in
    p_weight_sum=0.
    read(34) particle_id,particle_type,pos(1:3),vel(1:3),temp,frac_ion,star_age,lum,stellar_life_time,mass,frac_mass, h, dens
    nb_lya=0.d0

    pos=pos*3.085d+21
    h=h*3.085d+21
    if(particle_type==2) then
      ix_min=floor((pos(1)-h+box_size/2.)/box_size*array_size)+1
      ix_max=floor((pos(1)+h+box_size/2.)/box_size*array_size)+1
      iy_min=floor((pos(2)-h+box_size/2.)/box_size*array_size)+1
      iy_max=floor((pos(2)+h+box_size/2.)/box_size*array_size)+1
      iz_min=floor((pos(3)-h+box_size/2.)/box_size*array_size)+1
      iz_max=floor((pos(3)+h+box_size/2.)/box_size*array_size)+1
      do ix=ix_min,ix_max
        do iy=iy_min,iy_max
          do iz=iz_min,iz_max
            dx=pos(1)+box_size/2.-(ix-0.5)*box_size/array_size
            dy=pos(2)+box_size/2.-(iy-0.5)*box_size/array_size
            dz=pos(3)+box_size/2.-(iz-0.5)*box_size/array_size
            dist=sqrt(dx**2+ dy**2+ dz**2 )
            in_cell=(abs(dx)<clsize/2.).and.(abs(dy)<clsize/2.).and.  &
                 (abs(dz)<clsize/2.)
            ix_per=ix
            iy_per=iy
            iz_per=iz
            if(dist < h.or.in_cell) then
              if(ix_per < 1) ix_per=ix_per+array_size
              if(ix_per > array_size) ix_per=ix_per-array_size
              if(iy_per < 1) iy_per=iy_per+array_size
              if(iy_per > array_size) iy_per=iy_per-array_size
              if(iz_per < 1) iz_per=iz_per+array_size
              if(iz_per > array_size) iz_per=iz_per-array_size
              if(dist < h) then
                p_weight_sum=p_weight_sum+kernel(dist,real(h,8))*vol_unit
              else
                p_weight_sum=p_weight_sum+kernel(0.d+0,real(h,8))*vol_unit
              endif
            endif
          enddo
        enddo
      enddo
      do iz=iz_min,iz_max
        do iy=iy_min,iy_max
          do ix=ix_min,ix_max
            dx=pos(1)+box_size/2.-(ix-0.5)*box_size/array_size
            dy=pos(2)+box_size/2.-(iy-0.5)*box_size/array_size
            dz=pos(3)+box_size/2.-(iz-0.5)*box_size/array_size
            dist=sqrt(dx**2+ dy**2+ dz**2 )
            in_cell=(abs(dx)<clsize/2.).and.(abs(dy)<clsize/2.).and.  &
                  (abs(dz)<clsize/2.)
            ix_per=ix
            iy_per=iy
            iz_per=iz
            if(dist < h.or.in_cell) then
              if(ix_per < 1) ix_per=ix_per+array_size
              if(ix_per > array_size) ix_per=ix_per-array_size
              if(iy_per < 1) iy_per=iy_per+array_size
              if(iy_per > array_size) iy_per=iy_per-array_size
              if(iz_per < 1) iz_per=iz_per+array_size
              if(iz_per > array_size) iz_per=iz_per-array_size
              if(dist < h) then
                ker=kernel(dist,real(h,8))/p_weight_sum*vol_unit
              else
                ker=kernel(0.d+0,real(h,8))/p_weight_sum*vol_unit
              endif
              cell(ix_per,iy_per,iz_per)%vel(1)=cell(ix_per,iy_per,iz_per)%vel(1)+vel(1)*ker
              cell(ix_per,iy_per,iz_per)%vel(2)=cell(ix_per,iy_per,iz_per)%vel(2)+vel(2)*ker
              cell(ix_per,iy_per,iz_per)%vel(3)=cell(ix_per,iy_per,iz_per)%vel(3)+vel(3)*ker
              cell(ix_per,iy_per,iz_per)%T4=cell(ix_per,iy_per,iz_per)%T4+temp*ker
              cell(ix_per,iy_per,iz_per)%HI_number_density=cell(ix_per,iy_per,iz_per)%HI_number_density+(1.d+0-frac_ion)*ker
              cell(ix_per,iy_per,iz_per)%HII_number_density=cell(ix_per,iy_per,iz_per)%HII_number_density+frac_ion*ker
            ! Photons Lyman-alpha X
              diff_num(ix_per,iy_per,iz_per)=diff_num(ix_per,iy_per,iz_per)+nb_lya*ker
              weight(ix_per,iy_per,iz_per)=weight(ix_per,iy_per,iz_per)+ker
            endif
          enddo
        enddo
      enddo
    endif
    if(lum /=0.) then
      ix=floor((pos(1)+box_size/2.)/box_size*array_size)+1
      iy=floor((pos(2)+box_size/2.)/box_size*array_size)+1
      iz=floor((pos(3)+box_size/2.)/box_size*array_size)+1
      if(ix < 1) ix=ix+array_size
      if(ix > array_size) ix=ix-array_size
      if(iy < 1) iy=iy+array_size
      if(iy > array_size) iy=iy-array_size
      if(iz < 1) iz=iz+array_size
      if(iz > array_size) iz=iz-array_size
! Instantaneous luminosity extended over the whole inter-snapshot period. 
! No interpolation with future luminosity
      cell(ix,iy,iz)%luminosity=cell(ix,iy,iz)%luminosity+real(lum/luminosity_unit)*luminosity_ratio
    endif

    if(particle_type>0 .and. particle_type < 3) then
      compte = compte + 1
    endif

  enddo
enddo

i=compte
close(34)

CALL MPI_ALLREDUCE(int(i,KIND=8),compte_8,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierror)

do j=1,9
  do iz=1,array_size
    do iy=1,array_size
      do ix=1,array_size
        select case(j)
          case(1)
            work_array1(ix,iy)=cell(ix,iy,iz)%vel(1)
          case(2)
            work_array1(ix,iy)=cell(ix,iy,iz)%vel(2)
          case(3)
            work_array1(ix,iy)=cell(ix,iy,iz)%vel(3)
          case(4)
            work_array1(ix,iy)=cell(ix,iy,iz)%T4
          case(5)
            work_array1(ix,iy)=cell(ix,iy,iz)%HI_number_density
          case(6)
            work_array1(ix,iy)=cell(ix,iy,iz)%HII_number_density
          case(7)
             work_array1(ix,iy)=cell(ix,iy,iz)%luminosity
          case(8)
            work_array1(ix,iy)=diff_num(ix,iy,iz)
          case(9)
            work_array1(ix,iy)=weight(ix,iy,iz)
        end select
      enddo
    enddo
    CALL MPI_ALLREDUCE(work_array1,work_array2,array_size*array_size,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
    do iy=1,array_size
      do ix=1,array_size
        select case(j)
          case(1)
            cell(ix,iy,iz)%vel(1)=work_array2(ix,iy)
          case(2)
            cell(ix,iy,iz)%vel(2)=work_array2(ix,iy)
          case(3)
            cell(ix,iy,iz)%vel(3)=work_array2(ix,iy)
          case(4)
            cell(ix,iy,iz)%T4=work_array2(ix,iy)
          case(5)
            cell(ix,iy,iz)%HI_number_density=work_array2(ix,iy)
          case(6)
            cell(ix,iy,iz)%HII_number_density=work_array2(ix,iy)
          case(7)
            cell(ix,iy,iz)%luminosity=work_array2(ix,iy)
          case(8)
            diff_num(ix,iy,iz)=work_array2(ix,iy)
          case(9)
            weight(ix,iy,iz)=work_array2(ix,iy)
        end select
      enddo
    enddo
  enddo
enddo


IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
 print*,compte_8,'particules lues'
END IF

rho_crit=3.*(hubble*100.*1d+5/3.08d+24)**2/8/3.141592/6.67259d-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
part_mass=rho_crit*omega_b*(1-Y_p)*box_size**3/compte_8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'Nombre de cases vides:',count(weight==0)
  print*,'Recomputed particle mass (cgs):', part_mass
END IF


do ix=1,array_size
  do iy=1,array_size
    do iz=1,array_size
      ! Convert velocities to comoving cm/s
      cell(ix,iy,iz)%vel(1)=cell(ix,iy,iz)%vel(1)/weight(ix,iy,iz)*97.75d0*1.d5                 
      cell(ix,iy,iz)%vel(2)=cell(ix,iy,iz)%vel(2)/weight(ix,iy,iz)*97.75d0*1.d5                 
      cell(ix,iy,iz)%vel(3)=cell(ix,iy,iz)%vel(3)/weight(ix,iy,iz)*97.75d0*1.d5                 
      cell(ix,iy,iz)%T4=cell(ix,iy,iz)%T4/weight(ix,iy,iz)/1.d+4
      cell(ix,iy,iz)%HI_number_density=cell(ix,iy,iz)%HI_number_density*(part_mass/proton_mass/clsize**3)
      cell(ix,iy,iz)%HII_number_density=cell(ix,iy,iz)%HII_number_density*(part_mass/proton_mass/clsize**3)
     ! Do not divide by weight! Intensive quantity+ normalized kernel.
     ! In this model 10% of the total luminosity is in the form of X-rays. So 10% of the particle does not radiate in the Lyman band.
     ! Sunghye: 1.927334213369105E+046 pour UV, 5.677334371619719E+045 pour Lya
     ! cell(ix,iy,iz)%luminosity=cell(ix,iy,iz)%luminosity*(0.9*5.677334371619719d+45/luminosity_unit)
!      if(cell(ix,iy,iz)%luminosity.ne.0.0) then
!        write(44,666) num_input,ix,iy,iz,cell(ix,iy,iz)%luminosity
! 666  format (i2,1x,i3,1x,i3,1x,i3,1x,e26.15)
!      end if
    enddo
  enddo
enddo
IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'Average HI number density per commoving cm3', sum(dble(cell(:,:,:)%HI_number_density))/array_size**3
  print*,'Average T4 ', sum(dble(cell(:,:,:)%T4))/array_size**3,maxval(cell(:,:,:)%T4)
END IF



! call cpu_time(t2)
!print*,'Interpolation computed in', t2-t1, 'sec'
 
val=0.0
 call compute_t(tf,val)
UNIV_AGE=tf
IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'UNIV_AGE =', tf/86400./365./1d+9,' Gyr'
END IF
 call compute_t(tnow,exp_fact)
tnow=tf-tnow
write(NUMBER1,'(I3.3)') num_input+1
OPEN(34,FILE='../simu_GD/snapshots/rank00000_out'//trim(adjustl(NUMBER1)),FORM='unformatted',STATUS= 'OLD')
read(34) mxbody_data_in
read(34) exp_fact,numfich
close(34)

a_next=exp_fact 
call compute_t(t_next,a_next)
t_next=tf-t_next
LY_ALPHA_DT=t_next-tnow

call compute_t(t_next,dble(1./11.))
t_next=tf-t_next
call expansion_factor(t_next,a,adot)
Hz10=adot/a

IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'Current time (in Myr):',tnow/86400./365./1d+6
  print*,'Integration time until next output (in Myr):',LY_ALPHA_DT/86400./365./1d+6
END IF

!do ix=1,array_size
!  do iy=1,array_size
!    do iz=1,array_size
!      cell(ix,iy,iz)%luminosity=cell(ix,iy,iz)%luminosity*20.3*1d+6*365*86400/LY_ALPHA_DT
!    enddo
!  enddo
!enddo

deallocate(work_array1)
deallocate(work_array2)
deallocate(weight)

END SUBROUTINE LICORICE_NEW_INPUT


! --------------------------------------------------------------------
! INPARM: read in input parameters
!
! --------------------------------------------------------------------
SUBROUTINE INPARM
implicit none

IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'=============================================================='
  print*,'                 LECTURE FICHIER TREEPARS'
  print*,'=============================================================='
END IF

open(41,file='treepars',STATUS='OLD')
read(41,*) RESTART
read(41,*) FIRST_FILE
read(41,*) LAST_FILE
READ(41,*) BOX_SIZE 
 box_size=box_size*3.085d24/hubble
 CLSIZE=BOX_SIZE/ARRAY_SIZE
READ(41,'(a)') boundary_cond
READ(41,*) luminosity_ratio
READ(41,*) domain_nb_in_input_data
READ(41,*) nb_thread

IF((ierror==MPI_SUCCESS).and.(myrank==MASTER))then
  print*,'using',nb_thread,'threads'
  if(RESTART) then
    print*,'This is a restart'
  else
    print*,'This is not a restart'
  endif
  print*,'First file is number', first_file
  print*,'Last file is number', last_file
  print*,'Box size (comoving cm)',BOX_SIZE
  print*,'Cell size (comoving cm)',CLSIZE
  print*,'Lyman_alpha luminosity / ionizing continuum luminosity', luminosity_ratio 
  print*,'Number of domains in input data', domain_nb_in_input_data
END IF

 close(41)

END SUBROUTINE INPARM

! -----------------------------------------------------------------------
! KERNEL : kernel function : spherically symmetric spline kernel function
!          proposed by Monaghan and Lattanzio (1985).
! -----------------------------------------------------------------------

FUNCTION KERNEL(x,h0)
implicit none
REAL*8 KERNEL 
REAL*8 x,h0
REAL*8 R,denom

R=2.*x/h0
denom=pi*h0**3/8.


IF (R .LE. 1) THEN
  KERNEL = (1.-1.5*R**2+0.75*R**3)/denom
ELSE IF(R .LE. 2) THEN
  KERNEL = (0.25*(2.-R)**3)/denom
ELSE
  KERNEL = 0.
ENDIF
!print*,'gggg',kernel,x,h0,denom,R,pi

RETURN
END FUNCTION



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 29 April 2003
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE convert_big_small_endian_real (realIn,dim )
  IMPLICIT NONE
  integer, intent(in) :: dim
  REAL, dimension(dim),INTENT(INOUT)               :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
  REAL,dimension(dim) :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
  integer :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables (generic 32 bit INTEGER spaces):

  INTEGER                                       :: i_element
  INTEGER                                       :: i_element_br

  do i=1,dim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
  i_element = TRANSFER( realIn(i), 0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reverse order of 4 bytes in 32 bit INTEGER space:
  CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
  CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
  CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
  CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer reversed order bytes to 32 bit REAL space (realOut):
  realOut(i) = TRANSFER( i_element_br, 0.0 )
  realin(i)=realout(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE convert_big_small_endian_double(doublein)
  IMPLICIT none
  
  double precision, intent(inout) :: doublein
  double precision :: doubleout
  INTEGER, dimension(2) :: i_element
  INTEGER, dimension(2) :: i_element_br

  i_element = TRANSFER( doublein, (/0,0/) )

  CALL MVBITS( i_element(2), 24, 8, i_element_br(1), 0  )
  CALL MVBITS( i_element(2), 16, 8, i_element_br(1), 8  )
  CALL MVBITS( i_element(2),  8, 8, i_element_br(1), 16 )
  CALL MVBITS( i_element(2),  0, 8, i_element_br(1), 24 )
  CALL MVBITS( i_element(1), 24, 8, i_element_br(2), 0  )
  CALL MVBITS( i_element(1), 16, 8, i_element_br(2), 8  )
  CALL MVBITS( i_element(1),  8, 8, i_element_br(2), 16 )
  CALL MVBITS( i_element(1),  0, 8, i_element_br(2), 24 )

  doubleout= TRANSFER(i_element_br,0.0d+0)

  doublein=doubleout

  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  SUBROUTINE convert_big_small_endian_integer(intin)
  IMPLICIT NONE
  integer, intent(inout) :: intin
  integer            :: intout 
  CALL MVBITS( intin, 24, 8, intout, 0  )
  CALL MVBITS( intin, 16, 8, intout, 8  )
  CALL MVBITS( intin,  8, 8, intout, 16 )
  CALL MVBITS( intin,  0, 8, intout, 24 )

  intin=intout

  END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUBROUTINE MPI_SPLIT_ALLREDUCE(split_factor,work_array1,work_array2,array_size,MY_MPI_REAL,MY_MPI_SUM,MY_MPI_COMM_WORLD,ierror)

! real(KIND=4), dimension(array_size,array_size,array_size), INTENT(IN) :: work_array1
! real(KIND=4), dimension(array_size,array_size,array_size), INTENT(INOUT) :: work_array1
! integer, intent(in) :: split_fact, array_size, MY_MPI_REAL,MY_MPI_SUM,MY_MPI_COMM_WORLD,ierror

 
! END SUBROUTINE MPI_SPLIT_ALLREDUCE
!
END MODULE IO
