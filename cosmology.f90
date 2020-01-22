MODULE COSMOLOGY

USE VARS
USE OMP_LIB


integer, PARAMETER :: exp_fact_ntab=1000
real(KIND=8), SAVE, DIMENSION(0:exp_fact_ntab*2) :: a_tab,adot_tab

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute a et a dot at time t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine expansion_factor(t,a,adot) 
implicit none
real (KIND=8), intent(in) :: t
real(KIND=8), intent(out) :: a,adot

integer :: i
real(KIND=8) :: dt,a8,t8,h_0,da1,da2,da3,da4

h_0=hubble*100.*1d+5/3.08d+24
a8=a
t8=t

dt=(UNIV_AGE-t8)/1000000.
a8=1.
do i=1,1000000 
  da1=dt*a8*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a8**3)
  da2=dt*(a8+da1/2.) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a8+da1/2.)**3)
  da3=dt*(a8+da2/2.) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a8+da2/2.)**3)
  da4=dt*(a8+da3) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a8+da3/2.)**3)

  a8=a8-(da1/6.+da2/3.+da3/3.+da4/6.)
enddo

adot=a8*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a8**3)
a=a8


end subroutine expansion_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! increment a and adot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine increment_exp_fact(t_old,a_old,t_new,a_new,adot_new)
implicit none
real (KIND=8), intent(inout) :: t_old,a_old,t_new
real(KIND=8), intent(out) :: a_new,adot_new

integer :: i
real(KIND=8) :: dt,a08,t08,h_0,da01,da02,da03,da04,h_0_sq,a_sav

a_sav=a_old

h_0=hubble*100.*1d+5/3.08d+24
h_0_sq=H_0**2
a08=a_old
t08=t_old
i=0

dt=0.01
do while(t08 < t_new)
  da01=dt*a08*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/a08**3)
  da02=dt*(a08+da01/2.)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da01/2.)**3)
  da03=dt*(a08+da02/2.)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da02/2.)**3)
  da04=dt*(a08+da03)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da03/2.)**3)

  a08=a08+(da01/6.+da02/3.+da03/3.+da04/6.)
  t08=t08+dt
  i=i+1

!  if(mod(i,1000)==0) print*,'increment_exp_fact ALARM',i,t_old,a_old,t_new,t08,a08
enddo

a_new=a08-(da01/6.+da02/3.+da03/3.+da04/6.)*(t08-t_new)/dt
adot_new=a_new*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a_new**3)

end subroutine increment_exp_fact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute time form expansion factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_t(t1,a1)
implicit none
real(KIND=8) :: t1,a1
integer :: i
real(KIND=8) :: a,t,da1,da2,da3,da4,dt,h_0

h_0=hubble*100.*1d+5/3.08d+24

t=0.
dt=86400.*365.*10000.
a=1.
do while(a.gt.a1)
  da1=dt*a*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a**3)
  da2=dt*(a+da1/2.) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a+da1/2.)**3)
  da3=dt*(a+da2/2.) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a+da2/2.)**3)
  da4=dt*(a+da3) &
     *sqrt(H_0**2*omega_lambda+H_0**2*omega_0/(a+da3)**3)
  a=a-(da1/6.+da2/3.+da3/3.+da4/6.)
  t=t+dt
enddo
t1=t

end subroutine compute_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tabulate_exp_fact(t_ini,a_ini,t_end)

real(KIND=8), intent(IN):: t_ini,a_ini,t_end
integer :: i
real(KIND=8) :: dt,a08,h_0,da01,da02,da03,da04,h_0_sq


h_0=hubble*100.*1d+5/3.08d+24
h_0_sq=H_0**2
a_tab(0)=a_ini
adot_tab(0)=a_ini*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a_ini**3)
a08=a_ini

dt=(t_end-t_ini)/exp_fact_ntab

if(dt > 86400.*365*1d+6) then
  print*,'Inacurate computation of a_tab'
  STOP
endif
do i=1,exp_fact_ntab*2 
  da01=dt*a08*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/a08**3)
  da02=dt*(a08+da01/2.)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da01/2.)**3)
  da03=dt*(a08+da02/2.)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da02/2.)**3)
  da04=dt*(a08+da03)*sqrt(H_0_sq*omega_lambda+H_0_sq*omega_0/(a08+da03/2.)**3)

  a08=a08+(da01/6.+da02/3.+da03/3.+da04/6.)
  a_tab(i)=a08
  adot_tab(i)=a08*sqrt(H_0**2*omega_lambda+H_0**2*omega_0/a08**3)
enddo

end subroutine tabulate_exp_fact

END MODULE COSMOLOGY
