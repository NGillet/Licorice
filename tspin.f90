MODULE TSPIN
USE VARS


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the spin temperature, the Lyman-alpha
! coupling coefficient and the differential
! brightness temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TS

implicit none

integer :: ix,iy,iz,i,cpt_dtb
character(3) :: num1
real(KIND=8) :: nb_sec
real(KIND=8) :: tk,kappa,gammae,nHI,nHII,S_alpha_tilde,xi,xc,t_s,inv_T_c_eff,x_alpha_tilde,x_alpha,t_s_2,del,deltaT_b,deltaT_b_2
real(KIND=8) :: moy_ion,moy_ion2,moy_ion3,moy_dtb,moy_dtb2,moy_xalpha,moy_xalpha2,moy_ts,moy_ts2,moy_tk
REAL(KIND=4), allocatable :: dtb(:,:,:),xa(:,:,:),temp_spin(:,:,:)

!----------------------------------------------------------------------

!allocate(dtb(array_size,array_size,array_size))
allocate(xa(array_size,array_size,array_size))
!allocate(temp_spin(array_size,array_size,array_size))

write(num1,'(I3.3)') num_file
open(37,FILE='ts/ts_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
!open(38,FILE='ts/ts2_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
open(39,FILE='dTb/dTb_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
!open(40,FILE='dTb/dTb2_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
!open(42,FILE='ion_frac/ion_frac_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')
open(43,FILE='ly_alpha/ly_alpha_'//num1//'.dat',FORM='UNFORMATTED',STATUS='REPLACE')

print*,'Computes T_s and dT_b for file ',num_file

moy_ion       = 0.0
moy_ion2      = 0.0
moy_ion3      = 0.0
moy_dtb       = 0.0
moy_dtb2      = 0.0
moy_xalpha    = 0.0
moy_xalpha2   = 0.0
moy_ts        = 0.0
moy_ts2       = 0.0
moy_tk        = 0.0

cpt_dtb = 0

!dtb = 0.0
xa  = 0.0

nb_sec        = LY_ALPHA_DT

do ix=1,array_size
  do iy=1,array_size
    do iz=1,array_size
!======================================================================
!             Collisional coupling
!======================================================================
       tk     = cell(ix,iy,iz)%T4*1.0d4
       nHI    = cell(ix,iy,iz)%HI_number_density*exp_fact**(-3.0)                                                    ! comoving -> physical
       nHII   = cell(ix,iy,iz)%HII_number_density*exp_fact**(-3.0)
       kappa  = 3.1d-11*tk**(0.357)*exp(-32./tk)                                                                     ! Collisions with HI and with protons
       if (tk.le.1.0d4) then                                                                                         ! Collisions with electrons
         gammae = 10.0**(-9.607+0.5*log10(max(tk,1.0d0))*exp(-log10(max(tk,1.0d0))**(4.5)/1800.0))
       else
         gammae = 7.9043d-9
       endif
       xc            = 0.0682/2.85d-15/2.725*exp_fact*(nHI*kappa+nHII*(gammae+3.2*kappa))
!======================================================================
!             S_alpha_tilde
!======================================================================
!       x_alpha       = diff_num_MPI(ix,iy,iz)/(cell(ix,iy,iz)%HI_number_density*(clsize)**3)/nb_sec*4./27.*0.068/2.73*exp_fact/2.85d-15
       x_alpha       = diff_num(ix,iy,iz)/(cell(ix,iy,iz)%HI_number_density*(clsize)**3)/nb_sec*4./27.*0.068/2.73*exp_fact/2.85d-15

       xi            = (1.0d-7*3.0*nHI*1.215703398d-5**3.0*50.0d6/2.0/exp_fact_dot*exp_fact)**(1./3.)*tk**(-2./3.)   ! equations 41 and 35 Hirata
       t_s           = 2.725/exp_fact                                                                                ! initialization T_s = T_gamma
       do i=1,6                                                                                                      ! converges in a few steps
          S_alpha_tilde = (1.0-0.0631789/tk+0.115995/tk/tk-0.401403/t_s/tk+0.336463/t_s/tk/tk) &                      ! equation 40 Hirata
                        /(1.0+2.98394*xi+1.53583*xi*xi+3.85289*xi*xi*xi)
!======================================================================
!             Color temperature
!======================================================================
          inv_T_c_eff   = 1./tk+0.405535/tk*(1./t_s-1./tk)
!======================================================================
!             Spin temperature
!======================================================================
          x_alpha_tilde = x_alpha*S_alpha_tilde
          t_s           = (1.0+x_alpha_tilde+xc)/(exp_fact/2.725+x_alpha_tilde*inv_T_c_eff+xc/tk)
       enddo
!       temp_spin(ix,iy,iz) = t_s
!       write(37) real(t_s,8)
!       write(43) x_alpha_tilde*(inv_T_c_eff-1./t_s)/(1./tk-1./t_s)
       xa(ix,iy,iz) = real(x_alpha_tilde*(inv_T_c_eff-1./t_s)/(1./tk-1./t_s),4)
!======================================================================
!             Spin temperature without correction (for comparison)
!======================================================================
         t_s_2         = (1.0+x_alpha+xc)/(exp_fact/2.725+x_alpha/tk+xc/tk)
!       write(38) real(t_s_2,8)
!======================================================================
!             Computes the differential brightness temperature
!======================================================================
! Equation (8) Baek et al. Does not take PROPER VELOCITY GRADIENT into account
       del             = (nHI+nHII)/(3./8./3.141592/6.67259d-8*omega_b*(1-Y_p)*exp_fact**(-3.)*(hubble*100.*1.0d5/3.0856d24)**(2.)/proton_mass)
       deltaT_b        = 28.1*nHI/(nHI+nHII)*del*(1.0/10.0/exp_fact)**(0.5)*(1.0-2.725/exp_fact/t_s) &
                         *(omega_b*hubble/0.042/0.73)*(0.24/omega_0)**(0.5)*(1.0-Y_p)/0.752
       if (deltaT_b.gt.-500.0.or.deltaT_b.lt.500.0)then
       else
         print*,'deltaT_b=',deltaT_b,'t_s=',t_s,'x_alpha_tilde=',x_alpha_tilde,'inv_T_c_eff=',inv_T_c_eff
         print*,'xc=',xc,'kappa=',kappa,'gammae=',gammae,'tk',tk
       end if
!       write(39) real(deltaT_b,8)
!       dtb(ix,iy,iz) = real(deltaT_b,4)
       deltaT_b_2      = 28.1*nHI/(nHI+nHII)*del*(1.0/10.0/exp_fact)**(0.5)*(1.0-2.725/exp_fact/t_s_2)
!       write(40) real(deltaT_b_2,8)
!======================================================================
!             Computes the ionization fraction
!======================================================================
!       write(42) real(nHII/(nHI+nHII),8)
!======================================================================
!             Computes averages
!======================================================================
       moy_ion         = moy_ion + nHII
       moy_ion2        = moy_ion2 + nHI + nHII
       moy_ion3        = moy_ion3 + nHII/(nHI+nHII)
       moy_dtb         = moy_dtb + deltaT_b
       moy_dtb2        = moy_dtb + deltaT_b_2
       moy_xalpha      = moy_xalpha + x_alpha_tilde*(inv_T_c_eff-1./t_s)/(1./tk-1./t_s)
       moy_xalpha2     = moy_xalpha2 + x_alpha
       moy_ts          = moy_ts + t_s
       moy_ts2         = moy_ts2 + t_s_2
       moy_tk          = moy_tk + tk

       cpt_dtb = cpt_dtb+1
       if(cpt_dtb .eq. 0.25*array_size**3) print*,'Computing dTb: 1/4',myrank
       if(cpt_dtb .eq. 0.5*array_size**3) print*,'Computing dTb: 1/2',myrank
       if(cpt_dtb .eq. 0.75*array_size**3) print*,'Computing dTb: 3/4',myrank
       if(cpt_dtb .eq. array_size**3) print*,'Computing dTb: FINI',myrank
    enddo
  enddo
enddo

print*,'Writing the files...',myrank
!write(37) temp_spin
!write(39) dtb
write(43) xa

print*,'Average differential brightness temperature: ', moy_dtb/array_size**3
print*,'Mass averaged ionization fraction: ', moy_ion/moy_ion2
print*,'Volume averaged ionization fraction: ', moy_ion3/array_size**3
print*,'Average x_alpha: ', moy_xalpha/array_size**3
print*,'Average T_s: ',moy_ts/array_size**3
print*,'Average T_k: ',moy_tk/array_size**3
print*,'Domain',myrank,'Minimum x_alpha:',minval(xa)
print*,'Domain',myrank,'Maximum x_alpha:',maxval(xa)

write(50,119) 1.0/exp_fact-1.0,moy_dtb/array_size**3,moy_dtb2/array_size**3,moy_xalpha/array_size**3,moy_ts/array_size**3,moy_tk/array_size**3,(1.0/exp_fact)*2.725
 119  format (11(1x,e26.15))
!close(37)
!close(38)
close(39)
!close(40)
!close(42)
close(43)

!deallocate(dtb)
deallocate(xa)
!deallocate(temp_spin)

END SUBROUTINE TS


END MODULE TSPIN
