MODULE FUNCTIONS
!USE MKL_VSL_TYPE
!USE MKL_VSL
USE VARS


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random generator from numerical recipe, adapted for OpenMP parallelization.
! Each thread uses a different seed after init_rando is called
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION rando(idum)
USE OMP_LIB
implicit none
INTEGER, INTENT(IN)::idum
INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773
INTEGER, PARAMETER :: IR=2836,NDIV=1+(IM-1)/NTAB
REAL rando
REAL, PARAMETER :: AM=1./IM,EPS=1.2e-7,RNMX=1.-EPS
!Call with idum a negative integer to initialize. RNMX should approxima
!te the largest  oating value that is less than 1.

INTEGER j,k

if (idum .lt. 0 .or. iiy.eq.0) then
  idum_loc=max(-idum,1)
  do j=NTAB+8,1,-1
    k=idum_loc/IQ
    idum_loc=IA*(idum_loc-k*IQ)-IR*k
    if (idum_loc.lt.0) idum_loc=idum_loc+IM
    if (j.le.NTAB) iv(j)=idum_loc
  enddo
  iiy=iv(1)
endif
k=idum_loc/IQ
idum_loc=IA*(idum_loc-k*IQ)-IR*k
if (idum_loc.lt.0) idum_loc=idum_loc+IM
j=1+iiy/NDIV
iiy=iv(j)
iv(j)=idum_loc
rando=min(AM*iiy,RNMX)
return

end FUNCTION rando
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_rando()
USE OMP_LIB
REAL :: u1

!U1=rando(-3049*(omp_get_thread_num()+1))
U1=rando(-3049*(omp_get_thread_num()+myrank*omp_get_num_threads()))
!print*,'Resultat de init_rando :',U1,omp_get_thread_num()+myrank*omp_get_num_threads()

END SUBROUTINE init_rando

!------------------------------------------------------------
! Error function
!-----------------------------------------------------------

FUNCTION erfunc(x)
IMPLICIT NONE
real*8 erfunc
real*8 x
real*8 val
integer i

val=0.5*1*(x/500.)
do i=1,499
  val=val+exp(-(i*x/500.)**2)*(x/500.)
enddo
val=val+0.5*exp(-x*x)*(x/500.)

erfunc=2./sqrt(3.1415925389)*val
return

end FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Approx des fonctions erf et inverse de erf, erreur relative < 10^-4
!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION my_erf(x)
implicit none
real(KIND=8), intent(in) :: x
real(KIND=8) :: a,z

a=0.14d0
z=x*x

my_erf=sqrt(1.d0-exp(-z*(4.d0/pi+a*z)/(1+a*z)))
if(x < 0d0) my_erf=-my_erf
return

END FUNCTION my_erf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! approximate erfc, but exact inverse of my_erf to numerical truncation
! errors (which have been minimized).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION inv_erf(x)
implicit none
real(KIND=8), intent(in) :: x
real(KIND=8) :: a,b,c,z,sol1,sol2

a=0.14d0

z=log(1-x*x)
b=4d0/pi+a*z
c=z
! Solve a*x^4+b*x^2+c=0 using robust method if a*c << b^2
sol1=-1d0*sign(1d0,b)*(dabs(b)+sqrt(b*b-4.d0*a*c))/2d0/a
sol2=c/a/sol1
if(sol1>0d0 .and.sol2>0d0) print*,'NEW INV ERF ISSUE: two positive solutions'
if(sol1 >0d0) then
  inv_erf=sqrt(sol1)
else
  inv_erf=sqrt(sol2)
endif

!if(abs(x) > 1d-6) then
!  z=log(1-x*x)
!  inv_erf=sqrt( -2./pi/a - z/2. + sqrt( (2./pi/a+z/2.)**2 - 1./a*z ) )
!else
!  inv_erf=sqrt(pi)*abs(x)/2.
!endif
if(x < 0d0) inv_erf=-inv_erf

return

END FUNCTION inv_erf

END MODULE FUNCTIONS
