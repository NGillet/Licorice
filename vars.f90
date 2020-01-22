MODULE VARS
USE MPI

SAVE

! General parameters !
REAL(KIND=8),    PARAMETER :: PI=3.141592653589, CPULIMIT=500000
REAL(KIND=8),    PARAMETER :: inv_sqrt_pi=0.564189583547
REAL(KIND=8),    PARAMETER :: proton_mass=1.6726d-24 !cgs
REAL(KIND=8),    PARAMETER :: planck_cgs=6.62d-27
REAL(KIND=8),    PARAMETER :: luminosity_unit=1.d+45  !(To use single precision in the cell structure)

! Higher-order Lyman-series lines parameters
!REAL(KIND=8), dimension(23):: frec = (/0.0, 1.0, 0.0, 0.2609, 0.3078, 0.3259, 0.3353, 0.3410, 0.3448, 0.3476, 0.3496, 0.3512, 0.3524, 0.3535, 0.3543, &
!                                       0.3550, 0.3556, 0.3561, 0.3565, 0.3569, 0.3572, 0.3575, 0.3578/)        ! Pritchard & Furlanetto
REAL(KIND=8), dimension(10):: frec = (/0.0, 1.0, 0.0, 0.1305, 0.1539, 0.1630, 0.1677, 0.1705, 0.1724, 0.1738/)  ! /2 for n.gt.2 (injected photons)
real(KIND=8), dimension(10):: osc_strength = (/0.0d0, 4.162d-1, 7.9d-2, 2.9d-2, 1.39d-2, 7.8d-3, 4.82d-3, 3.19d-3, 2.22d-3, 1.6d-3/)
real(KIND=8), dimension(10):: natural = (/1.0d0, 9.971d7, 3.018d7, 1.288d7, 6.636d6, 3.857d6, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
real(KIND=8), dimension(10):: prob_casc = (/0.0, 0.0, 0.1183, 0.1610, 0.1822, 0.1947, 0.2028, 0.0, 0.0, 0.0/)

! Various array sizes !
INTEGER, PARAMETER :: ARRAY_SIZE=1024

! Cosmological variables
!real(KIND=8), parameter :: hubble=0.678,omega_0=0.308,omega_lambda=0.692,omega_b=0.0484,Y_p=0.0 !0.248
real(KIND=8), parameter :: hubble=0.678,omega_0=0.3175,omega_lambda=0.6825,omega_b=0.049,Y_p=0.0
REAL(KIND=8) :: UNIV_AGE
real(KIND=8) :: EXP_FACT,EXP_FACT_DOT
real(KIND=8) :: exp_fact_thread,exp_fact_dot_thread,inv_exp_fact_thread,tnow

! C2RAY parameter
INTEGER, PARAMETER :: cubeP3M_array_size=6144
REAL(KIND=8), PARAMETER :: C2ray_source_blackbody_temp=50000


!$OMP THREADPRIVATE(exp_fact_thread,exp_fact_dot_thread,inv_exp_fact_thread)

! Input parameters
CHARACTER(len=8)  :: BOUNDARY_COND
LOGICAL  :: RESTART
REAL(KIND=8) :: BOX_SIZE
REAL(KIND=8) :: luminosity_ratio
INTEGER:: domain_nb_in_input_data
INTEGER:: NB_THREAD

! Random generator variable
INTEGER :: ISEED=0
INTEGER, SAVE :: idum_loc,iiy
INTEGER, PARAMETER :: NTAB=32
INTEGER, DIMENSION(NTAB),SAVE :: iv

!$OMP THREADPRIVATE(idum_loc,iiy,iv)

INTEGER :: NUM_FILE,FIRST_FILE,LAST_FILE

! Radiative transfer variables
TYPE :: CELL_STRUCT
!  SEQUENCE
  REAL(KIND=4),dimension(3) :: vel
  REAL(KIND=4) :: HI_number_density
  REAL(KIND=4) :: HII_number_density
  REAL(KIND=4) :: luminosity
  REAL(KIND=4) :: T4
END TYPE

TYPE (CELL_STRUCT), dimension(ARRAY_SIZE,ARRAY_SIZE,ARRAY_SIZE) :: CELL

REAL(KIND=8) :: CLSIZE
INTEGER, PARAMETER :: NB_PHOT_PER_SNAPSHOT = 396000000                !13000000             !Valeur initiale du fichier de B. Semelin: 1600000000
INTEGER, PARAMETER :: NB_PHOT_PER_BATCH_PER_THREAD = 500000         !Valeur initiale du fichier de B. Semelin: 5000000
INTEGER, PARAMETER :: OLD_PHOT_ARRAY_SIZE=NB_PHOT_PER_BATCH_PER_THREAD
REAL(KIND=8) :: NB_PHOT_PER_SOURCE
REAL(KIND=8) ::LY_ALPHA_DT

REAL(KIND=8) :: PHOT_FREQ
REAL(kind=4), DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD,3) :: OLD_PHOT_POS_LT, OLD_PHOT_DIR_LT
REAL(kind=4), DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: OLD_PHOT_FREQ_LT,old_phot_target_opt_depth_lt
REAL(kind=4), DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: OLD_PHOT_time_lt
REAL(kind=8), DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: nb_physical_phot_per_old_phot_lt
REAL(kind=8), DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: diff_num_lt
INTEGER     , DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: diff_pos_lt
LOGICAL, DIMENSION(NB_PHOT_PER_BATCH_PER_THREAD) :: OLD_PHOT_ARRAY_IN_USE_LT
REAL(KIND=8) :: diff_num(ARRAY_SIZE,ARRAY_SIZE,ARRAY_SIZE)
INTEGER :: nb_inflight_phot_lt

!$OMP THREADPRIVATE(PHOT_FREQ,OLD_PHOT_ARRAY_IN_USE_LT,OLD_PHOT_time_lt,nb_physical_phot_per_old_phot_lt)
!$OMP THREADPRIVATE(OLD_PHOT_POS_LT,OLD_PHOT_DIR_LT,OLD_PHOT_FREQ_LT,old_phot_target_opt_depth_lt,nb_inflight_phot_lt)

REAL(KIND=8) :: T_INIT_BATCH,Hz10

! MPI variables
INTEGER :: ierror,size,myrank,somme
INTEGER :: PROVIDED
INTEGER, PARAMETER::MASTER=0
!REAL(KIND=8) :: diff_num_MPI(ARRAY_SIZE,ARRAY_SIZE,ARRAY_SIZE)


END MODULE VARS
