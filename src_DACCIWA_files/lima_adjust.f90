!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_LIMA_ADJUST
!     #######################
!
INTERFACE
!
      SUBROUTINE LIMA_ADJUST(KRR, KMI, TPFILE, HRAD,                           &
                             HTURBDIM, OCLOSE_OUT, OSUBG_COND, PTSTEP,         &
                             PRHODREF, PRHODJ, PEXNREF, PPABSM, PSIGS, PPABST, &
                             PRT, PRS, PSVT, PSVS,                             &
                             PTHS, PSRCS, PCLDFR, PTHT, PTHM, KTCOUNT          )

USE MODD_IO_ll, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
CHARACTER*4,              INTENT(IN)   :: HTURBDIM   ! Dimensionality of the
                                                     ! turbulence scheme
CHARACTER*4,              INTENT(IN)   :: HRAD       ! Radiation scheme name
LOGICAL,                  INTENT(IN)   :: OCLOSE_OUT ! Conditional closure of 
                                                     ! the OUTPUT FM-file
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid 
                                                     ! Condensation
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSM    ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT      ! Concentrations at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS      ! Concentration source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS      ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                     ! s'rc'/2Sigma_s2 at time t+1
                                                     ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT      ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM      ! Theta at time t-Dt
INTEGER,                  INTENT(IN)    :: KTCOUNT   ! Temporal loop counter
!
END SUBROUTINE LIMA_ADJUST
!
END INTERFACE
!
END MODULE MODI_LIMA_ADJUST
!
!     ##########################################################################
      SUBROUTINE LIMA_ADJUST(KRR, KMI, TPFILE, HRAD,                           &
                             HTURBDIM, OCLOSE_OUT, OSUBG_COND, PTSTEP,         &
                             PRHODREF, PRHODJ, PEXNREF, PPABSM, PSIGS, PPABST, &
                             PRT, PRS, PSVT, PSVS,                             &
                             PTHS, PSRCS, PCLDFR, PTHT, PTHM, KTCOUNT          )
!     ##########################################################################
!
!!****  *MIMA_ADJUST* -  compute the fast microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the fast microphysical sources
!!      through an explict scheme and a saturation ajustement procedure.
!!
!!
!!**  METHOD
!!    ------
!!      Reisin et al.,    1996 for the explicit scheme when ice is present
!!      Langlois, Tellus, 1973 for the implict adjustment for the cloud water
!!      (refer also to book 1 of the documentation).
!!
!!      Computations are done separately for three cases :
!!        - ri>0 and rc=0
!!        - rc>0 and ri=0
!!        - ri>0 and rc>0
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!      Module MODD_BUDGET:
!!         NBUMOD 
!!         CBUTYPE
!!         NBUPROCCTR 
!!         LBU_RTH    
!!         LBU_RRV  
!!         LBU_RRC  
!!      Module MODD_LES : NCTR_LES,LTURB_LES,NMODNBR_LES
!!                        XNA declaration (cloud fraction as global var)
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 1 and Book2 of documentation ( routine FAST_TERMS )
!!      Langlois, Tellus, 1973
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      JP Chaboureau *LA*   March 2014  fix the calculation of icy cloud fraction
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      M. Claeys  * LACy*   june 2019  add several shapes for ice crystals
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_CONF
USE MODD_CST
USE MODD_IO_ll,   ONLY: TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM
!
USE MODE_FIELD, ONLY : TFIELDDATA, TYPEREAL
USE MODE_FM
USE MODE_FMWRIT
!
USE MODI_BUDGET
USE MODI_CONDENS
USE MODI_LIMA_FUNCTIONS
USE MODI_COMPUTE_LBDA_SHAPE
USE MODI_LIMA_CHANGE_SHAPE
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
CHARACTER*4,              INTENT(IN)   :: HTURBDIM   ! Dimensionality of the
                                                     ! turbulence scheme
CHARACTER*4,              INTENT(IN)   :: HRAD       ! Radiation scheme name
LOGICAL,                  INTENT(IN)   :: OCLOSE_OUT ! Conditional closure of 
                                                     ! the OUTPUT FM-file
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid 
                                                     ! Condensation
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSM    ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT      ! Concentrations at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS      ! Concentration source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS      ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                     ! s'rc'/2Sigma_s2 at time t+1
                                                     ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT      ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM      ! Theta at time t-Dt
INTEGER,                  INTENT(IN)    :: KTCOUNT   ! Temporal loop counter
!
!
!*       0.2   Declarations of local variables :
!
! 3D Microphysical variables
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: ZRVT,        & ! Water vapor m.r. at t
                            ZRCT,        & ! Cloud water m.r. at t
                            ZRRT,        & ! Rain water m.r. at t
                            ZRIT,        & ! Cloud ice  m.r. at t
                            ZRST,        & ! Aggregate  m.r. at t
                            ZRGT,        & ! Graupel    m.r. at t
!
                            ZRVS,        & ! Water vapor m.r. source
                            ZRCS,        & ! Cloud water m.r. source
                            ZRRS,        & ! Rain water m.r. source
                            ZRIS,        & ! Cloud ice  m.r. source
                            ZRSS,        & ! Aggregate  m.r. source
                            ZRGS,        & ! Graupel    m.r. source
!
                            ZCCT,        & ! Cloud water conc. at t
!
                            ZCCS,        & ! Cloud water C. source
                            ZMAS           ! Mass of scavenged AP

REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                          :: ZCIT_TOT,   & ! Cloud ice   conc. at t
                             ZCIS_TOT      ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE &
                         :: ZNFS,        & ! Free      CCN C. source
                            ZNAS,        & ! Activated CCN C. source
                            ZIFS,        & ! Free      IFN C. source 
                            ZINS,        & ! Nucleated IFN C. source
                            ZNIS           ! Acti. IMM. nuclei C. source
!
!
!
REAL                     :: ZEPS         ! Mv/Md
REAL                     :: ZDT          ! Time increment (2*Delta t or Delta t if cold start)
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: ZEXNS,&      ! guess of the Exner function at t+1
                            ZT,   &      ! guess of the temperature at t+1
                            ZCPH, &      ! guess of the CPh for the mixing
                            ZW,   &
                            ZW1,  &
                            ZW2,  &
                            ZLV,  &      ! guess of the Lv at t+1
                            ZLS,  &      ! guess of the Ls at t+1
                            ZMASK, &
                            ZT0,  &      ! Temperature at t
                            ZTM          ! Temperature at t-1
LOGICAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: GMICRO, GMICRO_RI, GMICRO_RC ! Test where to compute cond/dep proc.
INTEGER                  :: IMICRO
REAL, DIMENSION(:), ALLOCATABLE &
                         :: ZRCT_V, ZRIT_V, ZRVS_V, ZRCS_V, ZRIS_V, ZTHS_V,              &
                            ZCCT_V, ZCCS_V,                                      &
                            ZRHODREF, ZZT, ZPRES, ZEXNREF, ZZCPH,            &
                            ZZW, ZLVFACT, ZLSFACT,                           &
                            ZRVSATW, ZRVSATI, ZRVSATW_PRIME, ZRVSATI_PRIME,  &
                            ZAW, ZAI, ZCJ, ZKA, ZDV, ZITW, ZITI, ZAWW, ZAIW, &
                            ZAWI, ZAII, ZFACT, ZDELTW,                       &
                            ZDELTI, ZDELT1, ZDELT2, ZCND, ZDEP
!
REAL, DIMENSION(:), ALLOCATABLE &
                         :: ZZT0, ZZTM   ! Temperature at t and t-1
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIT_TOT_V, ZCIS_TOT_V
REAL, DIMENSION(:),   ALLOCATABLE :: ZLBDAI ! slope parameter for the ice crystal distribution
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLBDAI_SHAPE_S, ZLBDAI_SHAPE_T
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCSHAPE_T_V
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCSHAPE_S_V
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCIS_SHAPE_V, ZCIT_SHAPE_V
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZONEOVER_VARS, ZONEOVER_VART  ! for optimization
! contient les concentrations par forme de cristaux
REAL, DIMENSION(:,:,:,:),  ALLOCATABLE :: ZCIS_SHAPE, &
                                          ZCIT_SHAPE
REAL, DIMENSION(:,:),      ALLOCATABLE :: ZRIS_SHAPE_V,  &
                                          ZRIT_SHAPE_V
! Rapport de concentration relative de chaque forme de cristaux
REAL, DIMENSION(:,:,:,:),  ALLOCATABLE :: ZCSHAPE_T, ZCSHAPE_S
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW_2D,   &
                                     ZITI_2D,  &
                                     ZDEP_2D,  & 
                                     ZCND_2D 
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZW_2D   
INTEGER, DIMENSION(:),  ALLOCATABLE :: ISHAPE_MAX_S, ISHAPE_MAX_T ! index of the dominant shape concentration 
!
INTEGER                  :: IRESP      ! Return code of FM routines
INTEGER                  :: IKB        ! K index value of the first inner mass point
INTEGER                  :: IKE        ! K index value of the last inner mass point
INTEGER                  :: IIB,IJB    ! Horz index values of the first inner mass points
INTEGER                  :: IIE,IJE    ! Horz index values of the last inner mass points
INTEGER                  :: JITER,ITERMAX  ! iterative loop for first order adjustment
INTEGER                  :: ILUOUT     ! Logical unit of output listing 
!
INTEGER                           :: ISIZE
REAL, DIMENSION(:), ALLOCATABLE   :: ZRTMIN
REAL, DIMENSION(:), ALLOCATABLE   :: ZCTMIN
!
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: JMOD, JMOD_IFN, JMOD_IMM
INTEGER                           ::  JSH     ! index for ice crystal shapes
!
INTEGER , DIMENSION(3) :: BV
TYPE(TFIELDDATA)  :: TZFIELD
!
CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: ZSHAPE_NAME_M, & ! Name of the shape corresponding to the
                                               ZSHAPE_NAME_0    ! temperature growth regime, at t-1 and t
CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: ZSHAPE_NAME_M_V, &
                                               ZSHAPE_NAME_0_V
LOGICAL, DIMENSION(:),   ALLOCATABLE :: GHABIT ! Test where to compute change of ice crystal habit
INTEGER                              :: IHABIT
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZCIS_SHAPE_VV
REAL,    DIMENSION(:),   ALLOCATABLE :: ZZTM_V, ZZT0_V
REAL,    DIMENSION(:),   ALLOCATABLE :: ZDEP_V
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZW_V
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ILUOUT = TLUOUT%NLU
!
IIB = 1 + JPHEXT
IIE = SIZE(PRHODJ,1) - JPHEXT
IJB = 1 + JPHEXT
IJE = SIZE(PRHODJ,2) - JPHEXT
IKB = 1 + JPVEXT
IKE = SIZE(PRHODJ,3) - JPVEXT
!
ZEPS= XMV / XMD
!
IF (OSUBG_COND) THEN
  ITERMAX=2
ELSE
  ITERMAX=1
END IF
!
ZDT = PTSTEP
!
ISIZE = SIZE(XRTMIN)
ALLOCATE(ZRTMIN(ISIZE))
ZRTMIN(:) = XRTMIN(:) / ZDT
ISIZE = SIZE(XCTMIN)
ALLOCATE(ZCTMIN(ISIZE))
ZCTMIN(:) = XCTMIN(:) / ZDT
!
! Prepare 3D water mixing ratios
ZRVT(:,:,:) = PRT(:,:,:,1)
ZRVS(:,:,:) = PRS(:,:,:,1)
!
ZRCT(:,:,:) = 0.
ZRCS(:,:,:) = 0.
ZRRT(:,:,:) = 0.
ZRRS(:,:,:) = 0.
ZRIT(:,:,:) = 0.
ZRIS(:,:,:) = 0.
ZRST(:,:,:) = 0.
ZRSS(:,:,:) = 0.
ZRGT(:,:,:) = 0.
ZRGS(:,:,:) = 0.
!
IF ( KRR .GE. 2 ) ZRCT(:,:,:) = PRT(:,:,:,2)
IF ( KRR .GE. 2 ) ZRCS(:,:,:) = PRS(:,:,:,2)
IF ( KRR .GE. 3 ) ZRRT(:,:,:) = PRT(:,:,:,3) 
IF ( KRR .GE. 3 ) ZRRS(:,:,:) = PRS(:,:,:,3)
IF ( KRR .GE. 4 ) ZRIT(:,:,:) = PRT(:,:,:,4)
IF ( KRR .GE. 4 ) ZRIS(:,:,:) = PRS(:,:,:,4) 
IF ( KRR .GE. 5 ) ZRST(:,:,:) = PRT(:,:,:,5) 
IF ( KRR .GE. 5 ) ZRSS(:,:,:) = PRS(:,:,:,5) 
IF ( KRR .GE. 6 ) ZRGT(:,:,:) = PRT(:,:,:,6)
IF ( KRR .GE. 6 ) ZRGS(:,:,:) = PRS(:,:,:,6)
!
! Prepare 3D number concentrations
ZCCT(:,:,:) = 0.
ZCCS(:,:,:) = 0.
ZCIT_TOT(:,:,:) = 0.
ZCIS_TOT(:,:,:) = 0.
IF (LCOLD .AND. LCRYSTAL_SHAPE) THEN
  ALLOCATE(ZCIT_SHAPE(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NB_CRYSTAL_SHAPE))
  ALLOCATE(ZCIS_SHAPE(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NB_CRYSTAL_SHAPE))
  ZCIT_SHAPE(:,:,:,:) = 0.
  ZCIS_SHAPE(:,:,:,:) = 0.
END IF
!
IF ( LWARM ) THEN
  ZCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC)
  ZCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
END IF
!
IF (LCOLD) THEN
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    ZCIT_TOT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
    ZCIS_TOT(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
  ELSE
    DO JSH = 1, NB_CRYSTAL_SHAPE
! concentration pour chaque forme
      ZCIT_SHAPE(:,:,:,JSH) = PSVT(:,:,:,NSV_LIMA_NI+JSH-1)
      ZCIS_SHAPE(:,:,:,JSH) = PSVS(:,:,:,NSV_LIMA_NI+JSH-1)
    END DO
! concentration totale
    ZCIS_TOT(:,:,:) = SUM(ZCIS_SHAPE, DIM=4)
    ZCIT_TOT(:,:,:) = SUM(ZCIT_SHAPE, DIM=4)
! calcul du rapport de concentration de chaque forme pour ZCIT et ZCIS
    ALLOCATE(ZCSHAPE_S(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NB_CRYSTAL_SHAPE))
    ALLOCATE(ZCSHAPE_T(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NB_CRYSTAL_SHAPE))
    ZCSHAPE_T (:,:,:,:) = 0.0
    ZCSHAPE_S (:,:,:,:) = 0.0
    ALLOCATE(ZONEOVER_VARS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)))
    ALLOCATE(ZONEOVER_VART(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)))
    ZONEOVER_VARS(:,:,:) = 0.
    ZONEOVER_VART(:,:,:) = 0.
    WHERE (ZCIS_TOT(:,:,:) .GT. 0.0) ZONEOVER_VARS(:,:,:) = 1.0 / ZCIS_TOT(:,:,:)
    WHERE (ZCIT_TOT(:,:,:) .GT. 0.0) ZONEOVER_VART(:,:,:) = 1.0 / ZCIT_TOT(:,:,:)
    DO JSH = 1, NB_CRYSTAL_SHAPE
      WHERE ((ZCIS_TOT(:,:,:) .GT. 0.0) .AND. (ZCIS_SHAPE(:,:,:,JSH) .GT. 0.0))
        ZCSHAPE_S(:,:,:,JSH) = MIN(ZCIS_SHAPE(:,:,:,JSH)*ZONEOVER_VARS(:,:,:),1.0)
      END WHERE
      WHERE ((ZCIT_TOT(:,:,:) .GT. 0.0) .AND. (ZCIT_SHAPE(:,:,:,JSH) .GT. 0.0))
        ZCSHAPE_T(:,:,:,JSH) = MIN(ZCIT_SHAPE(:,:,:,JSH)*ZONEOVER_VART(:,:,:),1.0)
      END WHERE
    END DO ! JSH
    DEALLOCATE(ZONEOVER_VARS)
    DEALLOCATE(ZONEOVER_VART)
  END IF ! LCRYSTAL_SHAPE
END IF ! LCOLD
!
!ALLOCATE(ZAUX(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NB_CRYSTAL_SHAPE))
!ZAUX(:,:,:,:) = ZCIS_SHAPE(:,:,:,:)
!
IF ( LSCAV .AND. LAERO_MASS ) ZMAS(:,:,:) = PSVS(:,:,:,NSV_LIMA_SCAVMASS)
! 
IF ( LWARM .AND. NMOD_CCN.GE.1 ) THEN
   ALLOCATE( ZNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( ZNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ZNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   ZNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
END IF
!
IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
   ALLOCATE( ZIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ALLOCATE( ZINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ZIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
   ZINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   ALLOCATE( ZNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
   ZNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
END IF
!
! Compute the temperature at t-1 et t
ZT0(:,:,:) = 0.
ZTM(:,:,:) = 0.
!

ZT0(:,:,:) = PTHT(:,:,:) * (PPABST(:,:,:) / XP00)**(XRD/XCPD)
!IF (KTCOUNT .GT. 1) THEN
!  ZTM(:,:,:) = PTHM(:,:,:) * (PPABSM(:,:,:) / XP00)**(XRD/XCPD)
!ELSE 
  ZTM(:,:,:) = ZT0(:,:,:)
!END IF
!
!-------------------------------------------------------------------------------
!
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!*       2.1    remove negative non-precipitating negative water
!               ------------------------------------------------
!
IF (ANY(ZRVS(:,:,:)+ZRCS(:,:,:)+ZRIS(:,:,:) < 0.) .AND. NVERB>5) THEN
  WRITE(ILUOUT,*) 'LIMA_ADJUST:  negative values of total water (reset to zero)'
  WRITE(ILUOUT,*) '  location of minimum ZRVS+ZRCS+ZRIS:',MINLOC(ZRVS+ZRCS+ZRIS)
  WRITE(ILUOUT,*) '  value of minimum    ZRVS+ZRCS+ZRIS:',MINVAL(ZRVS+ZRCS+ZRIS)
END IF
!
WHERE ( ZRVS(:,:,:)+ZRCS(:,:,:)+ZRIS(:,:,:) < 0.)
  ZRVS(:,:,:) = -  ZRCS(:,:,:) - ZRIS(:,:,:)
END WHERE
!
!*       2.2    estimate the Exner function at t+1
!
ZEXNS(:,:,:) = ( (2. * PPABST(:,:,:) - PPABSM(:,:,:)) / XP00 ) ** (XRD/XCPD)  
!
!    beginning of the iterative loop
!
DO JITER =1,ITERMAX
!
!*       2.3    compute the intermediate temperature at t+1, T*
!  
  ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
!
!*       2.4    compute the specific heat for moist air (Cph) at t+1
!
  ZCPH(:,:,:) = XCPD + XCPV  *ZDT*   ZRVS(:,:,:)                             &
                     + XCL   *ZDT* ( ZRCS(:,:,:) + ZRRS(:,:,:) )             &
                     + XCI   *ZDT* ( ZRIS(:,:,:) + ZRSS(:,:,:) + ZRGS(:,:,:) )
!
!*       2.5    compute the latent heat of vaporization Lv(T*) at t+1
!               and of sublimation Ls(T*) at t+1
!
  ZLV(:,:,:) = XLVTT + ( XCPV - XCL ) * ( ZT(:,:,:) -XTT )
  ZLS(:,:,:) = XLSTT + ( XCPV - XCI ) * ( ZT(:,:,:) -XTT )
!
!
!-------------------------------------------------------------------------------
!
!*       3.     FIRST ORDER SUBGRID CONDENSATION SCHEME
!               ---------------------------------------
!
  IF ( OSUBG_COND ) THEN
!
! not yet available
!
     STOP
  ELSE
!
!-------------------------------------------------------------------------------
!
!
!*       4.     FULLY EXPLICIT SCHEME FROM TZIVION et al. (1989)
!               -----------------------------------------------
! 
!*              select cases where r_i>0 and r_c=0
! 
    GMICRO(:,:,:) = .FALSE.
    GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                         &
                             (ZRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4) .AND.   &
                              ZCIS_TOT(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(4))    &
                 .AND. .NOT. (ZRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) .AND.   &
                              ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2)      )
    GMICRO_RI(:,:,:) = GMICRO(:,:,:)
    IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!++je++ 17/02/20 test //
!    IF( IMICRO >= 1 ) THEN
    IF (IMICRO >= 0) THEN
!--je--
      ALLOCATE(ZRIT_V(IMICRO))
      ALLOCATE(ZCIT_TOT_V(IMICRO))
      ALLOCATE(ZCIS_TOT_V(IMICRO))
!
      ALLOCATE(ZRVS_V(IMICRO))
      ALLOCATE(ZRIS_V(IMICRO))
      ALLOCATE(ZTHS_V(IMICRO))
!
      ALLOCATE(ZRHODREF(IMICRO))
      ALLOCATE(ZZT(IMICRO))
      ALLOCATE(ZZT0(IMICRO)) ! temperature at t
      ALLOCATE(ZZTM(IMICRO)) ! temperature at t+1
      ALLOCATE(ZPRES(IMICRO))
      ALLOCATE(ZEXNREF(IMICRO))
      ALLOCATE(ZZCPH(IMICRO))
      IF (LCRYSTAL_SHAPE) THEN
        ALLOCATE(ZCSHAPE_T_V(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCSHAPE_S_V(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCIS_SHAPE_V(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCIT_SHAPE_V(IMICRO,NB_CRYSTAL_SHAPE))
      END IF
      DO JL = 1, IMICRO
        ZRIT_V(JL) = ZRIT(I1(JL),I2(JL),I3(JL))
        ZCIT_TOT_V(JL) = ZCIT_TOT(I1(JL),I2(JL),I3(JL))
!
        ZRVS_V(JL) = ZRVS(I1(JL),I2(JL),I3(JL))
        ZRIS_V(JL) = ZRIS(I1(JL),I2(JL),I3(JL))
        ZCIS_TOT_V(JL) = ZCIS_TOT(I1(JL),I2(JL),I3(JL))
        ZTHS_V(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
        IF (LCRYSTAL_SHAPE) THEN
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCSHAPE_T_V(JL,JSH)  = ZCSHAPE_T(I1(JL),I2(JL),I3(JL),JSH)
            ZCSHAPE_S_V(JL,JSH)  = ZCSHAPE_S(I1(JL),I2(JL),I3(JL),JSH)
            ZCIS_SHAPE_V(JL,JSH) = ZCIS_SHAPE(I1(JL),I2(JL),I3(JL),JSH)
            ZCIT_SHAPE_V(JL,JSH) = ZCIT_SHAPE(I1(JL),I2(JL),I3(JL),JSH)
          END DO
        END IF
!
        ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
        ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
        ZPRES(JL) = 2.0*PPABST(I1(JL),I2(JL),I3(JL))-PPABSM(I1(JL),I2(JL),I3(JL))
        ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
        ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
!
        ZZT0(JL) = ZT0(I1(JL),I2(JL),I3(JL))
        ZZTM(JL) = ZTM(I1(JL),I2(JL),I3(JL))
      ENDDO
      ALLOCATE(ZZW(IMICRO))
      ALLOCATE(ZLSFACT(IMICRO))
      ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
      ALLOCATE(ZRVSATI(IMICRO))
      ALLOCATE(ZRVSATI_PRIME(IMICRO))
      ALLOCATE(ZDELTI(IMICRO))
      ALLOCATE(ZAI(IMICRO))
      ALLOCATE(ZCJ(IMICRO))
      ALLOCATE(ZKA(IMICRO))
      ALLOCATE(ZDV(IMICRO))
      ALLOCATE(ZITI(IMICRO))
!
      ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
      ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
      ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
      ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                          * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
!
      ZDELTI(:) = ZRVS_V(:)*ZDT - ZRVSATI(:)
      ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                     + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIT_TOT_V(:),XCTMIN(4))                       &
                              /(MAX(ZRIT_V(:),XRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
        ZITI(:) = ZCIT_TOT_V(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                        / (ZRVSATI(:)*ZAI(:))
      ELSE
        ALLOCATE(ISHAPE_MAX_S(IMICRO))
        ALLOCATE(ISHAPE_MAX_T(IMICRO))
        ALLOCATE(ZZW_2D(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZITI_2D(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZRIT_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZRIS_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZLBDAI_SHAPE_S(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZLBDAI_SHAPE_T(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZSHAPE_NAME_M(IMICRO))
        ALLOCATE(ZSHAPE_NAME_0(IMICRO))
!         
        ZRIT_SHAPE_V(:,:) = 0.0
        ZRIS_SHAPE_V(:,:) = 0.0
        ZZW_2D(:,:) = 0.0
        ZITI_2D(:,:) = 0.0
        ZLBDAI_SHAPE_S(:,:) = 1.E10
        ZLBDAI_SHAPE_T(:,:) = 1.E10
! Calcul de ZLBDAI pour retrouver le rapport de melange par forme
! Hyp : on considere la forme preponderante par point de grille pour le calcul de lambda
        ISHAPE_MAX_S(:) = MAXLOC(ZCSHAPE_S_V,DIM=2)
        ISHAPE_MAX_T(:) = MAXLOC(ZCSHAPE_T_V,DIM=2)
! calcul de ZLBDAI_S et T pour prendre en compte les valeurs differentes de ZCIS et ZCIT
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX_S, IMICRO, PTSTEP,       &
                                ZRIT_V, ZCIT_TOT_V, ZLBDAI_SHAPE_S, &
                                ZRIS_V, ZCIS_SHAPE_V, ZRIS_SHAPE_V)
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX_T, IMICRO, PTSTEP,       &
                                ZRIT_V, ZCIT_TOT_V, ZLBDAI_SHAPE_T, &
                                ZRIT_V, ZCIT_SHAPE_V, ZRIT_SHAPE_V)
! par forme : calcul de ZRIT_V en fonction des shapes (p110 sci doc) : ZRIT_SHAPE_V
! Rajout d'une condition pour le calcul de ZRIS_SHAPE_V
        DO JSH = 1, NB_CRYSTAL_SHAPE
          WHERE (ZRIS_SHAPE_V(:,JSH) > ZRTMIN(4) .AND. ZCIS_SHAPE_V(:,JSH) > ZCTMIN(4)) !++cb-- 2 mai 2018: utile ?
            ZZW_2D(:,JSH) = MIN(1.E8, &
                               (XLBI_SHAPE(JSH)* MAX(ZCIT_SHAPE_V(:,JSH),XCTMIN(4)) &
                             /(MAX(ZRIT_SHAPE_V(:,JSH),XRTMIN(4))))**XLBEXI_SHAPE(JSH))
            ZITI_2D(:,JSH) = ZCIT_SHAPE_V(:,JSH) * (X0DEPI_SHAPE(JSH)/ZZW_2D(:,JSH) + &
                             X2DEPI_SHAPE(JSH) * ZCJ(:) * ZCJ(:) / &
                             ZZW_2D(:,JSH)**(XDI_SHAPE(JSH)+2.0)) / (ZRVSATI(:) * ZAI(:))
          END WHERE
        END DO
      END IF  ! LCRYSTAL_SHAPE 
!
      ALLOCATE(ZAII(IMICRO))
      ALLOCATE(ZDEP(IMICRO))
      IF (LCRYSTAL_SHAPE) ALLOCATE(ZDEP_2D(IMICRO,NB_CRYSTAL_SHAPE)) 
      ZDEP(:) = 0.
!
      ZAII(:) = 1.0 + ZRVSATI_PRIME(:)*ZLSFACT(:)
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZZW(:)  = ZAII(:)*ZITI(:)*ZDT ! R*delta_T
        WHERE( ZZW(:)<1.0E-2 )
          ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - (ZZW(:)/2.0)*(1.0-ZZW(:)/3.0))
        ELSEWHERE          
          ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - EXP(-ZZW(:)))/ZZW(:)
        END WHERE
      ELSE 
        ZDEP_2D(:,:) = 0.0
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZZW_2D(:,JSH) = ZAII(:) * ZITI_2D(:,JSH) * ZDT ! R*delta_T
          WHERE( ZZW_2D(:,JSH) < 1.0E-2 )
            ZDEP_2D(:,JSH) = ZITI_2D(:,JSH) * ZDELTI(:) * &
                            (1.0 - (ZZW_2D(:,JSH) / 2.0) * &
                            (1.0 - ZZW_2D(:,JSH) / 3.0))
          ELSEWHERE         
            ZDEP_2D(:,JSH) = ZITI_2D(:,JSH) * ZDELTI(:) * &
                            (1.0 - EXP(-ZZW_2D(:,JSH))) / ZZW_2D(:,JSH)
          END WHERE
        END DO
        ZDEP(:) = SUM(ZDEP_2D, DIM=2)
      END IF ! LCRYSTAL_SHAPE
!
! Integration
!
      WHERE( ZDEP(:) < 0.0 )
         ZDEP(:) = MAX ( ZDEP(:), -ZRIS_V(:) )
      ELSEWHERE
         ZDEP(:) = MIN ( ZDEP(:),  ZRVS_V(:) )
!         ZDEP(:) = MIN ( ZDEP(:),  ZCIS(:)*5.E-10 ) !!!BVIE!!!
      END WHERE
      WHERE( ZRIS_V(:) < ZRTMIN(4) )
         ZDEP(:) = 0.0
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZDEP(:)
      ZRIS_V(:) = ZRIS_V(:) + ZDEP(:)
      ZTHS_V(:) = ZTHS_V(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
! Ice crystal can change form
!
      IF (LHABIT_CHANGE .AND. NB_CRYSTAL_SHAPE .GE. 3) THEN
        ZSHAPE_NAME_M(:) = FIND_SHAPE(ZZTM(:))
        ZSHAPE_NAME_0(:) = FIND_SHAPE(ZZT0(:))
        ALLOCATE(GHABIT(IMICRO))
        GHABIT(:) = .FALSE.
        GHABIT(:) = (((ZSHAPE_NAME_0(:) .EQ. "PLA"  .AND. ZCIS_SHAPE_V(:,2) .GT. ZCTMIN(4)          .AND. &
                                                          ZCIS_SHAPE_V(:,2) .LT. ZCIS_SHAPE_V(:,1)) .OR. &
                      (ZSHAPE_NAME_0(:) .EQ. "COL"  .AND. ZCIS_SHAPE_V(:,1) .GT. ZCTMIN(4)          .AND. &
                                                          ZCIS_SHAPE_V(:,1) .LT. ZCIS_SHAPE_V(:,2))) .AND.  &
                       ZDEP(1:IMICRO) .GT. 0. )
        IHABIT = COUNT(GHABIT)
!
        IF (IHABIT >= 1) THEN
          ALLOCATE(ZZTM_V(IHABIT))
          ALLOCATE(ZZT0_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_M_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_0_V(IHABIT))
          ALLOCATE(ZCIS_SHAPE_VV(IHABIT,NB_CRYSTAL_SHAPE))
          ALLOCATE(ZDEP_V(IHABIT))
! 
          ZZTM_V(:) = PACK(ZZTM(:),GHABIT(:))
          ZZT0_V(:) = PACK(ZZT0(:),GHABIT(:))
          ZSHAPE_NAME_M_V(:) = PACK(ZSHAPE_NAME_M(:),GHABIT(:))
          ZSHAPE_NAME_0_V(:) = PACK(ZSHAPE_NAME_0(:),GHABIT(:))
          ZDEP_V(:) = PACK(ZDEP(:),GHABIT(:))
!
!          ALLOCATE(ZDUM(IMICRO,NB_CRYSTAL_SHAPE))
!          ZDUM(:,:) = ZCIS_SHAPE_V(:,:)
!
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCIS_SHAPE_VV(:,JSH) = PACK(ZCIS_SHAPE_V(:,JSH),GHABIT(:))
          END DO
!
          CALL LIMA_CHANGE_SHAPE(ZZTM_V, ZZT0_V, ZSHAPE_NAME_M_V, &
                                 ZSHAPE_NAME_0_V, ZCIS_SHAPE_VV,ZCTMIN,ZDEP_V)
!
!   Unpack variables
          ALLOCATE(ZW_V(IMICRO, NB_CRYSTAL_SHAPE))
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZW_V(:,JSH) = ZCIS_SHAPE_V(:,JSH)
            ZCIS_SHAPE_V(:,JSH) = UNPACK(ZCIS_SHAPE_VV(:,JSH),MASK=GHABIT(:),FIELD=ZW_V(:,JSH))
          END DO
!
          DEALLOCATE(ZCIS_SHAPE_VV)
          DEALLOCATE(ZZTM_V)
          DEALLOCATE(ZZT0_V)
          DEALLOCATE(ZSHAPE_NAME_M_V)
          DEALLOCATE(ZSHAPE_NAME_0_V)
          DEALLOCATE(ZW_V)
          DEALLOCATE(ZDEP_V)
        END IF
!
        DEALLOCATE(GHABIT)
!
      END IF
!
!  Implicit ice crystal sublimation if ice saturated conditions are not met
!
      ZZT(:) = ( ZTHS_V(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
      WHERE( ZRVS_V(:)*ZDT<ZRVSATI(:) )
        ZZW(:)  = ZRVS_V(:) + ZRIS_V(:)
        ZRVS_V(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
        ZTHS_V(:) = ZTHS_V(:) + ( MAX( 0.0,ZZW(:)-ZRVS_V(:) )-ZRIS_V(:) ) &
                            * ZLSFACT(:) / ZEXNREF(:)
        ZRIS_V(:) = MAX( 0.0,ZZW(:)-ZRVS_V(:) )
      END WHERE
!
!
      ZW(:,:,:) = ZRVS(:,:,:)
      ZRVS(:,:,:) = UNPACK( ZRVS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = ZRIS(:,:,:)
      ZRIS(:,:,:) = UNPACK( ZRIS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = PTHS(:,:,:)
      PTHS(:,:,:) = UNPACK( ZTHS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
      IF (LCRYSTAL_SHAPE) THEN
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZW(:,:,:) = ZCIS_SHAPE(:,:,:,JSH)
          ZCIS_SHAPE(:,:,:,JSH) = UNPACK(ZCIS_SHAPE_V(:,JSH),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:))
        END DO
      END IF
!
      DEALLOCATE(ZRIT_V)
      DEALLOCATE(ZCIT_TOT_V)
      DEALLOCATE(ZRVS_V)
      DEALLOCATE(ZRIS_V)
      DEALLOCATE(ZCIS_TOT_V)
      DEALLOCATE(ZTHS_V)
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZZT)
      DEALLOCATE(ZZT0)
      DEALLOCATE(ZZTM)
      DEALLOCATE(ZPRES)
      DEALLOCATE(ZEXNREF)
      DEALLOCATE(ZZCPH)
      DEALLOCATE(ZZW)
      DEALLOCATE(ZLSFACT)
      DEALLOCATE(ZRVSATI)
      DEALLOCATE(ZRVSATI_PRIME)
      DEALLOCATE(ZDELTI)
      DEALLOCATE(ZAI)
      DEALLOCATE(ZCJ)
      DEALLOCATE(ZKA)
      DEALLOCATE(ZDV)
      DEALLOCATE(ZITI)
      DEALLOCATE(ZAII)
!
      IF (ALLOCATED(ZDEP))           DEALLOCATE(ZDEP)
      IF (ALLOCATED(ZCIT_SHAPE_V))   DEALLOCATE(ZCIT_SHAPE_V)
      IF (ALLOCATED(ZCIS_SHAPE_V))   DEALLOCATE(ZCIS_SHAPE_V)
      IF (ALLOCATED(ZCSHAPE_T_V))    DEALLOCATE(ZCSHAPE_T_V)
      IF (ALLOCATED(ZCSHAPE_S_V))    DEALLOCATE(ZCSHAPE_S_V)
      IF (ALLOCATED(ISHAPE_MAX_S))   DEALLOCATE(ISHAPE_MAX_S)
      IF (ALLOCATED(ISHAPE_MAX_T))   DEALLOCATE(ISHAPE_MAX_T)
      IF (ALLOCATED(ZLBDAI_SHAPE_S)) DEALLOCATE(ZLBDAI_SHAPE_S)
      IF (ALLOCATED(ZLBDAI_SHAPE_T)) DEALLOCATE(ZLBDAI_SHAPE_T)
      IF (ALLOCATED(ZSHAPE_NAME_M))  DEALLOCATE(ZSHAPE_NAME_M)
      IF (ALLOCATED(ZSHAPE_NAME_0))  DEALLOCATE(ZSHAPE_NAME_0)
      IF (ALLOCATED(ZZW_2D))         DEALLOCATE(ZZW_2D)
      IF (ALLOCATED(ZDEP_2D))        DEALLOCATE(ZDEP_2D)
      IF (ALLOCATED(ZITI_2D))        DEALLOCATE(ZITI_2D)
      IF (ALLOCATED(ZRIT_SHAPE_V))   DEALLOCATE(ZRIT_SHAPE_V)
      IF (ALLOCATED(ZRIS_SHAPE_V))   DEALLOCATE(ZRIS_SHAPE_V)
    END IF ! IMICRO
!
!
!-------------------------------------------------------------------------------
!
!
!*       5.     FULLY IMPLICIT CONDENSATION SCHEME
!               ---------------------------------
! 
!*              select cases where r_c>0 and r_i=0
! 
!
    GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                      & 
              .NOT. GMICRO_RI(IIB:IIE,IJB:IJE,IKB:IKE)                     &
              .AND. ( ZRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) .AND.        &
                      ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2)      )       &
        .AND. .NOT. ( ZRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4) .AND.        &
                      ZCIS_TOT(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(4)      )
    GMICRO_RC(:,:,:) = GMICRO(:,:,:)
    IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!++je++ 19/02/20 test //
!    IF( IMICRO >= 1 ) THEN
    IF ( IMICRO >= 0) THEN
!--je--
      ALLOCATE(ZRCT_V(IMICRO))
!
      ALLOCATE(ZRVS_V(IMICRO))
      ALLOCATE(ZRCS_V(IMICRO))
      ALLOCATE(ZTHS_V(IMICRO))
!
      ALLOCATE(ZRHODREF(IMICRO))
      ALLOCATE(ZZT(IMICRO))
      ALLOCATE(ZPRES(IMICRO))
      ALLOCATE(ZEXNREF(IMICRO))
      ALLOCATE(ZZCPH(IMICRO))
      DO JL = 1, IMICRO
        ZRCT_V(JL) = ZRCT(I1(JL),I2(JL),I3(JL))
!
        ZRVS_V(JL) = ZRVS(I1(JL),I2(JL),I3(JL))
        ZRCS_V(JL) = ZRCS(I1(JL),I2(JL),I3(JL))
        ZTHS_V(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
        ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
        ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
        ZPRES(JL) = 2.0*PPABST(I1(JL),I2(JL),I3(JL))-PPABSM(I1(JL),I2(JL),I3(JL))
        ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
        ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
      ENDDO
      ALLOCATE(ZZW(IMICRO))
      ALLOCATE(ZLVFACT(IMICRO))
      ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
      ALLOCATE(ZRVSATW(IMICRO))
      ALLOCATE(ZRVSATW_PRIME(IMICRO))
!
      ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
      ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw
      ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                         * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
      ALLOCATE(ZAWW(IMICRO))
      ALLOCATE(ZDELT1(IMICRO))
      ALLOCATE(ZDELT2(IMICRO))
      ALLOCATE(ZCND(IMICRO))
!
      ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZLVFACT(:)
      ZDELT2(:) = (ZRVSATW_PRIME(:)*ZLVFACT(:)/ZAWW(:)) *                     &
                  ( ((-2.*XBETAW+XGAMW*ZZT(:))/(XBETAW-XGAMW*ZZT(:))          &
                  + (XBETAW/ZZT(:)-XGAMW)*(1.0+2.0*ZRVSATW(:)/ZEPS))/ZZT(:) )
      ZDELT1(:) = (ZLVFACT(:)/ZAWW(:)) * ( ZRVSATW(:) - ZRVS_V(:)*ZDT )
      ZCND(:) = - ZDELT1(:)*( 1.0 + 0.5*ZDELT1(:)*ZDELT2(:) ) / (ZLVFACT(:)*ZDT)
!
! Integration
!
      WHERE( ZCND(:) < 0.0 )
        ZCND(:) = MAX ( ZCND(:), -ZRCS_V(:) )
      ELSEWHERE
        ZCND(:) = MIN ( ZCND(:),  ZRVS_V(:) )
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZCND(:)
      ZRCS_V(:) = ZRCS_V(:) + ZCND(:)
      ZTHS_V(:) = ZTHS_V(:) + ZCND(:) * ZLVFACT(:) / ZEXNREF(:)
!
      ZW(:,:,:) = ZRVS(:,:,:)
      ZRVS(:,:,:) = UNPACK( ZRVS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = ZRCS(:,:,:)
      ZRCS(:,:,:) = UNPACK( ZRCS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = PTHS(:,:,:)
      PTHS(:,:,:) = UNPACK( ZTHS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
      DEALLOCATE(ZRCT_V)
      DEALLOCATE(ZRVS_V)
      DEALLOCATE(ZRCS_V)
      DEALLOCATE(ZTHS_V)
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZZT)
      DEALLOCATE(ZPRES)
      DEALLOCATE(ZEXNREF)
      DEALLOCATE(ZZCPH)
      DEALLOCATE(ZZW)
      DEALLOCATE(ZLVFACT)
      DEALLOCATE(ZRVSATW)
      DEALLOCATE(ZRVSATW_PRIME)
      DEALLOCATE(ZAWW)
      DEALLOCATE(ZDELT1)
      DEALLOCATE(ZDELT2)
      DEALLOCATE(ZCND)
    END IF ! IMICRO
!
!
!-------------------------------------------------------------------------------
!
!
!*       6.     IMPLICIT-EXPLICIT SCHEME USING REISIN et al. (1996)
!               ---------------------------------------------------
! 
!*              select cases where r_i>0 and r_c>0 (supercooled water)
! 
!
    GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                &
               .NOT. GMICRO_RI(IIB:IIE,IJB:IJE,IKB:IKE)              &
         .AND. .NOT. GMICRO_RC(IIB:IIE,IJB:IJE,IKB:IKE)              &
               .AND. ( ZRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4) .AND. &
                       ZCIS_TOT(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(4) ) &
               .AND. ( ZRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) .AND. &
                       ZCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2) )
    IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!++je++ 19/02/20 test //
!    IF( IMICRO >= 1 ) THEN
    IF ( IMICRO >= 0) THEN
!--je--
      ALLOCATE(ZRCT_V(IMICRO))
      ALLOCATE(ZRIT_V(IMICRO))
      ALLOCATE(ZCCT_V(IMICRO))
      ALLOCATE(ZCIT_TOT_V(IMICRO))
!
      ALLOCATE(ZRVS_V(IMICRO))
      ALLOCATE(ZRCS_V(IMICRO))
      ALLOCATE(ZRIS_V(IMICRO))
      ALLOCATE(ZCCS_V(IMICRO))
      ALLOCATE(ZCIS_TOT_V(IMICRO))
      ALLOCATE(ZTHS_V(IMICRO))
!
      IF (LCRYSTAL_SHAPE) THEN
        ALLOCATE(ZCSHAPE_T_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCSHAPE_S_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCIS_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCIT_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
      END IF
!
      ALLOCATE(ZRHODREF(IMICRO))
      ALLOCATE(ZZT(IMICRO))
      ALLOCATE(ZZTM(IMICRO))
      ALLOCATE(ZZT0(IMICRO))
      ALLOCATE(ZPRES(IMICRO))
      ALLOCATE(ZEXNREF(IMICRO))
      ALLOCATE(ZZCPH(IMICRO))
      DO JL=1,IMICRO
        ZRCT_V(JL) = ZRCT(I1(JL),I2(JL),I3(JL))
        ZRIT_V(JL) = ZRIT(I1(JL),I2(JL),I3(JL))
        ZCCT_V(JL) = ZCCT(I1(JL),I2(JL),I3(JL))
        ZCIT_TOT_V(JL) = ZCIT_TOT(I1(JL),I2(JL),I3(JL))
!
        ZRVS_V(JL) = ZRVS(I1(JL),I2(JL),I3(JL))
        ZRCS_V(JL) = ZRCS(I1(JL),I2(JL),I3(JL))
        ZRIS_V(JL) = ZRIS(I1(JL),I2(JL),I3(JL))
        ZCCS_V(JL) = ZCCS(I1(JL),I2(JL),I3(JL))
        ZCIS_TOT_V(JL) = ZCIS_TOT(I1(JL),I2(JL),I3(JL))
        ZTHS_V(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
        IF (LCRYSTAL_SHAPE) THEN
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCSHAPE_T_V(JL,JSH) = ZCSHAPE_T(I1(JL),I2(JL),I3(JL),JSH)   
            ZCSHAPE_S_V(JL,JSH) = ZCSHAPE_S(I1(JL),I2(JL),I3(JL),JSH)
            ZCIS_SHAPE_V(JL,JSH) = ZCIS_SHAPE(I1(JL),I2(JL),I3(JL),JSH)   
            ZCIT_SHAPE_V(JL,JSH) = ZCIT_SHAPE(I1(JL),I2(JL),I3(JL),JSH)   
          END DO
        END IF
!
        ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
        ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
        ZZT0(JL) = ZT0(I1(JL),I2(JL),I3(JL))
        ZZTM(JL) = ZTM(I1(JL),I2(JL),I3(JL))
        ZPRES(JL) = 2.0*PPABST(I1(JL),I2(JL),I3(JL))-PPABSM(I1(JL),I2(JL),I3(JL))
        ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
        ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
      ENDDO
      ALLOCATE(ZZW(IMICRO))
      ALLOCATE(ZLVFACT(IMICRO))
      ALLOCATE(ZLSFACT(IMICRO))
      ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
      ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
      ALLOCATE(ZRVSATW(IMICRO))
      ALLOCATE(ZRVSATI(IMICRO))
      ALLOCATE(ZRVSATW_PRIME(IMICRO))
      ALLOCATE(ZRVSATI_PRIME(IMICRO))
      ALLOCATE(ZDELTW(IMICRO))
      ALLOCATE(ZDELTI(IMICRO))
      ALLOCATE(ZAW(IMICRO))
      ALLOCATE(ZAI(IMICRO))
      ALLOCATE(ZCJ(IMICRO))
      ALLOCATE(ZKA(IMICRO))
      ALLOCATE(ZDV(IMICRO))
      ALLOCATE(ZITW(IMICRO))
      ALLOCATE(ZITI(IMICRO))
!
      ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
      ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
      ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!*       6.2    implicit adjustment at water saturation
!
      ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
      ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw
      ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                         * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
      ZDELTW(:) = ABS( ZRVS_V(:)*ZDT - ZRVSATW(:) )
      ZAW(:) = ( XLSTT + (XCPV-XCL)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                     + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
      ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                         * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
      ZDELTI(:) = ABS( ZRVS_V(:)*ZDT - ZRVSATI(:) )
      ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                     + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!
      ZZW(:) = MIN(1.E8,( XLBC* MAX(ZCCT_V(:),XCTMIN(2))                       &
                              /(MAX(ZRCT_V(:),XRTMIN(2))) )**XLBEXC)
                                                                  ! Lbda_c
!
      ZITW(:) = ZCCT_V(:) * (X0CNDC/ZZW(:) + X2CNDC*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDC+2.0)) &
                        / (ZRVSATW(:)*ZAW(:))
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIT_TOT_V(:),XCTMIN(4))                       &
                          /(MAX(ZRIT_V(:),XRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
        ZITI(:) = ZCIT_TOT_V(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                          / (ZRVSATI(:)*ZAI(:))
!
      ELSE
!
        ALLOCATE (ZZW_2D(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE (ZITI_2D(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE (ZCND_2D(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE (ZDEP_2D(IMICRO,NB_CRYSTAL_SHAPE))
        ALLOCATE (ISHAPE_MAX_S(IMICRO))
        ALLOCATE (ISHAPE_MAX_T(IMICRO))
        ALLOCATE (ZRIT_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE (ZRIS_SHAPE_V(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE (ZLBDAI_SHAPE_T(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE (ZLBDAI_SHAPE_S(IMICRO, NB_CRYSTAL_SHAPE))
        ALLOCATE (ZSHAPE_NAME_M(IMICRO))
        ALLOCATE (ZSHAPE_NAME_0(IMICRO))
!
        ZRIT_SHAPE_V(:,:) = 0.0
        ZRIS_SHAPE_V(:,:) = 0.0
        ZZW_2D(:,:) = 0.0
        ZITI_2D(:,:) = 0.0
        ZLBDAI_SHAPE_S(:,:) = 1.E10
        ZLBDAI_SHAPE_T(:,:) = 1.E10
        ZITI(:) = 0.
!
! prevailing habit in each mesh
        ISHAPE_MAX_T(:) = MAXLOC(ZCSHAPE_T_V,DIM=2)
        ISHAPE_MAX_S(:) = MAXLOC(ZCSHAPE_S_V,DIM=2)
! compute lambda and the mixing ratio per habit
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX_S, IMICRO, PTSTEP,       &
                                ZRIT_V, ZCIT_TOT_V, ZLBDAI_SHAPE_S, &
                                ZRIS_V, ZCIS_SHAPE_V, ZRIS_SHAPE_V)
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX_T, IMICRO, PTSTEP,       &
                                ZRIT_V, ZCIS_TOT_V, ZLBDAI_SHAPE_T, &
                                ZRIT_V, ZCIT_SHAPE_V, ZRIT_SHAPE_V)
!
! on reprend les conditions de gmicro et on les adapte pour etre sur qu'on ne traite 
! qu'une fois chaque maille si une seule espece presente par maille
        DO JSH = 1, NB_CRYSTAL_SHAPE
          WHERE (ZRIS_SHAPE_V(:,JSH) > ZRTMIN(4) .AND. ZCIS_SHAPE_V(:,JSH) > ZCTMIN(4))
            ZZW_2D(:,JSH) = MIN(1.E8,( XLBI_SHAPE(JSH)* MAX(ZCIT_SHAPE_V(:,JSH),XCTMIN(4)) &
                          /(MAX(ZRIT_SHAPE_V(:,JSH),XRTMIN(4))))**XLBEXI_SHAPE(JSH))
            ZITI_2D(:,JSH) = ZCIT_SHAPE_V(:,JSH) * (X0DEPI_SHAPE(JSH) / ZZW_2D(:,JSH) + &
                             X2DEPI_SHAPE(JSH)  *ZCJ(:) * ZCJ(:) / &
                             ZZW_2D(:,JSH)**(XDI_SHAPE(JSH)+2.0)) / (ZRVSATI(:) * ZAI(:))
          END WHERE
        END DO
!
! to avoid double contribution from ice crystal, ziti is computed here
        ZITI(:) = SUM(ZITI_2D, DIM=2)
      END IF  ! lcrystal_shape
!
      ALLOCATE(ZAWW(IMICRO))
      ALLOCATE(ZAIW(IMICRO))
      ALLOCATE(ZAWI(IMICRO))
      ALLOCATE(ZAII(IMICRO))
!
      ALLOCATE(ZFACT(IMICRO))
      ALLOCATE(ZDELT1(IMICRO))
      ALLOCATE(ZDELT2(IMICRO))
!
      ZAII(:)  = ZITI(:)*ZDELTI(:)
      WHERE( ZAII(:)<1.0E-15 )
        ZFACT(:) = ZLVFACT(:)
      ELSEWHERE          
        ZFACT(:) = (ZLVFACT(:)*ZITW(:)*ZDELTW(:)+ZLSFACT(:)*ZITI(:)*ZDELTI(:)) &
                          / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
      END WHERE
      ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZFACT(:)
!
      ZDELT2(:) = (ZRVSATW_PRIME(:)*ZFACT(:)/ZAWW(:)) *                       &
                  ( ((-2.*XBETAW+XGAMW*ZZT(:))/(XBETAW-XGAMW*ZZT(:))          &
                    + (XBETAW/ZZT(:)-XGAMW)*(1.0+2.0*ZRVSATW(:)/ZEPS))/ZZT(:) )
      ZDELT1(:) = (ZFACT(:)/ZAWW(:)) * ( ZRVSATW(:) - ZRVS_V(:)*ZDT )
!
      ALLOCATE(ZCND(IMICRO))
      ALLOCATE(ZDEP(IMICRO))
      ZCND(:) = 0.0
      ZDEP(:) = 0.0
!
      ZZW(:) =  - ZDELT1(:)*( 1.0 + 0.5*ZDELT1(:)*ZDELT2(:) ) / (ZFACT(:)*ZDT) 
      WHERE( ZAII(:)<1.0E-15 )
        ZCND(:) = ZZW(:)
        ZDEP(:) = 0.0
      ELSEWHERE          
        ZCND(:) = ZZW(:)*ZITW(:)*ZDELTW(:) / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
        ZDEP(:) = ZZW(:)*ZITI(:)*ZDELTI(:) / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
      END WHERE
!
! Integration
!
      WHERE( ZCND(:) < 0.0 )
        ZCND(:) = MAX ( ZCND(:), -ZRCS_V(:) )
      ELSEWHERE 
        ZCND(:) = MIN ( ZCND(:),  ZRVS_V(:) )
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZCND(:)
      ZRCS_V(:) = ZRCS_V(:) + ZCND(:)
      ZTHS_V(:) = ZTHS_V(:) + ZCND(:) * ZLVFACT(:) / ZEXNREF(:)
!
      WHERE( ZDEP(:) < 0.0 )
        ZDEP(:) = MAX ( ZDEP(:), -ZRIS_V(:) )
      ELSEWHERE
        ZDEP(:) = MIN ( ZDEP(:),  ZRVS_V(:) )
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZDEP(:)
      ZRIS_V(:) = ZRIS_V(:) + ZDEP(:)
      ZTHS_V(:) = ZTHS_V(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
! Ice crystal can change form
!
      IF (LHABIT_CHANGE .AND. NB_CRYSTAL_SHAPE .GE. 3) THEN
        ZSHAPE_NAME_M(:) = FIND_SHAPE(ZZTM(:))
        ZSHAPE_NAME_0(:) = FIND_SHAPE(ZZT0(:))
!
        ALLOCATE(GHABIT(IMICRO))
        GHABIT(:) = .FALSE.
        GHABIT(:) = (((ZSHAPE_NAME_0(:) .EQ. "PLA" .AND. ZCIS_SHAPE_V(:,2) .GT. ZCTMIN(4) .AND. &
                                                         ZCIS_SHAPE_V(:,2) .LT. ZCIS_SHAPE_V(:,1)) .OR.  &
                      (ZSHAPE_NAME_0(:) .EQ. "COL" .AND. ZCIS_SHAPE_V(:,1) .GT. ZCTMIN(4) .AND. &
                                                         ZCIS_SHAPE_V(:,1) .LT. ZCIS_SHAPE_V(:,2))) .AND. &
                       ZDEP(:) .GT. 0.0) 
!
        IHABIT = COUNT(GHABIT)
!
        IF (IHABIT >= 1) THEN
          ALLOCATE(ZZTM_V(IHABIT))
          ALLOCATE(ZZT0_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_M_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_0_V(IHABIT))
          ALLOCATE(ZCIS_SHAPE_VV(IHABIT,NB_CRYSTAL_SHAPE))
          ALLOCATE(ZDEP_V(IHABIT))
!
!          ALLOCATE(ZDUM(IMICRO,NB_CRYSTAL_SHAPE))
!          ZDUM(:,:) = ZCIS_SHAPE_V(:,:)
!
          ZZTM_V(:) = PACK(ZZTM(:),GHABIT(:))
          ZZT0_V(:) = PACK(ZZT0(:),GHABIT(:))
          ZSHAPE_NAME_M_V(:) = PACK(ZSHAPE_NAME_M(:),GHABIT(:))
          ZSHAPE_NAME_0_V(:) = PACK(ZSHAPE_NAME_0(:),GHABIT(:))
          ZDEP_V(:) = PACK(ZDEP(:),GHABIT(:))
!
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCIS_SHAPE_VV(:,JSH) = PACK(ZCIS_SHAPE_V(:,JSH),GHABIT(:))
          END DO
!
          CALL LIMA_CHANGE_SHAPE(ZZTM_V, ZZT0_V, ZSHAPE_NAME_M_V, &
                                 ZSHAPE_NAME_0_V, ZCIS_SHAPE_VV,ZCTMIN,ZDEP_V)
!
!   Unpack variables
          ALLOCATE(ZW_V(IMICRO, NB_CRYSTAL_SHAPE))
          ZW_V(:,:) = 0.
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZW_V(:,JSH) = ZCIS_SHAPE_V(:,JSH)
            ZCIS_SHAPE_V(:,JSH) = UNPACK(ZCIS_SHAPE_VV(:,JSH),MASK=GHABIT(:),FIELD=ZW_V(:,JSH))
          END DO
!
          DEALLOCATE(ZCIS_SHAPE_VV)
          DEALLOCATE(ZZTM_V)
          DEALLOCATE(ZZT0_V)
          DEALLOCATE(ZSHAPE_NAME_M_V)
          DEALLOCATE(ZSHAPE_NAME_0_V)
          DEALLOCATE(ZW_V)
          DEALLOCATE(ZDEP_V)
        END IF ! ihabit
        DEALLOCATE(GHABIT)
      END IF ! lhabit
!
!
!*       6.3    explicit integration of the final eva/dep rates
!
      ZZT(:) = ( ZTHS_V(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
!
!  If Si < 0, implicit adjustment to Si=0 using ice only
!
      WHERE( ZRVS_V(:)*ZDT<ZRVSATI(:) )
         ZZW(:)  = ZRVS_V(:) + ZRIS_V(:)
         ZRVS_V(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
         ZTHS_V(:) = ZTHS_V(:) + ( MAX( 0.0,ZZW(:)-ZRVS_V(:) )-ZRIS_V(:) ) &
                             * ZLSFACT(:) / ZEXNREF(:)
         ZRIS_V(:) = MAX( 0.0,ZZW(:)-ZRVS_V(:) )
      END WHERE
!
!  Following the previous adjustment, the real procedure begins
!
      ZZT(:) = ( ZTHS_V(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
!
      ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
      ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
!
      ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
      ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
      ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
      ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
      ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw
      ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                         * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
      ZDELTW(:) = ZRVS_V(:)*ZDT - ZRVSATW(:)
      ZAW(:) = ( XLSTT + (XCPV-XCL)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                     + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
      ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                         * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
      ZDELTI(:) = ZRVS_V(:)*ZDT - ZRVSATI(:)
      ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                    + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!                    
      ZZW(:) = MIN(1.E8,( XLBC* MAX(ZCCS_V(:),ZCTMIN(2)) &
                              /(MAX(ZRCS_V(:),ZRTMIN(2))) )**XLBEXC)
                                                                  ! Lbda_c
      ZITW(:) = ZCCT_V(:) * (X0CNDC/ZZW(:) + X2CNDC*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDC+2.0)) &
                        / (ZRVSATW(:)*ZAW(:))
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIS_TOT_V(:),ZCTMIN(4)) &
                                /(MAX(ZRIS_V(:),ZRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
        ZITI(:) = ZCIT_TOT_V(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                       / (ZRVSATI(:)*ZAI(:))
      ELSE
        ZITI_2D(:,:) = 0.
        ZZW_2D(:,:)  = 0.
!
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX_S, IMICRO, PTSTEP,       &
                                ZRIT_V, ZCIT_TOT_V, ZLBDAI_SHAPE_S, &
                                ZRIS_V, ZCIS_SHAPE_V, ZRIS_SHAPE_V)
!
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZZW_2D(:,JSH)  = MIN(1.E8,( XLBI_SHAPE(JSH) * MAX(ZCIS_SHAPE_V(:,JSH),ZCTMIN(4)) / &
                               (MAX(ZRIS_SHAPE_V(:,JSH),ZRTMIN(4))) )**XLBEXI_SHAPE(JSH))
                                                                 ! Lbda_I
          ZITI_2D(:,JSH) = ZCIT_SHAPE_V(:,JSH) * (X0DEPI_SHAPE(JSH) / ZZW_2D(:,JSH) + &
                           X2DEPI_SHAPE(JSH) * ZCJ(:) * ZCJ(:) / &
                           ZZW_2D(:,JSH)**(XDI_SHAPE(JSH)+2.0)) / (ZRVSATI(:)*ZAI(:))
        END DO
        ZITI(:) = SUM(ZITI_2D, DIM=2)
      END IF ! lcrystal_shape
!
      ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZLVFACT(:)
      ZAIW(:) = 1.0 + ZRVSATI_PRIME(:)*ZLVFACT(:)
      ZAWI(:) = 1.0 + ZRVSATW_PRIME(:)*ZLSFACT(:)
      ZAII(:) = 1.0 + ZRVSATI_PRIME(:)*ZLSFACT(:)
!
      ZCND(:) = 0.0      
      ZDEP(:) = 0.0
      ZZW(:) = ZAWW(:)*ZITW(:) + ZAII(:)*ZITI(:) ! R
      WHERE( ZZW(:)<1.0E-2 )
        ZFACT(:) = ZDT*(0.5 - (ZZW(:)*ZDT)/6.0)
      ELSEWHERE          
        ZFACT(:) = (1.0/ZZW(:))*(1.0-(1.0-EXP(-ZZW(:)*ZDT))/(ZZW(:)*ZDT))
      END WHERE
      ZCND(:) = ZITW(:)*(ZDELTW(:)-( ZAWW(:)*ZITW(:)*ZDELTW(:)           &
                                   + ZAWI(:)*ZITI(:)*ZDELTI(:) )*ZFACT(:))
      ZDEP(:) = ZITI(:)*(ZDELTI(:)-( ZAIW(:)*ZITW(:)*ZDELTW(:)           &
                                   + ZAII(:)*ZITI(:)*ZDELTI(:) )*ZFACT(:))
!                    
! Integration        
!           
      WHERE( ZCND(:) < 0.0 )
        ZCND(:) = MAX ( ZCND(:), -ZRCS_V(:) )
      ELSEWHERE          
        ZCND(:) = MIN ( ZCND(:),  ZRVS_V(:) )
      END WHERE
      WHERE( ZRCS_V(:) < ZRTMIN(2) )
        ZCND(:) = 0.0
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZCND(:)
      ZRCS_V(:) = ZRCS_V(:) + ZCND(:)
      ZTHS_V(:) = ZTHS_V(:) + ZCND(:) * ZLVFACT(:) / ZEXNREF(:)
!
      WHERE( ZDEP(:) < 0.0 )
        ZDEP(:) = MAX ( ZDEP(:), -ZRIS_V(:) )
      ELSEWHERE          
        ZDEP(:) = MIN ( ZDEP(:),  ZRVS_V(:) )
      END WHERE
      WHERE( ZRIS_V(:) < ZRTMIN(4) )
        ZDEP(:) = 0.0
      END WHERE
      ZRVS_V(:) = ZRVS_V(:) - ZDEP(:)
      ZRIS_V(:) = ZRIS_V(:) + ZDEP(:)
      ZTHS_V(:) = ZTHS_V(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
!
! Ice crystals can change shape
!
      IF (LHABIT_CHANGE .AND. NB_CRYSTAL_SHAPE .GE. 3) THEN
        ALLOCATE(GHABIT(IMICRO))
        GHABIT(:) = .FALSE.
        GHABIT(1:IMICRO) = (((ZSHAPE_NAME_0(:) .EQ. "PLA" .AND.  ZCIS_SHAPE_V(:,2) .GT. ZCTMIN(4) .AND. &
                                                                 ZCIS_SHAPE_V(:,2) .LT. ZCIS_SHAPE_V(:,1)) .OR. &
                             (ZSHAPE_NAME_0(:) .EQ. "COL" .AND.  ZCIS_SHAPE_V(:,1) .GT. ZCTMIN(4) .AND. &
                                                                 ZCIS_SHAPE_V(:,1) .LT. ZCIS_SHAPE_V(:,2))) .AND. &
                             ZDEP(1:IMICRO) .GT. 0.0)
!
        IHABIT = COUNT( GHABIT)
!
        IF (IHABIT >= 1) THEN
          ALLOCATE(ZZTM_V(IHABIT))
          ALLOCATE(ZZT0_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_M_V(IHABIT))
          ALLOCATE(ZSHAPE_NAME_0_V(IHABIT))
          ALLOCATE(ZCIS_SHAPE_VV(IHABIT,NB_CRYSTAL_SHAPE))
          ALLOCATE(ZDEP_V(IHABIT))
!
!          ALLOCATE(ZDUM(IMICRO,NB_CRYSTAL_SHAPE))
!          ZDUM(:,:) = ZCIS_SHAPE_V(:,:)
! 
          ZZTM_V(:) = PACK(ZZTM(:),GHABIT(:))
          ZZT0_V(:) = PACK(ZZT0(:),GHABIT(:))
          ZSHAPE_NAME_M_V(:) = PACK(ZSHAPE_NAME_M(:),GHABIT(:))
          ZSHAPE_NAME_0_V(:) = PACK(ZSHAPE_NAME_0(:),GHABIT(:))
          ZDEP_V(:) = PACK(ZDEP(:),GHABIT(:))
!
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCIS_SHAPE_VV(:,JSH) = PACK(ZCIS_SHAPE_V(:,JSH),GHABIT(:))
          END DO
!
          CALL LIMA_CHANGE_SHAPE(ZZTM_V, ZZT0_V, ZSHAPE_NAME_M_V, &
                                 ZSHAPE_NAME_0_V, ZCIS_SHAPE_VV,ZCTMIN,ZDEP_V)
!
!   Unpack variables
          ALLOCATE(ZW_V(IMICRO, NB_CRYSTAL_SHAPE))
!
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZW_V(:,JSH) = ZCIS_SHAPE_V(:,JSH)
            ZCIS_SHAPE_V(:,JSH) = UNPACK(ZCIS_SHAPE_VV(:,JSH),MASK=GHABIT(:),FIELD=ZW_V(:,JSH))
          END DO
!
          DEALLOCATE(ZCIS_SHAPE_VV)
          DEALLOCATE(ZZTM_V)
          DEALLOCATE(ZZT0_V)
          DEALLOCATE(ZSHAPE_NAME_M_V)
          DEALLOCATE(ZSHAPE_NAME_0_V)
          DEALLOCATE(ZW_V)
          DEALLOCATE(ZDEP_V)
        END IF ! ihabit
        DEALLOCATE(GHABIT)
      END IF ! lhabit
!
!  Implicit ice crystal sublimation if ice saturated conditions are not met
!
      ZZT(:) = ( ZTHS_V(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
      ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
      ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph   
      ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
      ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
      WHERE( ZRVS_V(:)*ZDT<ZRVSATI(:) )
        ZZW(:)  = ZRVS_V(:) + ZRIS_V(:)
        ZRVS_V(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
        ZTHS_V(:) = ZTHS_V(:) + ( MAX( 0.0,ZZW(:)-ZRVS_V(:) )-ZRIS_V(:) ) &
                                              * ZLSFACT(:) / ZEXNREF(:)
        ZRIS_V(:) = MAX( 0.0,ZZW(:)-ZRVS_V(:) )
      END WHERE
!                    
      ZW(:,:,:) = ZRVS(:,:,:)
      ZRVS(:,:,:) = UNPACK( ZRVS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = ZRCS(:,:,:)
      ZRCS(:,:,:) = UNPACK( ZRCS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = ZRIS(:,:,:)
      ZRIS(:,:,:) = UNPACK( ZRIS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
      ZW(:,:,:) = PTHS(:,:,:)
      PTHS(:,:,:) = UNPACK( ZTHS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
      IF (LCRYSTAL_SHAPE) THEN
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZW(:,:,:) = ZCIS_SHAPE(:,:,:,JSH)
          ZCIS_SHAPE(:,:,:,JSH) = UNPACK(ZCIS_SHAPE_V(:,JSH),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
        END DO
      END IF
!
      DEALLOCATE(ZRCT_V)
      DEALLOCATE(ZRIT_V)
      DEALLOCATE(ZCCT_V)
      DEALLOCATE(ZCIT_TOT_V)
      DEALLOCATE(ZRVS_V)
      DEALLOCATE(ZRCS_V)
      DEALLOCATE(ZRIS_V)
      DEALLOCATE(ZCCS_V)
      DEALLOCATE(ZCIS_TOT_V)
      DEALLOCATE(ZTHS_V)
      IF (ALLOCATED (ZITI_2D))       DEALLOCATE(ZITI_2D)
      IF (ALLOCATED (ZZW_2D))        DEALLOCATE(ZZW_2D)
      IF (ALLOCATED(ZDEP_2D))        DEALLOCATE(ZDEP_2D)
      IF (ALLOCATED(ZCND_2D))        DEALLOCATE(ZCND_2D)
      IF (ALLOCATED(ZLBDAI_SHAPE_T)) DEALLOCATE(ZLBDAI_SHAPE_T)
      IF (ALLOCATED(ZLBDAI_SHAPE_S)) DEALLOCATE(ZLBDAI_SHAPE_S)
      IF (ALLOCATED(ISHAPE_MAX_T))   DEALLOCATE(ISHAPE_MAX_T)
      IF (ALLOCATED(ZSHAPE_NAME_M))  DEALLOCATE(ZSHAPE_NAME_M)
      IF (ALLOCATED(ZSHAPE_NAME_0))  DEALLOCATE(ZSHAPE_NAME_0)
      DEALLOCATE(ZRHODREF)
      DEALLOCATE(ZZT)
      DEALLOCATE(ZPRES)
      DEALLOCATE(ZEXNREF)
      DEALLOCATE(ZZCPH)
      DEALLOCATE(ZZW)
      DEALLOCATE(ZLVFACT)
      DEALLOCATE(ZLSFACT)
      DEALLOCATE(ZRVSATW)
      DEALLOCATE(ZRVSATI)
      DEALLOCATE(ZRVSATW_PRIME)
      DEALLOCATE(ZRVSATI_PRIME)
      DEALLOCATE(ZDELTW)
      DEALLOCATE(ZDELTI)
      DEALLOCATE(ZAW)
      DEALLOCATE(ZAI)
      DEALLOCATE(ZCJ)
      DEALLOCATE(ZKA)
      DEALLOCATE(ZDV)
      DEALLOCATE(ZITW)
      IF (ALLOCATED(ZITI)) DEALLOCATE(ZITI)
      DEALLOCATE(ZAWW)
      DEALLOCATE(ZAIW)
      DEALLOCATE(ZAWI)
      DEALLOCATE(ZAII)
      DEALLOCATE(ZFACT)
      DEALLOCATE(ZDELT1)
      DEALLOCATE(ZDELT2)
      DEALLOCATE(ZCND)
      IF (ALLOCATED(ZDEP)) DEALLOCATE(ZDEP)
    END IF ! IMICRO
!
  END IF ! OSUBG_COND
!
! call budget when ice crystals change shape
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_SV) THEN
      IF (LCRYSTAL_SHAPE .AND. LHABIT_CHANGE) THEN
        DO JSH = 1, NB_CRYSTAL_SHAPE
          CALL BUDGET (ZCIS_SHAPE(:,:,:,JSH) * PRHODJ(:,:,:),12+NSV_LIMA_NI+JSH-1,'IHAB_BU_RSV')
        END DO
      END IF
    END IF
  END IF
!
! full sublimation of the cloud ice crystals if there are few
!
  ZMASK(:,:,:) = 0.0
  ZW(:,:,:) = 0.
!
  IF (LCRYSTAL_SHAPE) ZCIS_TOT(:,:,:) = SUM(ZCIS_SHAPE, DIM=4)
!
  WHERE (ZRIS(:,:,:) <= ZRTMIN(4) .OR. ZCIS_TOT(:,:,:) <= ZCTMIN(4)) 
    ZRVS(:,:,:) = ZRVS(:,:,:) + ZRIS(:,:,:) 
    PTHS(:,:,:) = PTHS(:,:,:) - ZRIS(:,:,:)*ZLS(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
    ZRIS(:,:,:) = 0.0
    ZW(:,:,:)   = MAX(ZCIS_TOT(:,:,:),0.)
    ZCIS_TOT(:,:,:) = 0.0
  END WHERE
!
  IF (LCRYSTAL_SHAPE) THEN
    DO JSH = 1, NB_CRYSTAL_SHAPE
      WHERE (ZCIS_TOT(:,:,:) .EQ. 0.0)
        ZCIS_SHAPE(:,:,:,JSH) = 0.0
      END WHERE
    END DO
  END IF
!
  IF (LCOLD .AND. (NMOD_IFN .GE. 1 .OR. NMOD_IMM .GE. 1)) THEN
    ZW1(:,:,:) = 0.
    IF (NMOD_IFN .GE. 1) ZW1(:,:,:) = ZW1(:,:,:) + SUM(ZINS,DIM=4)
    IF (NMOD_IMM .GE. 1) ZW1(:,:,:) = ZW1(:,:,:) + SUM(ZNIS,DIM=4)
    ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
    ZW2(:,:,:) = 0.
    WHERE ( ZW(:,:,:) > 0. )
      ZMASK(:,:,:) = 1.0
      ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
    ENDWHERE
  END IF
!
  IF (LCOLD .AND. NMOD_IFN.GE.1) THEN
    DO JMOD_IFN = 1, NMOD_IFN
      ZIFS(:,:,:,JMOD_IFN) = ZIFS(:,:,:,JMOD_IFN) +                    &
           ZMASK(:,:,:) * ZINS(:,:,:,JMOD_IFN) * ZW2(:,:,:)
      ZINS(:,:,:,JMOD_IFN) = ZINS(:,:,:,JMOD_IFN) -                    &
           ZMASK(:,:,:) * ZINS(:,:,:,JMOD_IFN) * ZW2(:,:,:)
      ZINS(:,:,:,JMOD_IFN) = MAX( 0.0 , ZINS(:,:,:,JMOD_IFN) )
    ENDDO
  END IF
!
  IF (LCOLD .AND. NMOD_IMM.GE.1) THEN
    JMOD_IMM = 0
    DO JMOD = 1, NMOD_CCN
      IF (NIMM(JMOD) == 1) THEN 
        JMOD_IMM = JMOD_IMM + 1 
        IF (LWARM .AND. NMOD_CCN.GE.1 ) THEN
          ZNAS(:,:,:,JMOD)     = ZNAS(:,:,:,JMOD) +                     &
               ZMASK(:,:,:) * ZNIS(:,:,:,JMOD_IMM) * ZW2(:,:,:)
        END IF
        ZNIS(:,:,:,JMOD_IMM) = ZNIS(:,:,:,JMOD_IMM) -                 &
             ZMASK(:,:,:) * ZNIS(:,:,:,JMOD_IMM) * ZW2(:,:,:)
        ZNIS(:,:,:,JMOD_IMM) = MAX( 0.0 , ZNIS(:,:,:,JMOD_IMM) )
      END IF
    ENDDO
  END IF
!
! complete evaporation of the cloud droplets if there are few
!
  ZMASK(:,:,:) = 0.0
  ZW(:,:,:) = 0.
  WHERE (ZRCS(:,:,:) <= ZRTMIN(2) .OR. ZCCS(:,:,:) <= ZCTMIN(2)) 
    ZRVS(:,:,:) = ZRVS(:,:,:) + ZRCS(:,:,:) 
    PTHS(:,:,:) = PTHS(:,:,:) - ZRCS(:,:,:)*ZLV(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
    ZRCS(:,:,:) = 0.0
    ZW(:,:,:)   = MAX(ZCCS(:,:,:),0.)
    ZCCS(:,:,:) = 0.0
  END WHERE
!
  ZW1(:,:,:) = 0.
  IF (LWARM .AND. NMOD_CCN.GE.1) ZW1(:,:,:) = SUM(ZNAS,DIM=4)
  ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
  ZW2(:,:,:) = 0.
  WHERE ( ZW(:,:,:) > 0. )
    ZMASK(:,:,:) = 1.0
    ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
  ENDWHERE
!
  IF (LWARM .AND. NMOD_CCN.GE.1) THEN
    DO JMOD = 1, NMOD_CCN
      ZNFS(:,:,:,JMOD) = ZNFS(:,:,:,JMOD) +                           &
           ZMASK(:,:,:) * ZNAS(:,:,:,JMOD) * ZW2(:,:,:)
      ZNAS(:,:,:,JMOD) = ZNAS(:,:,:,JMOD) -                           &
           ZMASK(:,:,:) * ZNAS(:,:,:,JMOD) * ZW2(:,:,:)
      ZNAS(:,:,:,JMOD) = MAX( 0.0 , ZNAS(:,:,:,JMOD) )
    ENDDO
  END IF
!
  IF (LSCAV .AND. LAERO_MASS) ZMAS(:,:,:) = ZMAS(:,:,:) * (1-ZMASK(:,:,:))
!
!  end of the iterative loop
!
END DO
!
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
!
!*       5.2    compute the cloud fraction PCLDFR (binary !!!!!!!)
!
IF ( .NOT. OSUBG_COND ) THEN
!  WHERE (ZRCS(:,:,:) + ZRIS(:,:,:) + ZRSS(:,:,:) > 1.E-12 / ZDT)
  WHERE (ZRCS(:,:,:) + ZRIS(:,:,:)  > 1.E-12 / ZDT)
    ZW(:,:,:)  = 1.
  ELSEWHERE
    ZW(:,:,:)  = 0. 
  ENDWHERE
  IF ( SIZE(PSRCS,3) /= 0 ) THEN
    PSRCS(:,:,:) = ZW(:,:,:) 
  END IF
END IF
!
PCLDFR(:,:,:) = ZW(:,:,:)
!
IF ( OCLOSE_OUT ) THEN
  TZFIELD%CMNHNAME   = 'NEB'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'NEB'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_NEB'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZW)
END IF
!
!
!*       6.  SAVE CHANGES IN PRS AND PSVS
!            ----------------------------
!
!
! Prepare 3D water mixing ratios
PRS(:,:,:,1) = ZRVS(:,:,:)
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = ZRCS(:,:,:)
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = ZRRS(:,:,:)
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = ZRIS(:,:,:)
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = ZRSS(:,:,:)
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = ZRGS(:,:,:)
!
! Prepare 3D number concentrations
!
IF ( LWARM ) PSVS(:,:,:,NSV_LIMA_NC) = ZCCS(:,:,:)
!
IF ( LCOLD ) THEN
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    PSVS(:,:,:,NSV_LIMA_NI) = ZCIS_TOT(:,:,:)
  ELSE
    DO JSH = 1, NB_CRYSTAL_SHAPE
      PSVS(:,:,:,NSV_LIMA_NI+JSH-1) = ZCIS_SHAPE(:,:,:,JSH)
   END DO
  END IF
END IF
!
IF ( LSCAV .AND. LAERO_MASS ) PSVS(:,:,:,NSV_LIMA_SCAVMASS) = ZMAS(:,:,:)
! 
IF ( LWARM .AND. NMOD_CCN .GE. 1 ) THEN
  PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = ZNFS(:,:,:,:)
  PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = ZNAS(:,:,:,:)
END IF
!
IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
  PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = ZIFS(:,:,:,:)
  PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = ZINS(:,:,:,:)
END IF
!
IF ( LCOLD .AND. NMOD_IMM .GE. 1 ) THEN
  PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = ZNIS(:,:,:,:)
END IF
!
! write SSI in LFI
!
IF ( OCLOSE_OUT ) THEN
  ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
  ZW(:,:,:) = EXP( XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:) ) )
  ZW1(:,:,:)= 2.0*PPABST(:,:,:)-PPABSM(:,:,:)
  ZW(:,:,:) = ZRVT(:,:,:)*( ZW1(:,:,:)-ZW(:,:,:) ) / ( (XMV/XMD) * ZW(:,:,:) ) - 1.0

  TZFIELD%CMNHNAME   = 'SSI'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'SSI'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_SSI'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZW)
END IF
!
!
!*       7.  STORE THE BUDGET TERMS
!            ----------------------
!
!
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:) * PRHODJ(:,:,:),4,'CEDS_BU_RTH')
  IF (LBUDGET_RV) CALL BUDGET (ZRVS(:,:,:) * PRHODJ(:,:,:),6,'CEDS_BU_RRV')
  IF (LBUDGET_RC) CALL BUDGET (ZRCS(:,:,:) * PRHODJ(:,:,:),7,'CEDS_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET (ZRIS(:,:,:) * PRHODJ(:,:,:),9,'CEDS_BU_RRI')
  IF (LBUDGET_SV) THEN
    CALL BUDGET (ZCCS(:,:,:) * PRHODJ(:,:,:),12+NSV_LIMA_NC,'CEDS_BU_RSV') ! RCC
    IF (NMOD_CCN .GE. 1) THEN
      DO JL = 1, NMOD_CCN
        CALL BUDGET (ZNFS(:,:,:,JL)*PRHODJ(:,:,:),12+NSV_LIMA_CCN_FREE+JL-1,'CEDS_BU_RSV') ! RCC
      END DO
    END IF
    IF (NMOD_IFN .GE. 1) THEN
      DO JL = 1, NMOD_IFN
        CALL BUDGET (ZIFS(:,:,:,JL)*PRHODJ(:,:,:),12+NSV_LIMA_IFN_FREE+JL-1,'CEDS_BU_RSV') ! RCC
      END DO
    END IF
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET (ZCIS_TOT(:,:,:) * PRHODJ(:,:,:),12+NSV_LIMA_NI,'CEDS_BU_RSV')
    ELSE
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET (ZCIS_SHAPE(:,:,:,JSH) * PRHODJ(:,:,:),12+NSV_LIMA_NI+JSH-1,'CEDS_BU_RSV') 
      END DO
    END IF
  END IF
END IF
!
IF (ALLOCATED(ZNFS)) DEALLOCATE(ZNFS)
IF (ALLOCATED(ZNAS)) DEALLOCATE(ZNAS)
IF (ALLOCATED(ZIFS)) DEALLOCATE(ZIFS)
IF (ALLOCATED(ZINS)) DEALLOCATE(ZINS)
IF (ALLOCATED(ZNIS)) DEALLOCATE(ZNIS)
!
IF (ALLOCATED(ZCIS_SHAPE))   DEALLOCATE(ZCIS_SHAPE)
IF (ALLOCATED(ZCIT_SHAPE))   DEALLOCATE(ZCIT_SHAPE)
IF (ALLOCATED(ZRIS_SHAPE_V)) DEALLOCATE(ZRIS_SHAPE_V)
IF (ALLOCATED(ZRIT_SHAPE_V)) DEALLOCATE(ZRIT_SHAPE_V)
IF (ALLOCATED(ZCIS_SHAPE_V)) DEALLOCATE(ZCIS_SHAPE_V)
IF (ALLOCATED(ZCIT_SHAPE_V)) DEALLOCATE(ZCIT_SHAPE_V)
IF (ALLOCATED(ZCSHAPE_S))    DEALLOCATE(ZCSHAPE_S)
IF (ALLOCATED(ZCSHAPE_T))    DEALLOCATE(ZCSHAPE_T)
IF (ALLOCATED(ZCSHAPE_T_V))  DEALLOCATE(ZCSHAPE_T_V)
IF (ALLOCATED(ZCSHAPE_S_V))  DEALLOCATE(ZCSHAPE_S_V)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ADJUST
