!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
       MODULE MODI_LIMA_MIXED
!      ######################
!
INTERFACE
      SUBROUTINE LIMA_MIXED (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI, &
                             KRR, PZZ, PRHODJ,                    &
                             PRHODREF, PEXNREF, PPABST, PW_NU,    &
                             PTHM, PPABSM,                        &
                             PTHT, PRT, PSVT,                     &
                             PTHS, PRS, PSVS)
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                      ! integration for  ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM  ! abs. pressure at time t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Concentrations at t 

!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS    ! Concentrations source
!
END SUBROUTINE LIMA_MIXED
END INTERFACE
END MODULE MODI_LIMA_MIXED
!
!     #######################################################################
      SUBROUTINE LIMA_MIXED (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI, &
                             KRR, PZZ, PRHODJ,                    &
                             PRHODREF, PEXNREF, PPABST, PW_NU,    &
                             PTHM, PPABSM,                        &
                             PTHT, PRT, PSVT,                     &
                             PTHS, PRS, PSVS                      )
!     #######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    microphysical processes
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!      Most of the parameterizations come from the ICE3 scheme, described in
!!    the MESO-NH scientific documentation.
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      T. Hoarau  * LACy *  jui. 2016   add budgets for CIBU
!!      M. Claeys  * LACy *  mar. 2019   add ice crystal shapes
!!      J.-P. Pinty * LA *   jui. 2019   add RDSF
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,       ONLY : JPHEXT, JPVEXT
USE MODD_CST,              ONLY : XP00, XRD, XRV, XMV, XMD, XCPD, XCPV,       &
                                  XCL, XCI, XTT, XLSTT, XLVTT,                &
                                  XALPI, XBETAI, XGAMI
USE MODD_PARAM_LIMA,       ONLY : NMOD_IFN, XRTMIN, XCTMIN, LWARM, LCOLD,     &
                                  NMOD_CCN, NMOD_IMM, LRAIN, LSNOW, LHAIL,    &
                                  LCIBU, LCRYSTAL_SHAPE, NB_CRYSTAL_SHAPE,    &
                                  LRDSF
USE MODD_PARAM_LIMA_WARM,  ONLY : XLBC, XLBEXC, XLBR, XLBEXR
USE MODD_PARAM_LIMA_COLD,  ONLY : XLBI, XLBEXI, XLBS, XLBEXS, XSCFAC,         &
                                  XLBI_SHAPE, XLBEXI_SHAPE
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XLBH, XLBEXH
!USE MODD_BUDGET,           ONLY : LBU_ENABLE, NBUMOD
!
USE MODD_NSV
!
USE MODD_BUDGET
USE MODI_BUDGET
!
USE MODI_LIMA_FUNCTIONS,   ONLY : COUNTJV
USE MODI_LIMA_MIXED_SLOW_PROCESSES
USE MODI_LIMA_MIXED_FAST_PROCESSES
USE MODI_COMPUTE_LBDA_SHAPE
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                      ! integration for  ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM  ! abs. pressure at time t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Concentrations at t 

!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS    ! Concentrations source
!
!*       0.2   Declarations of local variables :
!
!3D microphysical variables
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))  &
                                    :: ZRVT,    & ! Water vapor m.r. at t 
                                       ZRCT,    & ! Cloud water m.r. at t 
                                       ZRRT,    & ! Rain water m.r. at t 
                                       ZRIT,    & ! Cloud ice m.r. at t 
                                       ZRST,    & ! Snow/aggregate m.r. at t 
                                       ZRGT,    & ! Graupel m.r. at t 
                                       ZRHT,    & ! Hail m.r. at t 
                                       !
                                       ZRVS,    & ! Water vapor m.r. source
                                       ZRCS,    & ! Cloud water m.r. source
                                       ZRRS,    & ! Rain water m.r. source
                                       ZRIS,    & ! Pristine ice m.r. source
                                       ZRSS,    & ! Snow/aggregate m.r. source
                                       ZRGS,    & ! Graupel m.r. source
                                       ZRHS,    & ! Hail m.r. source
                                       !
                                       ZCCT,    & ! Cloud water C. at t
                                       ZCRT,    & ! Rain water C. at t
                                       ZCIT,    & ! Ice crystal (tot) C. at t
                                       !
                                       ZCCS,    & ! Cloud water C. source
                                       ZCRS,    & ! Rain water C. source
                                       ZCIS       ! Ice crystal (tot) C. source
    
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZONEOVER_VAR ! for optimization
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCIS_SHAPE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCIT_SHAPE
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZNFS     ! CCN C. available source
                                                  !used as Free ice nuclei for
                                                  !HOMOGENEOUS nucleation of haze
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZNAS     ! Cloud  C. nuclei C. source
                                                  !used as Free ice nuclei for
                                                  !IMMERSION freezing
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZIFS     ! Free ice nuclei C. source 
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZINS     ! Activated ice nuclei C. source
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZNIS     ! Activated ice nuclei C. source
                                                  !for IMMERSION
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZNHS     ! Hom. freezing of CCN
!
! Replace PACK
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GMICRO
INTEGER :: IMICRO
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
! Packed microphysical variables
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVT_V    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCT_V    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRT_V    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIT_V    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRST_V    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGT_V    ! Graupel m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHT_V    ! Hail m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCCT_V    ! Cloud water conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCRT_V    ! Rain water conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCIT_V    ! Pristine ice conc. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCIS_V    ! Pristine ice conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCIT_SHAPE_V  ! Pristine ice conc. for each shape at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCIS_SHAPE_V  ! Pristine ice conc. for each shape source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCSHAPE_S_V    ! Vectorized
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLBDAI_SHAPE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCSHAPE_S  ! Conc. ratio for each ice crystal shape
!
INTEGER, DIMENSION (:), ALLOCATABLE :: ISHAPE_MAX ! index of the dominant shape concentration 
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVS_V    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCS_V    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRS_V    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIS_V    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRSS_V    ! Snow/aggregate m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGS_V    ! Graupel m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHS_V    ! Hail m.r. source
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZTHS_V    ! Theta source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCCS_V    ! Cloud water conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCRS_V    ! Rain water conc. source
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFS_V    ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINS_V    ! Nucleated Ice nuclei conc. source
!
! Other packed variables
REAL, DIMENSION(:), ALLOCATABLE &
                   :: ZRHODREF, & ! RHO Dry REFerence
                      ZRHODJ,   & ! RHO times Jacobian
                      ZZT,      & ! Temperature
                      ZPRES,    & ! Pressure
                      ZEXNREF,  & ! EXNer Pressure REFerence
                      ZZW,      & ! Work array
                      ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                      ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                      ZSSI,     & ! Supersaturation over ice
                      ZLBDAC,   & ! Slope parameter of the cloud droplet distr.
                      ZLBDAR,   & ! Slope parameter of the raindrop  distr.
                      ZLBDAI,   & ! Slope parameter of the ice crystal distr.
                      ZLBDAS,   & ! Slope parameter of the aggregate distr.
                      ZLBDAG,   & ! Slope parameter of the graupel   distr.
                      ZLBDAH,   & ! Slope parameter of the hail   distr.
                      ZAI,      & ! Thermodynamical function
                      ZCJ,      & ! used to compute the ventilation coefficient
                      ZKA,      & ! Thermal conductivity of the air
                      ZDV         ! Diffusivity of water vapor in the air
!
! 3D Temperature
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZT, ZW
!
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
INTEGER :: JMOD_IFN                            ! Loop index 
INTEGER :: JSH                                 ! Loop index for ice crystal shapes
!
!-------------------------------------------------------------------------------
!
!
!*       0.     3D MICROPHYSCAL VARIABLES
!	        -------------------------
!
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
ZRHT(:,:,:) = 0.
ZRHS(:,:,:) = 0.
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
IF ( KRR .GE. 7 ) ZRHT(:,:,:) = PRT(:,:,:,7)
IF ( KRR .GE. 7 ) ZRHS(:,:,:) = PRS(:,:,:,7)
!
! Prepare 3D number concentrations
ZCCT(:,:,:) = 0.
ZCRT(:,:,:) = 0.
ZCCS(:,:,:) = 0.
ZCRS(:,:,:) = 0.
ZCIT(:,:,:) = 0.
ZCIS(:,:,:) = 0.
!
IF ( LWARM ) ZCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC) 
IF ( LWARM .AND. LRAIN ) ZCRT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NR)
!
IF ( LWARM ) ZCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) ZCRS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NR)
!
IF ( LCOLD ) THEN
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    ZCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
    ZCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
!++cb++ 4/02/20 pour pouvoir passer ces variables en arguments sans pb
    ALLOCATE(ZCIS_SHAPE(0,0,0,0))
!--cb--
  ELSE
    ALLOCATE(ZCIT_SHAPE(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NB_CRYSTAL_SHAPE))
    ALLOCATE(ZCIS_SHAPE(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NB_CRYSTAL_SHAPE))
    ALLOCATE(ZCSHAPE_S(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NB_CRYSTAL_SHAPE))
    ALLOCATE(ZONEOVER_VAR(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)))
    ZCSHAPE_S(:,:,:,:) = 0.0
! Concentration per shape
    DO JSH = 1, NB_CRYSTAL_SHAPE
      ZCIT_SHAPE(:,:,:,JSH) = PSVT(:,:,:,NSV_LIMA_NI+JSH-1)
      ZCIS_SHAPE(:,:,:,JSH) = PSVS(:,:,:,NSV_LIMA_NI+JSH-1)
    END DO
! Total concentration
    ZCIT(:,:,:) = SUM(ZCIT_SHAPE, DIM=4)
    ZCIS(:,:,:) = SUM(ZCIS_SHAPE, DIM=4)
!
! Calcul du rapport de concentration de chaque forme par rapport a la
! concentration totale
    WHERE (ZCIS(:,:,:) .GT. 0.0) ZONEOVER_VAR(:,:,:) = 1.0 / ZCIS(:,:,:)
    DO JSH = 1, NB_CRYSTAL_SHAPE
      WHERE ((ZCIS(:,:,:) .GT. 0.0) .AND. (ZCIS_SHAPE(:,:,:,JSH) .GT. 0.0))
        ZCSHAPE_S(:,:,:,JSH) = MIN(ZCIS_SHAPE(:,:,:,JSH) * ZONEOVER_VAR(:,:,:), 1.0)
      END WHERE
    END DO
    DEALLOCATE(ZONEOVER_VAR)
  END IF
END IF
!
IF ( NMOD_CCN .GE. 1 ) THEN
   ALLOCATE( ZNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( ZNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ZNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   ZNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
ELSE
   ALLOCATE( ZNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( ZNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ZNFS(:,:,:,:) = 0.
   ZNAS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   ALLOCATE( ZIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ALLOCATE( ZINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ZIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
   ZINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
ELSE
   ALLOCATE( ZIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( ZINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ZIFS(:,:,:,:) = 0.
   ZINS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   ALLOCATE( ZNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
   ZNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
ELSE
   ALLOCATE( ZNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ZNIS(:,:,:,:) = 0.0
END IF
!
IF ( OHHONI ) THEN
   ALLOCATE( ZNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   ZNHS(:,:,:) = PSVS(:,:,:,NSV_LIMA_HOM_HAZE)
ELSE
   ALLOCATE( ZNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   ZNHS(:,:,:) = 0.0
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Pack variables, computations only where necessary
!	        -------------------------------------------------
!
! Physical domain
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
! Temperature
ZT(:,:,:)  = PTHT(:,:,:) * ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
! Looking for regions where computations are necessary
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) = ZRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR. &
                                  ZRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3) .OR. &
                                  ZRIT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(4) .OR. &
                                  ZRST(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(5) .OR. &
                                  ZRGT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(6) .OR. &
                                  ZRHT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(7)
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
!++je++ 19/02/20 test //
!IF( IMICRO >= 1 ) THEN
IF ( IMICRO >= 0 ) THEN
!--je--
!
   ALLOCATE(ZRVT_V(IMICRO)) 
   ALLOCATE(ZRCT_V(IMICRO))
   ALLOCATE(ZRRT_V(IMICRO))  
   ALLOCATE(ZRIT_V(IMICRO)) 
   ALLOCATE(ZRST_V(IMICRO)) 
   ALLOCATE(ZRGT_V(IMICRO))  
   ALLOCATE(ZRHT_V(IMICRO))  
   !
   ALLOCATE(ZCCT_V(IMICRO)) 
   ALLOCATE(ZCRT_V(IMICRO)) 
   !
   ALLOCATE(ZRVS_V(IMICRO))  
   ALLOCATE(ZRCS_V(IMICRO)) 
   ALLOCATE(ZRRS_V(IMICRO)) 
   ALLOCATE(ZRIS_V(IMICRO))
   ALLOCATE(ZRSS_V(IMICRO))
   ALLOCATE(ZRGS_V(IMICRO)) 
   ALLOCATE(ZRHS_V(IMICRO)) 
   ALLOCATE(ZTHS_V(IMICRO))
   !
   ALLOCATE(ZCCS_V(IMICRO)) 
   ALLOCATE(ZCRS_V(IMICRO)) 
   ALLOCATE(ZIFS_V(IMICRO,NMOD_IFN))
   ALLOCATE(ZINS_V(IMICRO,NMOD_IFN))
   !
   ALLOCATE(ZRHODREF(IMICRO)) 
   ALLOCATE(ZZT(IMICRO)) 
   ALLOCATE(ZPRES(IMICRO)) 
   ALLOCATE(ZEXNREF(IMICRO))
   ALLOCATE(ZCIT_V(IMICRO)) 
   ALLOCATE(ZCIS_V(IMICRO))
!
   IF (LCRYSTAL_SHAPE) THEN
     ALLOCATE(ZCIT_SHAPE_V(IMICRO,NB_CRYSTAL_SHAPE)) 
     ALLOCATE(ZCIS_SHAPE_V(IMICRO,NB_CRYSTAL_SHAPE)) 
     ALLOCATE(ZCSHAPE_S_V(IMICRO,NB_CRYSTAL_SHAPE))
   ELSE
     ALLOCATE(ZCIT_SHAPE_V(0,0))
     ALLOCATE(ZCIS_SHAPE_V(0,0))
     ALLOCATE(ZCSHAPE_S_V(0,0))
   END IF
!
   DO JL=1,IMICRO   
      ZRVT_V(JL) = ZRVT(I1(JL),I2(JL),I3(JL))
      ZRCT_V(JL) = ZRCT(I1(JL),I2(JL),I3(JL))
      ZRRT_V(JL) = ZRRT(I1(JL),I2(JL),I3(JL))
      ZRIT_V(JL) = ZRIT(I1(JL),I2(JL),I3(JL))
      ZRST_V(JL) = ZRST(I1(JL),I2(JL),I3(JL))
      ZRGT_V(JL) = ZRGT(I1(JL),I2(JL),I3(JL))
      ZRHT_V(JL) = ZRHT(I1(JL),I2(JL),I3(JL))
      !
      ZCCT_V(JL) = ZCCT(I1(JL),I2(JL),I3(JL))
      ZCRT_V(JL) = ZCRT(I1(JL),I2(JL),I3(JL))
      !
      ZRVS_V(JL) = ZRVS(I1(JL),I2(JL),I3(JL))
      ZRCS_V(JL) = ZRCS(I1(JL),I2(JL),I3(JL))
      ZRRS_V(JL) = ZRRS(I1(JL),I2(JL),I3(JL))
      ZRIS_V(JL) = ZRIS(I1(JL),I2(JL),I3(JL))
      ZRSS_V(JL) = ZRSS(I1(JL),I2(JL),I3(JL))
      ZRGS_V(JL) = ZRGS(I1(JL),I2(JL),I3(JL))
      ZRHS_V(JL) = ZRHS(I1(JL),I2(JL),I3(JL))
      ZTHS_V(JL) = PTHS(I1(JL),I2(JL),I3(JL))
      !
      ZCCS_V(JL) = ZCCS(I1(JL),I2(JL),I3(JL))
      ZCRS_V(JL) = ZCRS(I1(JL),I2(JL),I3(JL))
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZCIT_V(JL) = ZCIT(I1(JL),I2(JL),I3(JL))
        ZCIS_V(JL) = ZCIS(I1(JL),I2(JL),I3(JL))
      ELSE
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZCIT_SHAPE_V(JL,JSH) = ZCIT_SHAPE(I1(JL),I2(JL),I3(JL),JSH)
          ZCIS_SHAPE_V(JL,JSH) = ZCIS_SHAPE(I1(JL),I2(JL),I3(JL),JSH)
          ZCSHAPE_S_V(JL,JSH) = ZCSHAPE_S(I1(JL),I2(JL),I3(JL),JSH)
        END DO
!        ZCIT_V(:) = SUM(ZCIT_SHAPE_V, DIM=2)
!        ZCIS_V(:) = SUM(ZCIS_SHAPE_V, DIM=2)
      END IF
!
      DO JMOD_IFN = 1, NMOD_IFN
        ZIFS_V(JL,JMOD_IFN) = ZIFS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
        ZINS_V(JL,JMOD_IFN) = ZINS(I1(JL),I2(JL),I3(JL),JMOD_IFN)
      ENDDO
      !
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
   ENDDO
!
   IF (LCRYSTAL_SHAPE) THEN
     ZCIT_V(:) = SUM(ZCIT_SHAPE_V, DIM=2)
     ZCIS_V(:) = SUM(ZCIS_SHAPE_V, DIM=2)
   END IF
!
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
      ALLOCATE(ZRHODJ(IMICRO))
      ZRHODJ(:) = PACK( PRHODJ(:,:,:),MASK=GMICRO(:,:,:) )
   END IF
!
! Atmospheric parameters 
!
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ALLOCATE(ZLVFACT(IMICRO))
   ALLOCATE(ZSSI(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
!
   ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT_V(:)+XCL*(ZRCT_V(:)+ZRRT_V(:)) &
        +XCI*(ZRIT_V(:)+ZRST_V(:)+ZRGT_V(:)) )
!
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
   ZSSI(:) = ZRVT_V(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                       ! Supersaturation over ice
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
!
! Thermodynamical function ZAI = A_i(T,P)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                         + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
! ZCJ = c^prime_j (in the ventilation factor)
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!
! Particle distribution parameters
!
  ALLOCATE(ZLBDAC(IMICRO)) 
  ALLOCATE(ZLBDAR(IMICRO))
  ALLOCATE(ZLBDAI(IMICRO)) 
  IF (LCRYSTAL_SHAPE) ALLOCATE(ZLBDAI_SHAPE(IMICRO,NB_CRYSTAL_SHAPE))
  ALLOCATE(ZLBDAS(IMICRO))
  ALLOCATE(ZLBDAG(IMICRO))
  ALLOCATE(ZLBDAH(IMICRO))
  ZLBDAC(:)  = 1.E10
  WHERE (ZRCT_V(:)>XRTMIN(2) .AND. ZCCT_V(:)>XCTMIN(2))
    ZLBDAC(:) = ( XLBC*ZCCT_V(:) / ZRCT_V(:) )**XLBEXC
  END WHERE
  ZLBDAR(:)  = 1.E10
  WHERE (ZRRT_V(:)>XRTMIN(3) .AND. ZCRT_V(:)>XCTMIN(3))
    ZLBDAR(:) = ( XLBR*ZCRT_V(:) / ZRRT_V(:) )**XLBEXR
  END WHERE
  ZLBDAI(:)  = 1.E10
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    WHERE (ZRIT_V(:) > XRTMIN(4) .AND. ZCIT_V(:) > XCTMIN(4))
      ZLBDAI(:) = ( XLBI*ZCIT_V(:) / ZRIT_V(:) )**XLBEXI
    END WHERE
  ELSE
    ALLOCATE(ISHAPE_MAX(IMICRO)) 
    ZLBDAI_SHAPE(:,:)= 1.E10
    ISHAPE_MAX(:) = MAXLOC(ZCSHAPE_S_V,DIM=2)
    CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX, IMICRO, PTSTEP, &
                            ZRIT_V, ZCIT_V, ZLBDAI_SHAPE)
  END IF  ! LCRYSTAL_SHAPE
  ZLBDAS(:)  = 1.E10
  WHERE (ZRST_V(:)>XRTMIN(5) )
    ZLBDAS(:) = XLBS*( ZRHODREF(:)*ZRST_V(:) )**XLBEXS
  END WHERE
  ZLBDAG(:)  = 1.E10
  WHERE (ZRGT_V(:)>XRTMIN(6) )
    ZLBDAG(:) = XLBG*( ZRHODREF(:)*ZRGT_V(:) )**XLBEXG
  END WHERE
  ZLBDAH(:)  = 1.E10
  WHERE (ZRHT_V(:)>XRTMIN(7) )
    ZLBDAH(:) = XLBH*( ZRHODREF(:)*ZRHT_V(:) )**XLBEXH
  END WHERE
! 
!-------------------------------------------------------------------------------
!
!
!*       2.     Compute the slow processes involving cloud water and graupel
!	        ------------------------------------------------------------
!
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    CALL LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,            &
                                   ZLSFACT, ZLVFACT, ZAI, ZCJ,             &
                                   ZRGT_V, ZCIT_V, ZCIT_SHAPE_V,       &
                                   ZRVS_V, ZRCS_V, ZRIS_V, ZRGS_V, ZTHS_V, &
                                   ZCCS_V, ZCIS_V, ZCIS_SHAPE_V,       &
                                   ZIFS_V, ZINS_V,                         &
                                   ZLBDAI, ZLBDAG,                         &
                                   GMICRO, PRHODJ, KMI,                    &
                                   PTHS, ZRVS, ZRCS, ZRIS, ZRGS,           &
                                   ZCCS, ZCIS, ZCIS_SHAPE)
  ELSE
    CALL LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,            &
                                   ZLSFACT, ZLVFACT, ZAI, ZCJ,             &
                                   ZRGT_V, ZCIT_V, ZCIT_SHAPE_V,       &
                                   ZRVS_V, ZRCS_V, ZRIS_V, ZRGS_V, ZTHS_V, &
                                   ZCCS_V, ZCIS_V, ZCIS_SHAPE_V,       &
                                   ZIFS_V, ZINS_V,                         &
                                   ZLBDAI, ZLBDAG,                         &
                                   GMICRO, PRHODJ, KMI,                    &
                                   PTHS, ZRVS, ZRCS, ZRIS, ZRGS,           &
                                   ZCCS, ZCIS, ZCIS_SHAPE,             &
                                   ZLBDAI_SHAPE)
  END IF
! 
!-------------------------------------------------------------------------------
!
!
!        3.     Compute the fast RS and RG processes
!   	        ------------------------------------
!
  IF (LSNOW) THEN
    CALL LIMA_MIXED_FAST_PROCESSES(ZRHODREF, ZZT, ZPRES, PTSTEP,               &
                                   ZLSFACT, ZLVFACT, ZKA, ZDV, ZCJ,            &
                                   ZRVT_V, ZRCT_V, ZRRT_V, ZRIT_V, ZRST_V,     &
                                   ZRGT_V, ZRHT_V, ZCCT_V, ZCRT_V, ZCIT_V, &
                                   ZRCS_V, ZRRS_V, ZRIS_V, ZRSS_V, ZRGS_V,     &
                                   ZRHS_V, ZTHS_V, ZCCS_V, ZCRS_V, ZCIS_V, &
                                   ZCIS_SHAPE_V, ZCSHAPE_S_V,                  &
                                   ZLBDAC, ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH,     &
                                   PRHODJ, GMICRO, KMI, PTHS,                  &
                                   ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS,         &
                                   ZCCS, ZCRS, ZCIS, ZCIS_SHAPE )                     
  END IF
!
!-------------------------------------------------------------------------------
!
!        4.     Unpack variables
!   	        ----------------
!
!
  ZW(:,:,:) = ZRVS(:,:,:)
  ZRVS(:,:,:) = UNPACK( ZRVS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRCS(:,:,:)
  ZRCS(:,:,:) = UNPACK( ZRCS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRRS(:,:,:)
  ZRRS(:,:,:) = UNPACK( ZRRS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRIS(:,:,:)
  ZRIS(:,:,:) = UNPACK( ZRIS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRSS(:,:,:)
  ZRSS(:,:,:) = UNPACK( ZRSS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRGS(:,:,:)
  ZRGS(:,:,:) = UNPACK( ZRGS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZRHS(:,:,:)
  ZRHS(:,:,:) = UNPACK( ZRHS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  ZW(:,:,:) = PTHS(:,:,:)
  PTHS(:,:,:) = UNPACK( ZTHS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  ZW(:,:,:) = ZCCS(:,:,:)
  ZCCS(:,:,:) = UNPACK( ZCCS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZCRS(:,:,:)
  ZCRS(:,:,:) = UNPACK( ZCRS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = ZCIS(:,:,:)
  ZCIS(:,:,:) = UNPACK( ZCIS_V(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  IF (LCRYSTAL_SHAPE) THEN
    DO JSH = 1, NB_CRYSTAL_SHAPE
      ZW(:,:,:) = ZCIS_SHAPE(:,:,:,JSH)
      ZCIS_SHAPE(:,:,:,JSH) = UNPACK( ZCIS_SHAPE_V(:,JSH),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
    END DO
  END IF
!
  DO JMOD_IFN = 1, NMOD_IFN
    ZW(:,:,:) = ZIFS(:,:,:,JMOD_IFN)
    ZIFS(:,:,:,JMOD_IFN) = UNPACK( ZIFS_V(:,JMOD_IFN),MASK=GMICRO(:,:,:), &
                                                  FIELD=ZW(:,:,:) )
    ZW(:,:,:) = ZINS(:,:,:,JMOD_IFN)
    ZINS(:,:,:,JMOD_IFN) = UNPACK( ZINS_V(:,JMOD_IFN),MASK=GMICRO(:,:,:), &
                                                  FIELD=ZW(:,:,:) )
  ENDDO
!
  DEALLOCATE(ZRVT_V) 
  DEALLOCATE(ZRCT_V)
  DEALLOCATE(ZRRT_V)  
  DEALLOCATE(ZRIT_V) 
  DEALLOCATE(ZRST_V) 
  DEALLOCATE(ZRGT_V)
  DEALLOCATE(ZRHT_V)
!  
  DEALLOCATE(ZCCT_V) 
  DEALLOCATE(ZCRT_V)
! 
  DEALLOCATE(ZRVS_V)  
  DEALLOCATE(ZRCS_V) 
  DEALLOCATE(ZRRS_V) 
  DEALLOCATE(ZRIS_V)
  DEALLOCATE(ZRSS_V)
  DEALLOCATE(ZRGS_V) 
  DEALLOCATE(ZRHS_V) 
  DEALLOCATE(ZTHS_V)
!
  DEALLOCATE(ZCCS_V) 
  DEALLOCATE(ZCRS_V) 
  DEALLOCATE(ZIFS_V)
  DEALLOCATE(ZINS_V)
!
  DEALLOCATE(ZCIS_V)  
  DEALLOCATE(ZCIT_V)
  IF (ALLOCATED(ZCIS_SHAPE_V)) DEALLOCATE(ZCIS_SHAPE_V)  
  IF (ALLOCATED(ZCIT_SHAPE_V)) DEALLOCATE(ZCIT_SHAPE_V)
  IF (ALLOCATED(ISHAPE_MAX))   DEALLOCATE(ISHAPE_MAX)
  IF (ALLOCATED(ZCSHAPE_S_V))  DEALLOCATE(ZCSHAPE_S_V)
  IF (ALLOCATED(ZLBDAI_SHAPE)) DEALLOCATE(ZLBDAI_SHAPE)
!
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZT) 
  DEALLOCATE(ZPRES) 
  DEALLOCATE(ZEXNREF)
!
  DEALLOCATE(ZZW)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZLVFACT)
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZAI)
  DEALLOCATE(ZCJ)
  DEALLOCATE(ZKA)
  DEALLOCATE(ZDV)
!
  DEALLOCATE(ZLBDAC) 
  DEALLOCATE(ZLBDAR)
  DEALLOCATE(ZLBDAI) 
  DEALLOCATE(ZLBDAS)
  DEALLOCATE(ZLBDAG)
  DEALLOCATE(ZLBDAH)
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) DEALLOCATE(ZRHODJ)
!
!
ELSE
!
! Advance the budget calls
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) THEN
      ZW(:,:,:) = PTHS(:,:,:)*PRHODJ(:,:,:)
      IF (LSNOW) CALL BUDGET (ZW,4,'DEPG_BU_RTH')
      CALL BUDGET (ZW,4,'IMLT_BU_RTH')
      CALL BUDGET (ZW,4,'BERFI_BU_RTH')
      IF (LSNOW) CALL BUDGET (ZW,4,'RIM_BU_RTH')
      IF (LSNOW .AND. LRAIN) CALL BUDGET (ZW,4,'ACC_BU_RTH')
      IF (LSNOW) CALL BUDGET (ZW,4,'CFRZ_BU_RTH')
      IF (LSNOW) CALL BUDGET (ZW,4,'WETG_BU_RTH')
      IF (LSNOW) CALL BUDGET (ZW,4,'DRYG_BU_RTH')
      IF (LSNOW) CALL BUDGET (ZW,4,'GMLT_BU_RTH')
      IF (LHAIL) CALL BUDGET (ZW,4,'WETH_BU_RTH')
      IF (LHAIL) CALL BUDGET (ZW,4,'HMLT_BU_RTH')
    ENDIF
    IF (LBUDGET_RV) THEN
      ZW(:,:,:) = ZRVS(:,:,:)*PRHODJ(:,:,:)
      IF (LSNOW) CALL BUDGET (ZW,6,'DEPG_BU_RRV')
    ENDIF
    IF (LBUDGET_RC) THEN
      ZW(:,:,:) = ZRCS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,7,'IMLT_BU_RRC')
      CALL BUDGET (ZW,7,'BERFI_BU_RRC')
      IF (LSNOW) CALL BUDGET (ZW,7,'RIM_BU_RRC')
      IF (LSNOW) CALL BUDGET (ZW,7,'WETG_BU_RRC')
      IF (LSNOW) CALL BUDGET (ZW,7,'DRYG_BU_RRC')
      IF (LHAIL) CALL BUDGET (ZW,7,'WETH_BU_RRC')
    ENDIF
    IF (LBUDGET_RR .AND. LRAIN) THEN
      ZW(:,:,:) = ZRRS(:,:,:)*PRHODJ(:,:,:)
      IF (LSNOW .AND. LRAIN) CALL BUDGET (ZW,8,'ACC_BU_RRR')
      IF (LSNOW) CALL BUDGET (ZW,8,'CFRZ_BU_RRR')
      IF (LSNOW) CALL BUDGET (ZW,8,'WETG_BU_RRR')
      IF (LSNOW) CALL BUDGET (ZW,8,'DRYG_BU_RRR')
      IF (LSNOW) CALL BUDGET (ZW,8,'GMLT_BU_RRR')
      IF (LHAIL) CALL BUDGET (ZW,8,'WETH_BU_RRR')
      IF (LHAIL) CALL BUDGET (ZW,8,'HMLT_BU_RRR')
    ENDIF
    IF (LBUDGET_RI) THEN
      ZW(:,:,:) = ZRIS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,9,'IMLT_BU_RRI')
      CALL BUDGET (ZW,9,'BERFI_BU_RRI')
      IF (LSNOW) CALL BUDGET (ZW,9,'HMS_BU_RRI')
      IF (LCIBU) CALL BUDGET (ZW,9,'CIBU_BU_RRI')
      IF (LSNOW) CALL BUDGET (ZW,9,'CFRZ_BU_RRI')
      IF (LRDSF) CALL BUDGET (ZW,9,'RDSF_BU_RRI')
      IF (LSNOW) CALL BUDGET (ZW,9,'WETG_BU_RRI')
      IF (LSNOW) CALL BUDGET (ZW,9,'DRYG_BU_RRI')
      IF (LSNOW) CALL BUDGET (ZW,9,'HMG_BU_RRI')
      IF (LHAIL) CALL BUDGET (ZW,9,'WETH_BU_RRI')
    ENDIF
    IF (LBUDGET_RS .AND. LSNOW) THEN
      ZW(:,:,:) = ZRSS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,10,'RIM_BU_RRS')
      CALL BUDGET (ZW,10,'HMS_BU_RRS')
      IF (LCIBU) CALL BUDGET (ZW,10,'CIBU_BU_RRS')
      IF (LRAIN) CALL BUDGET (ZW,10,'ACC_BU_RRS')
      CALL BUDGET (ZW,10,'CMEL_BU_RRS')
      CALL BUDGET (ZW,10,'WETG_BU_RRS')
      CALL BUDGET (ZW,10,'DRYG_BU_RRS')
      IF (LHAIL) CALL BUDGET (ZW,10,'WETH_BU_RRS')
    ENDIF
    IF (LBUDGET_RG .AND. LSNOW) THEN
      ZW(:,:,:) = ZRGS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,11,'DEPG_BU_RRG')
      CALL BUDGET (ZW,11,'RIM_BU_RRG')
      IF (LRAIN) CALL BUDGET (ZW,11,'ACC_BU_RRG')
      CALL BUDGET (ZW,11,'CMEL_BU_RRG')
      CALL BUDGET (ZW,11,'CFRZ_BU_RRG')
      IF (LRDSF) CALL BUDGET (ZW,11,'RDSF_BU_RRG')
      CALL BUDGET (ZW,11,'WETG_BU_RRG')
      CALL BUDGET (ZW,11,'DRYG_BU_RRG')
      CALL BUDGET (ZW,11,'HMG_BU_RRG')
      CALL BUDGET (ZW,11,'GMLT_BU_RRG')
      IF (LHAIL) CALL BUDGET (ZW,11,'WETH_BU_RRG')
      IF (LHAIL) CALL BUDGET (ZW,11,'COHG_BU_RRG')
    ENDIF
    IF (LBUDGET_RH .AND. LHAIL) THEN
      ZW(:,:,:) = ZRHS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,12,'WETG_BU_RRH')
      IF (LHAIL) CALL BUDGET (ZW,12,'WETH_BU_RRH')
      IF (LHAIL) CALL BUDGET (ZW,12,'COHG_BU_RRH')
      IF (LHAIL) CALL BUDGET (ZW,12,'HMLT_BU_RRH')
    ENDIF
    IF (LBUDGET_SV) THEN
      ZW(:,:,:) = ZCCS(:,:,:)*PRHODJ(:,:,:)
      CALL BUDGET (ZW,12+NSV_LIMA_NC,'IMLT_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NC,'RIM_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NC,'WETG_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NC,'DRYG_BU_RSV')
      IF (LHAIL) CALL BUDGET (ZW,12+NSV_LIMA_NC,'WETH_BU_RSV')
!
      ZW(:,:,:) = ZCRS(:,:,:)*PRHODJ(:,:,:)
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NR,'ACC_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NR,'CFRZ_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NR,'WETG_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NR,'DRYG_BU_RSV')
      IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NR,'GMLT_BU_RSV')
      IF (LHAIL) CALL BUDGET (ZW,12+NSV_LIMA_NR,'WETH_BU_RSV')
      IF (LHAIL) CALL BUDGET (ZW,12+NSV_LIMA_NR,'HMLT_BU_RSV')
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN      
        ZW(:,:,:) = ZCIS(:,:,:)*PRHODJ(:,:,:)
        CALL BUDGET (ZW,12+NSV_LIMA_NI,'IMLT_BU_RSV')
        IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NI,'HMS_BU_RSV')
        IF (LCIBU) CALL BUDGET (ZW,12+NSV_LIMA_NI,'CIBU_BU_RSV')
        IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NI,'CFRZ_BU_RSV')
        IF (LRDSF) CALL BUDGET (ZW,12+NSV_LIMA_NI,'RDSF_BU_RSV')
        IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NI,'WETG_BU_RSV')
        IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NI,'DRYG_BU_RSV')
        IF (LSNOW) CALL BUDGET (ZW,12+NSV_LIMA_NI,'HMG_BU_RSV')
        IF (LHAIL) CALL BUDGET (ZW,12+NSV_LIMA_NI,'WETH_BU_RSV')

      ELSE 
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZW(:,:,:) = ZCIS_SHAPE(:,:,:,JSH)*PRHODJ(:,:,:)
          CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'IMLT_BU_RSV')
          IF (LSNOW) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'HMS_BU_RSV')
          IF (LCIBU) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'CIBU_BU_RSV')
          IF (LSNOW) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'CFRZ_BU_RSV')
          IF (LRDSF) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'RDSF_BU_RSV')
          IF (LSNOW) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'WETG_BU_RSV')
          IF (LSNOW) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'DRYG_BU_RSV')
          IF (LSNOW) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'HMG_BU_RSV')
          IF (LHAIL) CALL BUDGET (ZW(:,:,:),12+NSV_LIMA_NI+JSH-1,'WETH_BU_RSV')
        END DO
      END IF
    ENDIF ! BUDGET_SV
  ENDIF ! NBUMOD
!
END IF ! IMICRO >= 1
!
!------------------------------------------------------------------------------
!
!
!*       5.    REPORT 3D MICROPHYSICAL VARIABLES IN PRS AND PSVS
!              -------------------------------------------------
!
PRS(:,:,:,1) = ZRVS(:,:,:)
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = ZRCS(:,:,:)
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = ZRRS(:,:,:)
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = ZRIS(:,:,:)
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = ZRSS(:,:,:)
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = ZRGS(:,:,:)
IF ( KRR .GE. 7 ) PRS(:,:,:,7) = ZRHS(:,:,:)
!
! Prepare 3D number concentrations
!
PSVS(:,:,:,NSV_LIMA_NC) = ZCCS(:,:,:)
IF ( LRAIN ) PSVS(:,:,:,NSV_LIMA_NR) = ZCRS(:,:,:)
IF (.NOT. LCRYSTAL_SHAPE) THEN
  PSVS(:,:,:,NSV_LIMA_NI) = ZCIS(:,:,:)
ELSE 
  DO JSH = 1, NB_CRYSTAL_SHAPE
    PSVS(:,:,:,NSV_LIMA_NI+JSH-1) = ZCIS_SHAPE(:,:,:,JSH)
  END DO
END IF
!
IF ( NMOD_CCN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = ZNFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = ZNAS(:,:,:,:)
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = ZIFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = ZINS(:,:,:,:)
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = ZNIS(:,:,:,:)
END IF
!
IF (ALLOCATED(ZNFS)) DEALLOCATE(ZNFS)
IF (ALLOCATED(ZNAS)) DEALLOCATE(ZNAS)
IF (ALLOCATED(ZIFS)) DEALLOCATE(ZIFS)
IF (ALLOCATED(ZINS)) DEALLOCATE(ZINS)
IF (ALLOCATED(ZNIS)) DEALLOCATE(ZNIS)
IF (ALLOCATED(ZNHS)) DEALLOCATE(ZNHS)
IF (ALLOCATED(ZCIT_SHAPE))   DEALLOCATE(ZCIT_SHAPE)
IF (ALLOCATED(ZCIS_SHAPE))   DEALLOCATE(ZCIS_SHAPE)
IF (ALLOCATED(ZCSHAPE_S))    DEALLOCATE(ZCSHAPE_S)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED
