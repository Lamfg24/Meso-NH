!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #####################
       MODULE MODI_LIMA_COLD
!      #####################
!
INTERFACE
      SUBROUTINE LIMA_COLD (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,           &
                           KRR, PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,               &
                           PTHM, PPABSM,                                   &
                           PTHT, PRT, PSVT,                                &
                           PTHS, PRS, PSVS,                                &
                           PINPRS, PINPRG, PINPRH)
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                                   ! for ice sedimendation
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
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
END SUBROUTINE LIMA_COLD
END INTERFACE
END MODULE MODI_LIMA_COLD
!
!     ######################################################################
      SUBROUTINE LIMA_COLD (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,           &
                           KRR, PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,               &
                           PTHM, PPABSM,                                   &
                           PTHT, PRT, PSVT,                                &
                           PTHS, PRS, PSVS,                                &
                           PINPRS, PINPRG, PINPRH)
!     ######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase 
!!    microphysical sources involving only primary ice and snow, except for 
!!    the sedimentation which also includes graupelns, and the homogeneous
!!    freezing of CCNs, cloud droplets and raindrops.
!!
!!
!!**  METHOD
!!    ------
!!      The nucleation of IFN is parameterized following either Meyers (1992)
!!    or Phillips (2008, 2013).
!!
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!
!!    REFERENCES
!!    ----------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!      Phillips et al., 2008: An empirical parameterization of heterogeneous
!!        ice nucleation for multiple chemical species of aerosols, J. Atmos. Sci. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_NSV
USE MODD_PARAM_LIMA
!
USE MODD_BUDGET
USE MODI_BUDGET
!
USE MODI_LIMA_COLD_SEDIMENTATION
USE MODI_LIMA_MEYERS
USE MODI_LIMA_PHILLIPS
USE MODI_LIMA_COLD_HOM_NUCL
USE MODI_LIMA_COLD_SLOW_PROCESSES
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                                   ! for ice sedimendation
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
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))  &
                                    :: ZRVT,    & ! Water vapor m.r. at t 
                                       ZRCT,    & ! Cloud water m.r. at t 
                                       ZRRT,    & ! Rain water m.r. at t 
                                       ZRIT,    & ! Cloud ice m.r. at t 
                                       ZRST,    & ! Snow/aggregate m.r. at t 
                                       ZRGT,    & ! Graupel m.r. at t 
                                       ZRHT,    & ! Graupel m.r. at t 
                                       !
                                       ZRVS,    & ! Water vapor m.r. source
                                       ZRCS,    & ! Cloud water m.r. source
                                       ZRRS,    & ! Rain water m.r. source
                                       ZRIS,    & ! Pristine ice m.r. source
                                       ZRSS,    & ! Snow/aggregate m.r. source
                                       ZRGS,    & ! Graupel/hail m.r. source
                                       ZRHS,    & ! Graupel/hail m.r. source
                                       !
                                       ZCCT,    & ! Cloud water C. at t
                                       ZCRT,    & ! Rain water C. at t
                                       ZCCS,    & ! Cloud water C. source
                                       ZCRS,    & ! Rain water C. source
                                       ZCIT,&
                                       ZCIS
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCIT_SHAPE, & ! Ice crystal C. per shape at t
                                         ZCIS_SHAPE    ! Ice crystal C. per shape source
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
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZNHS     ! Haze homogeneous activation
!
INTEGER :: ISH ! index for ice crystal shapes
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
ZCIS(:,:,:) = 0.
ZCIT(:,:,:) = 0.
!
IF ( LWARM ) ZCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) ZCRT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NR)
!
IF ( LWARM ) ZCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) ZCRS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NR)
!
!print*, ' LC1 psvs ni= ', minval(psvs(:,:,:,nsv_lima_ni)), maxval(psvs(:,:,:,nsv_lima_ni))
!print*, ' LC1 psvs ni+1= ', minval(psvs(:,:,:,nsv_lima_ni+1)), maxval(psvs(:,:,:,nsv_lima_ni+1))
IF (LCOLD) THEN
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    ZCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
    ZCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
!++cb++ 4/02/20 rajout car variable passee en argument non optionnel dans subroutine ! ne plante pas !!!???
    ALLOCATE(ZCIT_SHAPE(0,0,0,0))
    ALLOCATE(ZCIS_SHAPE(0,0,0,0))
!--cb--
  ELSE
    ALLOCATE(ZCIT_SHAPE(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NB_CRYSTAL_SHAPE))
    ALLOCATE(ZCIS_SHAPE(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3),NB_CRYSTAL_SHAPE))
!
    ZCIT_SHAPE(:,:,:,:) = PSVT(:,:,:,NSV_LIMA_NI:NSV_LIMA_NI+NB_CRYSTAL_SHAPE-1)
    ZCIS_SHAPE(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_NI:NSV_LIMA_NI+NB_CRYSTAL_SHAPE-1)
    ZCIT(:,:,:) = SUM(ZCIT_SHAPE, DIM=4)
    ZCIS(:,:,:) = SUM(ZCIS_SHAPE, DIM=4)
!print*, 'LC1 zcis_tot ni= ', minval(zcis_tot(:,:,:)), maxval(zcis_tot(:,:,:))
!print*, ' LC1 zcis_shape ni= ', minval(zcis_shape(:,:,:,1)), maxval(zcis_shape(:,:,:,1))
!print*, ' LC1 zcis_shape ni+1= ', minval(zcis_shape(:,:,:,2)), maxval(zcis_shape(:,:,:,2))
!    END DO
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
!*       1.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
CALL LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,      &
                              PZZ, PRHODJ, PRHODREF,            &
                              ZRIT, ZCIT,                   &
                              ZRIS, ZRSS, ZRGS, ZRHS, ZCIS, &
                              PINPRS, PINPRG,                   &
                              PINPRH, ZCIS_SHAPE                )
!print*, 'LCsedim zcis_tot ni= ', minval(zcis_tot(:,:,:)), maxval(zcis_tot(:,:,:))
!print*, ' LCsedim zcis_shape ni= ', minval(zcis_shape(:,:,:,1)), maxval(zcis_shape(:,:,:,1))
!print*, ' LCsedim zcis_shape ni+1= ', minval(zcis_shape(:,:,:,2)), maxval(zcis_shape(:,:,:,2))

IF (LBU_ENABLE) THEN
  IF (LBUDGET_RI .AND. OSEDI) CALL BUDGET (ZRIS(:,:,:)*PRHODJ(:,:,:),9  ,'SEDI_BU_RRI')
  IF (LBUDGET_RS .AND. LSNOW) CALL BUDGET (ZRSS(:,:,:)*PRHODJ(:,:,:),10 ,'SEDI_BU_RRS')
  IF (LBUDGET_RG .AND. LSNOW) CALL BUDGET (ZRGS(:,:,:)*PRHODJ(:,:,:),11 ,'SEDI_BU_RRG')
  IF (LBUDGET_RH .AND. LHAIL) CALL BUDGET (ZRHS(:,:,:)*PRHODJ(:,:,:),12 ,'SEDI_BU_RRH')
  IF (LBUDGET_SV) THEN
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      IF (OSEDI) &
        CALL BUDGET (ZCIS(:,:,:)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'SEDI_BU_RSV') ! RCI
    ELSE
      IF (OSEDI) THEN
        DO ISH = 1, NB_CRYSTAL_SHAPE
          CALL BUDGET (ZCIS_SHAPE(:,:,:,ISH)*PRHODJ(:,:,:), &
                                12+NSV_LIMA_NI+ISH-1,'SEDI_BU_RSV') ! RCI
        END DO
      END IF
    END IF
  END IF
END IF
!-------------------------------------------------------------------------------
!
!
!               COMPUTE THE NUCLEATION PROCESS SOURCES
!   	        --------------------------------------
!
IF (LNUCL) THEN
!
   IF ( LMEYERS ) THEN
      ZIFS(:,:,:,:) = 0.0
      ZNIS(:,:,:,:) = 0.0
      CALL LIMA_MEYERS (OHHONI, PTSTEP, KMI,                            &
                        PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,         &
                        PTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZCCT, &
                        PTHS, ZRVS, ZRCS, ZRIS,                         &
                        ZCCS, ZCIS, ZINS )
   ELSE
      CALL LIMA_PHILLIPS (OHHONI, PTSTEP, KMI,                      &
                          PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,   &
                          PTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, &
                          PTHS, ZRVS, ZRCS, ZRIS,                   &
                          ZCIT, ZCCS, ZCIS, ZCIS_SHAPE,             &
                          ZNAS, ZIFS, ZINS, ZNIS  )
   END IF
!print*, 'LCphil zcis_tot ni= ', minval(zcis_tot(:,:,:)), maxval(zcis_tot(:,:,:))
!print*, ' LCphil zcis_shape ni= ', minval(zcis_shape(:,:,:,1)), maxval(zcis_shape(:,:,:,1))
!print*, ' LCphil zcis_shape ni+1= ', minval(zcis_shape(:,:,:,2)), maxval(zcis_shape(:,:,:,2))
!
   IF (LWARM .OR. (LHHONI .AND. NMOD_CCN.GE.1)) THEN
      CALL LIMA_COLD_HOM_NUCL (OHHONI, PTSTEP, KMI,                         &
                               PZZ, PRHODJ,                                 &
                               PRHODREF, PEXNREF, PPABST, PW_NU,            &
                               PTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,    &
                               PTHS, ZRVS, ZRCS, ZRRS, ZRIS, ZRGS,          &
                               ZCCT,                                        &
                               ZCCS, ZCRS, ZNFS,                            &
                               ZCIS, ZCIS_SHAPE, ZNIS, ZNHS )
   END IF
!print*, 'LChom zcis_tot ni= ', minval(zcis_tot(:,:,:)), maxval(zcis_tot(:,:,:))
!print*, ' LChom zcis_shape ni= ', minval(zcis_shape(:,:,:,1)), maxval(zcis_shape(:,:,:,1))
!print*, ' LChom zcis_shape ni+1= ', minval(zcis_shape(:,:,:,2)), maxval(zcis_shape(:,:,:,2))
!
END IF
!
!------------------------------------------------------------------------------
!
!
!*       4.    SLOW PROCESSES: depositions, aggregation
!              ----------------------------------------
!
IF (LSNOW) THEN
!
   CALL LIMA_COLD_SLOW_PROCESSES(PTSTEP, KMI, PZZ, PRHODJ,                 &
                                 PRHODREF, PEXNREF, PPABST,                &
                                 PTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, &
                                 PTHS, ZRVS, ZRIS, ZRSS,                   &
                                 ZCIT, ZCIS, ZCIS_SHAPE, ZCIT_SHAPE)
!print*, 'LCslow zcis_tot ni= ', minval(zcis_tot(:,:,:)), maxval(zcis_tot(:,:,:))
!print*, ' LCslow zcis_shape ni= ', minval(zcis_shape(:,:,:,1)), maxval(zcis_shape(:,:,:,1))
!print*, ' LCslow zcis_shape ni+1= ', minval(zcis_shape(:,:,:,2)), maxval(zcis_shape(:,:,:,2))
!
END IF
!
!------------------------------------------------------------------------------
!
!
!*       4.    REPORT 3D MICROPHYSICAL VARIABLES IN PRS AND PSVS
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
IF ( LWARM ) PSVS(:,:,:,NSV_LIMA_NC) = ZCCS(:,:,:)
IF ( LWARM .AND. LRAIN ) PSVS(:,:,:,NSV_LIMA_NR) = ZCRS(:,:,:)
IF (.NOT. LCRYSTAL_SHAPE) THEN 
  PSVS(:,:,:,NSV_LIMA_NI) = ZCIS(:,:,:)
ELSE
  DO ISH = 1, NB_CRYSTAL_SHAPE
    PSVS(:,:,:,NSV_LIMA_NI+ISH-1) = ZCIS_SHAPE(:,:,:,ISH)
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

IF ( OHHONI ) PSVS(:,:,:,NSV_LIMA_HOM_HAZE) = ZNHS(:,:,:)
!
IF (ALLOCATED(ZNFS)) DEALLOCATE(ZNFS)
IF (ALLOCATED(ZNAS)) DEALLOCATE(ZNAS)
IF (ALLOCATED(ZIFS)) DEALLOCATE(ZIFS)
IF (ALLOCATED(ZINS)) DEALLOCATE(ZINS)
IF (ALLOCATED(ZNIS)) DEALLOCATE(ZNIS)
IF (ALLOCATED(ZNHS)) DEALLOCATE(ZNHS)
IF (ALLOCATED(ZCIS_SHAPE)) DEALLOCATE(ZCIS_SHAPE)
IF (ALLOCATED(ZCIT_SHAPE)) DEALLOCATE(ZCIT_SHAPE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COLD
