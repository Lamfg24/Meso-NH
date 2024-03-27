!      #####################################
       MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!      #####################################
!
INTERFACE
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(PRHODREF, PZT, PSSI, PTSTEP,        &
                                           PLSFACT, PLVFACT, PAI, PCJ,         &
                                           PRGT, PCIT, PCIT_SHAPE,             &
                                           PRVS, PRCS, PRIS, PRGS, PTHS,       &
                                           PCCS, PCIS, PCIS_SHAPE,             &
                                           PIFS, PINS,                         &
                                           PLBDAI, PLBDAG,                     &
                                           OMICRO, PRHODJ, KMI, PTHS_3D,       &
                                           PRVS_3D, PRCS_3D, PRIS_3D, PRGS_3D, &
                                           PCCS_3D, PCIS_3D, PCIS_SHAPE_3D,&
                                           PLBDAI_SHAPE)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: PZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time-step
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT      ! Pristine ice conc. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PCIT_SHAPE      ! Pristine ice conc.
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS      ! Theta source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS      ! Pristine ice conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCIS_SHAPE      ! Pristine ice conc. at t
REAL, DIMENSION(:,:), INTENT(INOUT) :: PIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: PINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope parameter of the graupel distr.
!
! used for budget storage
LOGICAL, DIMENSION(:,:,:),   INTENT(IN) :: OMICRO 
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ
INTEGER,                     INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PTHS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRVS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRIS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRGS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCIS_3D
REAL,    DIMENSION(:,:,:,:), INTENT(IN) :: PCIS_SHAPE_3D      ! Pristine ice conc.
REAL,    DIMENSION(:,:),     INTENT(INOUT), OPTIONAL :: PLBDAI_SHAPE
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!
!     ##########################################################################
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(PRHODREF, PZT, PSSI, PTSTEP,        &
                                           PLSFACT, PLVFACT, PAI, PCJ,         &
                                           PRGT, PCIT, PCIT_SHAPE,             &
                                           PRVS, PRCS, PRIS, PRGS, PTHS,       &
                                           PCCS, PCIS, PCIS_SHAPE,             &
                                           PIFS, PINS,                         &
                                           PLBDAI, PLBDAG,                     &
                                           OMICRO, PRHODJ, KMI, PTHS_3D,       &
                                           PRVS_3D, PRCS_3D, PRIS_3D, PRGS_3D, &
                                           PCCS_3D, PCIS_3D, PCIS_SHAPE_3D,&
                                           PLBDAI_SHAPE)
!     ##########################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    slow processes : 
!!
!!      Deposition of water vapor on graupeln
!!      Cloud ice Melting
!!      Bergeron-Findeisen effect
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
!!      M. Claeys  * LACy *  mar. 2019   add ice crystal shapes
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT, XALPI, XBETAI, XGAMI,          &
                                       XALPW, XBETAW, XGAMW
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, NMOD_IFN, LSNOW,    &
                                  LCRYSTAL_SHAPE, NB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD,  ONLY : XDI, X0DEPI, X2DEPI, XSCFAC,        &
                                  XDI_SHAPE, X0DEPI_SHAPE, X2DEPI_SHAPE
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XLBDAG_MAX,           &
                                  X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG 
!
USE MODD_NSV
USE MODD_BUDGET
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: PZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time-step
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT      ! Pristine ice conc. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PCIT_SHAPE      ! Pristine ice conc.
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS      ! Theta source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS      ! Pristine ice conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCIS_SHAPE      ! Pristine ice conc. at t
REAL, DIMENSION(:,:), INTENT(INOUT) :: PIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: PINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope parameter of the graupel distr.
!
! used for budget storage
LOGICAL, DIMENSION(:,:,:),   INTENT(IN) :: OMICRO
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ
INTEGER,                     INTENT(IN) :: KMI
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PTHS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRVS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRIS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRGS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCIS_3D
REAL,    DIMENSION(:,:,:,:), INTENT(IN) :: PCIS_SHAPE_3D      ! Pristine ice conc.
REAL,    DIMENSION(:,:),     INTENT(INOUT), OPTIONAL :: PLBDAI_SHAPE
!
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PZT)) :: ZZW, ZMASK    ! Work vectors
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW_2D      ! 
!
INTEGER :: JMOD_IFN
INTEGER :: JSH
!
!-------------------------------------------------------------------------------
!
!*       1    Deposition of water vapor on r_g: RVDEPG
!        ---------------------------------------------
!
!
IF (LSNOW) THEN
  ZZW(:) = 0.0
  WHERE ( (PRGT(:)>XRTMIN(6)) .AND. (PRGS(:)>XRTMIN(6)/PTSTEP) )
!Correction BVIE RHODREF
!      ZZW(:) = ( PSSI(:)/(PRHODREF(:)*PAI(:)) ) *                               &
    ZZW(:) = ( PSSI(:)/(PAI(:)) ) *                               &
             ( X0DEPG*PLBDAG(:)**XEX0DEPG + X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG )
    ZZW(:) =         MIN( PRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                   - MIN( PRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
    PRGS(:) = PRGS(:) + ZZW(:)
    PRVS(:) = PRVS(:) - ZZW(:)
    PTHS(:) = PTHS(:) + ZZW(:)*PLSFACT(:)
  END WHERE
!
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET (                                             &
                  UNPACK(PTHS(:),MASK=OMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
                                                               4,'DEPG_BU_RTH')
    IF (LBUDGET_RV) CALL BUDGET (                                             &
                  UNPACK(PRVS(:),MASK=OMICRO(:,:,:),FIELD=PRVS_3D)*PRHODJ(:,:,:),&
                                                               6,'DEPG_BU_RRV')
    IF (LBUDGET_RG) CALL BUDGET (                                             &
                  UNPACK(PRGS(:),MASK=OMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:),&
                                                              11,'DEPG_BU_RRG')
  END IF
END IF
!
!
!*       2    cloud ice Melting: RIMLTC and CIMLTC
!        -----------------------------------------
!
!
ZMASK(:) = 1.0
IF (.NOT. LCRYSTAL_SHAPE) THEN
  WHERE( (PRIS(:)>XRTMIN(4)/PTSTEP) .AND. (PZT(:)>XTT) )
    PRCS(:) = PRCS(:) + PRIS(:)
    PTHS(:) = PTHS(:) - PRIS(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RIMLTC))
    PRIS(:) = 0.0
!
    PCCS(:) = PCCS(:) + PCIS(:)
    PCIS(:) = 0.0
    ZMASK(:)= 0.0
  END WHERE
ELSE
  WHERE( (PRIS(:)>XRTMIN(4)/PTSTEP) .AND. (PZT(:)>XTT) )
    PRCS(:) = PRCS(:) + PRIS(:)
    PTHS(:) = PTHS(:) - PRIS(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RIMLTC))
    PRIS(:) = 0.0
    PCCS(:) = PCCS(:) + SUM(PCIS_SHAPE,DIM=2)
    PCIS(:) = 0.0
    ZMASK(:)= 0.0
  END WHERE
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE (PCIS(:) .EQ. 0.)
      PCIS_SHAPE(:,JSH) = 0.0
    END WHERE
  END DO
END IF   ! LCRYSTAL_SHAPE

DO JMOD_IFN = 1,NMOD_IFN
! Correction BVIE aerosols not released but in droplets
!      PIFS(:,JMOD_IFN) = PIFS(:,JMOD_IFN) + PINS(:,JMOD_IFN)*(1.-ZMASK(:)) 
  PINS(:,JMOD_IFN) = PINS(:,JMOD_IFN) * ZMASK(:)
ENDDO
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                              &
                UNPACK(PTHS(:),MASK=OMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:), &
                                                             4,'IMLT_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET (                                              &
                UNPACK(PRCS(:),MASK=OMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
                                                             7,'IMLT_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET (                                              &
                UNPACK(PRIS(:),MASK=OMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                             9,'IMLT_BU_RRI')
  IF (LBUDGET_SV) THEN
    CALL BUDGET (UNPACK(PCCS(:),MASK=OMICRO(:,:,:),FIELD=PCCS_3D)*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NC,'IMLT_BU_RSV')
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET (UNPACK(PCIS(:),MASK=OMICRO(:,:,:),FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NI,'IMLT_BU_RSV')
    ELSE 
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET (UNPACK(PCIS_SHAPE(:,JSH),MASK=OMICRO(:,:,:),                         &
                                                FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                12+NSV_LIMA_NI+JSH-1,'IMLT_BU_RSV')
      END DO
    END IF  ! LCRYSTAL_SHAPE
  END IF
END IF
!
!
!*       3    Bergeron-Findeisen effect: RCBERI
!        --------------------------------------
!
!
IF (.NOT. LCRYSTAL_SHAPE) THEN
  ZZW(:) = 0.0
  WHERE( (PRCS(:)>XRTMIN(2)/PTSTEP) .AND. (PRIS(:)>XRTMIN(4)/PTSTEP) &
                                    .AND. (PCIT(:) > XCTMIN(4)) )
    ZZW(:) = EXP( (XALPW-XALPI) - (XBETAW-XBETAI)/PZT(:)          &
                                - (XGAMW-XGAMI)*ALOG(PZT(:)) ) -1.0 
                                  ! supersaturation of saturated water over ice
    ZZW(:) = MIN( PRCS(:),( ZZW(:) / PAI(:) ) * PCIT(:) *        &
                  ( X0DEPI/PLBDAI(:)+X2DEPI*PCJ(:)*PCJ(:)/PLBDAI(:)**(XDI+2.0) ) )
    PRCS(:) = PRCS(:) - ZZW(:)
    PRIS(:) = PRIS(:) + ZZW(:)
    PTHS(:) = PTHS(:) + ZZW(:) * (PLSFACT(:) - PLVFACT(:)) ! f(L_f*(RCBERI))
  END WHERE
ELSE
  ALLOCATE(ZZW_2D(SIZE(PZT),NB_CRYSTAL_SHAPE))
  ZZW(:) = 0.0
  ZZW_2D(:,:) = 0.0
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE((PRCS(:) > XRTMIN(2)/PTSTEP) .AND. (PRIS(:) > XRTMIN(4)/PTSTEP) &
                                       .AND. (PCIT_SHAPE(:,JSH) > XCTMIN(4)))
      ZZW(:) = EXP( (XALPW-XALPI) - (XBETAW - XBETAI) / PZT(:)          &
                                  - (XGAMW - XGAMI) * ALOG(PZT(:)) ) -1.0 
                            ! supersaturation of saturated water over ice
      ZZW_2D(:,JSH) = MIN( PRCS(:),                               &
                          (ZZW(:)/PAI(:))*PCIT_SHAPE(:,JSH) *   &
                          (X0DEPI_SHAPE(JSH)/PLBDAI_SHAPE(:,JSH)+ &
                          X2DEPI_SHAPE(JSH)*PCJ(:)*PCJ(:)/        &
                          PLBDAI_SHAPE(:,JSH)**(XDI_SHAPE(JSH)+2.0)) )
    END WHERE
  END DO ! JSH
  ZZW(:)  = 0.
  ZZW(:)  = MIN(PRCS, SUM(ZZW_2D, DIM=2))
  PRCS(:) = PRCS(:) - ZZW(:)
  PRIS(:) = PRIS(:) + ZZW(:)
  PTHS(:) = PTHS(:) + ZZW(:) * (PLSFACT(:) - PLVFACT(:)) ! f(L_f*(RCBERI))
END IF  ! LCRYSTAL_SHAPE
!
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                              &
                UNPACK(PTHS(:),MASK=OMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:), &
                                                            4,'BERFI_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET (                                              &
                UNPACK(PRCS(:),MASK=OMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
                                                            7,'BERFI_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET (                                              &
                UNPACK(PRIS(:),MASK=OMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                            9,'BERFI_BU_RRI')
END IF
!
IF (ALLOCATED(ZZW_2D)) DEALLOCATE(ZZW_2D)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
