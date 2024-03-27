!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###################################
       MODULE MODI_LIMA_COLD_SEDIMENTATION
!      ###################################
!
INTERFACE
      SUBROUTINE LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,  &
                                          PZZ, PRHODJ, PRHODREF,        &
                                          PRIT, PCIT,                   &
                                          PRIS, PRSS, PRGS, PRHS, PCIS, &
                                          PINPRS, PINPRG, PINPRH,       &
                                          PCIS_SHAPE )
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian (Budgets)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT       ! Cloud ice m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS       ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS       ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS       ! Graupel m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRHS       ! Hail m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT   ! Ice crystal C. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS   ! Ice crystal C. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCIS_SHAPE !  ice concentration per shape
!
      END SUBROUTINE LIMA_COLD_SEDIMENTATION
END INTERFACE
END MODULE MODI_LIMA_COLD_SEDIMENTATION
!
!
!     ###################################################################
      SUBROUTINE LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,  &
                                          PZZ, PRHODJ, PRHODREF,        &
                                          PRIT, PCIT,                   &
                                          PRIS, PRSS, PRGS, PRHS, PCIS, &
                                          PINPRS, PINPRG, PINPRH,       &
                                          PCIS_SHAPE                    )
!     ###################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sedimentation
!!    of primary ice, snow and graupel.
!!
!!    METHOD
!!    ------
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
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
!!      M. Claeys  * LACy *  mar. 2019   add ice crytal shapes
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_COLD,  ONLY : XLBEXI, XLBI, XDI,                 &
                                  XFSEDRI, XFSEDCI, XFSEDS, XEXSEDS, &
                                  XDI_SHAPE, XFSEDCI_SHAPE, XLBI_SHAPE, &
                                  XLBEXI_SHAPE, XFSEDRI_SHAPE, XBI_SHAPE
USE MODD_PARAM_LIMA_MIXED, ONLY : XFSEDG, XEXSEDG, XFSEDH, XEXSEDH
USE MODD_PARAM_LIMA,       ONLY : XCEXVT, XRTMIN, XCTMIN,           &
                                  LCRYSTAL_SHAPE, NB_CRYSTAL_SHAPE, &
                                  XALPHAI, XNUI
USE MODD_CST,              ONLY : XRHOLW
USE MODD_PARAMETERS,       ONLY : JPHEXT, JPVEXT
USE MODD_NSV
!
USE MODI_LIMA_FUNCTIONS,   ONLY : COUNTJV
USE MODI_GAMMA
USE MODI_COMPUTE_LBDA_SHAPE
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI      ! switch to activate the 
                                                      ! cloud ice sedimentation
INTEGER,                  INTENT(IN)    :: KSPLITG    ! Number of small time step 
                                                      ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step          
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ        ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian (Budgets)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT       ! Cloud ice m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS       ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS       ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS       ! Graupel m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRHS       ! Hail m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS   ! Ice crystal C. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT   ! Ice crystal C. at t
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCIS_SHAPE ! ice conc. per shape
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JL, JN                     ! Loop index
INTEGER :: JSH                            ! Loop index for ice crystal shapes
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
INTEGER :: ISEDIM                         ! Case number of sedimentation
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: GSEDIM      ! Test where to compute the SED processes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
                           :: ZW,       & ! Work array
                              ZWSEDR,   & ! Sedimentation of MMR
                              ZWSEDC      ! Sedimentation of number conc.
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZWSEDC_SHAPE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZWSEDR_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZRIS_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZZY_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZZW_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZZX_SHAPE
!
REAL, DIMENSION(:), ALLOCATABLE           &
                           :: ZRIS,       & ! Pristine ice m.r. source
                              ZCIS,       & ! Pristine ice conc. source
                              ZRSS,       & ! Snow/aggregate m.r. source
                              ZRGS,       & ! Graupel/hail m.r. source
                              ZRHS,       & ! Graupel/hail m.r. source
                              ZRIT,       & ! Pristine ice m.r. at t
                              ZCIT,       & ! Pristine ice conc. at t
                              ZRHODREF,   & ! RHO Dry REFerence
                              ZRHODJ,     & ! RHO times Jacobian
                              ZZW,        & ! Work array
                              ZZX,        & ! Work array
                              ZZY,        & ! Work array
                              ZLBDAI,     & ! Slope parameter of the ice crystal distr.
                              ZRTMIN
!
INTEGER, DIMENSION(SIZE(PRHODREF)) :: I1,I2,I3 ! Indexes for PACK replacement
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCSHAPE_S    ! Ratio of ice conc. per shape
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZCSHAPE_S_V
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZCIS_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZLBDAI_SHAPE
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZONEOVER_VAR  ! for optimization
INTEGER, DIMENSION (:),   ALLOCATABLE :: ISHAPE_MAX   ! index of the dominant shape conc.
!
REAL    :: ZTSPLITG                       ! Small time step for rain sedimentation
!
!-------------------------------------------------------------------------------
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
! Time splitting and ZRTMIN
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
!
ZTSPLITG= PTSTEP / FLOAT(KSPLITG)
!
PINPRS(:,:) = 0.
PINPRG(:,:) = 0.
PINPRH(:,:) = 0.
!
! ################################
! Compute the sedimentation fluxes
! ################################
!
DO JN = 1 , KSPLITG 
  ! Computation only where enough ice, snow, graupel or hail
  GSEDIM(:,:,:) = .FALSE.
  GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRSS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(5) &
                               .OR. PRGS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(6) &
                               .OR. PRHS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(7)
  IF( OSEDI ) THEN
     GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) &
                               .OR. PRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4)
  END IF
!
  ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  IF( ISEDIM >= 1 ) THEN
!
    IF (LCRYSTAL_SHAPE) THEN
      ALLOCATE(ZCSHAPE_S(SIZE(PCIS,1),SIZE(PCIS,2),SIZE(PCIS,3),NB_CRYSTAL_SHAPE))
      ALLOCATE(ZONEOVER_VAR(SIZE(PCIS,1),SIZE(PCIS,2),SIZE(PCIS,3)))
      ZONEOVER_VAR(:,:,:) = 0.
    END IF
!
    IF( JN==1 ) THEN
      IF( OSEDI ) THEN
        IF (.NOT. LCRYSTAL_SHAPE) THEN
          PCIS(:,:,:) = PCIS(:,:,:) * PTSTEP
        ELSE
          DO JSH = 1, NB_CRYSTAL_SHAPE
            PCIS_SHAPE(:,:,:,JSH) = PCIS_SHAPE(:,:,:,JSH) * PTSTEP
          END DO
          PCIS = SUM(PCIS_SHAPE,DIM=4)
          ZCSHAPE_S(:,:,:,:) = 0.0
          WHERE (PCIS(:,:,:) .GT. 0.0) ZONEOVER_VAR(:,:,:) = 1.0 / PCIS(:,:,:)
          DO JSH = 1, NB_CRYSTAL_SHAPE
! compute the ratio : shape/tot
            WHERE ((PCIS(:,:,:) .GT. 0.0) .AND. (PCIS_SHAPE(:,:,:,JSH) .GT. 0.0))
              ZCSHAPE_S(:,:,:,JSH) = MIN(PCIS_SHAPE(:,:,:,JSH)*ZONEOVER_VAR(:,:,:), 1.0)
            END WHERE
          END DO
        END IF
        PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
      END IF
!
      PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
      PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
      PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
      DO JK = IKB , IKE
        ZW(:,:,JK)=ZTSPLITG/(PZZ(:,:,JK+1)-PZZ(:,:,JK))
      END DO
    END IF
!
    ALLOCATE(ZRHODREF(ISEDIM))
    DO JL = 1,ISEDIM
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    END DO
!
    ALLOCATE(ZZW(ISEDIM)) ; ZZW(:) = 0.0
    ALLOCATE(ZZX(ISEDIM)) ; ZZX(:) = 0.0
    ALLOCATE(ZZY(ISEDIM)) ; ZZY(:) = 0.0
    IF (LCRYSTAL_SHAPE) THEN
      ALLOCATE(ZZW_SHAPE(ISEDIM, NB_CRYSTAL_SHAPE)) ; ZZW_SHAPE(:,:) = 0.0
      ALLOCATE(ZZX_SHAPE(ISEDIM, NB_CRYSTAL_SHAPE)) ; ZZX_SHAPE(:,:) = 0.0
      ALLOCATE(ZZY_SHAPE(ISEDIM, NB_CRYSTAL_SHAPE)) ; ZZY_SHAPE(:,:) = 0.0
    END IF
!
!*       2.21   for pristine ice
!
    IF( OSEDI .AND. MAXVAL(PRIS(:,:,:)) > ZRTMIN(4) ) THEN
      ALLOCATE(ZRIS(ISEDIM))
      ALLOCATE(ZRIT(ISEDIM))
      ALLOCATE(ZCIS(ISEDIM))  
      ALLOCATE(ZCIT(ISEDIM))
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ALLOCATE(ZLBDAI(ISEDIM))
      ELSE
        ALLOCATE(ZLBDAI_SHAPE(ISEDIM, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCSHAPE_S_V(ISEDIM, NB_CRYSTAL_SHAPE))
        ALLOCATE(ZCIS_SHAPE(ISEDIM,NB_CRYSTAL_SHAPE))
        ALLOCATE(ZRIS_SHAPE(ISEDIM, NB_CRYSTAL_SHAPE))  
        ALLOCATE(ZWSEDC_SHAPE(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3),NB_CRYSTAL_SHAPE))
        ALLOCATE(ZWSEDR_SHAPE(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3),NB_CRYSTAL_SHAPE))
        ZRIS_SHAPE(:,:) = 0.0
        ZWSEDC_SHAPE(:,:,:,:) = 0.0
        ZWSEDR_SHAPE(:,:,:,:) = 0.0
      END IF
!
      DO JL = 1, ISEDIM
        ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
        ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
        ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
        ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
        IF (LCRYSTAL_SHAPE) THEN
          DO JSH = 1, NB_CRYSTAL_SHAPE
            ZCSHAPE_S_V(JL,JSH)  = ZCSHAPE_S(I1(JL),I2(JL),I3(JL),JSH)
            ZCIS_SHAPE(JL,JSH) = PCIS_SHAPE(I1(JL),I2(JL),I3(JL),JSH)
          END DO
        END IF
      END DO
!
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        ZLBDAI(:) = 1.E10
        WHERE (ZRIT(:) > XRTMIN(4) .AND. ZCIT(:) > XCTMIN(4))
          ZLBDAI(:) = ( XLBI*ZCIT(:) / ZRIT(:) )**XLBEXI
          ZZY(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDAI(:)**(-XDI)
          ZZW(:) = XFSEDRI * ZRIS(:)       * ZZY(:) * ZRHODREF(:)
          ZZX(:) = XFSEDCI * ZCIS(:) * ZZY(:) * ZRHODREF(:)
        END WHERE
        ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
        ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRIS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
        ZWSEDC(:,:,:) = UNPACK( ZZX(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
        ZWSEDC(:,:,IKB:IKE) = MIN( ZWSEDC(:,:,IKB:IKE), PCIS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
        DO JK = IKB, IKE
          PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK) *    &
               (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
          PCIS(:,:,JK) = PCIS(:,:,JK) + ZW(:,:,JK) *    &
               (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))/PRHODREF(:,:,JK)
        END DO
!
        DEALLOCATE(ZLBDAI)
!
      ELSE
!
        ALLOCATE(ISHAPE_MAX(ISEDIM))
        ISHAPE_MAX(:) = MAXLOC(ZCSHAPE_S_V,DIM=2)
        ZLBDAI_SHAPE(:,:) = 1E10
        CALL COMPUTE_LBDA_SHAPE(ISHAPE_MAX, ISEDIM, PTSTEP,     &
                                ZRIT, ZCIT, ZLBDAI_SHAPE, &
                                ZRIS, ZCIS_SHAPE, ZRIS_SHAPE)
!
        DO JSH = 1, NB_CRYSTAL_SHAPE
          WHERE (ZRIS(:) > ZRTMIN(4))
            ZZY_SHAPE(:,JSH) = ZRHODREF(:)**(-XCEXVT) * &
                               ZLBDAI_SHAPE(:,JSH)**(-XDI_SHAPE(JSH))
            ZZW_SHAPE(:,JSH) = XFSEDRI_SHAPE(JSH) * ZRIS_SHAPE(:,JSH) * &
                               ZZY_SHAPE(:,JSH) * ZRHODREF(:)
            ZZX_SHAPE(:,JSH) = XFSEDCI_SHAPE(JSH) * ZCIS_SHAPE(:,JSH) * &
                               ZZY_SHAPE(:,JSH) * ZRHODREF(:)
          END WHERE
        END DO
!
        ZONEOVER_VAR(:,:,:) = 1.0 / PRHODREF(:,:,:)
        DO JSH = 1, NB_CRYSTAL_SHAPE
          ZWSEDR_SHAPE(:,:,:,JSH) = UNPACK( ZZW_SHAPE(:,JSH),MASK=GSEDIM(:,:,:),FIELD=0.0 )
          ZWSEDC_SHAPE(:,:,:,JSH) = UNPACK( ZZX_SHAPE(:,JSH),MASK=GSEDIM(:,:,:),FIELD=0.0 )
          ZWSEDC_SHAPE(:,:,IKB:IKE,JSH) = MIN( ZWSEDC_SHAPE(:,:,IKB:IKE,JSH),&
                            PCIS_SHAPE(:,:,IKB:IKE,JSH) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
        END DO
        ZWSEDR(:,:,:) = SUM(ZWSEDR_SHAPE, DIM=4)
        ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRIS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
        DO JK = IKB, IKE
          PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK) *    &
                        (ZWSEDR(:,:,JK+1) - ZWSEDR(:,:,JK)) * ZONEOVER_VAR(:,:,JK)
          DO JSH = 1, NB_CRYSTAL_SHAPE
            PCIS_SHAPE(:,:,JK,JSH) = PCIS_SHAPE(:,:,JK,JSH) + ZW(:,:,JK) *    &
                                    (ZWSEDC_SHAPE(:,:,JK+1,JSH) - ZWSEDC_SHAPE(:,:,JK,JSH)) * &
                                     ZONEOVER_VAR(:,:,JK)
          END DO
        END DO
        PCIS = SUM(PCIS_SHAPE,DIM=4)
!
        DEALLOCATE(ZLBDAI_SHAPE)
        DEALLOCATE(ZCSHAPE_S_V)
        DEALLOCATE(ISHAPE_MAX)
        DEALLOCATE(ZRIS_SHAPE)
        DEALLOCATE(ZCIS_SHAPE)
        DEALLOCATE(ZWSEDC_SHAPE)
        DEALLOCATE(ZWSEDR_SHAPE)
      END IF 
      DEALLOCATE(ZRIS)
      DEALLOCATE(ZRIT)
      DEALLOCATE(ZCIT)
      DEALLOCATE(ZCIS)
    END IF
    IF (ALLOCATED(ZCSHAPE_S)) DEALLOCATE(ZCSHAPE_S)
    IF (ALLOCATED(ZONEOVER_VAR)) DEALLOCATE(ZONEOVER_VAR)
!
    IF (ALLOCATED(ZZY_SHAPE)) DEALLOCATE(ZZY_SHAPE)
    IF (ALLOCATED(ZZW_SHAPE)) DEALLOCATE(ZZW_SHAPE)
    IF (ALLOCATED(ZZX_SHAPE)) DEALLOCATE(ZZX_SHAPE)
!
!
!*       2.22   for aggregates
!
    ZZW(:) = 0.
    IF( MAXVAL(PRSS(:,:,:))>XRTMIN(5) ) THEN
      ALLOCATE(ZRSS(ISEDIM)) 
      DO JL = 1,ISEDIM
        ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
      END DO
      WHERE( ZRSS(:)>XRTMIN(5) )
! Correction BVIE ZRHODREF
!        ZZW(:) = XFSEDS * ZRSS(:)**XEXSEDS * ZRHODREF(:)**(XEXSEDS-XCEXVT)
        ZZW(:) = XFSEDS * ZRSS(:)**XEXSEDS * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
      END WHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRSS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
      DO JK = IKB , IKE
        PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)* &
             (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
      END DO
      DEALLOCATE(ZRSS)
    ELSE
      ZWSEDR(:,:,IKB) = 0.0
    END IF
!    
    PINPRS(:,:) = PINPRS(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                          ! in m/s
!
!*       2.23   for graupeln
!
    ZZW(:) = 0.
    IF( MAXVAL(PRGS(:,:,:))>XRTMIN(6) ) THEN
      ALLOCATE(ZRGS(ISEDIM)) 
      DO JL = 1,ISEDIM
        ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
      END DO
      WHERE( ZRGS(:)>XRTMIN(6) )
! Correction BVIE ZRHODREF
!        ZZW(:) = XFSEDG * ZRGS(:)**XEXSEDG * ZRHODREF(:)**(XEXSEDG-XCEXVT)
        ZZW(:) = XFSEDG * ZRGS(:)**XEXSEDG * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
      END WHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRGS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
      DO JK = IKB , IKE
        PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)* &
             (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
      END DO
      DEALLOCATE(ZRGS)
    ELSE
      ZWSEDR(:,:,IKB) = 0.0
    END IF
!    
    PINPRG(:,:) = PINPRG(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                        ! in m/s
!
!*       2.23   for hail
!
    ZZW(:) = 0.
    IF( MAXVAL(PRHS(:,:,:))>XRTMIN(7) ) THEN
      ALLOCATE(ZRHS(ISEDIM)) 
      DO JL = 1,ISEDIM
        ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
      END DO
      WHERE( ZRHS(:)>XRTMIN(7) )
! Correction BVIE ZRHODREF
!        ZZW(:) = XFSEDH * ZRHS(:)**XEXSEDH * ZRHODREF(:)**(XEXSEDH-XCEXVT)
        ZZW(:) = XFSEDH * ZRHS(:)**XEXSEDH * ZRHODREF(:)**(-XCEXVT) * ZRHODREF(:)
      END WHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDR(:,:,IKB:IKE) = MIN( ZWSEDR(:,:,IKB:IKE), PRHS(:,:,IKB:IKE) * PRHODREF(:,:,IKB:IKE) / ZW(:,:,IKB:IKE) )
      DO JK = IKB , IKE
        PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)* &
             (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
      END DO
      DEALLOCATE(ZRHS)
    ELSE
      ZWSEDR(:,:,IKB) = 0.0
    END IF
!    
    PINPRH(:,:) = PINPRH(:,:) + ZWSEDR(:,:,IKB)/XRHOLW/KSPLITG                        ! in m/s
!
!*       2.24 End of sedimentation  
!
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZW)
    DEALLOCATE(ZZX)
    DEALLOCATE(ZZY)
!
    IF( JN==KSPLITG ) THEN
      IF( OSEDI ) THEN
        PRIS(:,:,:) = PRIS(:,:,:) / PTSTEP
        IF (.NOT. LCRYSTAL_SHAPE) THEN
          PCIS(:,:,:) = PCIS(:,:,:) / PTSTEP
        ELSE
          DO JSH = 1, NB_CRYSTAL_SHAPE
            PCIS_SHAPE(:,:,:,JSH) = PCIS_SHAPE(:,:,:,JSH) / PTSTEP
          END DO
          PCIS(:,:,:) = SUM(PCIS_SHAPE, DIM=4) 
        END IF
      END IF
      PRSS(:,:,:) = PRSS(:,:,:) / PTSTEP
      PRGS(:,:,:) = PRGS(:,:,:) / PTSTEP
      PRHS(:,:,:) = PRHS(:,:,:) / PTSTEP
    END IF
  END IF
END DO
!
DEALLOCATE(ZRTMIN)
!
!
END SUBROUTINE LIMA_COLD_SEDIMENTATION
!
!-------------------------------------------------------------------------------
