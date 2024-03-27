!     ##############################
      MODULE MODI_COMPUTE_LBDA_SHAPE
!     ##############################
!
INTERFACE
!
      SUBROUTINE COMPUTE_LBDA_SHAPE(KSHAPE_MAX, KLOOP, PTSTEP, &
                                    PRIT, PCIT,            &
                                    PLBDAI_SHAPE,              &
                                    PRI,PCI_SHAPE,PRI_SHAPE)
!
INTEGER, DIMENSION(:),   INTENT(IN)    :: KSHAPE_MAX
INTEGER,                 INTENT(IN)    :: KLOOP
REAL,                    INTENT(IN)    :: PTSTEP   ! Time step
REAL,    DIMENSION(:),   INTENT(IN)    :: PRIT     ! Cloud ice m.r. at t
!REAL,    DIMENSION(:),   INTENT(IN)    :: PRIS     ! Pristine ice m.r. source
REAL,    DIMENSION(:),   INTENT(IN)    :: PCIT ! total pristine ice conc. at t
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
REAL,    DIMENSION(:),   INTENT(IN),    OPTIONAL :: PRI       ! total pristine ice m.r. (source or at t)
REAL,    DIMENSION(:,:), INTENT(IN),    OPTIONAL :: PCI_SHAPE ! Pristine ice conc. (source or at t)
REAL,    DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: PRI_SHAPE ! pristine ice m.r. (source or at t)
!
      END SUBROUTINE COMPUTE_LBDA_SHAPE
END INTERFACE
END MODULE MODI_COMPUTE_LBDA_SHAPE
!
!     ###################################################################
      SUBROUTINE COMPUTE_LBDA_SHAPE (KSHAPE_MAX, KLOOP, PTSTEP, &
                                     PRIT, PCIT,            &
                                     PLBDAI_SHAPE,              &
                                     PRI,PCI_SHAPE,PRI_SHAPE)
!     ###################################################################
!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!      C. Barthe  * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/04/2018
!!      C. Barthe      26/11/2019 : supprime PRIS --> inutile
!!
!-----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_PARAM_LIMA, ONLY : XRTMIN, XCTMIN, NB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : XLBI_SHAPE, XLBEXI_SHAPE, XBI_SHAPE, &
                                 XFSEDRI_TOT_SHAPE
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
INTEGER, DIMENSION(:),   INTENT(IN)    :: KSHAPE_MAX
INTEGER,                 INTENT(IN)    :: KLOOP
REAL,                    INTENT(IN)    :: PTSTEP  ! Time step
REAL,    DIMENSION(:),   INTENT(IN)    :: PRIT     ! Cloud ice m.r. at t
!REAL,    DIMENSION(:),   INTENT(IN)    :: PRIS     ! Pristine ice m.r. source
REAL,    DIMENSION(:),   INTENT(IN)    :: PCIT ! total pristine ice conc. at t
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
REAL,    DIMENSION(:),   INTENT(IN),    OPTIONAL :: PRI       ! total pristine ice m.r. (source or at t)
REAL,    DIMENSION(:,:), INTENT(IN),    OPTIONAL :: PCI_SHAPE ! Pristine ice conc. (source or at t)
REAL,    DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: PRI_SHAPE ! pristine ice m.r. (source or at t)
!
!
!*      0.2   Declaration of local variables
!
INTEGER :: II, JSH, &    ! loop counter
           ICOUNT_SHAPE  ! count the number of points where a shape dominates
REAL    :: ZRTMIN
REAL, DIMENSION(:),   ALLOCATABLE :: ZLBDAI_AVG ! mean value of lambda per habit
REAL, DIMENSION(:),   ALLOCATABLE :: ZLBDAI
!REAL, DIMENSION(:),   ALLOCATABLE :: ZRTMIN, ZCTMIN
REAL, DIMENSION(:),   ALLOCATABLE :: ZAUX, ZSUM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRSHAPE  ! mass ratio for each habit
!
!-----------------------------------------------------------------------------
!
!ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
!ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZRTMIN = XRTMIN(4) / PTSTEP
!ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
!ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
ALLOCATE(ZLBDAI(KLOOP))
ZLBDAI(:) = 1.E10

!
!
!*      1.    COMPUTE LAMBDA FOR THE MAJOR HABIT PER MESH
!             -------------------------------------------
!
! calcul de ZLBDAI attache a l'espece dominante par maille 
WHERE (PRIT(:) > XRTMIN(4) .AND. PCIT(:) > XCTMIN(4))
  ZLBDAI(:) = (XLBI_SHAPE(KSHAPE_MAX(:)) * PCIT(:) / &
              PRIT(:))**XLBEXI_SHAPE(KSHAPE_MAX(:))
END WHERE
!
!
!*      2.    COMPUTE A MEAN LAMBDA PER HABIT
!             -------------------------------
!
ALLOCATE(ZLBDAI_AVG(NB_CRYSTAL_SHAPE))
ZLBDAI_AVG(:) = 0.0
!
! Note : to be used for non dominant habits in the mesh
DO JSH = 1, NB_CRYSTAL_SHAPE
  ICOUNT_SHAPE = COUNT(KSHAPE_MAX(:) .EQ. JSH .AND. ZLBDAI(:) .NE. 1.E10)
  IF (ICOUNT_SHAPE .GT. 0) THEN
    DO II = 1, KLOOP
      IF (KSHAPE_MAX(II) .EQ. JSH .AND. ZLBDAI(II) .NE. 1.E10) THEN
        ZLBDAI_AVG(JSH) = ZLBDAI_AVG(JSH) + ZLBDAI(II)
      END IF
    END DO
    ZLBDAI_AVG(JSH) = ZLBDAI_AVG(JSH) / ICOUNT_SHAPE
!
    WHERE (KSHAPE_MAX(:) .EQ. JSH)
      PLBDAI_SHAPE(:,JSH) = ZLBDAI(:)
    ELSEWHERE (KSHAPE_MAX(:) .NE. JSH .AND. ZLBDAI(:) .LT. 1.E10)
      PLBDAI_SHAPE(:,JSH) = ZLBDAI_AVG(JSH)
    END WHERE
  END IF
!
!
!*      3.    COMPUTE THE MIXING RATIO PER HABIT
!             ----------------------------------
!
  IF (PRESENT(PRI_SHAPE)) THEN
! marine : changement des conditions de calcul (comme lima_cold_sedim)
!    WHERE (PRI(:) > ZRTMIN(4) .AND. PCI_SHAPE(:,JSH) > ZCTMIN(4) .AND. &
!           JSH .EQ. KSHAPE_MAX(:))
!++cb++ 31/03/20
!    WHERE (PRI(:) > ZRTMIN .AND. JSH .EQ. KSHAPE_MAX(:))
    WHERE (PRI(:) > ZRTMIN)
      PRI_SHAPE(:,JSH) = XFSEDRI_TOT_SHAPE(JSH) * &
                         PLBDAI_SHAPE(:,JSH)**(-XBI_SHAPE(JSH)) * &
                         PCI_SHAPE(:,JSH)
!    ELSE WHERE (PRI(:) > ZRTMIN(4) .AND. PCI_SHAPE(:,JSH) > ZCTMIN(4) .AND. &
!               JSH .NE. KSHAPE_MAX(:))
!    ELSE WHERE (PRI(:) > ZRTMIN .AND. JSH .NE. KSHAPE_MAX(:))
!      PRI_SHAPE(:,JSH) = XFSEDRI_TOT_SHAPE(JSH) * &
!                         PLBDAI_SHAPE(:,JSH)**(-XBI_SHAPE(JSH)) * &
!                         PCI_SHAPE(:,JSH)
    END WHERE
  END IF
END DO
!
DEALLOCATE(ZLBDAI_AVG)
DEALLOCATE(ZLBDAI)
!
IF (PRESENT(PRI_SHAPE)) THEN
  ALLOCATE(ZAUX(KLOOP))
  ALLOCATE(ZSUM(KLOOP))
  ALLOCATE(ZRSHAPE(KLOOP,NB_CRYSTAL_SHAPE))
  ZRSHAPE(:,:) = 0.
!
! on s'assure que la somme des masses = la masse totale, sinon, on redistribue en fonction
! du pourcentage en masse de chaque forme
  ZSUM(:) = SUM(PRI_SHAPE,DIM=2)
  ZAUX(:) = PRI(:) - ZSUM(:) ! >0 si pas assez de masse / < 0 si masse en trop
! rapport en masse de chaque forme
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE ((PRI(:) .GT. 0.0) .AND. (PRI_SHAPE(:,JSH) .GT. 0.0)) ! ne faut-il pas plutot mettre .ne. 0
      ZRSHAPE(:,JSH) = MIN(PRI_SHAPE(:,JSH) / ZSUM(:), 1.0)
    END WHERE
!  END DO
!           
!  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE (ZAUX(:) .GT. 0.1 .OR. ZAUX(:) .LT. 0.1) ! erreur < 1%
      PRI_SHAPE(:,JSH) = PRI_SHAPE(:,JSH) + ZRSHAPE(:,JSH) * ZAUX(:)  
    END WHERE
    WHERE (PRI(:) .EQ. 0.)
      PRI_SHAPE(:,JSH) = 0.
    END WHERE
  END DO
!
  DEALLOCATE(ZAUX)
  DEALLOCATE(ZSUM)
  DEALLOCATE(ZRSHAPE)
END IF
!
END SUBROUTINE COMPUTE_LBDA_SHAPE
