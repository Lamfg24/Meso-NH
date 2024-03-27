!     ##############################
      MODULE MODI_LIMA_CHANGE_SHAPE
!     ##############################
!
INTERFACE
!
      SUBROUTINE LIMA_CHANGE_SHAPE(PZTM, PZTT,                   &
                                   HSHAPE_NAME_M, HSHAPE_NAME_T, &
                                   PCIS_SHAPE_V, PCTMIN, PDEP_V)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PZTM  ! Temp at t-1 (K)
REAL, DIMENSION(:),   INTENT(IN)    :: PZTT  ! Temp at t (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCIS_SHAPE_V  ! Nb conc. of ice crystal shapes
REAL, DIMENSION(:)                  :: PCTMIN
REAL, DIMENSION(:),   INTENT(IN)    :: PDEP_V  ! water vapor deposited on ice crystals
CHARACTER(LEN=3), DIMENSION(:), INTENT(IN) :: HSHAPE_NAME_M ! Name of the shape corresponding to 
                                                            ! the temperature growth regime, at t-1
CHARACTER(LEN=3), DIMENSION(:), INTENT(IN) :: HSHAPE_NAME_T ! Name of the shape corresponding to 
                                                            ! the temperature growth regime, at t
!
      END SUBROUTINE LIMA_CHANGE_SHAPE
END INTERFACE
END MODULE MODI_LIMA_CHANGE_SHAPE
!
!
!     #############################################################
      SUBROUTINE LIMA_CHANGE_SHAPE (PZTM, PZTT,                   &
                                    HSHAPE_NAME_M, HSHAPE_NAME_T, &
                                    PCIS_SHAPE_V, PCTMIN, PDEP_V)
!     #############################################################
!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the change of shape of small ice
!!    crystals, when plates are in a temperature range of column growth and vice
!!    versa. Plates and columns turn into irregulars definitely.
!! 
!!    QUESTIONS : -  Rajouter une valeur seuil sur la quantité de vapeur d'eau
!!                déposée sur les cristaux (PDEP en entrée de la routine pour
!!                l'instant mais qui ne sert pas pour l'instant)?
!!                - 
!!
!!    METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!      M. Claeys  * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/11/2018
!!
!-----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_PARAM_LIMA, ONLY : NB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : XLBI_SHAPE, XLBEXI_SHAPE, XBI_SHAPE!, &
                                 !XFSEDRI_TOT_SHAPE
USE MODD_CST, ONLY : XTT
USE MODD_BUDGET
USE MODI_BUDGET
!
IMPLICIT NONE
!
!
!*      0.1   Declaration of dummy arguments
!
REAL, DIMENSION(:),   INTENT(IN)    :: PZTM  ! Temp at t-1 (K)
REAL, DIMENSION(:),   INTENT(IN)    :: PZTT  ! Temp at t (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCIS_SHAPE_V  ! Nb conc. of ice crystal shapes
REAL, DIMENSION(:)                  :: PCTMIN
REAL, DIMENSION(:),   INTENT(IN)    :: PDEP_V  ! water vapor deposited on ice crystals
CHARACTER(LEN=3), DIMENSION(:), INTENT(IN) :: HSHAPE_NAME_M ! Name of the shape corresponding to 
                                                            ! the temperature growth regime, at t-1
CHARACTER(LEN=3), DIMENSION(:), INTENT(IN) :: HSHAPE_NAME_T ! Name of the shape corresponding to 
                                                            ! the temperature growth regime, at t
!
!
!*      0.2   Declaration of local variables
!
REAL, DIMENSION(SIZE(PZTT,1))         :: ZTT_CELSIUS ! temperature (C)
!REAL, DIMENSION(SIZE(PZTM,1))         :: ZTM_CELSIUS ! temperature (C)
REAL, DIMENSION(SIZE(PZTM,1))         :: ZFACTOR_PLA_IRR !  
REAL, DIMENSION(SIZE(PZTM,1))         :: ZFACTOR_COL_IRR !
!REAL, DIMENSION(SIZE(PCIS_SHAPE_V,1)) :: ZDN_DGAMMA_PLA_IRR
!REAL, DIMENSION(SIZE(PCIS_SHAPE_V,1)) :: ZDN_DGAMMA_COL_IRR
!REAL, DIMENSION(SIZE(PZTM,1))         :: ZDTEMP_DTIME
!
!INTEGER, DIMENSION(SIZE(PZTM,1))      :: IDX_CHEN_M
INTEGER, DIMENSION(SIZE(PZTT,1))      :: IDX_CHEN_T
INTEGER                               :: II  ! index

! Chen 90
REAL, DIMENSION(121)          :: ZTEMP_CHEN ! Lookup table T (Chen and Lamb, 1994)
REAL, DIMENSION(SIZE(PZTM,1)) :: ZGAM_CHEN_T
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZGAM_CHEN_M
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDGAMMA_DT
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDGAMMA
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDTEMP    ! temperature variation between t and t-dt
REAL, DIMENSION(121)          :: ZGAM_CHEN ! Lookup table gamma(T) (Chen and Lamb, 1994)
                                           ! Index_GAM_CHEN = MAX( 1,MIN( 121,INT((XTT-T)*4.0)+1 ) )
                                           ! XTT = 273.16 and T is the temperature
!REAL, DIMENSION(121)          :: ZDGAMDT_CHEN ! Lookup table of d(gamma(T))/dT
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDGAMDT_CHEN_T ! Lookup table of d(gamma(T))/dT
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZZDGAMDT_CHEN_T ! Lookup table of d(gamma(T))/dT
!
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDN_GAMMA_COL_IRR !
!REAL, DIMENSION(SIZE(PZTM,1)) :: ZDN_GAMMA_PLA_IRR !
!
REAL, DIMENSION(SIZE(PZTM,1)) :: ZDGAMMA_COL_IRR !
REAL, DIMENSION(SIZE(PZTM,1)) :: ZDGAMMA_PLA_IRR !
!
REAL, DIMENSION(SIZE(PZTM,1)) :: ZW1, ZW2
!
REAL                          :: ZGAMMA_COL, ZGAMMA_PLA, ZGAMMA_IRR
!REAL, DIMENSION(2)            :: ZDGAMMA_FIXE 
! 
!-----------------------------------------------------------------------------
!
!*       0.   PRELIMINARIES
!             -------------
!
! Change unit of temperature
!
ZTT_CELSIUS(:) = PZTT(:) - XTT
!
!
!*       1.   COMPUTE GAMMA(T) 
!             ----------------
!
! (marine : à mettre dans modd_lima_cold_mixed ?)
!
!*       1.1  Table for Gamma and temperature (Chen and Lamb, 1984)
!
ZTEMP_CHEN(:) = (/ &
& 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, &
& 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, &  
& 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, &
& 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75, 10.00, 10.25, 10.50, 10.75, & 
& 11.00, 11.25, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00, 13.25, 13.50, &
& 13.75, 14.00, 14.25, 14.50, 14.75, 15.00, 15.25, 15.50, 15.75, 16.00, 16.25, &  
& 16.50, 16.75, 17.00, 17.25, 17.50, 17.75, 18.00, 18.25, 18.50, 18.75, 19.00, &
& 19.25, 19.50, 19.75, 20.00, 20.25, 20.50, 20.75, 21.00, 21.25, 21.50, 21.75, & 
& 22.00, 22.25, 22.50, 22.75, 23.00, 23.25, 23.50, 23.75, 24.00, 24.25, 24.50, &
& 24.75, 25.00, 25.25, 25.50, 25.75, 26.00, 26.25, 26.50, 26.75, 27.00, 27.25, &  
& 27.50, 27.75, 28.00, 28.25, 28.50, 28.75, 29.00, 29.25, 29.50, 29.75, 30.00 /)
!
ZTEMP_CHEN(:) = -ZTEMP_CHEN(:)
!
ZGAM_CHEN(1:121) = (/ &
& 1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.83, 0.81, &
& 0.78, 0.76, 0.70, 0.54, 0.47, 0.52, 0.63, 0.81, 1.10, 1.48, &
& 1.91, 2.09, 2.29, 2.40, 2.45, 2.43, 2.37, 2.29, 2.14, 2.00, &
& 1.86, 1.74, 1.62, 1.51, 1.40, 1.29, 1.19, 1.10, 1.00, 0.92, &
& 0.85, 0.79, 0.72, 0.67, 0.62, 0.58, 0.54, 0.50, 0.47, 0.44, &
& 0.41, 0.38, 0.35, 0.33, 0.32, 0.30, 0.29, 0.29, 0.28, 0.28, &
& 0.28, 0.28, 0.28, 0.29, 0.29, 0.30, 0.31, 0.32, 0.33, 0.35, &
& 0.37, 0.39, 0.43, 0.46, 0.49, 0.52, 0.56, 0.61, 0.66, 0.72, &
& 0.79, 0.86, 0.95, 1.05, 1.15, 1.26, 1.38, 1.50, 1.60, 1.70, &
& 1.78, 1.84, 1.88, 1.91, 1.91, 1.88, 1.86, 1.84, 1.80, 1.74, &
& 1.70, 1.64, 1.58, 1.55, 1.51, 1.48, 1.45, 1.43, 1.41, 1.39, &
& 1.38, 1.36, 1.35, 1.34, 1.33, 1.32, 1.31, 1.30, 1.29, 1.29, &
& 1.28 /) ! Chen data from 0 to -30C and given each -0.25C
!
!
!*       1.2  Find indexes for Gamma
!
DO II = 1, SIZE(PZTT,1)
  IDX_CHEN_T(II) = MINLOC(ABS(ZTEMP_CHEN(:)-ZTT_CELSIUS(II)),1)
  IF (ZTT_CELSIUS(II) .GT. ZTEMP_CHEN(IDX_CHEN_T(II)) .AND. IDX_CHEN_T(II) .GT. 1) THEN
    IDX_CHEN_T(II) = IDX_CHEN_T(II) - 1
  END IF
END DO
!
! linear interpolation to find Gamma given the temperature
!
ZGAM_CHEN_T(:) = ZGAM_CHEN(IDX_CHEN_T(:))
!
WHERE (IDX_CHEN_T(:) .LT. 121 .AND. IDX_CHEN_T(:) .GT. 0)
  ZGAM_CHEN_T(:) = ZGAM_CHEN(IDX_CHEN_T(:)) +                                &
                              (ZTT_CELSIUS(:) - ZTEMP_CHEN(IDX_CHEN_T(:))) * &
                  (ZGAM_CHEN(IDX_CHEN_T(:)+1) - ZGAM_CHEN(IDX_CHEN_T(:))) /  &
                 (ZTEMP_CHEN(IDX_CHEN_T(:)+1) - ZTEMP_CHEN(IDX_CHEN_T(:)) )
END WHERE
!
WHERE (IDX_CHEN_T(:) .EQ. SIZE(ZTEMP_CHEN,1))
  ZGAM_CHEN_T(:) = ZGAM_CHEN(SIZE(ZTEMP_CHEN,1))
END WHERE
!
! Gamma values for plates, columns and irregulars (from Figure 3, Chen and Lamb, 1984)
ZGAMMA_PLA = 0.5
ZGAMMA_COL = 1.5
ZGAMMA_IRR = 1.0
!
!
!*       2.   RATE OF TRANSFERT BETWEEN SHAPES
!             --------------------------------
!
! proportion of colums to transfer to irregulars
ZFACTOR_COL_IRR(:) = 0.
WHERE (HSHAPE_NAME_T(:) == 'PLA' .AND. PCIS_SHAPE_V(:,2) .GT. PCTMIN(4))
  ZFACTOR_COL_IRR(:) = (ZGAM_CHEN_T(:) - ZGAMMA_IRR) / (ZGAMMA_PLA - ZGAMMA_IRR)
END WHERE
!
! the factor must be between 0 and 1
ZFACTOR_COL_IRR(:) = MAX(ZFACTOR_COL_IRR(:), 0.0) 
ZFACTOR_COL_IRR(:) = MIN(ZFACTOR_COL_IRR(:), 1.0)
!
! proportion of plates to transfer to irregulars
ZFACTOR_PLA_IRR(:) = 0.
WHERE (HSHAPE_NAME_T(:) == 'COL' .AND. PCIS_SHAPE_V(:,1) .GT. PCTMIN(4))
  ZFACTOR_PLA_IRR(:) = (ZGAM_CHEN_T(:) - ZGAMMA_IRR) / (ZGAMMA_COL - ZGAMMA_IRR)
END WHERE
!
! the factor must be between 0 and 1
ZFACTOR_PLA_IRR(:) = MAX(ZFACTOR_PLA_IRR(:), 0.0)
ZFACTOR_PLA_IRR(:) = MIN(ZFACTOR_PLA_IRR(:), 1.0)
!
!
!
!*       4.   CHANGE IN NUMBER CONCENTRATION
!             ------------------------------
! 
ZW1(:) = PCIS_SHAPE_V(:,1) * ZFACTOR_PLA_IRR(:)
ZW2(:) = PCIS_SHAPE_V(:,2) * ZFACTOR_COL_IRR(:)
PCIS_SHAPE_V(:,1) = PCIS_SHAPE_V(:,1) - ZW1(:)
PCIS_SHAPE_V(:,2) = PCIS_SHAPE_V(:,2) - ZW2(:)
PCIS_SHAPE_V(:,3) = PCIS_SHAPE_V(:,3) + ZW1(:) + ZW2(:)
!
END SUBROUTINE LIMA_CHANGE_SHAPE
