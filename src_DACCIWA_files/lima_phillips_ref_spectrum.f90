!      ######################################
       MODULE MODI_LIMA_PHILLIPS_REF_SPECTRUM
!      ######################################
!
INTERFACE
      SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM (PZT, PSI, PSI_W, PZY)
!
REAL, DIMENSION(:), INTENT(IN)    :: PZT    ! Temperature
REAL, DIMENSION(:), INTENT(IN)    :: PSI    ! Saturation over ice
REAL, DIMENSION(:), INTENT(IN)    :: PSI_W  ! Saturation over ice at water sat.
REAL, DIMENSION(:), INTENT(INOUT) :: PZY    ! Reference activity spectrum
!
END SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM
END INTERFACE
END MODULE MODI_LIMA_PHILLIPS_REF_SPECTRUM
!
!     ######################################################################
      SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM (PZT, PSI, PSI_W, PZY)
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the reference activation spectrum
!!    described by Phillips (2008)
!!
!!
!!    REFERENCE
!!    ---------
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
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XGAMMA, XRHO_CFDC
USE MODI_LIMA_FUNCTIONS,  ONLY : RECT, DELTA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:), INTENT(IN)    :: PZT    ! Temperature
REAL, DIMENSION(:), INTENT(IN)    :: PSI    ! Saturation over ice
REAL, DIMENSION(:), INTENT(IN)    :: PSI_W  ! Saturation over ice at water sat.
REAL, DIMENSION(:), INTENT(INOUT) :: PZY    ! Reference activity spectrum
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZMAX, &
                                     ZMOY, &
                                     ZZY1, &
                                     ZZY2, &
                                     Z1,   &
                                     Z2,   &
                                     ZSI2
!
REAL                              :: ZPSI  ! factor to match equations at -30C
!
!-------------------------------------------------------------------------------
!
ALLOCATE(ZMAX(SIZE(PZT))) ; ZMAX(:)= 0.0
ALLOCATE(ZMOY(SIZE(PZT))) ; ZMOY(:)= 0.0
ALLOCATE(ZZY1(SIZE(PZT))) ; ZZY1(:)= 0.0
ALLOCATE(ZZY2(SIZE(PZT))) ; ZZY2(:)= 0.0
ALLOCATE(Z1(SIZE(PZT)))   ; Z1(:)  = 0.0
ALLOCATE(Z2(SIZE(PZT)))   ; Z2(:)  = 0.0
ALLOCATE(ZSI2(SIZE(PZT))) ; ZSI2(:)= 0.0
!
PZY(:) = 0.0   
!
ZPSI   = 0.058707 * XGAMMA / XRHO_CFDC  ! Appendix A Phillips et al. (2008)
!
ZSI2(:) = MIN(PSI(:),PSI_W(:))
!
WHERE( PSI(:) > 1.0 )
!
!* T <= -35 C 
!
  PZY(:) = 1000. * XGAMMA / XRHO_CFDC                &  ! Eq. 1
         * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3 &  ! in
         * RECT(1.,0.,PZT(:),(XTT-80.),(XTT-35.))       ! Phillips et al. (2008)
!
!* -35 C < T <= -25 C (in Appendix A) 
!
  ZZY1(:) = 1000. * XGAMMA / XRHO_CFDC &
          * ( EXP(12.96*(MIN(ZSI2(:),7.)-1.1)) )**0.3
  ZZY2(:) = 1000. * ZPSI               &
          *   EXP(12.96*(MIN(ZSI2(:),7.)-1.0)-0.639)
!
!* -35 C < T <= -30 C
!
  ZMAX(:) = 1000. * XGAMMA / XRHO_CFDC         & ! Eq. 6
          * ( EXP(12.96*(PSI_W(:)-1.1)) )**0.3 &
          * RECT(1.,0.,PZT(:),(XTT-35.),(XTT-30.))
!
!* -30 C < T <= -25 C
!
  ZMAX(:) = ZMAX(:) + 1000. * ZPSI            & ! Eq. 7
          * EXP( 12.96*(PSI_W(:)-1.0)-0.639 ) &
          * RECT(1.,0.,PZT(:),(XTT-30.),(XTT-25.))
  Z1(:)   = MIN(ZZY1(:), ZMAX(:)) ! n~ in Appendix A
  Z2(:)   = MIN(ZZY2(:), ZMAX(:)) ! n^ in Appendix A
!
!* T > -25 C 
!
  PZY(:)  = PZY(:) + 1000. * ZPSI                    & ! Eq. 3
          * EXP( 12.96*(MIN(ZSI2(:),7.)-1.0)-0.639 ) &
          * RECT(1.,0.,PZT(:),(XTT-25.),(XTT-2.))
END WHERE
!
WHERE (Z2(:) > 0.0 .AND. Z1(:) > 0.0)
  ZMOY(:) = Z2(:) * (Z1(:) / Z2(:))**DELTA(1.,0.,PZT(:),(XTT-35.),(XTT-25.)) ! Eq. 5
  PZY(:)  = PZY(:) + MIN(ZMOY(:),ZMAX(:))  ! N_{IN,1,*} Eq. 4
END WHERE
!
DEALLOCATE(ZMAX)
DEALLOCATE(ZMOY)
DEALLOCATE(ZZY1)
DEALLOCATE(ZZY2)
DEALLOCATE(Z1)
DEALLOCATE(Z2)
!
END SUBROUTINE LIMA_PHILLIPS_REF_SPECTRUM
