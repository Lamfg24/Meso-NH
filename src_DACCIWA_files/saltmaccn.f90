!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/saltlfin.f90,v $ $Revision: 1.1.2.2.2.1.2.1 $
! MASDEV4_7 newsrc 2007/01/25 13:13:15
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_SALTMACC_n
!     ########################
!
INTERFACE
!
SUBROUTINE SALTMACC_n(PSV,PMASSMACC)
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),INTENT(INOUT) :: PSV
REAL,       DIMENSION(:,:,:,:),INTENT(IN) :: PMASSMACC
END SUBROUTINE SALTMACC_n
!
END INTERFACE
!
END MODULE MODI_SALTMACC_n
!
!
!     ############################################################
      SUBROUTINE SALTMACC_n(PSV,  PMASSMACC)
!     ############################################################
!
!!    PURPOSE
!!    -------
!!    Initialise le champs  de salts à partir des analyses MACC
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LACy)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SALT
USE MODD_NSV
USE MODD_CSTS_SALT
USE MODE_SALT_PSD
USE MODI_INIT_SALT
! 
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:,:,:), INTENT(INOUT) :: PSV
REAL,   DIMENSION(:,:,:,:),INTENT(IN) :: PMASSMACC ! macc salt concentration (kg.kg-1)
!
!
!*      0.2    declarations local variables
!
REAL   :: ZDEN2MOL, ZRHOI, ZMI, ZFAC, ZRGMIN
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZCTOTA
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM
REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  ::  ZMASS
INTEGER,DIMENSION(:),    ALLOCATABLE  :: IM0, IM3, IM6
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS, ZINISIGMA
INTEGER :: IKU, IMOMENTS
INTEGER :: JJ, JN, JK  ! loop counter
INTEGER :: IMODEIDX  ! index mode

REAL :: DELTA_1,DELTA_2,DELTA_3,DELTA_4,DELTA_5,DELTA_6,DELTA_7
REAL :: RATIO_1,RATIO_2,RATIO_3,RATIO_4,RATIO_5, RATIO_6,RATIO_7
REAL :: DELTA_MACC_1,DELTA_MACC_2,DELTA_MACC_3
REAL :: RAY_MACC_1,RAY_MACC_2,RAY_MACC_3,RAY_MACC_4
REAL :: RAY_2,RAY_3,RAY_4
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  ::  ZMASS_TEST
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               -----------------------------------
!
!        1.1    initialisation 
!
CALL INIT_SALT
IKU = SIZE(PSV,3)
!
ALLOCATE (IM0(NMODE_SLT))
ALLOCATE (IM3(NMODE_SLT))
ALLOCATE (IM6(NMODE_SLT))
ALLOCATE (ZCTOTA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT))
ALLOCATE (ZM(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT*3))
ALLOCATE (ZSIGMA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3)))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
ALLOCATE (ZINISIGMA(NMODE_SLT))
ALLOCATE (ZMMIN(NMODE_SLT*3))
ALLOCATE (ZMASS(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3),NMODE_SLT))
!
!
! Rayons des bins MACC

RAY_MACC_1 = 0.03
RAY_MACC_2 = 0.5
RAY_MACC_3 = 5
RAY_MACC_4 = 20

! Choix des diametres de separation (selon Ovadnevaite et al., 2014)

RAY_2 = 0.045
RAY_3 = 0.11
RAY_4 = 0.41

! Calcul des proportions 

! Calcul des écarts bin MACC

DELTA_MACC_1 = RAY_MACC_2 - RAY_MACC_1
DELTA_MACC_2 = RAY_MACC_3 - RAY_MACC_2
DELTA_MACC_3 = RAY_MACC_4 - RAY_MACC_3

! Calcul des ecarts par mode en fonction des rayons de separation
! puis calcul de la masse correspondante avec facteur correctif pour eviter
! la surestimation des concentrations en aerosols

DELTA_1 = RAY_2 - RAY_MACC_1
RATIO_1 = DELTA_1 / DELTA_MACC_1
ZMASS(:,:,:,2) = PMASSMACC(:,:,:,1) * RATIO_1 ! * 1E-2 ! Attribution Mode 2 ORILAM

DELTA_2 = RAY_3 - RAY_2
RATIO_2 = DELTA_2 / DELTA_MACC_1
ZMASS(:,:,:,3) = PMASSMACC(:,:,:,1) * RATIO_2 ! * 1E-2 ! Attribution Mode 3 ORILAM

DELTA_3 = RAY_4 - RAY_3
RATIO_3 = DELTA_3 / DELTA_MACC_1
ZMASS(:,:,:,4) = PMASSMACC(:,:,:,1) * RATIO_3 ! * 1E-1  ! Attribution Mode 4 ORILAM

DELTA_4 = RAY_MACC_2 - RAY_4
RATIO_4 = DELTA_4 / DELTA_MACC_1
ZMASS(:,:,:,5) = PMASSMACC(:,:,:,1) * RATIO_4 ! Attribution Mode 5 ORILAM

DELTA_5 = RAY_MACC_3 - RAY_MACC_2
RATIO_5 = DELTA_5 / DELTA_MACC_2
ZMASS(:,:,:,5) = (PMASSMACC(:,:,:,2) * RATIO_5) + ZMASS(:,:,:,5) ! Attribution Mode 5 bis ORILAM

DELTA_6 = 10 - RAY_MACC_3
RATIO_6 = DELTA_3 / DELTA_MACC_1
ZMASS(:,:,:,5) = (PMASSMACC(:,:,:,3) * RATIO_6) + ZMASS(:,:,:,5) ! Attribution Mode 5 ter ORILAM

ZMASS(:,:,:,5) = ZMASS(:,:,:,5) * 1E-1

! Hyp : the ultrafine mode is neglected for orilam-lima coupling
ZMASS(:,:,:,1) = PMASSMACC(:,:,:,1) * 1E-5 ! ultrafin mode
!
!========================================================
! Adjust the mass / SSA emissions after a few hours
ZMASS(:,:,:,1) = ZMASS(:,:,:,1) * 1.
ZMASS(:,:,:,2) = ZMASS(:,:,:,2) * 1.
ZMASS(:,:,:,3) = ZMASS(:,:,:,3) * 1.
ZMASS(:,:,:,4) = ZMASS(:,:,:,4) * 1.
ZMASS(:,:,:,5) = ZMASS(:,:,:,5) * 1.
!========================================================

DO JN = 1, NMODE_SLT
  IM0(JN) = 1 + (JN - 1) * 3
  IM3(JN) = 2 + (JN - 1) * 3
  IM6(JN) = 3 + (JN - 1) * 3
  !
  !Get the salt mode we are talking about, MODE 2 is treated first, then mode 3, then 1
  !This index is only needed to get the right radius out of the XINIRADIUS array and the
  !right XINISIG out of the XINISIG-array
  IMODEIDX = JPSALTORDER(JN)
  !
  !Convert initial mass median radius to number median radius
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF
  ZINISIGMA(JN)  = XINISIG_SLT(IMODEIDX)
  !
  ZMMIN(IM0(JN)) = XN0MIN_SLT(IMODEIDX)
  ZRGMIN   = ZINIRADIUS(JN)
  ZMMIN(IM3(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZMMIN(IM6(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**6)*EXP(18. * LOG(ZINISIGMA(JN))**2)

END DO

ZMASS(:,:,:,:) = MAX(ZMASS(:,:,:,:), 1E-40)
!
!
ZRHOI = XDENSITY_SALT !1.8e3 !++changed alfgr
ZMI   = XMOLARWEIGHT_SALT
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
ZFAC = (4. / 3.) * XPI * ZRHOI * 1.e-9

!
DO JN = 1, NMODE_SLT

!*       1.1    calculate moment 0 from ZMASS
!
  ZM(:,:,:,IM0(JN)) = ZMASS(:,:,:,JPSALTORDER(JN))       &![kg_{salt}/kg_{air}  
                    / XDENSITY_SALT                      &![kg__{salt}/m3_{salt}==>m3_{salt}/m3{air}
                    * (6.d0 / XPI)                       &
                    / (2.d0 * ZINIRADIUS(JN) * 1.d-6)**3 &![particle/m_salt^{-3}]==> particle/m3
                    * EXP(-4.5*(LOG(ZINISIGMA(JN)))**2) !Take into account distribution

  ZM(:,:,:,IM0(JN)) = MAX(ZMMIN(IM0(JN)), ZM(:,:,:,IM0(JN)))
!
!*       1.2    calculate moment 3 from m0,  RG and SIG 
!
  ZM(:,:,:,IM3(JN)) = ZM(:,:,:,IM0(JN)) * &
                      (ZINIRADIUS(JN)**3) * & 
                      EXP(4.5*LOG(ZINISIGMA(JN))**2) 

  ZM(:,:,:,IM3(JN)) = MAX(ZMMIN(IM3(JN)), ZM(:,:,:,IM3(JN)))
!
!*       1.3    calculate moment 6 from m0,  RG and SIG 
!
  ZM(:,:,:,IM6(JN))= ZM(:,:,:,IM0(JN)) * ((ZINIRADIUS(JN)**6)*&
                        EXP(18. * (LOG(ZINISIGMA(JN)))**2))
  ZM(:,:,:,IM6(JN)) = MAX(ZMMIN(IM6(JN)), ZM(:,:,:,IM6(JN)))
!
!*       1.4    output concentration
!
  IMOMENTS = INT(NSV_SLTEND - NSV_SLTBEG+1) / NMODE_SLT
  IF (IMOMENTS == 3) THEN
    PSV(:,:,:,1+(JN-1)*3) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO)
    PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,IM3(JN)) * XMD * XPI * 4. / 3. * ZRHOI / &
                            (ZMI*XM3TOUM3_SALT)

    PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,IM6(JN)) *  XMD / (XAVOGADRO*1.d-6)
  ELSE IF (IMOMENTS == 2) THEN
    PSV(:,:,:,1+(JN-1)*2) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO)
    PSV(:,:,:,2+(JN-1)*2) = ZM(:,:,:,IM3(JN)) * XMD * XPI * 4./3. * ZRHOI / &
                            (ZMI*XM3TOUM3_SALT)
  ELSE 
    PSV(:,:,:,JN) = ZM(:,:,:,IM3(JN)) * XMD * XPI * 4. / 3. * ZRHOI / &
                    (ZMI * XM3TOUM3_SALT)
  END IF
END DO

!
DEALLOCATE(ZMMIN)
DEALLOCATE(ZINISIGMA)
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZM)
DEALLOCATE(ZCTOTA)
DEALLOCATE(IM6)
DEALLOCATE(IM3)
DEALLOCATE(IM0)
DEALLOCATE(ZMASS)
!
!
END SUBROUTINE SALTMACC_n
