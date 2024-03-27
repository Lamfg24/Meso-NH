!        ############################
         MODULE MODI_LIMA_AEROOPT_GET
!        ############################

INTERFACE

   SUBROUTINE LIMA_AEROOPT_GET(PSVTA,PZZ,PRHODREFA,                    &
                               PPIZA_WVL,PCGA_WVL,PTAUREL_WVL,PTAU550, &
                               KSWB,PIR,PII                            )
   
    REAL,    DIMENSION(:,:,:,:), INTENT(IN)    :: PSVTA       ! [moments/molec_{air}] transported moments aerosol
    REAL,    DIMENSION(:,:,:),   INTENT(IN)    :: PZZ         ! [m] height of layers
    REAL,    DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREFA   ! [kg/m3] density of air
    REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PPIZA_WVL   ! [-] single scattering albedo aerosol layer for all SW wavelengths
    REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PCGA_WVL    ! [-] assymetry factor faerosol layer for all SW wavelengths
    REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PTAUREL_WVL ! [-] opt.depth/opt.depth(550) faerosol layer for all SW wvl 
    REAL,    DIMENSION(:,:,:),   INTENT(INOUT) :: PTAU550     ! [-] opt.depth at 550nm for aaerosol layer
    INTEGER,                     INTENT(IN)    :: KSWB        ! [nbr] number of shortwave wavelengths
    REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PIR         ! [-] Real part of the aerosol refractive index
    REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PII         ! [-] Imaginary part of the aerosol refractive index

  END SUBROUTINE LIMA_AEROOPT_GET
 !
END INTERFACE
 !
END MODULE MODI_LIMA_AEROOPT_GET

!#####################################################################
   SUBROUTINE LIMA_AEROOPT_GET(PSVTA,PZZ,PRHODREFA,                    &
                               PPIZA_WVL,PCGA_WVL,PTAUREL_WVL,PTAU550, &
                               KSWB,PIR,PII                            )
!#####################################################################
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!      Benjamin Aouizerats (CNRM/GMEI)
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CH_AEROSOL
USE MODD_DUST
USE MODD_PARAM_LIMA
USE MODD_NSV
IMPLICIT NONE 
!
! Arguments
!
REAL,    DIMENSION(:,:,:,:), INTENT(IN)    :: PSVTA       ! [moments/molec_{air}] transported moments aerosol
REAL,    DIMENSION(:,:,:),   INTENT(IN)    :: PZZ         ! [m] height of layers
REAL,    DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREFA   ! [kg/m3] density of air
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PPIZA_WVL   ! [-] single scattering albedo aerosol layer for all SW wavelengths
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PCGA_WVL    ! [-] assymetry factor faerosol layer for all SW wavelengths
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PTAUREL_WVL ! [-] opt.depth/opt.depth(550) faerosol layer for all SW wvl 
REAL,    DIMENSION(:,:,:),   INTENT(INOUT) :: PTAU550     ! [-] opt.depth at 550nm for aaerosol layer
INTEGER,                     INTENT(IN)    :: KSWB        ! [nbr] number of shortwave wavelengths
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PIR         ! [-] Real part of the aerosol refractive index
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PII         ! [-] Imaginary part of the aerosol refractive index
!
! Local variables
!
INTEGER :: NMODE_AER
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                           :: ZSVT
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                         :: ZMASS         ![kg/m3] mass of aerosol mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                           :: ZMASSeq       ![kg/m3] mass of aerosol mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                           :: ZRADIUS       ![um] number median radius ofaerosol mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                           :: ZSIGMA        ![-] dispersion coefficientaerosol mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)                           :: ZTAU550_MDE   ![-] opt.depth 550nm one mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                         :: ZTAU_WVL_MDE  ![-] opt.depth @ wvl, one mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                         :: ZPIZA_WVL_MDE ![-] single scattering albedo @ wvl, one mode
REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:)                         :: ZCGA_WVL_MDE  ![-] assymetry factor @ wvl, one mode
INTEGER                                                         :: JMDE          ![idx] counter for modes
INTEGER                                                         :: JWVL          ![idx] counter for wavelengths
COMPLEX, DIMENSION(6,6)                                         :: Ri            !Refactive  index
COMPLEX, DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: eps1          !Intermediate computations
COMPLEX, DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: eps2          !Intermediate computations
COMPLEX, DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: eps3          !Intermediate computations
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: f1            !Intermediate computations
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VOC           !Volume of aerosol
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VDDST
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VBC
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA1
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA2
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA3
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA4
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA5
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA6
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA7
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA8
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA9
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSOA10
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VAM
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VNI
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VH2O
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VSU
REAL,    DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3))   :: VEXTR
COMPLEX, DIMENSION(size(PSVTA,1),size(PSVTA,2),size(PSVTA,3),6) :: Req           !Equivalent refractive index
REAL, PARAMETER                                                 :: EPSILON=1.d-8 ![um] a small number used to avoid zero

INTEGER ::JJJ

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
! Computations for LIMA :
!     IFN modes may be mixed, and composed of up to 4 species (DM1, DM2, BC, O)
!     Dust are accounted for separately, -> BC and O
!     For the radiation scheme, IFN of the same kind are regrouped
!
!     If CCN type "M", -> sea salt, not considered here 
!     If CCN type "C", -> mix of Su, BC, OC 
!     IF CCN from MACC analyses, mode1=sea-salt, mode2=Su, mode3=BC/OC
!
NMODE_AER=0
IF (LWARM .AND. NMOD_CCN.GE.1) NMODE_AER = NMODE_AER + NMOD_CCN
IF (LCOLD .AND. NMOD_IFN.GE.1) NMODE_AER = NMODE_AER + 2 
JJJ=0
IF (LWARM .AND. NMOD_CCN.GE.1) JJJ = NMOD_CCN

! Allocate arrays which size depend on number of modes
!
ALLOCATE(ZMASS(SIZE(PSVTA,1),SIZE(PSVTA,2),SIZE(PSVTA,3),16,NMODE_AER))
ALLOCATE(ZMASSeq(SIZE(PSVTA,1),SIZE(PSVTA,2),SIZE(PSVTA,3),NMODE_AER))
ALLOCATE(ZRADIUS(SIZE(PSVTA,1),SIZE(PSVTA,2),SIZE(PSVTA,3),NMODE_AER))
ALLOCATE(ZSIGMA(SIZE(PSVTA,1),SIZE(PSVTA,2),SIZE(PSVTA,3),NMODE_AER))
ZMASS(:,:,:,:,:)=0.
ZMASSeq(:,:,:,:)=0.

ALLOCATE(ZSVT(SIZE(PSVTA,1),SIZE(PSVTA,2),SIZE(PSVTA,3),NMODE_AER))
ZSVT(:,:,:,:) = 0.
IF (LWARM .AND. NMOD_CCN.GE.1) THEN
   DO JMDE = 1, NMOD_CCN
      ZSVT(:,:,:,JMDE)=PSVTA(:,:,:,NSV_LIMA_CCN_FREE+JMDE-1)   
   END DO
END IF
IF (LCOLD .AND. NMOD_IFN.GE.1) THEN
   DO JMDE=1, NMOD_IFN
      ZSVT(:,:,:,JJJ+1) = ZSVT(:,:,:,JJJ+1) + PSVTA(:,:,:,NSV_LIMA_IFN_FREE+JMDE-1)*XFRAC(3,JMDE)
      ZSVT(:,:,:,JJJ+2) = ZSVT(:,:,:,JJJ+2) + PSVTA(:,:,:,NSV_LIMA_IFN_FREE+JMDE-1)*XFRAC(4,JMDE)
   END DO
END IF
! Conversion ZSVT en #/m3
!
DO JMDE = 1, NMODE_AER
   ZSVT(:,:,:,JMDE) = ZSVT(:,:,:,JMDE) * PRHODREFA(:,:,:)
END DO

ALLOCATE(ZTAU550_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),NMODE_AER)) 
ALLOCATE(ZTAU_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_AER))
ALLOCATE(ZPIZA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_AER))
ALLOCATE(ZCGA_WVL_MDE(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB,NMODE_AER))

! Prepare r and sigma arrays (NMOD_CCN, DM1, DM2, BC, O)
!
IF (LWARM .AND. NMOD_CCN.GE.1) THEN
   DO JMDE = 1, NMOD_CCN
      ZSIGMA(:,:,:,JMDE)=EXP(XLOGSIG_CCN(JMDE))
      ZRADIUS(:,:,:,JMDE)=XR_MEAN_CCN(JMDE)*1.E6
   END DO
END IF
IF (LCOLD .AND. NMOD_IFN.GE.1) THEN
   DO JMDE = 1, 2
      ZSIGMA(:,:,:,JJJ+JMDE)=XSIGMA_IFN(2+JMDE)
      ZRADIUS(:,:,:,JJJ+JMDE)=XMDIAM_IFN(2+JMDE)*1.E6/2.
   END DO
END IF

! Compute CCN aerosol mass mixing ratio (kg/m3) 
!
IF (LWARM .AND. NMOD_CCN.GE.1) THEN
   IF ( CCCN_MODES=='MACC' .OR. CINT_MIXING=='MACC' ) THEN
      ! Mode 1 CCN - Sea salt  -> not in here
      ! Mode 2 CCN - Sulfates
      ZMASS(:,:,:,1,2) = 3.1415/6. * XRHO_CCN(2) * ZSVT(:,:,:,2) &
           * EXP( 3.*LOG(2.*XR_MEAN_CCN(2)) + 4.5*XLOGSIG_CCN(2)**2 )
      ! Mode 3 CCN - BC + O
      ZMASS(:,:,:,6,3) = 1./2. * 3.1415/6. * XRHO_CCN(3) * ZSVT(:,:,:,3) &
           * EXP( 3.*LOG(2.*XR_MEAN_CCN(3)) + 4.5*XLOGSIG_CCN(3)**2 )
      ZMASS(:,:,:,5,3) = ZMASS(:,:,:,6,3)
   ELSE
      DO JMDE = 1, NMOD_CCN
         IF (HTYPE_CCN(JMDE) == 'C') THEN
           ! Mix of Su, BC and OC
           ZMASS(:,:,:,1,JMDE) = 1./3. * 3.1415/6. * XRHO_CCN(JMDE) * ZSVT(:,:,:,JMDE) &
              * EXP( 3.*LOG(2.*XR_MEAN_CCN(JMDE)) + 4.5*XLOGSIG_CCN(JMDE)**2 )
           ZMASS(:,:,:,5,JMDE) = ZMASS(:,:,:,1,JMDE)
           ZMASS(:,:,:,6,JMDE) = ZMASS(:,:,:,1,JMDE)
         END IF
!        IF (HTYPE_CCN(JMDE) == 'M') THEN
!      ! Sea salt -> not in here
!        END IF
      END DO
   END IF
END IF

! Compute IFN aerosol mass mixing ratio (kg/m3)
!
IF (LCOLD .AND. NMOD_IFN.GE.1) THEN
   ! Dust 1 -> not in here
   ! Dust 2 -> not in here
   ! BC
   ZMASS(:,:,:, 6,JJJ+1) = 3.1415/6. * XRHO_IFN(3) * ZSVT(:,:,:,1+JJJ) &
        * EXP( 3.*LOG(XMDIAM_IFN(3)) + 4.5*LOG(XSIGMA_IFN(3))**2 )
   ! O
   ZMASS(:,:,:, 5,JJJ+2) = 3.1415/6. * XRHO_IFN(4) * ZSVT(:,:,:,2+JJJ) &
        * EXP( 3.*LOG(XMDIAM_IFN(4)) + 4.5*LOG(XSIGMA_IFN(4))**2 )
END IF

! Convert the aerosol mass from kg/m3 to ug/m3
!
!ZMASS(:,:,:,:,:)=ZMASS(:,:,:,:,:)*1.E9   


  DO JMDE=1,NMODE_AER

      Ri(1,1)=CMPLX(1.80,-7.40E-1)
      Ri(1,2)=CMPLX(1.80,-7.40E-1)
      Ri(1,3)=CMPLX(1.83,-7.40E-1)
      Ri(1,4)=CMPLX(1.88,-6.90E-1)
      Ri(1,5)=CMPLX(1.97,-6.80E-1)
      Ri(1,6)=CMPLX(2.10,-7.20E-1)

      Ri(2,1)=CMPLX(1.45,-1.00E-3)
      Ri(2,2)=CMPLX(1.45,-1.00E-3)
      Ri(2,3)=CMPLX(1.45,-1.00E-3)
      Ri(2,4)=CMPLX(1.46,-1.00E-3)
      Ri(2,5)=CMPLX(1.49,-1.00E-3)
      Ri(2,6)=CMPLX(1.42,-1.26E-2)

      Ri(3,1)=CMPLX(1.36,-3.60E-8)
      Ri(3,2)=CMPLX(1.34,-3.00E-9)
      Ri(3,3)=CMPLX(1.33,-1.80E-8)
      Ri(3,4)=CMPLX(1.33,-5.75E-7)
      Ri(3,5)=CMPLX(1.31,-1.28E-4)
      Ri(3,6)=CMPLX(1.42,-2.54E-1)
      
      Ri(4,1)=CMPLX(1.52,-5.00E-4)
      Ri(4,2)=CMPLX(1.52,-5.00E-4)
      Ri(4,3)=CMPLX(1.52,-5.00E-4)
      Ri(4,4)=CMPLX(1.52,-5.00E-4)
      Ri(4,5)=CMPLX(1.51,-5.00E-4)
      Ri(4,6)=CMPLX(1.35,-1.40E-2)
      
      Ri(5,1)=CMPLX(1.53,-5.00E-3)
      Ri(5,2)=CMPLX(1.53,-5.00E-3)
      Ri(5,3)=CMPLX(1.53,-6.00E-3)
      Ri(5,4)=CMPLX(1.52,-1.30E-2)
      Ri(5,5)=CMPLX(1.52,-1.30E-2)
      Ri(5,6)=CMPLX(1.45,-5.00E-1)


      Ri(6,1)=CMPLX(1.448,-0.00292)
      Ri(6,2)=CMPLX(1.448,-0.00292)
      Ri(6,3)=CMPLX(1.4777,-0.01897)
      Ri(6,4)=CMPLX(1.44023,-0.00116)
      Ri(6,5)=CMPLX(1.41163,-0.00106)
      Ri(6,6)=CMPLX(1.41163,-0.00106)

      VDDST(:,:,:)=0.

      IF(.NOT.ALLOCATED(XFAC)) ALLOCATE(XFAC(16))
      XFAC(:)=1800.*1.E-9
      XFAC(4)=1000.*1.E-9

      VOC(:,:,:)=(ZMASS(:,:,:,5,JMDE))/XFAC(5)
      VH2O(:,:,:)=(ZMASS(:,:,:,4,JMDE))/XFAC(4)
      VAM(:,:,:)=(ZMASS(:,:,:,3,JMDE))/XFAC(3)
      VSU(:,:,:)=(ZMASS(:,:,:,1,JMDE))/XFAC(1)
      VNI(:,:,:)=(ZMASS(:,:,:,2,JMDE))/XFAC(2)
      VBC(:,:,:)=(ZMASS(:,:,:,6,JMDE))/XFAC(6)
      VSOA1(:,:,:)=(ZMASS(:,:,:,7,JMDE))/XFAC(7)
      VSOA2(:,:,:)=(ZMASS(:,:,:,8,JMDE))/XFAC(8)
      VSOA3(:,:,:)=(ZMASS(:,:,:,9,JMDE))/XFAC(9)
      VSOA4(:,:,:)=(ZMASS(:,:,:,10,JMDE))/XFAC(10)
      VSOA5(:,:,:)=(ZMASS(:,:,:,11,JMDE))/XFAC(11)
      VSOA6(:,:,:)=(ZMASS(:,:,:,12,JMDE))/XFAC(12)
      VSOA7(:,:,:)=(ZMASS(:,:,:,13,JMDE))/XFAC(13)
      VSOA8(:,:,:)=(ZMASS(:,:,:,14,JMDE))/XFAC(14)
      VSOA9(:,:,:)=(ZMASS(:,:,:,15,JMDE))/XFAC(15)
      VSOA10(:,:,:)=(ZMASS(:,:,:,16,JMDE))/XFAC(16)
      VSOA(:,:,:)=VSOA1(:,:,:)+VSOA2(:,:,:)+VSOA3(:,:,:)+VSOA4(:,:,:)+&
                  VSOA5(:,:,:)+VSOA6(:,:,:)+VSOA7(:,:,:)+VSOA8(:,:,:)+&
                  VSOA9(:,:,:)+VSOA10(:,:,:)

      VEXTR(:,:,:)=VSOA(:,:,:)+VH2O(:,:,:)+VAM(:,:,:)+VSU(:,:,:)+VNI(:,:,:)
     
     
     DO JWVL=1,KSWB                    !Number of SW wavelengths

     eps1(:,:,:)=CMPLX(0)
     WHERE ( (VBC(:,:,:)+VOC(:,:,:)).NE.0. )
     eps1(:,:,:)=CMPLX((Ri(1,JWVL)*VBC(:,:,:)+Ri(2,JWVL)*VOC(:,:,:)+VDDST(:,:,:)*Ri(6,JWVL))/(VBC(:,:,:)+VOC(:,:,:)))**2
     END WHERE

     Req(:,:,:,JWVL)=sqrt(CMPLX(eps1(:,:,:)))

     WHERE (VEXTR(:,:,:).NE.0. )
     eps2(:,:,:)=CMPLX((VSOA(:,:,:)*Ri(2,JWVL)+VH2O(:,:,:)*Ri(3,JWVL)+VAM(:,:,:)*Ri(4,JWVL)&
                 +VSU(:,:,:)*Ri(4,JWVL)+VNI(:,:,:)*Ri(5,JWVL))/&
                 (VSOA(:,:,:)+VH2O(:,:,:)+VAM(:,:,:)+VSU(:,:,:)+VNI(:,:,:)))**2
     f1(:,:,:)=(VOC(:,:,:)+VBC(:,:,:))/(VSOA(:,:,:)+VH2O(:,:,:)+VAM(:,:,:)+VSU(:,:,:)+VNI(:,:,:)+VOC(:,:,:)+VBC(:,:,:))
     eps3(:,:,:)=CMPLX(eps2(:,:,:)*(eps1(:,:,:)+2*eps2(:,:,:)+2*f1(:,:,:)*(eps1(:,:,:)-eps2(:,:,:)))/&
                      (eps1(:,:,:)+2*eps2(:,:,:)-f1(:,:,:)*(eps1(:,:,:)-eps2(:,:,:))))
     Req(:,:,:,JWVL)=sqrt(CMPLX(eps3(:,:,:)))
     ENDWHERE

     ENDDO   
        ZMASSeq(:,:,:,JMDE)=ZMASS(:,:,:,1,JMDE)+ZMASS(:,:,:,2,JMDE)+ZMASS(:,:,:,3,JMDE)&
                         +ZMASS(:,:,:,4,JMDE)+ZMASS(:,:,:,5,JMDE)+ZMASS(:,:,:,6,JMDE)+ZMASS(:,:,:,7,JMDE)&
                         +ZMASS(:,:,:,8,JMDE)+ZMASS(:,:,:,9,JMDE)+ZMASS(:,:,:,10,JMDE)+ZMASS(:,:,:,11,JMDE)&
                         +ZMASS(:,:,:,12,JMDE)+ZMASS(:,:,:,13,JMDE)+ZMASS(:,:,:,14,JMDE)+ZMASS(:,:,:,15,JMDE)&
                         +ZMASS(:,:,:,16,JMDE)    
     PII(:,:,:,:) = aimag(CMPLX(Req(:,:,:,:))) 
     PIR(:,:,:,:) = real(CMPLX(Req(:,:,:,:))) 
     !Get aerosol optical properties from look up tables


     CALL AEROOPT_LKT(                     &
          ZRADIUS(:,:,:,JMDE)              &  !I [um] number median radius for current mode
          ,ZSIGMA(:,:,:,JMDE)              &  !I [none] dispersion coefficient for current mode
          ,ZMASSeq(:,:,:, JMDE)            &  !I [kg/m3] Mass of aerosol for current mode
          ,ZTAU550_MDE(:,:,:,JMDE)         &  !O [-] optical depth at 550 nm wavelength
          ,ZTAU_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] opt.depth(lambda)/opt.depth(550nm)
          ,ZPIZA_WVL_MDE(:,:,:,:,JMDE)     &  !O [-] single scattering coefficient at any wavelength
          ,ZCGA_WVL_MDE(:,:,:,:,JMDE)      &  !O [-] assymetry factor at any wavelength
          ,PZZ(:,:,:)                      &  !I [m] height of layers
          ,KSWB                            &  !I [nbr] number of shortwave bands
          ,Req(:,:,:,:))

 
  ENDDO  !Loop on modes


  !Erase earlier value of optical depth at 550 nm
  PTAU550(:,:,:)=0.d0


  !Get total at 550 nm from all modes 
  DO JMDE=1,NMODE_AER
     PTAU550(:,:,:) =        &        !Aerosol optical depth at 550 nm for all aerosol
          PTAU550(:,:,:)     &        !Aerosol optical depth at 550 nm for all aerosol
          + ZTAU550_MDE(:,:,:,JMDE)   !Optical depth for one mode at 550 nm      
  ENDDO
  !Initialize output variables
  PTAUREL_WVL(:,:,:,:)=0.d0         !Initialize opt.depth at wvl=lambda
  PCGA_WVL(:,:,:,:)=0.d0           !Initialize assym.factor at wvl=lambda
  PPIZA_WVL(:,:,:,:)=0.d0          !Initialize single scattering albedo at wvl=lambda


  !Find the numerator in the expression for the average of the optical properties   
  DO JMDE=1,NMODE_AER             !Number of modes
     DO JWVL=1,KSWB                    !Number of SW wavelengths

        !Get sum of optical depth from all modes at wvl
        PTAUREL_WVL(:,:,:,JWVL)  =                &  !new opt.depth(lambda) / opt.depth(550)
             PTAUREL_WVL(:,:,:,JWVL)              &  !old sum for all modes at wvl=lambda
             +ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)       !optical depth for one mode at wvl=lambda
        
        !Get sum of all assymmetry factors  from all modes at wvl=lambda 
        PCGA_WVL(:,:,:,JWVL) =                     &  !New sum of assymetry factors
             PCGA_WVL(:,:,:,JWVL)                  &  !old sum of assymetry factors
             +ZCGA_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Assymetry factor for one mode and one wavelength
             *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)     &  !Optical depth of this wavelength and mode
             *ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)       !Fraction of radiation scattered

        !Get sum of single scattering albdedo at wvl=lambda
        PPIZA_WVL(:,:,:,JWVL)  =                  & !New sum of single scattering albedo
             PPIZA_WVL(:,:,:,JWVL)             & !Old sum of single scattering albedo
             +ZPIZA_WVL_MDE(:,:,:,JWVL,JMDE)   & !SSA for onen mode and one wavelength
             *ZTAU_WVL_MDE(:,:,:,JWVL,JMDE)        !Optical depth for this wavelength and mode
        
     ENDDO
  ENDDO
  !Compute the output values for aerosol optical properties
  DO JWVL=1,KSWB
     
     !Divide total single scattering albdeo by total optical depth at this wavelength
     !This is needed since we weight all single scattering alebdos by wavelengths just above
     PPIZA_WVL(:,:,:,JWVL) =               &     !The value we want is ....
          PPIZA_WVL(:,:,:,JWVL)            &     !..value weighted by optical depths of all wvl and modes
          /max(epsilon,PTAUREL_WVL(:,:,:,JWVL))            !..divided by the optical depth for all wvl 
     
     
     !Divide total assymetry factor by total optical depth at this wavelength
     !This is needed since we weight all assymetry factors by wavelengths just above
     PCGA_WVL(:,:,:,JWVL) =              &  !The value we want is ....
          PCGA_WVL(:,:,:,JWVL)           &  !..value weighted by optical depths of all wvl and modes
          /                              &
          (max(epsilon,                  &
          (PTAUREL_WVL(:,:,:,JWVL)       &  !..divided scattered fraction of by the optical depth
          *PPIZA_WVL(:,:,:,JWVL))))

     !Finally convert PTAUREL_WVL which was until now an optical depth to a fraction of optical depth
     PTAUREL_WVL(:,:,:,JWVL) =                    &
          PTAUREL_WVL(:,:,:,JWVL)                 &  !Opt.depth at lambda with contr. from all modes
          /max(epsilon,PTAU550(:,:,:))               !Optical depth at 550 contr. from all modes
     
  ENDDO !Loop on wavelenghts
     

  !DEALLOCATE local  arrays which size depend on number of modes
  DEALLOCATE(ZTAU550_MDE)
  DEALLOCATE(ZTAU_WVL_MDE)
  DEALLOCATE(ZPIZA_WVL_MDE)
  DEALLOCATE(ZCGA_WVL_MDE)
!
  DEALLOCATE(ZMASS)
  DEALLOCATE(ZMASSeq)
  DEALLOCATE(ZRADIUS)
  DEALLOCATE(ZSIGMA)
  DEALLOCATE(ZSVT)

 CONTAINS

       SUBROUTINE AEROOPT_LKT(              &
       PRG                             & !I [um] number median radius of aerosol mode
       ,PSIGMA                         & !I [-] lognormal dispersion coefficient
       ,PMASS                          & !I [kg/m3] Mass concentration of aerosol
       ,PTAU550                        & !O [optical depth at 550 nm
       ,PTAU_WVL                       & !O [-] opt.depth(lambda)/opt.depth(550nm)
       ,PPIZA_WVL                      & !O [-] single scattering coefficient at any wavelength
       ,PCGA_WVL                       & !O [-] assymetry factor at any wavelength
       ,PZZ                            & !I [m] height of layers
       ,KSWB                           & !I [nbr] number of short wave bands
       ,Req)

    !Purpose: Get optical properties of one aerosol mode from the mass concentration, 
    !dispersion coefficient and number median radius.

    !Use the module with the aerosol optical properties look up tables
    USE MODD_AEROSET  

    IMPLICIT NONE
    !INPUT
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PRG          !I [um] number median radius for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PSIGMA       !I [-] dispersion coefficient for one mode
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PMASS        !I [ug/m3] mass of aerosol
    REAL, DIMENSION(:,:,:), INTENT(IN)        :: PZZ          !I [m] height of layers
    COMPLEX, DIMENSION(:,:,:,:), INTENT(IN)        :: Req          !
    INTEGER, INTENT(IN)                       :: KSWB         !I [nbr] number of shortwave bands

    !OUTPUT
    REAL, DIMENSION(:,:,:), INTENT(OUT)       :: PTAU550      !O [-] optical depth at 550 nm
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PTAU_WVL     !O [-] optical depth at wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PPIZA_WVL    !O [-] single scattering albedo @ wvl
    REAL, DIMENSION(:,:,:,:), INTENT(OUT)     :: PCGA_WVL     !O [-] assymetry factor @ wvl

    !LOCALS
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3),KSWB) :: ZEXT_COEFF_WVL    ![m2/kg] Extinction coefficient at wvl
    REAL, DIMENSION(SIZE(PTAU550,1),SIZE(PTAU550,2),SIZE(PTAU550,3)     ) :: ZEXT_COEFF_550    ![m2/kg] Extinction coefficient at 550nm
    REAL                      :: FACT_SIGMA  ![-] factor needed to get right index in look up table for sigma
    REAL                      :: FACT_RADIUS ![-] factor needed to get right index in look up table for radius
    INTEGER                   :: WVL_IDX     ![idx] counter for wavelengths
    INTEGER                   :: JI, JJ, JK  ![idx] counters for lon, lat and lev
    INTEGER                   :: RG_IDX      ![idx] index for radius to get in look up table
    INTEGER                   :: SG_IDX      ![idx] index for sigma to get in look up table
    INTEGER                   :: JRAD      ![idx] index for sigma to get in look up table
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZRG     ![um] bounded value for number median radius
    REAL,DIMENSION(SIZE(PRG,1),SIZE(PRG,2),SIZE(PRG,3)) :: ZSIGMA  ![um] bounded value for sigma
    REAL, PARAMETER           :: EPSILON=1.d-8                     ![um] a small number used to avoid zero
    REAL                      :: ZRADIUS_LKT_MAX, ZRADIUS_LKT_MIN  ![um] values limited at midpoint values of bin
    REAL                      :: ZSIGMA_LKT_MAX, ZSIGMA_LKT_MIN    ![-] values limited at midpoint of bin
    INTEGER                   :: JKRAD,II,IS                             !Index valid for radiation code
    REAL,DIMENSION(300)         ::      RADI,DR,ND
    REAL                        :: LWC,PI,DELTA

    COMPLEX,DIMENSION(13,6)                   :: R
    REAL,dimension(13) ::Err
    real,dimension(1) ::pds2,pds1,temp1,temp2,lg1,lg2,lg3,pds3,pds4,pds5,pds6

    REAL::TAU1,TAU2,TAU3,TAU4,TAU5,TAU6,TAU7,TAU8,SSA1,SSA2,&
    SSA3,SSA4,SSA5,SSA6,SSA7,SSA8,GA1,GA2,GA3,GA4,GA5,GA6,GA7,GA8
    REAL,DIMENSION(13)::PTAU1,PTAU2,PTAU3,PTAU4,PTAU5,PTAU6,PTAU7,PTAU8,PSSA1,PSSA2,PSSA3,PSSA4,PSSA5,&
    PSSA6,PSSA7,PSSA8,PGA1,PGA2,PGA3,PGA4,PGA5,PGA6,PGA7,PGA8
    integer,dimension(1)::INDri,INDrr,INDsi    
    real,dimension (10)::SI ,Dsi  
    real,dimension (8)::RR   ,Drr
    real,dimension (6)::RI   ,Dri
    !Limit max and min values to be midpont of bin to avoid 0 or NMAX+1 values
    ZRADIUS_LKT_MAX=20.
    ZRADIUS_LKT_MIN=1.E-12
    ZSIGMA_LKT_MAX=2.55
    ZSIGMA_LKT_MIN=1.0


    !Remove unphysical values for rg
    ZRG(:,:,:) = min( max(ZRADIUS_LKT_MIN,PRG(:,:,:)), ZRADIUS_LKT_MAX)
    ZSIGMA(:,:,:) = min( max(ZSIGMA_LKT_MIN,PSIGMA(:,:,:)), ZSIGMA_LKT_MAX)

    !Initilalize arrays to make sure, they are intent(OUT), 
    !and may be initialized strangely by the computer
    PTAU550(:,:,:)=EPSILON 
    PTAU_WVL(:,:,:,:)=EPSILON     
    PPIZA_WVL(:,:,:,:)=EPSILON
    PCGA_WVL(:,:,:,:)=EPSILON

    PI=ACOS(-1.0)

      DO IS=1,10
        SI(IS)=0.85+IS*0.2
      ENDDO

      RR(1)=1.45
      RR(2)=1.50
      RR(3)=1.55
      RR(4)=1.60
      RR(5)=1.65
      RR(6)=1.70
      RR(7)=1.75
      RR(8)=1.80

      RI(1)=-0.001
      RI(2)=-0.006
      RI(3)=-0.008
      RI(4)=-0.020
      RI(5)=-0.100
      RI(6)=-0.400

    DO WVL_IDX = 1,KSWB
          
    
        DO JK=2,SIZE(PMASS,3)-1
          JKRAD = JK - 1  !Index in radiation code
          DO JJ=1,SIZE(PMASS,2)
             DO JI=1,SIZE(PMASS,1)

  
                  Drr(:)=abs((RR(:))-real(Req(JI,JJ,JK,WVL_IDX)))
                  Dri(:)=abs((RI(:))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                  Dsi(:)=abs((SI(:))-PSIGMA(JI,JJ,JK))

                  INDrr=minloc(Drr(:))
                  INDri=minloc(Dri(:))
                  INDsi=minloc(Dsi(:))

                  TAU1=0.
                  TAU2=0.
                  TAU3=0.
                  TAU4=0.
                  TAU5=0.
                  TAU6=0.
                  TAU7=0.
                  TAU8=0.

                  SSA1=0.
                  SSA2=0.
                  SSA3=0.
                  SSA4=0.
                  SSA5=0.
                  SSA6=0.
                  SSA7=0.
                  SSA8=0.

                  GA1=0.
                  GA2=0.
                  GA3=0.
                  GA4=0.
                  GA5=0.
                  GA6=0.
                  GA7=0.
                  GA8=0.

                  pds1=0.
                  pds2=0.
                  pds3=0.
                  pds4=0.
                  pds5=0.
                  pds6=0.

                  !SWITCH RR,RI
                  if(SI(INDsi(1)).gt.PSIGMA(JI,JJ,JK).and.INDsi(1).ne.1)  THEN
 
                     if( RR(INDrr(1)).gt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.1) THEN
                        
                         if (RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)-1<SI<INDsi(1)


                                

                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3



                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)-1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)

                                else 
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)

                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)-1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)-1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
  
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)+1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)+1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1)+1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                         else
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))

                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1

                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)-1,INDri(1),:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*TAU1+pds6(1)*TAU2))+&
                                        pds2(1)*((pds5(1)*TAU5+pds6(1)*TAU6))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*SSA1+pds6(1)*SSA2))+&
                                        pds2(1)*((pds5(1)*SSA5+pds6(1)*SSA6))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds5(1)*GA1+pds6(1)*GA2)+&
                                        pds2(1)*((pds5(1)*GA5+pds6(1)*GA6))
                         endif

                        
                  else if (RR(INDrr(1)).lt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.8) THEN
                         if( RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)-1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)-1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)-1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                        
                       
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)+1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)+1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1)+1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                         else
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))

                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1

                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif



                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif



                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif


                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1)+1,INDri(1),:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*TAU1+pds6(1)*TAU2))+&
                                        pds2(1)*((pds5(1)*TAU5+pds6(1)*TAU6))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*SSA1+pds6(1)*SSA2))+&
                                        pds2(1)*((pds5(1)*SSA5+pds6(1)*SSA6))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*GA1+pds6(1)*GA2))+&
                                        pds2(1)*((pds5(1)*GA5+pds6(1)*GA6))
                         endif

                  else
                         if( RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)-1<SI<INDsi(1)


                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif


                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))
                                        
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))
                       
                        
                         else if( RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))   

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))
                         else
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)-1<SI<INDsi(1)
                                 pds5=abs(SI(INDsi(1)-1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif


                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)-1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*TAU1+pds6(1)*TAU2
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*SSA1+pds6(1)*SSA2
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*GA1+pds6(1)*GA2
                         endif

                      endif

                else if (SI(INDsi(1)).lt.PSIGMA(JI,JJ,JK).and.INDsi(1).ne.10) THEN

                       if (RR(INDrr(1)).gt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.1)THEN
                        
                         if (RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)<SI<INDsi(1)+1

                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)-1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)-1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)-1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                                 

                         else if( RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)<SI<INDsi(1)+1
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)+1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)+1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1)+1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                         else
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)<SI<INDsi(1)+1
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif


                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)-1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif


                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*TAU1+pds6(1)*TAU2))+&
                                        pds2(1)*((pds5(1)*TAU5+pds6(1)*TAU6))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*SSA1+pds6(1)*SSA2))+&
                                        pds2(1)*((pds5(1)*SSA5+pds6(1)*SSA6))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*GA1+pds6(1)*GA2))+&
                                        pds2(1)*((pds5(1)*GA5+pds6(1)*GA6))
                         endif

                        
                  else if( RR(INDrr(1)).lt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.8) THEN
                         if( RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)<SI<INDsi(1)+1
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)-1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)-1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)-1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))
                        
                       
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)<SI<INDsi(1)+1
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                  PTAU8(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)+1,:)
                                if (PTAU8(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU8=PTAU8(1)*LOG(PRG(JI,JJ,JK))**5+PTAU8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(3)*LOG(PRG(JI,JJ,JK))**3+PTAU8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(5)*LOG(PRG(JI,JJ,JK))+PTAU8(6)
                                else
                                 TAU8=PTAU8(7)*LOG(PRG(JI,JJ,JK))**5+PTAU8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU8(9)*LOG(PRG(JI,JJ,JK))**3+PTAU8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU8(11)*LOG(PRG(JI,JJ,JK))+PTAU8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                  PSSA8(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)+1,:)
                                if (PSSA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA8=PSSA8(1)*LOG(PRG(JI,JJ,JK))**5+PSSA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(3)*LOG(PRG(JI,JJ,JK))**3+PSSA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(5)*LOG(PRG(JI,JJ,JK))+PSSA8(6)
                                else
                                 SSA8=PSSA8(7)*LOG(PRG(JI,JJ,JK))**5+PSSA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA8(9)*LOG(PRG(JI,JJ,JK))**3+PSSA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA8(11)*LOG(PRG(JI,JJ,JK))+PSSA8(12)
                                endif
                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                                  PGA8(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1)+1,:)
                                if (PGA8(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA8=PGA8(1)*LOG(PRG(JI,JJ,JK))**5+PGA8(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(3)*LOG(PRG(JI,JJ,JK))**3+PGA8(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(5)*LOG(PRG(JI,JJ,JK))+PGA8(6)
                                else
                                 GA8=PGA8(7)*LOG(PRG(JI,JJ,JK))**5+PGA8(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA8(9)*LOG(PRG(JI,JJ,JK))**3+PGA8(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA8(11)*LOG(PRG(JI,JJ,JK))+PGA8(12)
                                endif
                           
     
                                if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (SSA8.GT.1.0) SSA8=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
                                if (GA7.GT.1.0) GA7=1.0
                                if (GA8.GT.1.0) GA8=1.0
                                
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*TAU5+pds6(1)*TAU6)+pds4(1)*(pds5(1)*TAU7+pds6(1)*TAU8))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*SSA5+pds6(1)*SSA6)+pds4(1)*(pds5(1)*SSA7+pds6(1)*SSA8))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))+&
                                        pds2(1)*(pds3(1)*(pds5(1)*GA5+pds6(1)*GA6)+pds4(1)*(pds5(1)*GA7+pds6(1)*GA8))

                         else
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)<SI<INDsi(1)+1
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))

                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1

                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU6(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PTAU6(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU6=PTAU6(1)*LOG(PRG(JI,JJ,JK))**5+PTAU6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(3)*LOG(PRG(JI,JJ,JK))**3+PTAU6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(5)*LOG(PRG(JI,JJ,JK))+PTAU6(6)
                                else
                                 TAU6=PTAU6(7)*LOG(PRG(JI,JJ,JK))**5+PTAU6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU6(9)*LOG(PRG(JI,JJ,JK))**3+PTAU6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU6(11)*LOG(PRG(JI,JJ,JK))+PTAU6(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA6(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PSSA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA6=PSSA6(1)*LOG(PRG(JI,JJ,JK))**5+PSSA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(3)*LOG(PRG(JI,JJ,JK))**3+PSSA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(5)*LOG(PRG(JI,JJ,JK))+PSSA6(6)
                                else
                                 SSA6=PSSA6(7)*LOG(PRG(JI,JJ,JK))**5+PSSA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA6(9)*LOG(PRG(JI,JJ,JK))**3+PSSA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA6(11)*LOG(PRG(JI,JJ,JK))+PSSA6(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA6(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1)+1,INDri(1),:)
                                if (PGA6(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA6=PGA6(1)*LOG(PRG(JI,JJ,JK))**5+PGA6(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(3)*LOG(PRG(JI,JJ,JK))**3+PGA6(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(5)*LOG(PRG(JI,JJ,JK))+PGA6(6)
                                else
                                 GA6=PGA6(7)*LOG(PRG(JI,JJ,JK))**5+PGA6(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA6(9)*LOG(PRG(JI,JJ,JK))**3+PGA6(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA6(11)*LOG(PRG(JI,JJ,JK))+PGA6(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA6.GT.1.0) SSA6=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA6.GT.1.0) GA6=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*TAU1+pds6(1)*TAU2))+&
                                        pds2(1)*((pds5(1)*TAU5+pds6(1)*TAU6))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*SSA1+pds6(1)*SSA2))+&
                                        pds2(1)*((pds5(1)*SSA5+pds6(1)*SSA6))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*((pds5(1)*GA1+pds6(1)*GA2))+&
                                        pds2(1)*((pds5(1)*GA5+pds6(1)*GA6))
                         endif

                  else
                         if( RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)<SI<INDsi(1)+1

                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)-1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))
                        
                       
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)<SI<INDsi(1)+1

                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2
                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU4(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PTAU4(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU4=PTAU4(1)*LOG(PRG(JI,JJ,JK))**5+PTAU4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(3)*LOG(PRG(JI,JJ,JK))**3+PTAU4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(5)*LOG(PRG(JI,JJ,JK))+PTAU4(6)
                                else
                                 TAU4=PTAU4(7)*LOG(PRG(JI,JJ,JK))**5+PTAU4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU4(9)*LOG(PRG(JI,JJ,JK))**3+PTAU4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU4(11)*LOG(PRG(JI,JJ,JK))+PTAU4(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA4(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PSSA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA4=PSSA4(1)*LOG(PRG(JI,JJ,JK))**5+PSSA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(3)*LOG(PRG(JI,JJ,JK))**3+PSSA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(5)*LOG(PRG(JI,JJ,JK))+PSSA4(6)
                                else
                                 SSA4=PSSA4(7)*LOG(PRG(JI,JJ,JK))**5+PSSA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA4(9)*LOG(PRG(JI,JJ,JK))**3+PSSA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA4(11)*LOG(PRG(JI,JJ,JK))+PSSA4(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA4(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1)+1,:)
                                if (PGA4(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA4=PGA4(1)*LOG(PRG(JI,JJ,JK))**5+PGA4(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(3)*LOG(PRG(JI,JJ,JK))**3+PGA4(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(5)*LOG(PRG(JI,JJ,JK))+PGA4(6)
                                else
                                 GA4=PGA4(7)*LOG(PRG(JI,JJ,JK))**5+PGA4(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA4(9)*LOG(PRG(JI,JJ,JK))**3+PGA4(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA4(11)*LOG(PRG(JI,JJ,JK))+PGA4(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA4.GT.1.0) SSA4=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA4.GT.1.0) GA4=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*TAU1+pds6(1)*TAU2)+pds4(1)*(pds5(1)*TAU3+pds6(1)*TAU4))
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*SSA1+pds6(1)*SSA2)+pds4(1)*(pds5(1)*SSA3+pds6(1)*SSA4))
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=(pds3(1)*(pds5(1)*GA1+pds6(1)*GA2)+pds4(1)*(pds5(1)*GA3+pds6(1)*GA4))
                         else
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)<SI<INDsi(1)+1

                                 pds5=abs(SI(INDsi(1)+1)-PSIGMA(JI,JJ,JK))
                                 pds6=abs(SI(INDsi(1))-PSIGMA(JI,JJ,JK))                                

                                 lg3=pds5+pds6
                                 pds5=pds5/lg3
                                 pds6=pds6/lg3

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU2(:)=POLYTAU(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PTAU2(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU2=PTAU2(1)*LOG(PRG(JI,JJ,JK))**5+PTAU2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(3)*LOG(PRG(JI,JJ,JK))**3+PTAU2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(5)*LOG(PRG(JI,JJ,JK))+PTAU2(6)
                                else
                                 TAU2=PTAU2(7)*LOG(PRG(JI,JJ,JK))**5+PTAU2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU2(9)*LOG(PRG(JI,JJ,JK))**3+PTAU2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU2(11)*LOG(PRG(JI,JJ,JK))+PTAU2(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA2(:)=POLYSSA(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PSSA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA2=PSSA2(1)*LOG(PRG(JI,JJ,JK))**5+PSSA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(3)*LOG(PRG(JI,JJ,JK))**3+PSSA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(5)*LOG(PRG(JI,JJ,JK))+PSSA2(6)
                                else
                                 SSA2=PSSA2(7)*LOG(PRG(JI,JJ,JK))**5+PSSA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA2(9)*LOG(PRG(JI,JJ,JK))**3+PSSA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA2(11)*LOG(PRG(JI,JJ,JK))+PSSA2(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA2(:)=POLYG(WVL_IDX,INDsi(1)+1,INDrr(1),INDri(1),:)
                                if (PGA2(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA2=PGA2(1)*LOG(PRG(JI,JJ,JK))**5+PGA2(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(3)*LOG(PRG(JI,JJ,JK))**3+PGA2(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(5)*LOG(PRG(JI,JJ,JK))+PGA2(6)
                                else
                                 GA2=PGA2(7)*LOG(PRG(JI,JJ,JK))**5+PGA2(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA2(9)*LOG(PRG(JI,JJ,JK))**3+PGA2(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA2(11)*LOG(PRG(JI,JJ,JK))+PGA2(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA2.GT.1.0) SSA2=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA2.GT.1.0) GA2=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*TAU1+pds6(1)*TAU2
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*SSA1+pds6(1)*SSA2
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds5(1)*GA1+pds6(1)*GA2
                         endif

                      endif

                  else

                      if (RR(INDrr(1)).gt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.1) THEN
                        
                         if (RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)=SI

                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                          
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2


                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA7.GT.1.0) GA7=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*TAU1+pds4(1)*TAU3)+&
                                        pds2(1)*(pds3(1)*TAU5+pds4(1)*TAU7)
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*SSA1+pds4(1)*SSA3)+&
                                        pds2(1)*(pds3(1)*SSA5+pds4(1)*SSA7)
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*GA1+pds4(1)*GA3)+&
                                        pds2(1)*(pds3(1)*GA5+pds4(1)*GA7)                                    
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)=SI

                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                          
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1                               

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA7.GT.1.0) GA7=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*TAU1+pds4(1)*TAU3)+&
                                        pds2(1)*(pds3(1)*TAU5+pds4(1)*TAU7)
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*SSA1+pds4(1)*SSA3)+&
                                        pds2(1)*(pds3(1)*SSA5+pds4(1)*SSA7)
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*GA1+pds4(1)*GA3)+&
                                        pds2(1)*(pds3(1)*GA5+pds4(1)*GA7)
                         else
                                 !CAS INDrr(1)-1<RR<INDrr(1)
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)=SI
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                            
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif


                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA5.GT.1.0) GA5=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*TAU1+pds2(1)*TAU5
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*SSA1+pds2(1)*SSA5
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*GA1+pds2(1)*GA5
                         endif

                        
                  else if (RR(INDrr(1)).lt.real(Req(JI,JJ,JK,WVL_IDX)).and.INDrr(1).ne.8) THEN
                         if( RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)=SI
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                            
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)-1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA7.GT.1.0) GA7=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*TAU1+pds4(1)*TAU3)+&
                                        pds2(1)*(pds3(1)*TAU5+pds4(1)*TAU7)
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*SSA1+pds4(1)*SSA3)+&
                                        pds2(1)*(pds3(1)*SSA5+pds4(1)*SSA7)
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*GA1+pds4(1)*GA3)+&
                                  pds2(1)*(pds3(1)*GA5+pds4(1)*GA7)
                        
                       
                        
                         else if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)=SI
                                 pds1=abs(RR(INDrr(1)+1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds3=abs(RI(INDri(1)+1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1
                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                  PTAU7(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PTAU7(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU7=PTAU7(1)*LOG(PRG(JI,JJ,JK))**5+PTAU7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(3)*LOG(PRG(JI,JJ,JK))**3+PTAU7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(5)*LOG(PRG(JI,JJ,JK))+PTAU7(6)
                                else
                                 TAU7=PTAU7(7)*LOG(PRG(JI,JJ,JK))**5+PTAU7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU7(9)*LOG(PRG(JI,JJ,JK))**3+PTAU7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU7(11)*LOG(PRG(JI,JJ,JK))+PTAU7(12)
                                endif

                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                  PSSA7(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PSSA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA7=PSSA7(1)*LOG(PRG(JI,JJ,JK))**5+PSSA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(3)*LOG(PRG(JI,JJ,JK))**3+PSSA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(5)*LOG(PRG(JI,JJ,JK))+PSSA7(6)
                                else
                                 SSA7=PSSA7(7)*LOG(PRG(JI,JJ,JK))**5+PSSA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA7(9)*LOG(PRG(JI,JJ,JK))**3+PSSA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA7(11)*LOG(PRG(JI,JJ,JK))+PSSA7(12)
                                endif

                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)+1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                                  PGA7(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)+1,INDri(1)+1,:)
                                if (PGA7(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA7=PGA7(1)*LOG(PRG(JI,JJ,JK))**5+PGA7(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(3)*LOG(PRG(JI,JJ,JK))**3+PGA7(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(5)*LOG(PRG(JI,JJ,JK))+PGA7(6)
                                else
                                 GA7=PGA7(7)*LOG(PRG(JI,JJ,JK))**5+PGA7(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA7(9)*LOG(PRG(JI,JJ,JK))**3+PGA7(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA7(11)*LOG(PRG(JI,JJ,JK))+PGA7(12)
                                endif

                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (SSA7.GT.1.0) SSA7=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
                                if (GA5.GT.1.0) GA5=1.0
                                if (GA7.GT.1.0) GA7=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*TAU1+pds4(1)*TAU3)+&
                                        pds2(1)*(pds3(1)*TAU5+pds4(1)*TAU7)
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*SSA1+pds4(1)*SSA3)+&
                                        pds2(1)*(pds3(1)*SSA5+pds4(1)*SSA7)
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*(pds3(1)*GA1+pds4(1)*GA3)+&
                                        pds2(1)*(pds3(1)*GA5+pds4(1)*GA7)
                         else
                                 !CAS INDrr(1)<RR<INDrr(1)+1
                                 !CAS INDri(1)=RI
                                 !CAS INDsi(1)=SI
                                 pds1=abs(RR(INDrr(1)-1)-real(Req(JI,JJ,JK,WVL_IDX)))
                                 pds2=abs(RR(INDrr(1))-real(Req(JI,JJ,JK,WVL_IDX)))
                            
                                 lg1=pds2+pds1
                                 pds1=pds1/lg1
                                 pds2=pds2/lg1

                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif


                                  PTAU5(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PTAU5(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU5=PTAU5(1)*LOG(PRG(JI,JJ,JK))**5+PTAU5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(3)*LOG(PRG(JI,JJ,JK))**3+PTAU5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(5)*LOG(PRG(JI,JJ,JK))+PTAU5(6)
                                else
                                 TAU5=PTAU5(7)*LOG(PRG(JI,JJ,JK))**5+PTAU5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU5(9)*LOG(PRG(JI,JJ,JK))**3+PTAU5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU5(11)*LOG(PRG(JI,JJ,JK))+PTAU5(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA5(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PSSA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA5=PSSA5(1)*LOG(PRG(JI,JJ,JK))**5+PSSA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(3)*LOG(PRG(JI,JJ,JK))**3+PSSA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(5)*LOG(PRG(JI,JJ,JK))+PSSA5(6)
                                else
                                 SSA5=PSSA5(7)*LOG(PRG(JI,JJ,JK))**5+PSSA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA5(9)*LOG(PRG(JI,JJ,JK))**3+PSSA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA5(11)*LOG(PRG(JI,JJ,JK))+PSSA5(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA5(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1)-1,INDri(1),:)
                                if (PGA5(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA5=PGA5(1)*LOG(PRG(JI,JJ,JK))**5+PGA5(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(3)*LOG(PRG(JI,JJ,JK))**3+PGA5(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(5)*LOG(PRG(JI,JJ,JK))+PGA5(6)
                                else
                                 GA5=PGA5(7)*LOG(PRG(JI,JJ,JK))**5+PGA5(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA5(9)*LOG(PRG(JI,JJ,JK))**3+PGA5(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA5(11)*LOG(PRG(JI,JJ,JK))+PGA5(12)
                                endif

                            
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA5.GT.1.0) SSA5=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA5.GT.1.0) GA5=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*TAU1+pds2(1)*TAU5
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*SSA1+pds2(1)*SSA5
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds1(1)*GA1+pds2(1)*GA5
                         endif

                  else

                         if (RI(INDri(1)).lt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.1) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)-1<RI<INDri(1)
                                 !CAS INDsi(1)=SI

                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2


                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif


                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif


                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif


                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*TAU1+pds4(1)*TAU3
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*SSA1+pds4(1)*SSA3
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*GA1+pds4(1)*GA3
                        endif
                        
                       
                        
                        if (RI(INDri(1)).gt.aimag(Req(JI,JJ,JK,WVL_IDX)).and.INDri(1).ne.6) THEN
                                 !CAS INDrr(1)=RR
                                 !CAS INDri(1)<RI<INDri(1)+1
                                 !CAS INDsi(1)=SI

                                 pds3=abs(RI(INDri(1)-1)-aimag(Req(JI,JJ,JK,WVL_IDX)))
                                 pds4=abs(RI(INDri(1))-aimag(Req(JI,JJ,JK,WVL_IDX)))
                         

                                 lg2=pds3+pds4
                                 pds3=pds3/lg2
                                 pds4=pds4/lg2


                                 PTAU1(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PTAU1(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU1=PTAU1(1)*LOG(PRG(JI,JJ,JK))**5+PTAU1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(3)*LOG(PRG(JI,JJ,JK))**3+PTAU1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(5)*LOG(PRG(JI,JJ,JK))+PTAU1(6)
                                else
                                 TAU1=PTAU1(7)*LOG(PRG(JI,JJ,JK))**5+PTAU1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU1(9)*LOG(PRG(JI,JJ,JK))**3+PTAU1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU1(11)*LOG(PRG(JI,JJ,JK))+PTAU1(12)
                                endif

                                  PTAU3(:)=POLYTAU(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PTAU3(13).GE.PRG(JI,JJ,JK)) THEN
                                 TAU3=PTAU3(1)*LOG(PRG(JI,JJ,JK))**5+PTAU3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(3)*LOG(PRG(JI,JJ,JK))**3+PTAU3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(5)*LOG(PRG(JI,JJ,JK))+PTAU3(6)
                                else
                                 TAU3=PTAU3(7)*LOG(PRG(JI,JJ,JK))**5+PTAU3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PTAU3(9)*LOG(PRG(JI,JJ,JK))**3+PTAU3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PTAU3(11)*LOG(PRG(JI,JJ,JK))+PTAU3(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!SSA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                   PSSA1(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PSSA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA1=PSSA1(1)*LOG(PRG(JI,JJ,JK))**5+PSSA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(3)*LOG(PRG(JI,JJ,JK))**3+PSSA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(5)*LOG(PRG(JI,JJ,JK))+PSSA1(6)
                                else
                                 SSA1=PSSA1(7)*LOG(PRG(JI,JJ,JK))**5+PSSA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA1(9)*LOG(PRG(JI,JJ,JK))**3+PSSA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA1(11)*LOG(PRG(JI,JJ,JK))+PSSA1(12)
                                endif

                                  PSSA3(:)=POLYSSA(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PSSA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 SSA3=PSSA3(1)*LOG(PRG(JI,JJ,JK))**5+PSSA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(3)*LOG(PRG(JI,JJ,JK))**3+PSSA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(5)*LOG(PRG(JI,JJ,JK))+PSSA3(6)
                                else
                                 SSA3=PSSA3(7)*LOG(PRG(JI,JJ,JK))**5+PSSA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PSSA3(9)*LOG(PRG(JI,JJ,JK))**3+PSSA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PSSA3(11)*LOG(PRG(JI,JJ,JK))+PSSA3(12)
                                endif

                                !!!!!!!!!!!!!!!!!!!!!!!G!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                    PGA1(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1),:)
                                if (PGA1(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA1=PGA1(1)*LOG(PRG(JI,JJ,JK))**5+PGA1(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(3)*LOG(PRG(JI,JJ,JK))**3+PGA1(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(5)*LOG(PRG(JI,JJ,JK))+PGA1(6)
                                else
                                 GA1=PGA1(7)*LOG(PRG(JI,JJ,JK))**5+PGA1(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA1(9)*LOG(PRG(JI,JJ,JK))**3+PGA1(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA1(11)*LOG(PRG(JI,JJ,JK))+PGA1(12)
                                endif

                                  PGA3(:)=POLYG(WVL_IDX,INDsi(1),INDrr(1),INDri(1)-1,:)
                                if (PGA3(13).GE.PRG(JI,JJ,JK)) THEN
                                 GA3=PGA3(1)*LOG(PRG(JI,JJ,JK))**5+PGA3(2)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(3)*LOG(PRG(JI,JJ,JK))**3+PGA3(4)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(5)*LOG(PRG(JI,JJ,JK))+PGA3(6)
                                else
                                 GA3=PGA3(7)*LOG(PRG(JI,JJ,JK))**5+PGA3(8)*LOG(PRG(JI,JJ,JK))**4+&
                                 PGA3(9)*LOG(PRG(JI,JJ,JK))**3+PGA3(10)*LOG(PRG(JI,JJ,JK))**2+&
                                 PGA3(11)*LOG(PRG(JI,JJ,JK))+PGA3(12)
                                endif
                               if (SSA1.GT.1.0) SSA1=1.0
                                if (SSA3.GT.1.0) SSA3=1.0
                                if (GA1.GT.1.0) GA1=1.0
                                if (GA3.GT.1.0) GA3=1.0
        ZEXT_COEFF_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*TAU1+pds4(1)*TAU3
        PPIZA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*SSA1+pds4(1)*SSA3
        PCGA_WVL(JI,JJ,JKRAD,WVL_IDX)=pds3(1)*GA1+pds4(1)*GA3
                        
                        endif

                         
                       
          endif



 
             ENDDO
          ENDDO
       ENDDO
    ENDDO


    ZEXT_COEFF_550(:,:,:) = ZEXT_COEFF_WVL(:,:,:,3)

    !Get the optical depth of this mode using the looked up extinction coeffient
    DO JK=2,SIZE(PZZ,3)-1
       JKRAD = JK - 1           !Index in radiation code
       PTAU550(:,:,JKRAD) = ZEXT_COEFF_550(:,:,JKRAD) &
            * PMASS(:,:,JK)                    &
            * (PZZ(:,:,JK+1) - PZZ(:,:,JK))
    ENDDO




    !Get the optical depth of whatever wavelength using the looked up tables
    DO WVL_IDX=1,KSWB
       DO JK=2,SIZE(PZZ,3)-1
          JKRAD = JK -1
          PTAU_WVL(:,:,JKRAD,WVL_IDX) =              &
               PMASS(:,:,JK)                      & ![kg/m3] Mass in this mode
               *ZEXT_COEFF_WVL(:,:,JKRAD,WVL_IDX)    & ![m2/kg] mass exinction coefficient
               *(PZZ(:,:,JK+1) - PZZ(:,:,JK))       ![m] Height of layer
       ENDDO !Loop on levels
    ENDDO    !Loop on wavelengths

    !Avoid unphysical values (which might occur on grid edges) on grid edges
    PTAU550(:,:,:)=max(PTAU550(:,:,:),EPSILON)
    PTAU_WVL(:,:,:,:)=max(PTAU_WVL(:,:,:,:),EPSILON)
    PPIZA_WVL(:,:,:,:)=max(PPIZA_WVL(:,:,:,:),EPSILON)
    PCGA_WVL(:,:,:,:)=max(PCGA_WVL(:,:,:,:),EPSILON)

  END SUBROUTINE AEROOPT_LKT

END SUBROUTINE LIMA_AEROOPT_GET
