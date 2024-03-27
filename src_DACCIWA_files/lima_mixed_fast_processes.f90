!      #####################################
       MODULE MODI_LIMA_MIXED_FAST_PROCESSES
!      #####################################
!
INTERFACE
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (PRHODREF, PZT, PPRES, PTSTEP,           &
                                            PLSFACT, PLVFACT, PKA, PDV, PCJ,        &
                                            PRVT, PRCT, PRRT, PRIT, PRST, PRGT,     &
                                            PRHT, PCCT, PCRT, PCIT,                 &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,     &
                                            PTHS, PCCS, PCRS, PCIS, PCIS_SHAPE,     &
                                            PCSHAPE_S,                              &
                                            PLBDAC, PLBDAR, PLBDAS, PLBDAG, PLBDAH, &
                                            PRHODJ, GMICRO, KMI, PTHS_3D,           &
                                            PRCS_3D, PRRS_3D, PRIS_3D, PRSS_3D,     &
                                            PRGS_3D, PRHS_3D,                       &
                                            PCCS_3D, PCRS_3D, PCIS_3D,              &
                                            PCIS_SHAPE_3D)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: PZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: PDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRHS    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCRS    ! Rain water conc. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS    ! Pristine ice conc. source
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAH  ! Slope param of the hail distr.
!
! used for budget storage
LOGICAL, DIMENSION(:,:,:),   INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ
INTEGER,                     INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PTHS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRRS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRIS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRSS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRGS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PRHS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCCS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCRS_3D
REAL,    DIMENSION(:,:,:),   INTENT(IN) :: PCIS_3D
REAL,    DIMENSION(:,:,:,:), INTENT(IN) :: PCIS_SHAPE_3D
REAL,    DIMENSION(:,:),     INTENT(INOUT) :: PCIS_SHAPE   ! Pristine ice conc. per shape 
REAL,    DIMENSION(:,:),     INTENT(INOUT) :: PCSHAPE_S    ! Pristine ice conc. per shape 
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_FAST_PROCESSES
!
!     ###############################################################################
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (PRHODREF, PZT, PPRES, PTSTEP,           &
                                            PLSFACT, PLVFACT, PKA, PDV, PCJ,        &
                                            PRVT, PRCT, PRRT, PRIT, PRST, PRGT,     &
                                            PRHT, PCCT, PCRT, PCIT,             &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,     &
                                            PTHS, PCCS, PCRS, PCIS, PCIS_SHAPE, &
                                            PCSHAPE_S,                              &
                                            PLBDAC, PLBDAR, PLBDAS, PLBDAG, PLBDAH, &
                                            PRHODJ, GMICRO, KMI, PTHS_3D,           &
                                            PRCS_3D, PRRS_3D, PRIS_3D, PRSS_3D,     &
                                            PRGS_3D, PRHS_3D,                       &
                                            PCCS_3D, PCRS_3D, PCIS_3D,          &
                                            PCIS_SHAPE_3D)                       
!     ###############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    fast processes :
!!      
!!      - Fast RS processes :
!!          - Cloud droplet riming of the aggregates
!!          - Hallett-Mossop ice multiplication process due to snow riming
!!          - Rain accretion onto the aggregates
!!          - Conversion-Melting of the aggregates
!!
!!      - Fast RG processes :
!!          - Rain contact freezing
!!          - Wet/Dry growth of the graupel
!!          - Hallett-Mossop ice multiplication process due to graupel riming
!!          - Melting of the graupeln
!!
!!
!!**  METHOD
!!    ------
!!
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
!!      C. Barthe  * LACy *  jan. 2014    add budgets
!!      T. Hoarau  * LACy *  Jul. 2016    add CIBU process
!!      JP Pinty   * LA *    Feb. 2017    correction of HMG flag
!!      M. Claeys  * LACy *  mar. 2019    add ice crystal schapes
!!      JP Pinty   * LA*     jul. 2019    add RDSF process
!!      C. Barthe  * LACy *  jul. 2019    attribute shapes to secondary produced ice crystals
!!      C. Barthe  * LACy *  sep. 2019    change the transfert rate to each shape for dryg, wetg, cfrz and wetg
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM, ONLY : XDR
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
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: PDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRHS    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCRS    ! Rain water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS    ! Pristine ice conc. source
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PCIS_SHAPE_3D   ! 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCIS_SHAPE  ! Pristine ice conc. per shape 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PCSHAPE_S   ! Pristine ice conc. per shape 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAH  ! Slope param of the hail distr.
!
! used for budget storage
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRRS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRSS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCRS_3D
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS_3D
!
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PZT)) :: GRIM, GACC, GDRY, GWET, GHAIL ! Test where to compute
INTEGER :: IGRIM, IGACC, IGDRY, IGWET, IHAIL
INTEGER :: JJ
INTEGER :: JSH, JL 
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(PZT))  :: ZZW, ZZX      
REAL,    DIMENSION(SIZE(PZT))  :: ZRDRYG, ZRWETG   
REAL,    DIMENSION(SIZE(PZT),7)  :: ZZW1 
REAL :: NHAIL
REAL :: ZTHRH, ZTHRC
!
! Variables for CIBU
LOGICAL, DIMENSION(SIZE(PZT)) :: GCIBU ! Test where to compute collision process
LOGICAL, SAVE                 :: GFIRSTCALL = .TRUE. ! control switch for the first call
!
INTEGER                            :: ICIBU
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_S1,IVEC2_S2         ! Snow indice vector
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_G                   ! Graupel indice vector
INTEGER, PARAMETER                 :: I_SEED_PARAM = 26032012
INTEGER, DIMENSION(:), ALLOCATABLE :: I_SEED
INTEGER                            :: NI_SEED
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_S,ZVEC1_S1,ZVEC1_S2,  & ! Work vectors
                                      ZVEC1_S3,ZVEC1_S4,          &
                                      ZVEC1_S11,ZVEC1_S12,        & ! for snow
                                      ZVEC1_S21,ZVEC1_S22,        &
                                      ZVEC1_S31,ZVEC1_S32,        &
                                      ZVEC1_S41,ZVEC1_S42,        &
                                      ZVEC2_S1,ZVEC2_S2
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_G,ZVEC1_G1,ZVEC1_G2, & ! Work vectors
                                      ZVEC2_G                      ! for graupel
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_SNOW_1, & ! incomplete gamma function
                                      ZINTG_SNOW_2, & ! for snow
                                      ZINTG_SNOW_3, &
                                      ZINTG_SNOW_4
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_GRAUPEL_1, & ! incomplete gamma
                                      ZINTG_GRAUPEL_2    ! function for graupel
REAL,    DIMENSION(:), ALLOCATABLE :: ZNI_CIBU,ZRI_CIBU  ! CIBU rates
REAL,    DIMENSION(:), ALLOCATABLE :: ZFRAGMENTS, ZHARVEST, ZFRAG_CIBU
REAL                               :: ZFACT1_XNDEBRIS, ZFACT2_XNDEBRIS
!
LOGICAL, DIMENSION(SIZE(PZT))      :: GRDSF              ! Test where to compute collision process
INTEGER :: IRDSF
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R            ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R1           ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC2_R            ! Work vectors for rain                            
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_R            ! Rain indice vector
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_RAIN         ! incomplete gamma function for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZNI_RDSF,ZRI_RDSF  ! RDSF rates
!
REAL,    DIMENSION(:),   ALLOCATABLE :: ZAUX     ! used to distribute 
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZFACT    ! the total concentration in each shape
REAL,    DIMENSION(:),   ALLOCATABLE :: ZONEOVER_VAR ! for optimization
!
!
!-------------------------------------------------------------------------------
!
!                         #################
!                         FAST RS PROCESSES
!                         #################
!
IF (LSNOW) THEN
!
!
!*       1.1  Cloud droplet riming of the aggregates  
!        -------------------------------------------
!
!
  ZZW1(:,:) = 0.0
!
  GRIM(:) = (PRCT(:)>XRTMIN(2))        .AND. (PRST(:)>XRTMIN(5)) .AND. &
            (PRCS(:)>XRTMIN(2)/PTSTEP) .AND. (PZT(:)<XTT)
  IGRIM = COUNT( GRIM(:) )
!
  IF( IGRIM>0 ) THEN
!
!        1.1.0  allocations
!
    ALLOCATE(ZVEC1(IGRIM))
    ALLOCATE(ZVEC2(IGRIM))
    ALLOCATE(IVEC1(IGRIM))
    ALLOCATE(IVEC2(IGRIM))
!
!        1.1.1  select the PLBDAS
!
    ZVEC1(:) = PACK( PLBDAS(:),MASK=GRIM(:) )
!
!        1.1.2  find the next lower indice for the PLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete 
!               gamma function
!
    ZVEC2(1:IGRIM) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,           &
                         XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - FLOAT( IVEC2(1:IGRIM) )
!
!        1.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        1.1.4  riming of the small sized aggregates
!
    WHERE ( GRIM(:) )
      ZZW1(:,1) = MIN( PRCS(:),                         &
                       XCRIMSS * ZZW(:) * PRCT(:)       & ! RCRIMSS
                               *   PLBDAS(:)**XEXCRIMSS &
                               * PRHODREF(:)**(-XCEXVT) )
      PRCS(:) = PRCS(:) - ZZW1(:,1)
      PRSS(:) = PRSS(:) + ZZW1(:,1)
      PTHS(:) = PTHS(:) + ZZW1(:,1) * (PLSFACT(:) - PLVFACT(:)) ! f(L_f*(RCRIMSS))
!
      PCCS(:) = MAX( PCCS(:)-ZZW1(:,1)*(PCCT(:)/PRCT(:)),0.0 ) ! Lambda_c**3
    END WHERE
!
!        1.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        1.1.6  riming-conversion of the large sized aggregates into graupeln
!
!
    WHERE ( GRIM(:) .AND. (PRSS(:)>XRTMIN(5)/PTSTEP) )
      ZZW1(:,2) = MIN( PRCS(:),                         &
                       XCRIMSG * PRCT(:)                & ! RCRIMSG
                               *  PLBDAS(:)**XEXCRIMSG  &
                               * PRHODREF(:)**(-XCEXVT) &
               - ZZW1(:,1)              )
      ZZW1(:,3) = MIN( PRSS(:),                         &
                       XSRIMCG * PLBDAS(:)**XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(:) )/(PTSTEP*PRHODREF(:)))
      PRCS(:) = PRCS(:) - ZZW1(:,2)
      PRSS(:) = PRSS(:) - ZZW1(:,3)
      PRGS(:) = PRGS(:) + ZZW1(:,2) + ZZW1(:,3)
      PTHS(:) = PTHS(:) + ZZW1(:,2) * (PLSFACT(:) - PLVFACT(:)) ! f(L_f*(RCRIMSG))
!
      PCCS(:) = MAX( PCCS(:)-ZZW1(:,2)*(PCCT(:)/PRCT(:)),0.0 ) ! Lambda_c**3
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET (                                                  &
                   UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:), &
                                                                 4,'RIM_BU_RTH')
    IF (LBUDGET_RC) CALL BUDGET (                                                  &
                   UNPACK(PRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
                                                                 7,'RIM_BU_RRC')
    IF (LBUDGET_RS) CALL BUDGET (                                                  &
                   UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                                10,'RIM_BU_RRS')
    IF (LBUDGET_RG) CALL BUDGET (                                                  &
                   UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                                11,'RIM_BU_RRG')
    IF (LBUDGET_SV) THEN
      CALL BUDGET (UNPACK(PCCS(:),MASK=GMICRO(:,:,:),FIELD=PCCS_3D)*PRHODJ(:,:,:), &
                                                    12+NSV_LIMA_NC,'RIM_BU_RSV')
    END IF
  END IF
!
!
!*       1.2  Hallett-Mossop ice multiplication process due to snow riming  
!        -----------------------------------------------------------------
!
!
  GRIM(:) = (PZT(:)<XHMTMAX) .AND. (PZT(:)>XHMTMIN)                          &
                             .AND. (PRST(:)>XRTMIN(5)) .AND. (PRCT(:)>XRTMIN(2))
  IGRIM = COUNT( GRIM(:) )
  IF( IGRIM>0 ) THEN
    ALLOCATE(ZVEC1(IGRIM))
    ALLOCATE(ZVEC2(IGRIM))
    ALLOCATE(IVEC2(IGRIM))
!
    ZVEC1(:) = PACK( PLBDAC(:),MASK=GRIM(:) )
    ZVEC2(1:IGRIM) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,           &
                          XHMLINTP1 * LOG( ZVEC1(1:IGRIM) ) + XHMLINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - FLOAT( IVEC2(1:IGRIM) )
    ZVEC1(1:IGRIM) =   XGAMINC_HMC( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_HMC( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 ) ! Large droplets
!
    WHERE ( GRIM(:) .AND. ZZX(:)<0.99 )
      ZZW1(:,5) = (ZZW1(:,1)+ZZW1(:,2))*(PCCT(:)/PRCT(:))*(1.0-ZZX(:))* & 
                                                             XHM_FACTS* &
           MAX( 0.0, MIN( (PZT(:)-XHMTMIN)/3.0,(XHMTMAX-PZT(:))/2.0 ) ) ! CCHMSI
      PCIS(:) = PCIS(:) + ZZW1(:,5)
!
      ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMSI
      PRIS(:) = PRIS(:) + ZZW1(:,6)
      PRSS(:) = PRSS(:) - ZZW1(:,6)
    END WHERE
!
    IF (LCRYSTAL_SHAPE) THEN
! 2 shapes or more: secondary ice crystals are assumed to be columns (jsh=2)
      WHERE (GRIM(:) .AND. ZZX(:) < 0.99 .AND. ZZW1(:,5) .NE. 0.)
        PCIS_SHAPE(:,2) = PCIS_SHAPE(:,2) + ZZW1(:,5)
      END WHERE
    END IF
!
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_RI) CALL BUDGET (                                            &
                       UNPACK(PRIS(:),MASK=GMICRO,FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                                 9,'HMS_BU_RRI')
    IF (LBUDGET_RS) CALL BUDGET (                                            &
                       UNPACK(PRSS(:),MASK=GMICRO,FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                                10,'HMS_BU_RRS')
    IF (LBUDGET_SV)  THEN
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        CALL BUDGET  (UNPACK(PCIS(:),MASK=GMICRO,FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                 12+NSV_LIMA_NI,'HMS_BU_RSV')
      ELSE
        DO JSH = 1, NB_CRYSTAL_SHAPE
          CALL BUDGET  (UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO,FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                  12+NSV_LIMA_NI+ JSH-1,'HMS_BU_RSV')
        END DO
      END IF
    END IF
  END IF
!
!
!*      1.3  Ice multiplication process due to ice-ice collisions
!       ---------------------------------------------------------
!
  GCIBU(:) = LCIBU .AND. (PRST(:)>XRTMIN(5)) .AND. (PRGT(:)>XRTMIN(6))
  ICIBU    = COUNT( GCIBU(:) )
!
  IF (ICIBU > 0) THEN
!
!       1.3.0 randomization of XNDEBRIS_CIBU values
!
    IF (GFIRSTCALL) THEN
      CALL RANDOM_SEED(SIZE=NI_SEED) ! get size of seed
      ALLOCATE(I_SEED(NI_SEED))
      I_SEED(:) = I_SEED_PARAM !
      CALL RANDOM_SEED(PUT=I_SEED)
      GFIRSTCALL = .FALSE.
    END IF
!
    ALLOCATE(ZFRAGMENTS(ICIBU))
!
    IF (XNDEBRIS_CIBU >= 0.0) THEN
      ZFRAGMENTS(:) = XNDEBRIS_CIBU
    ELSE
!
! Mantissa gives the mean value (randomization around 10**MANTISSA)
! First digit after the comma provides the full range around 10**MANTISSA
!
      ALLOCATE(ZHARVEST(ICIBU))
!
      ZFACT1_XNDEBRIS = AINT(XNDEBRIS_CIBU)
      ZFACT2_XNDEBRIS = ABS(ANINT(10.0*(XNDEBRIS_CIBU - ZFACT1_XNDEBRIS)))
!
      CALL RANDOM_NUMBER(ZHARVEST(:))
!
      ZFRAGMENTS(:) = 10.0**(ZFACT2_XNDEBRIS*ZHARVEST(:) + ZFACT1_XNDEBRIS)
!
      DEALLOCATE(ZHARVEST)
!
! ZFRAGMENTS is a random variable containing the number of fragments per collision
! For XNDEBRIS_CIBU=-1.2345  => ZFRAGMENTS(:) = 10.0**(2.0*RANDOM_NUMBER(ZHARVEST(:)) - 1.0)
! and ZFRAGMENTS=[0.1, 10.0] centered around 1.0
!
    END IF
!
!       1.3.1 To compute the partial integration of snow gamma function
!
!       1.3.1.0 allocations

    ALLOCATE(ZVEC1_S(ICIBU))
    ALLOCATE(ZVEC1_S1(ICIBU))
    ALLOCATE(ZVEC1_S2(ICIBU))
    ALLOCATE(ZVEC1_S3(ICIBU))
    ALLOCATE(ZVEC1_S4(ICIBU))
    ALLOCATE(ZVEC1_S11(ICIBU))
    ALLOCATE(ZVEC1_S12(ICIBU))
    ALLOCATE(ZVEC1_S21(ICIBU))
    ALLOCATE(ZVEC1_S22(ICIBU))
    ALLOCATE(ZVEC1_S31(ICIBU))
    ALLOCATE(ZVEC1_S32(ICIBU))
    ALLOCATE(ZVEC1_S41(ICIBU))
    ALLOCATE(ZVEC1_S42(ICIBU))
    ALLOCATE(ZVEC2_S1(ICIBU))
    ALLOCATE(IVEC2_S1(ICIBU))
    ALLOCATE(ZVEC2_S2(ICIBU))
    ALLOCATE(IVEC2_S2(ICIBU))
!
!
!       1.3.1.1 select the PLBDAS
!
    ZVEC1_S(:) = PACK( PLBDAS(:),MASK=GCIBU(:) )

!
!       1.3.1.2 find the next lower indice for the PLBDAS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 1 (0.2 mm)
!
    ZVEC2_S1(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                        * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP1_S  ) )
    IVEC2_S1(1:ICIBU) = INT( ZVEC2_S1(1:ICIBU) )
    ZVEC2_S1(1:ICIBU) = ZVEC2_S1(1:ICIBU) - FLOAT( IVEC2_S1(1:ICIBU) )
!
!
!       1.3.1.3 find the next lower indice for the PLBDAS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 2 (1 mm)
!
    ZVEC2_S2(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                        * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP2_S  ) )
    IVEC2_S2(1:ICIBU) = INT( ZVEC2_S2(1:ICIBU) )
    ZVEC2_S2(1:ICIBU) = ZVEC2_S2(1:ICIBU) - FLOAT( IVEC2_S2(1:ICIBU) )
!
!
!       1.3.1.4 perform the linear interpolation of the
!               normalized "0"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
    ZVEC1_S11(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                       - XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
    ZVEC1_S12(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                       - XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! Computation of spectrum from 0.2 mm to 1 mm
    ZVEC1_S1(1:ICIBU) = ZVEC1_S12(1:ICIBU) - ZVEC1_S11(1:ICIBU)
!
!
!       1.3.1.5 perform the linear interpolation of the
!               normalized "XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
    ZVEC1_S21(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                       - XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
    ZVEC1_S22(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                       - XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
    ZVEC1_S2(1:ICIBU) = XMOMGS_CIBU_1 * (ZVEC1_S22(1:ICIBU) - ZVEC1_S21(1:ICIBU))
!
! For lower boundary (0.2 mm)
    ZVEC1_S31(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                       - XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
    ZVEC1_S32(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                       - XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
    ZVEC1_S3(1:ICIBU) = XMOMGS_CIBU_2 * (ZVEC1_S32(1:ICIBU) - ZVEC1_S31(1:ICIBU))
!
!
!       1.3.1.6 perform the linear interpolation of the
!               normalized "XBS+XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
    ZVEC1_S41(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                       - XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
    ZVEC1_S42(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                       - XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
    ZVEC1_S4(1:ICIBU) = XMOMGS_CIBU_3 * (ZVEC1_S42(1:ICIBU) - ZVEC1_S41(1:ICIBU))
!
    ALLOCATE(ZINTG_SNOW_1(SIZE(PZT)))
    ALLOCATE(ZINTG_SNOW_2(SIZE(PZT)))
    ALLOCATE(ZINTG_SNOW_3(SIZE(PZT)))
    ALLOCATE(ZINTG_SNOW_4(SIZE(PZT)))
!
    ZINTG_SNOW_1(:) = UNPACK ( VECTOR=ZVEC1_S1(:),MASK=GCIBU,FIELD=0.0 )
    ZINTG_SNOW_2(:) = UNPACK ( VECTOR=ZVEC1_S2(:),MASK=GCIBU,FIELD=0.0 )
    ZINTG_SNOW_3(:) = UNPACK ( VECTOR=ZVEC1_S3(:),MASK=GCIBU,FIELD=0.0 )
    ZINTG_SNOW_4(:) = UNPACK ( VECTOR=ZVEC1_S4(:),MASK=GCIBU,FIELD=0.0 )
!
!
!       1.3.2 Compute the partial integration of graupel gamma function
!
!       1.3.2.0 allocations
!
    ALLOCATE(ZVEC1_G(ICIBU))
    ALLOCATE(ZVEC1_G1(ICIBU))
    ALLOCATE(ZVEC1_G2(ICIBU))
    ALLOCATE(ZVEC2_G(ICIBU))
    ALLOCATE(IVEC2_G(ICIBU))
!
!
!       1.3.2.1 select the PLBDAG
!
    ZVEC1_G(:) = PACK( PLBDAG(:),MASK=GCIBU(:) )
!
!
!       1.3.2.2 find the next lower indice for the PLBDAG in the
!               geometrical set of Lbda_g used to tabulate some moments of the
!               incomplete gamma function, for the "2mm" boundary
!
    ZVEC2_G(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_G  &
                       * LOG( ZVEC1_G(1:ICIBU) ) + XCIBUINTP1_G  ) )
    IVEC2_G(1:ICIBU) = INT( ZVEC2_G(1:ICIBU) )
    ZVEC2_G(1:ICIBU) = ZVEC2_G(1:ICIBU) - FLOAT( IVEC2_G(1:ICIBU) )
!
!
!       1.3.2.3 perform the linear interpolation of the
!               normalized "2+XDG"-moment of the incomplete gamma function
!
    ZVEC1_G1(1:ICIBU) = XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                      - XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
    ZVEC1_G1(1:ICIBU) = XMOMGG_CIBU_1 * (1.0 - ZVEC1_G1(1:ICIBU))
!
!
!       1.3.2.4 perform the linear interpolation of the
!               normalized "2.0"-moment of the incomplete gamma function
!
    ZVEC1_G2(1:ICIBU) = XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                      - XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
    ZVEC1_G2(1:ICIBU) = XMOMGG_CIBU_2 * (1.0 - ZVEC1_G2(1:ICIBU))
!
!
    ALLOCATE(ZINTG_GRAUPEL_1(SIZE(PZT)))
    ALLOCATE(ZINTG_GRAUPEL_2(SIZE(PZT)))
!
    ZINTG_GRAUPEL_1(:) = UNPACK ( VECTOR=ZVEC1_G1(:),MASK=GCIBU,FIELD=0.0 )
    ZINTG_GRAUPEL_2(:) = UNPACK ( VECTOR=ZVEC1_G2(:),MASK=GCIBU,FIELD=0.0 )
!
!
!        1.3.3 To compute final "CIBU" contributions
!
    ALLOCATE(ZNI_CIBU(SIZE(PZT)))
    ALLOCATE(ZFRAG_CIBU(SIZE(PZT)))
!
    ZFRAG_CIBU(:) = UNPACK ( VECTOR=ZFRAGMENTS(:),MASK=GCIBU,FIELD=0.0 )
    ZNI_CIBU(:) = ZFRAG_CIBU(:) * (XFACTOR_CIBU_NI / (PRHODREF(:)**(XCEXVT-1.0))) * &
                  (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_1(:) *                     &
                   PLBDAS(:)**(XCXS) * PLBDAG(:)**(XCXG-(XDG+2.0))                  &
                 - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_2(:) *                     &
                   PLBDAS(:)**(XCXS-XDS) * PLBDAG(:)**(XCXG-2.0) )
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      PCIS(:) = PCIS(:) + ZNI_CIBU(:)
    ELSE
      IF (NB_CRYSTAL_SHAPE .EQ. 2) THEN
! only 2 shapes => the newly formed ice crystals are distributed in the 2 shapes
        PCIS_SHAPE(:,1) = PCIS_SHAPE(:,1) + ZNI_CIBU(:) * PCSHAPE_S(:,1)
        PCIS_SHAPE(:,2) = PCIS_SHAPE(:,2) + ZNI_CIBU(:) * PCSHAPE_S(:,2)
      ELSE IF (NB_CRYSTAL_SHAPE .GE. 3) THEN
! 3 shapes and more => the newly formed ice crystals are considered as irregular
        PCIS_SHAPE(:,3) = PCIS_SHAPE(:,3) + ZNI_CIBU(:)
      END IF
      PCIS(:) = PCIS(:) + ZNI_CIBU(:)
    END IF
!
    DEALLOCATE(ZFRAG_CIBU)
    DEALLOCATE(ZFRAGMENTS)
!
! Max value of rs removed by CIBU
    ALLOCATE(ZRI_CIBU(SIZE(PZT)))
    ZRI_CIBU(:) = (XFACTOR_CIBU_RI / (PRHODREF(:)**(XCEXVT+1.0))) *     &
                   (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_3(:) *        &
                    PLBDAS(:)**(XCXS-XBS) * PLBDAG(:)**(XCXG-(XDG+2.0)) &
                  - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_4(:) *        &
                    PLBDAS(:)**(XCXS-(XBS+XDS)) * PLBDAG(:)**(XCXG-2.0))
!
! The value of rs removed by CIBU is determined by the mean mass of pristine ice
    WHERE( PRIT(:)>XRTMIN(4) .AND. PCIT(:)>XCTMIN(4) )
      ZRI_CIBU(:) = MIN( ZRI_CIBU(:), PRSS(:), ZNI_CIBU(:)*PRIT(:)/PCIT(:) )
    ELSE WHERE
      ZRI_CIBU(:) = MIN( ZRI_CIBU(:), PRSS(:), MAX( ZNI_CIBU(:)*XMNU0,XRTMIN(4) ) )
    END WHERE
!
    PRIS(:) = PRIS(:) + ZRI_CIBU(:)   !
    PRSS(:) = PRSS(:) - ZRI_CIBU(:)   !
!
    DEALLOCATE(ZVEC1_S)
    DEALLOCATE(ZVEC1_S1)
    DEALLOCATE(ZVEC1_S2)
    DEALLOCATE(ZVEC1_S3)
    DEALLOCATE(ZVEC1_S4)
    DEALLOCATE(ZVEC1_S11)
    DEALLOCATE(ZVEC1_S12)
    DEALLOCATE(ZVEC1_S21)
    DEALLOCATE(ZVEC1_S22)
    DEALLOCATE(ZVEC1_S31)
    DEALLOCATE(ZVEC1_S32)
    DEALLOCATE(ZVEC1_S41)
    DEALLOCATE(ZVEC1_S42)
    DEALLOCATE(ZVEC2_S1)
    DEALLOCATE(IVEC2_S1)
    DEALLOCATE(ZVEC2_S2)
    DEALLOCATE(IVEC2_S2)
    DEALLOCATE(ZVEC1_G)
    DEALLOCATE(ZVEC1_G1)
    DEALLOCATE(ZVEC1_G2)
    DEALLOCATE(ZVEC2_G)
    DEALLOCATE(IVEC2_G)
    DEALLOCATE(ZINTG_SNOW_1)
    DEALLOCATE(ZINTG_SNOW_2)
    DEALLOCATE(ZINTG_SNOW_3)
    DEALLOCATE(ZINTG_SNOW_4)
    DEALLOCATE(ZINTG_GRAUPEL_1)
    DEALLOCATE(ZINTG_GRAUPEL_2)
    DEALLOCATE(ZNI_CIBU)
    DEALLOCATE(ZRI_CIBU)
  END IF
!
! Budget storage
!
  IF (LCIBU) THEN
    IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
      IF (LBUDGET_RI) CALL BUDGET (                                          &
                       UNPACK(PRIS(:),MASK=GMICRO,FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                                 9,'CIBU_BU_RRI')
      IF (LBUDGET_RS) CALL BUDGET (                                          &
                       UNPACK(PRSS(:),MASK=GMICRO,FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                                10,'CIBU_BU_RRS')
      IF (LBUDGET_SV) THEN
        IF (.NOT. LCRYSTAL_SHAPE) THEN
          CALL BUDGET  (UNPACK(PCIS(:),MASK=GMICRO,FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                   12+NSV_LIMA_NI,'CIBU_BU_RSV')
        ELSE
          DO JSH = 1, NB_CRYSTAL_SHAPE
            CALL BUDGET  (UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO,FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                    12+NSV_LIMA_NI+ JSH-1,'CIBU_BU_RSV')
          END DO
        END IF
      END IF
    END IF
  END IF
!
!
!*       1.3  Rain accretion onto the aggregates  
!        ---------------------------------------
!
!
  ZZW1(:,2:3) = 0.0
  GACC(:) = (PRRT(:)>XRTMIN(3)) .AND. (PRST(:)>XRTMIN(5)) .AND. (PRRS(:)>XRTMIN(3)/PTSTEP) .AND. (PZT(:)<XTT)
  IGACC = COUNT( GACC(:) )
!
  IF( IGACC>0 .AND. LRAIN) THEN
!
!        1.3.0  allocations
!
    ALLOCATE(ZVEC1(IGACC))
    ALLOCATE(ZVEC2(IGACC))
    ALLOCATE(ZVEC3(IGACC))
    ALLOCATE(IVEC1(IGACC))
    ALLOCATE(IVEC2(IGACC))
!
!        1.3.1  select the (PLBDAS,PLBDAR) couplet
!
    ZVEC1(:) = PACK( PLBDAS(:),MASK=GACC(:) )
    ZVEC2(:) = PACK( PLBDAR(:),MASK=GACC(:) )
!
!        1.3.2  find the next lower indice for the PLBDAS and for the PLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
    ZVEC1(1:IGACC) = MAX( 1.0001, MIN( FLOAT(NACCLBDAS)-0.0001,           &
                          XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - FLOAT( IVEC1(1:IGACC) )
!
    ZVEC2(1:IGACC) = MAX( 1.0001, MIN( FLOAT(NACCLBDAR)-0.0001,           &
                          XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - FLOAT( IVEC2(1:IGACC) )
!
!        1.3.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.3.4  raindrop accretion on the small sized aggregates
!
    WHERE ( GACC(:) )
      ZZW1(:,2) = PCRT(:) *                                           & !! coef of RRACCS
              XFRACCSS*( PLBDAS(:)**XCXS )*( PRHODREF(:)**(-XCEXVT-1.) ) &
         *( XLBRACCS1/((PLBDAS(:)**2)               ) +                  &
            XLBRACCS2/( PLBDAS(:)    * PLBDAR(:)    ) +                  &
            XLBRACCS3/(               (PLBDAR(:)**2)) )/PLBDAR(:)**3
      ZZW1(:,4) = MIN( PRRS(:),ZZW1(:,2)*ZZW(:) )           ! RRACCSS
      PRRS(:) = PRRS(:) - ZZW1(:,4)
      PRSS(:) = PRSS(:) + ZZW1(:,4)
      PTHS(:) = PTHS(:) + ZZW1(:,4)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RRACCSS))
!
      PCRS(:) = MAX( PCRS(:)-ZZW1(:,4)*(PCRT(:)/PRRT(:)),0.0 ) ! Lambda_r**3 
    END WHERE
!
!        1.3.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                         * ZVEC1(JJ) &
                 - (   XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) !! RRACCS
!
!        1.3.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * ZVEC2(JJ) &
                 - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.3.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
    WHERE ( GACC(:) .AND. (PRSS(:)>XRTMIN(5)/PTSTEP) )
      ZZW1(:,2) = MAX( MIN( PRRS(:),ZZW1(:,2)-ZZW1(:,4) ) , 0. )      ! RRACCSG
      ZZW1(:,3) = MIN( PRSS(:),XFSACCRG*ZZW(:)*                     & ! RSACCRG
            ( PLBDAS(:)**(XCXS-XBS) )*( PRHODREF(:)**(-XCEXVT-1.) ) &
           *( XLBSACCR1/((PLBDAR(:)**2)               ) +           &
              XLBSACCR2/( PLBDAR(:)    * PLBDAS(:)    ) +           &
              XLBSACCR3/(               (PLBDAS(:)**2)) ) )
      PRRS(:) = PRRS(:) - ZZW1(:,2)
      PRSS(:) = PRSS(:) - ZZW1(:,3)
      PRGS(:) = PRGS(:) + ZZW1(:,2)+ZZW1(:,3)
      PTHS(:) = PTHS(:) + ZZW1(:,2)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RRACCSG))
!
      PCRS(:) = MAX( PCRS(:)-ZZW1(:,2)*(PCRT(:)/PRRT(:)),0.0 ) ! Lambda_r**3 
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE .AND. LRAIN) THEN
    IF (LBUDGET_TH) CALL BUDGET (                                                 &
                   UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
                                                                 4,'ACC_BU_RTH')
    IF (LBUDGET_RR) CALL BUDGET (                                                  &
                   UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
                                                                 8,'ACC_BU_RRR')
    IF (LBUDGET_RS) CALL BUDGET (                                                  &
                   UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                                10,'ACC_BU_RRS')
    IF (LBUDGET_RG) CALL BUDGET (                                                  &
                   UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                                11,'ACC_BU_RRG')
    IF (LBUDGET_SV) THEN
      CALL BUDGET (UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
                                                    12+NSV_LIMA_NR,'ACC_BU_RSV')
    END IF
  END IF
!
!
!*       1.4  Conversion-Melting of the aggregates
!        -----------------------------------------
!
!
  ZZW(:) = 0.0
  WHERE( (PRST(:)>XRTMIN(5)) .AND. (PRSS(:)>XRTMIN(5)/PTSTEP) .AND. (PZT(:)>XTT) )
    ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
    ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
               ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                           *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RSMLT
!
    ZZW(:)  = MIN( PRSS(:), XFSCVMG*MAX( 0.0,( -ZZW(:) *             &
                           ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                             X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS ) -   &
                                     ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                              ( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                             ( PRHODREF(:)*XLMTT ) ) )
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
    PRSS(:) = PRSS(:) - ZZW(:)
    PRGS(:) = PRGS(:) + ZZW(:)
  END WHERE
!
! Budget storage
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_RS) CALL BUDGET (                                                     &
                         UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                                 10,'CMEL_BU_RRS')
    IF (LBUDGET_RG) CALL BUDGET (                                                     &
                         UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                                 11,'CMEL_BU_RRG')
  END IF
!
END IF ! LSNOW
!
!------------------------------------------------------------------------------
!
!                         #################
!                         FAST RG PROCESSES
!                         #################
!
!
!*       2.1  Rain contact freezing  
!        --------------------------
!
! Compute the proportion of each shape in the total concentration
! --> to be used to distribute the transfert rate in each shape
IF (LCRYSTAL_SHAPE) THEN
  ALLOCATE(ZFACT(SIZE(PZT),NB_CRYSTAL_SHAPE))
  ALLOCATE(ZAUX(SIZE(PZT)))
  ALLOCATE(ZONEOVER_VAR(SIZE(PZT)))
  ZFACT(:,:) = 0.
  ZONEOVER_VAR(:) = 0.
  WHERE (PCIS(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / PCIS(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE ((PCIS(:) .GT. 0.0) .AND. (PCIS_SHAPE(:,JSH) .GT. 0.0))
      ZFACT(:,JSH) = MIN(PCIS_SHAPE(:,JSH)*ZONEOVER_VAR(:), 1.0)
    END WHERE
  END DO
  ZAUX(:) = PCIS(:)
END IF
!
ZZW1(:,3:4) = 0.0
WHERE( (PRIT(:)>XRTMIN(4)) .AND. (PRRT(:)>XRTMIN(3)) .AND. &
       (PRIS(:)>XRTMIN(4)/PTSTEP) .AND. (PRRS(:)>XRTMIN(3)/PTSTEP) )
  ZZW1(:,3) = MIN( PRIS(:),XICFRR * PRIT(:) * PCRT(:)          & ! RICFRRG
                                  * PLBDAR(:)**XEXICFRR        &
                                  * PRHODREF(:)**(-XCEXVT-1.0) )
!
  ZZW1(:,4) = MIN( PRRS(:),XRCFRI * PCIT(:) * PCRT(:)          & ! RRCFRIG
                                 * PLBDAR(:)**XEXRCFRI        &
                                 * PRHODREF(:)**(-XCEXVT-2.0) )
!
  PRIS(:) = PRIS(:) - ZZW1(:,3)
  PRRS(:) = PRRS(:) - ZZW1(:,4)
  PRGS(:) = PRGS(:) + ZZW1(:,3) + ZZW1(:,4)
  PTHS(:) = PTHS(:) + ZZW1(:,4) * (PLSFACT(:) - PLVFACT(:)) ! f(L_f*RRCFRIG)
!
  PCRS(:)     = MAX( PCRS(:)-ZZW1(:,4)*(PCRT(:)/PRRT(:)),0.0 )     ! CRCFRIG
  PCIS(:) = MAX( PCIS(:)-ZZW1(:,3)*(PCIT(:)/PRIT(:)),0.0 )     ! CICFRRG
END WHERE
!
IF (LCRYSTAL_SHAPE) THEN
  ZAUX(:) = PCIS(:) - ZAUX(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE (ZAUX(:) .LT. 0.)
      PCIS_SHAPE(:,JSH) = PCIS_SHAPE(:,JSH) + ZAUX(:) * ZFACT(:,JSH)
    END WHERE
  END DO
END IF
!
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                                    &
                   UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:), &
                                                                4,'CFRZ_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET (                                                    &
                   UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
                                                                8,'CFRZ_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET (                                                    &
                   UNPACK(PRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                                9,'CFRZ_BU_RRI')
  IF (LBUDGET_RG) CALL BUDGET (                                                    &
                   UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                               11,'CFRZ_BU_RRG')
  IF (LBUDGET_SV) THEN
    CALL BUDGET ( UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
                                                  12+NSV_LIMA_NR,'CFRZ_BU_RSV')
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET ( UNPACK(PCIS(:),MASK=GMICRO(:,:,:),FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                  12+NSV_LIMA_NI,'CFRZ_BU_RSV')
    ELSE
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET ( UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO(:,:,:),FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                  12+NSV_LIMA_NI+JSH-1,'CFRZ_BU_RSV')
      END DO
    END IF
  END IF  ! LBUDGET_SV
END IF
!
!
!*       2.2  Ice multiplication process following rain contact freezing  
!        ---------------------------------------------------------------
!
GRDSF(:) = LRDSF .AND. (PRIT(:)>0.0) .AND. (PRRT(:)>0.0) .AND. &
                       (PRIS(:)>0.0) .AND. (PRRS(:)>0.0)
IRDSF    = COUNT( GRDSF(:) )
!
IF( IRDSF>0 ) THEN
!
  ALLOCATE(ZVEC1_R(IRDSF)) 
  ALLOCATE(ZVEC1_R1(IRDSF)) 
  ALLOCATE(ZVEC2_R(IRDSF)) 
  ALLOCATE(IVEC2_R(IRDSF))
!
!*       2.2.1  select the ZLBDAR
!
  ZVEC1_R(:) = PACK( PLBDAR(:),MASK=GRDSF(:) )

!*       2.2.2  find the next lower indice for the ZLBDAR in the
!               geometrical set of Lbda_r used to tabulate some moments of the
!               incomplete gamma function, for the lower boundary (0.1 mm)
!
  ZVEC2_R(1:IRDSF) = MAX( 1.00001, MIN( FLOAT(NGAMINC)-0.00001,XRDSFINTP_R  &
                      * LOG( ZVEC1_R(1:IRDSF) ) + XRDSFINTP1_R  ) )             
  IVEC2_R(1:IRDSF) = INT( ZVEC2_R(1:IRDSF) )
  ZVEC2_R(1:IRDSF) = ZVEC2_R(1:IRDSF) - FLOAT( IVEC2_R(1:IRDSF) )
!
!*       2.2.3  perform the linear interpolation of the 
!               normalized "2+XDR"-moment of the incomplete gamma function
!
  ZVEC1_R1(1:IRDSF) = XGAMINC_RDSF_R(IVEC2_R(1:IRDSF)+1)*ZVEC2_R(1:IRDSF)    &
                    - XGAMINC_RDSF_R(IVEC2_R(1:IRDSF))* (ZVEC2_R(1:IRDSF)-1.0)
!
!  From 0.1 mm to infinity we need  
  ZVEC1_R1(1:IRDSF) = XMOMGR_RDSF * (1.0 - ZVEC1_R1(1:IRDSF))
!
  ALLOCATE(ZINTG_RAIN(SIZE(PZT)))
  ZINTG_RAIN(:) = UNPACK ( VECTOR=ZVEC1_R1(:),MASK=GRDSF,FIELD=0.0 )
!
!*       2.2.4  To compute final "RDSF" contributions
!
  ALLOCATE(ZNI_RDSF(SIZE(PZT)))
  ZNI_RDSF(:) = (XFACTOR_RDSF_NI/(PRHODREF(:)**(XCEXVT-1.0)))*(  &
                PCIT(:)*PCRT(:)*ZINTG_RAIN(:)*PLBDAR(:)**(-(XDR+6.0)) )
!
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    PCIS(:) = PCIS(:) + ZNI_RDSF(:)
  ELSE
    IF (NB_CRYSTAL_SHAPE .EQ. 1) THEN
! only 1 shape => we put all the newly formed ice crystals in this shape
      PCIS_SHAPE(:,1) = PCIS_SHAPE(:,1) + ZNI_RDSF(:)
    ELSE IF (NB_CRYSTAL_SHAPE .EQ. 2) THEN
! only 2 shapes => the newly formed ice crystals are equally distributed in the 2 shapes
      PCIS_SHAPE(:,1) = PCIS_SHAPE(:,1) + ZNI_RDSF(:) * PCSHAPE_S(:,1)
      PCIS_SHAPE(:,2) = PCIS_SHAPE(:,2) + ZNI_RDSF(:) * PCSHAPE_S(:,2)
    ELSE IF (NB_CRYSTAL_SHAPE .GE. 3) THEN
! 3 shapes and more => the newly formed ice crystals are considered as irregular
      PCIS_SHAPE(:,3) = PCIS_SHAPE(:,3) + ZNI_RDSF(:)
    END IF
    PCIS(:) = PCIS(:) + ZNI_RDSF(:)
  END IF
!
! The value of rg removed by RDSF is determined by the mean mass of pristine ice
  ALLOCATE(ZRI_RDSF(SIZE(PZT)))
  ZRI_RDSF(:) = MIN( PRGS(:), MAX( ZNI_RDSF(:)*XMNU0,XRTMIN(5) ) )
!
  PRIS(:) = PRIS(:) + ZRI_RDSF(:) 
  PRGS(:) = PRGS(:) - ZRI_RDSF(:)
!
  DEALLOCATE(ZINTG_RAIN)
  DEALLOCATE(ZVEC1_R) 
  DEALLOCATE(ZVEC1_R1) 
  DEALLOCATE(ZVEC2_R) 
  DEALLOCATE(IVEC2_R)
  DEALLOCATE(ZNI_RDSF)
  DEALLOCATE(ZRI_RDSF)
ENDIF
!
! Budget storage
!
IF (LRDSF) THEN
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_RI) CALL BUDGET (                                          &
                     UNPACK(PRIS(:),MASK=GMICRO,FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                               9,'RDSF_BU_RRI')
    IF (LBUDGET_RG) CALL BUDGET (                                          &
                     UNPACK(PRGS(:),MASK=GMICRO,FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                              11,'RDSF_BU_RRG')
    IF (LBUDGET_SV) THEN
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        CALL BUDGET  (UNPACK(PCIS(:),MASK=GMICRO,FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                 12+NSV_LIMA_NI,'RDSF_BU_RSV')
      ELSE
        DO JSH = 1, NB_CRYSTAL_SHAPE
          CALL BUDGET  (UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO,FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                  12+NSV_LIMA_NI+ JSH-1,'RDSF_BU_RSV')
        END DO
      END IF
    END IF
  END IF
END IF
!
!
!*       2.3  Compute the Dry growth case
!        --------------------------------
!
!
ZZW1(:,:) = 0.0
WHERE( ((PRCT(:)>XRTMIN(2)) .AND. (PRGT(:)>XRTMIN(6)) .AND. (PRCS(:)>XRTMIN(2)/PTSTEP)) .OR. &
       ((PRIT(:)>XRTMIN(4)) .AND. (PRGT(:)>XRTMIN(6)) .AND. (PRIS(:)>XRTMIN(4)/PTSTEP))      )
  ZZW(:) = PLBDAG(:)**(XCXG-XDG-2.0) * PRHODREF(:)**(-XCEXVT)
  ZZW1(:,1) = MIN( PRCS(:),XFCDRYG * PRCT(:) * ZZW(:) )             ! RCDRYG
  ZZW1(:,2) = MIN( PRIS(:),XFIDRYG * EXP( XCOLEXIG*(PZT(:)-XTT) ) &
                                   * PRIT(:) * ZZW(:) )             ! RIDRYG
END WHERE
!
!*       2.3.1  accretion of aggregates on the graupeln
!        ----------------------------------------------
!
GDRY(:) = (PRST(:)>XRTMIN(5)) .AND. (PRGT(:)>XRTMIN(6)) .AND. (PRSS(:)>XRTMIN(5)/PTSTEP)
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.3.2  allocations
!
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(ZVEC3(IGDRY))
  ALLOCATE(IVEC1(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
!*       2.3.3  select the (PLBDAG,PLBDAS) couplet
!
  ZVEC1(:) = PACK( PLBDAG(:),MASK=GDRY(:) )
  ZVEC2(:) = PACK( PLBDAS(:),MASK=GDRY(:) )
!
!*       2.3.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
  ZVEC1(1:IGDRY) = MAX( 1.0001, MIN( FLOAT(NDRYLBDAG)-0.0001,           &
                        XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
  IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
  ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )
!
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( FLOAT(NDRYLBDAS)-0.0001,           &
                        XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2S ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       2.3.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
  DO JJ = 1,IGDRY
    ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                       * ZVEC1(JJ) &
               - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                       * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:) )
    ZZW1(:,3) = MIN( PRSS(:),XFSDRYG*ZZW(:)                         & ! RSDRYG
                                    * EXP( XCOLEXSG*(PZT(:)-XTT) )  &
                  *( PLBDAS(:)**(XCXS-XBS) )*( PLBDAG(:)**XCXG )    &
                  *( PRHODREF(:)**(-XCEXVT-1.) )                    &
                       *( XLBSDRYG1/( PLBDAG(:)**2              ) + &
                          XLBSDRYG2/( PLBDAG(:)   * PLBDAS(:)   ) + &
                          XLBSDRYG3/(               PLBDAS(:)**2) ) )
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC3)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
!*       2.3.6  accretion of raindrops on the graupeln
!        ---------------------------------------------
!
GDRY(:) = (PRRT(:)>XRTMIN(3)) .AND. (PRGT(:)>XRTMIN(6)) .AND. (PRRS(:)>XRTMIN(3))
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.3.7  allocations
!
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(ZVEC3(IGDRY))
  ALLOCATE(IVEC1(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
!*       2.3.8  select the (PLBDAG,PLBDAR) couplet
!
  ZVEC1(:) = PACK( PLBDAG(:),MASK=GDRY(:) )
  ZVEC2(:) = PACK( PLBDAR(:),MASK=GDRY(:) )
!
!*       2.3.9  find the next lower indice for the PLBDAG and for the PLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
  ZVEC1(1:IGDRY) = MAX( 1.0001, MIN( FLOAT(NDRYLBDAG)-0.0001,           &
                        XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
  IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
  ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )
!
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( FLOAT(NDRYLBDAR)-0.0001,           &
                        XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2R ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       2.3.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
  DO JJ = 1,IGDRY
    ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                       * ZVEC1(JJ) &
               - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                       * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:) )
    ZZW1(:,4) = MIN( PRRS(:),XFRDRYG*ZZW(:) * PCRT(:)                   & ! RRDRYG
                      *( PLBDAR(:)**(-3) )*( PLBDAG(:)**XCXG ) &
                              *( PRHODREF(:)**(-XCEXVT-1.) )   &
                  *( XLBRDRYG1/( PLBDAG(:)**2              ) + &
                     XLBRDRYG2/( PLBDAG(:)   * PLBDAR(:)   ) + &
                     XLBRDRYG3/(               PLBDAR(:)**2) ) )
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC3)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!
!*       2.4  Compute the Wet growth case
!        --------------------------------
!
!
ZZW(:) = 0.0
ZRWETG(:) = 0.0
WHERE( PRGT(:)>XRTMIN(6) )
  ZZW1(:,5) = MIN( PRIS(:),                                    &
              ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(PZT(:)-XTT)) ) ) ! RIWETG
  ZZW1(:,6) = MIN( PRSS(:),                                    &
              ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(PZT(:)-XTT)) ) ) ! RSWETG
!
  ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
  ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                  &
            ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT ))   &
                          *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RWETG
!
  ZRWETG(:)  = MAX( 0.0,                                               &
                  ( ZZW(:) * ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                               X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) +   &
                  ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                  ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PZT(:)))   ) ) / &
                                  ( PRHODREF(:)*(XLMTT-XCL*(XTT-PZT(:))) )   )
END WHERE
!
!
!*       2.5  Select Wet or Dry case
!        ---------------------------
!
!
! Wet case and partial conversion to hail
!
! Compute the proportion of each shape in the total concentration
! --> to be used to distribute the transfert rate in each shape
IF (LCRYSTAL_SHAPE) THEN
  ZFACT(:,:) = 0.
  ZONEOVER_VAR(:) = 0.
  WHERE (PCIS(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / PCIS(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE ((PCIS(:) .GT. 0.0) .AND. (PCIS_SHAPE(:,JSH) .GT. 0.0))
      ZFACT(:,JSH) = MIN(PCIS_SHAPE(:,JSH)*ZONEOVER_VAR(:), 1.0)
    END WHERE
  END DO
  ZAUX(:) = PCIS(:)
END IF
!
ZZW(:) = 0.0
NHAIL = 0.
IF (LHAIL) NHAIL = 1. 
WHERE( PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT                               &
                         .AND. ZRDRYG(:)>=ZRWETG(:) .AND. ZRWETG(:)>0.0 ) 
!   
  ZZW(:) = ZRWETG(:) - ZZW1(:,5) - ZZW1(:,6) ! RCWETG+RRWETG
!   
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!   
  ZZW1(:,7) = MAX( 0.0,MIN( ZZW(:),PRRS(:)+ZZW1(:,1) ) )
  ZZX(:)    = ZZW1(:,7) / ZZW(:)
  ZZW1(:,5) = ZZW1(:,5)*ZZX(:)
  ZZW1(:,6) = ZZW1(:,6)*ZZX(:)
  ZRWETG(:) = ZZW1(:,7) + ZZW1(:,5) + ZZW1(:,6)
!   
  PRCS(:) = PRCS(:) - ZZW1(:,1)
  PRIS(:) = PRIS(:) - ZZW1(:,5)
  PRSS(:) = PRSS(:) - ZZW1(:,6)
!
! assume a linear percent of conversion of graupel into hail
!
  PRGS(:) = PRGS(:) + ZRWETG(:)
  ZZW(:)  = PRGS(:) * ZRDRYG(:) * NHAIL / (ZRWETG(:) + ZRDRYG(:)) 
  PRGS(:) = PRGS(:) - ZZW(:)                        
  PRHS(:) = PRHS(:) + ZZW(:)
  PRRS(:) = MAX( 0.0,PRRS(:) - ZZW1(:,7) + ZZW1(:,1) )
  PTHS(:) = PTHS(:) + ZZW1(:,7) * (PLSFACT(:) - PLVFACT(:))
                                                ! f(L_f*(RCWETG+RRWETG))
!
  PCCS(:) = MAX( PCCS(:)-ZZW1(:,1)*(PCCT(:)/MAX(PRCT(:),XRTMIN(2))),0.0 )
  PCIS(:) = MAX( PCIS(:)-ZZW1(:,5)*(PCIT(:)/MAX(PRIT(:),XRTMIN(4))),0.0 )
  PCRS(:) = MAX( PCRS(:)-MAX( ZZW1(:,7)-ZZW1(:,1),0.0 )                 &
                        *(PCRT(:)/MAX(PRRT(:),XRTMIN(3))),0.0 )
END WHERE
!
IF (LCRYSTAL_SHAPE) THEN
  ZAUX(:) = PCIS(:) - ZAUX(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE (ZAUX(:) .LT. 0.)
      PCIS_SHAPE(:,JSH) = PCIS_SHAPE(:,JSH) + ZAUX(:) * ZFACT(:,JSH)
    END WHERE
  END DO
END IF
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                                 &
                 UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
                                                             4,'WETG_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET (                                                  &
                 UNPACK(PRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
                                                            7,'WETG_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET (                                                 &
                UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
                                                            8,'WETG_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET (                                                  &
                 UNPACK(PRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                             9,'WETG_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET (                                                  &
                 UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                            10,'WETG_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET (                                                  &
                 UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                            11,'WETG_BU_RRG')
  IF (LBUDGET_RH) CALL BUDGET (                                                  &
                 UNPACK(PRHS(:),MASK=GMICRO(:,:,:),FIELD=PRHS_3D)*PRHODJ(:,:,:), &
                                                            12,'WETG_BU_RRH')
  IF (LBUDGET_SV) THEN
    CALL BUDGET (UNPACK(PCCS(:),MASK=GMICRO(:,:,:),FIELD=PCCS_3D)*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NC,'WETG_BU_RSV')
    CALL BUDGET (UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NR,'WETG_BU_RSV')
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET (UNPACK(PCIS(:),MASK=GMICRO(:,:,:),FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NI,'WETG_BU_RSV')
    ELSE
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET (UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO(:,:,:),FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                            12+NSV_LIMA_NI+JSH-1,'WETG_BU_RSV')
      END DO
    END IF
  END IF  ! LBUDGET_SV
END IF
!
! Dry case
!
! Compute the proportion of each shape in the total concentration
! --> to be used to distribute the transfert rate in each shape
IF (LCRYSTAL_SHAPE) THEN
  ZFACT(:,:) = 0.
  ZONEOVER_VAR(:) = 0.
  WHERE (PCIS(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / PCIS(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE ((PCIS(:) .GT. 0.0) .AND. (PCIS_SHAPE(:,JSH) .GT. 0.0))
      ZFACT(:,JSH) = MIN(PCIS_SHAPE(:,JSH)*ZONEOVER_VAR(:), 1.0)
    END WHERE
  END DO
  ZAUX(:) = PCIS(:)
END IF
!
WHERE( PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT                              &
                         .AND. ZRDRYG(:)<ZRWETG(:) .AND. ZRDRYG(:)>0.0 ) ! case
  PRCS(:) = PRCS(:) - ZZW1(:,1)
  PRIS(:) = PRIS(:) - ZZW1(:,2)
  PRSS(:) = PRSS(:) - ZZW1(:,3)
  PRRS(:) = PRRS(:) - ZZW1(:,4)
  PRGS(:) = PRGS(:) + ZRDRYG(:)
  PTHS(:) = PTHS(:) + (ZZW1(:,1)+ZZW1(:,4))*(PLSFACT(:)-PLVFACT(:)) !
                                                                    ! f(L_f*(RCDRYG+RRDRYG))
!
  PCCS(:) = MAX( PCCS(:)-ZZW1(:,1)*(PCCT(:)/MAX(PRCT(:),XRTMIN(2))),0.0 )
  PCIS(:) = MAX( PCIS(:)-ZZW1(:,2)*(PCIT(:)/MAX(PRIT(:),XRTMIN(4))),0.0 )
  PCRS(:) = MAX( PCRS(:)-ZZW1(:,4)*(PCRT(:)/MAX(PRRT(:),XRTMIN(3))),0.0 ) 
                                                        ! Approximate rates
END WHERE
!
IF (LCRYSTAL_SHAPE) THEN
  ZAUX(:) = PCIS(:) - ZAUX(:)
  DO JSH = 1, NB_CRYSTAL_SHAPE
    WHERE (ZAUX(:) .LT. 0.)
      PCIS_SHAPE(:,JSH) = PCIS_SHAPE(:,JSH) + ZAUX(:) * ZFACT(:,JSH)
    END WHERE
  END DO
END IF
!
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                                   &
                   UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
                                                                4,'DRYG_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET (                                                    &
                   UNPACK(PRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
                                                                7,'DRYG_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET (                                                    &
                   UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
                                                                8,'DRYG_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET (                                                    &
                   UNPACK(PRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                                9,'DRYG_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET (                                                    &
                   UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
                                                               10,'DRYG_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET (                                                    &
                   UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                               11,'DRYG_BU_RRG')
  IF (LBUDGET_SV) THEN
    CALL BUDGET (  UNPACK(PCCS(:),MASK=GMICRO(:,:,:),FIELD=PCCS_3D)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NC,'DRYG_BU_RSV')
    CALL BUDGET (  UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NR,'DRYG_BU_RSV')
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET (  UNPACK(PCIS(:),MASK=GMICRO(:,:,:),FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NI,'DRYG_BU_RSV')
    ELSE
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET (  UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO(:,:,:),FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NI+JSH-1,'DRYG_BU_RSV')
      END DO
    END IF
  END IF  ! LBUDGET_SV
END IF
!
!
!*       2.6  Hallett-Mossop ice multiplication process due to graupel riming
!        --------------------------------------------------------------------
!
!
GDRY(:) = (PZT(:)<XHMTMAX) .AND. (PZT(:)>XHMTMIN)    .AND. (ZRDRYG(:)<ZRWETG(:))&
                           .AND. (PRGT(:)>XRTMIN(6)) .AND. (PRCT(:)>XRTMIN(2))
IGDRY = COUNT( GDRY(:) )
IF( IGDRY>0 ) THEN
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
  ZVEC1(:) = PACK( PLBDAC(:),MASK=GDRY(:) )
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,           &
                        XHMLINTP1 * LOG( ZVEC1(1:IGDRY) ) + XHMLINTP2 ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
  ZVEC1(1:IGDRY) =   XGAMINC_HMC( IVEC2(1:IGDRY)+1 )* ZVEC2(1:IGDRY)      &
                   - XGAMINC_HMC( IVEC2(1:IGDRY)   )*(ZVEC2(1:IGDRY) - 1.0)
  ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GDRY,FIELD=0.0 ) ! Large droplets
!
  WHERE ( GDRY(:) .AND. ZZX(:) < 0.99 ) ! Dry case
    ZZW1(:,5) = ZZW1(:,1)*(PCCT(:)/PRCT(:))*(1.0-ZZX(:))*XHM_FACTG*  &
         MAX( 0.0, MIN( (PZT(:)-XHMTMIN)/3.0,(XHMTMAX-PZT(:))/2.0 ) ) ! CCHMGI
    PCIS(:) = PCIS(:) + ZZW1(:,5)
    ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMGI
    PRIS(:) = PRIS(:) + ZZW1(:,6)
    PRGS(:) = PRGS(:) - ZZW1(:,6)
  END WHERE
! on traite HMG comme HMS --> regime de croissance colonnes => production de colonnes
  IF (LCRYSTAL_SHAPE .AND. NB_CRYSTAL_SHAPE .GE. 2) THEN
! 2 shapes or more: secondary ice crystals are assumed to be columns (jsh=2)
    WHERE (GDRY(:) .AND. ZZX(:) < 0.99 .AND. ZZW1(:,5) .NE. 0.)
      PCIS_SHAPE(:,2) = PCIS_SHAPE(:,2) + ZZW1(:,5)
    END WHERE
  END IF
!
  DEALLOCATE(IVEC2)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
! Budget storage
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_RI) CALL BUDGET (                                               &
                     UNPACK(PRIS(:),MASK=GMICRO,FIELD=PRIS_3D)*PRHODJ(:,:,:), &
                                                               9,'HMG_BU_RRI')
  IF (LBUDGET_RG) CALL BUDGET (                                               &
                     UNPACK(PRGS(:),MASK=GMICRO,FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                              11,'HMG_BU_RRG')
  IF (LBUDGET_SV) THEN
    IF (.NOT. LCRYSTAL_SHAPE) THEN
      CALL BUDGET ( UNPACK(PCIS(:),MASK=GMICRO,FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                                             12+NSV_LIMA_NI,'HMG_BU_RSV')
    ELSE
      DO JSH = 1, NB_CRYSTAL_SHAPE
        CALL BUDGET ( UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO,FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                                             12+NSV_LIMA_NI+JSH-1,'HMG_BU_RSV')
      END DO
    END IF
  END IF   !LBUDGET_SV
END IF
!
!
!*       2.7  Melting of the graupeln
!        ----------------------------
!
!
ZZW(:) = 0.0
WHERE( (PRGT(:)>XRTMIN(6)) .AND. (PRGS(:)>XRTMIN(6)/PTSTEP) .AND. (PZT(:)>XTT) )
  ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
  ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
             ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                         *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RGMLTR
!
  ZZW(:)  = MIN( PRGS(:), MAX( 0.0,( -ZZW(:) *                     &
                         ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                           X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) -   &
                                   ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                            ( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                           ( PRHODREF(:)*XLMTT ) ) )
  PRRS(:) = PRRS(:) + ZZW(:)
  PRGS(:) = PRGS(:) - ZZW(:)
  PTHS(:) = PTHS(:) - ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RGMLTR))
!
!   PCRS(:) = MAX( PCRS(:) + ZZW(:)*(XCCG*PLBDAG(:)**XCXG/PRGT(:)),0.0 )
  PCRS(:) = PCRS(:) + ZZW(:)*5.0E6  ! obtained after averaging
                                    ! Dshed=1mm and 500 microns
END WHERE
!
IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET (                                                   &
                   UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
                                                                4,'GMLT_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET (                                                    &
                   UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
                                                                8,'GMLT_BU_RRR')
  IF (LBUDGET_RG) CALL BUDGET (                                                    &
                   UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
                                                               11,'GMLT_BU_RRG')
  IF (LBUDGET_SV) THEN
    CALL BUDGET (  UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NR,'GMLT_BU_RSV')
  END IF
END IF
!
!
!------------------------------------------------------------------------------
!
!                         #################
!                         FAST RH PROCESSES
!                         #################
!
!
IF (LHAIL) THEN
!
  GHAIL(:) = PRHT(:)>XRTMIN(7)
  IHAIL = COUNT(GHAIL(:))
!
  IF( IHAIL>0 ) THEN
!
!*       3.1 Wet growth of hail 
!        ----------------------------
!
    ZZW1(:,:) = 0.0
    WHERE( GHAIL(:) .AND. ( (PRCT(:)>XRTMIN(2) .AND. PRCS(:)>XRTMIN(2)/PTSTEP) .OR. &
                            (PRIT(:)>XRTMIN(4) .AND. PRIS(:)>XRTMIN(4)/PTSTEP) )    )    
      ZZW(:) = PLBDAH(:)**(XCXH-XDH-2.0) * PRHODREF(:)**(-XCEXVT)
      ZZW1(:,1) = MIN( PRCS(:),XFWETH * PRCT(:) * ZZW(:) )             ! RCWETH
      ZZW1(:,2) = MIN( PRIS(:),XFWETH * PRIT(:) * ZZW(:) )             ! RIWETH
    END WHERE
!
!*       3.1.1  accretion of aggregates on the hailstones
!        ------------------------------------------------
!
    GWET(:) = GHAIL(:) .AND. (PRST(:)>XRTMIN(5) .AND. PRSS(:)>XRTMIN(5)/PTSTEP)
    IGWET = COUNT( GWET(:) )
!
    IF( IGWET>0 ) THEN
!
!*       3.1.2  allocations
!
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       3.1.3  select the (PLBDAH,PLBDAS) couplet
!
      ZVEC1(:) = PACK( PLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( PLBDAS(:),MASK=GWET(:) )
!
!*       3.1.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.0001, MIN( FLOAT(NWETLBDAH)-0.0001,           &
                            XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.0001, MIN( FLOAT(NWETLBDAS)-0.0001,           &
                            XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + XWETINTP2S ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       3.1.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                   - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,3) = MIN( PRSS(:),XFSWETH*ZZW(:)                       & ! RSWETH
                      *( PLBDAS(:)**(XCXS-XBS) )*( PLBDAH(:)**XCXH )  &
                         *( PRHODREF(:)**(-XCEXVT-1.) )               &
                         *( XLBSWETH1/( PLBDAH(:)**2              ) + &
                            XLBSWETH2/( PLBDAH(:)   * PLBDAS(:)   ) + &
                            XLBSWETH3/(               PLBDAS(:)**2) ) )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
!
!*       3.1.6  accretion of graupeln on the hailstones
!        ----------------------------------------------
!
    GWET(:) = GHAIL(:) .AND. (PRGT(:)>XRTMIN(6) .AND. PRGS(:)>XRTMIN(6)/PTSTEP)
    IGWET = COUNT( GWET(:) )
!
    IF( IGWET>0 ) THEN
!
!*       3.1.7  allocations
!
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       3.1.8  select the (PLBDAH,PLBDAG) couplet
!
      ZVEC1(:) = PACK( PLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( PLBDAG(:),MASK=GWET(:) )
!
!*       3.1.9  find the next lower indice for the PLBDAH and for the PLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.0001, MIN( FLOAT(NWETLBDAG)-0.0001,           &
                            XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.0001, MIN( FLOAT(NWETLBDAG)-0.0001,           &
                            XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + XWETINTP2G ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       3.1.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                  - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,5) = MAX(MIN( PRGS(:),XFGWETH*ZZW(:)                       & ! RGWETH
                     *( PLBDAG(:)**(XCXG-XBG) )*( PLBDAH(:)**XCXH )  &
                         *( PRHODREF(:)**(-XCEXVT-1.) )               &
                         *( XLBGWETH1/( PLBDAH(:)**2              ) + &
                            XLBGWETH2/( PLBDAH(:)   * PLBDAG(:)   ) + &
                            XLBGWETH3/(               PLBDAG(:)**2) ) ),0. )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF ! igwet
!
!*       3.2    compute the Wet growth of hail
!        -------------------------------------
!
    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. PZT(:)<XTT )
      ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
      ZZW(:) = PKA(:)*(XTT-PZT(:)) +                                 &
                ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                            *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RWETH
!
      ZZW(:)  =  MAX(0.,  ( ZZW(:) * ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                                X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) +   &
                   ( ZZW1(:,2)+ZZW1(:,3)+ZZW1(:,5) ) *                  &
                   ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PZT(:)))   ) ) / &
                         ( PRHODREF(:)*(XLMTT-XCL*(XTT-PZT(:))) ) )
!
      ZZW1(:,6) = MAX( ZZW(:) - ZZW1(:,2) - ZZW1(:,3) - ZZW1(:,5),0.) ! RCWETH+RRWETH
    END WHERE
!
! Compute the proportion of each shape in the total concentration
! --> to be used to distribute the transfert rate in each shape
    IF (LCRYSTAL_SHAPE) THEN
      ZFACT(:,:) = 0.
      WHERE (PCIS(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / PCIS(:)
      DO JSH = 1, NB_CRYSTAL_SHAPE
        WHERE ((PCIS(:) .GT. 0.0) .AND. (PCIS_SHAPE(:,JSH) .GT. 0.0))
          ZFACT(:,JSH) = MIN(PCIS_SHAPE(:,JSH)*ZONEOVER_VAR(:), 1.0)
        END WHERE
      END DO
      ZAUX(:) = PCIS(:)
    END IF
!
    WHERE ( GHAIL(:) .AND. PZT(:)<XTT  .AND. ZZW1(:,6)/=0.)
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
      ZZW1(:,4) = MAX( 0.0,MIN( ZZW1(:,6),PRRS(:)+ZZW1(:,1) ) )
      ZZX(:)    = ZZW1(:,4) / ZZW1(:,6)
      ZZW1(:,2) = ZZW1(:,2)*ZZX(:)
      ZZW1(:,3) = ZZW1(:,3)*ZZX(:)
      ZZW1(:,5) = ZZW1(:,5)*ZZX(:)
      ZZW(:)    = ZZW1(:,4) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5)
!
!*       3.2.1  integrate the Wet growth of hail
!
      PRCS(:) = PRCS(:) - ZZW1(:,1)
      PRIS(:) = PRIS(:) - ZZW1(:,2)
      PRSS(:) = PRSS(:) - ZZW1(:,3)
      PRGS(:) = PRGS(:) - ZZW1(:,5)
      PRHS(:) = PRHS(:) + ZZW(:)
      PRRS(:) = MAX( 0.0,PRRS(:) - ZZW1(:,4) + ZZW1(:,1) )
      PTHS(:) = PTHS(:) + ZZW1(:,4) * (PLSFACT(:) - PLVFACT(:)) 
                                  ! f(L_f*(RCWETH+RRWETH))
!
      PCCS(:) = MAX( PCCS(:)-ZZW1(:,1)*(PCCT(:)/MAX(PRCT(:),XRTMIN(2))),0.0 )
      PCIS(:) = MAX( PCIS(:)-ZZW1(:,2)*(PCIT(:)/MAX(PRIT(:),XRTMIN(4))),0.0 )
      PCRS(:) = MAX( PCRS(:)-MAX( ZZW1(:,4)-ZZW1(:,1),0.0 )                 &
                                        *(PCRT(:)/MAX(PRRT(:),XRTMIN(3))),0.0 )
    END WHERE
!
    IF (LCRYSTAL_SHAPE) THEN
      ZAUX(:) = PCIS(:) - ZAUX(:)
      DO JSH = 1, NB_CRYSTAL_SHAPE
        WHERE (ZAUX(:) .LT. 0.)
          PCIS_SHAPE(:,JSH) = PCIS_SHAPE(:,JSH) + ZAUX(:) * ZFACT(:,JSH)
        END WHERE
      END DO
    END IF
  END IF ! IHAIL>0
!
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET (                                        &
         UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:), &
         4,'WETH_BU_RTH')
    IF (LBUDGET_RC) CALL BUDGET (                                        &
         UNPACK(PRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS_3D)*PRHODJ(:,:,:), &
         7,'WETH_BU_RRC')
    IF (LBUDGET_RR) CALL BUDGET (                                        &
         UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
         8,'WETH_BU_RRR')
    IF (LBUDGET_RI) CALL BUDGET (                                        &
         UNPACK(PRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS_3D)*PRHODJ(:,:,:), &
         9,'WETH_BU_RRI')
    IF (LBUDGET_RS) CALL BUDGET (                                        &
         UNPACK(PRSS(:),MASK=GMICRO(:,:,:),FIELD=PRSS_3D)*PRHODJ(:,:,:), &
         10,'WETH_BU_RRS')
    IF (LBUDGET_RG) CALL BUDGET (                                        &
         UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
         11,'WETH_BU_RRG')
    IF (LBUDGET_RH) CALL BUDGET (                                        &
         UNPACK(PRHS(:),MASK=GMICRO(:,:,:),FIELD=PRHS_3D)*PRHODJ(:,:,:), &
         12,'WETH_BU_RRH')
    IF (LBUDGET_SV) THEN
      CALL BUDGET (UNPACK(PCCS(:),MASK=GMICRO(:,:,:),FIELD=PCCS_3D)*PRHODJ(:,:,:), &
           12+NSV_LIMA_NC,'WETH_BU_RSV')
      CALL BUDGET (UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
           12+NSV_LIMA_NR,'WETH_BU_RSV')
      IF (.NOT. LCRYSTAL_SHAPE) THEN
        CALL BUDGET (UNPACK(PCIS(:),MASK=GMICRO(:,:,:),FIELD=PCIS_3D)*PRHODJ(:,:,:), &
                        12+NSV_LIMA_NI,'WETH_BU_RSV')
      ELSE
        DO JSH = 1, NB_CRYSTAL_SHAPE
          CALL BUDGET (UNPACK(PCIS_SHAPE(:,JSH),MASK=GMICRO(:,:,:),FIELD=PCIS_SHAPE_3D(:,:,:,JSH))*PRHODJ(:,:,:), &
                   12+NSV_LIMA_NI+ JSH-1,'WETH_BU_RSV')
        END DO
      END IF  ! LCRYSTAL_SHAPE
    END IF    ! LBUDGET_SV
  END IF ! lbu_enable
!
!
! Partial reconversion of hail to graupel when rc and rh are small    
!
!
!*       3.3   Conversion of the hailstones into graupel
!        -----------------------------------------------
!
  IF ( IHAIL>0 ) THEN
    ZTHRH = 0.01E-3
    ZTHRC = 0.001E-3
    ZZW(:) = 0.0
!
    WHERE( PRHT(:)<ZTHRH .AND. PRCT(:)<ZTHRC .AND. PZT(:)<XTT )
      ZZW(:) = MIN( 1.0,MAX( 0.0,1.0-(PRCT(:)/ZTHRC) ) )
!
! assume a linear percent conversion rate of hail into graupel
!
      ZZW(:)  = PRHS(:) * ZZW(:)
      PRGS(:) = PRGS(:) + ZZW(:)   !   partial conversion
      PRHS(:) = PRHS(:) - ZZW(:)   ! of hail into graupel
!
    END WHERE
  END IF
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_RG) CALL BUDGET (                                    &
        UNPACK(PRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS_3D)*PRHODJ(:,:,:), &
        11,'COHG_BU_RRG')
    IF (LBUDGET_RH) CALL BUDGET (                                    &
        UNPACK(PRHS(:),MASK=GMICRO(:,:,:),FIELD=PRHS_3D)*PRHODJ(:,:,:), &
        12,'COHG_BU_RRH')
  END IF
!
!
!*       3.4    Melting of the hailstones
!
  IF ( IHAIL>0 ) THEN
    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. (PRHS(:)>XRTMIN(7)/PTSTEP) .AND. (PRHT(:)>XRTMIN(7)) .AND. (PZT(:)>XTT) )
      ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
      ZZW(:) = PKA(:)*(XTT-PZT(:)) +                              &
           ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
           *(XESTT-ZZW(:))/(XRV*PZT(:))         )
!
! compute RHMLTR
!
      ZZW(:)  = MIN( PRHS(:), MAX( 0.0,( -ZZW(:) *                     &
                             ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                               X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) -   &
                      ZZW1(:,6)*( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                               ( PRHODREF(:)*XLMTT ) ) )
      PRRS(:) = PRRS(:) + ZZW(:)
      PRHS(:) = PRHS(:) - ZZW(:)
      PTHS(:) = PTHS(:) - ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RHMLTR))
!
      PCRS(:) = MAX( PCRS(:) + ZZW(:)*(XCCH*PLBDAH(:)**XCXH/PRHT(:)),0.0 )
    END WHERE
  END IF ! ihail
!
  IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
    IF (LBUDGET_TH) CALL BUDGET (                                       &
         UNPACK(PTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS_3D)*PRHODJ(:,:,:),&
         4,'HMLT_BU_RTH')
    IF (LBUDGET_RR) CALL BUDGET (                                        &
         UNPACK(PRRS(:),MASK=GMICRO(:,:,:),FIELD=PRRS_3D)*PRHODJ(:,:,:), &
         8,'HMLT_BU_RRR')
    IF (LBUDGET_RH) CALL BUDGET (                                        &
         UNPACK(PRHS(:),MASK=GMICRO(:,:,:),FIELD=PRHS_3D)*PRHODJ(:,:,:), &
         12,'HMLT_BU_RRH')
    IF (LBUDGET_SV) THEN
       CALL BUDGET (  UNPACK(PCRS(:),MASK=GMICRO(:,:,:),FIELD=PCRS_3D)*PRHODJ(:,:,:), &
            12+NSV_LIMA_NR,'HMLT_BU_RSV')
    END IF
  END IF
END IF
!
IF (ALLOCATED(ZAUX))  DEALLOCATE(ZAUX)
IF (ALLOCATED(ZFACT)) DEALLOCATE(ZFACT)
IF (ALLOCATED(ZONEOVER_VAR)) DEALLOCATE(ZONEOVER_VAR)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
