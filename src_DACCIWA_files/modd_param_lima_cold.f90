!     ###########################
      MODULE MODD_PARAM_LIMA_COLD
!     ###########################
!
!!****  *MODD_PARAM_LIMA_COLD* - declaration of some descriptive parameters and
!!                               microphysical factors extensively used in 
!!                               the LIMA cold scheme.
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  *Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      T. Hoarau (LACy)     jui. 2016  add CIBU
!!      M. Claeys (LACy)     mar. 2019  add ice crystal shapes
!!      JP Pinty  (LA)       jul. 2019  add RSDF
!-------------------------------------------------------------------------------
!
IMPLICIT NONE 
!
!*       1.   DESCRIPTIVE PARAMETERS
!             ----------------------
!
!     Declaration of microphysical constants, including the descriptive
!     parameters for the raindrop and the ice crystal habits, and the 
!     parameters relevant of the dimensional distributions.
!
!         m(D)    = XAx * D**XBx      : Mass-MaxDim relationship
!         v(D)    = XCx * D**XDx      : Fallspeed-MaxDim relationship
!         N(Lbda) = XCCx * Lbda**XCXx : NumberConc-Slopeparam relationship
!         XF0x, XF1x, XF2x            : Ventilation factors
!         XC1x                        : Shape parameter for deposition
!
!              and
!
!         XALPHAx, XNUx                        : Generalized GAMMA law 
!         Lbda = XLBx * (r_x*rho_dref)**XLBEXx : Slope parameter of the 
!                                                distribution law
!
REAL,SAVE :: XLBEXI,XLBI              ! Prist. ice     distribution parameters
REAL,SAVE :: XLBEXS,XLBS              ! Snow/agg.      distribution parameters
!
REAL,SAVE :: XAI,XBI,XC_I,XDI         ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL,SAVE ::                           XF0IS,XF1IS    ! (large Di vent. coef.)
REAL,SAVE :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
!
REAL,SAVE :: XLBDAS_MAX               ! Max values allowed for the shape
                                      ! parameter of snow
!
REAL, DIMENSION(:), SAVE, ALLOCATABLE :: XLBEXI_SHAPE, & ! pristine ice distrib. parameters
                                         XLBI_SHAPE      ! when several shapes
REAL, DIMENSION(4), SAVE :: XAI_SHAPE,  XBI_SHAPE, & ! a and b parameters of the mass-diam relationship
                            XC_I_SHAPE, XDI_SHAPE, &
                            XC1I_SHAPE !relation masse-D et vitesse-D
!
CHARACTER(LEN=8),DIMENSION(5),PARAMETER &
                              :: CLIMA_COLD_NAMES=(/'CICE    ','CIFNFREE','CIFNNUCL', &
                                                        'CCNINIMM','CCCNNUCL'/)
                                 ! basenames of the SV articles stored
                                 ! in the binary files
                                 !with IF:Ice-nuclei Free (nonactivated IFN by Dep/Cond)
                                 !     IN:Ice-nuclei Nucleated (activated IFN by Dep/Cond)
                                 !     NI:Nuclei Immersed (activated IFN by Imm)
                                 !     HF:Homogeneous Freezing
CHARACTER(LEN=3),DIMENSION(5),PARAMETER &
                              :: CLIMA_COLD_CONC=(/'NI ','NIF','NIN','NNI','NNH'/)!for DIAG
CHARACTER(LEN=8),DIMENSION(4),PARAMETER &
                              :: CLIMA_SHAPE_NAMES=(/'CICE    ','CICECOL ','CICEIRR ','CICEDEN '/)
!                              :: CLIMA_SHAPE_NAMES=(/'PLA','COL','IRR','DEN'/)  
                                 ! basenames of the SV articles stored
                                 ! in the binary files
                                 !with PLA : plates
                                 !     COL : Columns
                                 !     IRR : Irregulars
                                 !     DEN : Dendrites
CHARACTER(LEN=5),DIMENSION(4),PARAMETER &
                              :: CLIMA_SHAPE_CONC=(/'N_PLA','N_COL','N_IRR','N_DEN'/)!for DIAG
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS
!             ---------------------
!
REAL,SAVE :: XFSEDRI,XFSEDCI,                  & ! Constants for sedimentation
             XFSEDS, XEXSEDS                     ! fluxes of ice and snow
REAL, DIMENSION(:),ALLOCATABLE, SAVE :: XFSEDCI_SHAPE  ! Constant per habit: N
REAL, DIMENSION(:),ALLOCATABLE, SAVE :: XFSEDRI_SHAPE  ! Constant per habit: mr
REAL, DIMENSION(:),ALLOCATABLE, SAVE :: XFSEDRI_TOT_SHAPE  ! Constant per habit: mr
!
REAL,SAVE :: XNUC_DEP,XEXSI_DEP,XEX_DEP,       & ! Constants for heterogeneous
             XNUC_CON,XEXTT_CON,XEX_CON,       & ! ice nucleation : DEP et CON
             XMNU0                               ! mass of nucleated ice crystal
!
REAL,SAVE :: XRHOI_HONH,XCEXP_DIFVAP_HONH,     & ! Constants for homogeneous
             XCOEF_DIFVAP_HONH,XRCOEF_HONH,    & ! haze freezing : HHONI
             XCRITSAT1_HONH,XCRITSAT2_HONH,    &
             XTMIN_HONH,XTMAX_HONH,            &
             XDLNJODT1_HONH,XDLNJODT2_HONH,    &
             XC1_HONH,XC2_HONH,XC3_HONH
!
REAL,SAVE :: XC_HONC,XR_HONC,                  & ! Constants for homogeneous
             XTEXP1_HONC,XTEXP2_HONC,          & ! droplet freezing : CHONI
             XTEXP3_HONC,XTEXP4_HONC,          &
             XTEXP5_HONC
!
REAL,SAVE :: XCSCNVI_MAX, XLBDASCNVI_MAX,      &
             XRHORSMIN,                        &
             XDSCNVI_LIM, XLBDASCNVI_LIM,      & ! Constants for snow
             XC0DEPSI,XC1DEPSI,                & ! sublimation conversion to
             XR0DEPSI,XR1DEPSI                   ! pristine ice : SCNVI
!
REAL, DIMENSION(:), ALLOCATABLE, SAVE          &
          :: XC0DEPSI_SHAPE,XC1DEPSI_SHAPE,    & ! sublimation conversion to
             XR0DEPSI_SHAPE,XR1DEPSI_SHAPE       ! pristine ice : SCNVI
!
REAL,SAVE :: XSCFAC,                           & ! Constants for the Bergeron
             X0DEPI,X2DEPI,                    & ! Findeisen process and
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS     ! deposition
!
REAL, DIMENSION(:), ALLOCATABLE, SAVE :: X0DEPI_SHAPE, & ! vapor deposition on ice
                                         X2DEPI_SHAPE
!
!
REAL,SAVE :: XDICNVS_LIM, XLBDAICNVS_LIM,      & ! Constants for pristine ice
             XC0DEPIS,XC1DEPIS,                & ! deposition conversion to
             XR0DEPIS,XR1DEPIS                   ! snow : ICNVS
!
REAL, DIMENSION(:), ALLOCATABLE, SAVE          &
          :: XC0DEPIS_SHAPE,XC1DEPIS_SHAPE,    & ! sublimation conversion to
             XR0DEPIS_SHAPE,XR1DEPIS_SHAPE       ! pristine ice : SCNVI
!
REAL,SAVE :: XCOLEXIS,                         & ! Constants for snow 
             XAGGS_CLARGE1,XAGGS_CLARGE2,      & ! aggregation : AGG
             XAGGS_RLARGE1,XAGGS_RLARGE2
!
REAL,DIMENSION(:), ALLOCATABLE, SAVE :: XAGGS_RLARGE1_SHAPE, & ! Constants for snow 
                                        XAGGS_RLARGE2_SHAPE   
!??????????????????
REAL,SAVE :: XKER_ZRNIC_A1,XKER_ZRNIC_A2         ! Long-Zrnic Kernels (ini_ice_coma)
!
REAL,SAVE :: XSELFI,XCOLEXII                     ! Constants for pristine ice
                                                 ! self-collection (ini_ice_coma)
!
REAL,SAVE :: XAUTO3, XAUTO4,                   & ! Constants for pristine ice
             XLAUTS,   XLAUTS_THRESHOLD,       & ! autoconversion : AUT
             XITAUTS, XITAUTS_THRESHOLD,       & ! (ini_ice_com) 
             XTEXAUTI
!
REAL,SAVE :: XCONCI_MAX                          ! Limitation of the pristine 
                                   ! ice concentration (init and grid-nesting) 
REAL,SAVE :: XFREFFI  ! Factor to compute the cloud ice effective radius
REAL, DIMENSION(:), ALLOCATABLE, SAVE :: XFREFFI_SHAPE !++cb--
!
! Constants for Ice-Ice Collision : CIBU
!
REAL, SAVE :: XDCSLIM_CIBU_MIN,                & ! aggregates min diam. : 0.2 mm
              XDCSLIM_CIBU_MAX,                & ! aggregates max diam. : 1.0 mm
              XDCGLIM_CIBU_MIN,                & ! graupel min diam. : 2 mm
              XGAMINC_BOUND_CIBU_SMIN,         & ! Min val. of Lbda_s*dlim
              XGAMINC_BOUND_CIBU_SMAX,         & ! Max val. of Lbda_s*dlim
              XGAMINC_BOUND_CIBU_GMIN,         & ! Min val. of Lbda_g*dlim
              XGAMINC_BOUND_CIBU_GMAX,         & ! Max val. of Lbda_g*dlim
              XCIBUINTP_S,XCIBUINTP1_S,        & !
              XCIBUINTP2_S,                    & !
              XCIBUINTP_G,XCIBUINTP1_G,        & !
              XFACTOR_CIBU_NI,XFACTOR_CIBU_RI, & ! Factor for final CIBU Eq.
              XMOMGG_CIBU_1,XMOMGG_CIBU_2,     & ! Moment computation
              XMOMGS_CIBU_1,XMOMGS_CIBU_2,     &
              XMOMGS_CIBU_3
!
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE        &
                       :: XGAMINC_CIBU_S,      & ! Tab.incomplete Gamma function
                          XGAMINC_CIBU_G         ! Tab.incomplete Gamma function
!
! Constants for Raindrop shattering : RDSF
!
REAL, SAVE :: XDCRLIM_RDSF_MIN,                & ! raindrops min diam. : 0.2 mm
              XGAMINC_BOUND_RDSF_RMIN,         & ! Min val. of Lbda_r*dlim
              XGAMINC_BOUND_RDSF_RMAX,         & ! Max val. of Lbda_r*dlim
              XRDSFINTP_R,XRDSFINTP1_R,        & !
              XFACTOR_RDSF_NI,                 & ! Factor for final RDSF Eq.
              XMOMGR_RDSF
!
REAL, DIMENSION(:), SAVE, ALLOCATABLE          &
                       :: XGAMINC_RDSF_R         ! Tab.incomplete Gamma function
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_COLD
