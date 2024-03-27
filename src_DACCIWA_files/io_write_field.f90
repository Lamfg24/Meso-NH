!MNH_LIC Copyright 2016-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Original version:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
MODULE MODE_IO_WRITE_FIELD
!
USE MODD_IO_ll, ONLY: TOUTBAK
USE MODE_FIELD
USE MODE_FMWRIT
!
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE IO_WRITE_FIELDLIST(TPOUTPUT)
!
USE MODE_MODELN_HANDLER, ONLY : GET_CURRENT_MODEL_INDEX
!
IMPLICIT NONE
!
TYPE(TOUTBAK),    INTENT(IN)  :: TPOUTPUT !Output structure
!
INTEGER :: IDX
INTEGER :: IMI
INTEGER :: JI
!
IMI = GET_CURRENT_MODEL_INDEX()
!
DO JI = 1,SIZE(TPOUTPUT%NFIELDLIST)
  IDX = TPOUTPUT%NFIELDLIST(JI)
  SELECT CASE (TFIELDLIST(IDX)%NDIMS)
    !
    !0D output
    !
    CASE (0)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !0D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D logical fields')
          END IF
        !
        !0D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D integer fields')
          END IF
        !
        !0D logical
        !
        CASE (TYPELOG)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_L0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_L0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_L0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D logical fields')
          END IF
        !
        !0D string
        !
        CASE (TYPECHAR)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_C0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_C0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_C0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D character fields')
          END IF
        !
        !0D date/time
        !
        CASE (TYPEDATE)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_T0D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_T0D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_T0D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_T0D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_T0D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 0D date/time fields')
          END IF
        !
        !0D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 0D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !1D output
    !
    CASE (1)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !1D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X1D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X1D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X1D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D real fields')
          END IF
!         !
!         !1D integer
!         !
!         CASE (TYPEINT)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D integer fields')
!           END IF
!         !
!         !1D logical
!         !
!         CASE (TYPELOG)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_L1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_L1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_L1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_L1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D logical fields')
!           END IF
!         !
!         !1D string
!         !
!         CASE (TYPECHAR)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_C1D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C1D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_C1D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_C1D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_C1D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 1D character fields')
!           END IF
        !
        !1D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 1D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !2D output
    !
    CASE (2)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !2D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X2D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X2D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X2D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X2D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X2D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 2D real fields')
          END IF
        !
        !2D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N2D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N2D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N2D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N2D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N2D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not allowed for 2D integer fields')
          END IF
        !
        !2D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 2D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !3D output
    !
    CASE (3)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !3D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X3D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X3D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X3D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 3D real fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X3D(IMI)%DATA)
          END IF
        !
        !3D integer
        !
        CASE (TYPEINT)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_N3D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N3D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_N3D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 3D integer fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_N3D(IMI)%DATA)
          END IF
        !
        !3D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 3D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
    !
    !4D output
    !
    CASE (4)
      SELECT CASE (TFIELDLIST(IDX)%NTYPE)
        !
        !4D real
        !
        CASE (TYPEREAL)
          IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X4D) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X4D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA) ) THEN
            PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X4D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
            STOP
          END IF
          IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
            CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA)
          ELSE
            CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 4D real fields')
            !PW: TODO?: add missing field in TFIELDLIST?
            !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X4D(IMI)%DATA)
          END IF
        !
        !4D other types
        !
        CASE DEFAULT
          PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 4D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
          STOP
      END SELECT
!     !
!     !5D output
!     !
!     CASE (5)
!       SELECT CASE (TFIELDLIST(IDX)%NTYPE)
!         !
!         !5D real
!         !
!         CASE (TYPEREAL)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X5D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X5D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X5D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 5D real fields')
!             !PW: TODO?: add missing field in TFIELDLIST?
!             !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X5D(IMI)%DATA)
!           END IF
!         !
!         !5D other types
!         !
!         CASE DEFAULT
!           PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 5D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!           STOP
!       END SELECT
!     !
!     !6D output
!     !
!     CASE (6)
!       SELECT CASE (TFIELDLIST(IDX)%NTYPE)
!         !
!         !6D real
!         !
!         CASE (TYPEREAL)
!           IF ( .NOT.ALLOCATED(TFIELDLIST(IDX)%TFIELD_X6D) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X6D is NOT allocated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( .NOT.ASSOCIATED(TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA) ) THEN
!             PRINT *,'FATAL: IO_WRITE_FIELDLIST: TFIELD_X6D%DATA is not associated for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!             STOP
!           END IF
!           IF ( TFIELDLIST(IDX)%CLBTYPE == 'NONE' ) THEN
!             CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TFIELDLIST(IDX),TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA)
!           ELSE
!             CALL PRINT_MSG(NVERB_ERROR,'IO','IO_WRITE_FIELDLIST','CLBTYPE/=NONE not (yet) allowed for 6D real fields')
!             !PW: TODO?: add missing field in TFIELDLIST?
!             !CALL IO_WRITE_FIELD_LB(TPOUTPUT%TFILE,TFIELDLIST(IDX),***,TFIELDLIST(IDX)%TFIELD_X6D(IMI)%DATA)
!           END IF
!         !
!         !6D other types
!         !
!         CASE DEFAULT
!           PRINT *,'FATAL: IO_WRITE_FIELDLIST: type not yet supported for 4D output of ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
!           STOP
!       END SELECT
    !
    !Other number of dimensions
    !
    CASE DEFAULT
      PRINT *,'FATAL: IO_WRITE_FIELDLIST: number of dimensions not yet supported for ',TRIM(TFIELDLIST(IDX)%CMNHNAME)
      STOP
  END SELECT
END DO
!
END SUBROUTINE IO_WRITE_FIELDLIST
!
! ##########################################################################
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE IO_WRITE_FIELD_USER(TPOUTPUT)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ##########################################################################
!
USE MODD_CONF_n,      ONLY : NRR,NRRL,NRRI
USE MODD_CST,         ONLY : XRV, XRD, XCPD, XG, XP00
USE MODD_DIM_n,       ONLY : NKMAX
USE MODD_PARAMETERS,  ONLY : JPVEXT
USE MODD_DYN_n,       ONLY : XTSTEP
USE MODD_FIELD_n,     ONLY : XUT, XVT, XWT, XRT, XTHT, XPABST, XCIT
USE MODD_REF_n,       ONLY : XRHODREF
USE MODD_PRECIP_n,    ONLY : XINPRR
USE MODD_GRID_n,      ONLY : XZZ, XZS
USE MODD_DIAG_IN_RUN
USE MODD_MNH_SURFEX_n
USE MODD_SURF_PAR,    ONLY : XUNDEF 
!
!USE MODI_RADAR_RAIN_ICE
USE MODI_UNPACK_SAME_RANK
! 
USE MODE_ll
!
IMPLICIT NONE
!
TYPE(TOUTBAK),    INTENT(IN)  :: TPOUTPUT !Output structure
TYPE(TFIELDDATA)              :: TZFIELD
!
INTEGER                             :: JLOOP
INTEGER                             :: IIU,IJU,IKU,IIB,IJB,IKB,IIE,IJE,IKE ! Arrays bounds
INTEGER                             :: IDIM1,IDIM2                         !
REAL                                :: ZGAMREF
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHETAV
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEMP, ZRARE, ZVDOP, ZRDZDR, ZRKDP
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZWORK2D1, ZWORK2D2
REAL, DIMENSION(:),     ALLOCATABLE :: ZWORK1D
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
CALL GET_INDICE_ll  (IIB,IJB,IIE,IJE)
!
IKU   = NKMAX+2*JPVEXT
IKB   = 1+JPVEXT
IKE   = IKU-JPVEXT
IDIM1 = IIE-IIB+1
IDIM2 = IJE-IJB+1
!
ALLOCATE(ZWORK2D1(IIU,IJU))
ALLOCATE(ZWORK2D2(IIU,IJU))
ALLOCATE(ZWORK1D(IDIM1*IDIM2))
ALLOCATE(ZTHETAV(IIU,IJU,IKU))
ALLOCATE(ZTEMP(IIU,IJU,IKU))
ALLOCATE(ZRARE(IIU,IJU,IKU)) 
ALLOCATE(ZVDOP(IIU,IJU,IKU))
ALLOCATE(ZRDZDR(IIU,IJU,IKU))
ALLOCATE(ZRKDP(IIU,IJU,IKU))
ZRARE = 0.0
!
!
!  1. u-component of wind at lowest physical level
!     --------------------------------------------
!
TZFIELD%CMNHNAME   = 'UTLOW'
TZFIELD%CSTDNAME   = 'x_wind'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_U component of wind at lowest physical level'
TZFIELD%NGRID      = 2
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XUT(:,:,IKB))
!
!
!  2. v-component of wind at lowest physical level
!     --------------------------------------------
!
TZFIELD%CMNHNAME   = 'VTLOW'
TZFIELD%CSTDNAME   = 'y_wind'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_V component of wind at lowest physical level'
TZFIELD%NGRID      = 3
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XVT(:,:,IKB))
!
!
!  3. w-component of wind at lowest physical level
!     --------------------------------------------
!
TZFIELD%CMNHNAME   = 'WTLOW'
TZFIELD%CSTDNAME   = 'z_wind'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_W component of wind at lowest physical level'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XWT(:,:,IKB))
!
!
!  4. potential temperature at lowest physical level
!     ----------------------------------------------
!
TZFIELD%CMNHNAME   = 'THTLOW'
TZFIELD%CSTDNAME   = 'air_potential_temperature'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'K'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_potential temperature at lowest physical level'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XTHT(:,:,IKB))
!
!
!  5. vapor mixing Ratio at lowest physical level
!     -------------------------------------------
!
TZFIELD%CMNHNAME   = 'RVTLOW'
TZFIELD%CSTDNAME   = 'specific_humidity'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'kg kg-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Z_Vapor mixing Ratio at lowest physical level'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XRT(:,:,IKB,1))
!
!
!  6. accumulated Precipitation Rain Rate during timestep
!     ---------------------------------------------------
!
TZFIELD%CMNHNAME   = 'ACPRRSTEP'
TZFIELD%CSTDNAME   = 'rainfall_amount'
TZFIELD%CLONGNAME  = ''
TZFIELD%CUNITS     = 'kg m-2'
TZFIELD%CDIR       = ''
TZFIELD%CCOMMENT   = 'X_Y_ACcumulated Precipitation Rain Rate during timestep'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
!XACPRR is multiplied by 1000. to convert from m to kg m-2 (water density is assumed to be 1000 kg m-3)
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XINPRR*XTSTEP*1.0E3)
!
!
!  7. mean sea level pressure (computation based on write_lfifm1_for_diag.f90)
!     ------------------------------------------------------------------------
!
IF(NRR > 0) THEN
  ! compute the ratio : 1 + total water mass / dry air mass
  ZTHETAV(:,:,:) = 1. + XRT(:,:,:,1)
  DO JLOOP = 2,1+NRRL+NRRI
    ZTHETAV(:,:,:) = ZTHETAV(:,:,:) + XRT(:,:,:,JLOOP)
  END DO
  ! compute the virtual potential temperature when water is present in any form
  ZTHETAV(:,:,:) = XTHT(:,:,:) * (1.+XRT(:,:,:,1)*XRV/XRD) / ZTHETAV(:,:,:)
ELSE
  ! compute the virtual potential temperature when water is absent
  ZTHETAV(:,:,:) = XTHT(:,:,:)
END IF
!
ZGAMREF=-6.5E-3
! Exner function at the first mass point
ZWORK2D1(:,:) = (XPABST(:,:,IKB) /XP00)**(XRD/XCPD)
! virtual temperature at the first mass point
ZWORK2D1(:,:) = ZWORK2D1(:,:) * ZTHETAV(:,:,IKB)
!  virtual temperature at ground level
ZWORK2D1(:,:) = ZWORK2D1(:,:) - ZGAMREF*((XZZ(:,:,IKB)+XZZ(:,:,IKB+1))/2.-XZS(:,:))
!  virtual temperature at sea level
ZWORK2D2(:,:) = ZWORK2D1(:,:) - ZGAMREF*XZS(:,:)
!  average underground virtual temperature
ZWORK2D2(:,:) = 0.5*(ZWORK2D1(:,:)+ZWORK2D2(:,:))
!  surface pressure
ZWORK2D1(:,:) = ( XPABST(:,:,IKB) + XPABST(:,:,IKB-1) )*.5
!  sea level pressure (hPa)
ZWORK2D2(:,:) = 1.E-2*ZWORK2D1(:,:)*EXP(XG*XZS(:,:)/(XRD*ZWORK2D2(:,:)))
!
TZFIELD%CMNHNAME   = 'MSLP'
TZFIELD%CSTDNAME   = 'air_pressure_at_sea_level'
TZFIELD%CLONGNAME  = 'MSLP'
TZFIELD%CUNITS     = 'hPa'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_Mean Sea Level Pressure'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D2)
!
!
!  8. ZON10M
!     ------
!
ZWORK2D1(:,:)             = XUNDEF
ZWORK2D1(IIB:IIE,IJB:IJE) = RESHAPE(YSURF_CUR%DU%XZON10M(:), (/IDIM1,IDIM2/))
!
!
TZFIELD%CMNHNAME   = 'ZON10M'
TZFIELD%CSTDNAME   = 'wind at 10 m'
TZFIELD%CLONGNAME  = 'ZON10M'
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_wind_at_10m'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D1(:,:))
!
!
!  8. MER10M
!     ------
!
ZWORK2D1(:,:)             = XUNDEF
ZWORK2D1(IIB:IIE,IJB:IJE) = RESHAPE(YSURF_CUR%DU%XMER10M(:), (/IDIM1,IDIM2/))
!
TZFIELD%CMNHNAME   = 'MER10M'
TZFIELD%CSTDNAME   = 'wind at 10 m'
TZFIELD%CLONGNAME  = 'MER10M'
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_wind_at_10m'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D1(:,:))
!

!
!  8. sensible heat flux
!     ------------------
!
!TZFIELD%CMNHNAME   = 'H'
!TZFIELD%CSTDNAME   = 'sensible_heat_flux'
!TZFIELD%CLONGNAME  = 'H'
!TZFIELD%CUNITS     = 'W m-2'
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_sensible heat flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_H(:,:))
!
!
!  9. sensible heat flux
!     ------------------
!
!TZFIELD%CMNHNAME   = 'LE'
!TZFIELD%CSTDNAME   = 'latent_heat_flux'
!TZFIELD%CLONGNAME  = 'H'
!TZFIELD%CUNITS     = 'W m-2'
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_sensible heat flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_LE(:,:))
!
!
!
!  10. sensible heat flux
!      ------------------
!
!TZFIELD%CMNHNAME   = 'SFTH'
!TZFIELD%CSTDNAME   = 'sensible_heat_flux'
!TZFIELD%CLONGNAME  = 'SFTH'
!TZFIELD%CUNITS     = ''
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_sensible heat flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_SFTH(:,:))
!
!
!  11. latent heat flux
!      ------------------
!
!TZFIELD%CMNHNAME   = 'SFTQ'
!TZFIELD%CSTDNAME   = 'latent_heat_flux'
!TZFIELD%CLONGNAME  = 'SFTQ'
!TZFIELD%CUNITS     = ''
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_sensible heat flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_SFTQ(:,:))
!
!  12. u-momentum flux
!      ------------------
!
!TZFIELD%CMNHNAME   = 'SFU'
!TZFIELD%CSTDNAME   = 'u_momentum_flux'
!TZFIELD%CLONGNAME  = 'SFU'
!TZFIELD%CUNITS     = ''
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_u_momentum_flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_SFU(:,:))
!
!
!  13. v-momentum flux
!      ------------------
!
!TZFIELD%CMNHNAME   = 'SFV'
!TZFIELD%CSTDNAME   = 'v_momentum_flux'
!TZFIELD%CLONGNAME  = 'SFV'
!TZFIELD%CUNITS     = ''
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_v_momentum flux'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,XCURRENT_SFV(:,:))
!
!  14. sea surface temperature
!      ------------------
!
ZWORK1D(:)    = 0.0
ZWORK2D1(:,:) = XUNDEF
!
IF (SIZE(YSURF_CUR%SM%S%XSST(:)) /= 0) CALL UNPACK_SAME_RANK(YSURF_CUR%U%NR_SEA(:),YSURF_CUR%SM%S%XSST(:),ZWORK1D(:),XUNDEF)
ZWORK2D1(IIB:IIE,IJB:IJE) = RESHAPE(ZWORK1D(:), (/IDIM1,IDIM2/))
!
TZFIELD%CMNHNAME   = 'SST'
TZFIELD%CSTDNAME   = 'sea_surface_temperature'
TZFIELD%CLONGNAME  = 'SST'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_sea_surface_temperature'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D1(:,:))
!
!  15. sea surface u-current
!      ------------------
!
ZWORK1D(:)=0.0
ZWORK2D1(:,:)=XUNDEF
!
IF (SIZE(YSURF_CUR%SM%S%XUMER(:)) /= 0) CALL UNPACK_SAME_RANK(YSURF_CUR%U%NR_SEA(:),YSURF_CUR%SM%S%XUMER(:),ZWORK1D(:),XUNDEF)
!
ZWORK2D1(IIB:IIE,IJB:IJE) = RESHAPE(ZWORK1D(:), (/IDIM1,IDIM2/))
!
TZFIELD%CMNHNAME   = 'UCUR'
TZFIELD%CSTDNAME   = 'sea_surface_u_current'
TZFIELD%CLONGNAME  = 'UCUR'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_sea_surface_u_current'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D1(:,:))
!
!  16. sea surface v-current
!      ------------------
!
ZWORK1D(:)=0.0
ZWORK2D1(:,:)=XUNDEF
!
IF (SIZE(YSURF_CUR%SM%S%XVMER(:)) /= 0) CALL UNPACK_SAME_RANK(YSURF_CUR%U%NR_SEA(:),YSURF_CUR%SM%S%XVMER(:),ZWORK1D(:),XUNDEF)
!
ZWORK2D1(IIB:IIE,IJB:IJE) = RESHAPE(ZWORK1D(:), (/IDIM1,IDIM2/))
!
TZFIELD%CMNHNAME   = 'VCUR'
TZFIELD%CSTDNAME   = 'sea_surface_v_current'
TZFIELD%CLONGNAME  = 'VCUR'
TZFIELD%CUNITS     = ''
TZFIELD%CDIR       = 'XY'
TZFIELD%CCOMMENT   = 'X_Y_sea_surface_v_current'
TZFIELD%NGRID      = 1
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 2
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZWORK2D1(:,:))
!
!
!  17. add radar reflectivity
!      ----------------------
!
!ZTEMP(:,:,:) = XTHT(:,:,:) * (XPABST(:,:,:)/XP00)**(XRD/XCPD)
!
!CALL RADAR_RAIN_ICE (XRT, XCIT, XRHODREF, ZTEMP, ZRARE, ZVDOP, ZRDZDR, ZRKDP)
!
!TZFIELD%CMNHNAME   = 'RARE_k23'
!TZFIELD%CSTDNAME   = 'radar_reflectivity_at_k23'
!TZFIELD%CLONGNAME  = 'RARE_k23'
!TZFIELD%CUNITS     = 'dBz'
!TZFIELD%CDIR       = 'XY'
!TZFIELD%CCOMMENT   = 'X_Y_radar_reflecticity_at_k23'
!TZFIELD%NGRID      = 1
!TZFIELD%NTYPE      = TYPEREAL
!TZFIELD%NDIMS      = 2
!TZFIELD%LTIMEDEP   = .TRUE.
!CALL IO_WRITE_FIELD(TPOUTPUT%TFILE,TZFIELD,ZRARE(:,:,23))
! 
!*    DEALLOCATE
!
!DEALLOCATE(ZTEMP)
!DEALLOCATE(ZRARE)
!DEALLOCATE(ZVDOP)
!DEALLOCATE(ZRDZDR)
!DEALLOCATE(ZRKDP)
DEALLOCATE(ZTHETAV)
DEALLOCATE(ZWORK2D1)
DEALLOCATE(ZWORK2D2)
DEALLOCATE(ZWORK1D)
!
END SUBROUTINE IO_WRITE_FIELD_USER
!
END MODULE MODE_IO_WRITE_FIELD
