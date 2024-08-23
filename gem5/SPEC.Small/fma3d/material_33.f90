!! MATERIAL(*)%PVAL(1:22) Usage
!!
!! 1:Density 2:Bulk_Ln 3:Bulk_Qd 4:HG_Viscosity 5:HG_Stiffness
!!
!! PVAL 1/2/3 6/7/8 10/20/30 11/21/31 22  25/35/45 32   33   36  17/27 38
!!                  40/50    41/51                                 37
!!  6:  K1    D1    E        E        E1Ln  E_aa  E          Ko        A0
!!  7:        D2    Nu       Nu       E1Qd  E_bb  Nu         Kinf      A1
!!  8:              Lambda   Lambda   E2Ln  E_cc  Lam        Kdec      A2
!!  9:              G        G        E2Qd  Nu_ba G     G    G0   Gftn G
!! 10:  Yield                Yield    Fric  Nu_ca Yield Bulk Ginf Bulk Bulk
!! 11:  Eplas                Eplas    Bulk  Nu_cb Ep_ki F_a  Gdec      Pftn
!! 12:  H                    H        Fcmp  G_ab  Ep_is F_b            2G
!! 13:  Beta                 Beta     Gmod  G_ac  pparm F-c            4G/3
!! 14:                                Gcof  G_bc  rparm C_r            Pfra
!! 15:                                                  C_p
!! 16:  p_exp                p_exp                      V_d
!! 17:  D_exp                D_exp                      V_w
!! 18:                                                  Ltyp
!! 19:                                                  Tmax
!!
      SUBROUTINE MATERIAL_33_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:58
!!
!! The following initialization operations are performed:
!!
!! 1. Based on the tension cutoff stress TMAX, tension positive,
!!    the tension cutoff parameter TCUT is found.
!! 2. Based on the initial value of X, XINT, the starting value
!!    of the state variable EL, ELSTRT, is found.
!! 3. Based on the constants for the failure surface CA, CB and
!!    CC, the intersection with the J1 axis, FCUT, is found.
!!
!! Secant iteration is used to find initial value of EL.
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
      DATA                                                                     &
     &          MAXIT /20/                                                     &
     &          CONV /1.0D-10/
!!
!! Statement function for exponential function with small/large argument.
!!
      EXPS(arg) = EXP(MIN(87.0D+0,MAX(-87.0D+0,arg)))
!!
!! Statement function for SJ2 failure envelope as a function of SJ1.
!!
      FAILURE(SJ1) = CA - CC * EXPS(CB*SJ1)
!!
!! Statement functions for cap behavior.
!! (CAPL is current value of CAP(EL).)
!! (  XL is current value of   X(EL).)
!!
      CAP(EL) = MIN(0.0D+0,EL)
      R(CAPL) = CR
      X(EL)   = EL - R(CAP(EL)) * FAILURE(EL)
      EVP(XL) = CW * (EXPS(CD*XL) - ONE)
      SJ2C(SJ1,XL,CAPL) =                                                      &
     &          SQRT((XL-CAPL)*(XL-CAPL)-(SJ1-CAPL)*(SJ1-CAPL))/R(CAPL)
!!
!! Statement functions for elastic moduli.
!! (EVPL is current value of EVP(XL).)
!!
      BMOD(SJ1,EVPL) = BULK
      SMOD(SJ2,EVPL) = SHEAR
!!
!! INITIALIZE MATERIAL CONSTANTS.
!!
      SHEAR = MATERIAL(MatID)%PVAL(9)
      BULK  = MATERIAL(MatID)%PVAL(10)
      CA    = MATERIAL(MatID)%PVAL(11)
      CB    = MATERIAL(MatID)%PVAL(12)
      CC    = MATERIAL(MatID)%PVAL(13)
      CR    = MATERIAL(MatID)%PVAL(14)
      XINT  = MATERIAL(MatID)%PVAL(15)
      CD    = MATERIAL(MatID)%PVAL(16)
      CW    = MATERIAL(MatID)%PVAL(17)
      LTYPE = MAX (1,MIN (2,NINT (MATERIAL(MatID)%PVAL(18))))
      TCUT  = 3.0D+0 * MATERIAL(MatID)%PVAL(19)
!!
      MATERIAL(MatID)%PVAL(8) = BULK - (2.0D+0/3.0D+0) * SHEAR
!!
      MATERIAL(MatID)%PVAL(18) = DBLE (LTYPE)
      MATERIAL(MatID)%PVAL(21) = TCUT
      GEOP = 0.0
!!
!! XINT is reset so that within the convergence criteria XINT is negative.
!!
      XINT = MIN (XINT,-CONV*BMOD(0.0D+0,0.0D+0))
      ELZ = XINT
      EL = XINT - 0.1D+0 * MAX (ABS (XINT),FAILURE(XINT))
      FZ = X(ELZ) - XINT
      COMP = CONV * FAILURE(CAP(EL) + CAP(ELZ))
      DO 100 IT = 1,MAXIT
        FL = X(EL)-XINT
        IF (ABS(FL).LT.COMP .OR. FL.EQ.FZ) GO TO 200
        ELN = EL - FL * (EL - ELZ) / (FL - FZ)
        ELZ = EL
        FZ = FL
        EL = ELN
  100   CONTINUE
  200   ELSTRT = EL
      MATERIAL(MatID)%PVAL(22) = ELSTRT
      FCUT = MAX (0.0D+0,ELSTRT)
      DEL = FAILURE(FCUT)
      IF (DEL .EQ. 0.0) GO TO 600
      DO 300 IT = 1,MAXIT
        EL = FCUT + DEL
        FZ = FAILURE(EL)
        IF (FZ .LT. 0.0) GO TO 400
        DEL = 10.0D+0 * DEL
        FCUT = EL
  300   CONTINUE
      GO TO 600
  400   CONTINUE
      DO 500 IT = 1,MAXIT
        FL = FAILURE(FCUT)
        IF (ABS(FL).LT.COMP .OR. FL.EQ.FZ) GO TO 600
        ELN = FCUT - FL * (FCUT - EL) / (FL - FZ)
        EL = FCUT
        FZ = FL
        FCUT = ELN
  500   CONTINUE
  600   MATERIAL(MatID)%PVAL(20) = FCUT
!!
      RETURN
!!
      ENTRY MATERIAL_33_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! Initialize stress point w/ starting value of cap parameter, and  tension
!! "debt" S1cut which must be overcome before material can reload in
!! compression after a tensile break.
!!
      Ioff = Isv - 1
      STATE_VARIABLES(Ioff+1) = MATERIAL(MatID)%PVAL(22)
      STATE_VARIABLES(Ioff+2) = 0.0
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_33                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 25-NOV-1991 08:50:49
!!
!! FINITE STRAIN ELASTIC-PLASTIC STRAIN HARDENING SOIL CAP
!!
!! This subroutine computes the current value of the stress STRESS(6), the
!! longitudinal sound speed CL, and the shear sound speed CS in the current
!! state, as a function of the stretching Dxx,...Dyz, the spin, the past
!! history data EL, and the old stress state.  This approach to plasicity is
!! based on hypoelastic concepts.
!!
!!      STATE_VARIABLES(1) = State variable, EL
!!      STATE_VARIABLES(2) = Tension "debt" for reloading after tension break.
!!      STATE_VARIABLES(3) = Inelastic process type.
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
!!      MATERIAL%BULK           ! Bulk modulus
!!      MATERIAL%SHEAR          ! Shear modulus
!!      MATERIAL%CAP75_A        ! Failure envelope constant
!!      MATERIAL%CAP75_B        ! Failure envelope constant
!!      MATERIAL%CAP75_C        ! Failure envelope constant
!!      MATERIAL%CAP75_R        ! Cap envelope constant
!!      MATERIAL%XINT           ! Cap envelope constant
!!      MATERIAL%CAP75_D        ! State function constant
!!      MATERIAL%CAP75_W        ! State function constant
!!      MATERIAL%LTYPE          ! State function constant
!!      MATERIAL%TMAX           ! Tensile cutoff, tension positive
!!      MATERIAL%FCUT           ! Failure envelope intersection w/ J1
!!      MATERIAL%TCUT           ! Tensile cutoff parameter, TCUT = 3*TMAX
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          Sinc(6),                                                       &
     &          Einc(6),                                                       &
     &          STRESS(6),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)

      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! This CAP75 common block has been retained in order to make as few changes
!! to CAP75 as possible.
!!
      COMMON /PROP/                                                            &
     &          TCUT,FCUT,BULK,SHEAR,CA,CB,CC,CD,CR,CW,CL2,CS2,LTYPE
!!
!! Initialize constants:
!!      P1 = Elastic Lame parameter lambda
!!      P2 = Elastic Lame parameter mu (shear modulus)
!!      P3 = Twice the shear modulus, 2mu
!!
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
      P3 = P2 + P2
      LTYPE = NINT (MATERIAL(MatID)%PVAL(18))
      FCUT  = MATERIAL(MatID)%PVAL(20)
      TCUT  = MATERIAL(MatID)%PVAL(21)
      BULK  = MATERIAL(MatID)%PVAL(10)
      SHEAR = P2
      CA = MATERIAL(MatID)%PVAL(11)
      CB = MATERIAL(MatID)%PVAL(12)
      CC = MATERIAL(MatID)%PVAL(13)
      CD = MATERIAL(MatID)%PVAL(16)
      CR = MATERIAL(MatID)%PVAL(14)
      CW = MATERIAL(MatID)%PVAL(17)
!!
!! Generate the required rotation operator (and new stretching components
!! if Ipolard equals one or two).
!!
      Igenrot = 1
      Ipolard = CONTROL%POLARD
!!
      CALL SYMMETRIC_TENSOR_ROTATION                                           &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
!!
!! Rotate stress from the configuration at time n to the configuration
!! at time n+1 (Ipolard=0).
!!
      IF (Ipolard .EQ. 0) CALL SYMMETRIC_TENSOR_ROTATION                       &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
!!
!! Save old stress state.
!!
      DO i = 1,6
        Sinc(i) = STRESS(i)
      ENDDO
!!
      SIGX = STRESS(1)
      SIGY = STRESS(2)
      SIGZ = STRESS(3)
      SXY  = STRESS(4)
      SXZ  = STRESS(5)
      SYZ  = STRESS(6)
!!
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D+0*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! Strain increment.
!!
      DEPSX = DTnext*Dxx
      DEPSY = DTnext*Dyy
      DEPSZ = DTnext*Dzz
      DEXY  = DTnext*Dxy
      DEXZ  = DTnext*Dxz
      DEYZ  = DTnext*Dyz
!!
!! Save strain increment.
!!
      Einc(1) = DEPSX
      Einc(2) = DEPSY
      Einc(3) = DEPSZ
      Einc(4) = DEXY
      Einc(5) = DEXZ
      Einc(6) = DEYZ
!!
!! Retrieve state variables and geostatic pressure.
!!
      EL    = STATE_VARIABLES(1)
      S1cut = STATE_VARIABLES(2)
      GEOP  = 0.0
!!
!! Update stress.
!!
      CALL CAP75                                                               &
     &          (                                                              &
     &          SIGX,SIGY,SIGZ,SXY,SXZ,SYZ,GEOP,DEPSX,DEPSY,                   &
     &          DEPSZ,DEXY,DEXZ,DEYZ,EL,S1cut,MTYPE,IT,NOCON                   &
     &          )
!!
!! Place new stresses in global storage.
!!
      STRESS(1) = SIGX
      STRESS(2) = SIGY
      STRESS(3) = SIGZ
      STRESS(4) = SXY
      STRESS(5) = SXZ
      STRESS(6) = SYZ
      STATE_VARIABLES(1) = EL
      STATE_VARIABLES(2) = S1cut
      STATE_VARIABLES(3) = DBLE(MTYPE)
!!
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D+0*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! If the end-of-the-interval Rashid approximate polar decomposition is
!! being used (Ipolard=1,2), rotate stress from the configuration at
!! time n to the configuration at time n+1. Note: the rotation operator R
!! has been saved from the generation entry (Igenrot=1).
!!
      IF (Ipolard .NE. 0) THEN
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
      ENDIF
!!
!! Sound speeds squared * RHO(t)
!!
      IF (CONTROL%SUBCYC .EQ. 0) THEN
        DO i = 1,6
          Sinc(i) = STRESS(i) - Sinc(i)
        ENDDO
        CALL LOADING_MODULI                                                    &
     &    (SOUND_SPEED%RCL2,SOUND_SPEED%RCS2,P1,P2,Sinc,Einc)
      ELSE
        SOUND_SPEED%RCL2 = P1 + P3
        SOUND_SPEED%RCS2 = P2
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE CAP75                                                         &
     &          (                                                              &
     &          SIGX,SIGY,SIGZ,SXY,SXZ,SYZ,GEOP,DEPSX,DEPSY,                   &
     &          DEPSZ,DEXY,DEXZ,DEYZ,EL,S1cut,MTYPE,IT,NOCON                   &
     &          )
!!
!!        The CAP75 model is a plasticity model defined by a non-
!!        softening convex yield surface and a plastic strain rate
!!        vector that is normal to the yield surface in stress space.
!!        The yield surface is defined by means of a failure envelope
!!        and a hardening cap.
!!
!!                                  By
!!                    Ivan S. Sandler AND David Rubin
!!               WEIDLINGER Associates, 333 Seventh Avenue
!!                    New York, New York, 10001, USA
!!
!!        I. S. Sandler and D. Rubin,"An Algorithm and a Modular
!!        Subroutine for the Cap Model," INTERNATIONAL JOURNAL FOR
!!        NUMERICAL AND ANALYTICAL METHODS IN GEOMECHANICS, Vol. 3,
!!        Pages 173-186 (1979).
!!
!!                              Adapted By
!!                            Samuel W. Key
!!
!!        INPUT:
!!        The input stresses and strain increments are 3-dimensional.
!!        GEOP, the geostatic presssure. (Compression is positive.)
!!        EL, the state variable at the beginning of the increment.
!!        S1cut, tension "debt" for reloading after a tension break.
!!
!!        OUTPUT:
!!        New stresses.
!!        EL, the state variable at the end of the increment.
!!        S1cut, tension "debt" for reloading after a tension break.
!!        MTYPE, indicates the form of material behavior which occurred:
!!           0 = Tension Cutoff, 1 = Elastic, 2 = Failure, 3 = Cap Action.
!!        IT, equals the number of iterations for the cap.
!!        NOCON, equals 1 when the cap iterations fail to converge within
!!           the maximum permitted, NIT.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      COMMON/PROP/TCUT,FCUT,BULK,SHEAR,CA,CB,CC,CD,CR,CW,CL2,CS2,LTYPE
      DATA        CONV/1.0D-3/, NIT/60/
!!
!! Statement function for exponential with small/large argument.
!!
      EXPS(arg) = EXP (MIN(87.0D+0,MAX(-87.0D+0,arg)))
!!
!! Statement function for SJ2 failure envelope as a function of SJ1.
!!
      FAILURE(SJ1) = CA-CC*EXPS(CB*SJ1)
!!
!! Statement functions for cap behavior.
!! (CAPL is current value of CAP(EL).)
!! (  XL is current value of   X(EL).)
!!
      CAP(EL) = MIN(0.0D+0,EL)
      R(CAPL) = CR
      X(EL)   = EL-R(CAP(EL))*FAILURE(EL)
      EVP(XL) = CW*(EXPS(CD*XL)-ONE)
      SJ2C(SJ1,XL,CAPL) =                                                      &
     & SQRT((XL-CAPL)*(XL-CAPL)-(SJ1-CAPL)*(SJ1-CAPL))/R(CAPL)
!!
!! Statement functions for elastic moduli.
!! (EVPL is current value of EVP(XL).)
!!
      BMOD(SJ1,EVPL) = BULK
      SMOD(SJ2,EVPL) = SHEAR
!!
      IT = 0
      NOCON = 0
      DEV = DEPSX+DEPSY+DEPSZ
      DEVB3 = DEV/3.0D+0
      DEXX = DEPSX-DEVB3
      DEYY = DEPSY-DEVB3
      PRESS = -(SIGX+SIGY+SIGZ)/3.0D+0
      SXX = SIGX+PRESS
      SYY = SIGY+PRESS
      SJ1I = -3.0D+0*(PRESS-GEOP)
      SJ2I = SQRT(SXX*SXX+SYY*SYY+SXX*SYY+SXY*SXY+SXZ*SXZ+SYZ*SYZ)
      CAPL = CAP(EL)
      XL = X(EL)
      EVPI = EVP(XL)
!!
!! Elastic material properties.
!!
      THREEK = 3.0D+0*BMOD(SJ1I,EVPI)
      TWOG = 2.0D+0*SMOD(SJ2I,EVPI)
!!
!! Elastic trial stresses.
!!
      SJ1 = THREEK*DEV+SJ1I
      SXX = SXX+TWOG*DEXX
      SYY = SYY+TWOG*DEYY
      SXY = SXY+TWOG*DEXY
      SXZ = SXZ+TWOG*DEXZ
      SYZ = SYZ+TWOG*DEYZ
      RATIO = ONE
      MTYPE = 1
!!
!! Tensile coding.
!!
      TENCUT = MIN(FCUT,TCUT+3.0D+0*GEOP)
      S1NEW = MIN(TENCUT,SJ1+MAX(0.0D+0,S1cut))
      S1cut = S1cut+SJ1-S1NEW
      SJ1 = S1NEW
      IF (SJ1 .LT. TENCUT) GO TO 10
      SJ1 = TENCUT
      RATIO = 0.0
      MTYPE = 0
      IF (LTYPE.EQ.2 .OR. EL.GE.0.0) GO TO 200
!!
!! Tension dilatancy coding.
!!
      ELL = MAX(0.0D+0,EL+CONV*FAILURE(EL))
      XLL = X(ELL)
      DENOM = EVP(XLL)-EVPI
      IF (DENOM .LE. 0.0) THEN
       EL = 0.0
      ELSE
       DEVP = DEV-(SJ1-SJ1I)/THREEK
       EL = EL+DEVP*(ELL-EL)/DENOM
       EL = MIN(0.0D+0,EL)
      ENDIF
      GO TO 200
!!
!! Check failure envelope.
!!
   10 CONTINUE
      SJ2 = SQRT(SXX*SXX+SYY*SYY+SXX*SYY+SXY*SXY+SXZ*SXZ+SYZ*SYZ)
      IF (SJ1 .LT. CAPL) GO TO 30
      TMISES = SJ2C(CAPL,XL,CAPL)
      FJ1 = FAILURE(SJ1)
      FF = SJ2-MIN(FJ1,TMISES)
      IF (FF .LE. 0.0) GO TO 200
!!
!! Failure surface calculation.
!!
      MTYPE = 2
      IF (FJ1 .LT. TMISES) THEN
       DFDJ1 = (FJ1-FAILURE(SJ1+CONV*SJ2))/(CONV*SJ2)
      ELSE
       DFDJ1 = 0.0
      ENDIF
      DEVP = 3.0D+0*DFDJ1*FF/(3.0D+0*THREEK*DFDJ1*DFDJ1+0.5D+0*TWOG)
      SJ1 = SJ1-THREEK*DEVP
!!
!! Dilatancy and corner coding.
!!
      IF (LTYPE.EQ.1 .AND. EL.LT.0.0 .AND. SJ1.GT.CAPL) THEN
      ELL = CAP(SJ1)
      XLL = X(ELL)
      IF (DEVP .GT. 0.0) THEN
        DEVPT = MAX(DEVP,EVP(XLL)-EVPI)
        EL = EL+(ELL-EL)*DEVP/DEVPT
      ENDIF
      ELSE
      SJ1 = MAX(SJ1,CAPL)
      ENDIF
      FJ1 = FAILURE(SJ1)
      RATIO = MIN(FJ1,TMISES)/SJ2
      GO TO 200
!!
!! Cap calculation.
!!
   30 CONTINUE
      IF (SJ1 .LT. XL) GO TO 40
      IF (SJ2 .LE. SJ2C(SJ1,XL,CAPL)) GO TO 200
   40 CONTINUE
      SJ1E = SJ1
      SJ2E = SJ2
      ELL = EL
      ELR = SJ1E
      IF (SJ1E .LE. XL) THEN
      FL = (EL-SJ1E)/(EL-XL)
      ELSE
      FL = 2.0D+0*SJ2E/(SJ2E+SJ2C(SJ1E,XL,CAPL))-ONE
      ENDIF
      XR = X(ELR)
      SJ1R = SJ1E-THREEK*(EVP(XR)-EVPI)
      FR = (XR-SJ1R)/(ELR-XR)
      COMP = CONV*FAILURE((FL*XR-FR*XL)/(FL-FR))
      IF (ABS(SJ1)+SJ2 .LT. COMP) GO TO 200
      MTYPE = 3
      FOLD = 0.0
      DO 190 IT = 1,NIT
      EL = (FL*ELR-FR*ELL)/(FL-FR)
      XL = X(EL)
      DEVP = EVP(XL)-EVPI
      SJ1 = SJ1E-THREEK*DEVP
      CAPL = CAP(EL)
      IF (SJ1 .LE. XL)   FC = (EL-SJ1)/(EL-XL)
      IF (SJ1 .GE. CAPL) FC = (XL-SJ1)/(CAPL-XL)
      IF (SJ1.GT.XL .AND. SJ1.LT.CAPL) THEN
        SJ2 = SJ2C(SJ1,XL,CAPL)
        DELJ1 = CONV*(XL-SJ1)
        IF (SJ1+DELJ1 .NE. SJ1) THEN
          DESP = (DEVP/6.0D+0)*(DELJ1/(SJ2-SJ2C(SJ1+DELJ1,XL,CAPL)))
        ELSE
          DESP = 0.0
        ENDIF
        SJ2TRY = SJ2+TWOG*DESP
        FC = (SJ2E-SJ2TRY)/(SJ2E+SJ2TRY)
        IF (ABS(SJ2E-SJ2TRY) .LE. COMP) GO TO 195
        IF (FC.GT.0.0 .AND. SJ1-CAPL.GE.DELJ1) GO TO 195
      ENDIF
      IF (FC .LE. 0.0) THEN
        ELR = EL
        FR = FC
        IF (FOLD .LT. 0.0) FL = 0.5D+0*FL
      ELSE
        ELL = EL
        FL = FC
        IF (FOLD .GT. 0.0) FR = 0.5D+0*FR
        FOLD = FC
      ENDIF
  190 CONTINUE
      NOCON = 1
      SJ1 = MAX(SJ1,XL)
      IF (SJ1 .GT. CAP(ELR)) SJ1 = CAPL
      SJ2 = MIN(SJ2E,SJ2C(SJ1,XL,CAPL))
!!
  195 CONTINUE
      IF (SJ2E .NE. 0.0) THEN
      RATIO = SJ2/SJ2E
      ELSE
      RATIO = 0.0
      ENDIF
!!
  200 CONTINUE
      SXX = SXX*RATIO
      SYY = SYY*RATIO
      SXY = SXY*RATIO
      SXZ = SXZ*RATIO
      SYZ = SYZ*RATIO
      PRESS = -SJ1/3.0D+0+GEOP
      SIGX = SXX-PRESS
      SIGY = SYY-PRESS
      SIGZ = -3.0D+0*PRESS-SIGX-SIGY
!!
!! Sound speeds squared * RHO(T).
!!
!!!    CL2 = (THREEK+TWOG+TWOG)/3.0D+0
!!!    CS2 = 0.5D+0*TWOG
!!
      RETURN
      END
!!_
      SUBROUTINE CAP75_ORG (SIGX,SIGY,SIGZ,SXY,SXZ,SYZ,GEOP,DEPSX,DEPSY,       &
     &                     DEPSZ,DEXY,DEXZ,DEYZ,EL,S1cut,MTYPE,IT,NOCON)
!!
!!        The CAP75 model is a plasticity model defined by a non-
!!        softening convex yield surface and a plastic strain rate
!!        vector that is normal to the yield surface in stress space.
!!        The yield surface is defined by means of a failure envelope
!!        and a hardening cap.
!!
!!                                  By
!!                    Ivan S. Sandler AND David Rubin
!!               WEIDLINGER Associates, 333 Seventh Avenue
!!                    New York, New York, 10001, USA
!!
!!        I. S. Sandler and D. Rubin,"An Algorithm and a Modular
!!        Subroutine for the Cap Model," INTERNATIONAL JOURNAL FOR
!!        NUMERICAL AND ANALYTICAL METHODS IN GEOMECHANICS, Vol. 3,
!!        Pages 173-186 (1979).
!!
!!                              Adapted By
!!                            Samuel W. Key
!!
!!        INPUT.
!!        The input stresses and strain increments are 3-dimensional.
!!        GEOP, the geostatic presssure. (Compression is positive.)
!!        EL, the state variable at the beginning of the increment.
!!        S1cut, tension "debt" for reloading after a tension break.
!!
!!        OUTPUT.
!!        New stresses.
!!        EL, the state variable at the end of the increment.
!!        S1cut, tension "debt" for reloading after a tension break.
!!        MTYPE, indicates the form of material behavior which occurred:
!!          0 = Tension Cutoff, 1 = Elastic, 2 = Failure, 3 = Cap Action.
!!        IT, equals the number of iterations for the cap.
!!        NOCON, equals 1 when the cap iterations fail to converge within
!!          the maximum permitted, NIT.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      COMMON/PROP/TCUT,FCUT,BULK,SHEAR,CA,CB,CC,CD,CR,CW,CL2,CS2,LTYPE
      DATA CONV/1.0D-3/, NIT/60/
!!
!! Statement function for exponential with small/large argument.
!!
      EXPS(Z) = EXP(MIN(87.0D+0,MAX(-87.0D+0,Z)))
!!
!! Statement function for SJ2 failure envelope as a function of SJ1.
!!
      FAILURE(SJ1) = CA-CC*EXPS(CB*SJ1)
!!
!! Cap statement functions.
!! (CAPL is current value of BIGL(EL).)
!! (XL   is current value of X(EL).   )
!!
      BIGL(EL) = MIN(0.0D+0,EL)
      R(CAPL) = CR
      X(EL) = EL-R(BIGL(EL))*FAILURE(EL)
      EVP(XL) = CW*(EXPS(CD*XL)-ONE)
      SJ2C(SJ1,XL,CAPL) =                                                      &
     & SQRT((XL-CAPL)*(XL-CAPL)-(SJ1-CAPL)*(SJ1-CAPL))/R(CAPL)
!!
!! Elastic moduli statement functions.
!! (EV is current value of EVP(XL).)
!!
      BMOD(SJ1,EV) = BULK
      SMOD(SJ2,EV) = SHEAR
!!
      IT = 0
      NOCON = 0
      DEV = DEPSX+DEPSY+DEPSZ
      DEVB3 = DEV/3.0D+0
      DEXX = DEPSX-DEVB3
      DEYY = DEPSY-DEVB3
      PRESS = -(SIGX+SIGY+SIGZ)/3.0D+0
      SXX = SIGX+PRESS
      SYY = SIGY+PRESS
      SJ1I = -3.0D+0*(PRESS-GEOP)
      SJ2I = SQRT(SXX*SXX+SYY*SYY+SXX*SYY+SXY*SXY+SXZ*SXZ+SYZ*SYZ)
      CAPL = BIGL(EL)
      XL = X(EL)
      EVPI = EVP(XL)
!!
!! Elastic material properties.
!!
      THREEK = 3.0D+0*BMOD(SJ1I,EVPI)
      TWOG = 2.0D+0*SMOD(SJ2I,EVPI)
!!
!! Elastic trial stresses.
!!
      SJ1 = THREEK*DEV+SJ1I
      SXX = SXX+TWOG*DEXX
      SYY = SYY+TWOG*DEYY
      SXY = SXY+TWOG*DEXY
      SXZ = SXZ+TWOG*DEXZ
      SYZ = SYZ+TWOG*DEYZ
      RATIO = ONE
      MTYPE = 1
!!
!! Tensile coding.
!!
      TENCUT = MIN(FCUT,TCUT+3.0D+0*GEOP)
      S1NEW = MIN(TENCUT,SJ1+MAX(0.0D+0,S1cut))
      S1cut = S1cut+SJ1-S1NEW
      SJ1 = S1NEW
      IF (SJ1 .LT. TENCUT) GO TO 10
      SJ1 = TENCUT
      RATIO = 0.0
      MTYPE = 0
      IF (LTYPE.EQ.2 .OR. EL.GE.0.0) GO TO 200
!!
!! Tension dilatancy coding.
!!
      ELL = MAX(0.0D+0,EL+CONV*FAILURE(EL))
      XLL = X(ELL)
      DENOM = EVP(XLL)-EVPI
      IF (DENOM .GT. 0.0) GO TO 5
      EL = 0.0
      GO TO 200
    5 DEVP = DEV-(SJ1-SJ1I)/THREEK
      EL = EL+DEVP*(ELL-EL)/DENOM
      EL = MIN(0.0D+0,EL)
      GO TO 200
!!
!! Check failure envelope.
!!
   10 CONTINUE
      SJ2 = SQRT(SXX*SXX+SYY*SYY+SXX*SYY+SXY*SXY+SXZ*SXZ+SYZ*SYZ)
      IF (SJ1 .LT. CAPL) GO TO 30
      TMISES = SJ2C(CAPL,XL,CAPL)
      FJ1 = FAILURE(SJ1)
      FF = SJ2-MIN(FJ1,TMISES)
      IF (FF .LE. 0.0) GO TO 200
!!
!! Failure surface calculation.
!!
      MTYPE = 2
      DFDJ1 = 0.0
      IF (FJ1.LT.TMISES) DFDJ1 = (FJ1-FAILURE(SJ1+CONV*SJ2))/(CONV*SJ2)
      DEVP = 3.0D+0*DFDJ1*FF/(3.0D+0*THREEK*DFDJ1*DFDJ1+0.5D+0*TWOG)
      SJ1 = SJ1-THREEK*DEVP
!!
!! Dilatancy and corner coding.
!!
      IF (LTYPE.EQ.1 .AND. EL.LT.0.0 .AND. SJ1.GT.CAPL) GO TO 60
      SJ1 = MAX(SJ1,CAPL)
      GO TO 70
   60 CONTINUE
      ELL = BIGL(SJ1)
      XLL = X(ELL)
      IF (DEVP .LE. 0.0) GO TO 70
      DEVPT = MAX(DEVP,EVP(XLL)-EVPI)
      EL = EL+(ELL-EL)*DEVP/DEVPT
   70 CONTINUE
      FJ1 = FAILURE(SJ1)
      RATIO = MIN(FJ1,TMISES)/SJ2
      GO TO 200
!!
!! Cap calculation.
!!
   30 CONTINUE
      IF (SJ1 .LT. XL) GO TO 40
      IF (SJ2 .LE. SJ2C(SJ1,XL,CAPL)) GO TO 200
   40 CONTINUE
      SJ1E = SJ1
      SJ2E = SJ2
      ELL = EL
      ELR = SJ1E
      IF (SJ1E .LE. XL) FL = (EL-SJ1E)/(EL-XL)
      IF (SJ1E .GT. XL) FL = 2.0D+0*SJ2E/(SJ2E+SJ2C(SJ1E,XL,CAPL))-ONE
      XR = X(ELR)
      SJ1R = SJ1E-THREEK*(EVP(XR)-EVPI)
      FR = (XR-SJ1R)/(ELR-XR)
      COMP = CONV*FAILURE((FL*XR-FR*XL)/(FL-FR))
      IF (ABS(SJ1)+SJ2 .LT. COMP) GO TO 200
      MTYPE = 3
      FOLD = 0.0
      DO 190 IT = 1,NIT
      EL = (FL*ELR-FR*ELL)/(FL-FR)
      XL = X(EL)
      DEVP = EVP(XL)-EVPI
      SJ1 = SJ1E-THREEK*DEVP
      CAPL = BIGL(EL)
      IF (SJ1 .LE. XL)   FC = (EL-SJ1)/(EL-XL)
      IF (SJ1 .GE. CAPL) FC = (XL-SJ1)/(CAPL-XL)
      IF (SJ1.LE.XL .OR. SJ1.GE.CAPL) GO TO 300
      SJ2 = SJ2C(SJ1,XL,CAPL)
      DELJ1 = CONV*(XL-SJ1)
      DESP = 0.0
      IF (SJ1+DELJ1 .NE. SJ1) DESP = (DEVP/6.0D+0)*(DELJ1/(SJ2-                   &
     &                                SJ2C(SJ1+DELJ1,XL,CAPL)))
      SJ2TRY = SJ2+TWOG*DESP
      FC = (SJ2E-SJ2TRY)/(SJ2E+SJ2TRY)
      IF (ABS(SJ2E-SJ2TRY) .LE. COMP) GO TO 195
      IF (FC.GT.0.0 .AND. SJ1-CAPL.GE.DELJ1) GO TO 195
  300 IF (FC .GT. 0.0) GO TO 320
      ELR = EL
      FR = FC
      IF (FOLD .LT. 0.0) FL = 0.5D+0*FL
      FOLD = FC
      GO TO 190
  320 CONTINUE
      ELL = EL
      FL = FC
      IF (FOLD .GT. 0.0) FR = 0.5D+0*FR
      FOLD = FC
  190 CONTINUE
      NOCON = 1
      SJ1 = MAX(SJ1,XL)
      IF (SJ1 .GT. BIGL(ELR)) SJ1 = CAPL
      SJ2 = MIN(SJ2E,SJ2C(SJ1,XL,CAPL))
  195 RATIO = 0.0
      IF (SJ2E .NE. 0.0) RATIO = SJ2/SJ2E
  200 CONTINUE
      SXX = SXX*RATIO
      SYY = SYY*RATIO
      SXY = SXY*RATIO
      SXZ = SXZ*RATIO
      SYZ = SYZ*RATIO
      PRESS = -SJ1/3.0D+0+GEOP
      SIGX = SXX-PRESS
      SIGY = SYY-PRESS
      SIGZ = -3.0D+0*PRESS-SIGX-SIGY
!!
!! Sound speeds squared * RHO(T).
!!
      CL2 = (THREEK+TWOG+TWOG)/3.0D+0
      CS2 = 0.5D+0*TWOG
!!
      RETURN
      END
