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
      SUBROUTINE MATERIAL_32_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:58
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Lmod,                                                          &
     &          STRESS(6)
      DATA                                                                     &
     &          TROOT /0.81649658092773/   ! SQRT (2/3)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lamé parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((1.0 + Prat) * (1.0 - 2.0*Prat))
      Gmod = Ymod / (2.0 + 2.0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the input linear isotropic hardening modulus Eiso into the plastic
!! hardening modulus Bi, and the input linear kinematic hardening modulus Ekin
!! into the plastic hardening modulus Bk.
!!
      Ekin = MATERIAL(MatID)%PVAL(11)
      Eiso = MATERIAL(MatID)%PVAL(12)

      Bkin = (Ymod * Ekin) / (Ymod - Ekin)
      Biso = (Ymod * Eiso) / (Ymod - Eiso)

      MATERIAL(MatID)%PVAL(11) = Bkin
      MATERIAL(MatID)%PVAL(12) = Biso
!!
      RETURN
!!
      ENTRY MATERIAL_32_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! Initialize stress point with starting value of yield surface radius.
!!
      Ioff = Isv - 1
      Yield = MATERIAL(MatID)%PVAL(10)
      STATE_VARIABLES(Ioff+7) = TROOT * Yield
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_32                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 26-MAY-1991 13:27:50
!!
!! FINITE STRAIN ELASTIC-PLASTIC, SMOOTH STRAIN HARDENING BEHAVIOR
!!
!! Module computes the current value of the stress STRESS(1:6),
!! the density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient), the past history
!! data STATE_VARIABLES(1:8), and the old stress state.
!!
!!      STATE_VARIABLES(1) = Yield surface center Sigma_XX
!!      STATE_VARIABLES(2) = Yield surface center Sigma_YY
!!      STATE_VARIABLES(3) = Yield surface center Sigma_ZZ
!!      STATE_VARIABLES(4) = Yield surface center Sigma_XY
!!      STATE_VARIABLES(5) = Yield surface center Sigma_XZ
!!      STATE_VARIABLES(6) = Yield surface center Sigma_YZ
!!      STATE_VARIABLES(7) = Yield surface radius SQRT(2/3) x Yield_Stress
!!      STATE_VARIABLES(8) = Accumulated effective plastic strain
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!        R. D. Krieg and S. W. Key, "On the Accurate Representation
!!        of Large Strain Non-Proportional Plastic Flow in Ductile
!!        Materials," CONSTITUTIVE EQUATIONS: Macro and Computational
!!        Aspects, edited by K. J. Willam, The American Society of
!!        Mechanical Engineers, New York, (1984).
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)

      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
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
!! Rotate the center of the yield surface STATE_VARIABLES(1:6).
!!
      IF (Ipolard .EQ. 0) CALL SYMMETRIC_TENSOR_ROTATION                       &
     &          (                                                              &
     &          STATE_VARIABLES,Igenrot,Ipolard,DTnext                         &
     &          )
!!
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! Update stress.
!!
      CALL MATERIAL_32_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES(1),   & ! Yield surface center
     &          STATE_VARIABLES(7),   & ! Yield surface radius
     &          STATE_VARIABLES(8),   & ! Effective plastic strain
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*DTnext) *                       &
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
!!
!! Rotate the center of the yield surface STATE_VARIABLES(1:6).
!!
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STATE_VARIABLES,Igenrot,Ipolard,DTnext                         &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_32_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,YSC,YSR,EPS,                                            &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,MatID               &
     &          )
!!
!! ELASTIC-PLASTIC WITH SMOOTH STRAIN HARDENING BEHAVIOR
!!
!! This subroutine is based on a numerical procedure by R. D. Krieg for
!! integrating the incremental equations of plasticity.  The numerical method
!! is a generalization of a method used in wave codes. The theory and further
!! references to it may be found in:
!!
!!        R. D. Krieg and S. W. Key, "On the Accurate Representation
!!        of Large Strain Non-Proportional Plastic Flow in Ductile
!!        Materials," CONSTITUTIVE EQUATIONS: Macro and Computational
!!        Aspects, edited by K. J. Willam, The American Society of
!!        Mechanical Engineers, New York, (1984).
!!
!!      MATERIAL%Ymod  = Young's modulus, (not used in calculation)
!!      MATERIAL%Prat  = Poisson's ratio, QR
!!      MATERIAL%Lmod  = Elastic Lame parameter Lambda, P1
!!      MATERIAL%Gmod  = Elastic Lame parameter µ (shear modulus), P2
!!      MATERIAL%Yield = Yield stress, QS
!!      MATERIAL%Ekin  = Kinematic hardening modulus, QA
!!      MATERIAL%Eiso  = Isotropic hardening modulus, QK
!!      MATERIAL%p     = Hardening parameter p, QP
!!      MATERIAL%r     = Hardening parameter r, QR
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0)) :: STRESS(6),YSC(6),Sinc(6),Einc(6)

      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: MAXSTEPS
      REAL(KIND(0D0)), SAVE :: THIRD,ROOT3,TROOT
!!
!! Obtain maximum number of integration substeps allowed.
!!
      IF (FIRST) THEN
        MAXSTEPS = NINT (PARAMVALUE%Max_Steps)
        THIRD = (1.0D+0 / 3.0D+0)
        ROOT3 = SQRT (1.0D+0 / 3.0D+0)
        TROOT = SQRT (2.0D+0 / 3.0D+0)
        FIRST = .FALSE.
      ENDIF
!!
!! Initialize constants:
!!      P1 = Elastic Lame parameter lambda
!!      P2 = Elastic Lame parameter µ (shear modulus)
!!      P3 = Twice the shear modulus, 2µ
!!      Blk = Elastic bulk modulus
!!      YSR = Yield surface radius
!!      YSC = Yield surface center
!!
      P1  = MATERIAL(MatID)%PVAL(8)
      P2  = MATERIAL(MatID)%PVAL(9)
      P3  = P2 + P2
      Blk = P1 + Third * P3
!!
      QA = MATERIAL(MatID)%PVAL(11)
      QK = MATERIAL(MatID)%PVAL(12)
      QP = MATERIAL(MatID)%PVAL(13)
      QR = MATERIAL(MatID)%PVAL(14)
!!
!! Save old stress state.
!!
      DO i = 1,6
        Sinc(i) = STRESS(i)
      ENDDO
!!
!! Compute strain increments.
!!
      Einc(1) = DTnext*Dxx
      Einc(2) = DTnext*Dyy
      Einc(3) = DTnext*Dzz
      Einc(4) = DTnext*Dxy
      Einc(5) = DTnext*Dxz
      Einc(6) = DTnext*Dyz
!!
!! Compute trial elastic step.
!!
      Dkk = Dxx + Dyy + Dzz
      Q1 = THIRD*Dkk
      Dxx = Dxx - Q1
      Dyy = Dyy - Q1
      Dzz = Dzz - Q1
      PRESS = THIRD*(STRESS(1) + STRESS(2) + STRESS(3))
      S1 = STRESS(1) + DTnext*(P3*Dxx) - PRESS
      S2 = STRESS(2) + DTnext*(P3*Dyy) - PRESS
      S3 = STRESS(3) + DTnext*(P3*Dzz) - PRESS
      S4 = STRESS(4) + DTnext*(P3*Dxy)
      S5 = STRESS(5) + DTnext*(P3*Dxz)
      S6 = STRESS(6) + DTnext*(P3*Dyz)
      PRESS = PRESS + DTnext*Blk*Dkk
      IF (YSR .EQ. 0.0) GO TO 400
!!
!! Evaluate yield criteria.
!!
      SX1 = S1 - YSC(1)
      SX2 = S2 - YSC(2)
      SX3 = S3 - YSC(3)
      SX4 = S4 - YSC(4)
      SX5 = S5 - YSC(5)
      SX6 = S6 - YSC(6)
      SXM2 = SX1*SX1 + SX2*SX2 + SX3*SX3                                       &
     &     + SX4*SX4 + SX4*SX4 + SX5*SX5                                       &
     &     + SX5*SX5 + SX6*SX6 + SX6*SX6
      SXM = SQRT(SXM2)
      PHI = SXM2 - YSR*YSR
!!
!! Check yield condition.
!!
      IF (PHI .LE. 0.0) GO TO 400
!!
!! Integrate plasticity equations.
!!
      XID = SX1*Dxx + SX2*Dyy + SX3*Dzz                                        &
     &    + SX4*Dxy + SX4*Dxy + SX5*Dxz                                        &
     &    + SX5*Dxz + SX6*Dyz + SX6*Dyz
      DD  = Dxx*Dxx + Dyy*Dyy + Dzz*Dzz                                        &
     &    + Dxy*Dxy + Dxy*Dxy + Dxz*Dxz                                        &
     &    + Dxz*Dxz + Dyz*Dyz + Dyz*Dyz
!!
!! Find plastic portion of time step and yield surface contact stress.
!!
      DTP = DTnext
      ARG = XID*XID - PHI*DD
      IF (ARG .GE. 0.0) DTP = (XID-SQRT(ARG))/(P3*DD)
      NSTEP = 0
      DTI = DTP
      DTR = 0.0
      CS1 = S1 - DTP*(P3*Dxx)
      CS2 = S2 - DTP*(P3*Dyy)
      CS3 = S3 - DTP*(P3*Dzz)
      CS4 = S4 - DTP*(P3*Dxy)
      CS5 = S5 - DTP*(P3*Dxz)
      CS6 = S6 - DTP*(P3*Dyz)
  100 CX1 = CS1 - YSC(1)
      CX2 = CS2 - YSC(2)
      CX3 = CS3 - YSC(3)
      CX4 = CS4 - YSC(4)
      CX5 = CS5 - YSC(5)
      CX6 = CS6 - YSC(6)
      CXM2 = CX1*CX1 + CX2*CX2 + CX3*CX3                                       &
     &     + CX4*CX4 + CX4*CX4 + CX5*CX5                                       &
     &     + CX5*CX5 + CX6*CX6 + CX6*CX6
      CXM = SQRT(CXM2)
      VC1 = CX1/CXM
      VC2 = CX2/CXM
      VC3 = CX3/CXM
      VC4 = CX4/CXM
      VC5 = CX5/CXM
      VC6 = CX6/CXM
      VE1 = SX1/SXM
      VE2 = SX2/SXM
      VE3 = SX3/SXM
      VE4 = SX4/SXM
      VE5 = SX5/SXM
      VE6 = SX6/SXM
  110 VCVE = VC1*VE1 + VC2*VE2 + VC3*VE3                                       &
     &     + VC4*VE4 + VC4*VE4 + VC5*VE5                                       &
     &     + VC5*VE5 + VC6*VE6 + VC6*VE6
      IF (VCVE .GE. 0.97) GO TO 120
      DTI = 0.5*DTI
      DTR = DTR + DTI
      VE1 = SX1 - DTR*(P3*Dxx)
      VE2 = SX2 - DTR*(P3*Dyy)
      VE3 = SX3 - DTR*(P3*Dzz)
      VE4 = SX4 - DTR*(P3*Dxy)
      VE5 = SX5 - DTR*(P3*Dxz)
      VE6 = SX6 - DTR*(P3*Dyz)
      VEM2 = VE1*VE1 + VE2*VE2 + VE3*VE3                                       &
     &     + VE4*VE4 + VE4*VE4 + VE5*VE5                                       &
     &     + VE5*VE5 + VE6*VE6 + VE6*VE6
      VEM = SQRT(VEM2)
      VE1 = VE1/VEM
      VE2 = VE2/VEM
      VE3 = VE3/VEM
      VE4 = VE4/VEM
      VE5 = VE5/VEM
      VE6 = VE6/VEM
      GO TO 110
  120 NSTEP = NSTEP + 1
      VA1 = VC1 + VE1
      VA2 = VC2 + VE2
      VA3 = VC3 + VE3
      VA4 = VC4 + VE4
      VA5 = VC5 + VE5
      VA6 = VC6 + VE6
      VAM2 = VA1*VA1 + VA2*VA2 + VA3*VA3                                       &
     &     + VA4*VA4 + VA4*VA4 + VA5*VA5                                       &
     &     + VA5*VA5 + VA6*VA6 + VA6*VA6
      VAM = SQRT(VAM2)
      VA1 = VA1/VAM
      VA2 = VA2/VAM
      VA3 = VA3/VAM
      VA4 = VA4/VAM
      VA5 = VA5/VAM
      VA6 = VA6/VAM
!!
      ETAZ = YSC(1)*VA1 + YSC(2)*VA2 + YSC(3)*VA3                              &
     &     + YSC(4)*VA4 + YSC(4)*VA4 + YSC(5)*VA5                              &
     &     + YSC(5)*VA5 + YSC(6)*VA6 + YSC(6)*VA6
      VD = VA1*Dxx + VA2*Dyy + VA3*Dzz                                         &
     &   + VA4*Dxy + VA4*Dxy + VA5*Dxz                                         &
     &   + VA5*Dxz + VA6*Dyz + VA6*Dyz
      Q1 = P3*VD*DTI
      Q2 = P3 + QK
      Q3 = EXP (QP*ETAZ/P2)
      Q4 = QR*P2 + QA*Q3
      QD = (P3 + QK)*Q3 + Q4
      ETA1 = ETAZ + (Q1*Q4)/QD
!!
!! Solve for ETA1.
!!
!! Case 1: QP = 0.0
!!
      IF (QP .EQ. 0.0) GO TO 300
!!
!! Case 2: QP non-zero, QA = 0.0
!!
      IF (QA .EQ. 0.0) THEN
        Q5 = Q2/(QR*QP)
        Q8 = P2/QP
        Q9 = -Q1 - ETAZ - Q5*Q3
        Q3 = EXP(QP*ETA1/P2)
        Z1 = Q3
        F1 = Q9 + ETA1 + Q5*Q3
        IK = 0
  210   IK = IK + 1
        Z2 = Z1 - F1 / (Q5 + Q8/Z1)
        F2 = Q9 + Q8*LOG(Z2) + Q5*Z2
  800   CONTINUE
        IF (ABS (F2/Q9) .GT. 1.0D-5 .AND. F2 .NE. F1) THEN
          F1 = F2
          Z1 = Z2
          IF (IK .LT. 10) GO TO 210
          Q6 = (P3 + QK)*Z2 + QR*P2
          ETA1 = Q8*LOG(Z2)
          GO TO 230
        ELSE
          Q6 = (P3 + QK)*Z2 + QR*P2
          ETA1 = Q8*LOG(Z2)
        ENDIF
!!
!! Case 3: QP and QA non-zero.
!!
      ELSE
        Q7 = P2*Q2/(QP*QA)
        Q8 = QP/P2
        Q9 = -Q1 - ETAZ - Q7*LOG(Q4)
        Q3 = EXP(Q8*ETA1)
        E1 = ETA1
        F1 = Q9 + ETA1 + Q7*LOG(QR*P2 + QA*Q3)
        IK = 0
  220   IK = IK + 1
        E2 = E1 - F1 / (1.0 + Q2*Q3/(QR*P2 + QA*Q3))
        Q3 = EXP(Q8*E2)
        F2 = Q9 + E2 + Q7*LOG(QR*P2 + QA*Q3)
        IF (ABS (F2/Q9) .GT. 1.0D-5 .AND. F2 .NE. F1) THEN
          F1 = F2
          E1 = E2
          IF (IK .LT. 10) GO TO 220
          Q6 = (P3 + QK + QA)*Q3 + QR*P2
          ETA1 = E2
          GO TO 230
        ELSE
          Q6 = (P3 + QK + QA)*Q3 + QR*P2
          ETA1 = E2
        ENDIF
      ENDIF
  230 CONTINUE
!!
!! Increment kinematic hardening.
!!
  300 YSC(1) = YSC(1) + (ETA1-ETAZ)*VA1
      YSC(2) = YSC(2) + (ETA1-ETAZ)*VA2
      YSC(3) = YSC(3) + (ETA1-ETAZ)*VA3
      YSC(4) = YSC(4) + (ETA1-ETAZ)*VA4
      YSC(5) = YSC(5) + (ETA1-ETAZ)*VA5
      YSC(6) = YSC(6) + (ETA1-ETAZ)*VA6
!!
!! Increment yield surface radius and effective plastic strain.
!!
      IF (QP .EQ. 0.0) THEN
        DEP = Q1/(P3 + QK + QA)
      ELSE IF (ABS (Q9) .LE. (P2*1.0D-3/QP)) THEN
        DEP = Q1/(P3 + QK + QA + QR*P2/Q3)
      ELSE
        Q11 = P2/(QP*Q9)
        DEP = (Q1/(P3 + QK + QA))*Q11*LOG(Q6/QD)
      ENDIF
      YSR = YSR + ROOT3*QK*DEP
      EPS = EPS + TROOT*DEP
!!
!! Compute subincrement deviatoric stress state.
!!
      S1 = VE1*YSR + YSC(1)
      S2 = VE2*YSR + YSC(2)
      S3 = VE3*YSR + YSC(3)
      S4 = VE4*YSR + YSC(4)
      S5 = VE5*YSR + YSC(5)
      S6 = VE6*YSR + YSC(6)
      IF (NSTEP .GT. MAXSTEPS) THEN
        WRITE (MSG1,'(I8)') MAXSTEPS
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MATERIAL_32_INTEGRATION.001.00'//                       &
     &          MSGL//'More Than'//MSG1//' Substeps Are Required'//            &
     &          MSGL//'To Increment Stress-Strain Model.'//                    &
     &          MSGL//'Please Examine Material Constants.'//                   &
     &          MSGL//'Maybe Increase Parameter "MAT32_STEPS."'                &
     &          )
      ENDIF
      IF (DTR .NE. 0.0) THEN
        DTI = DTR
        DTR = 0.0
        CS1 = S1
        CS2 = S2
        CS3 = S3
        CS4 = S4
        CS5 = S5
        CS6 = S6
        SX1 = S1 + DTI*(P3*Dxx) - YSC(1)
        SX2 = S2 + DTI*(P3*Dyy) - YSC(2)
        SX3 = S3 + DTI*(P3*Dzz) - YSC(3)
        SX4 = S4 + DTI*(P3*Dxy) - YSC(4)
        SX5 = S5 + DTI*(P3*Dxz) - YSC(5)
        SX6 = S6 + DTI*(P3*Dyz) - YSC(6)
        SXM2 = SX1*SX1 + SX2*SX2 + SX3*SX3                                     &
     &       + SX4*SX4 + SX4*SX4 + SX5*SX5                                     &
     &       + SX5*SX5 + SX6*SX6 + SX6*SX6
        SXM = SQRT(SXM2)
        GO TO 100
      ENDIF
!!
  400 CONTINUE
!!
!! Compute final stress state.
!!
      STRESS(1) = S1 + PRESS
      STRESS(2) = S2 + PRESS
      STRESS(3) = S3 + PRESS
      STRESS(4) = S4
      STRESS(5) = S5
      STRESS(6) = S6
      Dxx = Dxx + (THIRD*Dkk)
      Dyy = Dyy + (THIRD*Dkk)
      Dzz = Dzz + (THIRD*Dkk)
!!
!! Sound speeds squared * RHO(t).
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
