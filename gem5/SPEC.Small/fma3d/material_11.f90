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
      SUBROUTINE MATERIAL_11_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:54
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: MatID
      REAL(KIND(0D0)), INTENT(IN) :: STRESS(6)
      INTEGER,         INTENT(IN) :: SecID
      INTEGER,         INTENT(IN) :: Isv
      INTEGER,         INTENT(IN) :: Nsv
!!
!! Local variables.
      REAL(KIND(0D0)) :: Ymod,Ehrd,H
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_11_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = ONE
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_21_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:55
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: MatID
      REAL(KIND(0D0)), INTENT(IN) :: STRESS(3)
      INTEGER,         INTENT(IN) :: SecID
      INTEGER,         INTENT(IN) :: Isv
      INTEGER,         INTENT(IN) :: Nsv
!!
!! Local variables.
      REAL(KIND(0D0)) :: Ymod,Prat,Lmod,Gmod,Ehrd,H
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((ONE + Prat) * (ONE - 2.0D+0*Prat))
      Gmod = Ymod / (2.0D+0 + 2.0D+0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_21_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = ONE
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_31_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:57
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: MatID
      REAL(KIND(0D0)), INTENT(IN) :: STRESS(6)
      INTEGER,         INTENT(IN) :: LupID
      INTEGER,         INTENT(IN) :: Isv
      INTEGER,         INTENT(IN) :: Nsv
!!
!! Local variables.
      REAL(KIND(0D0)) :: Ymod,Prat,Lmod,Gmod,Ehrd,H
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((ONE + Prat) * (ONE - 2.0D+0*Prat))
      Gmod = Ymod / (2.0D+0 + 2.0D+0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_31_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = ONE
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_41_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:01
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: MatID
      INTEGER,         INTENT(IN) :: Ipts
      REAL(KIND(0D0)), INTENT(IN) :: STRESS(6)
      INTEGER,         INTENT(IN) :: SecID
      INTEGER,         INTENT(IN) :: Isv
      INTEGER,         INTENT(IN) :: Nsv
!!
!! Local variables.
      REAL(KIND(0D0)) :: Ymod,Prat,Lmod,Gmod,Ehrd,H
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((ONE + Prat) * (ONE - 2.0D+0*Prat))
      Gmod = Ymod / (2.0D+0 + 2.0D+0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_41_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = ONE
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_51_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:02
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: MatID
      INTEGER,         INTENT(IN) :: Ipts
      REAL(KIND(0D0)), INTENT(IN) :: STRESS(6,Ipts)
      INTEGER,         INTENT(IN) :: SecID
      INTEGER,         INTENT(IN) :: Isv
      INTEGER,         INTENT(IN) :: Nsv
!!
!! Local variables.
      REAL(KIND(0D0)) :: Ymod,Prat,Lmod,Gmod,Ehrd,H
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((ONE + Prat) * (ONE - 2.0D+0*Prat))
      Gmod = Ymod / (2.0D+0 + 2.0D+0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Convert the linear hardening modulus Ehrd into the plastic hardening
!! modulus H.
!!
      Ehrd = MATERIAL(MatID)%PVAL(11)

      H = (Ymod * Ehrd) / (Ymod - Ehrd)

      MATERIAL(MatID)%PVAL(12) = H
!!
      RETURN
!!
      ENTRY MATERIAL_51_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = ONE
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_11                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:28:00
!!
!! FINITE STRAIN UNIAXIAL-STRESS ELASTIC-PLASTIC BEHAVIOR.
!!
!! This subroutine computes the current value of the stress STRESS(1:1),
!! the longitudinal sound speed CL in the current state as a function of
!! the stretching Drr and the old stress state and the past history data
!! STATE_VARIABLES(1:2).
!!
!!      STATE_VARIABLES(1) = Yield surface center, Alpha.
!!      STATE_VARIABLES(2) = Effective plastic strain, Efps
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID

      COMMON /TRUSS/                                                           &
     &          Rx,Ry,Rz,Drr
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D+0*DTnext*Drr)*STRESS
!!
!! Update stress.
!!
      CALL MATERIAL_11_INTEGRATION                                             &
     &  (STRESS,STATE_VARIABLES(1),STATE_VARIABLES(2),DTnext,Drr,MatID)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D+0*DTnext*Drr)*STRESS
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_11_INTEGRATION                                       &
     &          (Sigma,Alpha,Efps,DTnext,Drr,MatID)
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:37:40
!!
!! UNIAXIAL-STRESS ELASTIC-PLASTIC STRAIN HARDENING BEHAVIOR.
!!
!! A material whose tensile stress-strain curve has a constant plastic
!! hardening modulus is described. The hardening is a linear combina-
!! tion of isotropic and kinematic which is controled by QB.
!!
!!              QE = Young's modulus
!!              QS = Initial yield stress
!!              QH = Hardening modulus
!!              QB = Isotropic-kinematic hardening parameter
!!                 = 0.0, Kinematic plastic hardening
!!                 = 1.0, Isotropic plastic hardening
!!                 = Fraction, Combined isotropic-kinematic hardening
!!              QP = Strain rate parameter p
!!              QD = Strain rate parameter D
!!
!!      R. P. Gëol and L .E. Malvern, Biaxial Plastic Simple Waves
!!      With Combined Kinematic And Isotropic Hardening, JOURNAL OF
!!      APPLIED MECHANICS, pages 1100-1106, (1970).
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: Sigma
      REAL(KIND(0D0)), INTENT(INOUT) :: Alpha
      REAL(KIND(0D0)), INTENT(INOUT) :: Efps
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      REAL(KIND(0D0)), INTENT(IN)    :: Drr
      INTEGER,         INTENT(IN)    :: MatID
!!
!! Local variables.
      REAL(KIND(0D0)) :: Kappa,Kinc
!!
!! Initialize constants.
!!
      QE = MATERIAL(MatID)%PVAL(6)
      QS = MATERIAL(MatID)%PVAL(10)
      QH = MATERIAL(MatID)%PVAL(12)
      QB = MATERIAL(MatID)%PVAL(13)
      QP = MATERIAL(MatID)%PVAL(16)
      QD = MATERIAL(MatID)%PVAL(17)
!!
!! Save old stress state.
!!
      Sold = Sigma
!!
!! Compute strain increment.
!!
      Einc = DTnext * Drr
!!
!! Compute trial elastic stress.
!!
      Sigma = Sigma + QE * Einc
!!
!! Test for zero yield stress, (Elastic behavior).
!!
      IF (QS .EQ. 0.0) GO TO 200
!!
!! Current radius of the yield surface.
!!
      Kappa = QS + (QB * QH * Efps)
!!
!! Find distance from center of the yield surface to stress state.
!!
      Xi = Sigma - Alpha
      Ximg = ABS (Xi)
!!
!! Check the yield condition.
!!
      IF (Ximg .LE. Kappa) GO TO 200
!!
!! Compute yield stress dynamic increment Sdyi.
!!
      IF (QP .NE. 0.0 .AND. QD .NE. 0.0) THEN
        Edot = Xi * Drr
        IF (Edot .GT. 0.0) THEN
          Edot = Edot/(Ximg*(ONE+(QH+QH)/(QE+QE+QE)))
          Kinc = QS * ( (Edot/QD) ** (ONE/QP) )
          Kappa = Kappa + Kinc
          IF (Ximg .LE. Kappa) GO TO 200
        ENDIF
      ENDIF
!!
!! Plastic loading condition. Calculate effective plastic strain
!! increment Epli, (constant plastic modulus QH).
!!
      IF (Xi .GE. Kappa) THEN
        dSigma = Xi - Kappa
      ELSE
        dSigma = Xi + Kappa
      ENDIF
      Epli = dSigma / (QE + QH)
!!
!! Find new values of stress Sigma, yield surface center Alpha, and
!! total effective plastic strain Efps.
!!
      Sigma = Sigma - QE * Epli
      Alpha = Alpha + (ONE - QB) * QH * Epli
      Efps  = Efps  + ABS (Epli)
!!
  200 CONTINUE
!!
!! Sound speeds squared * rho(t)
!!
      IF (CONTROL%SUBCYC .GT. 0 .OR. Einc .EQ. 0.0) THEN
        SOUND_SPEED%RCL2 = QE
        SOUND_SPEED%RCS2 = 0.0
      ELSE
        SOUND_SPEED%RCL2 = ABS ((Sigma - Sold) / Einc)
        SOUND_SPEED%RCS2 = 0.0
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_21                                                   &
     &          (STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 24-APR-1991 19:09:33
!!
!! FINITE STRAIN, PLANE-STRESS ELASTIC-PLASTIC STRAIN-HARDENING BEHAVIOR.
!!
!!      STATE_VARIABLES(1) = Yield surface center Sigma_RR
!!      STATE_VARIABLES(2) = Yield surface center Sigma_SS
!!      STATE_VARIABLES(3) = Yield surface center Sigma_RS
!!      STATE_VARIABLES(4) = Effective plastic strain
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(3)
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID
!!
!! Local variables.
      REAL(KIND(0D0)), SAVE :: LOCAL_STRESS(6)
      REAL(KIND(0D0)), SAVE :: YIELD_SURFACE_CENTER(4)
      REAL(KIND(0D0)), SAVE :: EFFECTIVE_PLASTIC_STRAIN

      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching, spin
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Gr,Gs,                  & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
!! Initialize material constants.
!!      P1 = Elastic Lame parameter, Lambda
!!      P2 = Elastic Lame parameter, Mu (shear modulus)
!!      P3 = Twice the shear modulus
!!
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
      P3 = P2 + P2
!!
!! Enforce plane stress condition by prescribing stretching in the normal
!! (or t) direction.  This development assumes that the initial value of
!! Dtt is zero.
!!
      Dtt = -P1*(Drr+Dss)/(P1+P3)
!!
!! Set transverse shear strains to zero.
!!
      Drt = 0.0
      Dst = 0.0
!!
!! "De-rotate" stress components to account for the fact that the
!! "co-rotational" axes have been redrawn from the last time step
!! along the same material fiber. The required rotational rate comes
!! from Vyx (=Dxy-Wxy). (The rotation based on Wxy contains the
!! finite rotations of the element relative to the co-rotational
!! coordinate axes, but does not account for the fact that the
!! co-rotational axes rotated with the element.)
!!
!! Rotate the stress from time n to time n+1. This is a proper orthognal
!! rotation using the rotation obtained from the polar decomposition of
!! the deformation gradient.
!!
      Xi = dBeta + DTnext * (Drs-Wrs)
      CB = COS (Xi)
      SB = SIN (Xi)
      S2B = SB*SB
      C2B = CB*CB
      SCB = SB*CB
      Q1 = (STRESS(3)*SCB) + (STRESS(3)*SCB)
      Srr = STRESS(1)*C2B + STRESS(2)*S2B + Q1
      Sss = STRESS(2)*C2B + STRESS(1)*S2B - Q1
      Srs = STRESS(3)*(C2B-S2B) + SCB*(STRESS(2)-STRESS(1))
!!
!! Put rotated stresses in temporary storage for use in MATERIAL_41_INTEGRATION.
!!
      LOCAL_STRESS(1) = Srr
      LOCAL_STRESS(2) = Sss
      LOCAL_STRESS(3) = 0.0
      LOCAL_STRESS(4) = Srs
      LOCAL_STRESS(5) = 0.0
      LOCAL_STRESS(6) = 0.0
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5D+0 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Rotate center of the yield surface from time n to time n+1. This is a
!! proper orthognal rotation using the rotation obtained from the polar
!! decomposition of the deformation gradient.
!!
      Q1 = (STATE_VARIABLES(3)*SCB) + (STATE_VARIABLES(3)*SCB)
      Srr = STATE_VARIABLES(1)*C2B + STATE_VARIABLES(2)*S2B + Q1
      Sss = STATE_VARIABLES(2)*C2B + STATE_VARIABLES(1)*S2B - Q1
      Srs = STATE_VARIABLES(3)*(C2B-S2B) +                                     &
     &  (STATE_VARIABLES(2)-STATE_VARIABLES(1))*SCB
!!
!! Put yield center in temporary storage for use in MATERIAL_41_INTEGRATION.
!!
      YIELD_SURFACE_CENTER(1) = Srr
      YIELD_SURFACE_CENTER(2) = Sss
      YIELD_SURFACE_CENTER(3) = 0.0
      YIELD_SURFACE_CENTER(4) = Srs
!!
!! Retrieve effective plastic strain.
!!
      EFFECTIVE_PLASTIC_STRAIN = STATE_VARIABLES(4)
!!
!! Update stress. Note that this membrane stress calculation uses the
!! plane stress shell module MATERIAL_41_INTEGRATION to insure identical
!! elastic-plastic behavior and to conserve on material model coding.
!!
      CALL MATERIAL_41_INTEGRATION                                             &
     &          (                                                              &
     &          LOCAL_STRESS,                                                  &
     &          YIELD_SURFACE_CENTER,                                          &
     &          EFFECTIVE_PLASTIC_STRAIN,                                      &
     &          DTnext,                                                        &
     &          Drr,Dss,Dtt,Drs,Drt,Dst,                                       &
     &          MatID,RCL2,RCS2                                                &
     &          )
!!
      SOUND_SPEED%RCL2 = RCL2
      SOUND_SPEED%RCS2 = RCS2
!!
!! Retrieve membrane stress components.
!!
      Srr = LOCAL_STRESS(1)
      Sss = LOCAL_STRESS(2)
      Srs = LOCAL_STRESS(4)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5D+0 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Return stresses to global storage.
!!
      STRESS(1) = Srr
      STRESS(2) = Sss
      STRESS(3) = Srs
!!
!! Return yield surface center to global storage.
!!
      STATE_VARIABLES(1) = YIELD_SURFACE_CENTER(1)
      STATE_VARIABLES(2) = YIELD_SURFACE_CENTER(2)
      STATE_VARIABLES(3) = YIELD_SURFACE_CENTER(4)
!!
!! Return effective plastic strain to global storage.
!!
      STATE_VARIABLES(4) = EFFECTIVE_PLASTIC_STRAIN
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_31                                                   &
     &          (STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 26-MAY-1991 13:27:50
!!
!! FINITE STRAIN ELASTIC-PLASTIC, LINEAR STRAIN HARDENING BEHAVIOR
!!
!! Module computes the current value of the stress STRESS(1:6),
!! the density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient), the past history
!! data STATE_VARIABLES(1:7), and the old stress state.
!!
!!      STATE_VARIABLES(1) = Yield surface center Sigma_XX
!!      STATE_VARIABLES(2) = Yield surface center Sigma_YY
!!      STATE_VARIABLES(3) = Yield surface center Sigma_ZZ
!!      STATE_VARIABLES(4) = Yield surface center Sigma_XY
!!      STATE_VARIABLES(5) = Yield surface center Sigma_XZ
!!      STATE_VARIABLES(6) = Yield surface center Sigma_YZ
!!      STATE_VARIABLES(7) = Accumulated effective plastic strain
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(6)
      REAL(KIND(0D0)), INTENT(INOUT) :: INTERNAL_ENERGY
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      INTEGER,         INTENT(IN)    :: MatID

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
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5D+0*DTnext) *                       &
     &          (                                                              &
     &           Dxx * STRESS(1) + Dyy * STRESS(2) + Dzz * STRESS(3)  +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) +        &
     &          (Dxy * STRESS(4) + Dxz * STRESS(5) + Dyz * STRESS(6)) -        &
     &          (Dxx+Dyy+Dzz) * INTERNAL_ENERGY                                &
     &          )
!!
!! Update stress.
!!
      CALL MATERIAL_31_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES(1),STATE_VARIABLES(7),                  &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,MatID               &
     &          )
!!
!! Internal energy increment from time n+1/2 to time n+1.
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
      SUBROUTINE MATERIAL_31_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,YLDC,EFPS,                                              &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,MatID               &
     &          )
!!
!! FINITE STRAIN, ELASTIC-PLASTIC, LINEAR STRAIN HARDENING BEHAVIOR.
!!
!! This model is based on a numerical procedure by R. D. Krieg for integrating
!! the incremental equations of plasticity.  The numerical method is based on
!! a generalization of the radial return method. material whose tensile
!! stress-strain curve has a constant plastic hardening modulus is described.
!! The hardening is a linear combination of isotropic and kinematic which is
!! controled by Beta.
!!
!!      R. P. Gëol and L .E. Malvern, Biaxial Plastic Simple Waves
!!      With Combined Kinematic And Isotropic Hardening, JOURNAL OF
!!      APPLIED MECHANICS, pages 1100-1106, (1970).
!!
!!      MATERIAL(MatID)%Ymod  = Young's modulus, (not used in calculation)
!!      MATERIAL(MatID)%Prat  = Poisson's ratio, QR
!!      MATERIAL(MatID)%Lmod  = Elastic Lame parameter Lambda, P1
!!      MATERIAL(MatID)%Gmod  = Elastic Lame parameter mu (shear modulus), P2
!!      MATERIAL(MatID)%Yield = Yield stress, QS
!!      MATERIAL(MatID)%H     = Hardening modulus, QH
!!      MATERIAL(MatID)%Beta  = Isotropic-kinematic hardening parameter, QB
!!                       0.0, Kinematic plastic hardening
!!                       1.0, Isotropic plastic hardening
!!                       Fraction, Combined isotropic-kinematic hardening
!!      MATERIAL(MatID)%p     = Strain rate exponent, QP
!!      MATERIAL(MatID)%D     = Strain rate scaling, QD
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(6)
      REAL(KIND(0D0)), INTENT(INOUT) :: YLDC(6)
      REAL(KIND(0D0)), INTENT(INOUT) :: EFPS
      REAL(KIND(0D0)), INTENT(IN)    :: DTnext
      REAL(KIND(0D0)), INTENT(IN)    :: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz
      INTEGER,         INTENT(IN)    :: MatID
!!
!! Local variables.
      REAL(KIND(0D0)), SAVE :: XI(6),Sinc(6),Einc(6)

      REAL(KIND(0D0)), PARAMETER :: Third = (1.0D+0 / 3.0D+0)
!!
!! Initialize constants:
!!      P1 = Elastic Lame parameter lambda
!!      P2 = Elastic Lame parameter mu (shear modulus)
!!      P3 = Twice the shear modulus, 2mu
!!
      Troot = SQRT (2.0D+0 / 3.0D+0)
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
      P3 = P2 + P2
!!
      QS = MATERIAL(MatID)%PVAL(10)
      QH = MATERIAL(MatID)%PVAL(12)
      QB = MATERIAL(MatID)%PVAL(13)
      QP = MATERIAL(MatID)%PVAL(16)
      QD = MATERIAL(MatID)%PVAL(17)
!!
!! AK equals the radius of yield surface (with hardening) times the
!! square root of 2/3.
!!
      AK = (QS + QB*QH*EFPS) * Troot
!!
!! Save old stress state.
!!
      DO i = 1,6
        Sinc(i) = STRESS(i)
      ENDDO
!!
!! Compute strain increments.
!!
      Einc(1) = DTnext * Dxx
      Einc(2) = DTnext * Dyy
      Einc(3) = DTnext * Dzz
      Einc(4) = DTnext * Dxy
      Einc(5) = DTnext * Dxz
      Einc(6) = DTnext * Dyz
!!
!! Compute trial elastic stresses.
!!
      DO i = 1,6
        STRESS(i) = STRESS(i) + P3 * Einc(i)
      ENDDO
      STRESS(1) = STRESS(1) + (DTnext * P1 * (Dxx+Dyy+Dzz))
      STRESS(2) = STRESS(2) + (DTnext * P1 * (Dxx+Dyy+Dzz))
      STRESS(3) = STRESS(3) + (DTnext * P1 * (Dxx+Dyy+Dzz))
!!
!! Test for zero yield stress, (Elastic response).
!!
      IF (AK .EQ. 0.0) GO TO 200
!!
!! Find distance from center of the yield surface to stress state.
!!
      DO I = 1,6
        XI(I) = STRESS(I) - YLDC(I)
      ENDDO
      XIM2 = Third *                                                           &
     &          (                                                              &
     &            (XI(1)-XI(2))*(XI(1)-XI(2))                                  &
     &          + (XI(1)-XI(3))*(XI(1)-XI(3))                                  &
     &          + (XI(2)-XI(3))*(XI(2)-XI(3))                                  &
     &          )                                                              &
     &          + (XI(4)*XI(4)) + (XI(4)*XI(4))                                &
     &          + (XI(5)*XI(5)) + (XI(5)*XI(5))                                &
     &          + (XI(6)*XI(6)) + (XI(6)*XI(6))
      XIM = SQRT (XIM2)
!!
!! Check the yield condition.
!!
      PHI = XIM2 - AK*AK
      IF (PHI .LE. 0.0) GO TO 200
!!
!! Compute dynamic yield stress.
!!
      IF (QP .NE. 0.0) THEN
        ED = Third *                                                           &
     &          (                                                              &
     &            (XI(1)-XI(2))*(Dxx-Dyy)                                      &
     &          + (XI(1)-XI(3))*(Dxx-Dzz)                                      &
     &          + (XI(2)-XI(3))*(Dyy-Dzz)                                      &
     &          )                                                              &
     &          + 4.0D+0 * (XI(4)*Dxy + XI(5)*Dxz + XI(6)*Dyz)
        IF (ED .GT. 0.0) THEN
          ED = ED / (XIM * (ONE + QH/(3.0D+0*P2)))
          ADD = QS*(ED/QD)**(ONE/QP)
          AK = AK + ADD*Troot
          PHI = XIM2 - AK*AK
          IF (PHI .LE. 0.0) GO TO 200
        ENDIF
      ENDIF
!!
!! Plastic loading condition. Calculate effective plastic strain (constant
!! plastic modulus).
!!
      C3 = 2.44948974278D+0*P2
      C4 = QH/(3.0D+0*P2)
      C5 = (ONE-AK/XIM)/(ONE+C4)
      EFPS = EFPS + C5*XIM/C3
      C6 = C5*C4*(ONE-QB)
!!
!! Find new values of stress and yield surface center.
!!
      XIA = Third*(XI(1)+XI(2)+XI(3))
      DO I = 1,3
        YLDC(I)   = YLDC(I)   + C6*(XI(I)-XIA)
        STRESS(I) = STRESS(I) - C5*(XI(I)-XIA)
      ENDDO
      DO I = 4,6
        YLDC(I)   = YLDC(I)   + C6*XI(I)
        STRESS(I) = STRESS(I) - C5*XI(I)
      ENDDO
!!
 200  CONTINUE
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
      SUBROUTINE MATERIAL_41                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,MatID,RCL2,RCS2                 &
     &          )
!!
!! Copyright (c) by KEY Associates; 21-MAY-1991 20:46:13
!!
!! FINITE STRAIN, PLANE-STRESS ELASTIC-PLASTIC STRAIN-HARDENING BEHAVIOR.
!!
!!      STATE_VARIABLES(1) = Yield surface center Stress_xx
!!      STATE_VARIABLES(2) = Yield surface center Stress_yy
!!      STATE_VARIABLES(3) = Yield surface center Stress_xy
!!      STATE_VARIABLES(4) = Effective plastic strain, EFPS
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(6)
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(4)
      REAL(KIND(0D0)), INTENT(IN   ) :: DTnext
      REAL(KIND(0D0)), INTENT(IN   ) :: Dxx,Dyy,Dxy,Dxz,Dyz,Wxy
      INTEGER,         INTENT(IN   ) :: MatID
      REAL(KIND(0D0)), INTENT(  OUT) :: RCL2,RCS2
!!
!! Local variables.
      REAL(KIND(0D0)) :: Dzz
      REAL(KIND(0D0)) :: YIELD_SURFACE_CENTER(4)
      REAL(KIND(0D0)) :: EFFECTIVE_PLASTIC_STRAIN
!!
      INTERFACE
        SUBROUTINE MATERIAL_41_INTEGRATION (STRESS,YLDC,EFPS,  &
     &    DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,MatID,RCL2,RCS2)

          REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(6)
          REAL(KIND(0D0)), INTENT(INOUT) :: YLDC(4)
          REAL(KIND(0D0)), INTENT(INOUT) :: EFPS
          REAL(KIND(0D0)), INTENT(IN   ) :: DTnext
          REAL(KIND(0D0)), INTENT(IN   ) :: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz
          INTEGER,         INTENT(IN   ) :: MatID
          REAL(KIND(0D0)), INTENT(  OUT) :: RCL2,RCS2

        END SUBROUTINE
      END INTERFACE
!!
!! Initialize material constants.
!!    P1 = Elastic Lame parameter, Lambda
!!    P2 = Elastic Lame parameter, Mu (shear modulus)
!!    P3 = Twice the shear modulus
!!
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
      P3 = P2 + P2
!!
!! Enforce plane stress condition by prescribing stretching in the normal
!! (or z) direction.  This development assumes that the initial value of
!! Dzz is zero.
!!
      Dzz = -P1*(Dxx+Dyy)/(P1+P3)
!!
!! "De-rotate" stress components to account for the fact that the
!! "co-rotational" axes have been redrawn from the last time step
!! along the same material fiber. The required rotational rate comes
!! from Vyx (=Dxy-Wxy). (The rotation based on Wxy contains the
!! finite rotations of the element relative to the co-rotational
!! coordinate axes, but does not account for the fact that the
!! co-rotational axes rotated with the element.)
!!
!! Rotation from time n to time n+1.
!!
      R1 = DTnext * (Wxy + (Dxy - Wxy))
      R2 = R1 + R1
      Sxx = STRESS(1) + (R2*STRESS(4))
      Syy = STRESS(2) - (R2*STRESS(4))
      Sxy = STRESS(4) + R1*(STRESS(2)-STRESS(1))
!!
!! Return stresses to global storage.
!!
      STRESS(1) = Sxx
      STRESS(2) = Syy
      STRESS(4) = Sxy
!!
!! Rotate center of the yield surface from time n to time n+1.
!!
      Sxx = STATE_VARIABLES(1) + (R2*STATE_VARIABLES(3))
      Syy = STATE_VARIABLES(2) - (R2*STATE_VARIABLES(3))
      Sxy = STATE_VARIABLES(3) + R1*(STATE_VARIABLES(2) - STATE_VARIABLES(1))
!!
!! Put yield center in temporary storage for use in MATERIAL_41_INTEGRATION.
!!
      YIELD_SURFACE_CENTER(1) = Sxx
      YIELD_SURFACE_CENTER(2) = Syy
      YIELD_SURFACE_CENTER(3) = 0.0
      YIELD_SURFACE_CENTER(4) = Sxy
!!
!! Retrieve effective plastic strain.
!!
      EFFECTIVE_PLASTIC_STRAIN = STATE_VARIABLES(4)
!!
!! Update stress.
!!
!!
      CALL MATERIAL_41_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,YIELD_SURFACE_CENTER,EFFECTIVE_PLASTIC_STRAIN,          &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,MatID,RCL2,RCS2                 &
     &          )
!!
!!
!! Return yield surface center to global storage.
!!
      STATE_VARIABLES(1) = YIELD_SURFACE_CENTER(1)
      STATE_VARIABLES(2) = YIELD_SURFACE_CENTER(2)
      STATE_VARIABLES(3) = YIELD_SURFACE_CENTER(4)
!!
!! Return effective plastic strain to global storage.
!!
      STATE_VARIABLES(4) = EFFECTIVE_PLASTIC_STRAIN
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_41_INTEGRATION ( STRESS,YLDC,EFPS,                   &
     &  DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,MatID,RCL2,RCS2 )
!!
!! Copyright (c) by KEY Associates; 25-APR-1991 20:30:30
!!
!! Plane Stress Elastic-Plastic Linear Strain Hardening Behavior.
!!
!! This module is based on a numerical procedure by Juan Simo for
!! integrating the incremental equations of plane stress plasticity.
!! The numerical method is a generalization of a method used in wave
!! codes. It is an implicit radial return algorithm using Newton's
!! Method to solve a nonlinear equation in the effective plastic strain
!! increment. A material whose tensile stress-strain curve has a
!! constant plastic hardening modulus is described.
!!
!!      J. C. Simo and S. Govindjee, Exact Closed-Form Solution
!!      of the Return Mapping Algorithm in Plane Stress Elasto-
!!      Viscoplasticity, ENGINEERING COMPUTATIONS, Vol. 5, No. 5,
!!      pages 254-258 (1988).
!!
!! A similar and somewhat earlier development by Ph. Jetteur can be
!! found in
!!
!!      Philippe Jetteur, Implicit Integration Algorithm for Elasto-
!!      Plasticity in Plane Stress, ENGINEERING COMPUTATIONS, Vol. 3,
!!      Number 3, pages 251-253, (1986).
!!
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!      *  NOTE: The plane-stress direction is along the z-axis.  *
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
!!      MATERIAL(MatID)%Ymod  = Young's modulus, (not used in calculation)
!!      MATERIAL(MatID)%Prat  = Poisson's ratio, QR
!!      MATERIAL(MatID)%Lmod  = Elastic Lame parameter Lambda, P1
!!      MATERIAL(MatID)%Gmod  = Elastic Lame parameter mu (shear modulus), P2
!!      MATERIAL(MatID)%Yield = Yield stress, QS
!!      MATERIAL(MatID)%H     = Hardening modulus, QH
!!      MATERIAL(MatID)%Beta  = Isotropic-kinematic hardening parameter, QB
!!      MATERIAL(MatID)%p     = Strain rate exponent, QP
!!      MATERIAL(MatID)%D     = Strain rate scaling, QD
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(INOUT) :: STRESS(6)
      REAL(KIND(0D0)), INTENT(INOUT) :: YLDC(4)
      REAL(KIND(0D0)), INTENT(INOUT) :: EFPS
      REAL(KIND(0D0)), INTENT(IN   ) :: DTnext
      REAL(KIND(0D0)), INTENT(IN   ) :: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz
      INTEGER,         INTENT(IN   ) :: MatID
      REAL(KIND(0D0)), INTENT(  OUT) :: RCL2,RCS2
!!
!! Local variables.
      INTEGER, SAVE   :: IFCONT = 0
      INTEGER, SAVE   :: IFPRNT = 1
      REAL(KIND(0D0)) :: Eta(4)
      REAL(KIND(0D0)) :: Einc(6)
      REAL(KIND(0D0)) :: Sinc(6)
      REAL(KIND(0D0)) :: P1,P2,P3
!!
      REAL(KIND(0D0)), PARAMETER :: ZERO = 0.0D+0
      REAL(KIND(0D0)), PARAMETER :: Epsilon = 1.0D-6
      INTEGER,         PARAMETER :: KMAX    = 8

      INTERFACE
        SUBROUTINE LOADING_MODULI (RCL2,RCS2,P1,P2,Sinc,Einc)

          REAL(KIND(0D0)), INTENT(OUT) :: RCL2
          REAL(KIND(0D0)), INTENT(OUT) :: RCS2
          REAL(KIND(0D0)), INTENT(IN)  :: P1
          REAL(KIND(0D0)), INTENT(IN)  :: P2
          REAL(KIND(0D0)), INTENT(IN)  :: Sinc(6)
          REAL(KIND(0D0)), INTENT(IN)  :: Einc(6)

        END SUBROUTINE
      END INTERFACE
!!
!! Initialize constants:
!!    P1 = Elastic Lame parameter lambda
!!    P2 = Elastic Lame parameter mu (shear modulus)
!!    P3 = Twice the shear modulus, 2mu
!!
      Third = (1.0D+0/3.0D+0)
      Troot = SQRT (2.0D+0/3.0D+0)
      Root2 = SQRT (2.0D+0)          ! Root2 = SQRT (2.0)
      Root6 = SQRT (6.0D+0)          ! Root6 = SQRT (6.0)
!!
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
      P3 = P2 + P2
!!
      QR = MATERIAL(MatID)%PVAL(7)
      QS = MATERIAL(MatID)%PVAL(10)
      QH = MATERIAL(MatID)%PVAL(12)
      QB = MATERIAL(MatID)%PVAL(13)
      QP = MATERIAL(MatID)%PVAL(16)
      QD = MATERIAL(MatID)%PVAL(17)
!!
!! AK equals the radius of yield surface (with hardening) times the
!! square root of 2/3.
!!
      AK = (QS + QB*QH*EFPS) * Troot
!!
!! Save old stress state. Note that the transverse shear stresses
!! Sxz and Syz are not used in the rho x c**2 calculation.
!!
      Sinc(1:6) = (/STRESS(1:4),ZERO,ZERO/)
!!
!! Compute strain increment.
!!
      Einc(1) = DTnext * Dxx
      Einc(2) = DTnext * Dyy
      Einc(3) = DTnext * Dzz
      Einc(4) = DTnext * Dxy
      Einc(5) = DTnext * Dxz
      Einc(6) = DTnext * Dyz
!!
!! Compute trial elastic stresses, including "transverse shears." Note
!! that this calculation assumes that Dzz has been pre-set to give Szz
!! a value of zero, the plane stress condition, ergo the actual formula
!! is commented-out in favor of the direct substitution of a zero value.
!!
      DO i = 1,6
        STRESS(i) = STRESS(i) + P3*Einc(i)
      ENDDO
      STRESS(1) = STRESS(1) + (DTnext * P1 * (Dxx+Dyy+Dzz))
      STRESS(2) = STRESS(2) + (DTnext * P1 * (Dxx+Dyy+Dzz))
!!    STRESS(3) = STRESS(3) + (DTnext * P1 * (Dxx+Dyy+Dzz))
      STRESS(3) = 0.0
!!
!! Eliminate the transverse shear strains from the rho x c**2 calculation.
!!
      Einc(5) = 0.0
      Einc(6) = 0.0
!!
!! Test for zero yield stress, (Elastic response).
!!
      IF (AK .EQ. 0.0) GO TO 200
!!
!! Compute magnitude of the deviatoric stress for use in yield calculation.
!!
      DO i = 1,4
        Eta(i) = STRESS(i) - YLDC(i)
      ENDDO
      Psi2 = Third*((Eta(1)-Eta(2)) * (Eta(1)-Eta(2))                        &
     &            + (Eta(1)*Eta(1)) + (Eta(2)*Eta(2)))                       &
     &            + (Eta(4)*Eta(4)) + (Eta(4)*Eta(4))
!!
!! Check the yield condition.
!!
      PHI = Psi2 - AK*AK
      IF (PHI .LE. 0.0) GO TO 200
!!
!! Compute dynamic yield stress.
!!
      IF (QP .NE. 0.0) THEN
        ED = 4.0D+0*(Eta(4)*Dxy) + Third *                                        &
     &    ((Eta(1)-Eta(2))*(Dxx-Dyy)+Eta(1)*(Dxx-Dzz)+Eta(2)*(Dyy-Dzz))
        IF (ED .GT. 0.0) THEN
          ED  = ED/(SQRT(Psi2)*(ONE+QH/(3.0D+0*P2)))
          ADD = QS*(ED/QD)**(ONE/QP)
          AK  = AK + ADD*Troot
          PHI = Psi2 - AK*AK
          IF (PHI .LE. 0.0) GO TO 200
        ENDIF
      ENDIF
!!
!! Plastic loading condition. Calculate discrete consistency parameter dL.
!! (Assumes a constant plastic modulus H and constant total strain rate.)
!!
      dL = 0.0
      KOUNT = 0
      A0 = (2.0D+0-QB-QB)*QH
      A1 = (Eta(1)+Eta(2))*(Eta(1)+Eta(2))
      A2 = Root6 * (P3*(ONE+QR)/(ONE-QR) + A0) * Third
      A3 = (Eta(1)-Eta(2))*(Eta(1)-Eta(2)) + 4.0D+0*Eta(4)*Eta(4)
      A4 = Root2 * (P3 + A0*Third)
      A5 = Troot*QB*QH
 100  CONTINUE
      B2 = Root6 + A2*dL
      C2 = B2 * B2
      B4 = Root2 + A4*dL
      C4 = B4 * B4
      F1 = A1/C2 + A3/C4
      A6 = SQRT (F1)
      A7 = AK + A5*(dL*A6)
      F2 = A7*A7
      IF ((F1-F2) .GT. Epsilon*F2) THEN
        KOUNT = KOUNT + 1
        DF1 = -2.0D+0 * ((A1/(B2*C2))*A2 + (A3/(B4*C4))*A4)
        DF2 =  2.0D+0 * A7*(A5*A6 + dL*A5*DF1/A6)
        dL = dL - (F1-F2)/(DF1-DF2)
        IF (KOUNT .LE. KMAX) GO TO 100
        IFCONT = IFCONT + 1
        IF (IFCONT .GE. IFPRNT) THEN
          IFPRNT = IFPRNT + IFPRNT
          WRITE (MSG1,'(I8)') KMAX
          WRITE (MSG2,'(I8)') IFCONT
          CALL USER_MESSAGE                                                    &
     &    (                                                                    &
     &    MSGL//'WARN'//                                                       &
     &    MSGL//'MATERIAL_41_INTEGRATION.001.00'//                             &
     &    MSGL//'Iteration For Consistency Parameter dL'//                     &
     &    MSGL//'Requires More Than'//MSG1//' Iterations To Converge.'//       &
     &    MSGL//'Number Of Failures To Converge:'//MSG2                        &
     &    )
        ENDIF
      ENDIF
      dEFPS = (dL*A6)
      EFPS = EFPS + dEFPS
      D2 = Root6 / B2
      D4 = Root2 / B4
      E1 = (Eta(1) + Eta(2)) * D2
      E2 = (Eta(1) - Eta(2)) * D4
      Eta(1) = 0.5D+0 * (E1+E2)
      Eta(2) = 0.5D+0 * (E1-E2)
      Eta(4) = Eta(4) * D4
      D5 = dL*Third*P3/(ONE-QR)
      STRESS(1) = STRESS(1) - D5*((2.0D+0-QR)*ETA(1)+(QR+QR-ONE)*ETA(2))
      STRESS(2) = STRESS(2) - D5*((2.0D+0-QR)*ETA(2)+(QR+QR-ONE)*ETA(1))
      STRESS(4) = STRESS(4) - D5*((3.0D+0-QR-QR-QR)*ETA(4))
      D5 = dL*Third*A0
      YLDC(1) = YLDC(1) + D5*ETA(1)
      YLDC(2) = YLDC(2) + D5*ETA(2)
      YLDC(4) = YLDC(4) + D5*ETA(4)
!!
 200  CONTINUE
!!
!! Sound speeds squared * rho(t)
!!
      IF (CONTROL%SUBCYC .EQ. 0) THEN
        DO i = 1,4
          Sinc(i) = STRESS(i) - Sinc(i)
        ENDDO
        CALL LOADING_MODULI ( RCL2,RCS2,P1,P2,Sinc,Einc )
        RCL2 = MIN (P1+P3, RCL2)
        RCS2 = MIN (P2   , RCS2)
      ELSE
        RCL2 = P1+P3
        RCS2 = P2
      ENDIF
!!
      RETURN
      END
