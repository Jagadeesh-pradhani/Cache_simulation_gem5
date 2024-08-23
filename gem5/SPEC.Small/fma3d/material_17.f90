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
      SUBROUTINE MATERIAL_17_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:26
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_17_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_27_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:28
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_27_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_37_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:29
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID,                                                         &
     &          LupID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_37_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_47_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:29
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_47_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_57_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 3-AUG-1993 19:22:33.31
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!!
!!      estimate Ymod (PVAL(9))
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_57_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_17                                                   &
     &          (STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:31
!!
!! FINITE STRAIN UNIAXIAL-STRESS NONLINEAR ELASTIC (RUBBER) BEHAVIOR
!!
!! This subroutine computes the current value of the stress STRESS(1:1),
!! the longitudinal sound speed CL in the current state as a function of
!! the stretching Drr and the old stress state.
!!
!!      STATE_VARIABLES(1) = Err, Axial Strain
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*),                                            &
     &          TABLE_LOOK_UP
!!
      COMMON /TRUSS/                                                           &
     &          Rx,Ry,Rz,Drr
!!
      DATA                                                                     &
     &          OneThird /0.33333333333333333/                                 &
     &          TwoThird /0.66666666666666667/
!!
!! Initialize material constants.
!!      P1 = Elastic Lame parameter, Lambda = k - 2µ/3
!!      P2 = Elastic Lame parameter, Mu (shear modulus)
!!      P3 = Twice the shear modulus
!!
      Bulk = MATERIAL(MatID)%PVAL(10)
      HstID = NINT (MATERIAL(MatID)%PVAL(9))
!!
!! Save old effective stress and strain.
!!
      Old_Effective_Stress = ABS (STRESS)
      Old_Effective_Strain = ABS (STATE_VARIABLES(1))
!!
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*STRESS
!!
!! Update axial strain.
!!
      Err = STATE_VARIABLES(1) + (DTnext*Drr)
!!
!! Return new strain to global storage.
!!
      STATE_VARIABLES(1) = Err
!!
!! Compute new effective strain.
!!
      Effective_Strain = ABS (Err)
!!
!! Look up new effective stress.
!!
      Effective_Stress =                                                       &
     &          TABLE_LOOK_UP (HstID,Effective_Strain)
!!
!! Compute new elastic stress.
!!
      STRESS = SIGN (Effective_Stress,Err)
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*STRESS
!!
!! Compute tangent shear modulus.
!!
      IF (ABS(Effective_Strain-Old_Effective_Strain) .GT. 1.0E-20) THEN
        P2 = OneThird * (Effective_Stress - Old_Effective_Stress)              &
     &                / (Effective_Strain - Old_Effective_Strain)
      ELSE
        P2 = 1.0E-20
      ENDIF
!!
!! Sound speed squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Bulk + TwoThird * (P2+P2)
      SOUND_SPEED%RCS2 = 0.0
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_27                                                   &
     &          (STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:32
!!
!! FINITE STRAIN PLANE-STRESS NONLINEAR ELASTIC (RUBBER) BEHAVIOR
!!
!! This subroutine computes the current value of the stress STRESS(1:4), the
!! longitudinal sound speed CL, and the shear sound speed CS in the current
!! state, as a function of the stretching Drr,Dss,Drs and the old stress
!! state.
!!
!!      STATE_VARIABLES(1) = Err, rr-component, total strain
!!      STATE_VARIABLES(2) = Ess, ss-component, total strain
!!      STATE_VARIABLES(3) = Ett, tt-component, total strain
!!      STATE_VARIABLES(4) = Ers, rs-component, total strain
!!      STATE_VARIABLES(5) = E-bar, effective strain, 2(E'ijE'ij)/3
!!      STATE_VARIABLES(6) = S-bar, effective stress, 3(S'ijS'ij)/2
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it maybe found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol. 32, pp 217-232, (1979).
!!
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!      *  NOTE: The plane-stress direction is along the t-axis.  *
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(3),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*),                                            &
     &          TABLE_LOOK_UP
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching, spin
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Gr,Gs,                  & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
      DATA                                                                     &
     &          OneThird /0.33333333333333333D0/                               &
     &          TwoThird /0.66666666666666667D0/
!!
!! Initialize material constants.
!!      P1 = Elastic Lame parameter, Lambda = k - 2µ/3
!!      P2 = Elastic Lame parameter, Mu (shear modulus)
!!      P3 = Twice the shear modulus
!!
      Bulk = MATERIAL(MatID)%PVAL(10)
      HstID = NINT (MATERIAL(MatID)%PVAL(9))
!!
!! Enforce plane stress condition, or in the case of rubber, incompressiblity
!! by prescribing stretching in the normal (or t) direction. This development
!! assumes that the initial value of Dtt is zero.
!!
      Dtt = -(Drr + Dss)
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
!! Internal energy increment from time n to time n+1/2.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Rotate the strain from time n to time n+1. This is a proper orthognal
!! rotation using the rotation obtained from the polar decomposition of
!! the deformation gradient.
!!
      Q1 = (STATE_VARIABLES(4)*SCB) + (STATE_VARIABLES(4)*SCB)
      Err = STATE_VARIABLES(1)*C2B + STATE_VARIABLES(2)*S2B + Q1
      Ess = STATE_VARIABLES(2)*C2B + STATE_VARIABLES(1)*S2B - Q1
      Ett = STATE_VARIABLES(3)
      Ers = STATE_VARIABLES(4)*(C2B-S2B) +                                     &
     &  SCB*(STATE_VARIABLES(2)-STATE_VARIABLES(1))
!!
!! Update strain.
!!
      Err = Err + DTnext * Drr
      Ess = Ess + DTnext * Dss
      Ett = Ett + DTnext * Dtt
      Ers = Ers + DTnext * Drs
!!
!! Return new strain to global storage.
!!
      STATE_VARIABLES(1) = Err
      STATE_VARIABLES(2) = Ess
      STATE_VARIABLES(3) = Ett
      STATE_VARIABLES(4) = Ers
!!
!! Compute new effective strain. (This strain has no volumetric component
!! and therefore is only a deviatoric strain.)
!!
      Effective_Strain = SQRT                                                  &
     &  (TwoThird*(Err*Err+Ess*Ess+Ett*Ett+(Ers*Ers)+(Ers*Ers)))
!!
!! Look up new effective stress.
!!
      Effective_Stress =                                                       &
     &          TABLE_LOOK_UP (HstID,Effective_Strain)
!!
!! Compute new elastic stress.
!!
      IF (Effective_Strain .GT. 1.0E-20) THEN
        Alpha = TwoThird * Effective_Stress / Effective_Strain
      ELSE
        Alpha = 0.0
      ENDIF
!!
      Srr = Alpha * ((Err + Ess) + Err)
      Sss = Alpha * ((Err + Ess) + Ess)
      Srs = Alpha * Ers
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Return new stress to global storage.
!!
      STRESS(1) = Srr
      STRESS(2) = Sss
      STRESS(3) = Srs
!!
!! Retrieve old effective stress and strain.
!!
      Old_Effective_Strain = STATE_VARIABLES(5)
      Old_Effective_Stress = STATE_VARIABLES(6)
!!
!! Update effective stress and strain.
!!
      STATE_VARIABLES(5) = Effective_Strain
      STATE_VARIABLES(6) = Effective_Stress
!!
!! Compute tangent shear modulus.
!!
      IF (ABS(Effective_Strain-Old_Effective_Strain) .GT. 1.0E-20) THEN
        P2 = OneThird * (Effective_Stress - Old_Effective_Stress)              &
     &                / (Effective_Strain - Old_Effective_Strain)
      ELSE
        P2 = 0.0
      ENDIF
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Bulk + TwoThird * (P2+P2)
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_37                                                   &
     &          (STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID)
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:33
!!
!! FINITE STRAIN HYPO-ELASTIC NONLINEAR ELASTIC (RUBBER) BEHAVIOR
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!!      STATE_VARIABLES(1) = Exx, xx-component, total strain
!!      STATE_VARIABLES(2) = Eyy, yy-component, total strain
!!      STATE_VARIABLES(3) = Ezz, zz-component, total strain
!!      STATE_VARIABLES(4) = Exy, xy-component, total strain
!!      STATE_VARIABLES(5) = Exz, xz-component, total strain
!!      STATE_VARIABLES(6) = Eyz, yz-component, total strain
!!      STATE_VARIABLES(7) = E-bar, effective strain, 2(E'ijE'ij)/3
!!      STATE_VARIABLES(8) = S-bar, effective stress, 3(S'ijS'ij)/2
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
      IF (Ipolard .EQ. 0) THEN
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STRESS,Igenrot,Ipolard,DTnext                                  &
     &          )
!!
!! Rotate the strain state STATE_VARIABLES(1:6).
!!
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STATE_VARIABLES(1),Igenrot,Ipolard,DTnext                      &
     &          )
      ENDIF
!!
      Exx = STATE_VARIABLES(1)
      Eyy = STATE_VARIABLES(2)
      Ezz = STATE_VARIABLES(3)
      Exy = STATE_VARIABLES(4)
      Exz = STATE_VARIABLES(5)
      Eyz = STATE_VARIABLES(6)
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
!! Update strain and stress.
!!
      CALL MATERIAL_37_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          Exx,Eyy,Ezz,Exy,Exz,Eyz,STATE_VARIABLES(7),                    &
     &          STATE_VARIABLES(8),MatID                                       &
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
      STATE_VARIABLES(1) = Exx
      STATE_VARIABLES(2) = Eyy
      STATE_VARIABLES(3) = Ezz
      STATE_VARIABLES(4) = Exy
      STATE_VARIABLES(5) = Exz
      STATE_VARIABLES(6) = Eyz
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
!! Rotate the strain state STATE_VARIABLES(1:6).
!!
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STATE_VARIABLES(1),Igenrot,Ipolard,DTnext                      &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_37_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,                                                        &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          Exx,Eyy,Ezz,Exy,Exz,Eyz,Effective_Strain,                      &
     &          Effective_Stress,MatID                                         &
     &          )
!!
!! Copyright (c) by KEY Associates; 2-MAR-1992 20:53:34
!!
!! ISOTROPIC NONLINEAR ELASTIC (RUBBER) BEHAVIOR
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          TABLE_LOOK_UP
      DATA                                                                     &
     &          OneThird /0.33333333333333333/                                 &
     &          TwoThird /0.66666666666666667/
!!
!! Initialize material constants.
!!      P1 = Elastic Lame parameter, Lambda = k - 2µ/3
!!      P2 = Elastic Lame parameter, Mu (shear modulus)
!!      P3 = Twice the shear modulus
!!
      Bulk = MATERIAL(MatID)%PVAL(10)
      HstID = NINT (MATERIAL(MatID)%PVAL(9))
!!
!! Save old effective stress and strain.
!!
      Old_Effective_Stress = Effective_Stress
      Old_Effective_Strain = Effective_Strain
!!
!! Update strain.
!!
      Exx = Exx + DTnext * Dxx
      Eyy = Eyy + DTnext * Dyy
      Ezz = Ezz + DTnext * Dzz
      Exy = Exy + DTnext * Dxy
      Exz = Exz + DTnext * Dxz
      Eyz = Eyz + DTnext * Dyz
!!
!! Compute new effective strain.
!!
      Hxx = Exx - (OneThird * (Exx + Eyy + Ezz))
      Hyy = Eyy - (OneThird * (Exx + Eyy + Ezz))
      Hzz = Ezz - (OneThird * (Exx + Eyy + Ezz))
      Effective_Strain = SQRT                                                  &
     &          (                                                              &
     &          TwoThird * (Hxx*Hxx + Hyy*Hyy + Hzz*Hzz +                      &
     &          (Exy*Exy) + (Exy*Exy) + (Exz*Exz) + (Exz*Exz) +                &
     &          (Eyz*Eyz) + (Eyz*Eyz))                                         &
     &          )
!!
!! Look up new effective stress.
!!
      Effective_Stress =                                                       &
     &          TABLE_LOOK_UP (HstID,Effective_Strain)
!!
!! Compute new elastic stress.
!!
      IF (Effective_Strain .GT. 1.0E-20) THEN
        Alpha = TwoThird * Effective_Stress / Effective_Strain
      ELSE
        Alpha = 0.0
      ENDIF
!!
      STRESS(1) = Alpha * Hxx + (Bulk * (Exx + Eyy + Ezz))
      STRESS(2) = Alpha * Hyy + (Bulk * (Exx + Eyy + Ezz))
      STRESS(3) = Alpha * Hzz + (Bulk * (Exx + Eyy + Ezz))
      STRESS(4) = Alpha * Exy
      STRESS(5) = Alpha * Exz
      STRESS(6) = Alpha * Eyz
!!
!! Compute tangent shear modulus.
!!
      IF (ABS(Effective_Strain-Old_Effective_Strain) .GT. 1.0E-20) THEN
        P2 = OneThird * (Effective_Stress - Old_Effective_Stress)              &
     &                / (Effective_Strain - Old_Effective_Strain)
      ELSE
        P2 = 0.0
      ENDIF
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Bulk + TwoThird * (P2+P2)
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_47                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,MatID                           &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-JUN-1992 13:07:03.56
!!
!! FINITE STRAIN PLANE-STRESS NONLINEAR ELASTIC (RUBBER) BEHAVIOR
!!
!! This subroutine computes the current value of the stress STRESS(1:6), the
!! longitudinal sound speed CL, and the shear sound speed CS in the current
!! state, as a function of the stretching Dxx,Dyy,Dxy,Dxz,Dyz, and the old
!! strain state.
!!
!!      STATE_VARIABLES(1) = Exx, xx-component, total strain
!!      STATE_VARIABLES(2) = Eyy, yy-component, total strain
!!      STATE_VARIABLES(3) = Ezz, zz-component, total strain
!!      STATE_VARIABLES(4) = Exy, xy-component, total strain
!!      STATE_VARIABLES(5) = Exz, xz-component, total strain
!!      STATE_VARIABLES(6) = Eyz, yz-component, total strain
!!      STATE_VARIABLES(7) = E-bar, effective strain, 2(E'ijE'ij)/3
!!      STATE_VARIABLES(8) = S-bar, effective stress, 3(S'ijS'ij)/2
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
!! Note: Until such time that a computation of the polar decomposition
!! of the deformation gradient, F = VR = RU is implemented to obtain
!! the orthogonal rotation R, the Jaumann stress flux based on the spin
!! W will be used. The Jaumann rate is good for shear strains up to 40%.
!!
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!      *  NOTE: The plane-stress direction is along the z-axis.  *
!!      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          STATE_VARIABLES(*),                                            &
     &          TABLE_LOOK_UP
!!
      REAL(KIND(0D0)), PARAMETER :: OneThird = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: TwoThird = (2.0D+0 / 3.0D+0)
!!
!! Initialize material constants.
!!      P1 = Elastic Lame parameter, Lambda = k - 2Mu/3
!!      P2 = Elastic Lame parameter, Mu (shear modulus)
!!      P3 = Twice the shear modulus
!!
      Bulk = MATERIAL(MatID)%PVAL(10)
      HstID = NINT (MATERIAL(MatID)%PVAL(9))
!!
!! Enforce plane stress condition, or in the case of rubber, incompressiblity
!! by prescribing stretching in the normal (or z) direction. This development
!! assumes that the initial value of Dzz is zero.
!!
      Dzz = -(Dxx + Dyy)
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
      Exx = STATE_VARIABLES(1) + (R2*STATE_VARIABLES(4))
      Eyy = STATE_VARIABLES(2) - (R2*STATE_VARIABLES(4))
      Ezz = STATE_VARIABLES(3)
      Exy = STATE_VARIABLES(4) +                                               &
     &  R1*(STATE_VARIABLES(2)-STATE_VARIABLES(1))
      Exz = STATE_VARIABLES(5)
      Eyz = STATE_VARIABLES(6)
!!
!! Update strain.
!!
      Exx = Exx + DTnext * Dxx
      Eyy = Eyy + DTnext * Dyy
      Ezz = Ezz + DTnext * Dzz
      Exy = Exy + DTnext * Dxy
      Exz = Exz + DTnext * Dxz
      Eyz = Eyz + DTnext * Dyz
!!
!! Return new strain to global storage.
!!
      STATE_VARIABLES(1) = Exx
      STATE_VARIABLES(2) = Eyy
      STATE_VARIABLES(3) = Ezz
      STATE_VARIABLES(4) = Exy
      STATE_VARIABLES(5) = Exz
      STATE_VARIABLES(6) = Eyz
!!
!! Compute new effective strain. (This strain has no volumetric component
!! and therefore is only a deviatoric strain.)
!!
      Effective_Strain = SQRT                                                  &
     &  (TwoThird * (Exx*Exx+Eyy*Eyy+Ezz*Ezz+(Exy*Exy)+(Exy*Exy)))
!!
!! Look up new effective stress.
!!
      Effective_Stress =                                                       &
     &          TABLE_LOOK_UP (HstID,Effective_Strain)
!!
!! Compute new elastic stress.
!!
      IF (Effective_Strain .GT. 1.0E-20) THEN
        Alpha = TwoThird * Effective_Stress / Effective_Strain
      ELSE
        Alpha = 0.0
      ENDIF
!!
      Sxx = Alpha * ((Exx + Eyy) + Exx)
      Syy = Alpha * ((Exx + Eyy) + Eyy)
      Szz = 0.0
      Sxy = Alpha * Exy
      Sxz = Alpha * Exz
      Syz = Alpha * Eyz
!!
!! Return new stresses to global storage.
!!
      STRESS(1) = Sxx
      STRESS(2) = Syy
      STRESS(4) = Sxy
      STRESS(5) = Sxz
      STRESS(6) = Syz
!!
!! Retrieve old effective stress and strain.
!!
      Old_Effective_Strain = STATE_VARIABLES(7)
      Old_Effective_Stress = STATE_VARIABLES(8)
!!
!! Update effective stress and strain.
!!
      STATE_VARIABLES(7) = Effective_Strain
      STATE_VARIABLES(8) = Effective_Stress
!!
!! Compute tangent shear modulus P2.
!!
      IF (ABS(Effective_Strain-Old_Effective_Strain) .GT. 1.0E-20) THEN
        P2 = OneThird * (Effective_Stress - Old_Effective_Stress)              &
     &                / (Effective_Strain - Old_Effective_Strain)
      ELSE
        P2 = 0.0
      ENDIF
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Bulk + TwoThird * (P2+P2)
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
