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
      SUBROUTINE MATERIAL_10_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:52
!!
      USE shared_common_data
      USE material_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER  SecID
      REAL(KIND(0D0))      STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_10_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_20_INIT (MatID)
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
      REAL(KIND(0D0))                                                          &
     &          Lmod,                                                          &
     &          STRESS(3)
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
      RETURN
!!
      ENTRY MATERIAL_20_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_30_INIT (MatID)
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
      REAL(KIND(0D0))                                                          &
     &          Lmod,                                                          &
     &          STRESS(6)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lamé parameters lambda (Lmod) and mu (Gmod).
!!
      write (IO_UNIT%LELO,*) ' Into MATER_30_INIT...,MatID:',MatID
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((1.0 + Prat) * (1.0 - 2.0*Prat))
      Gmod = Ymod / (2.0 + 2.0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
      RETURN
!!
      ENTRY MATERIAL_30_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_40_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:57
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
     &          Lmod,                                                          &
     &          STRESS(6,Ipts)
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
      RETURN
!!
      ENTRY MATERIAL_40_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_50_INIT (MatID)
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
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          MatID,                                                         &
     &          SecID
      REAL(KIND(0D0))                                                          &
     &          Lmod,                                                          &
     &          STRESS(6,Ipts)
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
      RETURN
!!
      ENTRY MATERIAL_50_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_10                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 16:14:25
!!
!! FINITE STRAIN UNIAXIAL-STRESS ELASTIC BEHAVIOR
!!
!! This subroutine computes the current value of the stress STRESS(1:1),
!! the longitudinal sound speed CL in the current state as a function of
!! the stretching Drr and the old stress state.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
      COMMON /TRUSS/                                                           &
     &          Rx,Ry,Rz,Drr
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*STRESS
!!
!! Update stress.
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      STRESS = STRESS + Ymod * (DTnext*Drr)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY + (0.5*(DTnext*Drr))*STRESS
!!
!! Sound speed squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Ymod
      SOUND_SPEED%RCS2 = 0.0
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_20                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 20-APR-1991 18:21:59
!!
!! FINITE STRAIN PLANE-STRESS ELASTIC BEHAVIOR
!!
!! This subroutine computes the current value of the stress STRESS(1:3), the
!! longitudinal sound speed CL, and the shear sound speed CS in the current
!! state, as a function of the stretching Dxx,Dyy,Dxy and the old stress
!! state.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it maybe found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol. 32, pp 217-232, (1979).
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
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(3),                                                     &
     &          INTERNAL_ENERGY,                                               &
     &          STATE_VARIABLES(*)
!!
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
!! Dzz is zero.
!!
      Dtt = -P1*(Drr+Dss)/(P1+P3)
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
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Compute elastic stress increments and add to old stress.
!!
      QA = (DTnext * P1 * (Drr+Dss+Dtt))
!!
      Srr = Srr + (DTnext * P3) * Drr + QA
      Sss = Sss + (DTnext * P3) * Dss + QA
      Srs = Srs + (DTnext * P3) * Drs
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Return new stresses to global storage.
!!
      STRESS(1) = Srr
      STRESS(2) = Sss
      STRESS(3) = Srs
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = P1 + P3
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_30                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 24-FEB-1991 12:22:38
!!
!! FINITE STRAIN HYPO-ELASTIC ISOTROPIC MATERIAL MODEL.
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!!      STATE_VARIABLES(-) = (not used)
!!
!! Note: This material model does not require any state variables. The
!! array STATE_VARIABLES is present to maintain the standard material
!! interface format.
!!
!! This approach to elasicity is based on hypo-elastic concepts. The theory
!! and further references to it may be found in:
!!
!!      J.K. Dienes, "On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies", ACTA MECHANICA, Vol.32, pp 217-232, (1979).
!!
      USE shared_common_data
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
      CALL MATERIAL_30_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,MatID,                                                  &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz                     &
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
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_30_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,MatID,                                                  &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz                     &
     &          )
!!
!! Copyright (c) by KEY Associates; 24-FEB-1991 20:40:14
!!
!! ISOTROPIC ELASTIC MATERIAL MODEL.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6)
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
!! Compute elastic stress increments and add to old stress.
!!
      QA = (DTnext * P1 * (Dxx+Dyy+Dzz))
!!
      STRESS(1) = STRESS(1) + (DTnext * P3) * Dxx + QA
      STRESS(2) = STRESS(2) + (DTnext * P3) * Dyy + QA
      STRESS(3) = STRESS(3) + (DTnext * P3) * Dzz + QA
      STRESS(4) = STRESS(4) + (DTnext * P3) * Dxy
      STRESS(5) = STRESS(5) + (DTnext * P3) * Dxz
      STRESS(6) = STRESS(6) + (DTnext * P3) * Dyz
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = P1 + P3
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_40                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,MatID                           &
     &          )
!!
!! Copyright (c) by KEY Associates; 26-MAR-1991 21:21:25
!!
!! FINITE STRAIN HYPO-ELASTIC ISOTROPIC MATERIAL MODEL.
!!
!! Module computes the current value of the stress STRESS(1:6,*),
!! the density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!! STATE_VARIABLES(1) = r-axis rotation rate from preceeding time step
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
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          STATE_VARIABLES(*)
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
      Sxz = STRESS(5)
      Syz = STRESS(6)
!!
!! Compute elastic stress increments and add to old stress.
!!
      QA = (DTnext * P1 * (Dxx+Dyy+Dzz))
!!
      Sxx = Sxx + (DTnext * P3) * Dxx + QA
      Syy = Syy + (DTnext * P3) * Dyy + QA
      Sxy = Sxy + (DTnext * P3) * Dxy
      Sxz = Sxz + (DTnext * P3) * Dxz
      Syz = Syz + (DTnext * P3) * Dyz
!!
!! Return new stresses to global storage.
!!
      STRESS(1) = Sxx
      STRESS(2) = Syy
      STRESS(4) = Sxy
      STRESS(5) = Sxz
      STRESS(6) = Syz
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = P1 + P3
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
