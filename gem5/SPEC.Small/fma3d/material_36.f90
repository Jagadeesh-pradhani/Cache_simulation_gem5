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
      SUBROUTINE MATERIAL_36_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:05:00
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
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_36_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_36                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:37:00
!!
!! FINITE STRAIN, LINEAR VISCOELASTIC MATERIAL MODEL.
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!!      STATE_VARIABLES( 1) = Integrated strain, xx-component
!!      STATE_VARIABLES( 2) = Integrated strain, yy-component
!!      STATE_VARIABLES( 3) = Integrated strain, zz-component
!!      STATE_VARIABLES( 4) = Integrated strain, xy-component
!!      STATE_VARIABLES( 5) = Integrated strain, xz-component
!!      STATE_VARIABLES( 6) = Integrated strain, yz-component
!!      STATE_VARIABLES( 7) = Viscoelastic state, xx-component
!!      STATE_VARIABLES( 8) = Viscoelastic state, yy-component
!!      STATE_VARIABLES( 9) = Viscoelastic state, zz-component
!!      STATE_VARIABLES(10) = Viscoelastic state, xy-component
!!      STATE_VARIABLES(11) = Viscoelastic state, xz-component
!!      STATE_VARIABLES(12) = Viscoelastic state, yz-component
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
!! Rotate the strain state STATE_VARIABLES(1:6).
!!
      IF (Ipolard .EQ. 0) CALL SYMMETRIC_TENSOR_ROTATION                       &
     &          (                                                              &
     &          STATE_VARIABLES(1),Igenrot,Ipolard,DTnext                      &
     &          )
!!
!! Rotate the viscoelastic state STATE_VARIABLES(7:12)
!!
      IF (Ipolard .EQ. 0) CALL SYMMETRIC_TENSOR_ROTATION                       &
     &          (                                                              &
     &          STATE_VARIABLES(7),Igenrot,Ipolard,DTnext                      &
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
!! Update strain, viscoelastic state, and stress.
!!
      CALL MATERIAL_36_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES(1),                                            &
     &          STATE_VARIABLES(7),                                            &
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
!! Rotate the strain state STATE_VARIABLES(1:6).
!!
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STATE_VARIABLES(1),Igenrot,Ipolard,DTnext                      &
     &          )
!!
!! Rotate the viscoelastic state STATE_VARIABLES(7:12)
!!
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &          (                                                              &
     &          STATE_VARIABLES(7),Igenrot,Ipolard,DTnext                      &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_36_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,STRAIN,VISCOUS,                                         &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:38:09
!!
!! LINEAR VISCOELASTIC MATERIAL MODEL.
!!
!! The numerical procedure for integrating the exponential relaxation
!! kernal comes from the following publication:
!!
!! L. R. Herrmann and F. E. Peterson, "A Numerical Procedure for Visco-
!! elastic Stress Analysis," Proceedings of the Seventh I.C.R.P.G.
!! Mechanical Behavior Working Group Meeting, C.P.I.A. Publication No.
!! 177, (October 1968).
!!
!!      Bulk constants:
!!              MATERIAL%Kzero, MATERIAL%Kinf, MATERIAL%Kbeta =
!!              K-zero, K-infinity, K-beta (exponential decay constant)
!!
!!      Deviatoric constants:
!!              MATERIAL%Gzero, MATERIAL%Ginf, MATERIAL%Gbeta =
!!              G-zero, G-infinity, G-beta (exponential decay constant)
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          STRESS(6),                                                     &
     &          STRAIN(6),                                                     &
     &          VISCOUS(6),                                                    &
     &          Sinc(6),Einc(6),                                               &
     &          D(6),Kz,Ki,Kd,Sz,Si,Sd
!!
!! Move the stretching Dmn (the symmetric part of the velocity gradient)
!! into the indexed array D(1:6) to facilitate the use of do-loops.
!!
      D(1) = Dxx
      D(2) = Dyy
      D(3) = Dzz
      D(4) = Dxy
      D(5) = Dxz
      D(6) = Dyz
!!
      DO i = 1,6
!!
!! Initialize Sinc with stress at time n.
!!
        Sinc(i) = STRESS(i)
!!
!! Compute strain increment Einc.
!!
        Einc(i) = (DTnext*D(i))
!!
!! Update total strain STRAIN. Note: For finite strains, the total strain
!! STRAIN can be obtained from the integral of the stretching Dmn only if
!! the total strain is simultaneously rotated as a part of the integration.
!! Here, the integration has been "split" into two sequential steps. First,
!! the total strain is rotated from time(n) to time(n+1). (The rotation step
!! is implemented in the companion calling module.) Second, the incremental
!! strain from time(n) to time(n+1) is added.
!!
        STRAIN(i) = STRAIN(i) + (DTnext*D(i))
      ENDDO
!!
!! Separate strain STRAIN, viscoelastic state VISCOUS, and stretching D into
!! bulk and deviatoric components.
!!
      Ekk  = STRAIN(1)  + STRAIN(2)  + STRAIN(3)
      Vkk  = VISCOUS(1) + VISCOUS(2) + VISCOUS(3)
      Dkk  = D(1)       + D(2)       + D(3)
      DO i = 1,3
        STRAIN(i)  = STRAIN(i)  - 0.33333333333*Ekk
        VISCOUS(i) = VISCOUS(i) - 0.33333333333*Vkk
        D(i)       = D(i)       - 0.33333333333*Dkk
      ENDDO
!!
!! Bulk response. Retrieve K-zero, K-infinity, and decay constant Beta.
!!
      Kz = MATERIAL(MatID)%PVAL(6)
      Ki = MATERIAL(MatID)%PVAL(7)
      Kd = Kz - Ki
      Beta = MATERIAL(MatID)%PVAL(8)
      IF (Beta .EQ. 0.0 .OR. Kd .EQ. 0.0) THEN
!!
!! Elastic bulk response.
!!
        Pressure  = Kz * Ekk
!!
      ELSE
!!
!! Viscoelastic bulk response (Pressure is positive in tension).
!!
        Q1 = EXP (-Beta*DTnext)
        Q2 = (1.0 - Q1) / Beta
        Vkk = (Q1*Vkk) + (Q2*Dkk)
        Pressure = (Ki*Ekk) + (Kd*Vkk)
!!
      ENDIF
!!
!! Shear response. Retrieve S-zero, S-infinity, and decay constant Beta.
!!
      Sz = MATERIAL(MatID)%PVAL(9)
      Si = MATERIAL(MatID)%PVAL(10)
      Sd = Sz - Si
      Beta = MATERIAL(MatID)%PVAL(11)
      IF (Beta .EQ. 0.0 .OR. Sd .EQ. 0.0) THEN
!!
!! Elastic shear response.
!!
        DO i = 1,6
          STRESS(i) = (Sz+Sz) * STRAIN(i)
        ENDDO
!!
      ELSE
!!
!! Viscoelastic shear response.
!!
        Q1 = EXP (-Beta*DTnext)
        Q2 = (1.0 - Q1) / Beta
        DO i = 1,6
          VISCOUS(i) = Q1*VISCOUS(i) + Q2*D(i)
        ENDDO
        DO i = 1,6
          STRESS(i) = (Si+Si)*STRAIN(i) + (Sd+Sd)*VISCOUS(i)
        ENDDO
!!
      ENDIF
!!
!! Combine bulk and deviatoric components of Strain STRAIN, Viscoelastic
!! state VISCOUS, and Stress STRESS.
!!
      DO i = 1,3
        STRAIN(i)  = STRAIN(i)  + 0.33333333333*Ekk
        VISCOUS(i) = VISCOUS(i) + 0.33333333333*Vkk
        STRESS(i)  = STRESS(i)  + Pressure
      ENDDO
!!
!! Estimate Lamé parameters for use in the event the strain
!! increment is zero.
!!
      P1 = Kz + 1.33333333333*Sz
      P2 = Sz
      P3 = P2 + P2
!!
!! Compute density times the longitudinal sound speed CL squared,
!! CL2, and density times the shear sound speed CS squared, CS2,
!! in the current state from stress and strain increments.
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
