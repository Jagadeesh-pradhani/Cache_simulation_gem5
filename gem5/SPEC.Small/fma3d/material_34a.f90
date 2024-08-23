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
      SUBROUTINE MATERIAL_34_INIT (MatID)
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
     &          RL(4),                                                         &
     &          STRESS(6)
!!
!! note: we must set number of state vars to 32. a lot of work
!! remains to be done defining input properties and handing them
!! off to the jointed rock model. the old method used EE.
!!
!! Temporary !!! until MATERIAL structure is set-up to carry material
!! properties
!!
      REAL(KIND(0D0))                                                          &
     &          EE(13,7)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Compute the Lame parameters lambda (Lmod) and mu (Gmod).
!!
      Ymod = MATERIAL(MatID)%PVAL(6)
      Prat = MATERIAL(MatID)%PVAL(7)

      Lmod = Prat * Ymod / ((1.0 + Prat) * (1.0 - 2.0*Prat))
      Gmod = Ymod / (2.0 + 2.0*Prat)

      MATERIAL(MatID)%PVAL(8) = Lmod
      MATERIAL(MatID)%PVAL(9) = Gmod
!!
!! Check the input data for valid ranges and supply default data for
!! unspecified input data.
!!
      NUMJT = 0  !  NINT(EE(1,3))
      DO I = 1,4
        RL(I) = 0.0  !  EE(I+1,3)
      ENDDO
!!
!!    CALL JNTDCK (NUMJT,RL,EE(1,4))
!!
!! Return block sliding length.
!!
      DO I = 1,4
        EE(I+1,3) = RL(I)
      ENDDO
!!
      RETURN
!!
      ENTRY MATERIAL_34_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! Initialize unit normals to joint planes, test initial stresses for
!! compatability with joint shear capabilty, and initialize joint gap
!! displacements.
!!
!! Define 3-D stress components and characteristic joint sliding length.
!!
      NUMJT = 0  !  NINT(EE(1,3))
      DO I = 1,4
        RL(I) = 0.0  !  EE(I+1,3)
      ENDDO
      Idim = 3
!!
!!    CALL JNTINI (NUMJT,RL,EE(1,4),STRESS,STATE_VARIABLES(Isv),Idim)
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_34                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 26-MAY-1991 13:27:50
!!
!! FINITE STRAIN ELASTIC MATRIX - DEFORMABLE JOINT BEHAVIOR.
!!
!! Module computes the current value of the stress STRESS(1:6),
!! the density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient), the past history
!! data STATE_VARIABLES(1:32) (= SVS), and the old stress state.
!!
!!   For i = 1,4
!!   SVS(8*i-7) = Ax, x-component of unit vector normal to joint surface.
!!   SVS(8*i-6) = Ay, y-component of unit vector normal to joint surface.
!!   SVS(8*i-5) = Az, z-component of unit vector normal to joint surface.
!!   SVS(8*i-4) = Accumulated shear displacement -- (Model = 1,2,3).
!!   SVS(8*i-3) = Accumulated normal displacement - (Model = 1,2,3).
!!   SVS(8*i-2) = Accumulated dilatation ---------- (Model = 1, , ).
!!   SVS(8*i-1) = Plastic displacement ------------ (Model =  , ,3).
!!   SVS(8*i  ) = Current yield stress ------------ (Model =  , ,3).
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
!! Temporary !!! until MATERIAL structure is set-up to carry material
!! properties
!!
      REAL(KIND(0D0))                                                          &
     &          EE(13,7)
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
!! Rotate unit vectors normal to joint planes.
!!
        NUMJT = 0  !  NINT(EE(1,3))
        IF (NUMJT .GT. 0) THEN
          dUxx = DTnext * Dxx
          dUyy = DTnext * Dyy
          dUzz = DTnext * Dzz
          dUxy = DTnext * (Dxy + Wxy)
          dUyx = DTnext * (Dxy - Wxy)
          dUxz = DTnext * (Dxz + Wxz)
          dUzx = DTnext * (Dxz - Wxz)
          dUyz = DTnext * (Dyz + Wyz)
          dUzy = DTnext * (Dyz - Wyz)
          DO i = 1,NUMJT
            Ax = STATE_VARIABLES(((8*i)-7))
            Ay = STATE_VARIABLES(((8*i)-6))
            Az = STATE_VARIABLES(((8*i)-5))
            dUax = Ax*dUxx + Ay*dUyx + Az*dUzx
            dUay = Ax*dUxy + Ay*dUyy + Az*dUzy
            dUaz = Ax*dUxz + Ay*dUyz + Az*dUzz
            dUaa = Ax*dUax + Ay*dUay + Az*dUaz
            Ax = Ax*(1.0 + dUaa) - dUax
            Ay = Ay*(1.0 + dUaa) - dUay
            Az = Az*(1.0 + dUaa) - dUaz
            Amiv = 1.0 / SQRT (Ax*Ax + Ay*Ay + Az*Az)
            STATE_VARIABLES(((8*i)-7)) = Ax * Amiv
            STATE_VARIABLES(((8*i)-6)) = Ay * Amiv
            STATE_VARIABLES(((8*i)-5)) = Az * Amiv
          ENDDO
        ENDIF
      ENDIF
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
      CALL MATERIAL_34_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
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
!! Rotate unit vectors normal to joint planes.
!!
        NUMJT = 0  !  NINT(EE(1,3))
        IF (NUMJT .GT. 0) THEN
          dUxx = DTnext * Dxx
          dUyy = DTnext * Dyy
          dUzz = DTnext * Dzz
          dUxy = DTnext * (Dxy + Wxy)
          dUyx = DTnext * (Dxy - Wxy)
          dUxz = DTnext * (Dxz + Wxz)
          dUzx = DTnext * (Dxz - Wxz)
          dUyz = DTnext * (Dyz + Wyz)
          dUzy = DTnext * (Dyz - Wyz)
          DO i = 1,NUMJT
            Ax = STATE_VARIABLES(((8*i)-7))
            Ay = STATE_VARIABLES(((8*i)-6))
            Az = STATE_VARIABLES(((8*i)-5))
            dUax = Ax*dUxx + Ay*dUyx + Az*dUzx
            dUay = Ax*dUxy + Ay*dUyy + Az*dUzy
            dUaz = Ax*dUxz + Ay*dUyz + Az*dUzz
            dUaa = Ax*dUax + Ay*dUay + Az*dUaz
            Ax = Ax*(1.0 + dUaa) - dUax
            Ay = Ay*(1.0 + dUaa) - dUay
            Az = Az*(1.0 + dUaa) - dUaz
            Amiv = 1.0 / SQRT (Ax*Ax + Ay*Ay + Az*Az)
            STATE_VARIABLES(((8*i)-7)) = Ax * Amiv
            STATE_VARIABLES(((8*i)-6)) = Ay * Amiv
            STATE_VARIABLES(((8*i)-5)) = Az * Amiv
          ENDDO
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_34_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Copyright (c) by KEY Associates; 23-JAN-1994 12:44:29.79
!!
!! Purpose: Interface to the three-dimensional material model JNTMAT.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          RL(4),                                                         &
     &          DTnext,                                                        &
     &          C(6,6),                                                        &
     &          D3D(6),                                                        &
     &          Einc(6),                                                       &
     &          Sinc(6),                                                       &
     &          STRESS(6),                                                     &
     &          STATE_VARIABLES(*)

!!
!! Temporary !!! until MATERIAL structure is set-up to carry material
!! properties
!!
      REAL(KIND(0D0))                                                          &
     &          EE(13,7)

!!
!! Initialize constants.
!!      YM = Young's modulus E
!!      PR = Poisson's ratio nu
!!      P1 = Elastic Lame parameter lambda
!!      P2 = Elastic Lame parameter mu, shear modulus
!!
      YM = MATERIAL(MatID)%PVAL(6)
      PR = MATERIAL(MatID)%PVAL(7)
      P1 = MATERIAL(MatID)%PVAL(8)
      P2 = MATERIAL(MatID)%PVAL(9)
!!
!! Define three-dimensional stretching components.
!!
      D3D(1) = Dxx
      D3D(2) = Dyy
      D3D(3) = Dzz
      D3D(4) = Dxy
      D3D(5) = Dxz
      D3D(6) = Dyz
!!
!! Save strain increment.
!!
      DO i = 1,6
        Einc(i) = DTnext * D3D(i)
      ENDDO
!!
!! Save old stress state.
!!
      DO i = 1,6
        Sinc(i) = STRESS(i)
      ENDDO
!!
!! Extract the number of joints and their characteristic length.
!!
      NUMJT = 0  !  NINT(EE(1,3))
      DO I = 1,4
        RL(I) = 0.0  !  EE(I+1,3)
      ENDDO
!!
!! Suppress tangent modulus calculation.
!!
      ICFLG = 0
!!
!! Convert I*2 to I*4
!!
      LELO = IO_UNIT%LELO
!!
!!    CALL JNTMAT
!!   &          (
!!   &          STRESS,D3D,DTnext,STATE_VARIABLES,YM,
!!   &          PR,NUMJT,RL,EE(1,4),C,ICFLG,LELO
!!   &          )
!!
!! Sound speeds squared * RHO(T)
!!
      DO i = 1,6
        Sinc(i) = STRESS(i) - Sinc(i)
      ENDDO
      CALL LOADING_MODULI (CL2,CS2,P1,P2,Sinc,Einc)
!!
      RETURN
      END
