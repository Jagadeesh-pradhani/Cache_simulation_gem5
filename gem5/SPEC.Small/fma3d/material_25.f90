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
      SUBROUTINE MATERIAL_25_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:59
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
!! Use triaxial orthotropic elastic material 35 material property
!! initialization
!!
      CALL MATERIAL_35_INIT (MatID)
!!
      RETURN
!!
      ENTRY MATERIAL_25_INI2 (STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! Initialize stress point w/ initial fiber orientation; global coordinates.
!!
      Ioff = Isv - 1
      STATE_VARIABLES(Ioff+1) = SECTION_2D(SecID)%Ax
      STATE_VARIABLES(Ioff+2) = SECTION_2D(SecID)%Ay
      STATE_VARIABLES(Ioff+3) = SECTION_2D(SecID)%Az
      STATE_VARIABLES(Ioff+4) = SECTION_2D(SecID)%Bx
      STATE_VARIABLES(Ioff+5) = SECTION_2D(SecID)%By
      STATE_VARIABLES(Ioff+6) = SECTION_2D(SecID)%Bz
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_35_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:59
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
!! Retrieve input Young's moduli and Poisson's ratios
!!
      Eaa = MATERIAL(MatID)%PVAL(6)
      Ebb = MATERIAL(MatID)%PVAL(7)
      Ecc = MATERIAL(MatID)%PVAL(8)
      Uba = MATERIAL(MatID)%PVAL(9)
      Uca = MATERIAL(MatID)%PVAl(10)
      Ucb = MATERIAL(MatID)%PVAL(11)
      Uab = Uba*(Eaa/Ebb)
      Uac = Uca*(Eaa/Ecc)
      Ubc = Ucb*(Ebb/Ecc)
!!
!! Check Poisson's ratios for consistency.
!!
      ERROR%COUNT = 0
      IF ((1.0-Uab-Uac) .LE. 0.0) THEN
        WRITE (MSG1,'(I8)') MATERIAL(MatID)%MatID
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'WARN'//                                                       &
     &    MSGL//'MATERIAL_35_INIT.001.01'//                                    &
     &    MSGL//'Material (MatID) Input Record ID:'//MSG1//                    &
     &    MSGL//'Orthotropic Poisson''s Ratios Fail Consistency:'//            &
     &    MSGL//'(1.0-Uab-Uac) Must Be Greater Than Zero.'                     &
     &    )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
      IF ((1.0-Uba-Ubc) .LE. 0.0) THEN
        WRITE (MSG1,'(I8)') MATERIAL(MatID)%MatID
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'WARN'//                                                       &
     &    MSGL//'MATERIAL_35_INIT.001.02'//                                    &
     &    MSGL//'Material (MatID) Input Record ID:'//MSG1//                    &
     &    MSGL//'Orthotropic Poisson''s Ratios Fail Consistency:'//            &
     &    MSGL//'(1.0-Uba-Ubc) Must Be Greater Than Zero.'                     &
     &    )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
      IF ((1.0-Ucb-Uca) .LE. 0.0) THEN
        WRITE (MSG1,'(I8)') MATERIAL(MatID)%MatID
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'WARN'//                                                       &
     &    MSGL//'MATERIAL_35_INIT.001.03'//                                    &
     &    MSGL//'Material (MatID) Input Record ID:'//MSG1//                    &
     &    MSGL//'Orthotropic Poisson''s Ratios Fail Consistency:'//            &
     &    MSGL//'(1.0-Ucb-Uca) Must Be Greater Than Zero.'                     &
     &    )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
      IF (ERROR%COUNT .GT. 0) THEN
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'MATERIAL_35_INIT.001.04'//                                    &
     &    MSGL//'Material (MatID) Input Record ID:'//MSG1//                    &
     &    MSGL//'One Or More Orthotropic Poisson''s Ratio '//                  &
     &    MSGL//'Pairs Have Failed The Consistency Check.'                     &
     &    )
      ENDIF
!!
!! Construct moduli matrix from input property data.
!!
      Det = 1. - Uba*Ucb*Uac - Uca*Uab*Ubc - Uac*Uca - Ubc*Ucb - Uab*Uba
      MATERIAL(MatID)%PVAL( 6) = Eaa * (1.0 + Ubc*Ucb) * (1.0 / Det)
      MATERIAL(MatID)%PVAL( 7) = Ebb * (1.0 + Uac*Uca) * (1.0 / Det)
      MATERIAL(MatID)%PVAL( 8) = Ecc * (1.0 + Uba*Uab) * (1.0 / Det)
      MATERIAL(MatID)%PVAL( 9) = Eaa * (Uba + Ubc*Uca) * (1.0 / Det)
      MATERIAL(MatID)%PVAL(10) = Eaa * (Uca + Uba*Ucb) * (1.0 / Det)
      MATERIAL(MatID)%PVAL(11) = Ebb * (Ucb + Uab*Uca) * (1.0 / Det)
!!
!! Estimate maximum uniaxial strain and shear stiffnesses.
!!
      RCL2 = MAX (MATERIAL(MatID)%PVAL( 6),                                    &
     &            MATERIAL(MatID)%PVAL( 7),                                    &
     &            MATERIAL(MatID)%PVAL( 8) )
      RCS2 = MAX (MATERIAL(MatID)%PVAL(12),                                    &
     &            MATERIAL(MatID)%PVAL(13),                                    &
     &            MATERIAL(MatID)%PVAL(14) )
      MATERIAL(MatID)%PVAL(15) = RCL2
      MATERIAL(MatID)%PVAL(16) = RCS2
!!
      RETURN
!!
      ENTRY MATERIAL_35_INI2 (STRESS,MatID,LupID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!!
      Ioff = Isv - 1
!!
!! Transform material axis vectors to global coordinates.
!!
      IF (LAYERING(LupID)%Isys .EQ. 0) THEN
        Ax = LAYERING(LupID)%Ax(1)
        Ay = LAYERING(LupID)%Ay(1)
        Az = LAYERING(LupID)%Az(1)
        Bx = LAYERING(LupID)%Bx(1)
        By = LAYERING(LupID)%By(1)
        Bz = LAYERING(LupID)%Bz(1)
        Cx = Ay*Bz - By*Az
        Cy = Az*Bx - Bz*Ax
        Cz = Ax*By - Bx*Ay
      ELSE
      ENDIF
!!
!!  Place material axis vectors in element state variable storage.
!!
      Qmag = 1.0 / SQRT (Ax*Ax + Ay*Ay + Az*Az)
      STATE_VARIABLES(Ioff+1) = Qmag * Ax
      STATE_VARIABLES(Ioff+2) = Qmag * Ay
      STATE_VARIABLES(Ioff+3) = Qmag * Az
      Qmag = 1.0 / SQRT (Bx*Bx + By*By + Bz*Bz)
      STATE_VARIABLES(Ioff+4) = Qmag * Bx
      STATE_VARIABLES(Ioff+5) = Qmag * By
      STATE_VARIABLES(Ioff+6) = Qmag * Bz
      Qmag = 1.0 / SQRT (Cx*Cx + Cy*Cy + Cz*Cz)
      STATE_VARIABLES(Ioff+7) = Qmag * Cx
      STATE_VARIABLES(Ioff+8) = Qmag * Cy
      STATE_VARIABLES(Ioff+9) = Qmag * Cz
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_45_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 21-FEB-1991 21:04:59
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
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! Use triaxial orthotropic elastic material 35 material property
!! initialization
!!
      CALL MATERIAL_35_INIT (MatID)
!!
      RETURN
!!
      ENTRY MATERIAL_45_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! Initialize stress point w/ initial fiber orientation; global coordinates.
!!
      Ioff = Isv - 1
      DO i = 0, (Ipts-1)*MATERIAL(MatID)%Nsv, MATERIAL(MatID)%Nsv
        STATE_VARIABLES(Ioff+i+1) = SECTION_2D(SecID)%Ax
        STATE_VARIABLES(Ioff+i+2) = SECTION_2D(SecID)%Ay
        STATE_VARIABLES(Ioff+i+3) = SECTION_2D(SecID)%Az
        STATE_VARIABLES(Ioff+i+4) = SECTION_2D(SecID)%Bx
        STATE_VARIABLES(Ioff+i+5) = SECTION_2D(SecID)%By
        STATE_VARIABLES(Ioff+i+6) = SECTION_2D(SecID)%Bz
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_55_INIT (MatID)
!!
!! Copyright (c) by KEY Associates; 3-AUG-1993 19:22:45.61
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
     &          STRESS(6,Ipts)
!!
!! INITIALIZE MATERIAL CONSTANTS.
!! (do nothing)
!!
      X = 1.0
      RETURN
!!
      ENTRY MATERIAL_55_INI2 (Ipts,STRESS,MatID,SecID,Isv,Nsv)
!!
!! INITIALIZE STATE VARIABLES.
!! (do nothing)
!!
      X = 1.0
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_25                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,INTERNAL_ENERGY,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:37:00
!!
!! FINITE STRAIN PLANE-STRESS ORTHOTROPIC MATERIAL MODEL.
!!
!!      STATE_VARIABLES(1) = Ar, r-component of material axis A.
!!      STATE_VARIABLES(2) = As, s-component of material axis A.
!!      STATE_VARIABLES(3) = Br, r-component of material axis B.
!!      STATE_VARIABLES(4) = Bs, s-component of material axis B.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching, spin
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Gr,Gs,                  & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(3),                                                     &
     &          VECTORS(9),                                                    &
     &          LOCAL_STRESS(6),                                               &
     &          STATE_VARIABLES(*)
      DATA                                                                     &
     &          VECTORS /8*0.0,1.0/
!!
!! Initialize material constants.
!!
      Cac = MATERIAL(MatID)%PVAL(10)
      Cbc = MATERIAL(MatID)%PVAL(11)
      Ccc = MATERIAL(MatID)%PVAL(8)
!!
!! Enforce plane stress condition by prescribing stretching in the normal
!! (or t) direction. This development assumes that the initial value of
!! Dzz is zero.
!!
      Dtt = -(Cac*Drr +Cbc*Dss)/Ccc
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
!! Put rotated stresses in temporary storage for use in MATERIAL_35_INTEGRATION.
!!
      LOCAL_STRESS(1) = Srr
      LOCAL_STRESS(2) = Sss
      LOCAL_STRESS(3) = 0.0
      LOCAL_STRESS(4) = Srs
      LOCAL_STRESS(5) = 0.0
      LOCAL_STRESS(6) = 0.0
!!
!! Rotate material directions A and B.
!!
      R1 = DTnext * Wrs
      Q1 = STATE_VARIABLES(1) + R1*STATE_VARIABLES(2)
      Q2 = STATE_VARIABLES(2) - R1*STATE_VARIABLES(1)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(1) = Q1 * Qim
      STATE_VARIABLES(2) = Q2 * Qim
      Q1 = STATE_VARIABLES(3) + R1*STATE_VARIABLES(4)
      Q2 = STATE_VARIABLES(4) - R1*STATE_VARIABLES(3)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(3) = Q1 * Qim
      STATE_VARIABLES(4) = Q2 * Qim
!!
!! Load material vectors into an array for use in the three-dimensional
!! module MATERIAL_35_INTEGRATION.
!!
      VECTORS(1) = STATE_VARIABLES(1)
      VECTORS(2) = STATE_VARIABLES(2)
      VECTORS(4) = STATE_VARIABLES(3)
      VECTORS(5) = STATE_VARIABLES(4)
!!
!! Internal energy.
!!
      INTERNAL_ENERGY = INTERNAL_ENERGY +                                      &
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Update stress. Note that this membrane stress calculation uses the
!! triaxial solid module MATERIAL_35_INTEGRATION to insure identical
!! orthotropic behavior and to conserve on material model coding.
!!
      CALL MATERIAL_35_INTEGRATION                                             &
     &          (                                                              &
     &          LOCAL_STRESS,VECTORS,                                          &
     &          DTnext,Drr,Dss,Dtt,Drs,Drt,Dst,Wrs,Wrt,Wst,                    &
     &          MatID                                                          &
     &          )
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
     &          (0.5 * DTnext) * (Drr*Srr + Dss*Sss + (Drs+Drs)*Srs)
!!
!! Return stresses to global storage.
!!
      STRESS(1) = Srr
      STRESS(2) = Sss
      STRESS(3) = Srs
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_35                                                   &
     &          (                                                              &
     &          STRESS,INTERNAL_ENERGY,STATE_VARIABLES,DTnext,MatID            &
     &          )
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:37:00
!!
!! FINITE STRAIN HYPO-ELASTIC ORTHOTROPIC MATERIAL MODEL.
!!
!! Module computes the current value of the stress STRESS(1:6), the
!! density times the longitudinal sound speed CL squared, and the
!! density times the shear sound speed CS squared in the current state,
!! as a function of the stretching Dxx,...,Dyz, the spin Wxy,Wxz,Wyz,
!! (the skew-symmetric part of the velocity gradient) and the old stress
!! state.
!!
!!      STATE_VARIABLES(1) = Ax, x-component of material axis A.
!!      STATE_VARIABLES(2) = Ay, y-component of material axis A.
!!      STATE_VARIABLES(3) = Az, z-component of material axis A.
!!      STATE_VARIABLES(4) = Bx, x-component of material axis B.
!!      STATE_VARIABLES(5) = By, y-component of material axis B.
!!      STATE_VARIABLES(6) = Bz, z-component of material axis B.
!!      STATE_VARIABLES(7) = Cx, x-component of material axis C.
!!      STATE_VARIABLES(8) = Cy, y-component of material axis C.
!!      STATE_VARIABLES(9) = Cz, z-component of material axis C.
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
      DTinc = DTnext
      dWxy = DTinc * Wxy
      dWxz = DTinc * Wxz
      dWyz = DTinc * Wyz
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
!! Rotate material-axis vectors from time n to time n+1. After the
!! rotation is applied, the material-axis vectors are re-normalized.
!!
        m = 0
        DO i = 1,3
          Qx =  dWxy*STATE_VARIABLES(m+2) + dWxz*STATE_VARIABLES(m+3)
          Qy = -dWxy*STATE_VARIABLES(m+1) + dWyz*STATE_VARIABLES(m+3)
          Qz = -dWxz*STATE_VARIABLES(m+1) - dWyz*STATE_VARIABLES(m+2)
          Qx = STATE_VARIABLES(m+1) - Qx
          Qy = STATE_VARIABLES(m+2) - Qy
          Qz = STATE_VARIABLES(m+3) - Qz
          Qmag = 1.0 / SQRT (Qx*Qx + Qy*Qy + Qz*Qz)
          STATE_VARIABLES(m+1) = Qx * Qmag
          STATE_VARIABLES(m+2) = Qy * Qmag
          STATE_VARIABLES(m+3) = Qz * Qmag
          m = m + 3
        ENDDO
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
      CALL MATERIAL_35_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,                                                        &
     &          STATE_VARIABLES,                                               &
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
!! Rotate material-axis vectors from time n to time n+1. After the
!! rotation is applied, the material-axis vectors are re-normalized.
!!
        m = 0
        DO i = 1,3
          Qx =  dWxy*STATE_VARIABLES(m+2) + dWxz*STATE_VARIABLES(m+3)
          Qy = -dWxy*STATE_VARIABLES(m+1) + dWyz*STATE_VARIABLES(m+3)
          Qz = -dWxz*STATE_VARIABLES(m+1) - dWyz*STATE_VARIABLES(m+2)
          Qx = STATE_VARIABLES(m+1) - Qx
          Qy = STATE_VARIABLES(m+2) - Qy
          Qz = STATE_VARIABLES(m+3) - Qz
          Qmag = 1.0 / SQRT (Qx*Qx + Qy*Qy + Qz*Qz)
          STATE_VARIABLES(m+1) = Qx * Qmag
          STATE_VARIABLES(m+2) = Qy * Qmag
          STATE_VARIABLES(m+3) = Qz * Qmag
          m = m + 3
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_35_INTEGRATION                                       &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:37:01
!!
!! ORTHOTROPIC ELASTIC MATERIAL MODEL.
!!
!!      STATE_VARIABLES(1) = Ax, x-component of material axis A.
!!      STATE_VARIABLES(2) = Ay, y-component of material axis A.
!!      STATE_VARIABLES(3) = Az, z-component of material axis A.
!!      STATE_VARIABLES(4) = Bx, x-component of material axis B.
!!      STATE_VARIABLES(5) = By, y-component of material axis B.
!!      STATE_VARIABLES(6) = Bz, z-component of material axis B.
!!      STATE_VARIABLES(7) = Cx, x-component of material axis C.
!!      STATE_VARIABLES(8) = Cy, y-component of material axis C.
!!      STATE_VARIABLES(9) = Cz, z-component of material axis C.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          STATE_VARIABLES(*),                                            &
     &          S(6),D(6)
!!
!! Retrieve stress from global storage. The vector storage sequence is
!! STRESS(1:6) = (Sxx,Syy,Szz,Sxy,Sxz,Syz)
!!
      Sxx = STRESS(1)
      Syy = STRESS(2)
      Szz = STRESS(3)
      Sxy = STRESS(4)
      Sxz = STRESS(5)
      Syz = STRESS(6)
!!
      Ax = STATE_VARIABLES(1)
      Ay = STATE_VARIABLES(2)
      Az = STATE_VARIABLES(3)
      Bx = STATE_VARIABLES(4)
      By = STATE_VARIABLES(5)
      Bz = STATE_VARIABLES(6)
      Cx = STATE_VARIABLES(7)
      Cy = STATE_VARIABLES(8)
      Cz = STATE_VARIABLES(9)
!!
!! Rotate stress from global axes to local material axes.
!!
      Sxa = Sxx*Ax + Sxy*Ay + Sxz*Az
      Sxb = Sxx*Bx + Sxy*By + Sxz*Bz
      Sxc = Sxx*Cx + Sxy*Cy + Sxz*Cz
      Sya = Sxy*Ax + Syy*Ay + Syz*Az
      Syb = Sxy*Bx + Syy*By + Syz*Bz
      Syc = Sxy*Cx + Syy*Cy + Syz*Cz
      Sza = Sxz*Ax + Syz*Ay + Szz*Az
      Szb = Sxz*Bx + Syz*By + Szz*Bz
      Szc = Sxz*Cx + Syz*Cy + Szz*Cz
!!
!! The local vector storage sequence is S(1:6) = (Saa,Sbb,Scc,Sab,Sac,Sbc)
!!
      S(1) = Ax*Sxa + Ay*Sya + Az*Sza
      S(4) = Ax*Sxb + Ay*Syb + Az*Szb
      S(5) = Ax*Sxc + Ay*Syc + Az*Szc
      S(2) = Bx*Sxb + By*Syb + Bz*Szb
      S(6) = Bx*Sxc + By*Syc + Bz*Szc
      S(3) = Cx*Sxc + Cy*Syc + Cz*Szc
!!
!! Rotate stretching from global axes to local material axes.
!!
      Dxa = Dxx*Ax + Dxy*Ay + Dxz*Az
      Dxb = Dxx*Bx + Dxy*By + Dxz*Bz
      Dxc = Dxx*Cx + Dxy*Cy + Dxz*Cz
      Dya = Dxy*Ax + Dyy*Ay + Dyz*Az
      Dyb = Dxy*Bx + Dyy*By + Dyz*Bz
      Dyc = Dxy*Cx + Dyy*Cy + Dyz*Cz
      Dza = Dxz*Ax + Dyz*Ay + Dzz*Az
      Dzb = Dxz*Bx + Dyz*By + Dzz*Bz
      Dzc = Dxz*Cx + Dyz*Cy + Dzz*Cz
!!
!! The local vector storage sequence is D(1:6) = (Daa,Dbb,Dcc,Dab,Dac,Dbc)
!!
      D(1) = DTnext * (Ax*Dxa + Ay*Dya + Az*Dza)
      D(4) = DTnext * (Ax*Dxb + Ay*Dyb + Az*Dzb)
      D(5) = DTnext * (Ax*Dxc + Ay*Dyc + Az*Dzc)
      D(2) = DTnext * (Bx*Dxb + By*Dyb + Bz*Dzb)
      D(6) = DTnext * (Bx*Dxc + By*Dyc + Bz*Dzc)
      D(3) = DTnext * (Cx*Dxc + Cy*Dyc + Cz*Dzc)
!!
!! Compute elastic stress increments. (The moduli matrix Cij was constructed
!! from the input quatities Eaa, Ebb, Ecc, nu-ab, nu-ac and nu-bc.)
!!
!!      Saa = Caa*Eaa + Cab*Ebb + Cac*Ecc
!!      Sbb = Cab*Eaa + Cbb*Ebb + Cbc*Ecc
!!      Scc = Cac*Eaa + Cbc*Ebb + Ccc*Ecc
!!      Sab = Gab*(2.0*Eab)
!!      Sac = Gac*(2.0*Eac)
!!      Sbc = Gbc*(2.0*Ebc)
!!
      S(1) = S(1) + MATERIAL(MatID)%PVAL( 6)*D(1)                              &
     &            + MATERIAL(MatID)%PVAL( 9)*D(2)                              &
     &            + MATERIAL(MatID)%PVAL(10)*D(3)
      S(2) = S(2) + MATERIAL(MatID)%PVAL( 9)*D(1)                              &
     &            + MATERIAL(MatID)%PVAL( 7)*D(2)                              &
     &            + MATERIAL(MatID)%PVAL(11)*D(3)
      S(3) = S(3) + MATERIAL(MatID)%PVAL(10)*D(1)                              &
     &            + MATERIAL(MatID)%PVAL(11)*D(2)                              &
     &            + MATERIAL(MatID)%PVAL( 8)*D(3)
      S(4) = S(4) + MATERIAL(MatID)%PVAL(12) * (D(4) + D(4))
      S(5) = S(5) + MATERIAL(MatID)%PVAL(13) * (D(5) + D(5))
      S(6) = S(6) + MATERIAL(MatID)%PVAL(14) * (D(6) + D(6))
!!
!! Evaluate stress for possible composite failure mechanisms.
!!
!!      (put failure criteria for composite material here)
!!
!! Return stress state to global coordinates.
!!
      Sax = S(1)*Ax + S(4)*Bx + S(5)*Cx
      Say = S(1)*Ay + S(4)*By + S(5)*Cy
      Saz = S(1)*Az + S(4)*Bz + S(5)*Cz
      Sbx = S(4)*Ax + S(2)*Bx + S(6)*Cx
      Sby = S(4)*Ay + S(2)*By + S(6)*Cy
      Sbz = S(4)*Az + S(2)*Bz + S(6)*Cz
      Scx = S(5)*Ax + S(6)*Bx + S(3)*Cx
      Scy = S(5)*Ay + S(6)*By + S(3)*Cy
      Scz = S(5)*Az + S(6)*Bz + S(3)*Cz
!!
!! The global storage sequence is STRESS(1:6) = (Sxx,Syy,Szz,Sxy,Sxz,Syz).
!!
      STRESS(1) = Ax*Sax + Bx*Sbx + Cx*Scx
      STRESS(4) = Ax*Say + Bx*Sby + Cx*Scy
      STRESS(5) = Ax*Saz + Bx*Sbz + Cx*Scz
      STRESS(2) = Ay*Say + By*Sby + Cy*Scy
      STRESS(6) = Ay*Saz + By*Sbz + Cy*Scz
      STRESS(3) = Az*Saz + Bz*Sbz + Cz*Scz
!!
!! Retrieve uniaxial strain and shear stiffnesses.
!!
      SOUND_SPEED%RCL2 = MATERIAL(MatID)%PVAL(15)
      SOUND_SPEED%RCS2 = MATERIAL(MatID)%PVAL(16)
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_45                                                   &
     &          (                                                              &
     &          STRESS,STATE_VARIABLES,                                        &
     &          DTnext,Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                                &
     &          MatID                                                          &
     &          )
!!
!!
!! Copyright (c) by KEY Associates; 28-NOV-1991 11:37:00
!!
!! FINITE STRAIN PLANE-STRESS ORTHOTROPIC MATERIAL MODEL.
!!
!!      STATE_VARIABLES(1) = Ax, x-component of material axis A.
!!      STATE_VARIABLES(2) = Ay, y-component of material axis A.
!!      STATE_VARIABLES(3) = Bx, x-component of material axis B.
!!      STATE_VARIABLES(4) = By, y-component of material axis B.
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          DTnext,                                                        &
     &          STRESS(6),                                                     &
     &          VECTORS(9),                                                    &
     &          STATE_VARIABLES(*)
      DATA                                                                     &
     &          VECTORS /8*0.0,1.0/
!!
!! Initialize material constants.
!!
      Cac = MATERIAL(MatID)%PVAL(10)
      Cbc = MATERIAL(MatID)%PVAL(11)
      Ccc = MATERIAL(MatID)%PVAL( 8)
!!
!! Enforce plane stress condition by prescribing stretching in the normal
!! (or z) direction.  This development assumes that the initial value of
!! Dzz is zero.
!!
      Dzz = -(Cac*Dxx +Cbc*Dyy)/Ccc
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
!! Rotate material directions A and B.
!!
      Q1 = STATE_VARIABLES(1) + R1*STATE_VARIABLES(2)
      Q2 = STATE_VARIABLES(2) - R1*STATE_VARIABLES(1)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(1) = Q1 * Qim
      STATE_VARIABLES(2) = Q2 * Qim
      Q1 = STATE_VARIABLES(3) + R1*STATE_VARIABLES(4)
      Q2 = STATE_VARIABLES(4) - R1*STATE_VARIABLES(3)
      Qim = 1.0 / SQRT (Q1*Q1 + Q2*Q2)
      STATE_VARIABLES(3) = Q1 * Qim
      STATE_VARIABLES(4) = Q2 * Qim
!!
!! Load material vectors into an array for use in the three-dimensional
!! module MATERIAL_35_INTEGRATION.
!!
      VECTORS(1) = STATE_VARIABLES(1)
      VECTORS(2) = STATE_VARIABLES(2)
      VECTORS(4) = STATE_VARIABLES(3)
      VECTORS(5) = STATE_VARIABLES(4)
!!
!! Update stress. Note that this shell stress calculation uses the
!! triaxial solid module MATERIAL_35_INTEGRATION to insure identical
!! orthotropic behavior and to conserve on material model coding.
!!
      CALL MATERIAL_35_INTEGRATION                                             &
     &          (                                                              &
     &          STRESS,VECTORS,                                                &
     &          DTnext,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,                    &
     &          MatID                                                          &
     &          )
!!
      RETURN
      END
