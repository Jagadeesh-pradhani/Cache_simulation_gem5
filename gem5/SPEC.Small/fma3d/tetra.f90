      SUBROUTINE TETRA_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:11:43
!!
      USE shared_common_data
      USE tetra_
      USE material_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: Isv,MatID
      LOGICAL :: FOUND
!!
!! Loop over all tetrahedral elements.
!!
      DO N = 1,NUMTX
!!
!! Gather element motion (position and velocity at time equal to zero).
!!
        DO i = 1,4
          EMOTION(i) = MOTION(TETRA(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        TETRA(N)%RES%Time   = 0.0
        TETRA(N)%RES%DTnext = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
        Isv = TETRA(N)%PAR%Isv
        MatID = TETRA(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies, only the
!! element volume calculation is required in order to construct the mass
!! matrix.)
!!
        CALL TETRA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Save initial element volume for later volume strain calculations.
!!
        TETRA(N)%PAR%Volume = TETRA(N)%RES%Volume
        IF (TETRA(N)%RES%Volume .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') TETRA(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'TETRA_INITIALIZATION.001.00'//                          &
     &          MSGL//'TXEL (4-Node Tetrahedron) Element ID:'//MSG1//          &
     &          MSGL//'Has A Zero Or Negative Volume.'                         &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL TETRA_MASS ( NEL,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (TETRA(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE( ELsize)
          DO i = 2,4
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &              + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                      &
     &              + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          TETRA(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Find initial sound speeds and stress.
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (30)
            CALL MATERIAL_30                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (31)
            CALL MATERIAL_31                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (32)
            CALL MATERIAL_32                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (33)
            CALL MATERIAL_33                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (34)
            CALL MATERIAL_34                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (35)
            CALL MATERIAL_35                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (36)
            CALL MATERIAL_36                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (37)
            CALL MATERIAL_37                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (38)
            CALL MATERIAL_38                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (39)
            CALL MATERIAL_39                                                   &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') TETRA(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'TETRA_INITIALIZATION.002.00'//                          &
     &          MSGL//'TXEL (4-Node tetrahedron) Element ID:'//MSG1//          &
     &          MSGL//'References An Unknown Material Model:'//MSG2            &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL TETRA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
          IF (NUMNR .NE. 0) THEN
            Itype = 2  !  Tetraahedron
            CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
          ENDIF
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTTtx = MAX (TIMSIM%DTTtx,TETRA(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = TETRA(N)%RES%DTelt .LT. TIMSIM%DTTet(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTTet(j + 1) = TIMSIM%DTTet(j)
              TIMSIM%TETRA(j + 1) = TIMSIM%TETRA(j)
            ENDDO
          ENDIF
          TIMSIM%DTTet(i) = TETRA(N)%RES%DTelt
          TIMSIM%TETRA(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE TETRA_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:11:49
!!
      USE shared_common_data
      USE tetra_
      USE material_
      USE node_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Isv,                                                           &
     &          MatID
!!
!! Loop over all tetrahdral elements.
!!
      DO N = 1,NUMTX
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,TETRA(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (TETRA(N)%PAR%ParID .GE. 0) THEN
!!
!! Count element execution
!!
            COUNTER%TETRA = COUNTER%TETRA + 1
!!
!! Gather element coordinates and motion.
!!
            DO i = 1,4
              EMOTION(i) = MOTION(TETRA(N)%PAR%IX(i))
            ENDDO
!!
!! Increment element clock.
!!
            TETRA(N)%RES%Time = TETRA(N)%RES%Time + TETRA(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            DO i = 1,4
              QA = NODE(TETRA(N)%PAR%IX(i))%Time - TETRA(N)%RES%Time
              EMOTION(i)%Ux =                                                  &
     &          EMOTION(i)%Ux - QA * EMOTION(i)%Vx
              EMOTION(i)%Uy =                                                  &
     &          EMOTION(i)%Uy - QA * EMOTION(i)%Vy
              EMOTION(i)%Uz =                                                  &
     &          EMOTION(i)%Uz - QA * EMOTION(i)%Vz
            ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
            Isv = TETRA(N)%PAR%Isv
            MatID = TETRA(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation evaluated at time n-1/2.
!!
            CALL TETRA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Incremental constitutive model evaluation.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (30)
              CALL MATERIAL_30                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (31)
              CALL MATERIAL_31                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (32)
              CALL MATERIAL_32                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (33)
              CALL MATERIAL_33                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (34)
              CALL MATERIAL_34                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (35)
              CALL MATERIAL_35                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (36)
              CALL MATERIAL_36                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (37)
              CALL MATERIAL_37                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (38)
              CALL MATERIAL_38                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (39)
              CALL MATERIAL_39                                                 &
     &          (                                                              &
     &          TETRA(N)%RES%STRESS,                                           &
     &          TETRA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          TETRA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Divergence operator evaluated at time n+1.
!!
!!!           IF (CONTROL%MIDINT .NE. 0) THEN
!!!             CALL TETRA_DIVERGENCE_OPERATOR ( NEL,MatID )
!!!           ENDIF
!!
!! Update stress divergence, viscosity stress and time step.
!!
            CALL TETRA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
            IF (NUMNR .NE. 0) THEN
              Itype = 2  !  Tetrahedron
              CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
            ENDIF
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTTtx = MAX (TIMSIM%DTTtx,TETRA(N)%RES%DTelt)
          IF (TETRA(N)%RES%DTelt .LT. TIMSIM%DTTet(1)) THEN
            TIMSIM%DTTet(1) = TETRA(N)%RES%DTelt
            TIMSIM%TETRA(1) = N
          ENDIF
!!
!! End of time-to-subcycle if-test.
!!
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE TETRA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:11:53
!!
      USE shared_common_data
      USE tetra_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          X(4),Y(4),Z(4)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Compute position of nodal points at time n+1 (whole-interval evaluation).
!!
      DO i = 1,4
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Construct useful differences.
!!
      X12 = X(1) - X(2)
      X13 = X(1) - X(3)
      X14 = X(1) - X(4)
      X24 = X(2) - X(4)
      X34 = X(3) - X(4)
      Y12 = Y(1) - Y(2)
      Y13 = Y(1) - Y(3)
      Y14 = Y(1) - Y(4)
      Y24 = Y(2) - Y(4)
      Y34 = Y(3) - Y(4)
      Z12 = Z(1) - Z(2)
      Z13 = Z(1) - Z(3)
      Z14 = Z(1) - Z(4)
      Z24 = Z(2) - Z(4)
      Z34 = Z(3) - Z(4)
!!
!! Gradient operators.
!!
      Bx(1) = 0.16666666667 * (Y34*Z24 - Y24*Z34)
      Bx(2) = 0.16666666667 * (Y13*Z14 - Y14*Z13)
      Bx(3) = 0.16666666667 * (Y14*Z12 - Y12*Z14)
      Bx(4) = 0.16666666667 * (Y12*Z13 - Y13*Z12)
!!
      By(1) = 0.16666666667 * (Z34*X24 - Z24*X34)
      By(2) = 0.16666666667 * (Z13*X14 - Z14*X13)
      By(3) = 0.16666666667 * (Z14*X12 - Z12*X14)
      By(4) = 0.16666666667 * (Z12*X13 - Z13*X12)
!!
      Bz(1) = 0.16666666667 * (X34*Y24 - X24*Y34)
      Bz(2) = 0.16666666667 * (X13*Y14 - X14*Y13)
      Bz(3) = 0.16666666667 * (X14*Y12 - X12*Y14)
      Bz(4) = 0.16666666667 * (X12*Y13 - X13*Y12)
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      TETRA(NEL)%RES%Volume = 0.0
      DO i = 1,4
        TETRA(NEL)%RES%Volume = TETRA(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0D+0 / TETRA(NEL)%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
!!
!! Velocity gradients, stretching, and spin.
!!
      Vxx = 0.0
      Vyx = 0.0
      Vzx = 0.0
      Vxy = 0.0
      Vyy = 0.0
      Vzy = 0.0
      Vxz = 0.0
      Vyz = 0.0
      Vzz = 0.0
      DO i = 1,4
        Vxx = Vxx + EMOTION(i)%Vx * Bx(i)
        Vyx = Vyx + EMOTION(i)%Vy * Bx(i)
        Vzx = Vzx + EMOTION(i)%Vz * Bx(i)
        Vxy = Vxy + EMOTION(i)%Vx * By(i)
        Vyy = Vyy + EMOTION(i)%Vy * By(i)
        Vzy = Vzy + EMOTION(i)%Vz * By(i)
        Vxz = Vxz + EMOTION(i)%Vx * Bz(i)
        Vyz = Vyz + EMOTION(i)%Vy * Bz(i)
        Vzz = Vzz + EMOTION(i)%Vz * Bz(i)
      ENDDO
      Vxx = Vxx * Vin
      Vyx = Vyx * Vin
      Vzx = Vzx * Vin
      Vxy = Vxy * Vin
      Vyy = Vyy * Vin
      Vzy = Vzy * Vin
      Vxz = Vxz * Vin
      Vyz = Vyz * Vin
      Vzz = Vzz * Vin
!!
      Dxx = Vxx
      Dyy = Vyy
      Dzz = Vzz
      Dxy = 0.5 * (Vxy + Vyx)
      Dxz = 0.5 * (Vxz + Vzx)
      Dyz = 0.5 * (Vyz + Vzy)
      Wxy = Vxy - Dxy
      Wxz = Vxz - Dxz
      Wyz = Vyz - Dyz
!!
      RETURN
      END
!!_
      SUBROUTINE TETRA_DIVERGENCE_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:11:57
!!
      USE shared_common_data
      USE tetra_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          X(4),Y(4),Z(4)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Current position of nodal points, time n+1.
!!
      DO i = 1,4
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Construct useful differences.
!!
      X12 = X(1) - X(2)
      X13 = X(1) - X(3)
      X14 = X(1) - X(4)
      X24 = X(2) - X(4)
      X34 = X(3) - X(4)
      Y12 = Y(1) - Y(2)
      Y13 = Y(1) - Y(3)
      Y14 = Y(1) - Y(4)
      Y24 = Y(2) - Y(4)
      Y34 = Y(3) - Y(4)
      Z12 = Z(1) - Z(2)
      Z13 = Z(1) - Z(3)
      Z14 = Z(1) - Z(4)
      Z24 = Z(2) - Z(4)
      Z34 = Z(3) - Z(4)
!!
!! Gradient operators.
!!
      Bx(1) = 0.16666666667 * (Y34*Z24 - Y24*Z34)
      Bx(2) = 0.16666666667 * (Y13*Z14 - Y14*Z13)
      Bx(3) = 0.16666666667 * (Y14*Z12 - Y12*Z14)
      Bx(4) = 0.16666666667 * (Y12*Z13 - Y13*Z12)
!!
      By(1) = 0.16666666667 * (Z34*X24 - Z24*X34)
      By(2) = 0.16666666667 * (Z13*X14 - Z14*X13)
      By(3) = 0.16666666667 * (Z14*X12 - Z12*X14)
      By(4) = 0.16666666667 * (Z12*X13 - Z13*X12)
!!
      Bz(1) = 0.16666666667 * (X34*Y24 - X24*Y34)
      Bz(2) = 0.16666666667 * (X13*Y14 - X14*Y13)
      Bz(3) = 0.16666666667 * (X14*Y12 - X12*Y14)
      Bz(4) = 0.16666666667 * (X12*Y13 - X13*Y12)
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      TETRA(NEL)%RES%Volume = 0.0
      DO i = 1,4
        TETRA(NEL)%RES%Volume = TETRA(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0D+0 / TETRA(NEL)%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
!!
      RETURN
      END
!!_
      SUBROUTINE TETRA_MASS ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:12:02
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one fourth at each node.
!!
      USE shared_common_data
      USE tetra_
      USE material_
      USE node_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Compute one fourth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = 0.25 * Density * TETRA(NEL)%RES%Volume
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,4
        NODE(TETRA(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(TETRA(NEL)%PAR%IX(i))%Mass + QMass
        MATERIAL(MatID)%Mass = MATERIAL(MatID)%Mass + QMass
        Px = EMOTION(i)%Px
        Py = EMOTION(i)%Py
        Pz = EMOTION(i)%Pz
        MATERIAL(MatID)%Xcm = MATERIAL(MatID)%Xcm + QMass * Px
        MATERIAL(MatID)%Ycm = MATERIAL(MatID)%Ycm + QMass * Py
        MATERIAL(MatID)%Zcm = MATERIAL(MatID)%Zcm + QMass * Pz
!!
!! Compute inertia tensor B wrt the origin from nodal point masses.
!!
        MATERIAL(MatID)%Bxx = MATERIAL(MatID)%Bxx + (Py*Py+Pz*Pz)*QMass
        MATERIAL(MatID)%Byy = MATERIAL(MatID)%Byy + (Px*Px+Pz*Pz)*QMass
        MATERIAL(MatID)%Bzz = MATERIAL(MatID)%Bzz + (Px*Px+Py*Py)*QMass
        MATERIAL(MatID)%Bxy = MATERIAL(MatID)%Bxy - Px*Py*QMass
        MATERIAL(MatID)%Bxz = MATERIAL(MatID)%Bxz - Px*Pz*QMass
        MATERIAL(MatID)%Byz = MATERIAL(MatID)%Byz - Py*Pz*QMass
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE TETRA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 18-NOV-1991 12:11:11
!!
      USE shared_common_data
      USE tetra_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Retrieve material properties.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      Bulk_Ln = MATERIAL(MatID)%PVAL(2)
      Bulk_Qd = MATERIAL(MatID)%PVAL(3)
      HG_Visc = MATERIAL(MatID)%PVAL(4)
!!
!! Find sound speed.
!!
      SOUND_SPEED%Density =                                                    &
     &          Density * TETRA(NEL)%PAR%Volume / TETRA(NEL)%RES%Volume
!!
      CSQ = MAX (SOUND_SPEED%RCL2,SOUND_SPEED%RCS2)                            &
     &    / SOUND_SPEED%Density
      IF (CSQ .GT. 0.0) THEN
        CX = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') TETRA(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'TETRA_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'TXEL (4-Node Tetrahedron) Element ID:'//MSG1//          &
     &          MSGL//'Sound Speed Imaginary (C**2 Is Zero Or Negative)'       &
     &          )
        CX = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Calculate generalized element dimension.
!!
      Dx = 1.01D+0 / Delta
!!
!! Artificial bulk viscosity pressure.
!!
      Dkk = Dxx + Dyy + Dzz
      IF (Dkk .LT. 0.0) THEN
        QC = Bulk_Qd * Bulk_Qd * Dx * ABS(Dkk) + Bulk_Ln * CX
        QP = SOUND_SPEED%Density * Dx * Dkk * QC
      ELSE
        QC = 0.0
        QP = 0.0
      ENDIF
!!
!! Critical time step calculation.
!!
      TETRA(NEL)%RES%DTelt = Dx / (QC + SQRT (QC*QC + CX*CX))
!!
!! Divergence of the (mean) stress.
!!
      Qxx = TETRA(NEL)%RES%Stress(1) + QP
      Qxy = TETRA(NEL)%RES%Stress(4)
      Qxz = TETRA(NEL)%RES%Stress(5)
      Qyy = TETRA(NEL)%RES%Stress(2) + QP
      Qyz = TETRA(NEL)%RES%Stress(6)
      Qzz = TETRA(NEL)%RES%Stress(3) + QP
!!
!! Accumulate element divergence results.
!!
      DO i = 1,4
        TETRA(NEL)%RES%Xint(i) = Bx(i)*Qxx + By(i)*Qxy + Bz(i)*Qxz
        TETRA(NEL)%RES%Yint(i) = Bx(i)*Qxy + By(i)*Qyy + Bz(i)*Qyz
        TETRA(NEL)%RES%Zint(i) = Bx(i)*Qxz + By(i)*Qyz + Bz(i)*Qzz
      ENDDO
!!
      RETURN
      END
