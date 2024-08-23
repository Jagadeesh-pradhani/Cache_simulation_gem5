      SUBROUTINE HEXAH_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 12-FEB-1991 20:22:44
!!
      USE shared_common_data
      USE hexah_
      USE motion_
      USE material_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: Isv,MatID
      LOGICAL :: FOUND
!!
!! Loop over all hexahedral elements.
!!
      DO N = 1,NUMHX
!!
!! Gather element motion (position and velocity at time equal to zero).
!!
        DO i = 1,8
          EMOTION(i) = MOTION(HEXAH(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        HEXAH(N)%RES%Time   = 0.0
        HEXAH(N)%RES%DTnext = 0.0
!!
!! Initialize hourglass restoring forces.
!!
        DO i = 1,4
          HEXAH(N)%RES%Px(i) = 0
          HEXAH(N)%RES%Py(i) = 0
          HEXAH(N)%RES%Pz(i) = 0
        ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
        Isv = HEXAH(N)%PAR%Isv
        MatID = HEXAH(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies, only the
!! element volume calculation is required in order to construct the mass
!! matrix.)
!!
        CALL HEXAH_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Save initial element volume for later volume strain calculations.
!!
        HEXAH(N)%PAR%Volume = HEXAH(N)%RES%Volume

        IF (HEXAH(N)%RES%Volume .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') HEXAH(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'HEXAH_INITIALIZATION.001.00'//                          &
     &          MSGL//'HXEL (8-Node Hexahedron) Element ID:'//MSG1//           &
     &          MSGL//'Has A Zero Or Negative Volume.'                         &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL HEXAH_MASS ( NEL,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (HEXAH(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,8
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &            + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                        &
     &            + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          HEXAH(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
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
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (31)
            CALL MATERIAL_31                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (32)
            CALL MATERIAL_32                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (33)
            CALL MATERIAL_33                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (34)
            CALL MATERIAL_34                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (35)
            CALL MATERIAL_35                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (36)
            CALL MATERIAL_36                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (37)
            CALL MATERIAL_37                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (38)
            CALL MATERIAL_38                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (39)
            CALL MATERIAL_39                                                   &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') HEXAH(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'HEXAH_INITIALIZATION.002.00'//                          &
     &          MSGL//'HXEL (8-Node Hexahedron) Element ID:'//MSG1//           &
     &          MSGL//'References An Unknown Material Model:'//MSG2            &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL HEXAH_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
          IF (NUMNR .NE. 0) THEN
            Itype = 0  !  Hexahedron
            CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
          ENDIF
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTHxx = MAX (TIMSIM%DTHxx,HEXAH(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = HEXAH(N)%RES%DTelt .LT. TIMSIM%DTHex(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTHex(j + 1) = TIMSIM%DTHex(j)
              TIMSIM%Hexah(j + 1) = TIMSIM%Hexah(j)
            ENDDO
          ENDIF
          TIMSIM%DTHex(i) = HEXAH(N)%RES%DTelt
          TIMSIM%Hexah(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE HEXAH_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE hexah_
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
!! Loop over all hexahedral elements.
!!
      DO N = 1,NUMHX
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,HEXAH(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (HEXAH(N)%PAR%ParID .GE. 0) THEN
!!
!! Count element execution.
!!
            COUNTER%HEXAH = COUNTER%HEXAH + 1
!!
!! Gather element coordinates and motion.
!!
            DO i = 1,8
              EMOTION(i) = MOTION(HEXAH(N)%PAR%IX(i))
            ENDDO
!!
!! Increment element clock.
!!
            HEXAH(N)%RES%Time = HEXAH(N)%RES%Time + HEXAH(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            DO i = 1,8
              QA = NODE(HEXAH(N)%PAR%IX(i))%Time - HEXAH(N)%RES%Time
              EMOTION(i)%Ux =                                                  &
     &        EMOTION(i)%Ux - QA * EMOTION(i)%Vx
              EMOTION(i)%Uy =                                                  &
     &        EMOTION(i)%Uy - QA * EMOTION(i)%Vy
              EMOTION(i)%Uz =                                                  &
     &        EMOTION(i)%Uz - QA * EMOTION(i)%Vz
            ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
            Isv = HEXAH(N)%PAR%Isv
            MatID = HEXAH(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation.
!!
            CALL HEXAH_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Incremental constitutive model evaluation.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (30)
              CALL MATERIAL_30                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (31)
              CALL MATERIAL_31                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (32)
              CALL MATERIAL_32                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (33)
              CALL MATERIAL_33                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (34)
              CALL MATERIAL_34                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (35)
              CALL MATERIAL_35                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (36)
              CALL MATERIAL_36                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (37)
              CALL MATERIAL_37                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (38)
              CALL MATERIAL_38                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (39)
              CALL MATERIAL_39                                                 &
     &          (                                                              &
     &          HEXAH(N)%RES%STRESS,                                           &
     &          HEXAH(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          HEXAH(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
            IF (MATERIAL(MatID)%PVAL(5) .GT. 0.0) THEN
              CALL HEXAH_HOURGLASS_FORCES ( NEL,MatID )
            ENDIF
!!
!! Divergence operator evaluated at time n+1.
!!
!!!           IF (CONTROL%MIDINT .NE. 0) THEN
!!!             CALL HEXAH_DIVERGENCE_OPERATOR ( NEL,MatID )
!!!           ENDIF
!!
!! Update stress divergence, viscosity stress and time step.
!!
            CALL HEXAH_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
            IF (NUMNR .NE. 0) THEN
              Itype = 0  !  Hexahedron
              CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
            ENDIF
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTHxx = MAX (TIMSIM%DTHxx,HEXAH(N)%RES%DTelt)
          IF (HEXAH(N)%RES%DTelt .LT. TIMSIM%DTHex(1)) THEN
            TIMSIM%DTHex(1) = HEXAH(N)%RES%DTelt
            TIMSIM%Hexah(1) = N
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
      SUBROUTINE HEXAH_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE hexah_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          I2(8),I3(8),I4(8),I5(8),I6(8),I8(8)
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,4,1,8,5,6,7/, I3 /3,4,1,2,7,8,5,6/,                    &
     &          I4 /4,1,2,3,6,7,8,5/, I5 /5,6,7,8,1,2,3,4/,                    &
     &          I6 /6,7,8,5,4,1,2,3/, I8 /8,5,6,7,2,3,4,1/
!!
!! Compute position of nodal points at time n+1 (whole-interval evaluation).
!!
      DO i = 1,8
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,8
        Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))               &
     &             +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))             &
     &             +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))             &
     &             +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                               &
     &             +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                               &
     &             +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))               &
     &             +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))             &
     &             +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))             &
     &             +Z(I3(i))*(X(I2(i))-X(I4(i)))                               &
     &             +Z(I6(i))*(X(I5(i))-X(I2(i)))                               &
     &             +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))               &
     &             +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))             &
     &             +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))             &
     &             +X(I3(i))*(Y(I2(i))-Y(I4(i)))                               &
     &             +X(I6(i))*(Y(I5(i))-Y(I2(i)))                               &
     &             +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      HEXAH(NEL)%RES%Volume = 0.0
      DO i = 1,8
        HEXAH(NEL)%RES%Volume = HEXAH(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0 / HEXAH(NEL)%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,8
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
      DO i = 1,8
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
!! Anti-hourglass gradients.
!!
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(4)) -(X(5)+X(6)) +(X(7)+X(8))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(4)) -(X(5)-X(6)) +(X(7)-X(8))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(4)) -(Y(5)+Y(6)) +(Y(7)+Y(8))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(4)) -(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(4)) -(Z(5)+Z(6)) +(Z(7)+Z(8))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(4)) -(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
!!
      DO i = 1,4
        Gx(i) = -(Vxx*Hx(i) + Vxy*Hy(i) + Vxz*Hz(i)) * Delta
        Gy(i) = -(Vyx*Hx(i) + Vyy*Hy(i) + Vyz*Hz(i)) * Delta
        Gz(i) = -(Vzx*Hx(i) + Vzy*Hy(i) + Vzz*Hz(i)) * Delta
      ENDDO
      Gx(1) = Gx(1) +                                                          &
     &          ( (EMOTION(1)%Vx + EMOTION(2)%Vx)                              &
     &          - (EMOTION(3)%Vx + EMOTION(4)%Vx)                              &
     &          - (EMOTION(5)%Vx + EMOTION(6)%Vx)                              &
     &          + (EMOTION(7)%Vx + EMOTION(8)%Vx) )*Delta
      Gx(2) = Gx(2) +                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          - (EMOTION(3)%Vx - EMOTION(4)%Vx)                              &
     &          - (EMOTION(5)%Vx - EMOTION(6)%Vx)                              &
     &          + (EMOTION(7)%Vx - EMOTION(8)%Vx) )*Delta
      Gx(3) = Gx(3) +                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          + (EMOTION(3)%Vx - EMOTION(4)%Vx)                              &
     &          + (EMOTION(5)%Vx - EMOTION(6)%Vx)                              &
     &          + (EMOTION(7)%Vx - EMOTION(8)%Vx) )*Delta
      Gx(4) = Gx(4) -                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          + (EMOTION(3)%Vx - EMOTION(4)%Vx)                              &
     &          - (EMOTION(5)%Vx - EMOTION(6)%Vx)                              &
     &          - (EMOTION(7)%Vx - EMOTION(8)%Vx) )*Delta
!!
      Gy(1) = Gy(1) +                                                          &
     &          ( (EMOTION(1)%Vy + EMOTION(2)%Vy)                              &
     &          - (EMOTION(3)%Vy + EMOTION(4)%Vy)                              &
     &          - (EMOTION(5)%Vy + EMOTION(6)%Vy)                              &
     &          + (EMOTION(7)%Vy + EMOTION(8)%Vy) )*Delta
      Gy(2) = Gy(2) +                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          - (EMOTION(3)%Vy - EMOTION(4)%Vy)                              &
     &          - (EMOTION(5)%Vy - EMOTION(6)%Vy)                              &
     &          + (EMOTION(7)%Vy - EMOTION(8)%Vy) )*Delta
      Gy(3) = Gy(3) +                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          + (EMOTION(3)%Vy - EMOTION(4)%Vy)                              &
     &          + (EMOTION(5)%Vy - EMOTION(6)%Vy)                              &
     &          + (EMOTION(7)%Vy - EMOTION(8)%Vy) )*Delta
      Gy(4) = Gy(4) -                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          + (EMOTION(3)%Vy - EMOTION(4)%Vy)                              &
     &          - (EMOTION(5)%Vy - EMOTION(6)%Vy)                              &
     &          - (EMOTION(7)%Vy - EMOTION(8)%Vy) )*Delta
!!
      Gz(1) = Gz(1) +                                                          &
     &          ( (EMOTION(1)%Vz + EMOTION(2)%Vz)                              &
     &          - (EMOTION(3)%Vz + EMOTION(4)%Vz)                              &
     &          - (EMOTION(5)%Vz + EMOTION(6)%Vz)                              &
     &          + (EMOTION(7)%Vz + EMOTION(8)%Vz) )*Delta
      Gz(2) = Gz(2) +                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          - (EMOTION(3)%Vz - EMOTION(4)%Vz)                              &
     &          - (EMOTION(5)%Vz - EMOTION(6)%Vz)                              &
     &          + (EMOTION(7)%Vz - EMOTION(8)%Vz) )*Delta
      Gz(3) = Gz(3) +                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          + (EMOTION(3)%Vz - EMOTION(4)%Vz)                              &
     &          + (EMOTION(5)%Vz - EMOTION(6)%Vz)                              &
     &          + (EMOTION(7)%Vz - EMOTION(8)%Vz) )*Delta
      Gz(4) = Gz(4) -                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          + (EMOTION(3)%Vz - EMOTION(4)%Vz)                              &
     &          - (EMOTION(5)%Vz - EMOTION(6)%Vz)                              &
     &          - (EMOTION(7)%Vz - EMOTION(8)%Vz) )*Delta
!!
      RETURN
      END
!!_
      SUBROUTINE HEXAH_DIVERGENCE_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE hexah_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          I2(8),I3(8),I4(8),I5(8),I6(8),I8(8)
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,4,1,8,5,6,7/, I3 /3,4,1,2,7,8,5,6/,                    &
     &          I4 /4,1,2,3,6,7,8,5/, I5 /5,6,7,8,1,2,3,4/,                    &
     &          I6 /6,7,8,5,4,1,2,3/, I8 /8,5,6,7,2,3,4,1/
!!
!! Current position of nodal points, time n+1.
!!
      DO i = 1,8
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,8
        Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))               &
     &             +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))             &
     &             +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))             &
     &             +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                               &
     &             +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                               &
     &             +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))               &
     &             +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))             &
     &             +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))             &
     &             +Z(I3(i))*(X(I2(i))-X(I4(i)))                               &
     &             +Z(I6(i))*(X(I5(i))-X(I2(i)))                               &
     &             +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))               &
     &             +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))             &
     &             +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))             &
     &             +X(I3(i))*(Y(I2(i))-Y(I4(i)))                               &
     &             +X(I6(i))*(Y(I5(i))-Y(I2(i)))                               &
     &             +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      HEXAH(NEL)%RES%Volume = 0.0
      DO i = 1,8
        HEXAH(NEL)%RES%Volume = HEXAH(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0 / HEXAH(NEL)%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,8
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
!!
!! Anti-hourglass divergence operator.
!!
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(4)) -(X(5)+X(6)) +(X(7)+X(8))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(4)) -(X(5)-X(6)) +(X(7)-X(8))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(4)) -(Y(5)+Y(6)) +(Y(7)+Y(8))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(4)) -(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(4)) -(Z(5)+Z(6)) +(Z(7)+Z(8))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(4)) -(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
!!
      RETURN
      END
!!_
      SUBROUTINE HEXAH_MASS ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 19:47:17
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one eighth at each node.
!!
      USE shared_common_data
      USE hexah_
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
!! Compute one eighth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = 0.125 * Density * HEXAH(NEL)%RES%Volume
!!
!! Accumulate mass at each nodal point and accumulate mass properties
!! for each material domain.
!!
      DO i = 1,8
        NODE(HEXAH(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(HEXAH(NEL)%PAR%IX(i))%Mass + QMass
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
      SUBROUTINE HEXAH_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 20:03:56
!!
      USE shared_common_data
      USE hexah_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Fx(8),Fy(8),Fz(8),PGx(4),PGy(4),PGz(4)
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
     &          Density * HEXAH(NEL)%PAR%Volume / HEXAH(NEL)%RES%Volume
!!
      CSQ = MAX(SOUND_SPEED%RCL2,SOUND_SPEED%RCS2)/SOUND_SPEED%Density
      IF (CSQ .GT. 0.0) THEN
        CX = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') HEXAH(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'HEXAH_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'HXEL (8-Node Hexahedron) Element ID:'//MSG1//           &
     &          MSGL//'Sound Speed Imaginary (C**2 Is Zero Or Negative)'       &
     &          )
        CX = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Calculate generalized element dimension.
!!
      Dx = 1.0 / Delta
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
!! Hourglass stiffness forces combined with hourglass viscosity forces.
!! Note: The inverse of the generalized element dimension used here is
!! based on the geometry at time n+1.
!!
      Qhg = SOUND_SPEED%Density * Dx * HG_Visc * CX
      DO i = 1,4
        PGx(i) = (HEXAH(NEL)%RES%Px(i) + Qhg*Gx(i)) * Delta
        PGy(i) = (HEXAH(NEL)%RES%Py(i) + Qhg*Gy(i)) * Delta
        PGz(i) = (HEXAH(NEL)%RES%Pz(i) + Qhg*Gz(i)) * Delta
      ENDDO
!!
!! Critical time step calculation.
!!
      HEXAH(NEL)%RES%DTelt = Dx / (QC + SQRT (QC*QC + CX*CX))
!!
!! Divergence of the hourglass forces.
!!
      Fx(1) = ( (PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(2) = ( (PGx(1)-PGx(2))-(PGx(3)-PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(3) = (-(PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(4) = (-(PGx(1)-PGx(2))-(PGx(3)-PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(5) = (-(PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(6) = (-(PGx(1)-PGx(2))-(PGx(3)+PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(7) = ( (PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * HEXAH(NEL)%RES%Volume
      Fx(8) = ( (PGx(1)-PGx(2))-(PGx(3)+PGx(4))) * HEXAH(NEL)%RES%Volume
!!
      Fy(1) = ( (PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(2) = ( (PGy(1)-PGy(2))-(PGy(3)-PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(3) = (-(PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(4) = (-(PGy(1)-PGy(2))-(PGy(3)-PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(5) = (-(PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(6) = (-(PGy(1)-PGy(2))-(PGy(3)+PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(7) = ( (PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * HEXAH(NEL)%RES%Volume
      Fy(8) = ( (PGy(1)-PGy(2))-(PGy(3)+PGy(4))) * HEXAH(NEL)%RES%Volume
!!
      Fz(1) = ( (PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(2) = ( (PGz(1)-PGz(2))-(PGz(3)-PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(3) = (-(PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(4) = (-(PGz(1)-PGz(2))-(PGz(3)-PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(5) = (-(PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(6) = (-(PGz(1)-PGz(2))-(PGz(3)+PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(7) = ( (PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * HEXAH(NEL)%RES%Volume
      Fz(8) = ( (PGz(1)-PGz(2))-(PGz(3)+PGz(4))) * HEXAH(NEL)%RES%Volume
!!
!! Divergence of the (mean) stress plus the orthogonalized correction Hi(1:4).
!!
      Qxx = HEXAH(NEL)%RES%Stress(1) + QP
      Qyx = HEXAH(NEL)%RES%Stress(4)
      Qzx = HEXAH(NEL)%RES%Stress(5)
      DO j = 1,4
        Qxx = Qxx - Hx(j)*PGx(j)
        Qyx = Qyx - Hy(j)*PGx(j)
        Qzx = Qzx - Hz(j)*PGx(j)
      ENDDO
!!
      Qxy = HEXAH(NEL)%RES%Stress(4)
      Qyy = HEXAH(NEL)%RES%Stress(2) + QP
      Qzy = HEXAH(NEL)%RES%Stress(6)
      DO j = 1,4
        Qxy = Qxy - Hx(j)*PGy(j)
        Qyy = Qyy - Hy(j)*PGy(j)
        Qzy = Qzy - Hz(j)*PGy(j)
      ENDDO
!!
      Qxz = HEXAH(NEL)%RES%Stress(5)
      Qyz = HEXAH(NEL)%RES%Stress(6)
      Qzz = HEXAH(NEL)%RES%Stress(3) + QP
      DO j = 1,4
        Qxz = Qxz - Hx(j)*PGz(j)
        Qyz = Qyz - Hy(j)*PGz(j)
        Qzz = Qzz - Hz(j)*PGz(j)
      ENDDO
!!
!! Accumulate element divergence results.
!!
      DO i = 1,8
        HEXAH(NEL)%RES%Xint(i) = Fx(i) + Bx(i)*Qxx+By(i)*Qyx+Bz(i)*Qzx
        HEXAH(NEL)%RES%Yint(i) = Fy(i) + Bx(i)*Qxy+By(i)*Qyy+Bz(i)*Qzy
        HEXAH(NEL)%RES%Zint(i) = Fz(i) + Bx(i)*Qxz+By(i)*Qyz+Bz(i)*Qzz
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE HEXAH_HOURGLASS_FORCES ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 26-MAY-1991 13:19:09
!!
!! Purpose: Increment stiffness based hourglass control forces.
!!
      USE shared_common_data
      USE hexah_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Q1(4),Q2(4),Q3(4)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Rotation of hour glass forces based on spin, matches Jaumann stress
!! flux.
!!
      IF (CONTROL%POLARD .EQ. 0) THEN
        dWxy = HEXAH(NEL)%RES%DTnext * Wxy
        dWxz = HEXAH(NEL)%RES%DTnext * Wxz
        dWyz = HEXAH(NEL)%RES%DTnext * Wyz
        DO i = 1,4
          Q1(i) =  dWxy*HEXAH(NEL)%RES%Py(i) + dWxz*HEXAH(NEL)%RES%Pz(i)
          Q2(i) = -dWxy*HEXAH(NEL)%RES%Px(i) + dWyz*HEXAH(NEL)%RES%Pz(i)
          Q3(i) = -dWxz*HEXAH(NEL)%RES%Px(i) - dWyz*HEXAH(NEL)%RES%Py(i)
        ENDDO
      ELSE
        DO i = 1,4
          Q1(i) = 0.0
          Q2(i) = 0.0
          Q3(i) = 0.0
        ENDDO
      ENDIF
!!
!! Increment elastic anti-hourglassing forces.
!!
        HG_Stiff = MATERIAL(MatID)%PVAL(5)
        Qe = HEXAH(NEL)%RES%DTnext * HG_Stiff * SOUND_SPEED%RCS2
        DO i = 1,4
          HEXAH(NEL)%RES%Px(i) = HEXAH(NEL)%RES%Px(i) + Qe*Gx(i) + Q1(i)
          HEXAH(NEL)%RES%Py(i) = HEXAH(NEL)%RES%Py(i) + Qe*Gy(i) + Q2(i)
          HEXAH(NEL)%RES%Pz(i) = HEXAH(NEL)%RES%Pz(i) + Qe*Gz(i) + Q3(i)
        ENDDO
!!
!! Rotation of hour glass forces based on rotation from polar decomposi-
!! tion of the deformation gradient, matches Green-McInnes stress flux.
!! This is a rotation forward from time n to the end of the interval at
!! time n+1 and thus, uses the inverse of the rotation R  generated in
!! SYMMETRIC_TENSOR_ROTATION That is, Rt * P is used in place of R * P,
!! the transformation we would have used if we had been able to construct
!! the "forward" deformation gradient in SYMMETRIC_TENSOR_ROTATION.
!!
      IF (CONTROL%POLARD .GT. 0) THEN
        DO i = 1,4
          Q1(i) = Rotation(1,1) * HEXAH(NEL)%RES%Px(i)                         &
     &          + Rotation(2,1) * HEXAH(NEL)%RES%Py(i)                         &
     &          + Rotation(3,1) * HEXAH(NEL)%RES%Pz(i)

          Q2(i) = Rotation(1,2) * HEXAH(NEL)%RES%Px(i)                         &
     &          + Rotation(2,2) * HEXAH(NEL)%RES%Py(i)                         &
     &          + Rotation(3,2) * HEXAH(NEL)%RES%Pz(i)

          Q3(i) = Rotation(1,3) * HEXAH(NEL)%RES%Px(i)                         &
     &          + Rotation(2,3) * HEXAH(NEL)%RES%Py(i)                         &
     &          + Rotation(3,3) * HEXAH(NEL)%RES%Pz(i)
        ENDDO
        DO i = 1,4
          HEXAH(NEL)%RES%Px(i) = Q1(i)
          HEXAH(NEL)%RES%Py(i) = Q2(i)
          HEXAH(NEL)%RES%Pz(i) = Q3(i)
        ENDDO
      ENDIF
!!
      RETURN
      END
