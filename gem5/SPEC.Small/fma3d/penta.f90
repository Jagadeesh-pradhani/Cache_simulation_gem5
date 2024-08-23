      SUBROUTINE PENTA_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 12-FEB-1991 20:22:44
!!
      USE shared_common_data
      USE penta_
      USE material_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: Isv,MatID
      LOGICAL :: FOUND
!!
!! Loop over all pentahedral elements.
!!
      DO N = 1,NUMPX
!!
!! Gather element motion (position and velocity at time equal to zero).
!!
        DO i = 1,6
          EMOTION(i) = MOTION(PENTA(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        PENTA(N)%RES%Time   = 0.0
        PENTA(N)%RES%DTnext = 0.0
!!
!! Initialize hourglass restoring forces.
!!
        DO i = 1,4
          PENTA(N)%RES%Px(i) = 0
          PENTA(N)%RES%Py(i) = 0
          PENTA(N)%RES%Pz(i) = 0
        ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
        Isv = PENTA(N)%PAR%Isv
        MatID = PENTA(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies, only the
!! element volume calculation is required in order to construct the mass
!! matrix.)
!!
        CALL PENTA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Save initial element volume for later volume strain calculations.
!!
        PENTA(N)%PAR%Volume = PENTA(N)%RES%Volume

        IF (PENTA(N)%RES%Volume .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') PENTA(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PENTA_INITIALIZATION.001.00'//                          &
     &          MSGL//'PXEL (6-Node Pentahedron) Element ID:'//MSG1//          &
     &          MSGL//'Has A Zero Or Negative Volume.'                         &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL PENTA_MASS ( NEL,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (PENTA(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,6
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &              + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                      &
     &              + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          PENTA(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
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
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (31)
            CALL MATERIAL_31                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (32)
            CALL MATERIAL_32                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (33)
            CALL MATERIAL_33                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (34)
            CALL MATERIAL_34                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (35)
            CALL MATERIAL_35                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (36)
            CALL MATERIAL_36                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (37)
            CALL MATERIAL_37                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (38)
            CALL MATERIAL_38                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (39)
            CALL MATERIAL_39                                                   &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') PENTA(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PENTA_INITIALIZATION.002.00'//                          &
     &          MSGL//'PXEL (6-Node pentahedron) Element ID:'//MSG1//          &
     &          MSGL//'References An Unknown Material Model:'//MSG2            &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL PENTA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
          IF (NUMNR .NE. 0) THEN
            Itype = 1  !  Pentahedron
            CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
          ENDIF
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTPnx = MAX (TIMSIM%DTPnx,PENTA(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = PENTA(N)%RES%DTelt .LT. TIMSIM%DTPen(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTPen(j + 1) = TIMSIM%DTPen(j)
              TIMSIM%PENTA(j + 1) = TIMSIM%PENTA(j)
            ENDDO
          ENDIF
          TIMSIM%DTPen(i) = PENTA(N)%RES%DTelt
          TIMSIM%PENTA(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PENTA_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE penta_
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
      REAL(KIND(0D0))                                                          &
     &          QA
!!
!! Loop over all pentahedral elements.
!!
      DO N = 1,NUMPX
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,PENTA(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (PENTA(N)%PAR%ParID .GE. 0) THEN
!!
!! Count element execution.
!!
            COUNTER%PENTA = COUNTER%PENTA + 1
!!
!! Gather element coordinates and motion.
!!
            DO i = 1,6
              EMOTION(i) = MOTION(PENTA(N)%PAR%IX(i))
            ENDDO
!!
!! Increment element clock.
!!
            PENTA(N)%RES%Time = PENTA(N)%RES%Time + PENTA(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            DO i = 1,6
              QA = NODE(PENTA(N)%PAR%IX(i))%Time - PENTA(N)%RES%Time
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
            Isv = PENTA(N)%PAR%Isv
            MatID = PENTA(N)%PAR%MatID
!!
!! Gradient operator, stretching, and rotation evaluated at time n-1/2.
!!
            CALL PENTA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Incremental constitutive model evaluation.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (30)
              CALL MATERIAL_30                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (31)
              CALL MATERIAL_31                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (32)
              CALL MATERIAL_32                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (33)
              CALL MATERIAL_33                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (34)
              CALL MATERIAL_34                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (35)
              CALL MATERIAL_35                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (36)
              CALL MATERIAL_36                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (37)
              CALL MATERIAL_37                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (38)
              CALL MATERIAL_38                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (39)
              CALL MATERIAL_39                                                 &
     &          (                                                              &
     &          PENTA(N)%RES%STRESS,                                           &
     &          PENTA(N)%RES%Int_Eng,                                          &
     &          STATE_VARIABLES(Isv),                                          &
     &          PENTA(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
            IF (MATERIAL(MatID)%PVAL(5) .GT. 0.0) THEN
              CALL PENTA_HOURGLASS_FORCES ( NEL,MatID )
            ENDIF
!!
!! Divergence operator evaluated at time n+1.
!!
!!!           IF (CONTROL%MIDINT .NE. 0) THEN
!!!             CALL PENTA_DIVERGENCE_OPERATOR ( NEL,MatID )
!!!           ENDIF
!!
!! Update stress divergence, viscosity stress and time step.
!!
            CALL PENTA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
            IF (NUMNR .NE. 0) THEN
              Itype = 1  !  Pentahedron
              CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
            ENDIF
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTPnx = MAX (TIMSIM%DTPnx,PENTA(N)%RES%DTelt)
          IF (PENTA(N)%RES%DTelt .LT. TIMSIM%DTPen(1)) THEN
            TIMSIM%DTPen(1) = PENTA(N)%RES%DTelt
            TIMSIM%PENTA(1) = N
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
      SUBROUTINE PENTA_GRADIENT_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE penta_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          I2(6),I3(6),I4(6),I5(6),I6(6)
      REAL(KIND(0D0))                                                          &
     &          X(6),Y(6),Z(6)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,1,6,4,5/, I3 /3,1,2,5,6,4/, I4 /4,5,6,1,2,3/,          &
     &          I5 /5,6,4,3,1,2/, I6 /6,4,5,2,3,1/
!!
!! Compute position of nodal points at time n+1 (whole-interval evaluation).
!!
      DO i = 1,6
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,6
        Bx(i) = ( Y(I2(i))*(Z(I4(i))-Z(I3(i))+Z(I5(i))-Z(I3(i)))               &
     &             +Y(I3(i))*(Z(I2(i))-Z(I4(i))+Z(I2(i))-Z(I6(i)))             &
     &             +Y(I4(i))*(Z(I3(i))-Z(I2(i))+Z(I6(i))-Z(I5(i)))             &
     &             +Y(I5(i))*(Z(I4(i))-Z(I2(i)))                               &
     &             +Y(I6(i))*(Z(I3(i))-Z(I4(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,6
        By(i) = ( Z(I2(i))*(X(I4(i))-X(I3(i))+X(I5(i))-X(I3(i)))               &
     &             +Z(I3(i))*(X(I2(i))-X(I4(i))+X(I2(i))-X(I6(i)))             &
     &             +Z(I4(i))*(X(I3(i))-X(I2(i))+X(I6(i))-X(I5(i)))             &
     &             +Z(I5(i))*(X(I4(i))-X(I2(i)))                               &
     &             +Z(I6(i))*(X(I3(i))-X(I4(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,6
        Bz(i) = ( X(I2(i))*(Y(I4(i))-Y(I3(i))+Y(I5(i))-Y(I3(i)))               &
     &             +X(I3(i))*(Y(I2(i))-Y(I4(i))+Y(I2(i))-Y(I6(i)))             &
     &             +X(I4(i))*(Y(I3(i))-Y(I2(i))+Y(I6(i))-Y(I5(i)))             &
     &             +X(I5(i))*(Y(I4(i))-Y(I2(i)))                               &
     &             +X(I6(i))*(Y(I3(i))-Y(I4(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      PENTA(NEL)%RES%Volume = 0.0
      DO i = 1,6
        PENTA(NEL)%RES%Volume = PENTA(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0 / PENTA(NEL)%RES%Volume
      Hin = 0.5 * Vin
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,6
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
      DO i = 1,6
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
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(1)) -(X(4)+X(5)) +(X(6)+X(4))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(1)) -(X(4)-X(5)) +(X(6)-X(4))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(1)) +(X(4)-X(5)) +(X(6)-X(4))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(1)) +(X(4)-X(5)) +(X(6)-X(4))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(1)) -(Y(4)+Y(5)) +(Y(6)+Y(4))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(1)) -(Y(4)-Y(5)) +(Y(6)-Y(4))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(1)) +(Y(4)-Y(5)) +(Y(6)-Y(4))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(1)) +(Y(4)-Y(5)) +(Y(6)-Y(4))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(1)) -(Z(4)+Z(5)) +(Z(6)+Z(4))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(1)) -(Z(4)-Z(5)) +(Z(6)-Z(4))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(1)) +(Z(4)-Z(5)) +(Z(6)-Z(4))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(1)) +(Z(4)-Z(5)) +(Z(6)-Z(4))
!!
      DO i = 1,4
        Gx(i) = -(Vxx*Hx(i) + Vxy*Hy(i) + Vxz*Hz(i)) * Delta
        Gy(i) = -(Vyx*Hx(i) + Vyy*Hy(i) + Vyz*Hz(i)) * Delta
        Gz(i) = -(Vzx*Hx(i) + Vzy*Hy(i) + Vzz*Hz(i)) * Delta
      ENDDO
      Gx(1) = Gx(1) +                                                          &
     &          ( (EMOTION(1)%Vx + EMOTION(2)%Vx)                              &
     &          - (EMOTION(3)%Vx + EMOTION(1)%Vx)                              &
     &          - (EMOTION(4)%Vx + EMOTION(5)%Vx)                              &
     &          + (EMOTION(6)%Vx + EMOTION(4)%Vx) )*Delta
      Gx(2) = Gx(2) +                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          - (EMOTION(3)%Vx - EMOTION(1)%Vx)                              &
     &          - (EMOTION(4)%Vx - EMOTION(5)%Vx)                              &
     &          + (EMOTION(6)%Vx - EMOTION(4)%Vx) )*Delta
      Gx(3) = Gx(3) +                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          + (EMOTION(3)%Vx - EMOTION(1)%Vx)                              &
     &          + (EMOTION(4)%Vx - EMOTION(5)%Vx)                              &
     &          + (EMOTION(6)%Vx - EMOTION(4)%Vx) )*Delta
      Gx(4) = Gx(4) -                                                          &
     &          ( (EMOTION(1)%Vx - EMOTION(2)%Vx)                              &
     &          + (EMOTION(3)%Vx - EMOTION(1)%Vx)                              &
     &          - (EMOTION(4)%Vx - EMOTION(5)%Vx)                              &
     &          - (EMOTION(6)%Vx - EMOTION(4)%Vx) )*Delta
!!
      Gy(1) = Gy(1) +                                                          &
     &          ( (EMOTION(1)%Vy + EMOTION(2)%Vy)                              &
     &          - (EMOTION(3)%Vy + EMOTION(1)%Vy)                              &
     &          - (EMOTION(4)%Vy + EMOTION(5)%Vy)                              &
     &          + (EMOTION(6)%Vy + EMOTION(4)%Vy) )*Delta
      Gy(2) = Gy(2) +                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          - (EMOTION(3)%Vy - EMOTION(1)%Vy)                              &
     &          - (EMOTION(4)%Vy - EMOTION(5)%Vy)                              &
     &          + (EMOTION(6)%Vy - EMOTION(4)%Vy) )*Delta
      Gy(3) = Gy(3) +                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          + (EMOTION(3)%Vy - EMOTION(1)%Vy)                              &
     &          + (EMOTION(4)%Vy - EMOTION(5)%Vy)                              &
     &          + (EMOTION(6)%Vy - EMOTION(4)%Vy) )*Delta
      Gy(4) = Gy(4) -                                                          &
     &          ( (EMOTION(1)%Vy - EMOTION(2)%Vy)                              &
     &          + (EMOTION(3)%Vy - EMOTION(1)%Vy)                              &
     &          - (EMOTION(4)%Vy - EMOTION(5)%Vy)                              &
     &          - (EMOTION(6)%Vy - EMOTION(4)%Vy) )*Delta
!!
      Gz(1) = Gz(1) +                                                          &
     &          ( (EMOTION(1)%Vz + EMOTION(2)%Vz)                              &
     &          - (EMOTION(3)%Vz + EMOTION(1)%Vz)                              &
     &          - (EMOTION(4)%Vz + EMOTION(5)%Vz)                              &
     &          + (EMOTION(6)%Vz + EMOTION(4)%Vz) )*Delta
      Gz(2) = Gz(2) +                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          - (EMOTION(3)%Vz - EMOTION(1)%Vz)                              &
     &          - (EMOTION(4)%Vz - EMOTION(5)%Vz)                              &
     &          + (EMOTION(6)%Vz - EMOTION(4)%Vz) )*Delta
      Gz(3) = Gz(3) +                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          + (EMOTION(3)%Vz - EMOTION(1)%Vz)                              &
     &          + (EMOTION(4)%Vz - EMOTION(5)%Vz)                              &
     &          + (EMOTION(6)%Vz - EMOTION(4)%Vz) )*Delta
      Gz(4) = Gz(4) -                                                          &
     &          ( (EMOTION(1)%Vz - EMOTION(2)%Vz)                              &
     &          + (EMOTION(3)%Vz - EMOTION(1)%Vz)                              &
     &          - (EMOTION(4)%Vz - EMOTION(5)%Vz)                              &
     &          - (EMOTION(6)%Vz - EMOTION(4)%Vz) )*Delta
!!
      RETURN
      END
!!_
      SUBROUTINE PENTA_DIVERGENCE_OPERATOR ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE penta_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          I2(6),I3(6),I4(6),I5(6),I6(6)
      REAL(KIND(0D0))                                                          &
     &          X(6),Y(6),Z(6)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,1,6,4,5/, I3 /3,1,2,5,6,4/, I4 /4,5,6,1,2,3/,          &
     &          I5 /5,6,4,3,1,2/, I6 /6,4,5,2,3,1/
!!
!! Current position of nodal points, time n+1.
!!
      DO i = 1,6
        X(i) = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i) = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i) = EMOTION(i)%Pz + EMOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,6
        Bx(i) = ( Y(I2(i))*(Z(I4(i))-Z(I3(i))+Z(I5(i))-Z(I3(i)))               &
     &             +Y(I3(i))*(Z(I2(i))-Z(I4(i))+Z(I2(i))-Z(I6(i)))             &
     &             +Y(I4(i))*(Z(I3(i))-Z(I2(i))+Z(I6(i))-Z(I5(i)))             &
     &             +Y(I5(i))*(Z(I4(i))-Z(I2(i)))                               &
     &             +Y(I6(i))*(Z(I3(i))-Z(I4(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,6
        By(i) = ( Z(I2(i))*(X(I4(i))-X(I3(i))+X(I5(i))-X(I3(i)))               &
     &             +Z(I3(i))*(X(I2(i))-X(I4(i))+X(I2(i))-X(I6(i)))             &
     &             +Z(I4(i))*(X(I3(i))-X(I2(i))+X(I6(i))-X(I5(i)))             &
     &             +Z(I5(i))*(X(I4(i))-X(I2(i)))                               &
     &             +Z(I6(i))*(X(I3(i))-X(I4(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,6
        Bz(i) = ( X(I2(i))*(Y(I4(i))-Y(I3(i))+Y(I5(i))-Y(I3(i)))               &
     &             +X(I3(i))*(Y(I2(i))-Y(I4(i))+Y(I2(i))-Y(I6(i)))             &
     &             +X(I4(i))*(Y(I3(i))-Y(I2(i))+Y(I6(i))-Y(I5(i)))             &
     &             +X(I5(i))*(Y(I4(i))-Y(I2(i)))                               &
     &             +X(I6(i))*(Y(I3(i))-Y(I4(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      PENTA(NEL)%RES%Volume = 0.0
      DO i = 1,6
        PENTA(NEL)%RES%Volume = PENTA(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0 / PENTA(NEL)%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,6
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
!!
!! Anti-hourglass divergence operator.
!!
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(1)) -(X(4)+X(5)) +(X(6)+X(4))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(1)) -(X(4)-X(5)) +(X(6)-X(4))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(1)) +(X(4)-X(5)) +(X(6)-X(4))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(1)) +(X(4)-X(5)) +(X(6)-X(4))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(1)) -(Y(4)+Y(5)) +(Y(6)+Y(4))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(1)) -(Y(4)-Y(5)) +(Y(6)-Y(4))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(1)) +(Y(4)-Y(5)) +(Y(6)-Y(4))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(1)) +(Y(4)-Y(5)) +(Y(6)-Y(4))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(1)) -(Z(4)+Z(5)) +(Z(6)+Z(4))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(1)) -(Z(4)-Z(5)) +(Z(6)-Z(4))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(1)) +(Z(4)-Z(5)) +(Z(6)-Z(4))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(1)) +(Z(4)-Z(5)) +(Z(6)-Z(4))
!!
      RETURN
      END
!!_
      SUBROUTINE PENTA_MASS ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 19:47:17
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one sixth at each node.
!!
      USE shared_common_data
      USE penta_
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
!! Compute one sixth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = 0.1666667 * Density * PENTA(NEL)%RES%Volume
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,6
        NODE(PENTA(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(PENTA(NEL)%PAR%IX(i))%Mass + QMass
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
      SUBROUTINE PENTA_STRESS_DIVERGENCE ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 20:03:56
!!
      USE shared_common_data
      USE penta_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Fx(6),Fy(6),Fz(6),PGx(4),PGy(4),PGz(4)
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
     &          Density * PENTA(NEL)%PAR%Volume / PENTA(NEL)%RES%Volume
!!
      CSQ = MAX(SOUND_SPEED%RCL2,SOUND_SPEED%RCS2) / SOUND_SPEED%Density
      IF (CSQ .GT. 0.0) THEN
        CX = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') PENTA(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PENTA_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'PXEL (6-Node Pentahedron) Element ID:'//MSG1//          &
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
        PGx(i) = (PENTA(NEL)%RES%Px(i) + Qhg*Gx(i)) * Delta
        PGy(i) = (PENTA(NEL)%RES%Py(i) + Qhg*Gy(i)) * Delta
        PGz(i) = (PENTA(NEL)%RES%Pz(i) + Qhg*Gz(i)) * Delta
      ENDDO
!!
!! Critical time step calculation.
!!
      PENTA(NEL)%RES%DTelt = Dx / (QC + SQRT (QC*QC + CX*CX))
!!
!! Divergence of the hourglass forces.
!!
      Fx(1) = ( (PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * PENTA(NEL)%RES%Volume
      Fx(2) = ( (PGx(1)-PGx(2))-(PGx(3)-PGx(4))) * PENTA(NEL)%RES%Volume
      Fx(3) = (-(PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * PENTA(NEL)%RES%Volume
      Fx(4) = (-(PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * PENTA(NEL)%RES%Volume
      Fx(5) = (-(PGx(1)-PGx(2))-(PGx(3)+PGx(4))) * PENTA(NEL)%RES%Volume
      Fx(6) = ( (PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * PENTA(NEL)%RES%Volume
!!
      Fy(1) = ( (PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * PENTA(NEL)%RES%Volume
      Fy(2) = ( (PGy(1)-PGy(2))-(PGy(3)-PGy(4))) * PENTA(NEL)%RES%Volume
      Fy(3) = (-(PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * PENTA(NEL)%RES%Volume
      Fy(4) = (-(PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * PENTA(NEL)%RES%Volume
      Fy(5) = (-(PGy(1)-PGy(2))-(PGy(3)+PGy(4))) * PENTA(NEL)%RES%Volume
      Fy(6) = ( (PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * PENTA(NEL)%RES%Volume
!!
      Fz(1) = ( (PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * PENTA(NEL)%RES%Volume
      Fz(2) = ( (PGz(1)-PGz(2))-(PGz(3)-PGz(4))) * PENTA(NEL)%RES%Volume
      Fz(3) = (-(PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * PENTA(NEL)%RES%Volume
      Fz(4) = (-(PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * PENTA(NEL)%RES%Volume
      Fz(5) = (-(PGz(1)-PGz(2))-(PGz(3)+PGz(4))) * PENTA(NEL)%RES%Volume
      Fz(6) = ( (PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * PENTA(NEL)%RES%Volume
!!
!! Divergence of the (mean) stress plus the orthogonalized correction Hi(1:4).
!!
      Qxx = PENTA(NEL)%RES%Stress(1) + QP
      Qyx = PENTA(NEL)%RES%Stress(4)
      Qzx = PENTA(NEL)%RES%Stress(5)
      DO j = 1,4
        Qxx = Qxx - Hx(j)*PGx(j)
        Qyx = Qyx - Hy(j)*PGx(j)
        Qzx = Qzx - Hz(j)*PGx(j)
      ENDDO
!!
      Qxy = PENTA(NEL)%RES%Stress(4)
      Qyy = PENTA(NEL)%RES%Stress(2) + QP
      Qzy = PENTA(NEL)%RES%Stress(6)
      DO j = 1,4
        Qxy = Qxy - Hx(j)*PGy(j)
        Qyy = Qyy - Hy(j)*PGy(j)
        Qzy = Qzy - Hz(j)*PGy(j)
      ENDDO
!!
      Qxz = PENTA(NEL)%RES%Stress(5)
      Qyz = PENTA(NEL)%RES%Stress(6)
      Qzz = PENTA(NEL)%RES%Stress(3) + QP
      DO j = 1,4
        Qxz = Qxz - Hx(j)*PGz(j)
        Qyz = Qyz - Hy(j)*PGz(j)
        Qzz = Qzz - Hz(j)*PGz(j)
      ENDDO
!!
!! Accumulate element divergence results.
!!
      DO i = 1,6
        PENTA(NEL)%RES%Xint(i) = Fx(i) + Bx(i)*Qxx+By(i)*Qyx+Bz(i)*Qzx
        PENTA(NEL)%RES%Yint(i) = Fy(i) + Bx(i)*Qxy+By(i)*Qyy+Bz(i)*Qzy
        PENTA(NEL)%RES%Zint(i) = Fz(i) + Bx(i)*Qxz+By(i)*Qyz+Bz(i)*Qzz
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PENTA_HOURGLASS_FORCES ( NEL,MatID )
!!
!! Copyright (c) by KEY Associates, 26-MAY-1991 13:19:09
!!
!! Purpose: Increment stiffness based hourglass control forces.
!!
      USE shared_common_data
      USE penta_
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
        dWxy = PENTA(NEL)%RES%DTnext * Wxy
        dWxz = PENTA(NEL)%RES%DTnext * Wxz
        dWyz = PENTA(NEL)%RES%DTnext * Wyz
        DO i = 1,4
          Q1(i) =  dWxy*PENTA(NEL)%RES%Py(i) + dWxz*PENTA(NEL)%RES%Pz(i)
          Q2(i) = -dWxy*PENTA(NEL)%RES%Px(i) + dWyz*PENTA(NEL)%RES%Pz(i)
          Q3(i) = -dWxz*PENTA(NEL)%RES%Px(i) - dWyz*PENTA(NEL)%RES%Py(i)
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
        Qe = PENTA(NEL)%RES%DTnext * HG_Stiff * SOUND_SPEED%RCS2
        DO i = 1,4
          PENTA(NEL)%RES%Px(i) = PENTA(NEL)%RES%Px(i) + Qe*Gx(i) + Q1(i)
          PENTA(NEL)%RES%Py(i) = PENTA(NEL)%RES%Py(i) + Qe*Gy(i) + Q2(i)
          PENTA(NEL)%RES%Pz(i) = PENTA(NEL)%RES%Pz(i) + Qe*Gz(i) + Q3(i)
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
          Q1(i) = Rotation(1,1) * PENTA(NEL)%RES%Px(i)                         &
     &          + Rotation(2,1) * PENTA(NEL)%RES%Py(i)                         &
     &          + Rotation(3,1) * PENTA(NEL)%RES%Pz(i)

          Q2(i) = Rotation(1,2) * PENTA(NEL)%RES%Px(i)                         &
     &          + Rotation(2,2) * PENTA(NEL)%RES%Py(i)                         &
     &          + Rotation(3,2) * PENTA(NEL)%RES%Pz(i)

          Q3(i) = Rotation(1,3) * PENTA(NEL)%RES%Px(i)                         &
     &          + Rotation(2,3) * PENTA(NEL)%RES%Py(i)                         &
     &          + Rotation(3,3) * PENTA(NEL)%RES%Pz(i)
        ENDDO
        DO i = 1,4
          PENTA(NEL)%RES%Px(i) = Q1(i)
          PENTA(NEL)%RES%Py(i) = Q2(i)
          PENTA(NEL)%RES%Pz(i) = Q3(i)
        ENDDO
      ENDIF
!!
      RETURN
      END
