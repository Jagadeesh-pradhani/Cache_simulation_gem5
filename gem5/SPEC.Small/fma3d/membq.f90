      SUBROUTINE MEMBQ_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:57
!!
      USE shared_common_data
      USE membq_
      USE material_
      USE node_
      USE motion_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: IX,Isv,MatID,SecID
      LOGICAL :: FOUND
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
      DO N = 1,NUMM4
!!
!! Gather element motion.
!!
        DO i = 1,4
          EMOTION(i) = MOTION(MEMBQ(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        MEMBQ(N)%RES%Time   = 0.0
        MEMBQ(N)%RES%DTnext = 0.0
!!
!! Initialize hourglass restoring forces.
!!
        MEMBQ(N)%RES%Pr = 0
        MEMBQ(N)%RES%Ps = 0
        MEMBQ(N)%RES%Pt = 0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve material ID and section ID from element data structure.
!!
        MatID = MEMBQ(N)%PAR%MatID
        SecID = MEMBQ(N)%PAR%SecID
        Isv   = MEMBQ(N)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies
!! only the element volume calculation is required.)
!!
        CALL MEMBQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Save initial element area for later area strain calculations.
!!
        MEMBQ(N)%PAR%Area = MEMBQ(N)%RES%Area

        IF (MEMBQ(N)%RES%Area .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') MEMBQ(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'MEMBQ_INITIALIZATION.001.00'//                          &
     &          MSGL//'M4EL (4-Node Membrane) Element ID:'//MSG1//             &
     &          MSGL//'Has A Zero Or Negative Area.'                           &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL MEMBQ_MASS ( NEL,SecID,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (MEMBQ(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,4
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &              + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                      &
     &              + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          MEMBQ(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Find initial sound speeds and stress.
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (20)
            CALL MATERIAL_20                                                   &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (21)
            CALL MATERIAL_21                                                   &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (22)
!!
!! Convert global fiber vectors to local element coordinates.
!!
            CALL GLOBAL_TO_LOCAL                                               &
     &        (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))
!!
            CALL MATERIAL_22                                                   &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (25)
!!
!! Convert global fiber vectors to local element coordinates.
!!
            CALL GLOBAL_TO_LOCAL                                               &
     &        (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))
!!
            CALL MATERIAL_25                                                   &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (27)
            CALL MATERIAL_27                                                   &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') MEMBQ(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBQ_INITIALIZATION.002.00'//                          &
     &          MSGL//'M4EL (4-Node Membrane) Element ID:'//MSG1//             &
     &          MSGL//'References An Unknown Material Model:'//MSG2            &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL MEMBQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTM4x = MAX (TIMSIM%DTM4x,MEMBQ(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = MEMBQ(N)%RES%DTelt .LT. TIMSIM%DTMb4(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTMb4(j + 1) = TIMSIM%DTMb4(j)
              TIMSIM%Memb4(j + 1) = TIMSIM%Memb4(j)
            ENDDO
          ENDIF
          TIMSIM%DTMb4(i) = MEMBQ(N)%RES%DTelt
          TIMSIM%Memb4(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBQ_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:59
!!
      USE shared_common_data
      USE membq_
      USE material_
      USE node_
      USE motion_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          IX,                                                            &
     &          Isv,                                                           &
     &          MatID,                                                         &
     &          SecID
!!
      DO N = 1,NUMM4
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,MEMBQ(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (MEMBQ(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            DO i = 1,4
              EMOTION(i) = MOTION(MEMBQ(N)%PAR%IX(i))
            ENDDO
!!
!! Count element execution.
!!
            COUNTER%MEMBQ = COUNTER%MEMBQ + 1
!!
!! Increment element clock.
!!
            MEMBQ(N)%RES%Time = MEMBQ(N)%RES%Time + MEMBQ(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            DO i = 1,4
              QA = NODE(MEMBQ(N)%PAR%IX(i))%Time - MEMBQ(N)%RES%Time
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
!! Retrieve material ID and section ID from element data structure.
!!
            MatID = MEMBQ(N)%PAR%MatID
            SecID = MEMBQ(N)%PAR%SecID
            Isv   = MEMBQ(N)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation.
!!
            CALL MEMBQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Update sound speeds and stress.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (20)
              CALL MATERIAL_20                                                 &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (21)
              CALL MATERIAL_21                                                 &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (22)
              CALL MATERIAL_22                                                 &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (25)
              CALL MATERIAL_25                                                 &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (27)
              CALL MATERIAL_27                                                 &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBQ(N)%RES%Int_Eng,                                          &
     &          MEMBQ(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
            IF (MATERIAL(MatID)%PVAL(5) .NE. 0.0) THEN
              CALL MEMBQ_HOURGLASS_FORCES ( NEL,SecID,MatID )
            ENDIF
!!
!! Construct the divergence operator (at the end of the interval).
!!
            CALL MEMBQ_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Compute stress divergence and time step.
!!
            CALL MEMBQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTM4x = MAX (TIMSIM%DTM4x,MEMBQ(N)%RES%DTelt)
          IF (MEMBQ(N)%RES%DTelt .LT. TIMSIM%DTMb4(1)) THEN
            TIMSIM%DTMb4(1) = MEMBQ(N)%RES%DTelt
            TIMSIM%Memb4(1) = N
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
      SUBROUTINE MEMBQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
      USE shared_common_data
      USE membq_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      REAL(KIND(0D0))                                                          &
     &          Vx(4),Vy(4),Vz(4),X(4),Y(4),Z(4),                              &
     &          Vr(4),Vs(4),Vt(4),R(4),S(4),T(4)
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
!! Mid-interval position of nodal points, and current translational velocities.
!!
      DT = 0.5 * MEMBQ(NEL)%RES%DTnext
      DO i = 1,4
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! CONSTRUCT LOCAL BASIS VECTORS.
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element sides. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag  = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (1.0 / Rmag)
      Ry = Ry * (1.0 / Rmag)
      Rz = Rz * (1.0 / Rmag)
      Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (1.0 / Smag)
      Sy = Sy * (1.0 / Smag)
      Sz = Sz * (1.0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') MEMBQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBQ_GRADIENT_OPERATOR.001.00'//                       &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
     &          )
      ENDIF
!!
      Tx = Tx * (1.0 / Tmag)
      Ty = Ty * (1.0 / Tmag)
      Tz = Tz * (1.0 / Tmag)
!!
!! Redefine S to be orthogonal to T and R.
!!
      Sx = Ty*Rz - Ry*Tz
      Sy = Tz*Rx - Rz*Tx
      Sz = Tx*Ry - Rx*Ty
!!
!! Transform position X,Y,Z, and translational velocity Vx,Vy,Vz to local
!! R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
      ENDDO
!!
!! MEMBRANE GRADIENT OPERATOR.
!!
      R13 = R(1) - R(3)
      R24 = R(2) - R(4)
      S13 = S(1) - S(3)
      S24 = S(2) - S(4)
!!
      Br(1) =  (0.5 * S24)
      Br(3) = -(0.5 * S24)
      Br(4) =  (0.5 * S13)
      Br(2) = -(0.5 * S13)
      Bs(3) =  (0.5 * R24)
      Bs(1) = -(0.5 * R24)
      Bs(2) =  (0.5 * R13)
      Bs(4) = -(0.5 * R13)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      MEMBQ(NEL)%RES%Area = 0.5 * (R13*S24 - S13*R24)
      Ain = 1.0 / MEMBQ(NEL)%RES%Area
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + Br(i)*Br(i) + Bs(i)*Bs(i)
      ENDDO
      Delta = Ain * SQRT(Dx)
!!
!! VELOCITY GRADIENTS, STRETCHING, AND SPIN.
!! Construct velocity gradients.
!!
      Vrr = 0.0
      Vsr = 0.0
      Vtr = 0.0
      Vrs = 0.0
      Vss = 0.0
      Vts = 0.0
      DO i = 1,4
        Vrr = Vrr + Vr(i)*Br(i)
        Vsr = Vsr + Vs(i)*Br(i)
        Vtr = Vtr + Vt(i)*Br(i)
        Vrs = Vrs + Vr(i)*Bs(i)
        Vss = Vss + Vs(i)*Bs(i)
        Vts = Vts + Vt(i)*Bs(i)
      ENDDO
      Vrr = Ain * Vrr
      Vsr = Ain * Vsr
      Vrs = Ain * Vrs
      Vss = Ain * Vss
!!
!! Stretching and spin components in local coordinates.
!!
      Drr = Vrr
      Dss = Vss
      Drs = 0.5 * (Vrs + Vsr)
      Wrs = 0.5 * (Vrs - Vsr)
!!
!! POLAR DECOMPOSITION OF THE DEFORMATION GRADIENT.
!! Compute the polar decomposition of the deformation gradient beteen
!! time (n) and time (n+1). Construct components of the deformation
!! gradient; compute the sum of F11 and F22, and the difference of F12
!! and F21.  If r(R,S,t) and s(R,S,t), then QSF = r,R + s,S  and
!! QDF = r,S - s,R. Hin = 1.0 / (2.0*Area(0)) cancels from the quotient
!! in the expression for the inverse tangent function, ATan.
!!
      DTnext = MEMBQ(NEL)%RES%DTnext
      dBeta = ATAN (DTnext*(Vrs-Vsr)/(2.0 + DTnext*(Vrr+Vss)))
      MEMBQ(NEL)%RES%Beta = MEMBQ(NEL)%RES%Beta + dBeta
!!
!! Anti-hourglass gradients.
!!
      Hr = Ain * (R(1) - R(2) + R(3) - R(4))
      Hs = Ain * (S(1) - S(2) + S(3) - S(4))
      Ht = Ain * (T(1) - T(2) + T(3) - T(4))
!!
      Gr = ( Vr(1)-Vr(2)+Vr(3)-Vr(4)-Vrr*Hr-Vrs*Hs )*Delta
      Gs = ( Vs(1)-Vs(2)+Vs(3)-Vs(4)-Vsr*Hr-Vss*Hs )*Delta
      Gt = ( Vt(1)-Vt(2)+Vt(3)-Vt(4)-Vtr*Hr-Vts*Hs )*Delta
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBQ_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
      USE shared_common_data
      USE membq_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      REAL(KIND(0D0))                                                          &
     &          Vx(4),Vy(4),Vz(4),X(4),Y(4),Z(4),                              &
     &          Vr(4),Vs(4),Vt(4),R(4),S(4),T(4)
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
!! Current position of nodal points, and current translational velocities.
!!
      DO i = 1,4
        X(i)  = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i)  = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i)  = EMOTION(i)%Pz + EMOTION(i)%Uz
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! CONSTRUCT LOCAL BASIS VECTORS.
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element sides. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag  = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (1.0 / Rmag)
      Ry = Ry * (1.0 / Rmag)
      Rz = Rz * (1.0 / Rmag)
      Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (1.0 / Smag)
      Sy = Sy * (1.0 / Smag)
      Sz = Sz * (1.0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') MEMBQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBQ_DIVERGENCE_OPERATOR.001.00'//                     &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
     &          )
      ENDIF
!!
      Tx = Tx * (1.0 / Tmag)
      Ty = Ty * (1.0 / Tmag)
      Tz = Tz * (1.0 / Tmag)
!!
!! Redefine S to be orthogonal to T and R.
!!
      Sx = Ty*Rz - Ry*Tz
      Sy = Tz*Rx - Rz*Tx
      Sz = Tx*Ry - Rx*Ty
!!
!! Transform position X,Y,Z, and translational velocity Vx,Vy,Vz to local
!! R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
      ENDDO
!!
!! MEMBRANE GRADIENT OPERATOR.
!!
      R13 = R(1) - R(3)
      R24 = R(2) - R(4)
      S13 = S(1) - S(3)
      S24 = S(2) - S(4)
!!
      Br(1) =  (0.5 * S24)
      Br(3) = -(0.5 * S24)
      Br(4) =  (0.5 * S13)
      Br(2) = -(0.5 * S13)
      Bs(3) =  (0.5 * R24)
      Bs(1) = -(0.5 * R24)
      Bs(2) =  (0.5 * R13)
      Bs(4) = -(0.5 * R13)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      MEMBQ(NEL)%RES%Area = 0.5 * (R13*S24 - S13*R24)
      Ain = 1.0 / MEMBQ(NEL)%RES%Area
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + Br(i)*Br(i) + Bs(i)*Bs(i)
      ENDDO
      Delta = Ain * SQRT(Dx)
!!
!! Anti-hourglass gradients.
!!
      Hr = Ain * (R(1) - R(2) + R(3) - R(4))
      Hs = Ain * (S(1) - S(2) + S(3) - S(4))
      Ht = Ain * (T(1) - T(2) + T(3) - T(4))
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBQ_MASS ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:58:00
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one fourth at each node and an isotropic inertia.
!!
      USE shared_common_data
      USE membq_
      USE section_2d_
      USE material_
      USE node_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
!! Compute one fourth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = Density*MEMBQ(NEL)%PAR%Area*SECTION_2D(SecID)%Thickness/4.
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,4
        NODE(MEMBQ(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(MEMBQ(NEL)%PAR%IX(i))%Mass + QMass
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
      SUBROUTINE MEMBQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:58:01
!!
      USE shared_common_data
      USE membq_
      USE section_2d_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      REAL(KIND(0D0))                                                          &
     &          Fr(4),Fs(4),Ft(4)
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
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
      CSQ = MAX (SOUND_SPEED%RCL2,SOUND_SPEED%RCS2) / Density
      IF (CSQ .GT. 0.0) THEN
        Cv = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') MEMBQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'MEMBQ_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'//                &
     &          MSGL//'M4EL (4-Node Membrane) Element ID:'//MSG1               &
     &          )
        Cv = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Calculate generalized element dimension.
!!
      Dx = 1.0 / Delta
!!
!! Artificial bulk viscosity pressure.
!!
      Dkk = Drr + Dss
      IF (Dkk .LT. 0.0) THEN
        QC = Bulk_Qd * Bulk_Qd * Dx * ABS(Dkk) + Bulk_Ln * Cv
        QP = Density * Dx * Dkk * QC
      ELSE
        QC = 0.0
        QP = 0.0
      ENDIF
!!
!! Critical time step calculation.
!!
      MEMBQ(NEL)%RES%DTelt = Dx / (QC + SQRT (QC*QC + Cv*Cv))
!!
!! Artificial hourglass viscosity.
!!
      Qhg = Density * Dx * HG_Visc * Cv
      PGr = (MEMBQ(NEL)%RES%Pr + Qhg*Gr) * Delta
      PGs = (MEMBQ(NEL)%RES%Ps + Qhg*Gs) * Delta
      PGt = (MEMBQ(NEL)%RES%Pt + Qhg*Gt) * Delta
!!
!! Compute current element thickness based on constant volume.
!!
      Thickness = SECTION_2D(SecID)%Thickness *                                &
     &          MEMBQ(NEL)%PAR%Area / MEMBQ(NEL)%RES%Area
!!
!! Divergence of the membrane stress resultants.
!!
      Vrr = Thickness * (MEMBQ(NEL)%RES%Stress(1) + QP - Hr*PGr)
      Vsr = Thickness * (MEMBQ(NEL)%RES%Stress(3)      - Hs*PGr)
      Vrs = Thickness * (MEMBQ(NEL)%RES%Stress(3)      - Hr*PGs)
      Vss = Thickness * (MEMBQ(NEL)%RES%Stress(2) + QP - Hs*PGs)
      Vrt = Thickness * (                              - Hr*PGt)
      Vst = Thickness * (                              - Hs*PGt)
!!
      DO i = 1,4
        Fr(i) = Br(i)*Vrr + Bs(i)*Vsr
        Fs(i) = Br(i)*Vrs + Bs(i)*Vss
        Ft(i) = Br(i)*Vrt + Bs(i)*Vst
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element force array.
!!
      DO i = 1,4
        MEMBQ(NEL)%RES%Xint(i) = Rx*Fr(i) + Sx*Fs(i)+ Tx*Ft(i)
        MEMBQ(NEL)%RES%Yint(i) = Ry*Fr(i) + Sy*Fs(i)+ Ty*Ft(i)
        MEMBQ(NEL)%RES%Zint(i) = Rz*Fr(i) + Sz*Fs(i)+ Tz*Ft(i)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBQ_HOURGLASS_FORCES ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:58:03
!!
!! Purpose: Increment stiffness based hourglass control forces.
!!
      USE shared_common_data
      USE membq_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching
     &          Delta,                        & ! Generalized element size
     &          dBeta,                        & ! Incremental rotation
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
      dWrs = MEMBQ(NEL)%RES%DTnext * Wrs
      HG_Stiff = MATERIAL(MatID)%PVAL(5)
      Qmod = MEMBQ(NEL)%RES%DTnext * HG_Stiff * SOUND_SPEED%RCS2
!!
      Qr = Qmod * Gr + dWrs * MEMBQ(NEL)%RES%Ps
      Qs = Qmod * Gs - dWrs * MEMBQ(NEL)%RES%Pr
      Qt = Qmod * Gt
!!
      MEMBQ(NEL)%RES%Pr = MEMBQ(NEL)%RES%Pr + Qr
      MEMBQ(NEL)%RES%Ps = MEMBQ(NEL)%RES%Ps + Qs
      MEMBQ(NEL)%RES%Pt = MEMBQ(NEL)%RES%Pt + Qt
!!
      RETURN
      END
