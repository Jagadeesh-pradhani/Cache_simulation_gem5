      SUBROUTINE MEMBT_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:01
!!
      USE shared_common_data
      USE membt_
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
      DO N = 1,NUMM3
!!
!! Gather element motion.
!!
        DO i = 1,3
          EMOTION(i) = MOTION(MEMBT(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        MEMBT(N)%RES%Time   = 0.0
        MEMBT(N)%RES%DTnext = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve material ID and section ID from element data structure.
!!
        MatID = MEMBT(N)%PAR%MatID
        SecID = MEMBT(N)%PAR%SecID
        Isv   = MEMBT(N)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies
!! only the element area calculation is required.)
!!
        CALL MEMBT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Save initial element area for later areal strain calculations.
!!
        MEMBT(N)%PAR%Area = MEMBT(N)%RES%Area

        IF (MEMBT(N)%RES%Area .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') MEMBT(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'MEMBT_INITIALIZATION.001.00'//                          &
     &          MSGL//'M3EL (3-Node Membrane) Element ID:'//MSG1//             &
     &          MSGL//'Has A Zero Or Negative Area.'                           &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL MEMBT_MASS ( NEL,SecID,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (MEMBT(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,3
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &              + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                      &
     &              + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          MEMBT(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
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
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (21)
            CALL MATERIAL_21                                                   &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
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
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
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
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (27)
            CALL MATERIAL_27                                                   &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') MEMBT(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBT_INITIALIZATION.002.00'//                          &
     &          MSGL//'M3EL (3-Node Membrane) Element ID:'//MSG1//             &
     &          MSGL//'References An Unknown Material Model:'//MSG2            &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL MEMBT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTM3x = MAX (TIMSIM%DTM3x,MEMBT(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = MEMBT(N)%RES%DTelt .LT. TIMSIM%DTMb3(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTMb3(j + 1) = TIMSIM%DTMb3(j)
              TIMSIM%Memb3(j + 1) = TIMSIM%Memb3(j)
            ENDDO
          ENDIF
          TIMSIM%DTMb3(i) = MEMBT(N)%RES%DTelt
          TIMSIM%Memb3(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBT_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:32
!!
      USE shared_common_data
      USE membt_
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
      DO N = 1,NUMM3
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,MEMBT(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (MEMBT(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            DO i = 1,3
              EMOTION(i) = MOTION(MEMBT(N)%PAR%IX(i))
            ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Count element execution.
!!
            COUNTER%MEMBT = COUNTER%MEMBT + 1
!!
!! Increment element clock.
!!
            MEMBT(N)%RES%Time = MEMBT(N)%RES%Time + MEMBT(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            DO i = 1,3
              QA = NODE(MEMBT(N)%PAR%IX(i))%Time - MEMBT(N)%RES%Time
              EMOTION(i)%Ux =                                                  &
     &          EMOTION(i)%Ux - QA * EMOTION(i)%Vx
              EMOTION(i)%Uy =                                                  &
     &          EMOTION(i)%Uy - QA * EMOTION(i)%Vy
              EMOTION(i)%Uz =                                                  &
     &          EMOTION(i)%Uz - QA * EMOTION(i)%Vz
            ENDDO
!!
!! Retrieve material ID and section ID from element data structure.
!!
            MatID = MEMBT(N)%PAR%MatID
            SecID = MEMBT(N)%PAR%SecID
            Isv   = MEMBT(N)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation.
!!
            CALL MEMBT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Update sound speeds and stress.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (20)
              CALL MATERIAL_20                                                 &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (21)
              CALL MATERIAL_21                                                 &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (22)
              CALL MATERIAL_22                                                 &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (25)
              CALL MATERIAL_25                                                 &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (27)
            CALL MATERIAL_27                                                   &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(Isv),                                          &
     &          MEMBT(N)%RES%Int_Eng,                                          &
     &          MEMBT(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Construct the divergence operator (at the end of the interval).
!!
            CALL MEMBT_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Compute stress divergence and time step.
!!
            CALL MEMBT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTM3x = MAX (TIMSIM%DTM3x,MEMBT(N)%RES%DTelt)
          IF (MEMBT(N)%RES%DTelt .LT. TIMSIM%DTMb3(1)) THEN
            TIMSIM%DTMb3(1) = MEMBT(N)%RES%DTelt
            TIMSIM%Memb3(1) = N
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
      SUBROUTINE MEMBT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
!! Purpose: Construct triangular membrane gradient operators
!!
      USE shared_common_data
      USE membt_
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
     &          X(3),Y(3),Z(3),Vx(3),Vy(3),Vz(3),                              &
     &          R(3),S(3),Vr(3),Vs(3)
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
      DT = 0.5 * MEMBT(NEL)%RES%DTnext
      DO i = 1,3
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! CONSTRUCT LOCAL BASIS VECTORS
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element sides. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
      Rx = X(2)-X(1)
      Ry = Y(2)-Y(1)
      Rz = Z(2)-Z(1)
      Qmag  = 1.0 / SQRT (Rx*Rx+Ry*Ry+Rz*Rz)
      Rx = Rx * Qmag
      Ry = Ry * Qmag
      Rz = Rz * Qmag
      Sx = X(3)-X(1)
      Sy = Y(3)-Y(1)
      Sz = Z(3)-Z(1)
      Qmag  = 1.0 / SQRT (Sx*Sx+Sy*Sy+Sz*Sz)
      Sx = Sx * Qmag
      Sy = Sy * Qmag
      Sz = Sz * Qmag
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') MEMBT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBT_GRADIENT_OPERATOR.001.00'//                       &
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
!! Transform position X,Y,Z, and translational velocity Vx,Vy,Vz, to local
!! coordinates.
!!
      DO i = 1,3
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
      ENDDO
!!
!! MEMBRANE GRADIENT OPERATOR.
!!
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S23 = S(2) - S(3)
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R23 = R(2) - R(3)
!!
      Br(1) =  0.5 * S23
      Br(2) = -0.5 * S13
      Br(3) =  0.5 * S12
      Bs(1) = -0.5 * R23
      Bs(2) =  0.5 * R13
      Bs(3) = -0.5 * R12
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      MEMBT(NEL)%RES%Area = 0.5 * (R12*S13 - R13*S12)
      Ain = 1.0 / MEMBT(NEL)%RES%Area
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,3
        Dx = Dx + Br(i)*Br(i) + Bs(i)*Bs(i)
      ENDDO
      Delta = Ain * SQRT(Dx)
!!
!! VELOCITY GRADIENTS, STRETCHING, AND SPIN.
!! Construct velocity gradients (still need to be divided by area).
!!
      Vrr = 0.0
      Vsr = 0.0
      Vrs = 0.0
      Vss = 0.0
      DO i = 1,3
        Vrr = Vrr + Vr(i)*Br(i)
        Vsr = Vsr + Vs(i)*Br(i)
        Vrs = Vrs + Vr(i)*Bs(i)
        Vss = Vss + Vs(i)*Bs(i)
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
      Wrs = Vrs - Drs
!!
!! POLAR DECOMPOSITION OF THE DEFORMATION GRADIENT.
!! Compute the polar decomposition of the deformation gradient beteen
!! time (n) and time (n+1). Construct components of the deformation
!! gradient; compute the sum of F11 and F22, and the difference of F12
!! and F21.  If r(R,S,t) and s(R,S,t), then QSF = r,R + s,S  and
!! QDF = r,S - s,R. Hin = 1.0 / (2.0*Area(0)) cancels from the quotient
!! in the expression for the inverse tangent function, ATan.
!!
      DTnext = MEMBT(NEL)%RES%DTnext
      dBeta = ATAN (DTnext*(Vrs-Vsr)/(2.0 + DTnext*(Vrr+Vss)))
      MEMBT(NEL)%RES%Beta = MEMBT(NEL)%RES%Beta + dBeta
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBT_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
!! Purpose: Construct triangular membrane gradient operators
!!
      USE shared_common_data
      USE membt_
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
     &          X(3),Y(3),Z(3),Vx(3),Vy(3),Vz(3),                              &
     &          R(3),S(3),Vr(3),Vs(3)
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
      DO i = 1,3
        X(i)  = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i)  = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i)  = EMOTION(i)%Pz + EMOTION(i)%Uz
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! CONSTRUCT LOCAL BASIS VECTORS
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element sides. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
      Rx = X(2)-X(1)
      Ry = Y(2)-Y(1)
      Rz = Z(2)-Z(1)
      Qmag  = 1.0 / SQRT (Rx*Rx+Ry*Ry+Rz*Rz)
      Rx = Rx * Qmag
      Ry = Ry * Qmag
      Rz = Rz * Qmag
      Sx = X(3)-X(1)
      Sy = Y(3)-Y(1)
      Sz = Z(3)-Z(1)
      Qmag  = 1.0 / SQRT (Sx*Sx+Sy*Sy+Sz*Sz)
      Sx = Sx * Qmag
      Sy = Sy * Qmag
      Sz = Sz * Qmag
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') MEMBT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'MEMBT_DIVERGENCE_OPERATOR.001.00'//                     &
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
!! Transform position X,Y,Z, and translational velocity Vx,Vy,Vz, to local
!! coordinates.
!!
      DO i = 1,3
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
      ENDDO
!!
!! MEMBRANE GRADIENT OPERATOR.
!!
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S23 = S(2) - S(3)
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R23 = R(2) - R(3)
!!
      Br(1) =  0.5 * S23
      Br(2) = -0.5 * S13
      Br(3) =  0.5 * S12
      Bs(1) = -0.5 * R23
      Bs(2) =  0.5 * R13
      Bs(3) = -0.5 * R12
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      MEMBT(NEL)%RES%Area = 0.5 * (R12*S13 - R13*S12)
      Ain = 1.0 / MEMBT(NEL)%RES%Area
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,3
        Dx = Dx + Br(i)*Br(i) + Bs(i)*Bs(i)
      ENDDO
      Delta = Ain * SQRT(Dx)
!!
      RETURN
      END
!!_
      SUBROUTINE MEMBT_MASS ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:35
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one third at each node.
!!
      USE shared_common_data
      USE membt_
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
!! Compute one third of the total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = Density*MEMBT(NEL)%PAR%Area*SECTION_2D(SecID)%Thickness/3.
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,3
        NODE(MEMBT(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(MEMBT(NEL)%PAR%IX(i))%Mass + QMass
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
      SUBROUTINE MEMBT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 27-NOV-1991 09:57:36
!!
      USE shared_common_data
      USE membt_
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
     &          Fr(3),Fs(3)
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
        WRITE (MSG1,'(I8)') MEMBT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'MEMBT_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'//                &
     &          MSGL//'M3EL (3-Node Membrane) Element ID:'//MSG1               &
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
      MEMBT(NEL)%RES%DTelt = Dx / (QC + SQRT (QC*QC + Cv*Cv))
!!
!! Compute current element thickness based on constant volume.
!!
      Thickness = SECTION_2D(SecID)%Thickness *                                &
     &          MEMBT(NEL)%PAR%Area / MEMBT(NEL)%RES%Area
!!
!! Divergence of the membrane stress resultants.
!!
      Vrr = Thickness * (MEMBT(NEL)%RES%Stress(1) + QP)
      Vss = Thickness * (MEMBT(NEL)%RES%Stress(2) + QP)
      Vrs = Thickness *  MEMBT(NEL)%RES%Stress(3)
!!
      DO i = 1,3
        Fr(i) = Br(i)*Vrr + Bs(i)*Vrs
        Fs(i) = Br(i)*Vrs + Bs(i)*Vss
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element force array.
!!
      DO i = 1,3
        MEMBT(NEL)%RES%Xint(i) = Rx*Fr(i) + Sx*Fs(i)
        MEMBT(NEL)%RES%Yint(i) = Ry*Fr(i) + Sy*Fs(i)
        MEMBT(NEL)%RES%Zint(i) = Rz*Fr(i) + Sz*Fs(i)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE GLOBAL_TO_LOCAL (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES)
!!
!! Copyright (c) by KEY Associates, 3-APR-1992 14:19:36
!!
!! Purpose: Convert global material direction vectors to local element
!! coordinates. The output direction vectors are unit magnitude.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN)    :: Rx,Ry,Rz
      REAL(KIND(0D0)), INTENT(IN)    :: Sx,Sy,Sz
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(6)
!!
!! Input is in terms of global components.
!!
!!      STATE_VARIABLES(1) = Ax
!!      STATE_VARIABLES(2) = Ay
!!      STATE_VARIABLES(3) = Az
!!      STATE_VARIABLES(4) = Bx
!!      STATE_VARIABLES(5) = By
!!      STATE_VARIABLES(6) = Bz
!!
!! Output is in terms of local element components.
!!
!!      STATE_VARIABLES(1) = Ar
!!      STATE_VARIABLES(2) = As
!!      STATE_VARIABLES(3) = Br
!!      STATE_VARIABLES(4) = Bs
!!      STATE_VARIABLES(5) = 0.0
!!      STATE_VARIABLES(6) = 0.0
!!
      Qr = STATE_VARIABLES(1)*Rx + STATE_VARIABLES(2)*Ry + STATE_VARIABLES(3)*Rz
      Qs = STATE_VARIABLES(1)*Sx + STATE_VARIABLES(2)*Sy + STATE_VARIABLES(3)*Sz
      Qm = 1.0D0 / SQRT (Qr*Qr + Qs*Qs)
      STATE_VARIABLES(1) = Qr * Qm
      STATE_VARIABLES(2) = Qs * Qm
      Qr = STATE_VARIABLES(4)*Rx + STATE_VARIABLES(5)*Ry + STATE_VARIABLES(6)*Rz
      Qs = STATE_VARIABLES(4)*Sx + STATE_VARIABLES(5)*Sy + STATE_VARIABLES(6)*Sz
      Qm = 1.0D0 / SQRT (Qr*Qr + Qs*Qs)
      STATE_VARIABLES(3) = Qr * Qm
      STATE_VARIABLES(4) = Qs * Qm
      STATE_VARIABLES(5) = 0.0
      STATE_VARIABLES(6) = 0.0
!!
      RETURN
      END
