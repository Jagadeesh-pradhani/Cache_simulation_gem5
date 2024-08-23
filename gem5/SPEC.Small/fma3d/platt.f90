      SUBROUTINE PLATT_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 15-JUN-1991 16:05:51
!!
      USE shared_common_data
      USE platt_
      USE material_
      USE node_
      USE motion_
      USE section_2d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: SecID,MatID
      LOGICAL :: FOUND
!!
      COMMON /PLATT/              &
     &  At(3),                    & ! mean averaging operators
     &  Br(3),Bs(3),              & ! mean gradient  operator
     &  Arr,Ass,Ars,Art,Ast,      & ! mean membrane  stress
     &  Brr,Bss,Brs,              & ! mean bending   stress
     &  Grr,Gss,Grs,Grt,Gst,Qrs,  & ! mean membrane  stretching
     &  Hrr,Hss,Hrs,        Urs,  & ! mean bending   stretching
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness
!!
      DO N = 1,NUMP3
!!
!! Gather element motion.
!!
        DO i = 1,6
          EMOTION(i) = MOTION(PLATT(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        PLATT(N)%RES%Time   = 0.0
        PLATT(N)%RES%DTnext = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
        MatID = PLATT(N)%PAR%MatID
        SecID = PLATT(N)%PAR%SecID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies
!! only the element volume calculation is required.)
!!
        CALL PLATT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Save initial element area for later area strain calculations.
!!
        PLATT(N)%PAR%Area = PLATT(N)%RES%Area

        IF (PLATT(N)%RES%Area .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') PLATT(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PLATT_INITIALIZATION.001.00'//                          &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Has A Zero Or Negative Area.'                           &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL PLATT_MASS ( NEL,SecID,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (PLATT(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,3
            Qdist = ABS (EMOTION(1)%Px - EMOTION(i)%Px)                        &
     &            + ABS (EMOTION(1)%Py - EMOTION(i)%Py)                        &
     &            + ABS (EMOTION(1)%Pz - EMOTION(i)%Pz)
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          PLATT(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Convert global material direction vectors to local element coordinates.
!!
          IF (MATERIAL(MatID)%Type .EQ. 45) THEN
            Isv = PLATT(N)%PAR%Isv
            Nsv = MATERIAL(MatID)%Nsv
            DO i = 1,Ipts_PLATT (NEL)
              CALL GLOBAL_TO_LOCAL                                             &
     &          (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))
              Isv = Isv + Nsv
            ENDDO
          ENDIF
!!
!! Find initial sound speeds and stress.
!!
          CALL PLATT_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Compute initial stress divergence and time step.
!!
          CALL PLATT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTP3x = MAX (TIMSIM%DTP3x,PLATT(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = PLATT(N)%RES%DTelt .LT. TIMSIM%DTPl3(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTPl3(j + 1) = TIMSIM%DTPl3(j)
              TIMSIM%Plat3(j + 1) = TIMSIM%Plat3(j)
            ENDDO
          ENDIF
          TIMSIM%DTPl3(i) = PLATT(N)%RES%DTelt
          TIMSIM%Plat3(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PLATT_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 15-JUN-1991 16:06:07
!!
      USE shared_common_data
      USE platt_
      USE material_
      USE node_
      USE motion_
      USE section_2d_
      USE stress_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: SecID,MatID
!!
      DO N = 1,NUMP3
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,PLATT(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (PLATT(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            DO i = 1,6
              EMOTION(i) = MOTION(PLATT(N)%PAR%IX(i))
            ENDDO
!!
!! Count element execution.
!!
            COUNTER%PLATT = COUNTER%PLATT + 1
!!
!! Increment element clock.
!!
            PLATT(N)%RES%Time = PLATT(N)%RES%Time + PLATT(N)%RES%DTnext
!!
!! Scale nodal positions to current element time. Note that rotational
!! degrees of freedom do not require any scaling since they are not used.
!!
            DO i = 1,3
              QA = NODE(PLATT(N)%PAR%IX(i))%Time - PLATT(N)%RES%Time
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
            MatID = PLATT(N)%PAR%MatID
            SecID = PLATT(N)%PAR%SecID
!!
!! Gradient operator, stretching, and rotation.
!!
            CALL PLATT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Find initial sound speeds and stress.
!!
            CALL PLATT_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Divergence operator.
!!
            CALL PLATT_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Compute initial stress divergence and time step.
!!
            CALL PLATT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTP3x = MAX (TIMSIM%DTP3x,PLATT(N)%RES%DTelt)
          IF (PLATT(N)%RES%DTelt .LT. TIMSIM%DTPl3(1)) THEN
            TIMSIM%DTPl3(1) = PLATT(N)%RES%DTelt
            TIMSIM%Plat3(1) = N
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
      SUBROUTINE PLATT_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platt_
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
     &          Vx(3),Vy(3),Vz(3),Wx(3),Wy(3),Wz(3),X(3),Y(3),Z(3),            &
     &          Vr(3),Vs(3),Vt(3),Wr(3),Ws(3),R(3),S(3),Theta_r(3),            &
     &          Theta_s(3)
!!
      COMMON /PLATT/              &
     &  At(3),                    & ! mean averaging operators
     &  Br(3),Bs(3),              & ! mean gradient  operator
     &  Arr,Ass,Ars,Art,Ast,      & ! mean membrane  stress
     &  Brr,Bss,Brs,              & ! mean bending   stress
     &  Grr,Gss,Grs,Grt,Gst,Qrs,  & ! mean membrane  stretching
     &  Hrr,Hss,Hrs,        Urs,  & ! mean bending   stretching
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness
!!
!! Mid-interval position of nodal points, and current translational and
!! rotational velocities.
!!
      DT = 0.5D0 * PLATT(NEL)%RES%DTnext
      DO i = 1,3
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
        Wx(i) = EMOTION(i+3)%Vx
        Wy(i) = EMOTION(i+3)%Vy
        Wz(i) = EMOTION(i+3)%Vz
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
      Qmag  = ONE / SQRT (Rx*Rx+Ry*Ry+Rz*Rz)
      Rx = Rx * Qmag
      Ry = Ry * Qmag
      Rz = Rz * Qmag
      Sx = X(3)-X(1)
      Sy = Y(3)-Y(1)
      Sz = Z(3)-Z(1)
      Qmag  = ONE / SQRT (Sx*Sx+Sy*Sy+Sz*Sz)
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
        WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATT_GRADIENT_OPERATOR.002.00'//                       &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Element Geometry Has An Undefined Normal.'              &
     &          )
      ENDIF
!!
      Tx = Tx * (ONE / Tmag)
      Ty = Ty * (ONE / Tmag)
      Tz = Tz * (ONE / Tmag)
!!
!! Redefine S to be orthogonal to T and R.
!!
      Sx = Ty*Rz - Ry*Tz
      Sy = Tz*Rx - Rz*Tx
      Sz = Tx*Ry - Rx*Ty
!!
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and
!! rotational velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,3
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
        Wr(i) = Rx*Wx(i) + Ry*Wy(i) + Rz*Wz(i)
        Ws(i) = Sx*Wx(i) + Sy*Wy(i) + Sz*Wz(i)
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      S23 = S(2) - S(3)
      S31 = S(3) - S(1)
      S12 = S(1) - S(2)
      R32 = R(3) - R(2)
      R13 = R(1) - R(3)
      R21 = R(2) - R(1)
!!
!! At(i); mean averaging operator (area weighting based on shape functions).
!!
      QA = 1.66666666666667D-01 * (R21*S31 - S12*R13)
      At(1) = QA
      At(2) = QA
      At(3) = QA
!!
!! Membrane gradient operator.
!!
      Br(1) = 0.5D0 * S23
      Br(2) = 0.5D0 * S31
      Br(3) = 0.5D0 * S12
      Bs(1) = 0.5D0 * R32
      Bs(2) = 0.5D0 * R13
      Bs(3) = 0.5D0 * R21
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATT(NEL)%RES%Area = 0.5D0 * (R21*S31 - S12*R13)
      Ain = ONE / PLATT(NEL)%RES%Area
      Hin = 0.5D0 * Ain
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,3
        Dx = Dx + Br(i)*Br(i) + Bs(i)*Bs(i)
      ENDDO
      Delta = Ain * SQRT(Dx)
!!
!! Convert global rotation rates to local theta's.
!!
      DO i = 1,3
        Theta_r(i) =  Ws(i)
        Theta_s(i) = -Wr(i)
      ENDDO
!!
!! I. MEMBRANE: mean membrane velocity gradients.
!!
      Vrr = 0.0
      Vsr = 0.0
      Vtr = 0.0
      Vrs = 0.0
      Vss = 0.0
      Vts = 0.0
      Vrt = 0.0
      Vst = 0.0
      DO i = 1,3
        Vrr = Vrr +  Vr(i)*Br(i)
        Vsr = Vsr +  Vs(i)*Br(i)
        Vtr = Vtr +  Vt(i)*Br(i)
        Vrs = Vrs +  Vr(i)*Bs(i)
        Vss = Vss +  Vs(i)*Bs(i)
        Vts = Vts +  Vt(i)*Bs(i)
        Vrt = Vrt +  Theta_r(i)*At(i)
        Vst = Vst +  Theta_s(i)*At(i)
      ENDDO
!!
!! Mean membrane stretching components.
!!
      Grr = Ain * Vrr
      Gss = Ain * Vss
      Grs = Hin * (Vrs + Vsr)
!!
!! Mean in-plane spin component.
!!
      Qrs = Hin * (Vrs - Vsr)
!!
!! Mean transverse-shear stretching components.
!!
      Grt = Hin * (Vrt + Vtr)
      Gst = Hin * (Vst + Vts)
!!
!! II. BENDING: Mean bending velocity gradients.
!!
      Vrr = 0.0
      Vsr = 0.0
      Vrs = 0.0
      Vss = 0.0
      DO i = 1,3
        Vrr = Vrr + Theta_r(i)*Br(i)
        Vsr = Vsr + Theta_s(i)*Br(i)
        Vrs = Vrs + Theta_r(i)*Bs(i)
        Vss = Vss + Theta_s(i)*Bs(i)
      ENDDO
!!
!! Mean bending stretching.
!!
      Hrr = Ain * Vrr
      Hss = Ain * Vss
      Hrs = Hin * (Vrs + Vsr)
      Urs = Hin * (Vrs - Vsr)
!!
      RETURN
      END
!!_
      SUBROUTINE PLATT_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Migrated by: S W Key, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platt_
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
     &          Vx(3),Vy(3),Vz(3),Wx(3),Wy(3),Wz(3),X(3),Y(3),Z(3),            &
     &          Vr(3),Vs(3),Vt(3),Wr(3),Ws(3),R(3),S(3),Theta_r(3),            &
     &          Theta_s(3)
!!
      COMMON /PLATT/              &
     &  At(3),                    & ! mean averaging operators
     &  Br(3),Bs(3),              & ! mean gradient  operator
     &  Arr,Ass,Ars,Art,Ast,      & ! mean membrane  stress
     &  Brr,Bss,Brs,              & ! mean bending   stress
     &  Grr,Gss,Grs,Grt,Gst,Qrs,  & ! mean membrane  stretching
     &  Hrr,Hss,Hrs,        Urs,  & ! mean bending   stretching
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness
!!
!! Current position of nodal points, and current translational and
!! rotational velocities.
!!
      DO i = 1,3
        X(i)  = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i)  = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i)  = EMOTION(i)%Pz + EMOTION(i)%Uz
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
        Wx(i) = EMOTION(i+3)%Vx
        Wy(i) = EMOTION(i+3)%Vy
        Wz(i) = EMOTION(i+3)%Vz
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
      Qmag  = ONE / SQRT (Rx*Rx+Ry*Ry+Rz*Rz)
      Rx = Rx * Qmag
      Ry = Ry * Qmag
      Rz = Rz * Qmag
      Sx = X(3)-X(1)
      Sy = Y(3)-Y(1)
      Sz = Z(3)-Z(1)
      Qmag  = ONE / SQRT (Sx*Sx+Sy*Sy+Sz*Sz)
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
        WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATT_DIVERGENCE_OPERATOR.002.00'//                     &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Element Geometry Has An Undefined Normal.'              &
     &          )
      ENDIF
!!
      Tx = Tx * (ONE / Tmag)
      Ty = Ty * (ONE / Tmag)
      Tz = Tz * (ONE / Tmag)
!!
!! Redefine S to be orthogonal to T and R.
!!
      Sx = Ty*Rz - Ry*Tz
      Sy = Tz*Rx - Rz*Tx
      Sz = Tx*Ry - Rx*Ty
!!
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and
!! rotational velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,3
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
        Wr(i) = Rx*Wx(i) + Ry*Wy(i) + Rz*Wz(i)
        Ws(i) = Sx*Wx(i) + Sy*Wy(i) + Sz*Wz(i)
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      S23 = S(2) - S(3)
      S31 = S(3) - S(1)
      S12 = S(1) - S(2)
      R32 = R(3) - R(2)
      R13 = R(1) - R(3)
      R21 = R(2) - R(1)
!!
!! At(i); mean averaging operator (area weighting based on shape functions).
!!
      QA = 1.66666666666667D-01 * (R21*S31 - S12*R13)
      At(1) = QA
      At(2) = QA
      At(3) = QA
!!
!! Membrane gradient operator.
!!
      Br(1) = 0.5D0 * S23
      Br(2) = 0.5D0 * S31
      Br(3) = 0.5D0 * S12
      Bs(1) = 0.5D0 * R32
      Bs(2) = 0.5D0 * R13
      Bs(3) = 0.5D0 * R21
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATT(NEL)%RES%Area = 0.5D0 * (R21*S31 - S12*R13)
      Ain = ONE / PLATT(NEL)%RES%Area
      Hin = 0.5D0 * Ain
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
      SUBROUTINE PLATT_MASS ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 15-JUN-1991 16:06:10
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one third at each node and an isotropic inertia.
!!
      USE shared_common_data
      USE platt_
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
!! Compute one third total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = Density*PLATT(NEL)%PAR%Area*SECTION_2D(SecID)%Thickness/3.0D0
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,3
        NODE(PLATT(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(PLATT(NEL)%PAR%IX(i))%Mass + QMass
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
!! Compute nodal isotropic inertia
!!
      RMass = QMass * (PLATT(NEL)%PAR%Area +                                   &
     &  SECTION_2D(SecID)%Thickness**2) / 12.0D0

      DO i = 4,6
        NODE(PLATT(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(PLATT(NEL)%PAR%IX(i))%Mass + RMass
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PLATT_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 15-JUN-1991 16:06:11
!!
      USE shared_common_data
      USE platt_
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
     &          Mr(3),Ms(3),Fr(3),Fs(3),Ft(3)
!!
      COMMON /PLATT/              &
     &  At(3),                    & ! mean averaging operators
     &  Br(3),Bs(3),              & ! mean gradient  operator
     &  Arr,Ass,Ars,Art,Ast,      & ! mean membrane  stress
     &  Brr,Bss,Brs,              & ! mean bending   stress
     &  Grr,Gss,Grs,Grt,Gst,Qrs,  & ! mean membrane  stretching
     &  Hrr,Hss,Hrs,        Urs,  & ! mean bending   stretching
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness
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
        WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PLATT_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'                  &
     &          )
        Cv = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Calculate generalized element dimension.
!!
      Dx = ONE / Delta
!!
!! Artificial bulk viscosity pressure.
!!
      Gkk = Grr + Gss
      IF (Gkk .LT. 0.0) THEN
        QC = Bulk_Qd * Bulk_Qd * Dx * ABS(Gkk) + Bulk_Ln * Cv
        QP = Density * Dx * Gkk * QC
      ELSE
        QC = 0.0
        QP = 0.0
      ENDIF
!!
!! Critical time step calculation. Three seperate conditions are considered:
!!      (1) a membrane critical time step DTmb
!!      (2) a transverse shear critical time step DTts
!!      (3) a bending critical time step DTbd
!!
      DTmb = Dx / (QC + SQRT (QC*QC + Cv*Cv))
      DTts = Thickness / (Factor*Cv)
      DTbd = (Dx*Dx) / ((Thickness+Thickness)*Cv)
      PLATT(NEL)%RES%DTelt = MIN (DTmb,DTts,DTbd)
!!
!! DIVERGENCE OF THE MEMBRANE STRESS RESULTANTS
!!
      Arr = Arr + (Thickness * QP)
      Ass = Ass + (Thickness * QP)
      DO i = 1,3
        Fr(i) = Arr*Br(i) + Ars*Bs(i)
        Fs(i) = Ars*Br(i) + Ass*Bs(i)
        Ft(i) = Art*Br(i) + Ast*Bs(i)
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
      DO i = 1,3
        PLATT(NEL)%RES%Xint(i) = Rx*Fr(i) + Sx*Fs(i) + Tx*Ft(i)
        PLATT(NEL)%RES%Yint(i) = Ry*Fr(i) + Sy*Fs(i) + Ty*Ft(i)
        PLATT(NEL)%RES%Zint(i) = Rz*Fr(i) + Sz*Fs(i) + Tz*Ft(i)
      ENDDO
!!
!! DIVERGENCE OF THE BENDING STRESS RESULTANTS
!! Includes area weighted terms in divergence calculation.
!!
      DO i = 1,3
        Fr(i) = Brr*Br(i) + Brs*Bs(i) + Art*At(i)
        Fs(i) = Brs*Br(i) + Bss*Bs(i) + Ast*At(i)
      ENDDO
      DO i = 1,3
        Mr(i) = -Fs(i)
        Ms(i) =  Fr(i)
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
      DO i = 1,3
        PLATT(NEL)%RES%Xint(i+3) = Rx*Mr(i) + Sx*Ms(i)
        PLATT(NEL)%RES%Yint(i+3) = Ry*Mr(i) + Sy*Ms(i)
        PLATT(NEL)%RES%Zint(i+3) = Rz*Mr(i) + Sz*Ms(i)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PLATT_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 15-JUN-1991 16:06:12
!!
      USE shared_common_data
      USE platt_
      USE material_
      USE section_2d_
      USE stress_
      USE state_variables_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      COMMON /PLATT/              &
     &  At(3),                    & ! mean averaging operators
     &  Br(3),Bs(3),              & ! mean gradient  operator
     &  Arr,Ass,Ars,Art,Ast,      & ! mean membrane  stress
     &  Brr,Bss,Brs,              & ! mean bending   stress
     &  Grr,Gss,Grs,Grt,Gst,Qrs,  & ! mean membrane  stretching
     &  Hrr,Hss,Hrs,        Urs,  & ! mean bending   stretching
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness
!!
      INTEGER, PARAMETER :: Ibgn(5) = (/1,3,6,10,15/)
      INTEGER, PARAMETER :: Iend(5) = (/2,5,9,14,20/)
      REAL(KIND(0D0))    :: Int_Eng,Beta(20),Awght(20)

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Non-dimensional thickness coordinates (trapezodial rule).
!!
      REAL(KIND(0D0)), PARAMETER :: ONE6TH = (1.0D0/6.0D0)
      REAL(KIND(0D0)), PARAMETER :: NONE6TH = (-ONE6TH)

      DATA (Beta(i),i= 1, 2) /-.50D+0, .500D+0/
      DATA (Beta(i),i= 3, 5) /-.50D+0, .000D+0, .500D+0/
      DATA (Beta(i),i= 6, 9) /-.50D+0, NONE6TH,  ONE6TH, .50D+0/
      DATA (Beta(i),i=10,14) /-.50D+0,-.250D+0, .000D+0, .25D+0, .50D+0/
      DATA (Beta(i),i=15,20) /-.50D+0,-.300D+0,-.100D+0, .10D+0, .30D+0, .50D+0/
!!
!! Membrane thickness integration weights.
!!
      REAL(KIND(0D0)), PARAMETER :: ONE3RD = (1.0D0/3.0D0)

      DATA (Awght(i),i= 1, 2) /.500D+0, .500D+0/
      DATA (Awght(i),i= 3, 5) /.250D+0, .500D+0, .250D+0/
      DATA (Awght(i),i= 6, 9) / ONE6TH,  ONE3RD,  ONE3RD,  ONE6TH/
      DATA (Awght(i),i=10,14) /.125D+0, .250D+0, .250D+0, .250D+0, .125D+0/
      DATA (Awght(i),i=15,20) /.100D+0, .200D+0, .200D+0, .200D+0, .200D+0, .100D+0/
!!
!! Define constants.
!!
      IF (FIRST) THEN
        SQRT6o1 = SQRT (6.0D+0)
        SQRT5o6 = SQRT (5.0D+0 / 6.0D+0)
        FIRST = .FALSE.
      ENDIF
!!
!! Compute current element thickness based on constant volume.
!!
      Thickness = SECTION_2D(SecID)%Thickness *                                &
     &  PLATT(NEL)%PAR%Area / PLATT(NEL)%RES%Area
!!
!! Transverse shear and thin shell limit correction factor.
!!
!!        1. Based on the work of Cowper the transverse shear strain
!!        energy should be scaled by a factor of 5/6'ths to obtain an
!!        optimal match with continuum theory.
!!
!!        G. R. Cowper, "On the Accuracy of Timoshenko's Beam Theory,"
!!        Journal of the Engineering Mechanics Division, Proceedings of
!!        the American Society of Civil Engineers, Vol. 94, No. EM6,
!!        December 1968.
!!
!!        2. Retention of the transverse shear stiffness drives the
!!        critical time step to extremely low values, as the thickness
!!        diminishes and thin shell behavior is approached. The correc-
!!        tion factor has a minimum psuedo thickness to limit the trans-
!!        verse shear strain energy, namely, 6.0*Thickness**2/Area. The
!!        factor allows sufficient transverse shear stiffness to obtain
!!        thin shell results without severe reduction in critical time
!!        step or locking in the element.
!!
!!        I. Fried, A. Johnson, and A. Tessler, "Minimal-Degree Thin
!!        Triangular Plate and Shell Bending Finite Elements of Order
!!        Two and Four," Computer Methods in Applied Mechanics and
!!        Engineering, Vol. 56, No. 3, July 1986.
!!
!! To maintain a self-adjoint statement of the problem, the transverse
!! shear portion of the gradient/divergence operator is scaled by a
!! combined factor which is the square root of the controling correction.
!! Here, SQRT5o6 = SQRT(5.0/6.0) and SQRT6o1 = SQRT(6.0).
!!
      Aspect_Ratio = Thickness / SQRT(PLATT(NEL)%RES%Area)
      Factor = MIN (SQRT5o6, SQRT6o1*Aspect_Ratio)
!!
!! Apply shear correction factor to transverse shear strains.
!!
      Grt = Factor * Grt
      Gst = Factor * Gst
!!
!! Initialize the mean stresses Arr, Ass, Ars, Art, and Ast.
!!
      Arr = 0.0
      Ass = 0.0
      Ars = 0.0
      Art = 0.0
      Ast = 0.0
!!
!! Initialize the moment stresses Brr, Bss, and Brs.
!!
      Brr = 0.0
      Bss = 0.0
      Brs = 0.0
!!
!! Initialize stress state pointer Ist, state variable pointer Isv and
!! state variable interval Nsv.
!!
      Ist = PLATT(NEL)%PAR%Ist
      Isv = PLATT(NEL)%PAR%Isv
      Nsv = MATERIAL(MatID)%Nsv
!!
!! Select between internal integration rule and user-supplied integration
!! rule. Calculate Xwght for use in constructing moment integration end-
!! point weighting Ewght.
!!
      IF (SECTION_2D(SecID)%Ipts .NE. 0) THEN
        Ione = Ibgn(SECTION_2D(SecID)%Ipts-1)
        Itwo = Iend(SECTION_2D(SecID)%Ipts-1)
        Xwght = 0.1666666666666667D0 * Thickness /                             &
     &    DBLE ((SECTION_2D(SecID)%Ipts-1)**2)
!!
!! Integrate through the thickness of the plate to obtain mean stress
!! and first moment of the stress.
!!
        DO i = Ione,Itwo
!!
!! Calculate the in-plane stretching and spin.
!!
          Eta = Thickness * (Beta(i) - SECTION_2D(SecID)%RefLoc)
          Vrr = Grr + Eta * Hrr
          Vss = Gss + Eta * Hss
          Vrs = Grs + Eta * Hrs
          Vrt = Grt
          Vst = Gst
          Xrs = Qrs + Eta * Urs
!!
!! Internal energy from time n to time n+1/2.
!!
          Int_Eng =  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                   &
     &            + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                  &
     &            + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                  &
     &            + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Calculate stresses. (Incremental constitutive evaluation.)
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (40)
            CALL MATERIAL_40                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (41)
            CALL MATERIAL_41                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID,RCL2,RCS2                                                &
     &          )
          CASE (45)
            CALL MATERIAL_45                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (47)
            CALL MATERIAL_47                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATT_STRESS_INTEGRATION.001.00'//                      &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Via a MATERIAL Record References An '                   &
     &              //'Invalid Material Type:'//MSG2                           &
     &          )
          END SELECT
!!
!! Internal energy from time n+1/2 to time n+1.
!!
          Int_Eng = Int_Eng                                                    &
     &            +  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                   &
     &            + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                  &
     &            + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                  &
     &            + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Accumulate layer internal energy density.
!!
          Int_Eng = (0.5D0*PLATT(NEL)%RES%DTnext) * AWght(i) * Int_Eng
!!
          PLATT(NEL)%RES%Int_Eng = PLATT(NEL)%RES%Int_Eng + Int_Eng
!!
!! Mean stress = Membrane stress resultant / Thickness.
!!
          Arr = Arr + Awght(i) * STRESS(1,Ist)
          Ass = Ass + Awght(i) * STRESS(2,Ist)
          Ars = Ars + Awght(i) * STRESS(4,Ist)
          Art = Art + Awght(i) * STRESS(5,Ist)
          Ast = Ast + Awght(i) * STRESS(6,Ist)
!!
!! Moment stress = Bending stress resultant / Thickness.
!!
          Ewght = 0.0
          IF (i .eq. Ione) Ewght = +Xwght
          IF (i .eq. Itwo) Ewght = -Xwght
!!
          Brr = Brr + (Eta * Awght(i) + Ewght) * STRESS(1,Ist)
          Bss = Bss + (Eta * Awght(i) + Ewght) * STRESS(2,Ist)
          Brs = Brs + (Eta * Awght(i) + Ewght) * STRESS(4,Ist)
!!
          Ist = Ist + 1
          Isv = Isv + Nsv
        ENDDO
!!
      ELSE IF (SECTION_2D(SecID)%Irule .NE. 0) THEN
        Ione=1
        Itwo=TABULATED_FUNCTION(SECTION_2D(SecID)%Irule)%Number_of_Pairs
!!
!! Integrate through the thickness of the plate to obtain mean stress
!! and first moment of the stress.
!!
        DO i = Ione,Itwo
          Zeta   = TABULATED_FUNCTION(SECTION_2D(SecID)%Irule)%X(i)
          Weight = TABULATED_FUNCTION(SECTION_2D(SecID)%Irule)%Y(i)
!!
!! Calculate the in-plane stretching and spin.
!!
          Eta = Thickness * (Zeta - SECTION_2D(SecID)%RefLoc)
          Vrr = Grr + Eta * Hrr
          Vss = Gss + Eta * Hss
          Vrs = Grs + Eta * Hrs
          Vrt = Grt
          Vst = Gst
          Xrs = Qrs + Eta * Urs
!!
!! Internal energy from time n to time n+1/2.
!!
          Int_Eng =  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                   &
     &            + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                  &
     &            + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                  &
     &            + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Calculate stresses. (Incremental constitutive evaluation.)
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (40)
            CALL MATERIAL_40                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (41)
            CALL MATERIAL_41                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID,RCL2,RCS2                                                &
     &          )
          CASE (45)
            CALL MATERIAL_45                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (47)
            CALL MATERIAL_47                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATT(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATT_STRESS_INTEGRATION.001.00'//                      &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Via a MATERIAL Record References An '                   &
     &              //'Invalid Material Type:'//MSG2                           &
     &          )
          END SELECT
!!
!! Internal energy from time n+1/2 to time n+1.
!!
          Int_Eng = Int_Eng                                                    &
     &            +  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                   &
     &            + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                  &
     &            + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                  &
     &            + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Accumulate layer internal energy density.
!!
          Int_Eng = (0.5D0*PLATT(NEL)%RES%DTnext) * Weight * Int_Eng
!!
          PLATT(NEL)%RES%Int_Eng = PLATT(NEL)%RES%Int_Eng + Int_Eng
!!
!! Mean stress = Membrane stress resultant / Thickness.
!!
          Arr = Arr + Weight * STRESS(1,Ist)
          Ass = Ass + Weight * STRESS(2,Ist)
          Ars = Ars + Weight * STRESS(4,Ist)
          Art = Art + Weight * STRESS(5,Ist)
          Ast = Ast + Weight * STRESS(6,Ist)
!!
!! Moment stress = Bending stress resultant / Thickness.
!!
          Brr = Brr + (Eta * Weight) * STRESS(1,Ist)
          Bss = Bss + (Eta * Weight) * STRESS(2,Ist)
          Brs = Brs + (Eta * Weight) * STRESS(4,Ist)
!!
          Ist = Ist + 1
          Isv = Isv + Nsv
        ENDDO
!!
      ELSE
        WRITE (MSG1,'(I8)') PLATT(NEL)%PAR%EleID
        WRITE (MSG2,'(I8)') SECTION_2D(SecID)%SecID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATT_STRESS_INTEGRATION.002.00'//                      &
     &          MSGL//'P3EL (3-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'References PSECTION ID:'//MSG2//                        &
     &          MSGL//'With Both Ipts And Irule Equal to Zero.'                &
     &          )
      ENDIF
!!
!! Convert mean stresses to stress resultants. Apply shear correction factor
!! to transverse-shear stresses.
!!
      Arr = Thickness * Arr
      Ass = Thickness * Ass
      Ars = Thickness * Ars
      Art = (Factor*Thickness) * Art
      Ast = (Factor*Thickness) * Ast
      Brr = Thickness * Brr
      Bss = Thickness * Bss
      Brs = Thickness * Brs
!!
      RETURN
      END
