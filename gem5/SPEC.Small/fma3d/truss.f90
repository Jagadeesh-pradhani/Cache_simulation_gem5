      SUBROUTINE TRUSS_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 18-MAY-1991 18:05:04
!!
      USE shared_common_data
      USE truss_
      USE material_
      USE node_
      USE motion_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: SecID,MatID
      LOGICAL :: FOUND
!!
      DO N = 1,NUMTR
!!
!! Gather element motion.
!!
        EMOTION(1) = MOTION(TRUSS(N)%PAR%IX(1))
        EMOTION(2) = MOTION(TRUSS(N)%PAR%IX(2))
!!
!! Initialize element clock and time step
!!
        TRUSS(N)%RES%Time   = 0.0
        TRUSS(N)%RES%DTnext = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element
!! data structure.
!!
        MatID = TRUSS(N)%PAR%MatID
        SecID = TRUSS(N)%PAR%SecID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies
!! only the element length calculation is required.)
!!
        CALL TRUSS_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Save initial element length for later axial strain calculation.
!!
        TRUSS(N)%PAR%Length = TRUSS(N)%RES%Length

        IF (TRUSS(N)%RES%Length .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') TRUSS(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'TRUSS_INITIALIZATION.001.00'//                          &
     &          MSGL//'TRUSS (2-Node Truss) Element ID:'//MSG1//               &
     &          MSGL//'Has A Zero Or Negative Length.'                         &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL TRUSS_MASS ( NEL,SecID,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (TRUSS(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = ABS (EMOTION(1)%Px - EMOTION(2)%Px)                         &
     &           + ABS (EMOTION(1)%Py - EMOTION(2)%Py)                         &
     &           + ABS (EMOTION(1)%Pz - EMOTION(2)%Pz)
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          TRUSS(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Find initial sound speeds and stress.
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (10)
            CALL MATERIAL_10                                                   &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (11)
            CALL MATERIAL_11                                                   &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE (17)
            CALL MATERIAL_17                                                   &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') TRUSS(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'TRUSS_INITIALIZITION.002.00'//                          &
     &          MSGL//'Truss (TRUSS) Element ID:'//MSG1//                      &
     &          MSGL//'References A Material ID:'//MSG2//                      &
     &          MSGL//'With An Invalid Material Type.'                         &
     &          )
          END SELECT
!!
!! Compute initial stress divergence and time step.
!!
          CALL TRUSS_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTTrx = MAX (TIMSIM%DTTrx,TRUSS(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = TRUSS(N)%RES%DTelt .LT. TIMSIM%DTTru(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTTru(j + 1) = TIMSIM%DTTru(j)
              TIMSIM%Truss(j + 1) = TIMSIM%Truss(j)
            ENDDO
          ENDIF
          TIMSIM%DTTru(i) = TRUSS(N)%RES%DTelt
          TIMSIM%Truss(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE TRUSS_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 18-MAY-1991 18:05:04
!!
      USE shared_common_data
      USE truss_
      USE material_
      USE node_
      USE motion_
      USE section_1d_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: SecID,MatID
!!
      DO N = 1,NUMTR
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,TRUSS(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (TRUSS(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            EMOTION(1) = MOTION(TRUSS(N)%PAR%IX(1))
            EMOTION(2) = MOTION(TRUSS(N)%PAR%IX(2))
!!
!! Count element execution.
!!
            COUNTER%TRUSS = COUNTER%TRUSS + 1
!!
!! Increment element clock.
!!
            TRUSS(N)%RES%Time = TRUSS(N)%RES%Time + TRUSS(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            QA = NODE(TRUSS(N)%PAR%IX(1))%Time - TRUSS(N)%RES%Time
            EMOTION(1)%Ux =                                                    &
     &      EMOTION(1)%Ux - QA * EMOTION(1)%Vx
            EMOTION(1)%Uy =                                                    &
     &      EMOTION(1)%Uy - QA * EMOTION(1)%Vy
            EMOTION(1)%Uz =                                                    &
     &      EMOTION(1)%Uz - QA * EMOTION(1)%Vz
!!
            QA = NODE(TRUSS(N)%PAR%IX(2))%Time - TRUSS(N)%RES%Time
            EMOTION(2)%Ux =                                                    &
     &      EMOTION(2)%Ux - QA * EMOTION(2)%Vx
            EMOTION(2)%Uy =                                                    &
     &      EMOTION(2)%Uy - QA * EMOTION(2)%Vy
            EMOTION(2)%Uz =                                                    &
     &      EMOTION(2)%Uz - QA * EMOTION(2)%Vz
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Retrieve state variable pointer and material ID from element
!! data structure.
!!
            MatID = TRUSS(N)%PAR%MatID
            SecID = TRUSS(N)%PAR%SecID
!!
!! Gradient operator, stretching, and rotation.
!!
            CALL TRUSS_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Update stress-strain model to obtain new stress.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (10)
              CALL MATERIAL_10                                                 &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (11)
              CALL MATERIAL_11                                                 &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            CASE (17)
              CALL MATERIAL_17                                                 &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress,                                           &
     &          STATE_VARIABLES(TRUSS(N)%PAR%Isv),                             &
     &          TRUSS(N)%RES%Int_Eng,                                          &
     &          TRUSS(N)%RES%DTnext,                                           &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Compute stress divergence and time step.
!!
            CALL TRUSS_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTTrx = MAX (TIMSIM%DTTrx,TRUSS(N)%RES%DTelt)
          IF (TRUSS(N)%RES%DTelt .LT. TIMSIM%DTTru(1)) THEN
            TIMSIM%DTTru(1) = TRUSS(N)%RES%DTelt
            TIMSIM%Truss(1) = N
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
      SUBROUTINE TRUSS_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 22-JUL-1989 18:53:25
!! Migrated by: S W Key, 18-MAY-1991 18:38:19
!!
!! Purpose: Calculate the gradient of the axial velocity.
!!
      USE shared_common_data
      USE truss_
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
     &          Vx(2),Vy(2),Vz(2),X(2),Y(2),Z(2)
!!
      COMMON /TRUSS/                                                           &
     &          Rx,Ry,Rz,Drr
!!
!! Current position of nodal points, and current translational velocities.
!!
      DO i = 1,2
        X(i)  = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i)  = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i)  = EMOTION(i)%Pz + EMOTION(i)%Uz
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! Define a basis vector R aligned with the truss.
!!
      Rx = X(2) - X(1)
      Ry = Y(2) - Y(1)
      Rz = Z(2) - Z(1)
      Tlength = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (ONE / Tlength)
      Ry = Ry * (ONE / Tlength)
      Rz = Rz * (ONE / Tlength)
      TRUSS(NEL)%RES%Length = Tlength
!!
!! Axial stretching using velocity Vx,Vy,Vz transformed to local R coordinate.
!!
      dVr = Rx*(Vx(2)-Vx(1)) + Ry*(Vy(2)-Vy(1)) + Rz*(Vz(2)-Vz(1))
      Drr = dVr * (ONE / Tlength)
!!
      RETURN
      END
!!_
      SUBROUTINE TRUSS_MASS ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 19-MAY-1991 10:45:07
!!
!! Purpose: For the truss element calculate the translational masses at the
!! nodes. The simplest of mass lumpings is used: one half the total mass at
!! each node.
!!
      USE shared_common_data
      USE truss_
      USE section_1d_
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
!! Compute one half total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = Density * TRUSS(NEL)%PAR%Length *                                &
     &  SECTION_1D(SecID)%Area / 2.0D0
!!
      DO i = 1,2
        NODE(TRUSS(NEL)%PAR%IX(i))%Mass =                                      &
     &    NODE(TRUSS(NEL)%PAR%IX(i))%Mass + QMass
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
      SUBROUTINE TRUSS_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 22-JUL-1989 19:05:24
!! Migrated by: S W Key, 19-MAY-1991 10:51:56
!!
!! Purpose: Calculate the critical time step and the divergence of the
!! stresses.
!!
      USE shared_common_data
      USE truss_
      USE section_1d_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL    ! Current element index
      INTEGER, INTENT(IN) :: SecID  ! Current element section index
      INTEGER, INTENT(IN) :: MatID  ! Current element material index
!!
      COMMON /TRUSS/                                                           &
     &          Rx,Ry,Rz,Drr
!!
!! Retrieve material properties.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      Bulk_Ln = MATERIAL(MatID)%PVAL(2)
      Bulk_Qd = MATERIAL(MatID)%PVAL(3)
!!
!! Find sound speed.
!!
      CSQ = SOUND_SPEED%RCL2 / Density
      IF (CSQ .GT. 0.0) THEN
        Cv = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') TRUSS(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'TRUSS_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'//                &
     &          MSGL//'Element ID:'//MSG1                                      &
     &          )
        Cv = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Artificial bulk viscosity.
!!
      IF (Drr .LT. 0.0) THEN
        QC = Bulk_Qd * Bulk_Qd * TRUSS(NEL)%PAR%Length * ABS(Drr)              &
     &     + Bulk_Ln * Cv
        QP = Density * TRUSS(NEL)%PAR%Length * Drr * QC
      ELSE
        QC = 0.0
        QP = 0.0
      ENDIF
!!
!! Critical time step calculation.
!!
      TRUSS(NEL)%RES%DTelt =                                                   &
     &  TRUSS(NEL)%RES%Length / (QC + SQRT (QC*QC + Cv*Cv))
!!
!! Accumulate element stress divergence results in local element arrays.
!!
      Area = SECTION_1D(SecID)%Area *                                          &
     &       TRUSS(NEL)%PAR%Length / TRUSS(NEL)%RES%Length
      Fr = -Area * (TRUSS(NEL)%RES%Stress + QP)
      Fx = Rx * Fr
      Fy = Ry * Fr
      Fz = Rz * Fr
      TRUSS(NEL)%RES%Xint(1) = +Fx
      TRUSS(NEL)%RES%Yint(1) = +Fy
      TRUSS(NEL)%RES%Zint(1) = +Fz
      TRUSS(NEL)%RES%Xint(2) = -Fx
      TRUSS(NEL)%RES%Yint(2) = -Fy
      TRUSS(NEL)%RES%Zint(2) = -Fz
!!
      RETURN
      END
