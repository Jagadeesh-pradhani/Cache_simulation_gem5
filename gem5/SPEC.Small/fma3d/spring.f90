      SUBROUTINE SPRING_INITIALIZATION
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
      USE shared_common_data
      USE spring_
      USE material_
      USE node_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: MatID
      LOGICAL :: FOUND
!!
      DO N = 1,NUMSP
!!
!! Gather element motion.
!!
        EMOTION(1) = MOTION(SPRING(N)%PAR%IX(1))
        EMOTION(2) = MOTION(SPRING(N)%PAR%IX(2))
!!
!! Initialize element clock and time step
!!
        SPRING(N)%RES%Time   = 0.0
        SPRING(N)%RES%DTnext = 0.0
!!
!! Initialize spring length change.
!!
        SPRING(N)%RES%Delta = 0.0
!!
!! For spring directions defined by element orientation, define
!! zero-length direction to be the x-axis.
!!
        IF (SPRING(N)%PAR%Idir .EQ. 0) THEN
          Rx = ONE
          Ry = 0.0
          Rz = 0.0
!!
!! For fixed direction springs define direction vector.
!!
        ELSE IF (SPRING(N)%PAR%Idir .EQ. 1) THEN
          Rx = ONE
          Ry = 0.0
          Rz = 0.0
        ELSE IF (SPRING(N)%PAR%Idir .EQ. 2) THEN
          Rx = 0.0
          Ry = ONE
          Rz = 0.0
        ELSE IF (SPRING(N)%PAR%Idir .EQ. 3) THEN
          Rx = 0.0
          Ry = 0.0
          Rz = ONE
        ELSE IF (SPRING(N)%PAR%Idir .EQ. 4) THEN
          Rx = SPRING(N)%PAR%Axis(1)
          Ry = SPRING(N)%PAR%Axis(2)
          Rz = SPRING(N)%PAR%Axis(3)
          SLength = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
          IF (SLength .GT. 1.0D-20) THEN
            Rx = Rx * (ONE / SLength)
            Ry = Ry * (ONE / SLength)
            Rz = Rz * (ONE / SLength)
          ELSE
            WRITE (MSG1,'(I8)') SPRING(N)%PAR%EleID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SPRING_INITIALIZATION.001.01'//                         &
     &          MSGL//'Spring (SPRING) Element ID:'//MSG1//                    &
     &          MSGL//'Spring Direction Idir=4 Has A'//                        &
     &          MSGL//'Zero-Length Direction Vector.'                          &
     &          )
          ENDIF
        ENDIF
        SPRING(N)%PAR%Axis(1) = Rx
        SPRING(N)%PAR%Axis(2) = Ry
        SPRING(N)%PAR%Axis(3) = Rz
        SPRING(N)%RES%Axis(1) = Rx
        SPRING(N)%RES%Axis(2) = Ry
        SPRING(N)%RES%Axis(3) = Rz
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve material ID from element data structure.
!!
        MatID = SPRING(N)%PAR%MatID
!!
!! Current spring length.
!!
        CALL SPRING_DISPLACEMENT ( NEL )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (SPRING(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = ABS (EMOTION(1)%Px - EMOTION(2)%Px)                         &
     &             + ABS (EMOTION(1)%Py - EMOTION(2)%Py)                       &
     &             + ABS (EMOTION(1)%Pz - EMOTION(2)%Pz)
          Density = MATERIAL(MatID)%PVAL(1)
          IF (Density .EQ. 0.0) THEN
            SPRING(N)%RES%DTelt = HUGE ( SPRING(N)%RES%DTelt )
          ELSE
            Ymod = MATERIAL(MatID)%PVAL(6)
            SPRING(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
          ENDIF
!!
        ELSE
!!
!! Find initial sound speeds and force.
!!
          MTRL_TYPE = MATERIAL(MatID)%Type
          SELECT CASE (MTRL_TYPE)
          CASE (1)
            CALL MATERIAL_S1                                                   &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (2)
            CALL MATERIAL_S2                                                   &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (3)
            CALL MATERIAL_S3                                                   &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (5)
            CALL MATERIAL_S5                                                   &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') SPRING(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SPRING_INITIALIZITION.001.02'//                         &
     &          MSGL//'Spring (SPRING) Element ID:'//MSG1//                    &
     &          MSGL//'References Material ID:'//MSG2//                        &
     &          MSGL//'With An Invalid Material Type.'                         &
     &          )
          END SELECT
!!
!! Compute initial force divergence and time step.
!!
          CALL SPRING_FORCE_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTSpx = MAX (TIMSIM%DTSpx,SPRING(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = SPRING(N)%RES%DTelt .LT. TIMSIM%DTSpr(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTSpr(j + 1) = TIMSIM%DTSpr(j)
              TIMSIM%Spring(j + 1) = TIMSIM%Spring(j)
            ENDDO
          ENDIF
          TIMSIM%DTSpr(i) = SPRING(N)%RES%DTelt
          TIMSIM%Spring(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE SPRING_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
      USE shared_common_data
      USE spring_
      USE material_
      USE node_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          MatID
!!
      DO N = 1,NUMSP
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,SPRING(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (SPRING(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            EMOTION(1) = MOTION(SPRING(N)%PAR%IX(1))
            EMOTION(2) = MOTION(SPRING(N)%PAR%IX(2))
!!
!! Count element execution.
!!
            COUNTER%SPRING = COUNTER%SPRING + 1
!!
!! Increment element clock.
!!
            SPRING(N)%RES%Time =                                               &
     &        SPRING(N)%RES%Time + SPRING(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            QA = NODE(SPRING(N)%PAR%IX(1))%Time - SPRING(N)%RES%Time
            EMOTION(1)%Ux =                                                    &
     &      EMOTION(1)%Ux - QA * EMOTION(1)%Vx
            EMOTION(1)%Uy =                                                    &
     &      EMOTION(1)%Uy - QA * EMOTION(1)%Vy
            EMOTION(1)%Uz =                                                    &
     &      EMOTION(1)%Uz - QA * EMOTION(1)%Vz
!!
            QA = NODE(SPRING(N)%PAR%IX(2))%Time - SPRING(N)%RES%Time
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
            MatID = SPRING(N)%PAR%MatID
!!
!! Computer current spring displacement.
!!
            CALL SPRING_DISPLACEMENT ( NEL )
!!
!! Update force-displacement model to obtain new force.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (1)
              CALL MATERIAL_S1                                                 &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (2)
              CALL MATERIAL_S2                                                 &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (3)
              CALL MATERIAL_S3                                                 &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (5)
              CALL MATERIAL_S5                                                 &
     &          (                                                              &
     &          SPRING(N)%RES%Force,                                           &
     &          STATE_VARIABLES(SPRING(N)%PAR%Isv),                            &
     &          SPRING(N)%RES%Int_Eng,                                         &
     &          SPRING(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Compute force divergence and time step.
!!
            CALL SPRING_FORCE_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTSpx = MAX (TIMSIM%DTSpx,SPRING(N)%RES%DTelt)
          IF (SPRING(N)%RES%DTelt .LT. TIMSIM%DTSpr(1)) THEN
            TIMSIM%DTSpr(1) = SPRING(N)%RES%DTelt
            TIMSIM%Spring(1) = N
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
      SUBROUTINE SPRING_DISPLACEMENT ( NEL )
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
!! Purpose: Calculate the gradient of the axial velocity.
!!
      USE shared_common_data
      USE spring_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Vx(2),Vy(2),Vz(2),SLength
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Current translational velocities.
!!
      Vx(1) = EMOTION(1)%Vx
      Vy(1) = EMOTION(1)%Vy
      Vz(1) = EMOTION(1)%Vz
      Vx(2) = EMOTION(2)%Vx
      Vy(2) = EMOTION(2)%Vy
      Vz(2) = EMOTION(2)%Vz
!!
!! Define a basis vector R aligned with the spring.
!!
      IF (SPRING(NEL)%PAR%Idir .EQ. 0) THEN
!!
!! Current position of nodal points.
!!
        X1 = EMOTION(1)%Px + EMOTION(1)%Ux
        Y1 = EMOTION(1)%Py + EMOTION(1)%Uy
        Z1 = EMOTION(1)%Pz + EMOTION(1)%Uz
        X2 = EMOTION(2)%Px + EMOTION(2)%Ux
        Y2 = EMOTION(2)%Py + EMOTION(2)%Uy
        Z2 = EMOTION(2)%Pz + EMOTION(2)%Uz

        Rx = X2 - X1
        Ry = Y2 - Y1
        Rz = Z2 - Z1
        SLength = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! If spring length is vanishingly small, use last recorded spring direction.
!!
        IF (SLength .LT. 1.0D-25) THEN
          Rx = SPRING(NEL)%RES%Axis(1)
          Ry = SPRING(NEL)%RES%Axis(2)
          Rz = SPRING(NEL)%RES%Axis(3)
          SLength = ONE
        ENDIF
        Rx = Rx * (ONE / SLength)
        Ry = Ry * (ONE / SLength)
        Rz = Rz * (ONE / SLength)
!!
!! Save current spring orientation (unit vector aligned with spring) for use
!! in the event that the spring length goes to zero in the future.
!!
        SPRING(NEL)%RES%Axis(1) = Rx
        SPRING(NEL)%RES%Axis(2) = Ry
        SPRING(NEL)%RES%Axis(3) = Rz
!!
!! This is a spring with a fixed orientation; retrieve orientation.
!!
      ELSE
        Rx = SPRING(NEL)%RES%Axis(1)
        Ry = SPRING(NEL)%RES%Axis(2)
        Rz = SPRING(NEL)%RES%Axis(3)
      ENDIF
!!
!! Relative axial velocity using Vx,Vy,Vz transformed to local R coordinate.
!!
      Drr = Rx*(Vx(2)-Vx(1)) + Ry*(Vy(2)-Vy(1)) + Rz*(Vz(2)-Vz(1))
!!
!! Update spring displacement.
!!
      SPRING(NEL)%RES%Delta =                                                  &
     &  SPRING(NEL)%RES%Delta + SPRING(NEL)%RES%DTnext * Drr
!!
      RETURN
      END
!!_
      SUBROUTINE SPRING_FORCE_DIVERGENCE ( NEL )
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
!! Purpose: Calculate the critical time step and the divergence of the
!! forces.
!!
      USE shared_common_data
      USE spring_
      USE material_
      USE node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Critical time step calculation. (This assumes that half the mass at each
!! end is "associated" with the spring which is similar to what occurs with
!! a truss element.) The two formulae used are (1) dTcrit < 2/Omega and (2)
!! Omega**2 = k/m. EFM is an "effective mass" that takes into account the
!! difference in size of the two masses at either end of the spring.
!!
      Spring_Constant = SOUND_SPEED%RCL2
      IF (Spring_Constant .LT. ONE) Spring_Constant = ONE
      SM1 = NODE(SPRING(NEL)%PAR%IX(1))%Mass
      SM2 = NODE(SPRING(NEL)%PAR%IX(2))%Mass
      EFM = (SM1 * SM2)/(SM1 + SM2)
      SPRING(NEL)%RES%DTelt = SQRT ((EFM + EFM) / Spring_Constant)
!!
!! Accumulate element force divergence results in local arrays.
!!
      Fr = SPRING(NEL)%RES%Force
      Fx = Rx * Fr
      Fy = Ry * Fr
      Fz = Rz * Fr
      SPRING(NEL)%RES%Xint(1) = -Fx
      SPRING(NEL)%RES%Yint(1) = -Fy
      SPRING(NEL)%RES%Zint(1) = -Fz
      SPRING(NEL)%RES%Xint(2) = +Fx
      SPRING(NEL)%RES%Yint(2) = +Fy
      SPRING(NEL)%RES%Zint(2) = +Fz
!!
      RETURN
      END
