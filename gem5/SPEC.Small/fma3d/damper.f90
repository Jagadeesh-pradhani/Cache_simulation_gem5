      SUBROUTINE DAMPER_INITIALIZATION
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
      USE shared_common_data
      USE damper_
      USE material_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: NEL,MatID
!!!      LOGICAL :: FOUND
!!
      DO N = 1,NUMDM
!!
!! Gather element motion.
!!
        EMOTION(1) = MOTION(DAMPER(N)%PAR%IX(1))
        EMOTION(2) = MOTION(DAMPER(N)%PAR%IX(2))
!!
!! Initialize element clock and time step
!!
        DAMPER(N)%RES%Time   = 0.0
        DAMPER(N)%RES%DTnext = 0.0
!!
!! For damper directions defined by element orientation, define
!! zero-length direction to be the x-axis.
!!
        IF (DAMPER(N)%PAR%Idir .EQ. 0) THEN
          Rx = ONE
          Ry = 0.0
          Rz = 0.0
!!
!! For fixed direction dampers, define direction vector.
!!
        ELSE IF (DAMPER(N)%PAR%Idir .EQ. 1) THEN
          Rx = ONE
          Ry = 0.0
          Rz = 0.0
        ELSE IF (DAMPER(N)%PAR%Idir .EQ. 2) THEN
          Rx = 0.0
          Ry = ONE
          Rz = 0.0
        ELSE IF (DAMPER(N)%PAR%Idir .EQ. 3) THEN
          Rx = 0.0
          Ry = 0.0
          Rz = ONE
        ELSE IF (DAMPER(N)%PAR%Idir .EQ. 4) THEN
          Rx = DAMPER(N)%PAR%Axis(1)
          Ry = DAMPER(N)%PAR%Axis(2)
          Rz = DAMPER(N)%PAR%Axis(3)
          SLength = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
          IF (SLength .GT. 1.0D-20) THEN
            Rx = Rx * (ONE / SLength)
            Ry = Ry * (ONE / SLength)
            Rz = Rz * (ONE / SLength)
          ELSE
            WRITE (MSG1,'(I8)') DAMPER(N)%PAR%EleID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'DAMPER_INITIALIZATION.001.01'//                         &
     &          MSGL//'Damper (DAMPER) Element ID:'//MSG1//                    &
     &          MSGL//'Damper Direction Idir=4 Has A'//                        &
     &          MSGL//'Zero-Length Direction Vector.'                          &
     &          )
          ENDIF
        ENDIF
        DAMPER(N)%PAR%Axis(1) = Rx
        DAMPER(N)%PAR%Axis(2) = Ry
        DAMPER(N)%PAR%Axis(3) = Rz
        DAMPER(N)%RES%Axis(1) = Rx
        DAMPER(N)%RES%Axis(2) = Ry
        DAMPER(N)%RES%Axis(3) = Rz
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve material ID from element data structure.
!!
        MatID = DAMPER(N)%PAR%MatID
!!
!! Initial damper velocity.
!!
        CALL DAMPER_VELOCITY ( NEL )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (DAMPER(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
!!!         ELsize = ABS (EMOTION(1)%Px - EMOTION(2)%Px)
!!!     2          + ABS (EMOTION(1)%Py - EMOTION(2)%Py)
!!!     2          + ABS (EMOTION(1)%Pz - EMOTION(2)%Pz)
!!!         Density = MATERIAL(MatID)%PVAL(1)
!!!         IF (Density .EQ. 0.0) THEN
!!!           DAMPER(N)%RES%DTelt = 1.0D+37
!!!         ELSE
!!!           Ymod = MATERIAL(MatID)%PVAL(6)
!!!           DAMPER(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!!         ENDIF
!!
        ELSE
!!
!! Find initial sound speeds and force.
!!
          Mtype = MATERIAL(MatID)%Type
          SELECT CASE (Mtype)
          CASE (6)
            CALL MATERIAL_D6                                                   &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (7)
            CALL MATERIAL_D7                                                   &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (8)
            CALL MATERIAL_D8                                                   &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE (9)
            CALL MATERIAL_D9                                                   &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') DAMPER(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'DAMPER_INITIALIZITION.001.00'//                         &
     &          MSGL//'Damper (DAMPER) Element ID:'//MSG1//                    &
     &          MSGL//'References Material ID:'//MSG2//                        &
     &          MSGL//'With An Invalid Material Type.'                         &
     &          )
          END SELECT
!!
!! Compute initial force divergence and time step.
!!
          CALL DAMPER_FORCE_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
!!!       TIMSIM%DTDmx = MAX (TIMSIM%DTDmx,DAMPER(N)%RES%DTelt)
!!!       i = 0
!!!       FOUND = .FALSE.
!!!       DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
!!!         i = i + 1
!!!         FOUND = DAMPER(N)%RES%DTelt .LT. TIMSIM%DTDmp(i)
!!!       ENDDO
!!!       IF (FOUND) THEN
!!!         IF (i .LT. NPNDT) THEN
!!!           DO j = NPNDT-1,i,-1
!!!             TIMSIM%DTDmp(j + 1) = TIMSIM%DTDmp(j)
!!!             TIMSIM%DAMPER(j + 1) = TIMSIM%DAMPER(j)
!!!           ENDDO
!!!         ENDIF
!!!         TIMSIM%DTDmp(i) = DAMPER(N)%RES%DTelt
!!!         TIMSIM%DAMPER(i) = N
!!!       ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE DAMPER_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
      USE shared_common_data
      USE damper_
      USE material_
      USE node_
      USE motion_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          NEL,MatID
!!
      DO N = 1,NUMDM
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,DAMPER(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (DAMPER(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            EMOTION(1) = MOTION(DAMPER(N)%PAR%IX(1))
            EMOTION(2) = MOTION(DAMPER(N)%PAR%IX(2))
!!
!! Count element execution.
!!
            COUNTER%DAMPER = COUNTER%DAMPER + 1
!!
!! Increment element clock.
!!
            DAMPER(N)%RES%Time =                                               &
     &        DAMPER(N)%RES%Time + DAMPER(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            QA = NODE(DAMPER(N)%PAR%IX(1))%Time - DAMPER(N)%RES%Time
            EMOTION(1)%Ux =                                                    &
     &        EMOTION(1)%Ux - QA * EMOTION(1)%Vx
            EMOTION(1)%Uy =                                                    &
     &        EMOTION(1)%Uy - QA * EMOTION(1)%Vy
            EMOTION(1)%Uz =                                                    &
     &        EMOTION(1)%Uz - QA * EMOTION(1)%Vz
!!
            QA = NODE(DAMPER(N)%PAR%IX(2))%Time - DAMPER(N)%RES%Time
            EMOTION(2)%Ux =                                                    &
     &        EMOTION(2)%Ux - QA * EMOTION(2)%Vx
            EMOTION(2)%Uy =                                                    &
     &        EMOTION(2)%Uy - QA * EMOTION(2)%Vy
            EMOTION(2)%Uz =                                                    &
     &        EMOTION(2)%Uz - QA * EMOTION(2)%Vz
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Retrieve state variable pointer and material ID from element
!! data structure.
!!
            MatID = DAMPER(N)%PAR%MatID
!!
!! Current damper velocity.
!!
            CALL DAMPER_VELOCITY ( NEL )
!!
!! Update force-velocity model to obtain new force.
!!
            MTRL_TYPE = MATERIAL(MatID)%Type
            SELECT CASE (MTRL_TYPE)
            CASE (6)
              CALL MATERIAL_D6                                                 &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (7)
              CALL MATERIAL_D7                                                 &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (8)
              CALL MATERIAL_D8                                                 &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            CASE (9)
              CALL MATERIAL_D9                                                 &
     &          (                                                              &
     &          DAMPER(N)%RES%Force,                                           &
     &          STATE_VARIABLES(DAMPER(N)%PAR%Isv),                            &
     &          DAMPER(N)%RES%Int_Eng,                                         &
     &          DAMPER(N)%RES%DTnext,                                          &
     &          MatID                                                          &
     &          )
            END SELECT
!!
!! Compute force divergence and time step.
!!
            CALL DAMPER_FORCE_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
!!!         TIMSIM%DTDmx = MAX (TIMSIM%DTDmx,DAMPER(N)%RES%DTelt)
!!!         IF (DAMPER(N)%RES%DTelt .LT. TIMSIM%DTDmp(1)) THEN
!!!           TIMSIM%DTDmp(1) = DAMPER(N)%RES%DTelt
!!!           TIMSIM%DAMPER(1) = N
!!!         ENDIF
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
      SUBROUTINE DAMPER_VELOCITY ( NEL )
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
!! Purpose: Calculate the gradient of the axial velocity.
!!
      USE shared_common_data
      USE damper_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0)) :: Vx(2),Vy(2),Vz(2)
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
!! Define a basis vector R aligned with the damper.
!!
      IF (DAMPER(NEL)%PAR%Idir .EQ. 0) THEN
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
!! If damper length is vanishingly small, use last recorded damper direction.
!!
        IF (SLength .LT. 1.0D-25) THEN
          Rx = DAMPER(NEL)%RES%Axis(1)
          Ry = DAMPER(NEL)%RES%Axis(2)
          Rz = DAMPER(NEL)%RES%Axis(3)
          SLength = ONE
        ENDIF
        Rx = Rx * (ONE / SLength)
        Ry = Ry * (ONE / SLength)
        Rz = Rz * (ONE / SLength)
!!
!! Save current damper orientation (unit vector aligned with damper) for use
!! in the event that the damper length goes to zero in the future.
!!
        DAMPER(NEL)%RES%Axis(1) = Rx
        DAMPER(NEL)%RES%Axis(2) = Ry
        DAMPER(NEL)%RES%Axis(3) = Rz
!!
!! This is a damper with a fixed orientation; retrieve orientation.
!!
      ELSE
        Rx = DAMPER(NEL)%RES%Axis(1)
        Ry = DAMPER(NEL)%RES%Axis(2)
        Rz = DAMPER(NEL)%RES%Axis(3)
      ENDIF
!!
!! Relative axial velocity using Vx,Vy,Vz transformed to local R coordinate.
!!
      Drr = Rx*(Vx(2)-Vx(1)) + Ry*(Vy(2)-Vy(1)) + Rz*(Vz(2)-Vz(1))
!!
!! Update damper velocity.
!!
      DAMPER(NEL)%RES%Delta = Drr
!!
      RETURN
      END
!!_
      SUBROUTINE DAMPER_FORCE_DIVERGENCE ( NEL )
!!
!! Copyright (c) by KEY Associates; 26-JUL-1991 18:50:59
!!
!! Purpose: Calculate the critical time step and the divergence of the
!! forces.
!!
!!!     USE shared_common_data
      USE damper_
!!!     USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
!! Retrieve material properties.
!!
!!!     MatID = DAMPER(NEL)%PAR%MatID
!!!     Density = MATERIAL(MatID)%PVAL(1)
!!
!! Critical time step calculation.
!!
!!!     IF (SOUND_SPEED%RCL2 .EQ. 0.0) THEN
!!!       DAMPER(NEL)%RES%DTelt = 1.0D-6
!!!     ELSE
!!!       DAMPER(NEL)%RES%DTelt = SQRT (SOUND_SPEED%RCL2/Density)
!!!     ENDIF
!!
!! Accumulate element force divergence results in local arrays.
!!
      Fr = DAMPER(NEL)%RES%Force
      Fx = Rx * Fr
      Fy = Ry * Fr
      Fz = Rz * Fr
      DAMPER(NEL)%RES%Xint(1) = -Fx
      DAMPER(NEL)%RES%Yint(1) = -Fy
      DAMPER(NEL)%RES%Zint(1) = -Fz
      DAMPER(NEL)%RES%Xint(2) = +Fx
      DAMPER(NEL)%RES%Yint(2) = +Fy
      DAMPER(NEL)%RES%Zint(2) = +Fz
!!
      RETURN
      END
