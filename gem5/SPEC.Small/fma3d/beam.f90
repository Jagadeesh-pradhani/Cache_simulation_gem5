      SUBROUTINE BEAM_INITIALIZATION
!!
!! Coded by: S W Key, 20-JUL-1989 21:03:00
!! Migrated by: S W Key, 4-MAY-1991 14:06:57
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:06:34.00
!!
!! Purpose: For beam elements evaluate mass terms, find initial critical
!! time step, and divergence of the initial stresses
!!
      USE shared_common_data
      USE beam_
      USE motion_
      USE material_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: MType,MatID
!!
      LOGICAL :: FOUND,Defined_Material
!!
      DO N = 1,NUMBM
!!
!! Gather element motion.
!!
        DO i = 1,4
          EMOTION(i) = MOTION(BEAM(N)%PAR%IX(i))
        ENDDO
!!
!! Initialize element clock and time step
!!
        BEAM(N)%RES%Time   = 0.0
        BEAM(N)%RES%DTnext = 0.0
!!
!! Initialize beam stresses and state variables.
!!
        DO i = 1,16
          BEAM(N)%RES%Axial(i) = 0.0
          BEAM(N)%RES%Shear(i) = 0.0
        ENDDO
        BEAM(N)%RES%Length  = 0.0
        BEAM(N)%RES%Int_Eng = 0.0
        BEAM(N)%RES%Trs     = 0.0
        BEAM(N)%RES%Trt     = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies,
!! only the element length calculation is required.)
!!
        CALL BEAM_GRADIENT_OPERATOR ( NEL )
!!
!! Save initial element length for later axial strain calculations.
!!
        BEAM(N)%PAR%Length = BEAM(N)%RES%Length
!!
!! Check for a zero-length beam.
!!
        IF (BEAM(N)%RES%Length .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') BEAM(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'BEAM_INITIALIZATION.001.00'//                           &
     &          MSGL//'BEAM (2-Node Beam) Element ID:'//MSG1//                 &
     &          MSGL//'Has A Zero Or Negative Length.'                         &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL BEAM_MASS ( NEL )
!!
!! Retrieve material ID from element data structure.
!!
        MatID = BEAM(N)%PAR%MatID
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (BEAM(N)%PAR%ParID .LT. 0) THEN
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
          IF (Density .EQ. 0.0) THEN
            BEAM(N)%RES%DTelt = HUGE ( BEAM(N)%RES%DTelt )
          ELSE
            Ymod = MATERIAL(MatID)%PVAL(6)
            BEAM(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
          ENDIF
!!
        ELSE
!!
!! Find initial sound speeds and stress.
!!
          Mtype = MATERIAL(MatID)%Type

          Defined_Material = Mtype .EQ. 50                                     &
     &                  .OR. Mtype .EQ. 51                                     &
     &                  .OR. Mtype .EQ. 55                                     &
     &                  .OR. Mtype .EQ. 57

          IF (Defined_Material) THEN
            CALL BEAM_STRESS_INTEGRATION                                       &
     &          ( NEL, STATE_VARIABLES(BEAM(N)%PAR%Isv) )
          ELSE
            WRITE (MSG1,'(I8)') BEAM(N)%PAR%EleID
            WRITE (MSG2,'(I8)') MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'BEAM_INITIALIZITION.002.00'//                           &
     &          MSGL//'Beam (BEAM) Element ID:'//MSG1//                        &
     &          MSGL//'References A Material ID:'//MSG2//                      &
     &          MSGL//'With An Invalid Material Type;'//                       &
     &          MSGL//'Type Not Equal To 50, 51, 55 or 57.'                    &
     &          )
          ENDIF
!!
!! Compute initial stress divergence and time step.
!!
          CALL BEAM_STRESS_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTBmx = MAX (TIMSIM%DTBmx,BEAM(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = BEAM(N)%RES%DTelt .LT. TIMSIM%DTBms(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTBms(j + 1) = TIMSIM%DTBms(j)
              TIMSIM%Beams(j + 1) = TIMSIM%Beams(j)
            ENDDO
          ENDIF
          TIMSIM%DTBms(i) = BEAM(N)%RES%DTelt
          TIMSIM%Beams(i) = N
        ENDIF
!!
!! End of element do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE BEAM_INTERNAL_FORCES
!!
!! Coded by: S W Key, 20-JUL-1989 21:03:00
!! Migrated by: S W Key, 4-MAY-1991 14:06:57
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:27:55.00
!!
!! Purpose: For beam elements integrate constitutive model, find critical
!! time step, and find the divergence of the stresses.
!!
      USE shared_common_data
      USE beam_
      USE node_
      USE motion_
      USE material_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: MatID
!!
!! Loop over all beam elements.
!!
      DO N = 1,NUMBM
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,BEAM(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (BEAM(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            DO i = 1,4
              EMOTION(i) = MOTION(BEAM(N)%PAR%IX(i))
            ENDDO
!!
!! Count element execution.
!!
            COUNTER%BEAM = COUNTER%BEAM + 1
!!
!! Increment element clock.
!!
            BEAM(N)%RES%Time = BEAM(N)%RES%Time + BEAM(N)%RES%DTnext
!!
!! Scale nodal positions to current element time. Note that rotational
!! degrees of freedom do not require any scaling since they are not used.
!!
            QA = NODE(BEAM(N)%PAR%IX(1))%Time - BEAM(N)%RES%Time
            EMOTION(1)%Ux =                                                    &
     &      EMOTION(1)%Ux - QA * EMOTION(1)%Vx
            EMOTION(1)%Uy =                                                    &
     &      EMOTION(1)%Uy - QA * EMOTION(1)%Vy
            EMOTION(1)%Uz =                                                    &
     &      EMOTION(1)%Uz - QA * EMOTION(1)%Vz
!!
            QA = NODE(BEAM(N)%PAR%IX(2))%Time - BEAM(N)%RES%Time
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
!! Gradient operator, stretching, and rotation.
!!
            CALL BEAM_GRADIENT_OPERATOR ( NEL )
!!
!! Update stress-strain model to obtain new stress.
!!
            CALL BEAM_STRESS_INTEGRATION                                       &
     &          ( NEL, STATE_VARIABLES(BEAM(N)%PAR%Isv) )
!!
!! Compute stress divergence and time step.
!!
            CALL BEAM_STRESS_DIVERGENCE ( NEL )
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTBmx = MAX (TIMSIM%DTBmx,BEAM(N)%RES%DTelt)
          IF (BEAM(N)%RES%DTelt .LT. TIMSIM%DTBms(1)) THEN
            TIMSIM%DTBms(1) = BEAM(N)%RES%DTelt
            TIMSIM%Beams(1) = N
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
      SUBROUTINE BEAM_GRADIENT_OPERATOR ( NEL )
!!
!! Coded by: S W Key, 22-JUL-1989 18:53:25
!! Migrated by: S W Key, 4-MAY-1991 14:40:24
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:15:43.00
!!
!! Purpose: Calculate the mean gradients of the velocity and the rotational
!! rate, update the "z-axis" vector which defines the orientation of the beam
!! cross section, calculate the length and cross-sectional area of the beam,
!! and the transverse shear correction factor.
!!
      USE shared_common_data
      USE beam_
      USE motion_
      USE section_1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL  ! Current beam index.
!!
!! Local variables.
      INTEGER :: MatID,SecID
      REAL(KIND(0D0))                                                          &
     &          Vx(2),Vy(2),Vz(2),X(2),Y(2),Z(2),Wx(2),Wy(2),Wz(2)
!!
      COMMON /BEAM/                                                            &
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Vrr,Vsr,Vtr,Wrr,Wsr,Wtr,            &
     &          Ws,Wt,S_Factor,T_Factor
!!
      REAL(KIND(0D0)), PARAMETER :: C2F = 5.0000000000000000D-1  !  1/2!
      REAL(KIND(0D0)), PARAMETER :: C3F = 1.6666666666666667D-1  !  1/3!
      REAL(KIND(0D0)), PARAMETER :: C4F = 4.1666666666666667D-2  !  1/4!
      REAL(KIND(0D0)), PARAMETER :: C5F = 8.3333333333333333D-3  !  1/5!
      REAL(KIND(0D0)), PARAMETER :: C6F = 1.3888888888888889D-3  !  1/6!

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
        SQRT6o1 = SQRT (6.0D+0)
        SQRT5o6 = SQRT (5.0D+0 / 6.0D+0)
        FIRST = .FALSE.
     ENDIF
!!
!! Current position of nodal points, and current translational and
!! rotational velocities.
!!
      DO i = 1,2
        X(i)  = EMOTION(i)%Px + EMOTION(i)%Ux
        Y(i)  = EMOTION(i)%Py + EMOTION(i)%Uy
        Z(i)  = EMOTION(i)%Pz + EMOTION(i)%Uz
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
        Wx(i) = EMOTION(i+2)%Vx
        Wy(i) = EMOTION(i+2)%Vy
        Wz(i) = EMOTION(i+2)%Vz
      ENDDO
!!
!! Define a basis vector R aligned with the beam. Note: The beam
!! coordinate system is defined by the orthonormal basis system R,S,T
!! to distinguish it from the global X,Y,Z-coordinates. The "z-axis"
!! vector which tracks the orientation of the beam cross section is the
!! T basis vector. The "z-axis" vector defines the "up" direction for
!! the beam cross section. The "z-axis" vector gets its name from the
!! special case when the local and global coordinate systems coincide.
!!
      Rx = X(2) - X(1)
      Ry = Y(2) - Y(1)
      Rz = Z(2) - Z(1)
      Blength = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (ONE / Blength)
      Ry = Ry * (ONE / Blength)
      Rz = Rz * (ONE / Blength)
      BEAM(NEL)%RES%Length = Blength
!!
!! Retrieve components of "z-axis" vector.
!!
      Tx = BEAM(NEL)%RES%Zaxis(1)
      Ty = BEAM(NEL)%RES%Zaxis(2)
      Tz = BEAM(NEL)%RES%Zaxis(3)
!!
!! Mean rotation rate.
!!
      Ox = 0.5 * (Wx(2) + Wx(1))
      Oy = 0.5 * (Wy(2) + Wy(1))
      Oz = 0.5 * (Wz(2) + Wz(1))
      Wr = Ox*Rx + Oy*Ry + Oz*Rz
!!
!! Update orientation of "z-axis" vector. First, reorthogonalize and
!! normalize "z-axis" vector. Second, rotate "z-axis" vector (which
!! defines the "up" direction of the beam cross section) using the mean
!! rotation rate about the beam axis Wr = W dot R.
!!
      Rmag = Tx*Rx + Ty*Ry + Tz*Rz
      Tx = Tx - Rmag * Rx
      Ty = Ty - Rmag * Ry
      Tz = Tz - Rmag * Rz
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
      Tx = Tx * (ONE / Tmag)
      Ty = Ty * (ONE / Tmag)
      Tz = Tz * (ONE / Tmag)
!!
      IF (ABS(BEAM(NEL)%RES%DTnext*Wr) .GT. 1.0D-18) THEN
!!
!! Compute angle of rotation and evaluate the trigonometric functions
!! SinPhi = Sin(Phi) and CosPm1 = Cos(Phi) - 1.
!!
        Phi = BEAM(NEL)%RES%DTnext*Wr
        Ph2 = Phi*Phi
        CosPm1 = -Ph2*(C2F - Ph2*(C4F - C6F*Ph2))
        SinPhi =  Phi*(ONE - Ph2*(C3F - C5F*Ph2))
!!
!! Compute S' = T x R. The vector S' is orthogonal to the plane defined
!! by T and R. (This is a temporary construction for use in rotating T.)
!!
        Sx = Ty*Rz - Tz*Ry
        Sy = Tz*Rx - Tx*Rz
        Sz = Tx*Ry - Ty*Rx
!!
!! Compute incremental rotation of the "z-axis" unit vector (Tx,Ty,Tz).
!! The rotation used is an incremental, proper orthogonal transformation.
!!
        Tx = Tx + Tx*CosPm1 - Sx*SinPhi
        Ty = Ty + Ty*CosPm1 - Sy*SinPhi
        Tz = Tz + Tz*CosPm1 - Sz*SinPhi
!!
      ENDIF
!!
!! Update "z-axis" unit vector.
!!
      BEAM(NEL)%RES%Zaxis(1) = Tx
      BEAM(NEL)%RES%Zaxis(2) = Ty
      BEAM(NEL)%RES%Zaxis(3) = Tz
!!
!! Construct the "y-axis" unit vector S orthogonal to R and T; S = T x R.
!!
      Sx = Ty*Rz - Tz*Ry
      Sy = Tz*Rx - Tx*Rz
      Sz = Tx*Ry - Ty*Rx
!!
!! Axial gradients using relative velocity dVx,dVy,dVz transformed to
!! local r,s,t-coordinates.
!!
      dVx = Vx(2) - Vx(1)
      dVy = Vy(2) - Vy(1)
      dVz = Vz(2) - Vz(1)
      Vrr = (Rx*dVx + Ry*dVy + Rz*dVz) * (ONE / Blength)
      Vsr = (Sx*dVx + Sy*dVy + Sz*dVz) * (ONE / Blength)
      Vtr = (Tx*dVx + Ty*dVy + Tz*dVz) * (ONE / Blength)
!!
!! Axial gradients using rotational rate Wx,Wy,Wz transformed to local
!! r,s,t-coordinates.
!!
      dWr = Rx*(Wx(2)-Wx(1)) + Ry*(Wy(2)-Wy(1)) + Rz*(Wz(2)-Wz(1))
      dWs = Sx*(Wx(2)-Wx(1)) + Sy*(Wy(2)-Wy(1)) + Sz*(Wz(2)-Wz(1))
      dWt = Tx*(Wx(2)-Wx(1)) + Ty*(Wy(2)-Wy(1)) + Tz*(Wz(2)-Wz(1))
      Wrr = dWr * (ONE / Blength)
      Wsr = dWs * (ONE / Blength)
      Wtr = dWt * (ONE / Blength)
!!
!! Calculate average rotations in local r,s,t-coordinates.
!!
      Ws = Ox*Sx + Oy*Sy + Oz*Sz
      Wt = Ox*Tx + Oy*Ty + Oz*Tz
!!
!! Transverse shear and slender beam limit correction factor.
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
!!        verse shear strain energy, namely, 6.0*(thickness/length)**2.
!!        The factor allows sufficient transverse shear stiffness to
!!        obtain thin shell results without severe reduction in critical
!!        time step or locking in the element. Two factors are used; one
!!        in the "y-axis" direction and one in the "z-axis" direction.
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
      SecID = BEAM(NEL)%PAR%SecID
      Aw = SECTION_1D(SecID)%Width / Blength
      S_Factor = MIN (SQRT5o6, SQRT6o1*Aw)
      Ah = SECTION_1D(SecID)%Height / Blength
      T_Factor = MIN (SQRT5o6, SQRT6o1*Ah)
!!
      RETURN
      END
!!_
      SUBROUTINE BEAM_MASS ( NEL )
!!
!! Coded by: S W Key, 4-MAY-1991 15:30:13
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:27:31.00
!!
!! Purpose: For the elastic beam element calculate the translational and
!! rotational masses at the nodes. The simplest of mass lumpings is used:
!! one half at each node and 1/12-th times the length squared times the
!! half mass.
!!
      USE shared_common_data
      USE beam_
      USE node_
      USE motion_
      USE material_
      USE section_1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL
!!
!! Local variables.
      INTEGER :: MatID,SecID
!!
!! Retrieve current beam element's material ID and section ID.
!!
      MatID = BEAM(NEL)%PAR%MatID
      SecID = BEAM(NEL)%PAR%SecID
!!
!! Compute one half total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      Qmass = 0.5*Density*BEAM(NEL)%PAR%Length*SECTION_1D(SecID)%Area
!!
      DO i = 1,2
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
!!
!! Distribute mass to the element's nodal points.
!!
        NODE(BEAM(NEL)%PAR%IX(i))%Mass =                                       &
     &    NODE(BEAM(NEL)%PAR%IX(i))%Mass + Qmass
      ENDDO
!!
!! Compute nodal isotropic inertia
!!
      Rmass = Qmass * BEAM(NEL)%PAR%Length**2 / 12.0
!!
      DO i = 3,4
        NODE(BEAM(NEL)%PAR%IX(i))%Mass =                                       &
     &    NODE(BEAM(NEL)%PAR%IX(i))%Mass + Rmass
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE BEAM_STRESS_DIVERGENCE ( NEL )
!!
!! Coded by: S W Key, 22-JUL-1989 19:05:24
!! Migrated by: S W Key, 4-MAY-1991 15:39:50
!! Copyright (c) by KEY Associates;  8-AUG-1993 13:59:43.23
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:28:35.00
!!
!! Purpose: Calculate the critical time step and the divergence of the
!! stresses.
!!
      USE shared_common_data
      USE beam_
      USE material_
      USE section_1d_
      USE mean_stress_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL
!!
!! Local variables.
      INTEGER :: MatID,SecID
!!
      COMMON /BEAM/                                                            &
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Vrr,Vsr,Vtr,Wrr,Wsr,Wtr,            &
     &          Ws,Wt,S_Factor,T_Factor
!!
!! Retrieve material properties.
!!
      SecID = BEAM(NEL)%PAR%SecID
      MatID = BEAM(NEL)%PAR%MatID
      Density = MATERIAL(MatID)%PVAL(1)
      Bulk_Ln = MATERIAL(MatID)%PVAL(2)
      Bulk_Qd = MATERIAL(MatID)%PVAL(3)
      HG_Visc = MATERIAL(MatID)%PVAL(4)
!!
!! Check for stiffness-only element. (This functionality is provided
!! only because structural analysis systems based on implicit time
!! integration algorithms such as NASTRAN allow these abberations.)
!!
      IF (Density .EQ. 0.0) THEN
        QP = 0.0
        BEAM(NEL)%RES%DTelt = HUGE ( BEAM(NEL)%RES%DTelt )
      ELSE
!!
!! Find sound speed.
!!
        CSQ = SOUND_SPEED%RCL2 / Density
        IF (CSQ .GT. 0.0) THEN
          Cv = SQRT(CSQ)
          Cs = SQRT(SOUND_SPEED%RCS2 / Density)
        ELSE
          WRITE (MSG1,'(I8)') BEAM(NEL)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'BEAM_STRESS_DIVERGENCE.001.00'//                        &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'//                &
     &          MSGL//'Element ID:'//MSG1                                      &
     &          )
          Cv = TIMSIM%DTlast * 1.0D-6
          Cs = TIMSIM%DTlast * 1.0D-6
        ENDIF
!!
!! Artificial bulk viscosity.
!!
        IF (Vrr .LT. 0.0) THEN
          Beam_Length = BEAM(NEL)%PAR%Length
          QC = Bulk_Qd * Bulk_Qd * Beam_Length * ABS(Vrr) + Bulk_Ln * Cv
          QP = Density * Beam_Length * Vrr * QC
        ELSE
          QC = 0.0
          QP = 0.0
        ENDIF
!!
!! Critical time step calculation. Three seperate conditions are considered:
!! (1) an axial critical time step DTmb, (2) a transverse shear critical
!! time step DTts, and (3) a bending critical time step DTbd.
!!
        DTmb = BEAM(NEL)%RES%Length / (QC + SQRT (QC*QC + Cv*Cv))
        DTts = BEAM(NEL)%RES%Length / (MAX (S_Factor,T_Factor) * Cs)
        DTbd = BEAM(NEL)%RES%Length * BEAM(NEL)%RES%Length * 0.5 /             &
     &     (MAX (SECTION_1D(SecID)%Height,SECTION_1D(SecID)%Width) * Cv)
        BEAM(NEL)%RES%DTelt = MIN (DTmb,DTts,DTbd)
      ENDIF
!!
!! Accumulate element stress divergence results in local arrays.
!!
      Volume = SECTION_1D(SecID)%Area * BEAM(NEL)%PAR%Length
      Area = Volume / BEAM(NEL)%RES%Length
      Fr = -Area * (MEAN%Nrr + Qp)
      Fs = -Area * (MEAN%Nrs + BEAM(NEL)%RES%Trs * S_Factor)
      Ft = -Area * (MEAN%Nrt + BEAM(NEL)%RES%Trt * T_Factor)
      Fx = Rx * Fr + Sx * Fs + Tx * Ft
      Fy = Ry * Fr + Sy * Fs + Ty * Ft
      Fz = Rz * Fr + Sz * Fs + Tz * Ft
      BEAM(NEL)%RES%Xint(1) =  Fx
      BEAM(NEL)%RES%Yint(1) =  Fy
      BEAM(NEL)%RES%Zint(1) =  Fz
      BEAM(NEL)%RES%Xint(2) = -Fx
      BEAM(NEL)%RES%Yint(2) = -Fy
      BEAM(NEL)%RES%Zint(2) = -Fz
      Fr = Area * MEAN%Mrr
      Fs = Area * MEAN%Mrs
      Ft = Area * MEAN%Mrt
      Fx = Rx * Fr + Sx * Fs + Tx * Ft
      Fy = Ry * Fr + Sy * Fs + Ty * Ft
      Fz = Rz * Fr + Sz * Fs + Tz * Ft
      Srs = (0.5 * Volume) * (MEAN%Nrs + BEAM(NEL)%RES%Trs * S_Factor)
      Srt = (0.5 * Volume) * (MEAN%Nrt + BEAM(NEL)%RES%Trt * T_Factor)
      S1 = Sx * Srt - Tx * Srs
      S2 = Sy * Srt - Ty * Srs
      S3 = Sz * Srt - Tz * Srs
      BEAM(NEL)%RES%Xint(3) = S1 + Fx
      BEAM(NEL)%RES%Yint(3) = S2 + Fy
      BEAM(NEL)%RES%Zint(3) = S3 + Fz
      BEAM(NEL)%RES%Xint(4) = S1 - Fx
      BEAM(NEL)%RES%Yint(4) = S2 - Fy
      BEAM(NEL)%RES%Zint(4) = S3 - Fz
!!
      RETURN
      END
!!_
      SUBROUTINE BEAM_STRESS_INTEGRATION ( NEL,STATE_VARIABLES )
!!
!! Coded by: S W Key, 2-APR-1989 13:39:23
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:28:56.00
!!
!! Purpose: Increment beam stresses using the current values of the axial
!! velocity gradient and gradients in the nodal rates of rotation, compute
!! the "extreme" fiber strain rates, the "extreme" fiber stresses
!! BEAM(NEL)%RES%Axial(1:Npts) & .Shear(1:Npts) for post-processing, and the
!! "mean" stresses Nrr, Nrs, Nrt, Mrr, Mrs, and Mrt for use in the stress
!! divergence calculation.
!!
      USE shared_common_data
      USE beam_
      USE material_
      USE section_1d_
      USE mean_stress_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN)    :: NEL
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
!!
!! Local variables. 
      INTEGER :: MatID,SecID
!!
      COMMON /BEAM/                                                            &
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Vrr,Vsr,Vtr,Wrr,Wsr,Wtr,            &
     &          Ws,Wt,S_Factor,T_Factor
!!
      REAL(KIND(0D0))         &
     &          Spos(16,6),   & ! Section non-dimensional s-coordinate
     &          Tpos(16,6),   & ! Section non-dimensional t-coordinate
     &          Ds(16,6),     & ! Circumferential vector s-component
     &          Dt(16,6),     & ! Circumferential vector t-component
     &          Awgt(16,5),   & ! Nondimensional integration weights
     &          Bwgt(16,5),   & ! Nondimensional integration weights
     &          S(16),        & ! Section extreme fiber s-coordinates
     &          T(16),        & ! Section extreme fiber t-coordinates
     &          Cs(16),       & ! Circumferential vector, s-component
     &          Ct(16),       & ! Circumferential vector, t-component
     &          Drr(16),      & ! Axial component of the stretching
     &          Drs(16),      & ! Shearing component of the stretching
     &          Drt(16),      & ! Shearing component of the stretching
     &          Drc(16),      & ! Circumferential component of the stretching
     &          Srr(16),      & ! Axial component of the stress
     &          Srs(16),      & ! Shear component of the stress
     &          Srt(16),      & ! Shear component of the stress
     &          Awt(16),      & ! Weights in numeric integration for resultants
     &          Bwt(16),      & ! Weights in numeric integration for resultants
     &          Int_Eng
      INTEGER                 &
     &          Iel(16),      & ! El. integration segments; closed, hollow sec'
     &          Ipl(32),      & ! Pl. integration segments; closed, hollow sec'
     &          Ipr(2,16)       ! Integration segments; closed, hollow sections
      LOGICAL                                                                  &
     &          Elastic_Material

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Define constants for use in DATA statements.
!!
      PARAMETER (UP =  1.0D+0, PP =  0.923879532D+0, RP =  0.707106781D+0)
      PARAMETER (FP =  0.5D+0, ZP =  0.0D+0,         VP =  0.65465367D+0 )
      PARAMETER (QP =  0.25D+0,GP =  0.44721360D+0,  DP =  0.382683432D+0)

      PARAMETER (UN = -1.0D+0, PN = -0.923879532D+0, RN = -0.707106781D+0)
      PARAMETER (FN = -0.5D+0, ZN = -0.0D+0,         VN = -0.65465367D+0 )
      PARAMETER (QN = -0.25D+0,GN = -0.44721360D+0,  DN = -0.382683432D+0)
!!
!! I. DEFINE NON-DIMENSIONAL STRESS EVALUATION POINTS.
!!
!! Solid eliptic cross section          (Isec = 1, Twall .eq. 0.0).
!!
      DATA (Spos(i,4),i=1, 5) /UP, ZP, UN, ZP, ZP/
      DATA (Tpos(i,4),i=1, 5) /ZP, UP, ZP, UN, ZP/
      DATA (Spos(i,4),i=6, 9) /FP, FN, FN, FP/
      DATA (Tpos(i,4),i=6, 9) /FP, FP, FN, FN/
!!
!! Hollow eliptic cross section                 (Isec = 2, Twall .ne. 0.0).
!!
      DATA (Spos(i,1),i=1, 8) /UP, ZP, UN, ZP, RP, RN, RN, RP/
      DATA (Tpos(i,1),i=1, 8) /ZP, UP, ZP, UN, RP, RP, RN, RN/
      DATA (Spos(i,1),i=9,16) /PP, DP, DN, PN, PN, DN, DP, PP/
      DATA (Tpos(i,1),i=9,16) /DP, PP, PP, DP, DN, PN, PN, DN/
!!
!! Solid rectangular cross section      (Isec = 3, Twall .eq. 0.0).
!!
      DATA (Spos(i,5),i=1, 9) /UP,UN,UN,UP,ZP,UN,ZP,UP,ZP/  !  3x3 Lobatto
      DATA (Tpos(i,5),i=1, 9) /UP,UP,UN,UN,UP,ZP,UN,ZP,ZP/
      DATA (Spos(i,6),i=1, 8) /UP,UN,UN,UP,GP,GN,UN,UN/     !  4x4 Lobatto
      DATA (Tpos(i,6),i=1, 8) /UP,UP,UN,UN,UP,UP,GP,GN/
      DATA (Spos(i,6),i=9,16) /GN,GP,UP,UP,GP,GN,GN,GP/
      DATA (Tpos(i,6),i=9,16) /UN,UN,GN,GP,GP,GP,GN,GN/
!!
!! Hollow rectangular cross section     (Isec = 4, Twall .ne. 0.0).
!!
      DATA (Spos(i,2),i=1, 8) /UP, UN, UN, UP, ZP, UN, ZP, UP/
      DATA (Tpos(i,2),i=1, 8) /UP, UP, UN, UN, UP, ZP, UN, ZP/
      DATA (Spos(i,2),i=9,16) /UP, QP, QN, UN, UN, QN, QP, UP/
      DATA (Tpos(i,2),i=9,16) /QP, UP, UP, QP, QN, UN, UN, QN/
!!
!! I-beam cross section                         (Isec = 5, Twall .ne. 0.0).
!!
      DATA (Spos(i,3),i= 1, 9) /UP, UN, UN, UP, ZP, ZP, ZP, ZP, ZP/
      DATA (Tpos(i,3),i= 1, 9) /UP, UP, UN, UN, UP, UN, UP, ZP, UN/
      DATA (Spos(i,3),i=10,15) /VP, VN, VN, VP, ZP, ZP/
      DATA (Tpos(i,3),i=10,15) /UP, UP, UN, UN, VP, VN/
!!
!! II. DEFINE "CIRCUMFERENTIAL" SHEAR STRAIN RATE DIRECTIONS.
!!
!! Solid eliptic cross section          (Isec = 1, Twall .eq. 0.0).
!!
      DATA (Ds(i,4),i=1, 5) /ZP, UN, ZP, UP, ZP/
      DATA (Dt(i,4),i=1, 5) /UP, ZP, UN, ZP, ZP/
      DATA (Ds(i,4),i=6, 9) /UN, UN, UP, UP/
      DATA (Dt(i,4),i=6, 9) /UP, UN, UN, UP/
!!
!! Hollow eliptic cross section                 (Isec = 2, Twall .ne. 0.0).
!!
      DATA (Ds(i,1),i=1, 8) /ZP, UN, ZP, UP, RN, RN, RP, RP/
      DATA (Dt(i,1),i=1, 8) /UP, ZP, UN, ZP, RP, RN, RN, RP/
      DATA (Ds(i,1),i=9,16) /DN, PN, PN, DN, DP, PP, PP, DP/
      DATA (Dt(i,1),i=9,16) /PP, DP, DN, PN, PN, DN, DP, PP/
!!
!! Solid rectangular cross section      (Isec = 3, Twall .eq. 0.0).
!!
      DATA (Ds(i,5),i=1, 9) /ZP,ZP,ZP,ZP,UN,ZP,UP,ZP,ZP/  !  3x3 Lobatto
      DATA (Dt(i,5),i=1, 9) /ZP,ZP,ZP,ZP,ZP,UN,ZP,UP,ZP/
      DATA (Ds(i,6),i=1, 8) /ZP,ZP,ZP,ZP,UN,UN,ZP,ZP/     !  4x4 Lobatto
      DATA (Dt(i,6),i=1, 8) /ZP,ZP,ZP,ZP,ZP,ZP,UN,UN/
      DATA (Ds(i,6),i=9,16) /UP,UP,ZP,ZP,UN,UN,UP,UP/
      DATA (Dt(i,6),i=9,16) /ZP,ZP,UP,UP,UP,UN,UN,UP/
!!
!! Hollow rectangular cross section     (Isec = 4, Twall .ne. 0.0).
!!
      DATA (Ds(i,2),i=1, 8) /FN, FN, FP, FP, UN, ZP, UP, ZP/
      DATA (Dt(i,2),i=1, 8) /FP, FN, FN, FP, ZP, UN, ZP, UP/
      DATA (Ds(i,2),i=9,16) /ZP, UN, UN, ZP, ZP, UP, UP, ZP/
      DATA (Dt(i,2),i=9,16) /UP, ZP, ZP, UN, UN, ZP, ZP, UP/
!!
!! I-beam cross section                         (Isec = 5, Twall .ne. 0.0).
!!
      DATA (Ds(i,3),i= 1, 9) /ZP, ZP, ZP, ZP, UN, UP, ZP, ZP, ZP/
      DATA (Dt(i,3),i= 1, 9) /ZP, ZP, ZP, ZP, ZP, ZP, ZP, UP, ZP/
      DATA (Ds(i,3),i=10,15) /UN, UN, UP, UP, ZP, ZP/
      DATA (Dt(i,3),i=10,15) /ZP, ZP, ZP, ZP, UP, UP/
!!
!! III. DEFINE NON-DIMENSIONAL NUMERICAL INTEGRATION WEIGHTS.
!!
!! Solid eliptic cross section          (Isec = 1, Twall .eq. 0.0).
!!
      DATA (Awgt(i,4),i=1, 9) /4*0.25D+0, 5*ONE/
!!
!! Hollow eliptic cross section                 (Isec = 2, Twall .ne. 0.0).
!!
      DATA (Awgt(i,1),i=1,16) /16*ONE/
!!
      DATA (Iel(i),i= 1,16)/1,5,5,2,2,6,6,3,3,7,7,4,4,8,8,1/
      DATA (Ipl(i),i= 1,16)/1,9,9,5,5,10,10,2,2,11,11,6,6,12,12,3/
      DATA (Ipl(i),i=17,32)/3,13,13,7,7,14,14,4,4,15,15,8,8,16,16,1/
!!
!! Solid rectangular cross section      (Isec = 3, Twall .eq. 0.0).
!!
      DATA (Awgt(i,5),i=1, 9) /4*ONE, 4*4.0D+0,   16.0D+0/  !  3x3 Lobatto
      DATA (Bwgt(i,5),i=1,16) /4*ONE, 8*5.0D+0, 4*25.0D+0/  !  4x4 Lobatto
!!
!! Hollow rectangular cross section     (Isec = 4, Twall .ne. 0.0).
!!
      DATA (Awgt(i,2),i=1,16) /16*ONE/
!!
!! I-beam cross section                         (Isec = 5, Twall .ne. 0.0).
!!
      PARAMETER (W31 = 1.0D+0, W32 = 4.0D+0)
      DATA (Awgt(i,3),i= 1, 9) /4*W31, 2*W32, W31, W32, W31/
!!
      PARAMETER (W51 = 0.3D+0)
      PARAMETER (W52 = 1.633333333333333333D+0);
      PARAMETER (W53 = 2.133333333333333333D+0);
      DATA (Bwgt(i,3),i= 1,15) /4*W51, 2*W53, W51, W53, W51, 6*W52/
!!
!! Normalize the circumferential shear vectors.
!!
      IF (FIRST) THEN
        DO Nvec = 1,6
          DO i = 1,16
            Qmag = SQRT (Ds(i,Nvec)*Ds(i,Nvec) + Dt(i,Nvec)*Dt(i,Nvec))
            IF (Qmag .GT. 0.0) THEN
              Ds(i,Nvec) = Ds(i,Nvec) * (ONE / Qmag)
              Dt(i,Nvec) = Dt(i,Nvec) * (ONE / Qmag)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Retrieve cross section geometric data.
!!
      SecID   = BEAM(NEL)%PAR%SecID
      Isec    = SECTION_1D(SecID)%Section
      NPLoc   = SECTION_1D(SecID)%NPLoc
      Width   = SECTION_1D(SecID)%Width
      Height  = SECTION_1D(SecID)%Height
      Twall   = SECTION_1D(SecID)%Twall
      Tflange = SECTION_1D(SecID)%Tflange
      Srefloc = SECTION_1D(SecID)%Yrefloc
      Trefloc = SECTION_1D(SecID)%Zrefloc
!!
!! Retrieve cross section material type (constitutive model)
!!
      MatID = BEAM(NEL)%PAR%MatID
      Mtype = MATERIAL(MatID)%Type
      Elastic_Material = .NOT.(Mtype .EQ. 51)
!!
!! Select stress evaluation points based on cross section geometry and
!! elastic versus inelastic material constitutive model. Compute off-sets
!! Sos and Tos to locate evaluation points in the middle of thin-wall
!! sections. Set up weights to compute stress resultants.
!!
!! Solid eliptic cross section          (Isec = 1, Twall .eq. 0.0).
!!
      IF (Isec .EQ. 1) THEN
        IF (Elastic_Material) THEN
          Npts = 5
          Nsec = 4
          Nvec = 4
          Nmom = 2
          DO i = 1,Npts
            Awt(i) = (Awgt(i,4)* (ONE/2.0D+0))
            Bwt(i) = (Awgt(i,4)* (ONE/2.0D+0))
          ENDDO
        ELSE
          Npts = 9
          Nsec = 4
          Nvec = 4
          Nmom = 2
          DO i = 1,Npts
            Awt(i) = (Awgt(i,4)* (ONE/6.0D+0))
            Bwt(i) = (Awgt(i,4)* (ONE/6.0D+0))
          ENDDO
        ENDIF
        Sos = 0.0
        Tos = 0.0
      ENDIF
!!
!! Hollow eliptic cross section                 (Isec = 2, Twall .ne. 0.0).
!!
      IF (Isec .EQ. 2) THEN
        IF (Elastic_Material) THEN
          Npts = 8
          Nsec = 1
          Nvec = 1
          Nmom = 1
          DO i = 1,Npts
            Awt(i) = (Awgt(i,1)* (ONE/ 8.0D+0))
            Bwt(i) = (Awgt(i,1)* (ONE/48.0D+0))
          ENDDO
          DO i = 1,2*Npts
            Ipr(i,1) = Iel(i)
          ENDDO
        ELSE
          Npts = 16
          Nsec = 1
          Nvec = 1
          Nmom = 1
          DO i = 1,Npts
            Awt(i) = (Awgt(i,1)* (ONE/16.0D+0))
            Bwt(i) = (Awgt(i,1)* (ONE/96.0D+0))
          ENDDO
          DO i = 1,2*Npts
            Ipr(i,1) = Ipl(i)
          ENDDO
        ENDIF
        Sos = (0.5 * Twall)
        Tos = (0.5 * Twall)
      ENDIF
!!
!! Solid rectangular cross section      (Isec = 3, Twall .eq. 0.0).
!!
      IF (Isec .EQ. 3) THEN
        IF (Elastic_Material) THEN
          Npts = 9
          Nsec = 5
          Nvec = 5
          Nmom = 2
          DO i = 1,Npts
            Awt(i) = (Awgt(i,5)* (ONE/36.0D+0))
            Bwt(i) = (Awgt(i,5)* (ONE/36.0D+0))
          ENDDO
        ELSE
          Npts = 16
          Nsec = 6
          Nvec = 6
          Nmom = 2
          DO i = 1,Npts
            Awt(i) = (Bwgt(i,5)* (ONE/144.0D+0))
            Bwt(i) = (Bwgt(i,5)* (ONE/144.0D+0))
          ENDDO
        ENDIF
        Sos = 0.0
        Tos = 0.0
      ENDIF
!!
!! Hollow rectangular cross section     (Isec = 4, Twall .ne. 0.0).
!!
      IF (Isec .EQ. 4) THEN
        IF (Elastic_Material) THEN
          Npts = 8
          Nsec = 2
          Nvec = 2
          Nmom = 1
          DO i = 1,Npts
            Awt(i) = (Awgt(i,2)* (ONE/ 8.0D+0))
            Bwt(i) = (Awgt(i,2)* (ONE/48.0D+0))
          ENDDO
          DO i = 1,2*Npts
            Ipr(i,1) = Iel(i)
          ENDDO
        ELSE
          Npts = 16
          Nsec = 2
          Nvec = 2
          Nmom = 1
          DO i = 1,Npts
            Awt(i) = (Awgt(i,2)* (ONE/16.0D+0))
            Bwt(i) = (Awgt(i,2)* (ONE/96.0D+0))
          ENDDO
          DO i = 1,2*Npts
            Ipr(i,1) = Ipl(i)
          ENDDO
        ENDIF
        Sos = (0.5 * Twall)
        Tos = (0.5 * Twall)
      ENDIF
!!
!! I-beam cross section                         (Isec = 5, Twall .ne. 0.0).
!!
      IF (Isec .EQ. 5) THEN
        Area = SECTION_1D(SecID)%Area
        IF (Elastic_Material) THEN
          Npts = 9
          Nsec = 3
          Nvec = 3
          Nmom = 2
          DO i = 1,Npts-3
            Awt(i) = Awgt(i,3) * ((Width*Tflange)/(6.0*Area))
          ENDDO
          Awt(7) = Awgt(7,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt(8) = Awgt(8,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt(9) = Awgt(9,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          DO i = 1,Npts
            Bwt(i) = Awt(i)
          ENDDO
        ELSE
          Npts = 15
          Nsec = 3
          Nvec = 3
          Nmom = 2
          DO i = 1,Npts-2
            Awt(i) = Bwgt(i,3)* ((Width*Tflange)/(6.0*Area))
          ENDDO
          Awt( 7) = Bwgt( 7,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt( 8) = Bwgt( 8,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt( 9) = Bwgt( 9,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt(14) = Bwgt(14,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          Awt(15) = Bwgt(15,3) * ((Height-2.0*Tflange)*Twall/(6.0*Area))
          DO i = 1,Npts
            Bwt(i) = Awt(i)
          ENDDO
        ENDIF
        Sos = 0.0
        Tos = (0.5 * Tflange)
      ENDIF
!!
!! Locate beam reference axis for off-set cross sections.
!! (NPLoc = 1, Right; NPLoc = 2, Top; NPLoc = 3, Left; NPLoc = 4, Bottom)
!!
      IF (NPLoc .EQ. 0) THEN
        Sra = 0.0
        Tra = 0.0
      ELSE IF (NPLoc .EQ. 1) THEN
        Sra = (0.5 * Width)
        Tra = 0.0
      ELSE IF (NPLoc .EQ. 2) THEN
        Sra = 0.0
        Tra = (0.5 * Height)
      ELSE IF (NPLoc .EQ. 3) THEN
        Sra =-(0.5 * Width)
        Tra = 0.0
      ELSE IF (NPLoc .EQ. 4) THEN
        Sra = 0.0
        Tra =-(0.5 * Height)
      ELSE IF (NPLoc .EQ. 5) THEN
        Sra = Srefloc
        Tra = Trefloc
      ENDIF
!!
!! Compute cross section positions and "circumferential" shear directions.
!!
      DO i = 1,Npts
        S(i) = (0.5 * Width  - Sos)
        T(i) = (0.5 * Height - Tos)
      ENDDO
      IF (Isec .EQ. 5) THEN
        T( 7) = (0.5 * Height - Tflange)
        T( 8) = (0.5 * Height - Tflange)
        T( 9) = (0.5 * Height - Tflange)
        T(14) = (0.5 * Height - Tflange)
        T(15) = (0.5 * Height - Tflange)
      ENDIF
      DO i = 1,Npts
        S(i)  = S(i) * Spos(i,Nsec) - Sra
        T(i)  = T(i) * Tpos(i,Nsec) - Tra
        Cs(i) = Ds(i,Nvec)
        Ct(i) = Dt(i,Nvec)
      ENDDO
!!
!! Compute extensional stretching rates. The stretching rate given by
!! Drr(1:Npts) is the extensional stretching evaluated at specific points
!! in the cross section for use in evaluating the mean-stress integrals.
!!
      DO i = 1,Npts
        Drr(i) = Vrr + T(i) * Wsr - S(i) * Wtr
      ENDDO
!!
!! Compute shearing rates. The shearing rates Drs(1:Npts) and Drt(1:Npts) are
!! the shearing components of the stretching evaluated at specific locations
!! within the cross section for use in evaluating the mean stress integrals.
!! Drs(1:Npts) and Drt(1:Npts) are defined in terms of the twisting Wrr. Ert
!! and Ers are the transverse shear stretchings; they are kept separate from
!! the twisting since for thin beams they are important numerical parameters
!! and are treated differently.
!!
      Ers = -(0.5 * (Wt - Vsr) * S_Factor)
      Ert =  (0.5 * (Ws + Vtr) * T_Factor)
      DO i = 1,Npts
        Drs(i) = -T(i) * (0.5 * Wrr)
        Drt(i) = +S(i) * (0.5 * Wrr)
      ENDDO
!!
!! Rotate the shearing components of the stretching to obtain the
!! "circumferential" component Drc(1:Npts). (The "radial" component
!! Dab(1:Npts) is assumed to vanish.)
!!
      DO i = 1,Npts
        Drc(i) = Drt(i) * Ct(i) + Drs(i) * Cs(i)
      ENDDO
!!
!! Internal energy increment from time n to time n+1/2.
!!
      Int_Eng = 0.0
      DO i = 1,Npts
        Int_Eng = Int_Eng                                                      &
     &          + Awt(i) * (Drr(i)*BEAM(NEL)%RES%Axial(i)                      &
     &          + (Drc(i)+Drc(i))*BEAM(NEL)%RES%Shear(i))
      ENDDO
!!
!! Calculate new stresses; incremental constitutive evaluation.
!!
      SELECT CASE (Mtype)
        CASE (50)
          CALL MATERIAL_50 ( NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts )
        CASE (51)
          CALL MATERIAL_51 ( NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts )
        CASE (55)
          CALL MATERIAL_55 ( NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts )
        CASE (57)
          CALL MATERIAL_57 ( NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts )
      END SELECT
!!
!! Internal energy increment from time n+1/2 to time n+1.
!!
      DO i = 1,Npts
        Int_Eng = Int_Eng                                                      &
     &          + Awt(i) * (Drr(i)*BEAM(NEL)%RES%Axial(i)                      &
     &          + (Drc(i)+Drc(i))*BEAM(NEL)%RES%Shear(i))
      ENDDO
      BEAM(NEL)%RES%Int_Eng =                                                  &
     &  BEAM(NEL)%RES%Int_Eng + (0.5*BEAM(NEL)%RES%DTnext*Int_Eng)
!!
!! Obtain axial and shear stresses in r,s,t-coordinates.
!!
      DO i = 1,Npts
        Srr(i) = BEAM(NEL)%RES%Axial(i)
        Srs(i) = BEAM(NEL)%RES%Shear(i) * Cs(i)
        Srt(i) = BEAM(NEL)%RES%Shear(i) * Ct(i)
      ENDDO
!!
!! Evaluate Mean_Stresses (= Stress_Resultants / Area)
!!
      MEAN%Nrr = 0.0
      MEAN%Nrs = 0.0
      MEAN%Nrt = 0.0
      MEAN%Mrr = 0.0
      MEAN%Mrs = 0.0
      MEAN%Mrt = 0.0
!!
!! Integrate over section to obtain mean normal resultants.
!!
      DO i = 1,Npts
        MEAN%Nrr = MEAN%Nrr + Awt(i) * Srr(i)
        MEAN%Nrs = MEAN%Nrs + Awt(i) * Srs(i)
        MEAN%Nrt = MEAN%Nrt + Awt(i) * Srt(i)
      ENDDO
!!
!! Integrate over section to obtain mean moment resultants.
!!
      IF (Nmom .EQ. 1) THEN
        DO i = 1,Npts
          m = Ipr(1,i)
          n = Ipr(2,i)
          MEAN%Mrr = MEAN%Mrr + (Bwt(m) * Srs(m)) * (T(m)+T(m)+T(n))
          MEAN%Mrr = MEAN%Mrr + (Bwt(n) * Srs(n)) * (T(n)+T(n)+T(m))
          MEAN%Mrr = MEAN%Mrr - (Bwt(m) * Srt(m)) * (S(m)+S(m)+S(n))
          MEAN%Mrr = MEAN%Mrr - (Bwt(n) * Srt(n)) * (S(n)+S(n)+S(m))
          MEAN%Mrs = MEAN%Mrs - (Bwt(m) * Srr(m)) * (T(m)+T(m)+T(n))
          MEAN%Mrs = MEAN%Mrs - (Bwt(n) * Srr(n)) * (T(n)+T(n)+T(m))
          MEAN%Mrt = MEAN%Mrt + (Bwt(m) * Srr(m)) * (S(m)+S(m)+S(n))
          MEAN%Mrt = MEAN%Mrt + (Bwt(n) * Srr(n)) * (S(n)+S(n)+S(m))
        ENDDO
      ELSE IF (Nmom .EQ. 2) THEN
        DO i = 1,Npts
          MEAN%Mrr = MEAN%Mrr + Bwt(i) * (Srs(i) * T(i) - Srt(i) * S(i))
          MEAN%Mrs = MEAN%Mrs - Bwt(i) *  Srr(i) * T(i)
          MEAN%Mrt = MEAN%Mrt + Bwt(i) *  Srr(i) * S(i)
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_50 (NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts)
!!
!! Coded by: S W Key, 10-MAY-1991 18:36:45
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:29:16.00
!!
!! FINITE STRAIN UNIAXIAL-STRESS ELASTIC BEHAVIOR.
!!
!! Purpose: Compute the current fiber stresses, BEAM(NEL)%Axial(1:Npts)
!! and BEAM(NEL)%Shear(1:Npts). This approach to elasicity is based on
!! hypoelastic concepts. The theory and further references to it may be
!! found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol 32, pp 217-232, (1979).
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE beam_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN)    :: NEL
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(IN)    :: Drr(16),Drc(16),Ers,Ert
      INTEGER,         INTENT(IN)    :: Npts
!!
!! Retrieve elastic constants for beam material model.
!! YM = Young's modulus; SM = Shear modulus.
!!
      MatID = BEAM(NEL)%PAR%MatID
      YM = MATERIAL(MatID)%PVAL(6)
      SM = MATERIAL(MatID)%PVAL(9)
!!
!! Define (rho c**2) for longitudinal (RCL2) and shear waves (RCS2).
!!
      SOUND_SPEED%RCL2 = YM
      SOUND_SPEED%RCS2 = SM
!!
!! Integrate elastic constitutive model.
!!
      G2 = SM + SM
      DO i = 1,Npts
        BEAM(NEL)%RES%Axial(i) =                                               &
     &    BEAM(NEL)%RES%Axial(i) + (YM * BEAM(NEL)%RES%DTnext) * Drr(i)
        BEAM(NEL)%RES%Shear(i) =                                               &
     &    BEAM(NEL)%RES%Shear(i) + (G2 * BEAM(NEL)%RES%DTnext) * Drc(i)
      ENDDO
!!
!! Transverse shear stresses.
!!
      BEAM(NEL)%RES%Trs =                                                      &
     &          BEAM(NEL)%RES%Trs + (G2 * BEAM(NEL)%RES%DTnext) * Ers
      BEAM(NEL)%RES%Trt =                                                      &
     &          BEAM(NEL)%RES%Trt + (G2 * BEAM(NEL)%RES%DTnext) * Ert
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_51 (NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts)
!!
!! Coded by: S W Key, 10-MAY-1991 18:39:08
!! Copyright (c) by KEY Associates; 25-SEP-1993 21:42:53.27
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:29:34.00
!!
!! FINITE STRAIN UNIAXIAL-STRESS ELASTIC-PLASTIC BEHAVIOR.
!!
!! Purpose: Compute the current value of the fiber stresses, BEAM(NEL)%Axial(1:N
!! and BEAM(NEL)%Shear(1:Npts). This approach to plasticity is based on hypoelas
!! concepts. The theory and further references to it may be found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol 32, pp 217-232, (1979).
!!
!! The implicit radial return algorithm comes from:
!!
!!      J. C. Simo and S. Govindjee, Exact Closed-Form Solution
!!      of the Return Mapping Algorithm in Plane Stress Elasto-
!!      Viscoplasticity, ENGINEERING COMPUTATIONS, Vol 5, No 5,
!!      pp 254-258 (1988).
!!
!! A similar and somewhat earlier development by Ph. Jetteur can be
!! found in
!!
!!      Philippe Jetteur, Implicit Integration Algorithm for Elasto-
!!      Plasticity in Plane Stress, ENGINEERING COMPUTATIONS, Vol 3,
!!      No 3, pp 251-253 (1986).
!!
!!      STATE_VARIABLES(1:16) = Effective plastic strain
!!
      USE shared_common_data
      USE beam_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN)    :: NEL
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(IN)    :: Drr(16),Drc(16),Ers,Ert
      INTEGER,         INTENT(IN)    :: Npts
!!
!! Local variables.
!!
!!      DATA    Root2   /1.414213562/           ! Root2 = SQRT (2.0)
!!      DATA    Root6   /2.449489743/           ! Root6 = SQRT (6.0)
!!      DATA    Third   /0.33333333333333/      ! Third = 1.0 / 3.0
!!      DATA    Troot   /0.81649658092773/      ! Troot = SQRT (2.0/3.0)
!!
      REAL(KIND(0D0)), PARAMETER :: Third   = 1.0D+0/3.0D+0
      REAL(KIND(0D0)), PARAMETER :: Epsilon = 1.0D-6
      INTEGER,         PARAMETER :: KMAX    = 8
!!
      Troot = SQRT (2.0D+0 / 3.0D+0)
      Root2 = SQRT (2.0D+0)                      ! Root2 = SQRT (2.0)
      Root6 = SQRT (6.0D+0)                      ! Root6 = SQRT (6.0)
!!
!! Retrieve elastic constants for beam material model.
!! YM = Young's modulus; SM = Shear modulus.
!!
      MatID = BEAM(NEL)%PAR%MatID
      YM = MATERIAL(MatID)%PVAL(6)
      SM = MATERIAL(MatID)%PVAL(9)
!!
!! Define (rho c**2) for longitudinal (RCL2) and shear waves (RCS2).
!!
      SOUND_SPEED%RCL2 = YM
      SOUND_SPEED%RCS2 = SM
!!
!! Integrate elastic constitutive model to obtain trial elastic stress.
!!
      G2 = SM + SM
      DO i = 1,Npts
        BEAM(NEL)%RES%Axial(i) =                                               &
     &    BEAM(NEL)%RES%Axial(i) + (YM * BEAM(NEL)%RES%DTnext) * Drr(i)
        BEAM(NEL)%RES%Shear(i) =                                               &
     &    BEAM(NEL)%RES%Shear(i) + (G2 * BEAM(NEL)%RES%DTnext) * Drc(i)
      ENDDO
!!
!! Retrieve plastic material properties.
!!
      QR = MATERIAL(MatID)%PVAL( 7)   ! Poisson's ratio, nu
      QS = MATERIAL(MatID)%PVAL(10)   ! Yield stress
      QH = MATERIAL(MatID)%PVAL(12)   ! Plastic hardening modulus, H
      QB = MATERIAL(MatID)%PVAL(13)   ! Kinematic/isotrpic hard'g ratio, Beta
      QP = MATERIAL(MatID)%PVAL(16)   ! Strain rate parameter, p
      QD = MATERIAL(MatID)%PVAL(17)   ! Strain rate parameter, D
!!
!! Force beta to be 1.0 (isotropic hardening) until such time as kinematic
!! hardening has been developed.
!!
      QB = ONE
!!
!! Loop over all integration points in the beam cross section.
!!
      Isv = 1
      Nsv = MATERIAL(MatID)%Nsv
      DO 200 i = 1,Npts
        Axial = BEAM(NEL)%RES%Axial(i)
        Shear = BEAM(NEL)%RES%Shear(i)
        Eplas = STATE_VARIABLES(Isv)
!!
!! AK = Radius of yield surface x Troot.
!!
        AK = (QS + QB*QH*Eplas) * Troot
!!
!! Test for zero yield stress, (Elastic response).
!!
        IF (AK .EQ. 0.0) GO TO 200
!!
!! Compute magnitude of the deviatoric stress for use in yield calculation.
!!
        Psi = 2.0*(Third*(Axial*Axial) + (Shear*Shear))
!!
!! Check the yield condition.
!!
        IF ((Psi - AK*AK) .LE. 0.0) GO TO 200
!!
!! Compute dynamic yield stress.
!!
        IF (QP .NE. 0.0) THEN
          ED = 4.0*(Shear*Drc(i)) + 2.0*Third*(Axial*Drr(i))
          IF (ED .GT. 0.0) THEN
            ED  = ED/(SQRT(Psi)*(ONE+QH/(3.0D+0*SM)))
            ADD = QS*(ED/QD)**(ONE/QP)
            AK  = AK + ADD*Troot
            IF ((Psi - AK*AK) .LE. 0.0) GO TO 200
          ENDIF
        ENDIF
!!
!! Plastic loading condition. Calculate discrete consistency parameter dL.
!! (Assumes a constant plastic modulus H and constant total strain rate.)
!!
        dL = 0.0
        KOUNT = 0
        A0 = (2.0-QB-QB)*QH
        A1 = Axial*Axial
        A2 = Root6 * (G2*(ONE+QR)/(ONE-QR) + A0) * Third
        A3 = A1 + 4.0*(Shear*Shear)
        A4 = Root2 * (G2 + A0*Third)
        A5 = Troot*QB*QH
 100      CONTINUE
        B2 = Root6 + A2*dL
        C2 = B2 * B2
        B4 = Root2 + A4*dL
        C4 = B4 * B4
        F1 = A1/C2 + A3/C4
        A6 = SQRT (F1)
        A7 = AK + A5*(dL*A6)
        F2 = A7*A7
        IF ((F1-F2) .GT. Epsilon*F2) THEN
          KOUNT = KOUNT + 1
          DF1 = -2.0 * ((A1/(B2*C2))*A2 + (A3/(B4*C4))*A4)
          DF2 =  2.0 * A7*(A5*A6 + dL*A5*DF1/A6)
          dL = dL - (F1-F2)/(DF1-DF2)
          IF (KOUNT .LE. KMAX) GO TO 100
          WRITE (MSG1,'(I8)') BEAM(NEL)%PAR%EleID
          WRITE (MSG2,'(I8)') KMAX
          CALL USER_MESSAGE                                                    &
     &      (                                                                  &
     &      MSGL//'WARN'//                                                     &
     &      MSGL//'MATERIAL_51; Beam (BEAM) Element ID:'//MSG1//               &
     &      MSGL//'Iteration For Consistency Parameter dL'//                   &
     &      MSGL//'Requires More Than'//MSG2//' Iterations To Converge.'       &
     &      )
        ENDIF
        D2 = Root6 / B2
        D4 = Root2 / B4
        BEAM(NEL)%RES%Axial(i) = Axial * 0.5 * (D2 + D4)
        BEAM(NEL)%RES%Shear(i) = Shear * D4
        Eplas = Eplas + dL*A6
!!
!! End of loop on each integration point in the beam cross section.
!!
        STATE_VARIABLES(Isv) = Eplas
        Isv = Isv + Nsv
 200    ENDDO
!!
!! Transverse shear stresses.
!!
      BEAM(NEL)%RES%Trs =                                                      &
     &          BEAM(NEL)%RES%Trs + (G2 * BEAM(NEL)%RES%DTnext) * Ers
      BEAM(NEL)%RES%Trt =                                                      &
     &          BEAM(NEL)%RES%Trt + (G2 * BEAM(NEL)%RES%DTnext) * Ert
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_55 (NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts)
!!
!! Copyright (c) by KEY Associates;  3-AUG-1993 13:05:50.81
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:30:20.00
!!
!! FINITE STRAIN UNIAXIAL-STRESS ORTHOTROPIC ELASTIC BEHAVIOR.
!!
!! Purpose: Compute the current fiber stresses, BEAM(NEL)%Axial(1:Npts)
!! and BEAM(NEL)%Shear(1:Npts). This approach to elasicity is based on
!! hypoelastic concepts. The theory and further references to it may be
!! found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol 32, pp 217-232, (1979).
!!
!!      STATE_VARIABLES(-) = (not used)
!!
      USE shared_common_data
      USE beam_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN)    :: NEL
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(*)
      REAL(KIND(0D0)), INTENT(IN)    :: Drr(16),Drc(16),Ers,Ert
      INTEGER,         INTENT(IN)    :: Npts
!!
!! Retrieve elastic constants for orthotropic beam material model.
!! YM = Young's modulus in the a-direction (along the beam axis), Eaa;
!! SM = Shear modulus in the ac-/ab-plane (shear modulus for torsion),
!!    = (Gab + Gac)/2
!!
      MatID = BEAM(NEL)%PAR%MatID
      YM = MATERIAL(MatID)%PVAL(6)
      SM = 0.5 * (MATERIAL(MatID)%PVAL(12) + MATERIAL(MatID)%PVAL(13))
!!
!! Define (rho c**2) for longitudinal (RCL2) and shear waves (RCS2).
!!
      SOUND_SPEED%RCL2 = YM
      SOUND_SPEED%RCS2 = SM
!!
!! Integrate elastic constitutive model.
!!
      G2 = SM + SM
      DO i = 1,Npts
        BEAM(NEL)%RES%Axial(i) =                                               &
     &    BEAM(NEL)%RES%Axial(i) + (YM * BEAM(NEL)%RES%DTnext) * Drr(i)
        BEAM(NEL)%RES%Shear(i) =                                               &
     &    BEAM(NEL)%RES%Shear(i) + (G2 * BEAM(NEL)%RES%DTnext) * Drc(i)
      ENDDO
!!
!! Transverse shear stresses.
!!
      BEAM(NEL)%RES%Trs =                                                      &
     &          BEAM(NEL)%RES%Trs + (G2 * BEAM(NEL)%RES%DTnext) * Ers
      BEAM(NEL)%RES%Trt =                                                      &
     &          BEAM(NEL)%RES%Trt + (G2 * BEAM(NEL)%RES%DTnext) * Ert
!!
      RETURN
      END
!!_
      SUBROUTINE MATERIAL_57 (NEL,STATE_VARIABLES,Drr,Drc,Ers,Ert,Npts)
!!
!! Copyright (c) by KEY Associates;  3-AUG-1993 13:21:58.16
!! Copyright (c) by KEY Associates; 10-OCT-1997 09:30:42.00
!!
!! FINITE STRAIN UNIAXIAL-STRESS RUBBER ELASTIC BEHAVIOR.
!!
!! Purpose: Compute the current value of the fiber stresses, BEAM(NEL)%Axial(1:N
!! and BEAM(NEL)%Shear(1:Npts). This approach to elasicity is based on hypoelast
!! concepts. The theory and further references to it may be found in:
!!
!!      J.K. Dienes, On the Analysis of Rotation and Stress Rate in
!!      Deforming Bodies, ACTA MECHANICA, Vol 32, pp 217-232, (1979).
!!
!!      STATE_VARIABLES(1:16,1) = Err, rr-component, total strain
!!      STATE_VARIABLES(1:16,2) = Erc, rc-component, total strain
!!      STATE_VARIABLES(1:16,3) = E-bar, effective strain, 2(E'ijE'ij)/3
!!      STATE_VARIABLES(1:16,4) = S-bar, effective stress, 3(S'ijS'ij)/2
!!
      USE shared_common_data
      USE beam_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN)    :: NEL
      REAL(KIND(0D0)), INTENT(INOUT) :: STATE_VARIABLES(1:16,*)
      REAL(KIND(0D0)), INTENT(IN)    :: Drr(16),Drc(16),Ers,Ert
      INTEGER,         INTENT(IN)    :: Npts
!!
!! Local variables.
!!
      INTEGER :: HstID
      REAL(KIND(0D0)) :: TABLE_LOOK_UP
!!
      DATA                                                                     &
     &          OneThird /0.33333333333333333D+0/                               &
     &          TwoThird /0.66666666666666667D+0/
!!
!! Retrieve elastic constants for beam material model.
!! Bulk = Bulk modulus; HstID = Eff. stress vs. eff. strain tabulated ftn.
!!
      MatID = BEAM(NEL)%PAR%MatID
      Bulk = MATERIAL(MatID)%PVAL(10)
      HstID = NINT (MATERIAL(MatID)%PVAL(9))
      P2 = 0.0
!!
!! Loop over all beam cross section stress evaluation points.
!!
      DO i = 1,Npts
!!
!! Update strain.
!!
        Err = STATE_VARIABLES(i,1) + BEAM(NEL)%RES%DTnext * Drr(i)
        Erc = STATE_VARIABLES(i,2) + BEAM(NEL)%RES%DTnext * Drc(i)

        STATE_VARIABLES(i,1) = Err
        STATE_VARIABLES(i,2) = Erc
!!
!! Retrieve old effective stress and strain.
!!
        Old_Effective_Strain = STATE_VARIABLES(i,3)
        Old_Effective_Stress = STATE_VARIABLES(i,4)
!!
!! Compute new effective strain. (This strain has no volumetric component
!! and therefore is only a deviatoric strain.)
!!
        Effective_Strain =                                                     &
     &          SQRT (Err*Err + TwoThird * ((Erc*Erc) + (Erc*Erc)))
!!
!! Look up new effective stress.
!!
        Effective_Stress =                                                     &
     &          TABLE_LOOK_UP ( HstID,Effective_Strain )
!!
!! Update effective stress and strain.
!!
        STATE_VARIABLES(i,3) = Effective_Strain
        STATE_VARIABLES(i,4) = Effective_Stress
!!
!! Compute new elastic stress.
!!
        IF (Effective_Strain .GT. 1.0D-20) THEN
          Alpha = TwoThird * Effective_Stress / Effective_Strain
        ELSE
          Alpha = 0.0
        ENDIF

        BEAM(NEL)%RES%Axial(i) = Alpha * Err * 1.5D0
        BEAM(NEL)%RES%Shear(i) = Alpha * Erc
!!
!! Compute tangent shear modulus.
!!
        IF (ABS(Effective_Strain-Old_Effective_Strain).GT.1.0D-20) THEN
          Q2 = OneThird * (Effective_Stress - Old_Effective_Stress)            &
     &                    / (Effective_Strain - Old_Effective_Strain)
        ELSE
          Q2 = 1.0D-3 * Bulk
        ENDIF
        P2 = MAX (P2, Q2)

      ENDDO
!!
!! Transverse shear stresses.
!!
      G2 = P2+P2
      BEAM(NEL)%RES%Trs =                                                      &
     &          BEAM(NEL)%RES%Trs + (G2 * BEAM(NEL)%RES%DTnext) * Ers
      BEAM(NEL)%RES%Trt =                                                      &
     &          BEAM(NEL)%RES%Trt + (G2 * BEAM(NEL)%RES%DTnext) * Ert
!!
!! Sound speeds squared * RHO(T)
!!
      SOUND_SPEED%RCL2 = Bulk + TwoThird * G2
      SOUND_SPEED%RCS2 = P2
!!
      RETURN
      END
