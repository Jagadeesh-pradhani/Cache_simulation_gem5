      SUBROUTINE PLATQ_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
      USE stress_
      USE state_variables_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Local variables.
      INTEGER :: NEL,SecID,MatID
      INTEGER :: ERRORCOUNT
      LOGICAL :: FOUND
!!
      REAL(KIND(0D0)) :: Nr,Ns,Nt
!!
      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
      TIME_STEP_MIN = HUGE(TIME_STEP_MIN)
      TIME_STEP_MAX = 0.0D0
	ERRORCOUNT=ERROR%COUNT
!!
!$OMP PARALLEL DO &
!$OMP   DEFAULT(PRIVATE), SHARED(PLATQ,MOTION,MATERIAL,STATE_VARIABLES), &
!$OMP   SHARED(CONTROL,TIMSIM,NODE,SECTION_2D,TABULATED_FUNCTION,STRESS),&
!$OMP   SHARED(NUMP4) REDUCTION(+:ERRORCOUNT),                          &
!$OMP   REDUCTION(MIN:TIME_STEP_MIN),                                   &
!$OMP   REDUCTION(MAX:TIME_STEP_MAX)
!!
      DO N = 1,NUMP4
!!
!! Gather element motion.
!!
        DO i = 1,8
          IX = PLATQ(N)%PAR%IX(i)
          Px(i) = MOTION(IX)%Px
          Py(i) = MOTION(IX)%Py
          Pz(i) = MOTION(IX)%Pz
          Ux(i) = MOTION(IX)%Ux
          Uy(i) = MOTION(IX)%Uy
          Uz(i) = MOTION(IX)%Uz
          Vx(i) = MOTION(IX)%Vx
          Vy(i) = MOTION(IX)%Vy
          Vz(i) = MOTION(IX)%Vz
        ENDDO
!!
!! Initialize element clock and time step
!!
        PLATQ(N)%RES%Time   = 0.0
        PLATQ(N)%RES%DTnext = 0.0
!!
!! Initialize hourglass stiffness forces.
!!
        PLATQ(N)%RES%Pr(1) = 0.0
        PLATQ(N)%RES%Pr(2) = 0.0
        PLATQ(N)%RES%Ps(1) = 0.0
        PLATQ(N)%RES%Ps(2) = 0.0
        PLATQ(N)%RES%Pt(1) = 0.0
        PLATQ(N)%RES%Pt(2) = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
        MatID = PLATQ(N)%PAR%MatID
        SecID = PLATQ(N)%PAR%SecID
!!
!! Gradient operator, stretching, and rotation. (For rigid bodies
!! only the element volume calculation is required.)
!!
        IF (CONTROL%BTPLTQ .EQ. 0) THEN
          CALL KHPLQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
        ELSE
          CALL BTPLQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
        ENDIF
!!
!! Save initial element area for later area strain calculations.
!!
        PLATQ(N)%PAR%Area = PLATQ(N)%RES%Area

        IF (PLATQ(N)%RES%Area .LE. 0.0) THEN
          WRITE (MSG1,'(I8)') PLATQ(N)%PAR%EleID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PLATQ_INITIALIZATION.001.00'//                          &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Has A Zero Or Negative Area.'                           &
     &          )
          ERRORCOUNT = ERRORCOUNT + 1
        ENDIF
!!
!! Compute element mass matrix and store in global mass matrix.
!!
        CALL PLATQ_MASS ( NEL,SecID,MatID )
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (PLATQ(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present. ELsize is an estimate of the distance between the
!! first nodal point and it's closest neighbor in the element using taxicab
!! geometry.
!!
          ELsize = HUGE ( ELsize )
          DO i = 2,4
            Qdist = ABS (Px(1) - Px(i)) + ABS (Py(1) - Py(i)) + ABS (Pz(1) - Pz(i))
            ELsize = MIN (ELsize,Qdist)
          ENDDO
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          PLATQ(N)%RES%DTelt = Elsize / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Convert global material direction vectors to local element coordinates.
!!
          IF (MATERIAL(MatID)%Type .EQ. 45) THEN
            Isv = PLATQ(N)%PAR%Isv
            Nsv = MATERIAL(MatID)%Nsv
            DO i = 1,Ipts_PLATQ (NEL)
              CALL GLOBAL_TO_LOCAL                                             &
     &          (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))
              Isv = Isv + Nsv
            ENDDO
          ENDIF
!!
!! Find initial sound speeds and stress.
!!
          CALL PLATQ_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Compute initial stress divergence and time step.
!!
          IF (CONTROL%BTPLTQ .EQ. 0) THEN
            CALL KHPLQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
          ELSE
            CALL BTPLQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
          ENDIF
!!
!! End of rigid body element if-test.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
!SPEC_omp2001        TIMSIM%DTP4x = MAX (TIMSIM%DTP4x,PLATQ(N)%RES%DTelt)
!SPEC_omp2001        i = 0
!SPEC_omp2001        FOUND = .FALSE.
!SPEC_omp2001        DO WHILE (.NOT.FOUND .and. i.LT.NPNDT)
!SPEC_omp2001          i = i + 1
!SPEC_omp2001          FOUND = PLATQ(N)%RES%DTelt .LT. TIMSIM%DTPl4(i)
!SPEC_omp2001        ENDDO
!SPEC_omp2001        IF (FOUND) THEN
!SPEC_omp2001          TIMSIM%DTPl4(i:NPNDT) = (/PLATQ(N)%RES%DTelt,TIMSIM%DTPl4(i:NPNDT-1)/)
!SPEC_omp2001          TIMSIM%Plat4(i:NPNDT) = (/N                 ,TIMSIM%Plat4(i:NPNDT-1)/)
!SPEC_omp2001        ENDIF
!!
!! End of element do-loop.
!!
        TIME_STEP_MIN = MIN (TIME_STEP_MIN, PLATQ(N)%RES%DTelt)
        TIME_STEP_MAX = MAX (TIME_STEP_MAX, PLATQ(N)%RES%DTelt)
      ENDDO
!!
!$OMP END PARALLEL DO
!!
      ERROR%COUNT = ERRORCOUNT
!!
!! Write out element time step values for compare using diff.
!!
!!!      do n = 1,nump4
!!!        write (1,*) n,PLATQ(n)%RES%DTelt
!!!      enddo
!!
      TIMSIM%DTP4x    = TIME_STEP_MAX
      TIMSIM%DTPl4(1) = TIME_STEP_MIN
      TIMSIM%Plat4(1) = 1 ! (Any element number will do.)
!!
!!!      write (3,*) TIME_STEP_MIN, TIME_STEP_MAX
!!
!END SPEC_omp2001
!!
      RETURN
      END
!!_
      SUBROUTINE PLATQ_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
      USE stress_
      USE state_variables_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Local variables.
      INTEGER :: NEL,SecID,MatID
!!
      REAL(KIND(0D0)) :: Nr,Ns,Nt
!!
      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
      TIME_STEP_MIN = HUGE(TIME_STEP_MIN)
      TIME_STEP_MAX = 0.0D0
!!
!$OMP PARALLEL DO                                                        &
!$OMP   DEFAULT(PRIVATE), SHARED(PLATQ,MOTION,MATERIAL,STATE_VARIABLES), &
!$OMP   SHARED(CONTROL,TIMSIM,NODE,SECTION_2D,TABULATED_FUNCTION,STRESS),&
!$OMP   SHARED(NUMP4,COUNTER),                                   &
!$OMP   REDUCTION(MIN:TIME_STEP_MIN),                                   &
!$OMP   REDUCTION(MAX:TIME_STEP_MAX)
!!
      DO N = 1,NUMP4
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,PLATQ(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (PLATQ(N)%PAR%ParID .GE. 0) THEN
!!
!! Gather element motion.
!!
            DO i = 1,8
              IX = PLATQ(N)%PAR%IX(i)
              Px(i) = MOTION(IX)%Px
              Py(i) = MOTION(IX)%Py
              Pz(i) = MOTION(IX)%Pz
              Ux(i) = MOTION(IX)%Ux
              Uy(i) = MOTION(IX)%Uy
              Uz(i) = MOTION(IX)%Uz
              Vx(i) = MOTION(IX)%Vx
              Vy(i) = MOTION(IX)%Vy
              Vz(i) = MOTION(IX)%Vz
            ENDDO
!!
!! Count element execution.
!!
!SPEC_omp2001            COUNTER%PLATQ = COUNTER%PLATQ + 1
!!
!! Increment element clock.
!!
            PLATQ(N)%RES%Time = PLATQ(N)%RES%Time + PLATQ(N)%RES%DTnext
!!
!! Scale nodal positions to current element time. Note that rotational
!! degrees of freedom do not require any scaling since they are not used.
!!
            DO i = 1,4
              QA = NODE(PLATQ(N)%PAR%IX(i))%Time - PLATQ(N)%RES%Time
              Ux(i) = Ux(i) - QA * Vx(i)
              Uy(i) = Uy(i) - QA * Vy(i)
              Uz(i) = Uz(i) - QA * Vz(i)
            ENDDO
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Retrieve state variable pointer and material ID from element data structure.
!!
            MatID = PLATQ(N)%PAR%MatID
            SecID = PLATQ(N)%PAR%SecID
!!
!! Distinguish between the Key-Hoff and Belytschko-Tsay shell formulations.
!!
            IF (CONTROL%BTPLTQ .EQ. 0) THEN
!!
!! Gradient operator, stretching, and rotation.
!!
              CALL KHPLQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Find initial sound speeds and stress.
!!
              CALL PLATQ_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Divergence operator (at the end of the interval).
!!
              CALL KHPLQ_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Compute initial stress divergence and time step.
!!
              CALL KHPLQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
            ELSE
!!
!! Gradient operator, stretching, and rotation.
!!
              CALL BTPLQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Find initial sound speeds and stress.
!!
              CALL PLATQ_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
              IF (MATERIAL(MatID)%PVAL(5) .GT. 0.0) THEN
                CALL PLATQ_HOURGLASS_FORCES ( NEL,SecID,MatID )
              ENDIF
!!
!! Divergence operator (at the end of the interval).
!!
              CALL BTPLQ_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Compute initial stress divergence and time step.
!!
              CALL BTPLQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
            ENDIF
!!
!! End of rigid body element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
!SPEC_omp2001         TIMSIM%DTP4x = MAX (TIMSIM%DTP4x,PLATQ(N)%RES%DTelt)
!SPEC_omp2001         IF (PLATQ(N)%RES%DTelt .LT. TIMSIM%DTPl4(1)) THEN
!SPEC_omp2001           TIMSIM%DTPl4(1) = PLATQ(N)%RES%DTelt
!SPEC_omp2001           TIMSIM%Plat4(1) = N
!SPEC_omp2001         ENDIF
!!
!! End of time-to-subcycle if-test.
!!
        ENDIF
!!
!! End of element do-loop.
!!
        TIME_STEP_MIN = MIN (TIME_STEP_MIN, PLATQ(N)%RES%DTelt)
        TIME_STEP_MAX = MAX (TIME_STEP_MAX, PLATQ(N)%RES%DTelt)
      ENDDO
!!
!$OMP END PARALLEL DO
!!
!! Write out element time step values for compare using diff.
!!
!!!      do n = 1,nump4
!!!        write (2,*) TIMSIM%Step,n,PLATQ(n)%RES%DTelt
!!!      enddo
!!
!BEGIN SPEC_omp2001
!!
      TIMSIM%DTP4x    = TIME_STEP_MAX
      TIMSIM%DTPl4(1) = TIME_STEP_MIN
      TIMSIM%Plat4(1) = 1 ! (Any element number will do.)
!!
!!!      write (4,*) TIMSIM%Step,TIME_STEP_MIN,TIME_STEP_MAX
!!!      if (TIMSIM%Step .ge. 10) stop ' You now have time step arrays.'
!!
!END SPEC_omp2001
!!
      RETURN
      END
!!_
      SUBROUTINE PLATQ_MASS ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 15:28:02
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one fourth at each node and an isotropic inertia.
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
      REAL(KIND(0D0)) :: Nr,Ns,Nt
!!
      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Compute one fourth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = Density*PLATQ(NEL)%PAR%Area*SECTION_2D(SecID)%Thickness/4.0D+0
!!
!! Accumulate mass at each nodal point.
!!
!!
!$OMP CRITICAL (PLATQ_MASS_VALUES)
      DO i = 1,4
        NODE(PLATQ(NEL)%PAR%IX(i))%Mass = NODE(PLATQ(NEL)%PAR%IX(i))%Mass + QMass
        MATERIAL(MatID)%Mass = MATERIAL(MatID)%Mass + QMass
        MATERIAL(MatID)%Xcm  = MATERIAL(MatID)%Xcm  + QMass * Px(I)
        MATERIAL(MatID)%Ycm  = MATERIAL(MatID)%Ycm  + QMass * Py(I)
        MATERIAL(MatID)%Zcm  = MATERIAL(MatID)%Zcm  + QMass * Pz(I)
!!
!! Compute inertia tensor B wrt the origin from nodal point masses.
!!
        MATERIAL(MatID)%Bxx = MATERIAL(MatID)%Bxx + (Py(I)*Py(I)+Pz(I)*Pz(I))*QMass
        MATERIAL(MatID)%Byy = MATERIAL(MatID)%Byy + (Px(I)*Px(I)+Pz(I)*Pz(I))*QMass
        MATERIAL(MatID)%Bzz = MATERIAL(MatID)%Bzz + (Px(I)*Px(I)+Py(I)*Py(I))*QMass
        MATERIAL(MatID)%Bxy = MATERIAL(MatID)%Bxy - Px(I)*Py(I)*QMass
        MATERIAL(MatID)%Bxz = MATERIAL(MatID)%Bxz - Px(I)*Pz(I)*QMass
        MATERIAL(MatID)%Byz = MATERIAL(MatID)%Byz - Py(I)*Pz(I)*QMass
      ENDDO
!!
!!
!! Compute nodal isotropic inertia
!!
      RMass = QMass * (PLATQ(NEL)%PAR%Area + SECTION_2D(SecID)%Thickness**2) / 12.0D+0
!!
!!
      NODE(PLATQ(NEL)%PAR%IX(5))%Mass = NODE(PLATQ(NEL)%PAR%IX(5))%Mass + RMass
      NODE(PLATQ(NEL)%PAR%IX(6))%Mass = NODE(PLATQ(NEL)%PAR%IX(6))%Mass + RMass
      NODE(PLATQ(NEL)%PAR%IX(7))%Mass = NODE(PLATQ(NEL)%PAR%IX(7))%Mass + RMass
      NODE(PLATQ(NEL)%PAR%IX(8))%Mass = NODE(PLATQ(NEL)%PAR%IX(8))%Mass + RMass
!$OMP END CRITICAL (PLATQ_MASS_VALUES)
!!
!!
      RETURN
      END
!!_
      SUBROUTINE KHPLQ_GRADIENT_OPERATOR ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) ::                                                       &
     &                            Wx(4),Wy(4),Wz(4),X(4),Y(4),Z(4),            &
     &          Vr(4),Vs(4),Vt(4),Wr(4),Ws(4),Wt(4),R(4),S(4),T(4),            &
     &          Theta_r(4),Theta_s(4),Theta_t(4),Pr(4),Ps(4),Pt(4)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Mid-interval position of nodal points, and current translational and
!! rotational velocities.
!!
      DT = 0.5D+0 * PLATQ(NEL)%RES%DTnext
      DO i = 1,4
        X(i)  = Px(i) + (Ux(i) - DT * Vx(i))
        Y(i)  = Py(i) + (Uy(i) - DT * Vy(i))
        Z(i)  = Pz(i) + (Uz(i) - DT * Vz(i))
        Wx(i) = Vx(i+4)
        Wy(i) = Vy(i+4)
        Wz(i) = Vz(i+4)
      ENDDO
!!
!! Construct local basis vectors at the center of the element.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag = 0.5D+0 * SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (0.5D+0 / Rmag)
      Ry = Ry * (0.5D+0 / Rmag)
      Rz = Rz * (0.5D+0 / Rmag)
      Smag = 0.5D+0 * SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (0.5D+0 / Smag)
      Sy = Sy * (0.5D+0 / Smag)
      Sz = Sz * (0.5D+0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATQ_GRADIENT_OPERATOR.002.00'//                       &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and rotational
!! velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
        Wr(i) = Rx*Wx(i) + Ry*Wy(i) + Rz*Wz(i)
        Ws(i) = Sx*Wx(i) + Sy*Wy(i) + Sz*Wz(i)
        Wt(i) = Tx*Wx(i) + Ty*Wy(i) + Tz*Wz(i)
      ENDDO
!!
!! Define a set of element edge vectors P (Pr,Ps,Pt) for use in constructing
!! a unit normal at each corner node.
!!
      i = 4
      j = 1
      DO k = 1,4
        Pr(i) = R(j) - R(i)
        Ps(i) = S(j) - S(i)
        Pt(i) = T(j) - T(i)
        j = i
        i = i - 1
      ENDDO
!!
!! Define the unit vector N (Nr,Ns,Nt) normal to the element at each corner
!! node. N = (P- x P+)/||P- x P+||
!!
      j = 4
      DO i = 1,4
        Qr = Ps(j)*Pt(i) - Ps(i)*Pt(j)
        Qs = Pt(j)*Pr(i) - Pt(i)*Pr(j)
        Qt = Pr(j)*Ps(i) - Pr(i)*Ps(j)
        Qmag  = ONE / SQRT (Qr*Qr + Qs*Qs + Qt*Qt)
        Nr(i) = Qr * Qmag
        Ns(i) = Qs * Qmag
        Nt(i) = Qt * Qmag
        j = i
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R14 = R(1) - R(4)
      R23 = R(2) - R(3)
      R24 = R(2) - R(4)
      R34 = R(3) - R(4)
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S14 = S(1) - S(4)
      S23 = S(2) - S(3)
      S24 = S(2) - S(4)
      S34 = S(3) - S(4)
!!
!! A0t(i); mean averaging operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.25000000000000000D-01 * (R13*S24 - S13*R24)
      QB = 4.16666666666666667D-02 * (R34*S12 - S34*R12)
      QC = 4.16666666666666667D-02 * (R23*S14 - S23*R14)
      A0t(1) = QA - (QB + QC)
      A0t(2) = QA + (QB - QC)
      A0t(3) = QA + (QB + QC)
      A0t(4) = QA - (QB - QC)
!!
!! A1t(i) and A2t(i); linear averaging operators (area weighting based on
!! bilinear shape functions times Xi1 and Xi2, respectively).
!!
      QD = 4.16666666666666667D-02 * (R23*S12 - S23*R12)
      QE = 4.16666666666666667D-02 * (R12*S24 - S12*R24)
      QF = 4.16666666666666667D-02 * (R13*S34 - S13*R34)
      QG = 4.16666666666666667D-02 * (R34*S23 - S34*R23)
      QH = 0.333333333333333333D+0 * QC
      QI = 0.333333333333333333D+0 * QB
      A1t(1) = ( QD - QH)
      A1t(2) = ( QE + QH)
      A1t(3) = ( QF - QH)
      A1t(4) = ( QG + QH)
      A2t(1) = (-QF - QI)
      A2t(2) = ( QG + QI)
      A2t(3) = (-QD - QI)
      A2t(4) = ( QE + QI)
!!
!! B0r(i) and B0s(i); mean membrane gradient operators.
!!
      B0r(1) =  (0.5D+0 * S24)
      B0r(3) = -(0.5D+0 * S24)
      B0r(4) =  (0.5D+0 * S13)
      B0r(2) = -(0.5D+0 * S13)
      B0s(3) =  (0.5D+0 * R24)
      B0s(1) = -(0.5D+0 * R24)
      B0s(2) =  (0.5D+0 * R13)
      B0s(4) = -(0.5D+0 * R13)
!!
!! B1s(i); Xi1-linear membrane gradient operators.
!!
      B1s(1) =  (8.3333333333333333D-02 * R34)
      B1s(2) = -(8.3333333333333333D-02 * R34)
      B1s(3) = -(8.3333333333333333D-02 * R12)
      B1s(4) =  (8.3333333333333333D-02 * R12)
!!
!! B2r(i); Xi2-linear membrane gradient operators.
!!
      B2r(1) = -(8.3333333333333333D-02 * S23)
      B2r(4) =  (8.3333333333333333D-02 * S23)
      B2r(2) =  (8.3333333333333333D-02 * S14)
      B2r(3) = -(8.3333333333333333D-02 * S14)
!!
      X12 = Nr(1) - Nr(2)
      X13 = Nr(1) - Nr(3)
      X14 = Nr(1) - Nr(4)
      X23 = Nr(2) - Nr(3)
      X24 = Nr(2) - Nr(4)
      X34 = Nr(3) - Nr(4)
      Y12 = Ns(1) - Ns(2)
      Y13 = Ns(1) - Ns(3)
      Y14 = Ns(1) - Ns(4)
      Y23 = Ns(2) - Ns(3)
      Y24 = Ns(2) - Ns(4)
      Y34 = Ns(3) - Ns(4)
!!
!! C0r(i) & C0s(i); mean coupling operator (membrane-bending coupling)
!!
      C0r(1) =  (0.5D+0 * Y24)
      C0r(3) = -(0.5D+0 * Y24)
      C0r(4) =  (0.5D+0 * Y13)
      C0r(2) = -(0.5D+0 * Y13)
      C0s(3) =  (0.5D+0 * X24)
      C0s(1) = -(0.5D+0 * X24)
      C0s(2) =  (0.5D+0 * X13)
      C0s(4) = -(0.5D+0 * X13)
!!
!! C0t(i); mean coupling operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.2500000000000000D-1 * (X13*S24 - S13*X24)
      QB = 4.1666666666666667D-2 * (X34*S12 - S34*X12)
      QC = 4.1666666666666667D-2 * (X23*S14 - S23*X14)
      C0t(1) = QA - (QB + QC)
      C0t(2) = QA + (QB - QC)
      C0t(3) = QA + (QB + QC)
      C0t(4) = QA - (QB - QC)
!!
      QA = 1.2500000000000000D-1 * (R13*Y24 - Y13*R24)
      QB = 4.1666666666666667D-2 * (R34*Y12 - Y34*R12)
      QC = 4.1666666666666667D-2 * (R23*Y14 - Y23*R14)
      C0t(1) = C0t(1) + QA - (QB + QC)
      C0t(2) = C0t(2) + QA + (QB - QC)
      C0t(3) = C0t(3) + QA + (QB + QC)
      C0t(4) = C0t(4) + QA - (QB - QC)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATQ(NEL)%RES%Area = 0.5D+0 * (R13*S24 - S13*R24)
      Ain = ONE / PLATQ(NEL)%RES%Area
      Hin = 0.5D+0 * Ain
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + B0r(i)*B0r(i) + B0s(i)*B0s(i)
      ENDDO
      Delta = Ain * SQRT (Dx+Dx)
!!
!! Convert global rotation rates to local theta's.
!!
      DO i = 1,4
        Theta_r(i) = -Ns(i)*Wt(i) + Nt(i)*Ws(i)
        Theta_s(i) = -Nt(i)*Wr(i) + Nr(i)*Wt(i)
        Theta_t(i) = -Nr(i)*Ws(i) + Ns(i)*Wr(i)
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
      DO i = 1,4
        Vrr = Vrr +  Vr(i)*B0r(i)
        Vsr = Vsr +  Vs(i)*B0r(i)
        Vtr = Vtr +  Vt(i)*B0r(i)
        Vrs = Vrs +  Vr(i)*B0s(i)
        Vss = Vss +  Vs(i)*B0s(i)
        Vts = Vts +  Vt(i)*B0s(i)
        Vrt = Vrt +  Theta_r(i)*A0t(i)
        Vst = Vst +  Theta_s(i)*A0t(i)
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
!! Xi1-linear transverse-shear stretching components.
!!
      Vts = 0.0
      Vst = 0.0
      DO i = 1,4
        Vts = Vts +  Vt(i)*B1s(i)
        Vst = Vst +  Theta_s(i)*A1t(i)
      ENDDO
!!
      Est = Hin * (Vst + Vts)
!!
!! Xi2-linear transverse-shear stretching components.
!!
      Vtr = 0.0
      Vrt = 0.0
      DO i = 1,4
        Vtr = Vtr +  Vt(i)*B2r(i)
        Vrt = Vrt +  Theta_r(i)*A2t(i)
      ENDDO
!!
      Frt = Hin * (Vrt + Vtr)
!!
!! II. BENDING: Mean bending velocity gradients.
!!
      Vrr = 0.0
      Vsr = 0.0
      Vtr = 0.0
      Vrs = 0.0
      Vss = 0.0
      Vts = 0.0
      Vrt = 0.0
      Vst = 0.0
      DO i = 1,4
        Vrr = Vrr + Vr(i)*C0r(i) + Theta_r(i)*B0r(i)
        Vsr = Vsr + Vs(i)*C0r(i) + Theta_s(i)*B0r(i)
        Vtr = Vtr + Vt(i)*C0r(i) + Theta_t(i)*B0r(i)
        Vrs = Vrs + Vr(i)*C0s(i) + Theta_r(i)*B0s(i)
        Vss = Vss + Vs(i)*C0s(i) + Theta_s(i)*B0s(i)
        Vts = Vts + Vt(i)*C0s(i) + Theta_t(i)*B0s(i)
        Vrt = Vrt + Theta_r(i)*C0t(i)
        Vst = Vst + Theta_s(i)*C0t(i)
      ENDDO
!!
!! Mean bending stretching.
!!
      Hrr = Ain * Vrr
      Hss = Ain * Vss
      Hrs = Hin * (Vrs + Vsr)
      Urs = Hin * (Vrs - Vsr)
!!
!! Mean transverse-shear stretching components.
!!
      Hrt = Hin * (Vrt + Vtr)
      Hst = Hin * (Vst + Vts)
!!
!! III. HOURGLASS CONTROL: Anti-hourglass gradients.
!!
      RETURN
      END
!!_
      SUBROUTINE KHPLQ_DIVERGENCE_OPERATOR ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) ::                                                       &
     &          X(4),Y(4),Z(4),R(4),S(4),T(4),Theta_r(4),Theta_s(4),           &
     &          Theta_t(4),Pr(4),Ps(4),Pt(4)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Current position of nodal points, and current translational and
!! rotational velocities.
!!
      DO i = 1,4
        X(i)  = Px(i) + Ux(i)
        Y(i)  = Py(i) + Uy(i)
        Z(i)  = Pz(i) + Uz(i)
      ENDDO
!!
!! Construct local basis vectors at the center of the element.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag = 0.5D+0 * SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (0.5D+0 / Rmag)
      Ry = Ry * (0.5D+0 / Rmag)
      Rz = Rz * (0.5D+0 / Rmag)
      Smag = 0.5D+0 * SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (0.5D+0 / Smag)
      Sy = Sy * (0.5D+0 / Smag)
      Sz = Sz * (0.5D+0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATQ_DIVERGENCE_OPERATOR.002.00'//                     &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and rotational
!! velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
      ENDDO
!!
!! Define a set of element edge vectors P (Pr,Ps,Pt) for use in constructing
!! a unit normal at each corner node.
!!
      i = 4
      j = 1
      DO k = 1,4
        Pr(i) = R(j) - R(i)
        Ps(i) = S(j) - S(i)
        Pt(i) = T(j) - T(i)
        j = i
        i = i - 1
      ENDDO
!!
!! Define the unit vector N (Nr,Ns,Nt) normal to the element at each corner
!! node. N = (P- x P+)/||P- x P+||
!!
      j = 4
      DO i = 1,4
        Qr = Ps(j)*Pt(i) - Ps(i)*Pt(j)
        Qs = Pt(j)*Pr(i) - Pt(i)*Pr(j)
        Qt = Pr(j)*Ps(i) - Pr(i)*Ps(j)
        Qmag  = ONE / SQRT (Qr*Qr + Qs*Qs + Qt*Qt)
        Nr(i) = Qr * Qmag
        Ns(i) = Qs * Qmag
        Nt(i) = Qt * Qmag
        j = i
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R14 = R(1) - R(4)
      R23 = R(2) - R(3)
      R24 = R(2) - R(4)
      R34 = R(3) - R(4)
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S14 = S(1) - S(4)
      S23 = S(2) - S(3)
      S24 = S(2) - S(4)
      S34 = S(3) - S(4)
!!
!! A0t(i); mean averaging operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.2500000000000000D-01 * (R13*S24 - S13*R24)
      QB = 4.1666666666666667D-02 * (R34*S12 - S34*R12)
      QC = 4.1666666666666667D-02 * (R23*S14 - S23*R14)
      A0t(1) = QA - (QB + QC)
      A0t(2) = QA + (QB - QC)
      A0t(3) = QA + (QB + QC)
      A0t(4) = QA - (QB - QC)
!!
!! A1t(i) and A2t(i); linear averaging operators (area weighting based on
!! bilinear shape functions times Xi1 and Xi2, respectively).
!!
      QD = 4.1666666666666667D-02 * (R23*S12 - S23*R12)
      QE = 4.1666666666666667D-02 * (R12*S24 - S12*R24)
      QF = 4.1666666666666667D-02 * (R13*S34 - S13*R34)
      QG = 4.1666666666666667D-02 * (R34*S23 - S34*R23)
      QH = 0.33333333333333333D+0 * QC
      QI = 0.33333333333333333D+0 * QB
      A1t(1) = ( QD - QH)
      A1t(2) = ( QE + QH)
      A1t(3) = ( QF - QH)
      A1t(4) = ( QG + QH)
      A2t(1) = (-QF - QI)
      A2t(2) = ( QG + QI)
      A2t(3) = (-QD - QI)
      A2t(4) = ( QE + QI)
!!
!! B0r(i) and B0s(i); mean membrane gradient operators.
!!
      B0r(1) =  (0.5D+0 * S24)
      B0r(3) = -(0.5D+0 * S24)
      B0r(4) =  (0.5D+0 * S13)
      B0r(2) = -(0.5D+0 * S13)
      B0s(3) =  (0.5D+0 * R24)
      B0s(1) = -(0.5D+0 * R24)
      B0s(2) =  (0.5D+0 * R13)
      B0s(4) = -(0.5D+0 * R13)
!!
!! B1s(i); Xi1-linear membrane gradient operators.
!!
      B1s(1) =  (8.333333333333333D-02 * R34)
      B1s(2) = -(8.333333333333333D-02 * R34)
      B1s(3) = -(8.333333333333333D-02 * R12)
      B1s(4) =  (8.333333333333333D-02 * R12)
!!
!! B2r(i); Xi2-linear membrane gradient operators.
!!
      B2r(1) = -(8.333333333333333D-02 * S23)
      B2r(4) =  (8.333333333333333D-02 * S23)
      B2r(2) =  (8.333333333333333D-02 * S14)
      B2r(3) = -(8.333333333333333D-02 * S14)
!!
      X12 = Nr(1) - Nr(2)
      X13 = Nr(1) - Nr(3)
      X14 = Nr(1) - Nr(4)
      X23 = Nr(2) - Nr(3)
      X24 = Nr(2) - Nr(4)
      X34 = Nr(3) - Nr(4)
      Y12 = Ns(1) - Ns(2)
      Y13 = Ns(1) - Ns(3)
      Y14 = Ns(1) - Ns(4)
      Y23 = Ns(2) - Ns(3)
      Y24 = Ns(2) - Ns(4)
      Y34 = Ns(3) - Ns(4)
!!
!! C0r(i) & C0s(i); mean coupling operator (membrane-bending coupling)
!!
      C0r(1) =  (0.5D+0 * Y24)
      C0r(3) = -(0.5D+0 * Y24)
      C0r(4) =  (0.5D+0 * Y13)
      C0r(2) = -(0.5D+0 * Y13)
      C0s(3) =  (0.5D+0 * X24)
      C0s(1) = -(0.5D+0 * X24)
      C0s(2) =  (0.5D+0 * X13)
      C0s(4) = -(0.5D+0 * X13)
!!
!! C0t(i); mean coupling operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.2500000000000000D-1 * (X13*S24 - S13*X24)
      QB = 4.1666666666666667D-2 * (X34*S12 - S34*X12)
      QC = 4.1666666666666667D-2 * (X23*S14 - S23*X14)
      C0t(1) = QA - (QB + QC)
      C0t(2) = QA + (QB - QC)
      C0t(3) = QA + (QB + QC)
      C0t(4) = QA - (QB - QC)
!!
      QA = 1.2500000000000000D-1 * (R13*Y24 - Y13*R24)
      QB = 4.1666666666666667D-2 * (R34*Y12 - Y34*R12)
      QC = 4.1666666666666667D-2 * (R23*Y14 - Y23*R14)
      C0t(1) = C0t(1) + QA - (QB + QC)
      C0t(2) = C0t(2) + QA + (QB - QC)
      C0t(3) = C0t(3) + QA + (QB + QC)
      C0t(4) = C0t(4) + QA - (QB - QC)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATQ(NEL)%RES%Area = 0.5D+0 * (R13*S24 - S13*R24)
      Ain = ONE / PLATQ(NEL)%RES%Area
      Hin = 0.5D+0 * Ain
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + B0r(i)*B0r(i) + B0s(i)*B0s(i)
      ENDDO
      Delta = Ain * SQRT (Dx+Dx)
!!
!!
!! III. HOURGLASS CONTROL: Anti-hourglass gradients.
!!
!!
      RETURN
      END
!!_
      SUBROUTINE KHPLQ_STRESS_DIVERGENCE ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 15:50:32
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) :: Mr,Ms,Mt,Fr,Fs,Ft

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
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
      CSQ = MAX (RCL2,RCS2) / Density
      IF (CSQ .GT. 0.0) THEN
        Cv = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PLATQ_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
      PLATQ(NEL)%RES%DTelt = MIN (DTmb,DTts,DTbd)
!!
!! DIVERGENCE OF THE MEMBRANE STRESS RESULTANTS
!! Stress divergence without hourglass control.
!!
      Arr = Arr + (Thickness * QP)
      Ass = Ass + (Thickness * QP)
      DO i = 1,4
        Fr = Arr*B0r(i) + Ars*B0s(i) + Brr*C0r(i) + Brs*C0s(i)
        Fs = Ars*B0r(i) + Ass*B0s(i) + Brs*C0r(i) + Bss*C0s(i)
        Ft = Art*B0r(i) + Ast*B0s(i) + Brt*C0r(i) + Bst*C0s(i)              &
     &     + Drt*B2r(i) + Cst*B1s(i)
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
        PLATQ(NEL)%RES%Xint(i) = Rx*Fr + Sx*Fs + Tx*Ft
        PLATQ(NEL)%RES%Yint(i) = Ry*Fr + Sy*Fs + Ty*Ft
        PLATQ(NEL)%RES%Zint(i) = Rz*Fr + Sz*Fs + Tz*Ft
      ENDDO
!!
!! DIVERGENCE OF THE BENDING STRESS RESULTANTS
!! Includes area weighted terms in divergence calculation.
!!
      DO i = 1,4
        Fr = Brr*B0r(i)+Brs*B0s(i)+Brt*C0t(i)+Art*A0t(i)+Drt*A2t(i)
        Fs = Brs*B0r(i)+Bss*B0s(i)+Bst*C0t(i)+Ast*A0t(i)+Cst*A1t(i)
        Ft = Brt*B0r(i)+Bst*B0s(i)

        Mr = -Nt(i)*Fs + Ns(i)*Ft
        Ms = -Nr(i)*Ft + Nt(i)*Fr
        Mt = -Ns(i)*Fr + Nr(i)*Fs
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
        PLATQ(NEL)%RES%Xint(i+4) = Rx*Mr + Sx*Ms + Tx*Mt
        PLATQ(NEL)%RES%Yint(i+4) = Ry*Mr + Sy*Ms + Ty*Mt
        PLATQ(NEL)%RES%Zint(i+4) = Rz*Mr + Sz*Ms + Tz*Mt
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE BTPLQ_GRADIENT_OPERATOR (NEL,SecID,MatID)
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) ::                                                       &
     &                            Wx(4),Wy(4),Wz(4),X(4),Y(4),Z(4),            &
     &          Vr(4),Vs(4),Vt(4),Wr(4),Ws(4),Wt(4),R(4),S(4),T(4),            &
     &          Theta_r(4),Theta_s(4),Theta_t(4)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Mid-interval position of nodal points, and current translational and
!! rotational velocities.
!!
      DT = 0.5D+0 * PLATQ(NEL)%RES%DTnext
      DO i = 1,4
        X(i)  = Px(i) + (Ux(i) - DT * Vx(i))
        Y(i)  = Py(i) + (Uy(i) - DT * Vy(i))
        Z(i)  = Pz(i) + (Uz(i) - DT * Vz(i))
!!      Vx(i) = Vx(i)
!!      Vy(i) = Vy(i)
!!      Vz(i) = Vz(i)
        Wx(i) = Vx(i+4)
        Wy(i) = Vy(i+4)
        Wz(i) = Vz(i+4)
      ENDDO
!!
!! Construct local basis vectors at the center of the element.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag = 0.5D+0 * SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (0.5D+0 / Rmag)
      Ry = Ry * (0.5D+0 / Rmag)
      Rz = Rz * (0.5D+0 / Rmag)
      Smag = 0.5D+0 * SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (0.5D+0 / Smag)
      Sy = Sy * (0.5D+0 / Smag)
      Sz = Sz * (0.5D+0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'BTPLQ_GRADIENT_OPERATOR.001.00'//                       &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and rotational
!! velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
        Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
        Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
        Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
        Wr(i) = Rx*Wx(i) + Ry*Wy(i) + Rz*Wz(i)
        Ws(i) = Sx*Wx(i) + Sy*Wy(i) + Sz*Wz(i)
        Wt(i) = Tx*Wx(i) + Ty*Wy(i) + Tz*Wz(i)
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R14 = R(1) - R(4)
      R23 = R(2) - R(3)
      R24 = R(2) - R(4)
      R34 = R(3) - R(4)
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S14 = S(1) - S(4)
      S23 = S(2) - S(3)
      S24 = S(2) - S(4)
      S34 = S(3) - S(4)
!!
!! A0t(i); mean averaging operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.25000000000000D-01 * (R13*S24 - S13*R24)
      QB = 4.16666666666667D-02 * (R34*S12 - S34*R12)
      QC = 4.16666666666667D-02 * (R23*S14 - S23*R14)
      A0t(1) = QA - (QB + QC)
      A0t(2) = QA + (QB - QC)
      A0t(3) = QA + (QB + QC)
      A0t(4) = QA - (QB - QC)
!!
!! B0r(i) and B0s(i); mean membrane gradient operators.
!!
      B0r(1) =  (0.5D+0 * S24)
      B0r(3) = -(0.5D+0 * S24)
      B0r(4) =  (0.5D+0 * S13)
      B0r(2) = -(0.5D+0 * S13)
      B0s(3) =  (0.5D+0 * R24)
      B0s(1) = -(0.5D+0 * R24)
      B0s(2) =  (0.5D+0 * R13)
      B0s(4) = -(0.5D+0 * R13)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATQ(NEL)%RES%Area = 0.5D+0 * (R13*S24 - S13*R24)
      Ain = ONE / PLATQ(NEL)%RES%Area
      Hin = 0.5D+0 * Ain
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + B0r(i)*B0r(i) + B0s(i)*B0s(i)
      ENDDO
      Delta = Ain * SQRT (Dx+Dx)
!!
!! I. MEMBRANE: Velocity gradients.
!!
      Vrr  = 0.0
      Vsr  = 0.0
      Vtr  = 0.0
      Vrs  = 0.0
      Vss  = 0.0
      Vts  = 0.0
      Wrav = 0.0
      Wsav = 0.0
      DO i = 1,4
        Vrr  = Vrr  + Vr(i)*B0r(i)
        Vsr  = Vsr  + Vs(i)*B0r(i)
        Vtr  = Vtr  + Vt(i)*B0r(i)
        Vrs  = Vrs  + Vr(i)*B0s(i)
        Vss  = Vss  + Vs(i)*B0s(i)
        Vts  = Vts  + Vt(i)*B0s(i)
        Wrav = Wrav + Wr(i)*A0t(i)
        Wsav = Wsav + Ws(i)*A0t(i)
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
      Grt = Hin * ( Wsav + Vtr)
      Gst = Hin * (-Wrav + Vts)
!!
!! II. BENDING: Rotation gradient.
!!
      Wrr = 0.0
      Wsr = 0.0
      Wtr = 0.0
      Wrs = 0.0
      Wss = 0.0
      Wts = 0.0
      DO i = 1,4
        Wrr = Wrr + Wr(i)*B0r(i)
        Wsr = Wsr + Ws(i)*B0r(i)
        Wtr = Wtr + Wt(i)*B0r(i)
        Wrs = Wrs + Wr(i)*B0s(i)
        Wss = Wss + Ws(i)*B0s(i)
        Wts = Wts + Wt(i)*B0s(i)
      ENDDO
!!
!! Mean bending stretching.
!!
      Hrr =  Ain * Wsr
      Hss = -Ain * Wrs
      Hrs =  Hin * (Wss - Wrr)
      Urs =  Hin * (Wss + Wrr)
!!
!! Zero out higher order strain terms.
!!
      Hrt = 0.0
      Hst = 0.0
      Frt = 0.0
      Est = 0.0
!!
!! III. Anti-hourglass gradients.
!!
      Hr = Ain * (R(1) - R(2) + R(3) - R(4))
      Hs = Ain * (S(1) - S(2) + S(3) - S(4))
      Ht = Ain * (T(1) - T(2) + T(3) - T(4))
!!
      Rr = -Wsav * Ht
      Rs =  Wrav * Ht
!!
      Gr(1) = ( Vr(1)-Vr(2)+Vr(3)-Vr(4) -Vrr*Hr-Vrs*Hs + Rr )*Delta
      Gs(1) = ( Vs(1)-Vs(2)+Vs(3)-Vs(4) -Vsr*Hr-Vss*Hs + Rs )*Delta
      Gt(1) = ( Vt(1)-Vt(2)+Vt(3)-Vt(4) -Vtr*Hr-Vts*Hs )*Delta
!!
      Gr(2) = ( Wr(1)-Wr(2)+Wr(3)-Wr(4) -Wrr*Hr-Wrs*Hs )*Delta
      Gs(2) = ( Ws(1)-Ws(2)+Ws(3)-Ws(4) -Wsr*Hr-Wss*Hs )*Delta
      Gt(2) = ( Wt(1)-Wt(2)+Wt(3)-Wt(4) -Wtr*Hr-Wts*Hs )*Delta
!!
      RETURN
      END
!!_
      SUBROUTINE BTPLQ_DIVERGENCE_OPERATOR (NEL,SecID,MatID)
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) ::                                                       &
     &                            Wx(4),Wy(4),Wz(4),X(4),Y(4),Z(4),            &
     &          Vr(4),Vs(4),Vt(4),Wr(4),Ws(4),Wt(4),R(4),S(4),T(4),            &
     &          Theta_r(4),Theta_s(4),Theta_t(4)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Current position of nodal points, and current translational and
!! rotational velocities.
!!
      DO i = 1,4
        X(i)  = Px(i) + Ux(i)
        Y(i)  = Py(i) + Uy(i)
        Z(i)  = Pz(i) + Uz(i)
      ENDDO
!!
!! Construct local basis vectors at the center of the element.
!!
      Rx = (X(3)-X(1)) + (X(2)-X(4))
      Ry = (Y(3)-Y(1)) + (Y(2)-Y(4))
      Rz = (Z(3)-Z(1)) + (Z(2)-Z(4))
      Sx = (X(3)-X(1)) - (X(2)-X(4))
      Sy = (Y(3)-Y(1)) - (Y(2)-Y(4))
      Sz = (Z(3)-Z(1)) - (Z(2)-Z(4))
      Rmag = 0.5D+0 * SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
      Rx = Rx * (0.5D+0 / Rmag)
      Ry = Ry * (0.5D+0 / Rmag)
      Rz = Rz * (0.5D+0 / Rmag)
      Smag = 0.5D+0 * SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (0.5D+0 / Smag)
      Sy = Sy * (0.5D+0 / Smag)
      Sz = Sz * (0.5D+0 / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'BTPLQ_DIVERGENCE_OPERATOR.001.00'//                     &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
!! Transform position X,Y,Z, translational velocity Vx,Vy,Vz, and rotational
!! velocity Wx,Wy,Wz to local R,S,T-coordinate system.
!!
      DO i = 1,4
        R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
        S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
        T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
      ENDDO
!!
!! Construct useful differences in local coordinates.
!!
      R12 = R(1) - R(2)
      R13 = R(1) - R(3)
      R14 = R(1) - R(4)
      R23 = R(2) - R(3)
      R24 = R(2) - R(4)
      R34 = R(3) - R(4)
      S12 = S(1) - S(2)
      S13 = S(1) - S(3)
      S14 = S(1) - S(4)
      S23 = S(2) - S(3)
      S24 = S(2) - S(4)
      S34 = S(3) - S(4)
!!
!! A0t(i); mean averaging operator (area weighting based on bilinear shape
!! functions).
!!
      QA = 1.25000000000000D-01 * (R13*S24 - S13*R24)
      QB = 4.16666666666667D-02 * (R34*S12 - S34*R12)
      QC = 4.16666666666667D-02 * (R23*S14 - S23*R14)
      A0t(1) = QA - (QB + QC)
      A0t(2) = QA + (QB - QC)
      A0t(3) = QA + (QB + QC)
      A0t(4) = QA - (QB - QC)
!!
!! B0r(i) and B0s(i); mean membrane gradient operators.
!!
      B0r(1) =  (0.5D+0 * S24)
      B0r(3) = -(0.5D+0 * S24)
      B0r(4) =  (0.5D+0 * S13)
      B0r(2) = -(0.5D+0 * S13)
      B0s(3) =  (0.5D+0 * R24)
      B0s(1) = -(0.5D+0 * R24)
      B0s(2) =  (0.5D+0 * R13)
      B0s(4) = -(0.5D+0 * R13)
!!
!! Calculate current element area; Ain = 1.0/Area
!!
      PLATQ(NEL)%RES%Area = 0.5D+0 * (R13*S24 - S13*R24)
      Ain = ONE / PLATQ(NEL)%RES%Area
      Hin = 0.5D+0 * Ain
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,4
        Dx = Dx + B0r(i)*B0r(i) + B0s(i)*B0s(i)
      ENDDO
      Delta = Ain * SQRT (Dx+Dx)
!!
!! III. Anti-hourglass gradients.
!!
      Hr = Ain * (R(1) - R(2) + R(3) - R(4))
      Hs = Ain * (S(1) - S(2) + S(3) - S(4))
      Ht = Ain * (T(1) - T(2) + T(3) - T(4))
!!
      RETURN
      END
!!_
      SUBROUTINE BTPLQ_STRESS_DIVERGENCE (NEL,SecID,MatID)
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 15:50:32
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
      REAL(KIND(0D0)) :: Fr(4),Fs(4),Ft(4)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
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
      CSQ = MAX (RCL2,RCS2) / Density
      IF (CSQ .GT. 0.0) THEN
        Cv = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'BTPLQ_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
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
      PLATQ(NEL)%RES%DTelt = MIN (DTmb,DTts,DTbd)
!!
!! I. MEMBRANE: Artificial velocity hourglass viscosity.
!!
      Volume = Thickness * PLATQ(NEL)%RES%Area
      Vin = ONE / Volume
!!
      Qhg = Density * Dx * HG_Visc * Cv
      PGr = Volume * (PLATQ(NEL)%RES%Pr(1) + Qhg*Gr(1)) * Delta
      PGs = Volume * (PLATQ(NEL)%RES%Ps(1) + Qhg*Gs(1)) * Delta
      PGt = Volume * (PLATQ(NEL)%RES%Pt(1) + Qhg*Gt(1)) * Delta
!!
!! Stress divergence with hourglass control.
!!
      Vrr = Arr - Hr*PGr + (Thickness * QP)
      Vrs = Ars - Hs*PGr
!!
      Vsr = Ars - Hr*PGs
      Vss = Ass - Hs*PGs + (Thickness * QP)
!!
      Vtr = Art - Hr*PGt
      Vts = Ast - Hs*PGt
!!
      Fr(1) =  PGr
      Fr(2) = -PGr
      Fr(3) =  PGr
      Fr(4) = -PGr
!!
      Fs(1) =  PGs
      Fs(2) = -PGs
      Fs(3) =  PGs
      Fs(4) = -PGs
!!
      Ft(1) =  PGt
      Ft(2) = -PGt
      Ft(3) =  PGt
      Ft(4) = -PGt
!!
      DO i = 1,4
        Fr(i) = Fr(i) + B0r(i)*Vrr + B0s(i)*Vrs
        Fs(i) = Fs(i) + B0r(i)*Vsr + B0s(i)*Vss
        Ft(i) = Ft(i) + B0r(i)*Vtr + B0s(i)*Vts
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
      DO i = 1,4
        PLATQ(NEL)%RES%Xint(i) = Rx*Fr(i) + Sx*Fs(i) + Tx*Ft(i)
        PLATQ(NEL)%RES%Yint(i) = Ry*Fr(i) + Sy*Fs(i) + Ty*Ft(i)
        PLATQ(NEL)%RES%Zint(i) = Rz*Fr(i) + Sz*Fs(i) + Tz*Ft(i)
      ENDDO
!!
!! II. BENDING: Area weighted terms in divergence calculation.
!!
      Qr = -Ast + PGs * (Ht * Vin)
      Qs =  Art + PGr * (Ht * Vin)
!!
!! Artificial rotation hourglass viscosity.
!!
      PGr = Volume * (PLATQ(NEL)%RES%Pr(2) + Qhg*Gr(2)) * Delta
      PGs = Volume * (PLATQ(NEL)%RES%Ps(2) + Qhg*Gs(2)) * Delta
      PGt = Volume * (PLATQ(NEL)%RES%Pt(2) + Qhg*Gt(2)) * Delta
!!
!! Bending stress divergence.
!!
      Vrr = -Brs - Hr*PGr
      Vrs = -Bss - Hs*PGr

      Vsr =  Brr - Hr*PGs
      Vss =  Brs - Hs*PGs
!!
      Vtr =  Brt - Hr*PGt
      Vts =  Bst - Hs*PGt
!!
      Fr(1) =  PGr
      Fr(2) = -PGr
      Fr(3) =  PGr
      Fr(4) = -PGr
!!
      Fs(1) =  PGs
      Fs(2) = -PGs
      Fs(3) =  PGs
      Fs(4) = -PGs
!!
      Ft(1) =  PGt
      Ft(2) = -PGt
      Ft(3) =  PGt
      Ft(4) = -PGt
!!
      DO i = 1,4
        Fr(i) = Fr(i) + B0r(i)*Vrr + B0s(i)*Vrs + Qr*A0t(i)
        Fs(i) = Fs(i) + B0r(i)*Vsr + B0s(i)*Vss + Qs*A0t(i)
        Ft(i) = Ft(i) + B0r(i)*Vtr + B0s(i)*Vts
      ENDDO
!!
!! Transform internal forces to global coordinates, and accumulate element
!! divergence results in local element array.
!!
      DO i = 1,4
        PLATQ(NEL)%RES%Xint(i+4) = Rx*Fr(i) + Sx*Fs(i) + Tx*Ft(i)
        PLATQ(NEL)%RES%Yint(i+4) = Ry*Fr(i) + Sy*Fs(i) + Ty*Ft(i)
        PLATQ(NEL)%RES%Zint(i+4) = Rz*Fr(i) + Sz*Fs(i) + Tz*Ft(i)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE PLATQ_STRESS_INTEGRATION ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 17-MAR-1991 12:27:16
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)), SAVE :: Beta(20),Awght(20),SQRT6o1,SQRT5o6
      REAL(KIND(0D0))       :: Int_Eng

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
      INTEGER, PARAMETER :: Ibgn(5) = (/1,3,6,10,15/)
      INTEGER, PARAMETER :: Iend(5) = (/2,5,9,14,20/)

      REAL(KIND(0D0)), PARAMETER :: POS3RD = (+1.0D0/3.0D0)
      REAL(KIND(0D0)), PARAMETER :: POS6TH = (+1.0D0/6.0D0)
      REAL(KIND(0D0)), PARAMETER :: NEG6TH = (-1.0D0/6.0D0)

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Non-dimensional thickness coordinates (trapezodial rule).
!!
      DATA (Beta(i),i= 1, 2) /-.50D+0, .500D+0/
      DATA (Beta(i),i= 3, 5) /-.50D+0, .000D+0, .500D+0/
      DATA (Beta(i),i= 6, 9) /-.50D+0,  NEG6TH,  POS6TH, .50D+0/
      DATA (Beta(i),i=10,14) /-.50D+0,-.250D+0, .000D+0, .25D+0, .50D+0/
      DATA (Beta(i),i=15,20) /-.50D+0,-.300D+0,-.100D+0, .10D+0, .30D+0, .50D+0/
!!
!! Membrane thickness integration weights.
!!
      DATA (Awght(i),i= 1, 2) /.500D+0, .500D+0/
      DATA (Awght(i),i= 3, 5) /.250D+0, .500D+0, .250D+0/
      DATA (Awght(i),i= 6, 9) / POS6TH,  POS3RD,  POS3RD,  POS6TH/
      DATA (Awght(i),i=10,14) /.125D+0, .250D+0, .250D+0, .250D+0, .125D+0/
      DATA (Awght(i),i=15,20) /.100D+0, .200D+0, .200D+0, .200D+0, .200D+0, .100D+0/
!!
!! Define constants.
!!
      IF (FIRST) THEN
        SQRT6o1 = SQRT (6.0D+0/1.0D+0)
        SQRT5o6 = SQRT (5.0D+0/6.0D+0)
        FIRST = .FALSE.
      ENDIF
!!
!! Compute current element thickness based on constant volume.
!!
      Thickness = SECTION_2D(SecID)%Thickness *                                &
     &  PLATQ(NEL)%PAR%Area / PLATQ(NEL)%RES%Area
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
      Aspect_Ratio = Thickness / Rmag
      R_Factor = MIN (SQRT5o6, SQRT6o1*Aspect_Ratio)
      Aspect_Ratio = Thickness / Smag
      S_Factor = MIN (SQRT5o6, SQRT6o1*Aspect_Ratio)
      Factor = MAX (R_Factor,S_Factor)
!!
!! Apply shear correction factors to transverse shear strains.
!!
      Grt = R_Factor * Grt
      Gst = S_Factor * Gst
      Hrt = R_Factor * Hrt
      Hst = S_Factor * Hst
      Est = S_Factor * Est
      Frt = R_Factor * Frt
!!
!! Initialize the mean stresses Arr, Ass, Ars, Art, and Ast.
!!
      Arr = 0.0
      Ass = 0.0
      Ars = 0.0
      Art = 0.0
      Ast = 0.0
!!
!! Initialize the moment stresses Brr, Bss, Brs, Brt, and Bst.
!!
      Brr = 0.0
      Bss = 0.0
      Brs = 0.0
      Brt = 0.0
      Bst = 0.0
!!
!! Initialize stress state pointer Ist, state variable pointer Isv and
!! state variable interval Nsv.
!!
      Ist  = PLATQ(NEL)%PAR%Ist
      Isv  = PLATQ(NEL)%PAR%Isv
      Nsv  = MATERIAL(MatID)%Nsv
!!
!! Select between internal integration rule and user-supplied integration
!! rule. Calculate Xwght for use in constructing moment integration end-
!! point weighting Ewght.
!!
      IF (SECTION_2D(SecID)%Ipts .NE. 0) THEN
        Ione = Ibgn(SECTION_2D(SecID)%Ipts-1)
        Itwo = Iend(SECTION_2D(SecID)%Ipts-1)
        Xwght = POS6TH * Thickness / DBLE ((SECTION_2D(SecID)%Ipts-1)**2)
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
          Vrt = Grt + Eta * Hrt
          Vst = Gst + Eta * Hst
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
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (41)
            CALL MATERIAL_41                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID,RCL2,RCS2                                                &
     &          )
          CASE (45)
            CALL MATERIAL_45                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (47)
            CALL MATERIAL_47                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATQ_STRESS_INTEGRATION.001.00'//                      &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Via a MATERIAL Record References An '                   &
     &              //'Invalid Material Type:'//MSG2                           &
     &          )
          END SELECT
!!
!! Internal energy from time n+1/2 to time n+1.
!!
          Int_Eng = Int_Eng                                                    &
     &              +  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                 &
     &              + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                &
     &              + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                &
     &              + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Accumulate layer internal energy density.
!!
          Int_Eng = (0.5D+0*PLATQ(NEL)%RES%DTnext) * AWght(i) * Int_Eng
!!
          PLATQ(NEL)%RES%Int_Eng = PLATQ(NEL)%RES%Int_Eng + Int_Eng
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
          Brt = Brt + (Eta * Awght(i) + Ewght) * STRESS(5,Ist)
          Bst = Bst + (Eta * Awght(i) + Ewght) * STRESS(6,Ist)
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
          Vrt = Grt + Eta * Hrt
          Vst = Gst + Eta * Hst
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
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (41)
            CALL MATERIAL_41                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID,RCL2,RCS2                                                &
     &          )
          CASE (45)
            CALL MATERIAL_45                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE (47)
            CALL MATERIAL_47                                                   &
     &          (                                                              &
     &          STRESS(1,Ist),                                                 &
     &          STATE_VARIABLES(Isv),                                          &
     &          PLATQ(NEL)%RES%DTnext,                                         &
     &          Vrr,Vss,Vrs,Vrt,Vst,Xrs,                                       &
     &          MatID                                                          &
     &          )
          CASE DEFAULT
            WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATQ_STRESS_INTEGRATION.001.00'//                      &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'Via a MATERIAL Record References An '                   &
     &              //'Invalid Material Type:'//MSG2                           &
     &          )
          END SELECT
!!
!! Internal energy from time n+1/2 to time n+1.
!!
          Int_Eng = Int_Eng                                                    &
     &              +  Vrr*STRESS(1,Ist)  +  Vss*STRESS(2,Ist)                 &
     &              + (Vrs*STRESS(4,Ist)) + (Vrs*STRESS(4,Ist))                &
     &              + (Vrt*STRESS(5,Ist)) + (Vrt*STRESS(5,Ist))                &
     &              + (Vst*STRESS(6,Ist)) + (Vst*STRESS(6,Ist))
!!
!! Accumulate layer internal energy density.
!!
          Int_Eng = (0.5D+0*PLATQ(NEL)%RES%DTnext) * Weight * Int_Eng
!!
          PLATQ(NEL)%RES%Int_Eng = PLATQ(NEL)%RES%Int_Eng + Int_Eng
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
          Brt = Brt + (Eta * Weight) * STRESS(5,Ist)
          Bst = Bst + (Eta * Weight) * STRESS(6,Ist)
!!
          Ist = Ist + 1
          Isv = Isv + Nsv
        ENDDO
!!
      ELSE
        WRITE (MSG1,'(I8)') PLATQ(NEL)%PAR%EleID
        WRITE (MSG2,'(I8)') SECTION_2D(SecID)%SecID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PLATQ_STRESS_INTEGRATION.002.00'//                      &
     &          MSGL//'P4EL (4-Node Plate) Element ID:'//MSG1//                &
     &          MSGL//'References PSECTION ID:'//MSG2//                        &
     &          MSGL//'With Both Ipts And Irule Equal to Zero.'                &
     &          )
      ENDIF
!!
!! Linearly varying components of the transverse shear. Provides coupling
!! between rotations and reference surface warping to generate the twisting
!! moment Brs = Mrs/Thickness.
!!
      Qmod = 12.0D+0 * PLATQ(NEL)%RES%DTnext * RCS2
      Cst = PLATQ(NEL)%RES%Shear(1) + Qmod * Est
      Drt = PLATQ(NEL)%RES%Shear(2) + Qmod * Frt
!!
      PLATQ(NEL)%RES%Shear(1) = Cst
      PLATQ(NEL)%RES%Shear(2) = Drt
!!
!! Convert mean stresses to stress resultants. Apply shear correction factors
!! to transverse-shear stresses.
!!
      Arr = Thickness * Arr
      Ass = Thickness * Ass
      Ars = Thickness * Ars
      Art = (R_Factor*Thickness) * Art
      Ast = (S_Factor*Thickness) * Ast
      Brr = Thickness * Brr
      Bss = Thickness * Bss
      Brs = Thickness * Brs
      Brt = (R_Factor*Thickness) * Brt
      Bst = (S_Factor*Thickness) * Bst
      Cst = (S_Factor*Thickness) * Cst
      Drt = (R_Factor*Thickness) * Drt
!!
      RETURN
      END
!!_
      SUBROUTINE PLATQ_HOURGLASS_FORCES ( NEL,SecID,MatID )
!!
!! Copyright (c) by KEY Associates, 26-MAY-1991 13:19:09
!!
!! Purpose: Increment stiffness based hourglass control forces.
!!
      USE shared_common_data
      USE platq_
      USE material_
      USE section_2d_
      USE node_
      USE motion_
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
!! Local variables.
      REAL(KIND(0D0)) :: Q1(2),Q2(2)

      REAL(KIND(0D0)) :: Nr,Ns,Nt

      COMMON /PLATQ_COMMON/      &
     &  A0t(4),A1t(4),A2t(4),    & ! mean, xi1, xi2 averaging operators
     &  B0r(4),B0s(4),           & ! mean           gradient  operators
     &  B1s(4),                  & ! xi1-linear     gradient  operator
     &  B2r(4),                  & ! xi2-linear     gradient  operator
     &  C0r(4),C0s(4),C0t(4),    & ! mean           coupling  operators
     &  Arr,Ass,Ars,Art,Ast,     & ! mean           membrane  stresses
     &  Brr,Bss,Brs,Brt,Bst,     & ! mean           bending   stresses
     &  Cst,                     & ! xi1-linear     membrane  shear
     &  Drt,                     & ! xi2-linear     membrane  shear
     &  Grr,Gss,Grs,Grt,Gst,Qrs, & ! mean           membrane  stretching
     &  Hrr,Hss,Hrs,Hrt,Hst,Urs, & ! mean           bending   stretching
     &  Est,                     & ! xi1            membrane  stretching
     &  Frt,                     & ! xi2            membrane  stretching
     &  Nr(4),Ns(4),Nt(4),       & ! unit vectors normal to corners
     &  Delta,Factor,Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz,Thickness,Rmag,Smag,   &
     &  Hr,Hs,Ht,Gr(2),Gs(2),Gt(2),Px(8),Py(8),Pz(8),Ux(8),Uy(8),Uz(8),&
     &  Vx(8),Vy(8),Vz(8),RCL2,RCS2
!!
!$OMP THREADPRIVATE (/PLATQ_COMMON/)
!!
!! Retrieve material properties.
!!
      HG_Stiff = MATERIAL(MatID)%PVAL(5)
!!
!! Rotate and increment elastic anti-hourglassing forces.
!!
      dQrs = 0.5D+0 * PLATQ(NEL)%RES%DTnext * Qrs
      DO i = 1,2
        Q1(i) =  dQrs * PLATQ(NEL)%RES%Ps(i)
        Q2(i) = -dQrs * PLATQ(NEL)%RES%Pr(i)
      ENDDO
      Qe = PLATQ(NEL)%RES%DTnext * HG_Stiff * RCS2
      DO i = 1,2
        PLATQ(NEL)%RES%Pr(i) = PLATQ(NEL)%RES%Pr(i) + Qe * Gr(i) + Q1(i)
        PLATQ(NEL)%RES%Ps(i) = PLATQ(NEL)%RES%Ps(i) + Qe * Gs(i) + Q2(i)
        PLATQ(NEL)%RES%Pt(i) = PLATQ(NEL)%RES%Pt(i) + Qe * Gt(i)
      ENDDO
      DO i = 1,2
        Q1(i) =  dQrs * PLATQ(NEL)%RES%Ps(i)
        Q2(i) = -dQrs * PLATQ(NEL)%RES%Pr(i)
      ENDDO
      DO i = 1,2
        PLATQ(NEL)%RES%Pr(i) = PLATQ(NEL)%RES%Pr(i) + Q1(i)
        PLATQ(NEL)%RES%Ps(i) = PLATQ(NEL)%RES%Ps(i) + Q2(i)
      ENDDO
!!
      RETURN
      END
