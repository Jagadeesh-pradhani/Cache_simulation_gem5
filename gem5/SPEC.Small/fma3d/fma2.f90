      SUBROUTINE SOLVE
!!
!! Copyright (c) by KEY Associates, 9-FEB-1991 21:03:29
!!
!! Purpose: Integrate forward in time to obtain the transient response.
!! This module is entered after all the input data has been assembled.
!! The solution is controlled from this module. The central difference
!! time integrator is contained in this module along with all of the
!! calls to routines which calculate the forces acting on the "nodal
!! masses."
!!
      USE shared_common_data
!!
!! The complete simulation data set.
!!
      USE indx_;          USE node_;          USE tabulated_function_;
      USE beam_;          USE coord_;         USE sliding_interface_;
      USE force_;         USE hexah_;         USE nonreflecting_bc_;
      USE penta_;         USE tetra_;         USE nodal_point_mass_;
      USE lsold_;         USE membt_;         USE constrained_node_;
      USE membq_;         USE truss_;         USE displacement_bc_;
      USE platt_;         USE platq_;         USE rigid_body_mass_;
      USE motion_;        USE stress_;        USE enumerated_sets_;
      USE spring_;        USE damper_;        USE contact_surface_;
      USE segment_;       USE tied_bc_;       USE state_variables_;
      USE results_;       USE gauge1d_;       USE rigid_wall_bc_;
      USE gauge2d_;       USE gauge3d_;       USE contact_node_;
      USE node_set_;      USE force_bc_;      USE sliding_node_;
      USE material_;      USE layering_;      USE segment_set_;
      USE massprop_;      USE spring_bc_;     USE element_set_;
      USE damper_bc_;     USE spot_weld_;     USE periodic_bc_;
      USE qa_record_;     USE nrbc_data_;     USE pressure_bc_;
      USE plate_pair_;    USE section_2d_;    USE section_1d_;
      USE rigid_body_;    USE body_force_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          ERRMSG*36               ! ERRor_MeSsaGe tag
      LOGICAL                         &
     &          REZONING,             & ! Controls use of rezoning             &
     &          SUBCYCLING              ! Controls use of subcycling
!!
      LOGICAL :: GO_TO_401 = .FALSE.
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SOLVE.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Record time spent processing input data.
!!
!SPEC_CPU2000      CALL TIMER (1)
!!
!! Initialize time information with controling tabulated function index.
!!
      TIMSIM%DTcntrl = CONTROL%DTTABF
!!
!! Check for time integration subcycling and rezoning.
!!
      SUBCYCLING = (CONTROL%SUBCYC .NE. 0)
      REZONING   = (CONTROL%REZONE .NE. 0)
!!
!! If the ending simulation time was not defined with a STOP card, scan
!! print and results requests for ending times. Use the maximum found.
!!
      IF (TIMSIM%Stop .EQ. 0.0) THEN
        TIMSIM%Stop = PRINT%End
        DO i = 1,NUMRF
          TIMSIM%Stop = MAX (TIMSIM%Stop, RESULTS(i)%End)
        ENDDO
      ENDIF
!!
!! Check restart status.
!!
      IF (CONTROL%RDSTAR .EQ. 1) THEN
!!
!! Initialize Next-Time-to-Print value.
!!
        IF (TIMSIM%Total .GT. PRINT%End) THEN
          PRINT%Time = 2.0 * TIMSIM%Stop
        ELSE
          PRINT%Time = PRINT%Begin
          IF (PRINT%Delta .GT. 0.0) THEN
            DO WHILE (PRINT%Time .LE. TIMSIM%Total)
              PRINT%Time = PRINT%Time + PRINT%Delta
            ENDDO
          ENDIF
        ENDIF
!!
!! Initialize Next-Time-to-Write-Plotting_Database values.
!!
        IF (NUMRF .GT. 0) THEN
          DO i = 1,NUMRF
            IF (TIMSIM%Total .GT. RESULTS(i)%End) THEN
              RESULTS(i)%Time = 2.0 * TIMSIM%Stop
            ELSE
              RESULTS(i)%Time = RESULTS(i)%Begin
              IF (RESULTS(i)%Delta .GT. 0.0) THEN
                DO WHILE (RESULTS(i)%Time .LE. TIMSIM%Total)
                  RESULTS(i)%Time = RESULTS(i)%Time + RESULTS(i)%Delta
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
!!
!! Jump to "resume" address in time integration loop.
!!
        IF (CONTROL%RZSTAR .EQ. 1) THEN
          GO TO 200
        ELSE
          GO_TO_401 = .TRUE.
          GO TO 400
        ENDIF
!!
      ENDIF
!!
!! *************************************************************************
!! *** PART I. INITIALIZATION AND ACCELERATIONS FOR TIME EQUAL TO ZERO. ****
!! *************************************************************************
!!
!! Initialize time information.
!!
      TIMSIM%Step   = 0
      TIMSIM%Cycle  = 0
      TIMSIM%Total  = 0.0
      TIMSIM%DTlast = 0.0
      TIMSIM%DTnext = 0.0
      CALL TIME_INITIALIZATION
!!
!!
      DO N = 1,NUMRT
        NODE(N)%Time   = 0.0
        NODE(N)%DTlast = 0.0
        NODE(N)%DTnext = 0.0
      ENDDO
!!
!!
!! Initialize time-to-print with time-to-begin-printing.
!!
      IF (PRINT%End .EQ. 0.0) THEN
        PRINT%Time = 2.0D+0 * TIMSIM%Stop
      ELSE
        PRINT%Time = PRINT%Begin
      ENDIF
!!
!! Initialize time-to-write-Plotting_Database with time-to-begin-writing.
!!
      IF (NUMRF .GT. 0) THEN
        DO i = 1,NUMRF
          RESULTS(i)%Time = RESULTS(i)%Begin
        ENDDO
      ENDIF
!!
!! Obtain initial stresses.
!!
      IF (CONTROL%RDSTAR .EQ. 2) CALL READ_INITIAL_CONDITIONS
!!
!! Initialize material constants and auxillary storage.
!!
      CALL INITIALIZE_MATERIALS
!!
!! CONCENTRATED NODAL MASSES
!! Initialize mass vector with concentrated masses.
!!
      IF (NUMCM .GT. 0) CALL INIT_CONCENTRATED_MASSES
!!
!! INTERNAL FORCES.
!! Initialize internal and external forces. Internal forces come from
!! the divergence of the stresses. External forces come from applied
!! loads, sliding interfaces, et cetera.
!!
!!
      DO N = 1,NUMRT
        FORCE(N) = force_type (0,0,0,0,0,0)
      ENDDO
!!
!!
!! Initial internal forces from the divergence of the initial stress state.
!!
      CALL ELEMENT_INITIALIZATION
!!
!! Record time spent setting up initial accelerations (1 of 2 locations).
!!
!SPEC_CPU2000      CALL TIMER (2)
!!
!! EXTERNAL FORCES.
!! Compute the initial body forces, surface tractions, concentrated forces,
!! and active boundary conditions (spring, damper, and nonreflecting BC's).
!!
      CALL EXTERNAL_FORCES
!!
!! Initialize strain gauges.
!!
      CALL STRAIN_GAUGE_INITIALIZATION
!!
!! Initilize mass and inertia for rigid body domains whose properties are
!! based on distributed nodal masses. (Nodal masses were computed during
!! element initialization above.)
!!
      IF (NUMRB .GT. 0) CALL INITIALIZE_RIGID_BODY_MASS
!!
!! Initialize effective mass for constrained nodal point triples.
!!
      IF (NUMNC .GT. 0) CALL BUILD_EFFECTIVE_MASS
!!
!! Report the simulations mass properties.
!!
!SPEC_CPU2000      CALL TOTAL_MASS_REPORT ( IO_UNIT%LELO )
!!
!! TIME STEP OPERATIONS AND GROUPING FOR SUBCYCLING
!! Minimum critical time step for central difference integration algorithm.
!!
      CALL NEXT_TIME_STEP ( TIMSIM%DTnext )
!!
!! Minimum acceptable time step. (Time step at which calculations will
!! be terminated in the event the integration time step TIMSIM%DTnext
!! drops too far.)
!!
      TIMSIM%DTlim = PARAMVALUE%DTratio * TIMSIM%DTnext
!!
!! Report the smallest element critical time steps from each element group.
!!
!!
!! Partition elements into subcycling groups.
!!
      CALL SUBCYCLING_PARTITION ( SUBCYCLING )
!!
!! Initialize subcycling counter.
!!
      TIMSIM%Cycle = TIMSIM%Cymax
!!
!! The following on/off-function is used to make nodal do-loop if-tests
!! more efficient. The nodal point value of NODE(*)%On is either .TRUE.
!! (on) or .FALSE. (off). Put in the form of a floating point number
!! (1.0 or 0.0) it can be used as a multiplier in nodal calculations.
!! For example,
!!
!!      NODE(N)%xxx = FLOAT (MAX (0, 1 - MOD (TIMSIM%Cycle,NODE(N)%ISI)))
!!
!!
!!
!!
!! INVERT MASS MATRIX
!! Check for non-zero nodal masses.
!!
      ERROR%COUNT = 0
      DO N = 1,NUMRT
        NODE(N)%On = .TRUE.  !  (MOD (TIMSIM%Cycle,NODE(N)%ISI) .EQ. 0)
        IF (NODE(N)%Mass .GT. 0.0) THEN
          NODE(N)%Minv = ONE / NODE(N)%Mass
        ELSE
          WRITE (MSG1,'(I8)') NODE(N)%ID
          IF(N .LE. NUMNP) THEN
            IF (NODE(N)%Mass .EQ. 0.0) THEN
              ERRMSG = 'Has A Zero Translational Mass.'
            ELSE
              ERRMSG = 'Has A Negative Translational Mass.'
            ENDIF
          ELSE
            IF (NODE(N)%Mass .EQ. 0.0) THEN
              ERRMSG = 'Has A Zero Rotational Inertia.'
            ELSE
              ERRMSG = 'Has A Negative Rotational Inertia.'
            ENDIF
          ENDIF
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SOLVE.001.01'//                                         &
     &          MSGL//'Nodal Point ID:'//MSG1//                                &
     &          MSGL//ERRMSG                                                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SOLVE.001.02'//                                         &
     &          MSGL//'Total Number Of Zero Masses/Inertias:'//MSG1            &
     &          )
      ENDIF
!!
!! INITIALIZE ENERGY BALANCE
!! Compute initial kinetic energy, and initialize internal energy and
!! external energy.
!!
      ENGK = 0.0
!!
!!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N) REDUCTION(+:ENGK)
      DO N = 1,NUMRT
        ENGK = ENGK + NODE(N)%Mass *                                           &
     &        (                                                                &
     &        MOTION(N)%Vx * MOTION(N)%Vx +                                    &
     &        MOTION(N)%Vy * MOTION(N)%Vy +                                    &
     &        MOTION(N)%Vz * MOTION(N)%Vz                                      &
     &        )
      ENDDO
!!
!$OMP END PARALLEL DO  
      ENERGY%Kinetic   = (0.5D+0 * ENGK)
      ENERGY%External  = (0.5D+0 * ENGK)
      ENERGY%Internal  = 0.0
      ENERGY%Sliding   = 0.0
      ENERGY%Contact   = 0.0
      ENERGY%Hourglass = 0.0
      ENERGY%Bulk_Vis  = 0.0
!!
!! Record time spent setting up initial accelerations (2 of 2 locations).
!!
!SPEC_CPU2000      CALL TIMER (2)
!!
!! *************************************************************************
!! ******************** PART II. TIME INTEGRATION LOOP. ********************
!! *************************************************************************
!!
!!         T O P   O F   T I M E   I N T E G R A T I O N   L O O P.
!!
 200    CONTINUE
!!
!! Compute accelerations.
!!
!!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N)
      DO N = 1,NUMRT
        MOTION(N)%Ax = NODE(N)%Minv * (FORCE(N)%Xext-FORCE(N)%Xint)
        MOTION(N)%Ay = NODE(N)%Minv * (FORCE(N)%Yext-FORCE(N)%Yint)
        MOTION(N)%Az = NODE(N)%Minv * (FORCE(N)%Zext-FORCE(N)%Zint)
      ENDDO
!$OMP END PARALLEL DO
!!
!!
!! Recompute the acceleration of nodal points with concentrated inertia.
!!
      IF (NUMCM .GT. 0) CALL UPDATE1_CONCENTRATED_MASSES
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000      CALL TIMER (3)
!!
!! KINEMATIC CONSTRAINTS
!! Modify nodal point accelerations (and, therefore, motion) due to
!! explicit kinematic constraints.
!!
      CALL KINEMATIC_CONSTRAINTS
!!
!! SUBCYCLING, REZONING, AND RESULTS OUTPUT MANAGEMENT
!! Check for completion of subcycling, that is, for all nodal/element
!! clocks back in synchronization. (Note that TIMSIM%Cycle is initialized
!! to Cymax to cause this section to be executed for TIMSIM%Step equal to
!! zero.)
!!
 400  CONTINUE ! A stop along the way to the real restart-resume location.
!!
      IF (TIMSIM%Cycle.EQ.TIMSIM%Cymax .OR. GO_TO_401) THEN
!!
!! Now, off to the real restart resume location -- all this to avoid
!! a compiler warning message even though we know exactly what we are
!! doing and exactly where we want to go.
!!
        IF (GO_TO_401) GO TO 401
!!
!! For time equal to zero the following operations were performed in Part I.
!!
        IF (TIMSIM%Step .GT. 0) THEN
!!
!! Check for elements requiring rezoning.
!!
          IF (REZONING) THEN
            IF (MOD(TIMSIM%Step,REZONE%INTERVAL) .EQ. 0) THEN
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000              CALL TIMER (3)
!!!           CALL H_ADAPTIVE_REMESHING
!!
!! Record time spent in evaluating rezoning criteria.
!!
!SPEC_CPU2000              CALL TIMER (26)
            ENDIF
          ENDIF
!!
!! Time step for next integration step. (This time step will be used for the
!! next round of subcycling. It will be reset when this section is re-entered
!! at the completion of the next round of subcycling.)
!!
          TIMSIM%DTlast = TIMSIM%DTnext
          CALL NEXT_TIME_STEP ( TIMSIM%DTnext )
!!
!! Re-partition elements into new subcycling groups as required.
!!
          CALL SUBCYCLING_PARTITION ( SUBCYCLING )
!!
!! Reset time step minimums and maximums to extreme values.
!!
          CALL TIME_INITIALIZATION
!!
        ENDIF
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000        CALL TIMER (3)
!!
!! A calculation which is restarted picks up here (CONTROL%RDSTAR .EQ. 1).
!! (GO_TO_401 is set to false to make the hop to this location a one-time
!! event.)
!!
 401    CONTINUE; GO_TO_401 = .FALSE.
!!
!! Digress for processing all output requests.
!!
        CALL OUTPUT_PROCESSING
!!
!! Check to see if this is a "data check."
!!
        IF (CONTROL%DCHECK .NE. 0) RETURN
!!
!! Check to see if the simulation is completed, if not, resume calculation.
!!
        IF (TIMSIM%Total .GE. TIMSIM%Stop) RETURN
!!
!! Increment global step counter and reset subcycling counter to zero.
!!
        TIMSIM%Step = TIMSIM%Step + 1
        TIMSIM%Cycle = 0
!!
      ENDIF
!!
!!            S T A R T   O F   N E W   T I M E   S T E P
!!
!! ENERGY BALANCE CALCULATION.
!! External Energy (ENG1): Contribution to External_Energy from external
!! forces. Includes surface traction work, sliding interface work (In a
!! "good" calculation the sliding interface forces are self-equilibrating
!! and thus do not generate energy.), and work from kinematic boundary
!! conditions which involve movement.
!!
!! Internal Energy (ENG2): Contribution to Internal_Energy from internal
!! forces. Internal forces come from the divergence of the stress field
!! within each element.
!!
      ENG1 = 0.0
      ENG2 = 0.0
      IF (SUBCYCLING) THEN
!!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,DTaver) REDUCTION(+:ENG1,ENG2)
!!
        DO N = 1,NUMRT
          IF (NODE(N)%On) THEN
            ENG1 = ENG1 + NODE(N)%DTlast *                                     &
     &          (                                                              &
     &          FORCE(N)%Xext * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yext * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zext * MOTION(N)%Vz                                   &
     &          )
            ENG2 = ENG2 + NODE(N)%DTlast *                                     &
     &          (                                                              &
     &          FORCE(N)%Xint * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yint * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zint * MOTION(N)%Vz                                   &
     &          )
!!
!! CENTRAL DIFFERENCE INTEGRATION
!! Update velocities and displacements. (A central difference integration of
!! angular velocities to get angles gives nonsense, and, thus, in principle
!! the loop limit should be NUMNP. However, kinematic boundary conditions
!! that fix the direction about which rotations can occur result in meaningful
!! angles when integrated.)
!!
            DTaver = 0.5D+0 * (NODE(N)%DTlast + NODE(N)%DTnext)
            MOTION(N)%Vx = MOTION(N)%Vx + DTaver * MOTION(N)%Ax
            MOTION(N)%Vy = MOTION(N)%Vy + DTaver * MOTION(N)%Ay
            MOTION(N)%Vz = MOTION(N)%Vz + DTaver * MOTION(N)%Az
            MOTION(N)%Ux = MOTION(N)%Ux + NODE(N)%DTnext * MOTION(N)%Vx
            MOTION(N)%Uy = MOTION(N)%Uy + NODE(N)%DTnext * MOTION(N)%Vy
            MOTION(N)%Uz = MOTION(N)%Uz + NODE(N)%DTnext * MOTION(N)%Vz
            NODE(N)%Time = NODE(N)%Time + NODE(N)%DTnext
            NODE(N)%DTlast = NODE(N)%DTnext

            ENG1 = ENG1 + NODE(N)%DTnext *                                     &
     &          (                                                              &
     &          FORCE(N)%Xext * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yext * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zext * MOTION(N)%Vz                                   &
     &          )
            ENG2 = ENG2 + NODE(N)%DTnext *                                     &
     &          (                                                              &
     &          FORCE(N)%Xint * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yint * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zint * MOTION(N)%Vz                                   &
     &          )
          ENDIF
        ENDDO
!!
!$OMP END PARALLEL DO
!!
      ELSE
!!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N) REDUCTION(+:ENG1,ENG2)
!!
        DO N = 1,NUMRT
          ENG1 = ENG1 + TIMSIM%DTlast *                                        &
     &          (                                                              &
     &          FORCE(N)%Xext * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yext * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zext * MOTION(N)%Vz                                   &
     &          )
          ENG2 = ENG2 + TIMSIM%DTlast *                                        &
     &          (                                                              &
     &          FORCE(N)%Xint * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yint * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zint * MOTION(N)%Vz                                   &
     &          )

          DTaver = 0.5D+0 * (TIMSIM%DTlast + TIMSIM%DTnext)

          MOTION(N)%Vx = MOTION(N)%Vx + DTaver * MOTION(N)%Ax
          MOTION(N)%Vy = MOTION(N)%Vy + DTaver * MOTION(N)%Ay
          MOTION(N)%Vz = MOTION(N)%Vz + DTaver * MOTION(N)%Az
          MOTION(N)%Ux = MOTION(N)%Ux + TIMSIM%DTnext * MOTION(N)%Vx
          MOTION(N)%Uy = MOTION(N)%Uy + TIMSIM%DTnext * MOTION(N)%Vy
          MOTION(N)%Uz = MOTION(N)%Uz + TIMSIM%DTnext * MOTION(N)%Vz
          NODE(N)%Time = NODE(N)%Time + TIMSIM%DTnext

          ENG1 = ENG1 + TIMSIM%DTnext *                                        &
     &          (                                                              &
     &          FORCE(N)%Xext * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yext * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zext * MOTION(N)%Vz                                   &
     &          )
          ENG2 = ENG2 + TIMSIM%DTnext *                                        &
     &          (                                                              &
     &          FORCE(N)%Xint * MOTION(N)%Vx +                                 &
     &          FORCE(N)%Yint * MOTION(N)%Vy +                                 &
     &          FORCE(N)%Zint * MOTION(N)%Vz                                   &
     &          )
        ENDDO
!!
!$OMP END PARALLEL DO
!!
      ENDIF
!!
      ENERGY%External = ENERGY%External + 0.5D+0 * ENG1
      ENERGY%Internal = ENERGY%Internal + 0.5D+0 * ENG2
!!
!! Integrate rotational motion of nodal points with concentrated inertia.
!! (The inertia tensor is also rotated from the configuration at time n-1
!! to the configuration at time n.)
!!
      IF (NUMCM .GT. 0) CALL UPDATE2_CONCENTRATED_MASSES
!!
!! Increment simulation time and subcycle step counter.
!!
      TIMSIM%Cycle = TIMSIM%Cycle + 1
      TIMSIM%Total = TIMSIM%Total + TIMSIM%DTnext
!!
!! Reset nodal point on/off-function for use in nodal if-tests.
!!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(N)
      IF (SUBCYCLING) THEN
!!
!$OMP DO
!!
        DO N = 1,NUMRT
          NODE(N)%On = (MOD (TIMSIM%Cycle,NODE(N)%ISI) .EQ. 0)
        ENDDO
!!
!$OMP END DO NOWAIT
!!
      ENDIF
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000      CALL TIMER (3)
!!
!! INTERNAL FORCES
!! Clear internal and external forces. Internal forces come from the
!! divergence of the stresses. External forces come from applied loads,
!! sliding interfaces, et cetera.
!!
!$OMP DO
!!
      DO N = 1,NUMRT
        FORCE(N) = force_type (0,0,0,0,0,0)
      ENDDO
!!
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!!
!! Compute new internal forces and local element critical time steps for next
!! time integration step.
!!
      CALL INTERNAL_FORCES
!!
!! EXTERNAL FORCES.
!! Compute new body forces, surface tractions, concentrated forces, and
!! active boundary conditions (spring, damper, and nonreflecting BC's).
!!
      CALL EXTERNAL_FORCES
!!
!! Update strain gauges. Note: Timing calls are contained within the
!! following module.
!!
      IF ((NUMG1+NUMG2+NUMG3) .GT. 0) CALL STRAIN_GAUGE_INTEGRATION
!!
!! Check for time step drop.
!!
      IF (TIMSIM%DTmin .LT. TIMSIM%DTlim) THEN
        WRITE (MSG1,'(I8)') TIMSIM%Ncrit
        WRITE (MSGF,'(1PE12.4)') TIMSIM%DTmin
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SOLVE.002.00'//                                         &
     &          MSGL//'Time Step Has Dropped Several Orders'//                 &
     &          MSGL//'In Magnitude From First Time Step.'//                   &
     &          MSGL//'Current Time Step:'//MSGF//                             &
     &          MSGL//'Critical Element ID:'//MSG1//TIMSIM%Source              &
     &          )
      ENDIF
!!
!!       B O T T O M   O F   T I M E   I N T E G R A T I O N   L O O P.
!!
      GO TO 200
!!
      END
!!_
      SUBROUTINE READ_INITIAL_CONDITIONS
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 15:30:17
!!
!! Purpose: Read a priviously written restart file and extract as initial
!! conditions the values found there. No check is made to make sure all
!! values are used and no check is made to be sure all variables of the
!! current model are initialized by data from the restart file. This method
!! provides the maximum user flexibility, but is also very dangerous.
!!
      USE shared_common_data
!!
!! The complete simulation data set.
!!
      USE indx_;          USE node_;          USE tabulated_function_;
      USE beam_;          USE coord_;         USE sliding_interface_;
      USE force_;         USE hexah_;         USE nonreflecting_bc_;
      USE penta_;         USE tetra_;         USE nodal_point_mass_;
      USE lsold_;         USE membt_;         USE constrained_node_;
      USE membq_;         USE truss_;         USE displacement_bc_;
      USE platt_;         USE platq_;         USE rigid_body_mass_;
      USE motion_;        USE stress_;        USE enumerated_sets_;
      USE spring_;        USE damper_;        USE contact_surface_;
      USE segment_;       USE tied_bc_;       USE state_variables_;
      USE results_;       USE gauge1d_;       USE rigid_wall_bc_;
      USE gauge2d_;       USE gauge3d_;       USE contact_node_;
      USE node_set_;      USE force_bc_;      USE sliding_node_;
      USE material_;      USE layering_;      USE segment_set_;
      USE massprop_;      USE spring_bc_;     USE element_set_;
      USE damper_bc_;     USE spot_weld_;     USE periodic_bc_;
      USE qa_record_;     USE nrbc_data_;     USE pressure_bc_;
      USE plate_pair_;    USE section_2d_;    USE section_1d_;
      USE rigid_body_;    USE body_force_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Define local records used to buffer-in initial data.
!!
      TYPE (gauge1d_type)     :: Ith_GAUGE1D  ! 1-D strain gauge
      TYPE (gauge2d_type)     :: Ith_GAUGE2D  ! 2-D strain gauge
      TYPE (gauge3d_type)     :: Ith_GAUGE3D  ! 3-D strain gauge
      TYPE (node_type)        :: Ith_NODE     ! Nodal point data
      TYPE (motion_type)      :: Ith_MOTION   ! Nodal point motion
      TYPE (force_type)       :: Ith_FORCE    ! Nodal point forces
      TYPE (hexah_type)       :: Ith_HEXAH    ! Solid hexahedron
      TYPE (penta_type)       :: Ith_PENTA    ! Solid pentahedron
      TYPE (tetra_type)       :: Ith_TETRA    ! Solid tetrahedron
      TYPE (lsold_type)       :: Ith_LSOLD    ! Layered solid
      TYPE (hexah_type)       :: Ith_LSHEX    ! Hexahedron for layered solid
      TYPE (membq_type)       :: Ith_MEMBQ    ! Membrane quadrilateral
      TYPE (membq_type)       :: Ith_LSMBQ    ! Membrane for layered solid
      TYPE (membt_type)       :: Ith_MEMBT    ! Membrane triangle
      TYPE (truss_type)       :: Ith_TRUSS    ! Truss (axial force only)
      TYPE (platq_type)       :: Ith_PLATQ    ! Plate quadrilateral
      TYPE (platt_type)       :: Ith_PLATT    ! Plate triangle
      TYPE (beam_type)        :: Ith_BEAM     ! Beam (axial/bending/torsion)
      TYPE (spring_type)      :: Ith_SPRING   ! Spring, translational/rotational
      TYPE (damper_type)      :: Ith_DAMPER   ! Damper, translational/rotational
      TYPE (spring_bc_type)   :: Ith_SPRING_BC! Spring BC, translational/rotatio
      TYPE (damper_bc_type)   :: Ith_DAMPER_BC! Damper BC, translational/rotatio
      TYPE (TIMESIMULATION)   :: RESTART_TIME ! Used to get TIME record
!!
      INTEGER, PARAMETER :: MAXSV=320   ! 320 = 20 x 16 beam locations
!!
      CHARACTER                                                                &
     &          File_Name*80            ! Holds full file name from INQUIRE
      INTEGER                         &
     &          ID,                   & ! Restart record value of element ID   &
     &          Ipt,                  & ! Plate integration station read       &
     &          Isv,                  & ! Starting location in state variable  &
     &          Nsv,                  & ! Number of state variables expected   &
     &          Msv                     ! Number of state variables read.
      REAL(KIND(0D0))                 &
     &          Ith_STRESS(6),        & ! Shell stresses, i-th layer           &
     &          Ith_STATE_VS(MAXSV)     ! State Variables, i-th element
      LOGICAL                                                                  &
     &          FOUND,                                                         &
     &          EXSTAT,                                                        &
     &          IOERROR,                                                       &
     &          ID_FOUND
!!
!! Inquire about file's existence and obtain full file name.
!!
      IOERROR = .TRUE.
      INQUIRE                                                                  &
     &          (                                                              &
     &          FILE        = 'fmardi',                                        &
     &          EXIST       = EXSTAT,                                          &
     &          NAME        = File_Name,                                       &
     &          ERR         = 100                                              &
     &          )
      IOERROR = .FALSE.
!!
!! Warning for failed INQUIRE operation.
!!
 100    IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_INITIAL_CONDITIONS.001.01'//                       &
     &          MSGL//'Error Executing INQUIRE On: fmardi'                     &
     &          )
!!
!! Warning for non-existant file.
!!
      ELSE IF (.NOT.EXSTAT) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_INITIAL_CONDITIONS.001.02'//                       &
     &          MSGL//'INQUIRE Unable To Find: fmardi'                         &
     &          )
      ENDIF
!!
!! Open restart data input file.
!!
!SPEC_CPU2000      IOERROR = .TRUE.
!SPEC_CPU2000      OPEN
!SPEC_CPU2000     &          (
!SPEC_CPU2000     &          UNIT   =  IO_UNIT%LRDI,
!SPEC_CPU2000     &          FILE   =  File_Name,
!SPEC_CPU2000     &          STATUS = 'OLD',
!SPEC_CPU2000     &          FORM   = 'UNFORMATTED',
!SPEC_CPU2000     &          ERR    =  200
!SPEC_CPU2000     &          )
!SPEC_CPU2000      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
!SPEC_CPU2000 200    IF (IOERROR) THEN
!SPEC_CPU2000        CALL USER_MESSAGE
!SPEC_CPU2000     &          (
!SPEC_CPU2000     &          MSGL//'FATAL'//
!SPEC_CPU2000     &          MSGL//'READ_INITIAL_CONDITIONS.002.00'//
!SPEC_CPU2000     &          MSGL//'Unable To Execute OPEN On: '//TRIM(File_Name)
!SPEC_CPU2000     &          )
!SPEC_CPU2000      ELSE
!!
!! Read title and program identification information.
!!
!SPEC_CPU2000        READ (IO_UNIT%LRDI)
!SPEC_CPU2000     &          JOB_ID_RECORD%RESTART%TITLE,
!SPEC_CPU2000     &          JOB_ID_RECORD%RESTART%DATEE,
!SPEC_CPU2000     &          JOB_ID_RECORD%RESTART%TIMEE
!!
!! Read TIME record, but skip over energy balance data, operation
!! counts, rezone data, and subprocess time accumulations.
!!
!SPEC_CPU2000        READ (IO_UNIT%LRDI) RESTART_TIME
!!
!! Read counters in restart data file.
!!
!SPEC_CPU2000        READ (IO_UNIT%LRDI)
!SPEC_CPU2000     &          IRSIF,IRSQA,IRSNP,IRSEL,IRSHX,IRSPX,IRSTX,IRSLS,IRSLX,
!SPEC_CPU2000     &          IRSLM,IRSM4,IRSM3,IRSTR,IRSP4,IRSP3,IRSBM,IRSSP,IRSDM,
!SPEC_CPU2000     &          IRSSG,IRSDC,IRSTC,IRSSW,IRSWC,IRSBF,IRSPC,IRSFC,IRSSC,
!SPEC_CPU2000     &          IRSVC,IRSCC,IRSNR,IRSND,IRSFS,IRSIT,IRSSI,IRSSN,IRSCE,
!SPEC_CPU2000     &          IRSCN,IRSCX,IRSNS,IRSNE,IRSES,IRSEE,IRSSS,IRSSE,IRSMT,
!SPEC_CPU2000     &          IRSLU,IRSTF,IRSFP,IRSRF,IRSPV,IRSS1,IRSS2,IRSG1,IRSG2,
!SPEC_CPU2000     &          IRSG3,IRSMP,IRSRM,IRSCM,IRSIC,IRSRB,IRSRT,IRSST,IRSAX,
!SPEC_CPU2000     &          IRSPP,IRSNC,IRSRE,IRSID
!SPEC_CPU2000!!
!SPEC_CPU2000!! Skip over QA data records, not used in the "changed" restart initializa
!!
!SPEC_CPU2000        DO i = 1,IRSQA
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! QA_RECORD(i)
!SPEC_CPU2000        ENDDO
!!
!! The following strain gauge data will be used to initialize the current
!! simulation.
!!
!SPEC_CPU2000        DO i = 1,IRSG1
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_GAUGE1D
!SPEC_CPU2000          ID = Ith_GAUGE1D%PAR%GauID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMG1)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. GAUGE1D(n)%PAR%GauID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            GAUGE1D(n)%RES = Ith_GAUGE1D%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSG2
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_GAUGE2D
!SPEC_CPU2000          ID = Ith_GAUGE2D%PAR%GauID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMG2)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. GAUGE2D(n)%PAR%GauID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            GAUGE2D(n)%RES = Ith_GAUGE2D%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSG3
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_GAUGE3D
!SPEC_CPU2000          ID = Ith_GAUGE3D%PAR%GauID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMG3)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. GAUGE3D(n)%PAR%GauID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            GAUGE3D(n)%RES = Ith_GAUGE3D%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!! Skip over data records not used in the "changed" restart initialization.
!!
!SPEC_CPU2000        DO i = 1,IRSMT
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! MATERIAL(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSLU
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! LAYERING(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSS2
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SECTION_2D(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSS1
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SECTION_1D(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSRB
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! RIGID_BODY(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSRM
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! RIGID_BODY_MASS(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSCM
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NODAL_POINT_MASS(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSDC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! DISPLACEMENT_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSTC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! TIED_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSSW
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SPOT_WELD(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSWC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! RIGID_WALL_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSBF
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! BODY_FORCE(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSPC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! PRESSURE_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSFC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! FORCE_BC(i)
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSSC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_SPRING_BC
!SPEC_CPU2000          ID = Ith_SPRING_BC%SprID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMSC)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. SPRING_BC(n)%SprID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            SPRING_BC(n)%RES = Ith_SPRING_BC%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSVC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_DAMPER_BC
!SPEC_CPU2000          ID = Ith_DAMPER_BC%DprID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMVC)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. DAMPER_BC(n)%DprID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            DAMPER_BC(n)%RES = Ith_DAMPER_BC%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSCC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! PERIODIC_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSNR
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NONREFLECTING_BC(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSSI
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SLIDING_INTERFACE(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSTF
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! TABULATED_FUNCTION(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSNS
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NODE_SET(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSES
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! ELEMENT_SET(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSSS
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SEGMENT_SET(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSSG
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! SEGMENT(i)
!SPEC_CPU2000        ENDDO
!!
!! The following data will be used to initialize the current simulation.
!!
!SPEC_CPU2000        DO i = 1,IRSRT
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_NODE,Ith_MOTION,Ith_FORCE
!SPEC_CPU2000          ID = Ith_NODE%ID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMNP)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. NODE(n)%ID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            IF (i .LE. IRSNP) THEN
!SPEC_CPU2000              MOTION(n)%Ux = Ith_MOTION%Ux
!SPEC_CPU2000              MOTION(n)%Uy = Ith_MOTION%Uy
!SPEC_CPU2000              MOTION(n)%Uz = Ith_MOTION%Uz
!SPEC_CPU2000              MOTION(n)%Vx = Ith_MOTION%Vx
!SPEC_CPU2000              MOTION(n)%Vy = Ith_MOTION%Vy
!SPEC_CPU2000              MOTION(n)%Vz = Ith_MOTION%Vz
!SPEC_CPU2000            ELSE IF (NODE(n)%IRT .GT. 0) THEN
!SPEC_CPU2000              n = NODE(n)%IRT
!SPEC_CPU2000              MOTION(n)%Ux = Ith_MOTION%Ux
!SPEC_CPU2000              MOTION(n)%Uy = Ith_MOTION%Uy
!SPEC_CPU2000              MOTION(n)%Uz = Ith_MOTION%Uz
!SPEC_CPU2000              MOTION(n)%Vx = Ith_MOTION%Vx
!SPEC_CPU2000              MOTION(n)%Vy = Ith_MOTION%Vy
!SPEC_CPU2000              MOTION(n)%Vz = Ith_MOTION%Vz
!SPEC_CPU2000            ENDIF
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSHX
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_HEXAH
!SPEC_CPU2000          ID = Ith_HEXAH%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMHX)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. HEXAH(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            HEXAH(n)%RES = Ith_HEXAH%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSPX
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_PENTA
!SPEC_CPU2000          ID = Ith_PENTA%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMPX)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. PENTA(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            PENTA(n)%RES = Ith_PENTA%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSTX
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_TETRA
!SPEC_CPU2000          ID = Ith_TETRA%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMTX)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. TETRA(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            TETRA(n)%RES = Ith_TETRA%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSLS
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_LSOLD
!SPEC_CPU2000          ID = Ith_LSOLD%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMLS)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. LSOLD(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            LSOLD(n)%RES = Ith_LSOLD%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSLX
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_LSHEX
!SPEC_CPU2000          ID = Ith_LSHEX%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMLX)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. LSHEX(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            LSHEX(n)%RES = Ith_LSHEX%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSLM
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_LSMBQ
!SPEC_CPU2000          ID = Ith_LSMBQ%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMLM)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. LSMBQ(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            LSMBQ(n)%RES = Ith_LSMBQ%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSM3
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_MEMBT
!SPEC_CPU2000          ID = Ith_MEMBT%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMM3)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. MEMBT(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            MEMBT(n)%RES = Ith_MEMBT%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSM4
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_MEMBQ
!SPEC_CPU2000          ID = Ith_MEMBQ%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMM4)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. MEMBQ(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            MEMBQ(n)%RES = Ith_MEMBQ%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSTR
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_TRUSS
!SPEC_CPU2000          ID = Ith_TRUSS%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMTR)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. TRUSS(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            TRUSS(n)%RES = Ith_TRUSS%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSP3
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_PLATT
!SPEC_CPU2000          ID = Ith_PLATT%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMP3)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. PLATT(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            PLATT(n)%RES = Ith_PLATT%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSP4
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_PLATQ
!SPEC_CPU2000          ID = Ith_PLATQ%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMP4)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. PLATQ(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            PLATQ(n)%RES = Ith_PLATQ%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSBM
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_BEAM
!SPEC_CPU2000          ID = Ith_BEAM%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMBM)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000!SPEC_CPU2000            FOUND = (ID .EQ. BEAM(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000           BEAM(n)%RES = Ith_BEAM%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000!SPEC_CPU2000        DO i = 1,IRSSP
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_SPRING
!SPEC_CPU2000          ID = Ith_SPRING%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMSP)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. SPRING(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            SPRING(n)%RES = Ith_SPRING%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        DO i = 1,IRSDM
!SPEC_CPU2000          READ (IO_UNIT%LRDI) Ith_DAMPER
!SPEC_CPU2000          ID = Ith_DAMPER%PAR%EleID
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMDM)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. DAMPER(n)%PAR%EleID)
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          IF (FOUND) THEN
!SPEC_CPU2000            DAMPER(n)%RES = Ith_DAMPER%RES
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000!! Skip over data records not used in the "changed" restart initialization
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSNE
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NNPSETS(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSEE
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NELSETS(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSSE
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! NSGSETS(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSPP
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! PLATE_PAIR(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000        DO i = 1,IRSNC
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! CONSTRAINED_NODE(i)
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000!! Skip over bulk write of plate stresses.
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSST
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ! (STRESS(k,i), k=1,6
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000        DO i = 1,IRSST
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ID,Ipt,Ith_STRESS
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          FOUND = .FALSE.
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMP3)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. PLATT(n)%PAR%EleID)
!SPEC_CPU2000            IF (FOUND) THEN
!SPEC_CPU2000              Ist = PLATT(n)%PAR%Ist - 1
!SPEC_CPU2000              Ipt = MIN (Ipt,Ipts_PLATT(n))
!SPEC_CPU2000              DO k = 1,6
!SPEC_CPU2000                STRESS(k,Ist+Ipt) = Ith_STRESS(k)
!SPEC_CPU2000              ENDDO
!SPEC_CPU2000            ENDIF
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000          n = 0
!SPEC_CPU2000          DO WHILE (.NOT.FOUND .AND. n .LT. NUMP4)
!SPEC_CPU2000            n = n + 1
!SPEC_CPU2000            FOUND = (ID .EQ. PLATQ(n)%PAR%EleID)
!SPEC_CPU2000            IF (FOUND) THEN
!SPEC_CPU2000              Ist = PLATQ(n)%PAR%Ist - 1
!SPEC_CPU2000              Ipt = MIN (Ipt,Ipts_PLATQ(n))
!SPEC_CPU2000              DO k = 1,6
!SPEC_CPU2000                STRESS(k,Ist+Ipt) = Ith_STRESS(k)
!SPEC_CPU2000              ENDDO
!SPEC_CPU2000            ENDIF
!SPEC_CPU2000          ENDDO
!SPEC_CPU2000        ENDDO
!SPEC_CPU2000!!
!SPEC_CPU2000!! Skip over bulk write of state variables.
!SPEC_CPU2000!!
!SPEC_CPU2000        READ (IO_UNIT%LRDI) ! (STATE_VARIABLES(i), i=1,NUMAX)
!SPEC_CPU2000!!
!! Read state variables element-by-element.
!!
!SPEC_CPU2000        NUMBER_OF_ELEMENTS_WRITTEN = IRSHX+IRSPX+IRSTX+IRSLX+IRSLM+
!SPEC_CPU2000     &  IRSM3+IRSM4+IRSTR+IRSP3+IRSP4+IRSBM+IRSSP+IRSDM+IRSSC+IRSVC
!SPEC_CPU2000
!SPEC_CPU2000        DO i = 1,NUMBER_OF_ELEMENTS_WRITTEN
!SPEC_CPU2000          READ (IO_UNIT%LRDI) ID,Msv,(Ith_STATE_VS(m), m=1,MAXSV)
!SPEC_CPU2000          Nsv = Msv
!SPEC_CPU2000          IF (ID_FOUND(ID,Isv,Nsv)) THEN
!SPEC_CPU2000            Isv = Isv - 1
!SPEC_CPU2000            DO n = 1,MIN(Msv,Nsv)
!SPEC_CPU2000              STATE_VARIABLES(Isv+n) = Ith_STATE_VS(n)
!SPEC_CPU2000            ENDDO
!SPEC_CPU2000          ENDIF
!SPEC_CPU2000        ENDDO
!!
!SPEC_CPU2000        CLOSE (UNIT=IO_UNIT%LRDI, STATUS='KEEP')
!!
!! Inform user that a restart file has been read.
!!
!SPEC_CPU2000        WRITE (MSGF,'(1PE12.4)') RESTART_TIME%Total
!SPEC_CPU2000        WRITE (MSG1,'(I8)') RESTART_TIME%Step
!SPEC_CPU2000        CALL USER_MESSAGE
!SPEC_CPU2000     &       (
!SPEC_CPU2000     &       MSGL//'INFORM'//
!SPEC_CPU2000     &       MSGL//'READ_INITIAL_CONDITIONS.003.00'//
!SPEC_CPU2000     &       MSGL//'Restart File Read: '//TRIM(File_Name)//
!SPEC_CPU2000     &       MSGL//'Restart File Title: '//JOB_ID_RECORD%RESTART%TITLE//
!SPEC_CPU2000     &       MSGL//'Restart Data Simulation Time:'//MSGF//
!SPEC_CPU2000     &       MSGL//'Restart Data Simulation Step:'//MSG1
!SPEC_CPU2000     &       )
!SPEC_CPU2000      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION ID_FOUND (ID,Isv,Nsv)
!!
!! Copyright (c) by KEY Associates; 25-APR-1992 09:57:15.48
!!
      USE shared_common_data
      USE material_
      USE section_2d_
      USE hexah_
      USE penta_
      USE tetra_
      USE lsold_
      USE membt_
      USE membq_
      USE truss_
      USE platt_
      USE platq_
      USE beam_
      USE spring_
      USE damper_
      USE spring_bc_
      USE damper_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)    :: ID
      INTEGER, INTENT(OUT)   :: Isv
      INTEGER, INTENT(INOUT) :: Nsv
!!
!! Local variables.
      LOGICAL :: FOUND
!!
      Isv = 0
      IF (Nsv .EQ. 0) THEN
        ID_FOUND = .FALSE.
      ELSE
        FOUND = .FALSE.
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMHX)
          n = n + 1
          FOUND = (ID .EQ. HEXAH(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = HEXAH(n)%PAR%Isv
            Nsv = MATERIAL(HEXAH(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMPX)
          n = n + 1
          FOUND = (ID .EQ. PENTA(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = PENTA(n)%PAR%Isv
            Nsv = MATERIAL(PENTA(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMTX)
          n = n + 1
          FOUND = (ID .EQ. TETRA(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = TETRA(n)%PAR%Isv
            Nsv = MATERIAL(TETRA(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMLX)
          n = n + 1
          FOUND = (ID .EQ. LSHEX(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = LSHEX(n)%PAR%Isv
            Nsv = MATERIAL(LSHEX(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMLM)
          n = n + 1
          FOUND = (ID .EQ. LSMBQ(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = LSMBQ(n)%PAR%Isv
            Nsv = MATERIAL(LSMBQ(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMM3)
          n = n + 1
          FOUND = (ID .EQ. MEMBT(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = MEMBT(n)%PAR%Isv
            Nsv = MATERIAL(MEMBT(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMM4)
          n = n + 1
          FOUND = (ID .EQ. MEMBQ(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = MEMBQ(n)%PAR%Isv
            Nsv = MATERIAL(MEMBQ(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMTR)
          n = n + 1
          FOUND = (ID .EQ. TRUSS(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = TRUSS(n)%PAR%Isv
            Nsv = MATERIAL(TRUSS(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMP3)
          n = n + 1
          FOUND = (ID .EQ. PLATT(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = PLATT(n)%PAR%Isv
            Nsv = MATERIAL(PLATT(n)%PAR%MatID)%Nsv
            Nsv = Nsv * Ipts_PLATT(n)
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMP4)
          n = n + 1
          FOUND = (ID .EQ. PLATQ(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = PLATQ(n)%PAR%Isv
            Nsv = MATERIAL(PLATQ(n)%PAR%MatID)%Nsv
            Nsv = Nsv * Ipts_PLATQ(n)
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMBM)
          n = n + 1
          FOUND = (ID .EQ. BEAM(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = BEAM(n)%PAR%Isv
            Nsv = 16 * MATERIAL(BEAM(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMSP)
          n = n + 1
          FOUND = (ID .EQ. SPRING(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = SPRING(n)%PAR%Isv
            Nsv = MATERIAL(SPRING(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMDM)
          n = n + 1
          FOUND = (ID .EQ. DAMPER(n)%PAR%EleID)
          IF (FOUND) THEN
            Isv = DAMPER(n)%PAR%Isv
            Nsv = MATERIAL(DAMPER(n)%PAR%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMSC)
          n = n + 1
          FOUND = (ID .EQ. SPRING_BC(n)%SprID)
          IF (FOUND) THEN
            Isv = SPRING_BC(n)%Isv
            Nsv = MATERIAL(SPRING_BC(n)%MatID)%Nsv
          ENDIF
        ENDDO
        n = 0
        DO WHILE (.NOT.FOUND .AND. n .LT. NUMVC)
          n = n + 1
          FOUND = (ID .EQ. DAMPER_BC(n)%DprID)
          IF (FOUND) THEN
            Isv = DAMPER_BC(n)%Isv
            Nsv = MATERIAL(DAMPER_BC(n)%MatID)%Nsv
          ENDIF
        ENDDO
!!
        ID_FOUND = FOUND
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE ELEMENT_INITIALIZATION
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 20:51:22.00
!!
!! Purpose: Compute divergence of initial stress state, compute
!! terms in mass matrix, and find sound speed.
!!
      USE shared_common_data
!!
!! INTERNAL FORCES
!! Look for zero or negative element volumes/areas/lengths.
!!
      ERROR%COUNT = 0
!!
!! Initialize 8-node hexahedron solid elements.
!!
      IF (NUMHX .GT. 0) CALL HEXAH_INITIALIZATION
!!
!! Initialize 6-node pentahedron solid elements.
!!
      IF (NUMPX .GT. 0) CALL PENTA_INITIALIZATION
!!
!! Initialize 4-node tetrahedron solid elements.
!!
      IF (NUMTX .GT. 0) CALL TETRA_INITIALIZATION
!!
!! Initialize 8-node hexahedron layered solid elements.
!!
      IF (NUMLS .GT. 0) CALL LSOLD_INITIALIZATION
!!
!! Initialize 3-node membrane elements.
!!
      IF (NUMM3 .GT. 0) CALL MEMBT_INITIALIZATION
!!
!! Initialize 4-node membrane elements.
!!
      IF (NUMM4 .GT. 0) CALL MEMBQ_INITIALIZATION
!!
!! Initialize 2-node truss elements.
!!
      IF (NUMTR .GT. 0) CALL TRUSS_INITIALIZATION
!!
!! Initialize 3-node plate elements.
!!
      IF (NUMP3 .GT. 0) CALL PLATT_INITIALIZATION
!!
!! Initialize 4-node plate elements.
!!
      IF (NUMP4 .GT. 0) CALL PLATQ_INITIALIZATION
!!
!! Initialize 2-node beam elements.
!!
      IF (NUMBM .GT. 0) CALL BEAM_INITIALIZATION
!!
!! Initialize 2-node spring elements.
!!
      IF (NUMSP .GT. 0) CALL SPRING_INITIALIZATION
!!
!! Initialize 2-node damper elements.
!!
      IF (NUMDM .GT. 0) CALL DAMPER_INITIALIZATION
!!
!! Check for non-zero count of zero or negative element volumes.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &      (                                                                  &
     &      MSGL//'FATAL'//                                                    &
     &      MSGL//'ELEMENT_INITIALIZATION.001.00'//                            &
     &      MSGL//'Total Number Of Zero/Negative Element Volumes:'//MSG1       &
     &      )
      ENDIF
!!
!! Scatter element internal forces to nodes.
!!
      CALL SCATTER_ELEMENT_NODAL_FORCES
!!
      RETURN
      END
!!_
      SUBROUTINE INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 20:51:30.00
!!
!! Purpose: Compute divergence of current stress state and find sound
!! speeds for next time step.
!!
      USE shared_common_data
!!
!! INTERNAL FORCES
!! Calculate forces from 8-node hexahedron solid elements.
!!
      IF (NUMHX .GT. 0) THEN
        CALL HEXAH_INTERNAL_FORCES
!!
!! Record time spent in solid elements.
!!
!SPEC_CPU2000        CALL TIMER (10)
      ENDIF
!!
!! Calculate forces from 6-node pentahedron solid elements.
!!
      IF (NUMPX .GT. 0) THEN
        CALL PENTA_INTERNAL_FORCES
!!
!! Record time spent in solid elements.
!!
!SPEC_CPU2000        CALL TIMER (11)
      ENDIF
!!
!! Calculate forces from 4-node tetrahedron solid elements.
!!
      IF (NUMTX .GT. 0) THEN
        CALL TETRA_INTERNAL_FORCES
!!
!! Record time spent in solid elements.
!!
!SPEC_CPU2000        CALL TIMER (12)
      ENDIF
!!
!! Calculate forces from 8-node hexahedron layered solid elements.
!!
      IF (NUMLS .GT. 0) THEN
        CALL LSOLD_INTERNAL_FORCES
!!
!! Record time spent in layered solid elements.
!!
!SPEC_CPU2000        CALL TIMER (13)
      ENDIF
!!
!! Calculate forces from 3-node membrane elements.
!!
      IF (NUMM3 .GT. 0) THEN
        CALL MEMBT_INTERNAL_FORCES
!!
!! Record time spent in 3-node membrane elements.
!!
!SPEC_CPU2000        CALL TIMER (14)
      ENDIF
!!
!! Calculate forces from 4-node membrane elements.
!!
      IF (NUMM4 .GT. 0) THEN
        CALL MEMBQ_INTERNAL_FORCES
!!
!! Record time spent in 4-node membrane elements.
!!
!SPEC_CPU2000        CALL TIMER (15)
      ENDIF
!!
!! Calculate forces from 2-node truss elements.
!!
      IF (NUMTR .GT. 0) THEN
        CALL TRUSS_INTERNAL_FORCES
!!
!! Record time spent in 2-node truss elements.
!!
!SPEC_CPU2000        CALL TIMER (16)
      ENDIF
!!
!! Calculate forces from 3-node plate elements.
!!
      IF (NUMP3 .GT. 0) THEN
        CALL PLATT_INTERNAL_FORCES
!!
!! Record time spent in 3-node plate elements.
!!
!SPEC_CPU2000        CALL TIMER (17)
      ENDIF
!!
!! Calculate forces from 4-node plate elements.
!!
      IF (NUMP4 .GT. 0) THEN
        CALL PLATQ_INTERNAL_FORCES
!!
!! Record time spent in 4-node plate elements.
!!
!SPEC_CPU2000        CALL TIMER (18)
      ENDIF
!!
!! Calculate forces from 2-node beam elements.
!!
      IF (NUMBM .GT. 0) THEN
        CALL BEAM_INTERNAL_FORCES
!!
!! Record time spent in 2-node beam elements.
!!
!SPEC_CPU2000        CALL TIMER (19)
      ENDIF
!!
!! Calculate forces from 2-node spring elements.
!!
      IF (NUMSP .GT. 0) THEN
        CALL SPRING_INTERNAL_FORCES
!!
!! Record time spent in 2-node spring elements.
!!
!SPEC_CPU2000        CALL TIMER (20)
      ENDIF
!!
!! Calculate forces from 2-node damper elements.
!!
      IF (NUMDM .GT. 0) THEN
        CALL DAMPER_INTERNAL_FORCES
!!
!! Record time spent in 2-node damper elements.
!!
!SPEC_CPU2000        CALL TIMER (21)
      ENDIF
!!
!! Scatter element internal forces to nodes.
!!
      CALL SCATTER_ELEMENT_NODAL_FORCES
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000      CALL TIMER (3)
!!
      RETURN
      END
!!_
      SUBROUTINE EXTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 20:51:36.00
!!
!! Purpose: Evaluate external forcing functions and active boundary
!! conditions acting of the nodes (and element "edges").
!!
      USE shared_common_data
!!
!! EXTERNAL FORCES.
!! Compute body forces acting on nodes.
!!
      IF (NUMBF .GT. 0 .AND. NUMTF .NE. 0) CALL IMPOSE_BODY_FORCE
!!
!! Compute surface tractions.
!!
      IF (NUMPC .NE. 0 .AND. NUMTF .NE. 0) CALL IMPOSE_PRESSURE_BC
!!
!! Compute concentrated forces.
!!
      IF (NUMFC .NE. 0 .AND. NUMTF .NE. 0) CALL IMPOSE_FORCE_BC
!!
!! ACTIVE BOUNDARY CONDITIONS
!! Compute boundary spring forces.
!!
      IF (NUMSC .NE. 0) CALL IMPOSE_SPRING_BC
!!
!! Compute boundary damper forces.
!!
      IF (NUMVC .NE. 0) CALL IMPOSE_DAMPER_BC
!!
!! Record time spent in applying traction boundary conditions.
!!
!SPEC_CPU2000      CALL TIMER (9)
!!
!! Non-reflecting boundary conditions. (This is an external force on the
!! mesh. However, it is an "active" boundary condition and requires the
!! current element density and sound speed to compute the appropriate
!! boundary tractions. Thus, it must be called after the internal forces
!! have been evaluated.)
!!
      IF (NUMNR .GT. 0) THEN
        CALL IMPOSE_NONREFLECTING_BC
!!
!! Record time spent in applying non-reflecting boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (22)
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE KINEMATIC_CONSTRAINTS
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 21:36:06.00
!!
!! Purpose: Impose kinematic constraints on nodal point motion. Note
!! that the order in which the kinematic boundary conditions are applied
!! is important in the event that a given degree of freedom at a given
!! nodal point has two or more kinematic constraints specified. The
!! last kinematic constraint controls the node's motion. Here, rigid
!! body motion is considered primary and is thus applied last. The
!! order is as follows:
!!
!!      1. Mid-Side Constraint 1. Introduce influence of mid-side
!!         nodal point on motion of bounding nodal points.
!!
!!      2. Sliding Interface Contact. Look ahead to check for
!!         contacts; introduce forces to mitigate mesh overlap.
!!
!!      3. Tied Boundary Condtions. For nodal points that have a
!!         "common" motion, make accelerations identical.
!!
!!      4. Spot Weld Constraints. For nodes or elements "connected"
!!         by spot welds, make their motion identical.
!!
!!      5. Displacement Boundary Conditions. For individual nodal
!!         points with specified motion, modify accelerations to
!!         match the required motion.
!!
!!      6. Periodic Boundary Conditions. For node sets representing
!!         the cyclically repeating surfaces, make the motions match.
!!
!!      7. Rigid Wall Boundary Condition. For those nodes contacting
!!         a rigid wall, modify the motion to preclude penetration.
!!
!!      8. Mid-Side Constraint 2. Make the motion of the mid-side
!!         nodal point conform to the bounding nodal points.
!!
!!      9. Rigid Body Motion. For those nodal points comprising a
!!         rigid body, modify the accelerations to obtain rigid
!!         body motion.
!!
      USE shared_common_data
      USE node_
      USE force_
      USE motion_
!!
!! Modify the external forces and accelerations of nodal points
!! influenced by midside constraints.
!!
      IF (NUMNC .GT. 0) THEN
        CALL IMPOSE_MIDSIDE_CONSTRAINTS1
!!
!! Record time spent in computing midside nodal point constraints.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Check for sliding interface contacts. This is a "look-ahead" calculation.
!! If a nodal point in the next time step crosses into the contact-element
!! which it is currently opposite, reaction forces are added which will cause
!! the nodal point in question to end up located on the contact element. If
!! the nodal point leaves the contact element, no action is taken.
!!
      IF (NUMSI .NE. 0) THEN
        CALL SLIDING_INTERFACE_CONTACT
!!
!! Recompute accelerations to take into account any sliding interface
!! interactions which may have occurred.
!!
        DO N = 1,NUMRT
          MOTION(N)%Ax = NODE(N)%Minv * (FORCE(N)%Xext-FORCE(N)%Xint)
          MOTION(N)%Ay = NODE(N)%Minv * (FORCE(N)%Yext-FORCE(N)%Yint)
          MOTION(N)%Az = NODE(N)%Minv * (FORCE(N)%Zext-FORCE(N)%Zint)
        ENDDO
!!
!! Recompute the acceleration of nodal points with concentrated inertia.
!!
        IF (NUMCM .GT. 0) CALL UPDATE1_CONCENTRATED_MASSES
!!
!! Record time spent in sliding interface calculations.
!!
!SPEC_CPU2000        CALL TIMER (4)
      ENDIF
!!
!! Apply tied boundary conditions.
!!
      IF (NUMTC .GT. 0) THEN
        CALL IMPOSE_TIED_BC
!!
!! Record time spent in applying kinematic boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Apply spot weld constraints.
!!
      IF (NUMSW .GT. 0) THEN
        CALL IMPOSE_SPOT_WELD
!!
!! Record time spent in applying kinematic boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Apply displacement boundary conditions.
!!
      IF (NUMDC .GT. 0) THEN
        CALL IMPOSE_DISPLACEMENT_BC
!!
!! Record time spent in applying kinematic boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Apply periodic boundary conditions.
!!
      IF (NUMCC .GT. 0) THEN
        CALL IMPOSE_PERIODIC_BC
!!
!! Record time spent in applying kinematic boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Apply rigid wall boundary conditions. This is a "look-ahead" calculation.
!! If a nodal point in the next time step crosses into the rigid wall exclu-
!! sion domain, reaction forces are added which will cause the nodal point in
!! question to end up located on the rigid wall. If the nodal point "leaves"
!! the rigid wall, no action is taken.
!!
      IF (NUMWC .NE. 0) THEN
        CALL IMPOSE_RIGID_WALL_BC
!!
!! Record time spent in applying rigid wall boundary conditions.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Insure that the motion of mid-side nodal points conforms to the motion
!! of the bounding nodes.
!!
      IF (NUMNC .GT. 0) THEN
        CALL IMPOSE_MIDSIDE_CONSTRAINTS2
!!
!! Record time spent in computing mid-side nodal point constraints.
!!
!SPEC_CPU2000        CALL TIMER (5)
      ENDIF
!!
!! Compute motion of "rigid body" regions.
!!
      IF (NUMRB .GT. 0) THEN
        CALL IMPOSE_RIGID_BODY_CONSTRAINTS
!!
!! Record time spent in computing rigid body motion.
!!
!SPEC_CPU2000        CALL TIMER (6)
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE OUTPUT_PROCESSING
!!
!! Copyright (c) by KEY Associates, 9-NOV-1991 10:51:29
!!
!! Purpose: Process all output requests during main time integration loop.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE results_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Local variables.
      CHARACTER(8)  :: SRQTYP  ! Status_ReQuest_TYPe
      INTEGER, SAVE :: NEBHDG  ! Energy balance heading flag
      INTEGER, SAVE :: SRQINT  ! Status request interval
      INTEGER, SAVE :: EBLINT  ! Energy balance interval
      LOGICAL, SAVE :: STATUS_REQUESTED ! Flag indicating a *.SRQ file exists

      LOGICAL, SAVE :: FIRST = .TRUE.

      IF (FIRST) THEN
        NEBHDG = 1
        SRQINT = MAX (1,NINT (PARAMVALUE%STATUS ))
        EBLINT = MAX (1,NINT (PARAMVALUE%ENG_BAL))
        STATUS_REQUESTED = .FALSE.
        FIRST = .FALSE.
      ENDIF
!!
      ITSTEP = TIMSIM%Step
!!
!! OUTPUT REQUEST 1.
!! Check for existence of a user status request.
!!
      IF (MOD(ITSTEP,SRQINT) .EQ. 0) THEN
        INQUIRE (FILE='fmasri', EXIST=STATUS_REQUESTED)
      ENDIF
!!
!! Process user status request.
!!
      IF (STATUS_REQUESTED) THEN
!!
!! Open status request file (*.SRI) and read keyword.
!!
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSRI,                                        &
     &          FILE   = 'fmasri',                                             &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  100                                                  &
     &          )
        READ (IO_UNIT%LSRI,'(A)') SRQTYP
        CLOSE (UNIT=IO_UNIT%LSRI, STATUS='KEEP')
        CALL CUPPER (SRQTYP)
!!
!! Open status request return message file (*.SRO) for reply.
!!
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSRO,                                        &
     &          FILE   = 'fmasro',                                             &
     &          STATUS = 'UNKNOWN',                                            &
     &          FORM   = 'FORMATTED',                                          &
     &          POSITION   = 'APPEND',                                         &
     &          ERR    =  100                                                  &
     &          )
!!
!! Echo status request keyword and progress to this point.
!!
        WRITE (IO_UNIT%LSRO,'(A)')                                             &
     &          ' FMA-3D> Status Request Keyword: '//SRQTYP
        WRITE (IO_UNIT%LSRO,'(A,F5.2,A)')                                      &
     &          ' FMA-3D> Percent Complete:       ',                           &
     &          (TIMSIM%Total/TIMSIM%Stop)*100.0D+0,'%'
!!
!! Interogate status request type SRQTYP.
!!
        IF (INDEX(SRQTYP,'PROGRESS') .NE. 0) THEN
        ELSE IF (INDEX(SRQTYP,'KILL') .NE. 0) THEN
          TIMSIM%Stop = TIMSIM%Total - TIMSIM%DTnext
          WRITE (IO_UNIT%LSRO,'(A)')                                           &
     &      ' FMA-3D> Calculation Will Be Terminated'
          WRITE (IO_UNIT%LSRO,'(A)') ' FMA-3D> At Next Opportunity.'
        ELSE
          WRITE (IO_UNIT%LSRO,'(A)') ' FMA-3D> Unrecognized Keyword.'
        ENDIF

        CLOSE (UNIT=IO_UNIT%LSRO, STATUS='KEEP')
!!
!! Error exit in the event the OPEN's fail.
!!
 100    CONTINUE
        STATUS_REQUESTED = .FALSE.
      ENDIF
!!
!! OUTPUT REQUEST 2.
!! Print current time, last time step, and energy balance.
!!
      IF (MOD(ITSTEP,EBLINT) .EQ. 0) THEN
        IF (NEBHDG .NE. 0) THEN
          WRITE (IO_UNIT%LELO,420)
          NEBHDG = 0
        ENDIF
        ENGK = 0.0
        DO N = 1,NUMRT
          ENGK = ENGK + NODE(N)%Mass *                                         &
     &        (                                                                &
     &        MOTION(N)%Vx * MOTION(N)%Vx +                                    &
     &        MOTION(N)%Vy * MOTION(N)%Vy +                                    &
     &        MOTION(N)%Vz * MOTION(N)%Vz                                      &
     &        )
        ENDDO
        ENERGY%Kinetic = 0.5D+0 * ENGK
        IF (ENERGY%External .EQ. 0.0) THEN
          ENERGY%Balance = 0.0
        ELSE
          ENERGY%Balance = 100.0D+0 * (ENERGY%External-ENERGY%Kinetic-         &
     &      ENERGY%Internal)/ENERGY%External
        ENDIF
        WRITE (IO_UNIT%LELO,440)                                               &
     &          TIMSIM%Step,                                                   &
     &          TIMSIM%Total,                                                  &
     &          TIMSIM%DTlast,                                                 &
     &          TIMSIM%DTmin,                                                  &
     &          TIMSIM%Ncrit,                                                  &
     &          TIMSIM%DTmax,                                                  &
     &          TIMSIM%Source,                                                 &
     &          ENERGY%External,                                               &
     &          ENERGY%Kinetic,                                                &
     &          ENERGY%Internal,                                               &
     &          ENERGY%Balance

      ENDIF
!!
!! OUTPUT REQUEST 3.
!! Write restart file.
!!
      IF (CONTROL%WRSTAR .GT. 0) THEN
        IF (MOD(TIMSIM%Step,CONTROL%WRSTAR) .EQ. 0) THEN
          CALL WRITE_RESTART_DATA
        ENDIF
      ENDIF
!!
!! Record time spent in solution management.
!!
!SPEC_CPU2000      CALL TIMER (3)
!!
!! OUTPUT REQUEST 5.
!! Write Plotting_Database files.
!!
      IF (NUMRF .GT. 0) THEN
        DO i = 1,NUMRF
          IF (TIMSIM%Total .GE. RESULTS(i)%Time) THEN
!!
!! Reset Next-Time-to-Write-Plotting_Database value.
!!
            IF (TIMSIM%Total .GT. RESULTS(i)%End) THEN
              RESULTS(i)%Time = 2.0 * TIMSIM%Stop
            ELSE IF (RESULTS(i)%Delta .GT. 0.0) THEN
              DO WHILE (RESULTS(i)%Time .LE. TIMSIM%Total)
                RESULTS(i)%Time = RESULTS(i)%Time + RESULTS(i)%Delta
              ENDDO
            ENDIF
!!
!! Mass property calculations.
!!
            IF (NUMMP .GT. 0) CALL MASS_PROPERTY_CALCULATION

            CALL WRITE_TO_PLOTTING_DATABASE

          ENDIF
        ENDDO
!!
!! Record time spent in creating the Plotting_Database files.
!!
!SPEC_CPU2000        CALL TIMER (7)
      ENDIF
!!
!! OUTPUT REQUEST 6.
!! Print accelerations,velocities,displacements and stresses.
!!
      IF (TIMSIM%Total .GE. PRINT%Time) THEN
!!
!! Reset Next-Time-to-Print value.
!!
        IF (TIMSIM%Total .GE. PRINT%End) THEN
          PRINT%Time = 2.0 * TIMSIM%Stop
        ELSE IF (PRINT%Delta .GT. 0.0) THEN
          DO WHILE (PRINT%Time .LE. TIMSIM%Total)
            PRINT%Time = PRINT%Time + PRINT%Delta
          ENDDO
        ENDIF
!!
!! Mass property calculations.
!!
        IF (NUMMP .GT. 0) CALL MASS_PROPERTY_CALCULATION
!!
        CALL PRINTER_OUTPUT
!!
!! Record time spent creating printed output.
!!
!SPEC_CPU2000        CALL TIMER (8)
      ENDIF
!!
!! FORMAT statements.
!!
 420    FORMAT                                                                 &
     &    (                                                                    &
     &    //' ',2X,'Cycle No.   Time',5X,'Dt_inc',4X,'Dt_min',                 &
     &    2X,'Crit_El  Dt_max',6X,'Critical Element Type',7X,                  &
     &    'Ext_Eng',5X,'K.E.',4X,'Int_Eng',2X,'% Error'                        &
     &    )
 440    FORMAT                                                                 &
     &    (                                                                    &
     &    ' ',I8,2X,3(1PE10.2),I8,1X,E10.2,1X,A30,E9.2,2E10.2,1X,0PF8.2        &
     &    )
!!
      RETURN
      END
!!_
      SUBROUTINE WRITE_RESTART_DATA
!!
!! Copyright (c) by KEY Associates, 19-FEB-1992 19:33:13
!!
!! Purpose: Write all of the data needed to restart the analysis.
!!
      USE shared_common_data
!!
!! The complete simulation data set.
!!
      USE indx_;          USE node_;          USE tabulated_function_;
      USE beam_;          USE coord_;         USE sliding_interface_;
      USE force_;         USE hexah_;         USE nonreflecting_bc_;
      USE penta_;         USE tetra_;         USE nodal_point_mass_;
      USE lsold_;         USE membt_;         USE constrained_node_;
      USE membq_;         USE truss_;         USE displacement_bc_;
      USE platt_;         USE platq_;         USE rigid_body_mass_;
      USE motion_;        USE stress_;        USE enumerated_sets_;
      USE spring_;        USE damper_;        USE contact_surface_;
      USE segment_;       USE tied_bc_;       USE state_variables_;
      USE results_;       USE gauge1d_;       USE rigid_wall_bc_;
      USE gauge2d_;       USE gauge3d_;       USE contact_node_;
      USE node_set_;      USE force_bc_;      USE sliding_node_;
      USE material_;      USE layering_;      USE segment_set_;
      USE massprop_;      USE spring_bc_;     USE element_set_;
      USE damper_bc_;     USE spot_weld_;     USE periodic_bc_;
      USE qa_record_;     USE nrbc_data_;     USE pressure_bc_;
      USE plate_pair_;    USE section_2d_;    USE section_1d_;
      USE rigid_body_;    USE body_force_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, PARAMETER :: MAXSV = 320 ! 320 = 16 x 20, Beam SV's
!!
      CHARACTER                                                                &
     &          CRDO*4                  ! Buffer for unique file extension
      LOGICAL                                                                  &
     &          IOERROR
      INTEGER                                                                  &
     &          IRDO                    ! Counter for unique file extension

!!
!! Build unique extension for next restart file: fmardo.nnn
!!
      IRDO = JOB_ID_RECORD%CURRENT%SEQUENCE_NUMBER + 1
      WRITE (CRDO,'(I4)') IRDO
      CRDO(1:1) = '.'
      JOB_ID_RECORD%CURRENT%SEQUENCE_NUMBER = IRDO
!!
!! Open restart data output file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LRDO,                                        &
     &          FILE   = 'fmardo'//CRDO,                                       &
     &          STATUS = 'UNKNOWN',                                            &
     &          FORM   = 'UNFORMATTED',                                        &
     &          ERR    =  100                                                  &
     &          )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 100    IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'WRITE_RESTART_DATA.001.01'//                            &
     &          MSGL//'Unable To Execute OPEN On: '//'fmardo'//CRDO            &
     &          )
      ELSE
!!
!! Rewind restart file; put new data in place of old.
!!
        REWIND (IO_UNIT%LRDO)
!!
!! Write title and program identification information.
!!
        WRITE (IO_UNIT%LRDO)                                                   &
     &          JOB_ID_RECORD%CURRENT%TITLE,                                   &
     &          JOB_ID_RECORD%CURRENT%DATEE,                                   &
     &          JOB_ID_RECORD%CURRENT%TIMEE,                                   &
     &          JOB_ID_RECORD%CURRENT%SEQUENCE_NUMBER,                         &
     &          JOB_ID_RECORD%PROGRAM,                                         &
     &          JOB_ID_RECORD%VERSION
!!
!! Write time and energy balance data.
!!
!SPEC_CPU2000        WRITE (IO_UNIT%LRDO)
!SPEC_CPU2000     &          ENERGY,COUNTER,REZONE,SUBPROCESS_TIME
!!
!! Write (all) counters to restart data file.
!!
        WRITE (IO_UNIT%LRDO)                                                   &
     &          NUMIF,NUMQA,NUMNP,NUMEL,NUMHX,NUMPX,NUMTX,NUMLS,NUMLX,         &
     &          NUMLM,NUMM4,NUMM3,NUMTR,NUMP4,NUMP3,NUMBM,NUMSP,NUMDM,         &
     &          NUMSG,NUMDC,NUMTC,NUMSW,NUMWC,NUMBF,NUMPC,NUMFC,NUMSC,         &
     &          NUMVC,NUMCC,NUMNR,NUMND,NUMFS,NUMIT,NUMSI,NUMSN,NUMCE,         &
     &          NUMCN,NUMCX,NUMNS,NUMNE,NUMES,NUMEE,NUMSS,NUMSE,NUMMT,         &
     &          NUMLU,NUMTF,NUMFP,NUMRF,NUMPV,NUMS1,NUMS2,NUMG1,NUMG2,         &
     &          NUMG3,NUMMP,NUMRM,NUMCM,NUMIC,NUMRB,NUMRT,NUMST,NUMAX,         &
     &          NUMPP,NUMNC,MAXRE,MXEID
!!
!! Write all data structures needed to restart simulation.
!!
        DO i = 1,NUMQA
          WRITE (IO_UNIT%LRDO) QA_RECORD(i)
        ENDDO
        DO i = 1,NUMG1
          WRITE (IO_UNIT%LRDO) GAUGE1D(i)
        ENDDO
        DO i = 1,NUMG2
          WRITE (IO_UNIT%LRDO) GAUGE2D(i)
        ENDDO
        DO i = 1,NUMG3
          WRITE (IO_UNIT%LRDO) GAUGE3D(i)
        ENDDO
        DO i = 1,NUMMT
          WRITE (IO_UNIT%LRDO) MATERIAL(i)
        ENDDO
        DO i = 1,NUMLU
          WRITE (IO_UNIT%LRDO) LAYERING(i)
        ENDDO
        DO i = 1,NUMS2
          WRITE (IO_UNIT%LRDO) SECTION_2D(i)
        ENDDO
        DO i = 1,NUMS1
          WRITE (IO_UNIT%LRDO) SECTION_1D(i)
        ENDDO
        DO i = 1,NUMRB
          WRITE (IO_UNIT%LRDO) RIGID_BODY(i)
        ENDDO
        DO i = 1,NUMRM
          WRITE (IO_UNIT%LRDO) RIGID_BODY_MASS(i)
        ENDDO
        DO i = 1,NUMCM
          WRITE (IO_UNIT%LRDO) NODAL_POINT_MASS(i)
        ENDDO
        DO i = 1,NUMDC
          WRITE (IO_UNIT%LRDO) DISPLACEMENT_BC(i)
        ENDDO
        DO i = 1,NUMTC
          WRITE (IO_UNIT%LRDO) TIED_BC(i)
        ENDDO
        DO i = 1,NUMSW
          WRITE (IO_UNIT%LRDO) SPOT_WELD(i)
        ENDDO
        DO i = 1,NUMWC
          WRITE (IO_UNIT%LRDO) RIGID_WALL_BC(i)
        ENDDO
        DO i = 1,NUMBF
          WRITE (IO_UNIT%LRDO) BODY_FORCE(i)
        ENDDO
        DO i = 1,NUMPC
          WRITE (IO_UNIT%LRDO) PRESSURE_BC(i)
        ENDDO
        DO i = 1,NUMFC
          WRITE (IO_UNIT%LRDO) FORCE_BC(i)
        ENDDO
        DO i = 1,NUMSC
          WRITE (IO_UNIT%LRDO) SPRING_BC(i)
        ENDDO
        DO i = 1,NUMVC
          WRITE (IO_UNIT%LRDO) DAMPER_BC(i)
        ENDDO
        DO i = 1,NUMCC
          WRITE (IO_UNIT%LRDO) PERIODIC_BC(i)
        ENDDO
        DO i = 1,NUMNR
          WRITE (IO_UNIT%LRDO) NONREFLECTING_BC(i)
        ENDDO
        DO i = 1,NUMSI
          WRITE (IO_UNIT%LRDO) SLIDING_INTERFACE(i)
        ENDDO
        DO i = 1,NUMTF
          WRITE (IO_UNIT%LRDO) TABULATED_FUNCTION(i)
        ENDDO
        DO i = 1,NUMNS
          WRITE (IO_UNIT%LRDO) NODE_SET(i)
        ENDDO
        DO i = 1,NUMES
          WRITE (IO_UNIT%LRDO) ELEMENT_SET(i)
        ENDDO
        DO i = 1,NUMSS
          WRITE (IO_UNIT%LRDO) SEGMENT_SET(i)
        ENDDO
        DO i = 1,NUMSG
          WRITE (IO_UNIT%LRDO) SEGMENT(i)
        ENDDO
        DO i = 1,NUMRT
          WRITE (IO_UNIT%LRDO) NODE(i),MOTION(i),FORCE(i)
        ENDDO
        DO i = 1,NUMHX
          WRITE (IO_UNIT%LRDO) HEXAH(i)
        ENDDO
        DO i = 1,NUMPX
          WRITE (IO_UNIT%LRDO) PENTA(i)
        ENDDO
        DO i = 1,NUMTX
          WRITE (IO_UNIT%LRDO) TETRA(i)
        ENDDO
        DO i = 1,NUMLS
          WRITE (IO_UNIT%LRDO) LSOLD(i)
        ENDDO
        DO i = 1,NUMLX
          WRITE (IO_UNIT%LRDO) LSHEX(i)
        ENDDO
        DO i = 1,NUMLM
          WRITE (IO_UNIT%LRDO) LSMBQ(i)
        ENDDO
        DO i = 1,NUMM3
          WRITE (IO_UNIT%LRDO) MEMBT(i)
        ENDDO
        DO i = 1,NUMM4
          WRITE (IO_UNIT%LRDO) MEMBQ(i)
        ENDDO
        DO i = 1,NUMTR
          WRITE (IO_UNIT%LRDO) TRUSS(i)
        ENDDO
        DO i = 1,NUMP3
          WRITE (IO_UNIT%LRDO) PLATT(i)
        ENDDO
        DO i = 1,NUMP4
          WRITE (IO_UNIT%LRDO) PLATQ(i)
        ENDDO
        DO i = 1,NUMBM
          WRITE (IO_UNIT%LRDO) BEAM(i)
        ENDDO
        DO i = 1,NUMSP
          WRITE (IO_UNIT%LRDO) SPRING(i)
        ENDDO
        DO i = 1,NUMDM
          WRITE (IO_UNIT%LRDO) DAMPER(i)
        ENDDO
        DO i = 1,NUMNE
          WRITE (IO_UNIT%LRDO) NNPSETS(i)
        ENDDO
        DO i = 1,NUMEE
          WRITE (IO_UNIT%LRDO) NELSETS(i)
        ENDDO
        DO i = 1,NUMSE
          WRITE (IO_UNIT%LRDO) NSGSETS(i)
        ENDDO
        DO i = 1,NUMPP
          WRITE (IO_UNIT%LRDO) PLATE_PAIR(i)
        ENDDO
        DO i = 1,NUMNC
          WRITE (IO_UNIT%LRDO) CONSTRAINED_NODE(i)
        ENDDO
!!
!! Write plate stresses without reference to which element they are related.
!!
        DO i = 1,NUMST
          WRITE (IO_UNIT%LRDO) (STRESS(k,i), k=1,6)
        ENDDO
!!
!! This second, laborious approach to writing the plate stresses elemen-by-
!! element is done solely for the purposes of being able to do a "changed"
!! restart.
!!
        DO N = 1,NUMP3
          NPL = N
          Ist = PLATT(N)%PAR%Ist - 1
          DO Ipt = 1,Ipts_PLATT(NPL)
            Ist = Ist + 1
            WRITE (IO_UNIT%LRDO)                                               &
     &        PLATT(N)%PAR%EleID,Ipt,(STRESS(k,Ist),k=1,6)
          ENDDO
        ENDDO
        DO N = 1,NUMP4
          NPL = N
          Ist = PLATQ(N)%PAR%Ist - 1
          DO Ipt = 1,Ipts_PLATQ(NPL)
            Ist = Ist + 1
            WRITE (IO_UNIT%LRDO)                                               &
     &        PLATQ(N)%PAR%EleID,Ipt,(STRESS(k,Ist),k=1,6)
          ENDDO
        ENDDO
!!
!! Write state_variables "in bulk" for use with an "unchanged" restart.
!!
        WRITE (IO_UNIT%LRDO) (STATE_VARIABLES(i), i=1,NUMAX)
!!
!! This second, laborious approach to writing the state variables element-by-
!! element is done solely for the purposes of being able to do a "changed"
!! restart. STATE_VAR_BUFFER returns a "large" fixed-length buffer, The state
!! variables, if any, occupy the first Nsv locations.
!!
        DO N = 1,NUMHX
          Isv = HEXAH(N)%PAR%Isv - 1
          Nsv = MATERIAL(HEXAH(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) HEXAH(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMPX
          Isv = PENTA(N)%PAR%Isv - 1
          Nsv = MATERIAL(PENTA(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) PENTA(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMTX
          Isv = TETRA(N)%PAR%Isv - 1
          Nsv = MATERIAL(TETRA(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) TETRA(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        DO N = 1,NUMLX
          Isv = LSHEX(N)%PAR%Isv - 1
          Nsv = MATERIAL(LSHEX(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) LSHEX(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMLM
          Isv = LSMBQ(N)%PAR%Isv - 1
          Nsv = MATERIAL(LSMBQ(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) LSMBQ(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        DO N = 1,NUMM3
          Isv = MEMBT(N)%PAR%Isv - 1
          Nsv = MATERIAL(MEMBT(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) MEMBT(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMM4
          Isv = MEMBQ(N)%PAR%Isv - 1
          Nsv = MATERIAL(MEMBQ(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) MEMBQ(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMTR
          Isv = TRUSS(N)%PAR%Isv - 1
          Nsv = MATERIAL(TRUSS(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) TRUSS(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        DO N = 1,NUMP3
          NPL = N
          Isv = PLATT(N)%PAR%Isv - 1
          Nsv = MATERIAL(PLATT(N)%PAR%MatID)%Nsv
          Nsv = Nsv * Ipts_PLATT(NPL)
          WRITE (IO_UNIT%LRDO) PLATT(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMP4
          NPL = N
          Isv = PLATQ(N)%PAR%Isv - 1
          Nsv = MATERIAL(PLATQ(N)%PAR%MatID)%Nsv
          Nsv = Nsv * Ipts_PLATQ(NPL)
          WRITE (IO_UNIT%LRDO) PLATQ(N)%PAR%EleID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMBM
          Isv = BEAM(N)%PAR%Isv - 1
          Nsv = 16 * MATERIAL(BEAM(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) BEAM(N)%PAR%EleID,Nsv,                          &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        DO N = 1,NUMSP
          Isv = SPRING(N)%PAR%Isv - 1
          Nsv = MATERIAL(SPRING(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) SPRING(N)%PAR%EleID,Nsv,                        &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMDM
          Isv = DAMPER(N)%PAR%Isv - 1
          Nsv = MATERIAL(DAMPER(N)%PAR%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) DAMPER(N)%PAR%EleID,Nsv,                        &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        DO N = 1,NUMSC
          Isv = SPRING_BC(N)%Isv - 1
          Nsv = MATERIAL(SPRING_BC(N)%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) SPRING_BC(N)%SprID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
        DO N = 1,NUMVC
          Isv = DAMPER_BC(N)%Isv - 1
          Nsv = MATERIAL(DAMPER_BC(N)%MatID)%Nsv
          WRITE (IO_UNIT%LRDO) DAMPER_BC(N)%DprID,Nsv,                         &
     &    (STATE_VARIABLES(Isv+i),i=1,Nsv),(0.0,i=Nsv+1,MAXSV)
        ENDDO
!!
        WRITE (IO_UNIT%LRDO) (NRBC_DATA(i),       i=1,NUMND)
        WRITE (IO_UNIT%LRDO) (SLIDING_NODE(i),    i=1,NUMSN)
        WRITE (IO_UNIT%LRDO) (CONTACT_SURFACE(i), i=1,NUMCE)
        WRITE (IO_UNIT%LRDO) (CONTACT_NODE(i),    i=1,NUMCN)
!!
        CLOSE (UNIT=IO_UNIT%LRDO, STATUS='KEEP')
!!
!! Inform user that a restart file has been written.
!!
        WRITE (MSGF,'(1PE12.4)') TIMSIM%Total
        WRITE (MSG1,'(I8)')      TIMSIM%Step
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'WRITE_RESTART_DATA.001.00'//                            &
     &          MSGL//'Restart File Written.'//                                &
     &          MSGL//'Simulation Time:'//MSGF//                               &
     &          MSGL//'Simulation Step:'//MSG1//                               &
     &          MSGL//'Output Sequence:'//CRDO                                 &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE MASS_PROPERTY_CALCULATION
!!
!! Copyright (c) by KEY Associates, 8-DEC-1991 11:59:54
!!
!! Purpose: For either a deformable or rigid body domain calculate the current
!! mass property characteristics: mass, center of mass, inertia, momemtum, KE,
!! angular momentum and angular KE.
!!
      USE shared_common_data
      USE massprop_
      USE node_
      USE motion_
      USE node_set_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID
      LOGICAL                                                                  &
     &          NEXT_NP_ID
!!
!! Provide for one-time call from TOTAL_MASS_REPORT to this module.
!!
      IF (NUMMP .EQ. 0) THEN
        NMPbgn = 0
        NMPend = 0
      ELSE
        NMPbgn = 1
        NMPend = NUMMP
      ENDIF
!!
!! Find total mass and center of mass for each "mass property" domain. Compute
!! velocity of the center of mass based on conservation of linear momentum.
!!
      DO NMP = NMPbgn,NMPend
        SetID = MASSPROP(NMP)%SetID
        Ams = 0.0
        Qcx = 0.0
        Qcy = 0.0
        Qcz = 0.0
        Xmv = 0.0
        Ymv = 0.0
        Zmv = 0.0
        Qke = 0.0
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Compute mass and first moment of mass.
!!
          Ams = Ams + NODE(N)%Mass
          Qcx = Qcx + NODE(N)%Mass * (MOTION(N)%Px + MOTION(N)%Ux)
          Qcy = Qcy + NODE(N)%Mass * (MOTION(N)%Py + MOTION(N)%Uy)
          Qcz = Qcz + NODE(N)%Mass * (MOTION(N)%Pz + MOTION(N)%Uz)
!!
!! Compute linear momentum.
!!
          Xmv = Xmv + NODE(N)%Mass * MOTION(N)%Vx
          Ymv = Ymv + NODE(N)%Mass * MOTION(N)%Vy
          Zmv = Zmv + NODE(N)%Mass * MOTION(N)%Vz
!!
!! Compute kinetic energy.
!!
          Qke = Qke + NODE(N)%Mass *                                           &
     &          (                                                              &
     &          MOTION(N)%Vx * MOTION(N)%Vx +                                  &
     &          MOTION(N)%Vy * MOTION(N)%Vy +                                  &
     &          MOTION(N)%Vz * MOTION(N)%Vz                                    &
     &          )
        ENDDO
!!
!! Store total mass of domain.
!!
        MASSPROP(NMP)%Mass = Ams
!!
!! Compute and store center of mass coordinates.
!!
        Xcm = Qcx * (ONE / Ams)
        Ycm = Qcy * (ONE / Ams)
        Zcm = Qcz * (ONE / Ams)
        MASSPROP(NMP)%Xcm  = Xcm
        MASSPROP(NMP)%Ycm  = Ycm
        MASSPROP(NMP)%Zcm  = Zcm
!!
!! Store momemtum.
!!
        MASSPROP(NMP)%Xmv  = Xmv
        MASSPROP(NMP)%Ymv  = Ymv
        MASSPROP(NMP)%Zmv  = Zmv
!!
!! Compute the velocity of the center of mass based on the conservation of
!! linear momentum.
!!
        Vxcm = Xmv / Ams
        Vycm = Ymv / Ams
        Vzcm = Zmv / Ams
        MASSPROP(NMP)%Vxcm = Vxcm
        MASSPROP(NMP)%Vycm = Vycm
        MASSPROP(NMP)%Vzcm = Vzcm
!!
!! Store kinetic energy of the "mass property" domain.
!!
        MASSPROP(NMP)%KE = 0.5D+0 * Qke
!!
!! Compute rotational characteristics wrt to the point specified based on
!! conservation of angular momemtum.
!!
        IF (MASSPROP(NMP)%Irot .GT. 0) THEN
!!
!! Base rotational properties on current center of mass.
!!
          IF (MASSPROP(NMP)%Irot .EQ. 1) THEN
            Xnert  = Xcm
            Ynert  = Ycm
            Znert  = Zcm
            Vxnert = Vxcm
            Vynert = Vycm
            Vznert = Vzcm
!!
!! Base rotational properties on origin of spatial coordinates.
!!
          ELSE IF (MASSPROP(NMP)%Irot .EQ. 2) THEN
            Xnert  = 0.0
            Ynert  = 0.0
            Znert  = 0.0
            Vxnert = 0.0
            Vynert = 0.0
            Vznert = 0.0
!!
!! Base rotational properties on point specified by the user.
!!
          ELSE IF (MASSPROP(NMP)%Irot .EQ. 3) THEN
            Xnert  = MASSPROP(NMP)%Xnert
            Ynert  = MASSPROP(NMP)%Ynert
            Znert  = MASSPROP(NMP)%Znert
            Vxnert = MASSPROP(NMP)%Vxnert
            Vynert = MASSPROP(NMP)%Vynert
            Vznert = MASSPROP(NMP)%Vznert
          ENDIF
!!
!! Save data on point with respect to which rotational properties are
!! computed.
!!
          MASSPROP(NMP)%Xnert  = Xnert
          MASSPROP(NMP)%Ynert  = Ynert
          MASSPROP(NMP)%Znert  = Znert
          MASSPROP(NMP)%Vxnert = Vxnert
          MASSPROP(NMP)%Vynert = Vynert
          MASSPROP(NMP)%Vznert = Vznert
!!
!! Compute rotational properties.
!!
          Oxmv = 0.0
          Oymv = 0.0
          Ozmv = 0.0
          Bxx = 0.0
          Byy = 0.0
          Bzz = 0.0
          Bxy = 0.0
          Bxz = 0.0
          Byz = 0.0
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Compute angular momentum
!!
            Pxcm = MOTION(N)%Px + MOTION(N)%Ux - Xnert
            Pycm = MOTION(N)%Py + MOTION(N)%Uy - Ynert
            Pzcm = MOTION(N)%Pz + MOTION(N)%Uz - Znert
            dVxc = MOTION(N)%Vx - Vxnert
            dVyc = MOTION(N)%Vy - Vynert
            dVzc = MOTION(N)%Vz - Vznert
            Oxmv  = Oxmv + NODE(N)%Mass * (Pycm*dVzc - Pzcm*dVyc)
            Oymv  = Oymv + NODE(N)%Mass * (Pzcm*dVxc - Pxcm*dVzc)
            Ozmv  = Ozmv + NODE(N)%Mass * (Pxcm*dVyc - Pycm*dVxc)
!!
!! Compute inertia tensor B from nodal point masses.
!!
            RSQ = Pxcm*Pxcm + Pycm*Pycm + Pzcm*Pzcm
            Bxx = Bxx + (RSQ - Pxcm*Pxcm) * NODE(N)%Mass
            Byy = Byy + (RSQ - Pycm*Pycm) * NODE(N)%Mass
            Bzz = Bzz + (RSQ - Pzcm*Pzcm) * NODE(N)%Mass
            Bxy = Bxy - Pxcm*Pycm * NODE(N)%Mass
            Bxz = Bxz - Pxcm*Pzcm * NODE(N)%Mass
            Byz = Byz - Pycm*Pzcm * NODE(N)%Mass
          ENDDO
          MASSPROP(NMP)%B(1) = Bxx
          MASSPROP(NMP)%B(2) = Byy
          MASSPROP(NMP)%B(3) = Bzz
          MASSPROP(NMP)%B(4) = Bxy
          MASSPROP(NMP)%B(5) = Bxz
          MASSPROP(NMP)%B(6) = Byz
!!
!! Invert the inertia tensor B.
!!
          Det = (Bxx*Byy)*Bzz                                                  &
     &        + (Bxy*Byz)*(Bxz+Bxz)                                            &
     &        - (Bxz*Bxz)*Byy                                                  &
     &        - (Byz*Byz)*Bxx                                                  &
     &        - (Bxy*Bxy)*Bzz
          IF (ABS (Det) .LT. 1.0D-30) THEN
            WRITE (MSG1,'(I8)') MASSPROP(NMP)%MPID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'MASS_PROPERTY_CALCULATION.001.00'//                     &
     &          MSGL//'Inertia Tensor For MASSPROP ID:'//MSG1//                &
     &          MSGL//'Is Singular. Do all nodes lie on a line?'               &
     &          )
            Cxx = 0.0
            Cyy = 0.0
            Czz = 0.0
            Cxy = 0.0
            Cxz = 0.0
            Cyz = 0.0
          ELSE
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
          ENDIF
!!
!! Store angular momentum and angular velocity.
!!
          MASSPROP(NMP)%Oxmv = Oxmv
          MASSPROP(NMP)%Oymv = Oymv
          MASSPROP(NMP)%Ozmv = Ozmv
          Ox = Cxx*Oxmv + Cxy*Oymv + Cxz*Ozmv
          Oy = Cxy*Oxmv + Cyy*Oymv + Cyz*Ozmv
          Oz = Cxz*Oxmv + Cyz*Oymv + Czz*Ozmv
!!
!! Convert angular momentum into magnitude and direction.
!!
          Omega = SQRT (Ox*Ox + Oy*Oy + Oz*Oz)
          IF (Omega .GT. 1.0D-27) THEN
            MASSPROP(NMP)%Omega = Omega
            MASSPROP(NMP)%Ax = Ox * (ONE / Omega)
            MASSPROP(NMP)%Ay = Oy * (ONE / Omega)
            MASSPROP(NMP)%Az = Oz * (ONE / Omega)
          ELSE
            MASSPROP(NMP)%Omega = 0.0
            MASSPROP(NMP)%Ax = 0.0
            MASSPROP(NMP)%Ay = 0.0
            MASSPROP(NMP)%Az = 0.0
          ENDIF
        ENDIF
!!
!! End of mass property do-loop
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE NEXT_TIME_STEP ( DTnext )
!!
!! Copyright (c) by KEY Associates, 10-FEB-1991 21:33:25
!!
      USE shared_common_data
      USE hexah_
      USE penta_
      USE tetra_
      USE lsold_
      USE membt_
      USE membq_
      USE truss_
      USE platt_
      USE platq_
      USE beam_
      USE spring_
      USE damper_
      USE spring_bc_
      USE damper_bc_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument
      REAL(KIND(0D0)), INTENT(OUT) :: DTnext  ! Next time step to use
!!
!! Local variables.
      INTEGER :: HstID
      REAL(KIND(0D0)) :: TABLE_LOOK_UP
!!
      TIMSIM%DTmax  = MAX                                                      &
     &          (                                                              &
     &          TIMSIM%DTHxx,                                                  &
     &          TIMSIM%DTPnx,                                                  &
     &          TIMSIM%DTTtx,                                                  &
     &          TIMSIM%DTLSx,                                                  &
     &          TIMSIM%DTM3x,                                                  &
     &          TIMSIM%DTM4x,                                                  &
     &          TIMSIM%DTTrx,                                                  &
     &          TIMSIM%DTP3x,                                                  &
     &          TIMSIM%DTP4x,                                                  &
     &          TIMSIM%DTBmx,                                                  &
     &          TIMSIM%DTSpx,                                                  &
     &          TIMSIM%DTDmx,                                                  &
     &          TIMSIM%DTSBx,                                                  &
     &          TIMSIM%DTDBx                                                   &
     &          )

      TIMSIM%DTmin = HUGE( TIMSIM%DTmin )

      IF (NUMHX .GT. 0) THEN
        IF (TIMSIM%DTHex(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTHex(1)
          TIMSIM%Ncrit = HEXAH(TIMSIM%Hexah(1))%PAR%EleID
          TIMSIM%Source = '      3-D Hexahedron'
        ENDIF
      ENDIF

      IF (NUMPX .GT. 0) THEN
        IF (TIMSIM%DTPen(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTPen(1)
          TIMSIM%Ncrit = PENTA(TIMSIM%Penta(1))%PAR%EleID
          TIMSIM%Source = '      3-D Pentrahedron'
        ENDIF
      ENDIF

      IF (NUMTX .GT. 0) THEN
        IF (TIMSIM%DTTet(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTTet(1)
          TIMSIM%Ncrit = TETRA(TIMSIM%Tetra(1))%PAR%EleID
          TIMSIM%Source = '      3-D Tetrahedron'
        ENDIF
      ENDIF

      IF (NUMLS .GT. 0) THEN
        IF (TIMSIM%DTLYS(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTLYS(1)
          TIMSIM%Ncrit = LSOLD(TIMSIM%LSold(1))%PAR%EleID
          TIMSIM%Source = '     3-D Layered Solid'
        ENDIF
      ENDIF

      IF (NUMM3 .GT. 0) THEN
        IF (TIMSIM%DTMb3(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTMb3(1)
          TIMSIM%Ncrit = MEMBT(TIMSIM%Memb3(1))%PAR%EleID
          TIMSIM%Source = '    Triangular Membrane'
        ENDIF
      ENDIF

      IF (NUMM4 .GT. 0) THEN
        IF (TIMSIM%DTMb4(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTMb4(1)
          TIMSIM%Ncrit = MEMBQ(TIMSIM%Memb4(1))%PAR%EleID
          TIMSIM%Source = '   Quadrilateral Membrane'
        ENDIF
      ENDIF

      IF (NUMTR .GT. 0) THEN
        IF (TIMSIM%DTTru(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTTru(1)
          TIMSIM%Ncrit = TRUSS(TIMSIM%Truss(1))%PAR%EleID
          TIMSIM%Source = '         1-D Truss'
        ENDIF
      ENDIF

      IF (NUMP3 .GT. 0) THEN
        IF (TIMSIM%DTPl3(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTPl3(1)
          TIMSIM%Ncrit = PLATT(TIMSIM%Plat3(1))%PAR%EleID
          TIMSIM%Source = '      Triangular Plate'
        ENDIF
      ENDIF

      IF (NUMP4 .GT. 0) THEN
        IF (TIMSIM%DTPl4(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTPl4(1)
          TIMSIM%Ncrit = PLATQ(TIMSIM%Plat4(1))%PAR%EleID
          TIMSIM%Source = '     Quadrilateral Plate'
        ENDIF
      ENDIF

      IF (NUMBM .GT. 0) THEN
        IF (TIMSIM%DTBms(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTBms(1)
          TIMSIM%Ncrit = BEAM(TIMSIM%Beams(1))%PAR%EleID
          TIMSIM%Source = '           Beam'
        ENDIF
      ENDIF

      IF (NUMSP .GT. 0) THEN
        IF (TIMSIM%DTSpr(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTSpr(1)
          TIMSIM%Ncrit = SPRING(TIMSIM%Spring(1))%PAR%EleID
          TIMSIM%Source = '        1-D Spring'
        ENDIF
      ENDIF

!!!   IF (NUMDM .GT. 0) THEN  !  (not yet implemented)
!!!     IF (TIMSIM%DTDmp(1) .LT. TIMSIM%DTmin) THEN
!!!       TIMSIM%DTmin = TIMSIM%DTDmp(1)
!!!       TIMSIM%Ncrit = DAMPER(TIMSIM%Damper(1))%PAR%EleID
!!!       TIMSIM%Source = '        1-D Damper'
!!!     ENDIF
!!!   ENDIF

      IF (NUMSC .GT. 0) THEN
        IF (TIMSIM%DTSBC(1) .LT. TIMSIM%DTmin) THEN
          TIMSIM%DTmin = TIMSIM%DTSBC(1)
          TIMSIM%Ncrit = SPRING_BC(TIMSIM%SprBC(1))%SprID
          TIMSIM%Source = '         Spring BC'
        ENDIF
      ENDIF

!!!   IF (NUMVC .GT. 0) THEN   !  (not yet implemented)
!!!     IF (TIMSIM%DTDBC(1) .LT. TIMSIM%DTmin) THEN
!!!       TIMSIM%DTmin = TIMSIM%DTDBC(1)
!!!       TIMSIM%Ncrit = DAMPER_BC(TIMSIM%DmpBC(1))%DprID
!!!       TIMSIM%Source = '         Damper BC'
!!!     ENDIF
!!!   ENDIF
!!
!BEGIN SPEC_omp2001
!!
      TIMSIM%Ncrit = 0
      TIMSIM%DTmin = REAL(REAL(TIMSIM%DTmin,KIND(0E0)),KIND(0D0))
!!
!END SPEC_omp2001
!!
!! Check to see if there is a tabulated function capping the time step.
!!
      HstID = TIMSIM%DTcntrl
      IF (HstID .EQ. 0) THEN
        DTnext = PARAMVALUE%DTscale * TIMSIM%DTmin
      ELSE
        PaDTmin = PARAMVALUE%DTscale * TIMSIM%DTmin
        TaTotal = TABLE_LOOK_UP (HstID,TIMSIM%Total)
        DTnext = MIN ( PaDTmin, TaTotal )
!!
!! Check to see if current simulation time TIMSIM%Total exceeds last tabulated
!! function entry. If it does, terminate use of the tabulated function by
!! re-setting the pointer TIMSIM%DTcntrl to zero.
!!
        NTP = TABULATED_FUNCTION(HstID)%Number_of_Pairs
        IF (TABULATED_FUNCTION(HstID)%X(NTP) .LT. TIMSIM%Total) THEN
          TIMSIM%DTcntrl = 0
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0)) FUNCTION TABLE_LOOK_UP ( Table_ID,X_Value )
!!
!! Copyright (c) by KEY Associates, 10-FEB-1991 22:11:56
!!
!! Locate the position of X_Value in the argument table TABLE(ID).X(*) and
!! interpolate or extrapolate for the Y_Value from the table TABLE(ID).Y(*).
!!
      USE shared_common_data
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: Table_ID
      REAL(KIND(0D0)), INTENT(IN) :: X_Value
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Pre-compute the individual segment slopes. This will occur the first
!! time this function is accessed.
!!
      IF (FIRST) THEN
        ERROR%COUNT = 0
        DO NTF = 1,NUMTF
          WRITE (MSG1,'(I8)') TABULATED_FUNCTION(NTF)%TFID
          DO i = 1,TABULATED_FUNCTION(NTF)%Number_of_Pairs-1

            Delta_Y = TABULATED_FUNCTION(NTF)%Y(i+1) -                         &
     &                TABULATED_FUNCTION(NTF)%Y(  i)
            Delta_X = TABULATED_FUNCTION(NTF)%X(i+1) -                         &
     &                TABULATED_FUNCTION(NTF)%X(  i)

            IF (Delta_X .LE. 1.0D-35) THEN

              WRITE (MSG2,'(I8)') i
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'TABLE_LOOK_UP.001.01'//                                 &
     &          MSGL//'For Tabulated Function ID:'//MSG1//                     &
     &          MSGL//'Delta X For Interval Number:'//MSG2//                   &
     &          MSGL//'Is Negative, Zero Or Too Close To Zero.'                &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1

            ELSE

              IF (ABS(Delta_Y) .LT. 1.0D-35) THEN
                Exponent = 0.0
              ELSE
                Exponent = LOG10(ABS(Delta_Y)) - LOG10(ABS(Delta_X))
              ENDIF

              IF (Exponent .GT. 35.0) THEN

                WRITE (MSG2,'(I8)') i
                CALL USER_MESSAGE                                              &
     &                  (                                                      &
     &                  MSGL//'WARN'//                                         &
     &                  MSGL//'TABLE_LOOK_UP.001.02'//                         &
     &                  MSGL//'For Tabulated Function ID:'//MSG1//             &
     &                  MSGL//'Slope For Interval Number:'//MSG2//             &
     &                  MSGL//'Is Too Close To Infinity.'                      &
     &                  )
                ERROR%COUNT = ERROR%COUNT + 1

              ELSE

                TABULATED_FUNCTION(NTF)%SLOPE(i) = Delta_Y / Delta_X

              ENDIF
            ENDIF
          ENDDO

          TABULATED_FUNCTION(NTF)%Last_Location = 1

        ENDDO

        IF (ERROR%COUNT .GT. 0) THEN
          WRITE (MSG1,'(I8)') ERROR%COUNT
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'TABLE_LOOK_UP.001.03'//                                 &
     &          MSGL//'Total Number Of Near Infinite Slopes,'//                &
     &          MSGL//'Or Negative/Zero X-Axis Intervals:'//MSG1               &
     &          )
        ENDIF
        FIRST = .FALSE.
      ENDIF
!!
!! Interogate table. Start at the point where the last value was found
!! provided the x-value is greater than the last location accessed. Other-
!! wise reset the table search to start with the first interval.
!!
      TABLE_LOOK_UP = 0.0
      IF (Table_ID .GE. 1 .AND. Table_ID .LE. NUMTF) THEN

        NTP = TABULATED_FUNCTION(Table_ID)%Number_of_Pairs
        IF (NTP .EQ. 1) THEN

          TABLE_LOOK_UP = TABULATED_FUNCTION(Table_ID)%Y(1)

        ELSE

          i = TABULATED_FUNCTION(Table_ID)%Last_location
          IF (TABULATED_FUNCTION(Table_ID)%X(i) .GT. X_Value) i = 1
          DO WHILE                                                             &
     &    (TABULATED_FUNCTION(Table_ID)%X(i).LT.X_Value .AND. i.LT.NTP)
            i = i + 1
          ENDDO
          m = MAX (1,i-1)

          TABLE_LOOK_UP = TABULATED_FUNCTION(Table_ID)%Y(m) +                  &
     &                    TABULATED_FUNCTION(Table_ID)%SLOPE(m) *              &
     &                    (X_Value-TABULATED_FUNCTION(Table_ID)%X(m))
          TABULATED_FUNCTION(Table_ID)%Last_Location = m

        ENDIF

      ELSE IF (Table_ID .NE. 0) THEN
!!
!! This module is frequently called with Table_ID equal to zero, thus,
!! this section is only needed (can only be used) for ID < 0 or for
!! ID > NUMTF.
!!
        WRITE (MSG1,'(I8)') Table_ID
        WRITE (MSG2,'(I8)') NUMTF
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'TABLE_LOOK_UP.002.00'//                                 &
     &          MSGL//'Module Accessed With A Tabulated '//                    &
     &          MSGL//'Function Internal Index:'//MSG1//                       &
     &          MSGL//'Smaller Than 0 Or Greater Than NUMTF:'//MSG2            &
     &          )

      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE TIME_INITIALIZATION
!!
!! Copyright (c) by KEY Associates, 12-FEB-1991 19:58:58
!!
!! Set time data min's and max's to limit values, clear critical elements.
!!
      USE shared_common_data
!!
      TIMSIM%DTHxx  = 0.0
      TIMSIM%DTPnx  = 0.0
      TIMSIM%DTTtx  = 0.0
      TIMSIM%DTLSx  = 0.0
      TIMSIM%DTM3x  = 0.0
      TIMSIM%DTM4x  = 0.0
      TIMSIM%DTTrx  = 0.0
      TIMSIM%DTP3x  = 0.0
      TIMSIM%DTP4x  = 0.0
      TIMSIM%DTBmx  = 0.0
      TIMSIM%DTSpx  = 0.0
      TIMSIM%DTDmx  = 0.0
      TIMSIM%DTSBx  = 0.0
      TIMSIM%DTDBx  = 0.0
!!!   TIMSIM%Ncrit  = 0
      DO i = 1,NPNDT
        TIMSIM%DTHex(i)  = HUGE( TIMSIM%DTHex(i) )
        TIMSIM%DTPen(i)  = HUGE( TIMSIM%DTPen(i) )
        TIMSIM%DTTet(i)  = HUGE( TIMSIM%DTTet(i) )
        TIMSIM%DTLYS(i)  = HUGE( TIMSIM%DTLYS(i) )
        TIMSIM%DTMb3(i)  = HUGE( TIMSIM%DTMb3(i) )
        TIMSIM%DTMb4(i)  = HUGE( TIMSIM%DTMb4(i) )
        TIMSIM%DTTru(i)  = HUGE( TIMSIM%DTTru(i) )
        TIMSIM%DTPl3(i)  = HUGE( TIMSIM%DTPl3(i) )
        TIMSIM%DTPl4(i)  = HUGE( TIMSIM%DTPl4(i) )
        TIMSIM%DTBms(i)  = HUGE( TIMSIM%DTBms(i) )
        TIMSIM%DTSpr(i)  = HUGE( TIMSIM%DTSpr(i) )
        TIMSIM%DTDmp(i)  = HUGE( TIMSIM%DTDmp(i) )
        TIMSIM%DTSBC(i)  = HUGE( TIMSIM%DTSBC(i) )
        TIMSIM%DTDBC(i)  = HUGE( TIMSIM%DTDBC(i) )
        TIMSIM%Hexah(i)  = 0
        TIMSIM%Penta(i)  = 0
        TIMSIM%Tetra(i)  = 0
        TIMSIM%LSold(i)  = 0
        TIMSIM%Memb3(i)  = 0
        TIMSIM%Memb4(i)  = 0
        TIMSIM%Truss(i)  = 0
        TIMSIM%Plat3(i)  = 0
        TIMSIM%Plat4(i)  = 0
        TIMSIM%Beams(i)  = 0
        TIMSIM%Spring(i) = 0
        TIMSIM%Damper(i) = 0
        TIMSIM%SprBC(i)  = 0
        TIMSIM%DmpBC(i)  = 0
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE xSCATTER_ELEMENT_NODAL_FORCES
!!
!! Copyright (c) by KEY Associates, 3-NOV-1991 21:15:04
!!
      USE shared_common_data
      USE hexah_
      USE penta_
      USE tetra_
      USE lsold_
      USE membt_
      USE membq_
      USE truss_
      USE platt_
      USE platq_
      USE beam_
      USE spring_
      USE damper_
      USE force_
!!
!! Scatter element internal forces to nodes.
!!
      DO N = 1,NUMHX
        DO i = 1,8
          IX = HEXAH(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + HEXAH(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + HEXAH(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + HEXAH(N)%RES%Zint(i)
        ENDDO
      ENDDO
      DO N = 1,NUMPX
        DO i = 1,6
          IX = PENTA(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + PENTA(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + PENTA(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + PENTA(N)%RES%Zint(i)
        ENDDO
      ENDDO
      DO N = 1,NUMTX
        DO i = 1,4
          IX = TETRA(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + TETRA(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + TETRA(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + TETRA(N)%RES%Zint(i)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMLS
        DO i = 1,8
          IX = LSOLD(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + LSOLD(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + LSOLD(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + LSOLD(N)%RES%Zint(i)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMM3
        IX = MEMBT(N)%PAR%IX(1)
        FORCE(IX)%Xint = FORCE(IX)%Xint + MEMBT(N)%RES%Xint(1)
        FORCE(IX)%Yint = FORCE(IX)%Yint + MEMBT(N)%RES%Yint(1)
        FORCE(IX)%Zint = FORCE(IX)%Zint + MEMBT(N)%RES%Zint(1)
        IX = MEMBT(N)%PAR%IX(2)
        FORCE(IX)%Xint = FORCE(IX)%Xint + MEMBT(N)%RES%Xint(2)
        FORCE(IX)%Yint = FORCE(IX)%Yint + MEMBT(N)%RES%Yint(2)
        FORCE(IX)%Zint = FORCE(IX)%Zint + MEMBT(N)%RES%Zint(2)
        IX = MEMBT(N)%PAR%IX(3)
        FORCE(IX)%Xint = FORCE(IX)%Xint + MEMBT(N)%RES%Xint(3)
        FORCE(IX)%Yint = FORCE(IX)%Yint + MEMBT(N)%RES%Yint(3)
        FORCE(IX)%Zint = FORCE(IX)%Zint + MEMBT(N)%RES%Zint(3)
      ENDDO
      DO N = 1,NUMM4
        DO i = 1,4
          IX = MEMBQ(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + MEMBQ(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + MEMBQ(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + MEMBQ(N)%RES%Zint(i)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMTR
        IX = TRUSS(N)%PAR%IX(1)
        FORCE(IX)%Xint = FORCE(IX)%Xint + TRUSS(N)%RES%Xint(1)
        FORCE(IX)%Yint = FORCE(IX)%Yint + TRUSS(N)%RES%Yint(1)
        FORCE(IX)%Zint = FORCE(IX)%Zint + TRUSS(N)%RES%Zint(1)
        IX = TRUSS(N)%PAR%IX(2)
        FORCE(IX)%Xint = FORCE(IX)%Xint + TRUSS(N)%RES%Xint(2)
        FORCE(IX)%Yint = FORCE(IX)%Yint + TRUSS(N)%RES%Yint(2)
        FORCE(IX)%Zint = FORCE(IX)%Zint + TRUSS(N)%RES%Zint(2)
      ENDDO
!!
      DO N = 1,NUMP3
        DO i = 1,6
          IX = PLATT(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + PLATT(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + PLATT(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + PLATT(N)%RES%Zint(i)
        ENDDO
      ENDDO
      DO N = 1,NUMP4
        DO i = 1,8
          IX = PLATQ(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + PLATQ(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + PLATQ(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + PLATQ(N)%RES%Zint(i)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMBM
        DO i = 1,4
          IX = BEAM(N)%PAR%IX(i)
          FORCE(IX)%Xint = FORCE(IX)%Xint + BEAM(N)%RES%Xint(i)
          FORCE(IX)%Yint = FORCE(IX)%Yint + BEAM(N)%RES%Yint(i)
          FORCE(IX)%Zint = FORCE(IX)%Zint + BEAM(N)%RES%Zint(i)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMSP
        IX = SPRING(N)%PAR%IX(1)
        FORCE(IX)%Xint = FORCE(IX)%Xint + SPRING(N)%RES%Xint(1)
        FORCE(IX)%Yint = FORCE(IX)%Yint + SPRING(N)%RES%Yint(1)
        FORCE(IX)%Zint = FORCE(IX)%Zint + SPRING(N)%RES%Zint(1)
        IX = SPRING(N)%PAR%IX(2)
        FORCE(IX)%Xint = FORCE(IX)%Xint + SPRING(N)%RES%Xint(2)
        FORCE(IX)%Yint = FORCE(IX)%Yint + SPRING(N)%RES%Yint(2)
        FORCE(IX)%Zint = FORCE(IX)%Zint + SPRING(N)%RES%Zint(2)
      ENDDO
      DO N = 1,NUMDM
        IX = DAMPER(N)%PAR%IX(1)
        FORCE(IX)%Xint = FORCE(IX)%Xint + DAMPER(N)%RES%Xint(1)
        FORCE(IX)%Yint = FORCE(IX)%Yint + DAMPER(N)%RES%Yint(1)
        FORCE(IX)%Zint = FORCE(IX)%Zint + DAMPER(N)%RES%Zint(1)
        IX = DAMPER(N)%PAR%IX(2)
        FORCE(IX)%Xint = FORCE(IX)%Xint + DAMPER(N)%RES%Xint(2)
        FORCE(IX)%Yint = FORCE(IX)%Yint + DAMPER(N)%RES%Yint(2)
        FORCE(IX)%Zint = FORCE(IX)%Zint + DAMPER(N)%RES%Zint(2)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE SCATTER_ELEMENT_NODAL_FORCES
!!
!! Copyright (c) by KEY Associates, 29-JUL-2000 13:37:29
!!
      USE shared_common_data
      USE hexah_
      USE penta_
      USE tetra_
      USE lsold_
      USE membt_
      USE membq_
      USE truss_
      USE platt_
      USE platq_
      USE beam_
      USE spring_
      USE damper_
      USE force_
!!
!! Local variables.
      INTEGER       :: IX(8)
      INTEGER, SAVE :: NGRPS
      LOGICAL       :: UNCOUPLED
      LOGICAL, SAVE :: FIRST = .TRUE.

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPGRP
!!
!! Define parallel-safe accumulation groups.
!!
      IF (FIRST) THEN

        NBGN_HEXAH  = 1
        NBGN_PENTA  = 1
        NBGN_TETRA  = 1
        NBGN_LSOLD  = 1
        NBGN_MEMBT  = 1
        NBGN_MEMBQ  = 1
        NBGN_TRUSS  = 1
        NBGN_PLATT  = 1
        NBGN_PLATQ  = 1
        NBGN_BEAM   = 1
        NBGN_SPRING = 1
        NBGN_DAMPER = 1
!!
!! Clear the group number from the element data. (Needed in case this
!! is a restart run.) Note: we are counting on a clever Fortran-90
!! compiler to clear all entries, that is, HEXAH(1:NUMHX)%PAR%IGR = 0,
!! et cetera
!!

        DO N = 1,NUMHX
          HEXAH(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMPX
          PENTA(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMTX
          TETRA(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMLS
          LSOLD(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMM3
          MEMBT(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMM4
          MEMBQ(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMTR
          TRUSS(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMP3
          PLATT(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMP4
          PLATQ(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMBM
          BEAM(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMSP
          SPRING(N)%PAR%IGR = 0
        ENDDO
        DO N = 1,NUMDM
          DAMPER(N)%PAR%IGR = 0
        ENDDO
!!
!! Allocate scratch storage.
!!
        ALLOCATE (NPGRP(1:NUMNP), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NPGRP',NUMNP)
        NPGRP = 0
!!
        NGRPS = 0
 100    CONTINUE
        NLAST = NGRPS
        NEXTG = NGRPS + 1

        FIRST = .TRUE.
        NEXT_NBGN = NUMHX + 1
        DO n = NBGN_HEXAH,NUMHX
          IF (HEXAH(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:8) = HEXAH(n)%PAR%IX(1:8)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,8)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:8)) = (/(NGRPS,i=1,8)/)
              HEXAH(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_HEXAH = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMPX + 1
        DO n = NBGN_PENTA,NUMPX
          IF (PENTA(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:6) = PENTA(n)%PAR%IX(1:6)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,6)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:6)) = (/(NGRPS,i=1,6)/)
              PENTA(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_PENTA = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMTX + 1
        DO n = NBGN_TETRA,NUMTX
          IF (TETRA(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:4) = TETRA(n)%PAR%IX(1:4)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,4)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:4)) = (/(NGRPS,i=1,4)/)
              TETRA(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_TETRA = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMLS + 1
        DO n = NBGN_LSOLD,NUMLS
          IF (LSOLD(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:8) = LSOLD(n)%PAR%IX(1:8)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,8)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:8)) = (/(NGRPS,i=1,8)/)
              LSOLD(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_LSOLD = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMM3 + 1
        DO n = NBGN_MEMBT,NUMM3
          IF (MEMBT(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:3) = MEMBT(n)%PAR%IX(1:3)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,3)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:3)) = (/(NGRPS,i=1,3)/)
              MEMBT(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_MEMBT = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMM4 + 1
        DO n = NBGN_MEMBQ,NUMM4
          IF (MEMBQ(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:4) = MEMBQ(n)%PAR%IX(1:4)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,4)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:4)) = (/(NGRPS,i=1,4)/)
              MEMBQ(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_MEMBQ = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMTR + 1
        DO n = NBGN_TRUSS,NUMTR
          IF (TRUSS(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:2) = TRUSS(n)%PAR%IX(1:2)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,2)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:2)) = (/(NGRPS,i=1,2)/)
              TRUSS(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_TRUSS = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMP3 + 1
        DO n = NBGN_PLATT,NUMP3
          IF (PLATT(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:3) = PLATT(n)%PAR%IX(1:3)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,3)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:3)) = (/(NGRPS,i=1,3)/)
              PLATT(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_PLATT = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMP4 + 1
        DO n = NBGN_PLATQ,NUMP4
          IF (PLATQ(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:4) = PLATQ(n)%PAR%IX(1:4)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,4)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:4)) = (/(NGRPS,i=1,4)/)
              PLATQ(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_PLATQ = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMBM + 1
        DO n = NBGN_BEAM,NUMBM
          IF (BEAM(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:2) = BEAM(n)%PAR%IX(1:2)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,2)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:2)) = (/(NGRPS,i=1,2)/)
              BEAM(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_BEAM = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMSP + 1
        DO n = NBGN_SPRING,NUMSP
          IF (SPRING(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:2) = SPRING(n)%PAR%IX(1:2)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,2)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:2)) = (/(NGRPS,i=1,2)/)
              SPRING(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_SPRING = NEXT_NBGN

        FIRST = .TRUE.
        NEXT_NBGN = NUMDM + 1
        DO n = NBGN_DAMPER,NUMDM
          IF (DAMPER(n)%PAR%IGR .EQ. 0) THEN
            IF (NGRPS .LT. NEXTG) NGRPS = NEXTG
            IX(1:2) = DAMPER(n)%PAR%IX(1:2)
            UNCOUPLED = ALL ((/(NPGRP(IX(i)) .LT. NGRPS,i=1,2)/))
            IF (UNCOUPLED) THEN
              NPGRP(IX(1:2)) = (/(NGRPS,i=1,2)/)
              DAMPER(n)%PAR%IGR = NGRPS
            ELSE
              IF (FIRST) THEN
                NEXT_NBGN = n
                FIRST = .FALSE.
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        NBGN_DAMPER = NEXT_NBGN

        IF (NLAST .LT. NGRPS) GO TO 100
!!
!! Deallocate scratch storage.
!!
        DEALLOCATE (NPGRP)

        WRITE (MSG1,'(I8)') NGRPS
        CALL USER_MESSAGE                                        &
     &    (                                                      &
     &    MSGL//'INFORM'//                                       &
     &    MSGL//'SCATTER_ELEMENT_NODAL_FORCES.001.00'//          &
     &    MSGL//'Number Of Parallel-Safe Scatter Groups:'//MSG1  &
     &    )

        FIRST = .FALSE.
      ENDIF
!!
!! Scatter element internal forces to nodes.
!!
      IF (NUMHX .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMHX
            IF (HEXAH(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:8) = HEXAH(N)%PAR%IX(1:8)
              FORCE(IX(1:8))%Xint = FORCE(IX(1:8))%Xint + HEXAH(N)%RES%Xint
              FORCE(IX(1:8))%Yint = FORCE(IX(1:8))%Yint + HEXAH(N)%RES%Yint
              FORCE(IX(1:8))%Zint = FORCE(IX(1:8))%Zint + HEXAH(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF (NUMPX .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMPX
            IF (PENTA(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:6) = PENTA(N)%PAR%IX(1:6)
              FORCE(IX(1:6))%Xint = FORCE(IX(1:6))%Xint + PENTA(N)%RES%Xint 
              FORCE(IX(1:6))%Yint = FORCE(IX(1:6))%Yint + PENTA(N)%RES%Yint
              FORCE(IX(1:6))%Zint = FORCE(IX(1:6))%Zint + PENTA(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF (NUMTX .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMTX
            IF (TETRA(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:4) = TETRA(N)%PAR%IX(1:4)
              FORCE(IX(1:4))%Xint = FORCE(IX(1:4))%Xint + TETRA(N)%RES%Xint 
              FORCE(IX(1:4))%Yint = FORCE(IX(1:4))%Yint + TETRA(N)%RES%Yint
              FORCE(IX(1:4))%Zint = FORCE(IX(1:4))%Zint + TETRA(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMLS .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMLS
            IF (LSOLD(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:8) = LSOLD(N)%PAR%IX(1:8)
              FORCE(IX(1:8))%Xint = FORCE(IX(1:8))%Xint + LSOLD(N)%RES%Xint 
              FORCE(IX(1:8))%Yint = FORCE(IX(1:8))%Yint + LSOLD(N)%RES%Yint
              FORCE(IX(1:8))%Zint = FORCE(IX(1:8))%Zint + LSOLD(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMM3 .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMM3
            IF (MEMBT(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:3) = MEMBT(N)%PAR%IX(1:3)
              FORCE(IX(1:3))%Xint = FORCE(IX(1:3))%Xint + MEMBT(N)%RES%Xint
              FORCE(IX(1:3))%Yint = FORCE(IX(1:3))%Yint + MEMBT(N)%RES%Yint
              FORCE(IX(1:3))%Zint = FORCE(IX(1:3))%Zint + MEMBT(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF (NUMM4 .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMM4
            IF (MEMBQ(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:4) = MEMBQ(N)%PAR%IX(1:4)
              FORCE(IX(1:4))%Xint = FORCE(IX(1:4))%Xint + MEMBQ(N)%RES%Xint 
              FORCE(IX(1:4))%Yint = FORCE(IX(1:4))%Yint + MEMBQ(N)%RES%Yint
              FORCE(IX(1:4))%Zint = FORCE(IX(1:4))%Zint + MEMBQ(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMTR .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMTR
            IF (TRUSS(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:2) = TRUSS(N)%PAR%IX(1:2)
              FORCE(IX(1:2))%Xint = FORCE(IX(1:2))%Xint + TRUSS(N)%RES%Xint
              FORCE(IX(1:2))%Yint = FORCE(IX(1:2))%Yint + TRUSS(N)%RES%Yint
              FORCE(IX(1:2))%Zint = FORCE(IX(1:2))%Zint + TRUSS(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMP3 .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMP3
            IF (PLATT(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:6) = PLATT(N)%PAR%IX(1:6)
              FORCE(IX(1:6))%Xint = FORCE(IX(1:6))%Xint + PLATT(N)%RES%Xint 
              FORCE(IX(1:6))%Yint = FORCE(IX(1:6))%Yint + PLATT(N)%RES%Yint
              FORCE(IX(1:6))%Zint = FORCE(IX(1:6))%Zint + PLATT(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMP4 .GT. 0) THEN
        DO NGR = 1,NGRPS
          CURRENT_GROUP = NGR
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,IX)
          DO N = 1,NUMP4
            IF (PLATQ(N)%PAR%IGR .EQ. CURRENT_GROUP) THEN
              IX(1:8) = PLATQ(N)%PAR%IX(1:8)
              FORCE(IX(1:8))%Xint = FORCE(IX(1:8))%Xint + PLATQ(N)%RES%Xint 
              FORCE(IX(1:8))%Yint = FORCE(IX(1:8))%Yint + PLATQ(N)%RES%Yint
              FORCE(IX(1:8))%Zint = FORCE(IX(1:8))%Zint + PLATQ(N)%RES%Zint
            ENDIF
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
      ENDIF
!!
      IF (NUMBM .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMBM
            IF (BEAM(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:4) = BEAM(N)%PAR%IX(1:4)
              FORCE(IX(1:4))%Xint = FORCE(IX(1:4))%Xint + BEAM(N)%RES%Xint 
              FORCE(IX(1:4))%Yint = FORCE(IX(1:4))%Yint + BEAM(N)%RES%Yint
              FORCE(IX(1:4))%Zint = FORCE(IX(1:4))%Zint + BEAM(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      IF (NUMSP .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMSP
            IF (SPRING(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:2) = SPRING(N)%PAR%IX(1:2)
              FORCE(IX(1:2))%Xint = FORCE(IX(1:2))%Xint + SPRING(N)%RES%Xint 
              FORCE(IX(1:2))%Yint = FORCE(IX(1:2))%Yint + SPRING(N)%RES%Yint
              FORCE(IX(1:2))%Zint = FORCE(IX(1:2))%Zint + SPRING(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF (NUMDM .GT. 0) THEN
        DO NGR = 1,NGRPS
          DO N = 1,NUMDM
            IF (DAMPER(N)%PAR%IGR .EQ. NGR) THEN
              IX(1:2) = DAMPER(N)%PAR%IX(1:2)
              FORCE(IX(1:2))%Xint = FORCE(IX(1:2))%Xint + DAMPER(N)%RES%Xint 
              FORCE(IX(1:2))%Yint = FORCE(IX(1:2))%Yint + DAMPER(N)%RES%Yint
              FORCE(IX(1:2))%Zint = FORCE(IX(1:2))%Zint + DAMPER(N)%RES%Zint
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_RIGID_BODY_MASS
!!
!! Copyright (c) by KEY Associates;  7-AUG-1993 14:12:09.74
!!
!! Purpose: For those rigid body domains with mass and inertia based on
!! the nodal point masses, initialize the total mass and inertia.
!!
      USE shared_common_data
      USE rigid_body_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! The material in this module draws heavily on the article:
!!
!!      J. L. Synge, "Classical Dynamics," Encylopedia of Physics,
!!      Volume III, Part 1, Springer-Verlag, Berlin, 1960.
!!
      ERROR%COUNT = 0
      DO i = 1,NUMRB

        IF (RIGID_BODY(i)%Prop .EQ. 0) THEN
!!
!! Calculate mass and inertia from finite element discretization. Find total
!! mass and center of mass for each rigid body domain. Compute velocity of
!! the center of mass based on conservation of linear momentum.
!!
          Ams = 0.0
          Qcx = 0.0
          Qcy = 0.0
          Qcz = 0.0
          Vmx = 0.0
          Vmy = 0.0
          Vmz = 0.0
          M = RIGID_BODY(i)%FirstNP
          DO WHILE (M .GT. 0)
            Ams = Ams + NODE(M)%Mass
            Qcx = Qcx + NODE(M)%Mass * MOTION(M)%Px
            Qcy = Qcy + NODE(M)%Mass * MOTION(M)%Py
            Qcz = Qcz + NODE(M)%Mass * MOTION(M)%Pz
            Vmx = Vmx + NODE(M)%Mass * MOTION(M)%Vx
            Vmy = Vmy + NODE(M)%Mass * MOTION(M)%Vy
            Vmz = Vmz + NODE(M)%Mass * MOTION(M)%Vz
            M = NODE(M)%IRB
          ENDDO
          IF (ABS (Ams) .LT. 1.0D-20) THEN
            WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'INITIALIZE_RIGID_BODY_MASS.001.01'//                    &
     &          MSGL//'Total Mass For Rigid Body ID:'//MSG1//                  &
     &          MSGL//'Is Zero. Cannot Proceed Beyond Initialization.'         &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            RIGID_BODY(i)%Mass = Ams
            Xcm  = Qcx * (ONE / Ams)
            Ycm  = Qcy * (ONE / Ams)
            Zcm  = Qcz * (ONE / Ams)
            Vxcm = Vmx * (ONE / Ams)
            Vycm = Vmy * (ONE / Ams)
            Vzcm = Vmz * (ONE / Ams)
          ENDIF
!!
!! Compute angular velocity of the center of mass based on conservation
!! of angular momemtum. (To get the final result the inertia tensor B
!! must also be computed.)
!!
          Omx = 0.0
          Omy = 0.0
          Omz = 0.0
          Bxx = 0.0
          Byy = 0.0
          Bzz = 0.0
          Bxy = 0.0
          Bxz = 0.0
          Byz = 0.0
          M = RIGID_BODY(i)%FirstNP
          DO WHILE (M .GT. 0)
!!
!! Compute angular momentum
!!
            Pxcm = MOTION(M)%Px - Xcm
            Pycm = MOTION(M)%Py - Ycm
            Pzcm = MOTION(M)%Pz - Zcm
            dVxc = MOTION(M)%Vx - Vxcm
            dVyc = MOTION(M)%Vy - Vycm
            dVzc = MOTION(M)%Vz - Vzcm
            Omx  = Omx - NODE(M)%Mass * (Pycm*dVzc - Pzcm*dVyc)
            Omy  = Omy - NODE(M)%Mass * (Pzcm*dVxc - Pxcm*dVzc)
            Omz  = Omz - NODE(M)%Mass * (Pxcm*dVyc - Pycm*dVxc)
!!
!! Compute inertia tensor B from nodal point masses.
!!
            RSQ = Pxcm*Pxcm + Pycm*Pycm + Pzcm*Pzcm
            Bxx = Bxx + (RSQ - Pxcm*Pxcm) * NODE(M)%Mass
            Byy = Byy + (RSQ - Pycm*Pycm) * NODE(M)%Mass
            Bzz = Bzz + (RSQ - Pzcm*Pzcm) * NODE(M)%Mass
            Bxy = Bxy - Pxcm*Pycm * NODE(M)%Mass
            Bxz = Bxz - Pxcm*Pzcm * NODE(M)%Mass
            Byz = Byz - Pycm*Pzcm * NODE(M)%Mass
            M = NODE(M)%IRB
          ENDDO
          RIGID_BODY(i)%B(1,1) = Bxx
          RIGID_BODY(i)%B(1,2) = Bxy
          RIGID_BODY(i)%B(1,3) = Bxz
          RIGID_BODY(i)%B(2,1) = Bxy
          RIGID_BODY(i)%B(2,2) = Byy
          RIGID_BODY(i)%B(2,3) = Byz
          RIGID_BODY(i)%B(3,1) = Bxz
          RIGID_BODY(i)%B(3,2) = Byz
          RIGID_BODY(i)%B(3,3) = Bzz
!!
!! Invert the inertia tensor B.
!!
          Det = (Bxx*Byy)*Bzz                                                  &
     &        + (Bxy*Byz)*(Bxz+Bxz)                                            &
     &        - (Bxz*Bxz)*Byy                                                  &
     &        - (Byz*Byz)*Bxx                                                  &
     &        - (Bxy*Bxy)*Bzz

          IF (ABS (Det) .LE. 1.0D-25) THEN
            WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'INITIALIZE_RIGID_BODY_MASS.001.02'//                    &
     &          MSGL//'Inertia Tensor For Rigid Body ID:'//MSG1//              &
     &          MSGL//'Is Singular. Do all nodes lie on a line?'               &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
          ENDIF
!!
!! Store initial position, initial velocity, and initial angular velocity
!! of the center of mass.
!!
          RIGID_BODY(i)%Px = Xcm
          RIGID_BODY(i)%Py = Ycm
          RIGID_BODY(i)%Pz = Zcm
          RIGID_BODY(i)%Vx = Vxcm
          RIGID_BODY(i)%Vy = Vycm
          RIGID_BODY(i)%Vz = Vzcm
          RIGID_BODY(i)%Ox = Cxx*Omx + Cxy*Omy + Cxz*Omz
          RIGID_BODY(i)%Oy = Cxy*Omx + Cyy*Omy + Cyz*Omz
          RIGID_BODY(i)%Oz = Cxz*Omx + Cyz*Omy + Czz*Omz
        ENDIF
      ENDDO

      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_RIGID_BODY_MASS.002.00'//                    &
     &          MSGL//'Total Number Of Zero Masses/Inertias:'//MSG1            &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_RIGID_BODY_CONSTRAINTS
!!
!! Copyright (c) by KEY Associates, 17-APR-1991 19:16:13
!!
!! Purpose: Modify the components of the acceleration vector Ax,Ay,Az to
!! effect rigid body motion of the nodal points comprising rigid body
!! domains.
!!
      USE shared_common_data
      USE rigid_body_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0)) :: dU(3),A(3,3),C(3),PRX(3,3),RTX(3,3)

      REAL(KIND(0D0)), PARAMETER :: C2F = (1.0D+0 /   2.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C3F = (1.0D+0 /   6.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C4F = (1.0D+0 /  24.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C5F = (1.0D+0 / 120.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C6F = (1.0D+0 / 720.0D+0)
!!
!! The material in this module draws heavily on the article:
!!
!!      J. L. Synge, "Classical Dynamics," Encylopedia of Physics,
!!      Volume III, Part 1, Springer-Verlag, Berlin, 1960.
!!
!! COMPUTE THE RIGID BODY MOTION OF EACH RIGID BODY DOMAIN.
!! Once the motion has been incremented, the global acceleration arrays are
!! modified so that when the central difference time integration in module
!! SOLVE is invoked the correct positions will be obtained.
!!
      DTaver = 0.5D+0 * (TIMSIM%DTlast + TIMSIM%DTnext)
      D1  = ONE / TIMSIM%DTnext
      D2  = ONE / DTaver
!!
      DO i = 1,NUMRB
!!
!! MOTION OF THE CENTER OF MASS.
!! Retrieve position of the center of mass.
!!
        Xcm = RIGID_BODY(i)%Px + RIGID_BODY(i)%Ux
        Ycm = RIGID_BODY(i)%Py + RIGID_BODY(i)%Uy
        Zcm = RIGID_BODY(i)%Pz + RIGID_BODY(i)%Uz
!!
!! Retrieve velocity of the center of mass.
!!
        Vxcm = RIGID_BODY(i)%Vx
        Vycm = RIGID_BODY(i)%Vy
        Vzcm = RIGID_BODY(i)%Vz
!!
!! Retrieve angular velocity about the center of mass.
!!
        Oxcm = RIGID_BODY(i)%Ox
        Oycm = RIGID_BODY(i)%Oy
        Ozcm = RIGID_BODY(i)%Oz
!!
!! Compute the total force acting on the rigid body.
!!
        Fxcm = 0.0
        Fycm = 0.0
        Fzcm = 0.0
        M = RIGID_BODY(i)%FirstNP
        DO WHILE (M .GT. 0)
          Fxcm = Fxcm + FORCE(M)%Xext - FORCE(M)%Xint
          Fycm = Fycm + FORCE(M)%Yext - FORCE(M)%Yint
          Fzcm = Fzcm + FORCE(M)%Zext - FORCE(M)%Zint
          M = NODE(M)%IRB
        ENDDO
        RIGID_BODY(i)%Force(1) = Fxcm
        RIGID_BODY(i)%Force(2) = Fycm
        RIGID_BODY(i)%Force(3) = Fzcm
!!
!! Compute linear acceleration of the center of mass.
!!
        Qms = ONE / RIGID_BODY(i)%Mass
        Axcm = Fxcm * Qms
        Aycm = Fycm * Qms
        Azcm = Fzcm * Qms
!!
!! Update acceleration of the center of mass
!!
        RIGID_BODY(i)%Ax = Axcm
        RIGID_BODY(i)%Ay = Aycm
        RIGID_BODY(i)%Az = Azcm
!!
!! Update velocity of the center of mass
!!
        RIGID_BODY(i)%Vx = (Vxcm + DTaver * Axcm)
        RIGID_BODY(i)%Vy = (Vycm + DTaver * Aycm)
        RIGID_BODY(i)%Vz = (Vzcm + DTaver * Azcm)
!!
!! Compute incremental displacement of the center of mass.
!!
        Uxcm = TIMSIM%DTnext * (Vxcm + DTaver * Axcm)
        Uycm = TIMSIM%DTnext * (Vycm + DTaver * Aycm)
        Uzcm = TIMSIM%DTnext * (Vzcm + DTaver * Azcm)
!!
!! Update position of the center of mass.
!!
        RIGID_BODY(i)%Ux = RIGID_BODY(i)%Ux + Uxcm
        RIGID_BODY(i)%Uy = RIGID_BODY(i)%Uy + Uycm
        RIGID_BODY(i)%Uz = RIGID_BODY(i)%Uz + Uzcm
!!
!! Compute torque of applied forces about the center of mass.
!! T = Sum(1:N){(Xn-Xcm) x Fn}
!!
        Tx = 0.0
        Ty = 0.0
        Tz = 0.0
        M = RIGID_BODY(i)%FirstNP
        DO WHILE (M .GT. 0)
          Pxcm = MOTION(M)%Px + MOTION(M)%Ux - Xcm
          Pycm = MOTION(M)%Py + MOTION(M)%Uy - Ycm
          Pzcm = MOTION(M)%Pz + MOTION(M)%Uz - Zcm
          Fx = FORCE(M)%Xext - FORCE(M)%Xint
          Fy = FORCE(M)%Yext - FORCE(M)%Yint
          Fz = FORCE(M)%Zext - FORCE(M)%Zint
          Tx = Tx + Pycm*Fz - Pzcm*Fy
          Ty = Ty + Pzcm*Fx - Pxcm*Fz
          Tz = Tz + Pxcm*Fy - Pycm*Fx
          IF (NODE(M)%IRT .GT. 0) THEN
            Tx = Tx + FORCE(NODE(M)%IRT)%Xext - FORCE(NODE(M)%IRT)%Xint
            Ty = Ty + FORCE(NODE(M)%IRT)%Yext - FORCE(NODE(M)%IRT)%Yint
            Tz = Tz + FORCE(NODE(M)%IRT)%Zext - FORCE(NODE(M)%IRT)%Zint
          ENDIF
          M = NODE(M)%IRB
        ENDDO
!!
!! Update torque applied to center of mass.
!!
        RIGID_BODY(i)%Torque(1) = Tx
        RIGID_BODY(i)%Torque(2) = Ty
        RIGID_BODY(i)%Torque(3) = Tz
!!
!! Invert the inertia tensor B.
!!
        Bxx = RIGID_BODY(i)%B(1,1)
        Byy = RIGID_BODY(i)%B(2,2)
        Bzz = RIGID_BODY(i)%B(3,3)
        Bxy = RIGID_BODY(i)%B(1,2)
        Bxz = RIGID_BODY(i)%B(1,3)
        Byz = RIGID_BODY(i)%B(2,3)
        Det = (Bxx*Byy)*Bzz                                                    &
     &      + (Bxy*Byz)*(Bxz+Bxz)                                              &
     &      - (Bxz*Bxz)*Byy                                                    &
     &      - (Byz*Byz)*Bxx                                                    &
     &      - (Bxy*Bxy)*Bzz

        IF (ABS (Det) .LE. 1.0D-25) THEN
          WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_RIGID_BODY_CONSTRAINTS.001.00'//                 &
     &          MSGL//'Inertia Tensor For Rigid Body ID:'//MSG1//              &
     &          MSGL//'Has Become Singular. Should Have Been'//                &
     &          MSGL//'Detected Earlier In INITIALIZE_RIGID_BODY_MASS.'        &
     &          )
        ENDIF

        Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
        Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
        Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
        Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
        Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
        Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
!!
!! Using Euler's equation, compute new angular velocity rate about center
!! of mass; dO/dt = Binv (T - OxB*O)
!!
        Q1 = Bxx*Oxcm + Bxy*Oycm + Bxz*Ozcm
        Q2 = Bxy*Oxcm + Byy*Oycm + Byz*Ozcm
        Q3 = Bxz*Oxcm + Byz*Oycm + Bzz*Ozcm
        Tx = Tx + Q3*Oycm - Q2*Ozcm
        Ty = Ty + Q1*Ozcm - Q3*Oxcm
        Tz = Tz + Q2*Oxcm - Q1*Oycm
!!
!! Compute new angular velocity of the center of mass.
!!
        Oxcm = Oxcm + DTaver * (Cxx*Tx + Cxy*Ty + Cxz*Tz)
        Oycm = Oycm + DTaver * (Cxy*Tx + Cyy*Ty + Cyz*Tz)
        Ozcm = Ozcm + DTaver * (Cxz*Tx + Cyz*Ty + Czz*Tz)
!!
!! Update angular velocity of the center of mass
!!
        RIGID_BODY(i)%Ox = Oxcm
        RIGID_BODY(i)%Oy = Oycm
        RIGID_BODY(i)%Oz = Ozcm
!!
!! RIGID BODY MOTION.
!! Construct rigid body motion at each point in the rigid body domain.
!!
        Qmag = SQRT (Oxcm*Oxcm + Oycm*Oycm + Ozcm*Ozcm)
        IF (Qmag .GT. 1.0D-20) THEN
!!
!! Construct unit vector C aligned with angular velocity "vector."
!!
          C(1) = Oxcm * (ONE / Qmag)
          C(2) = Oycm * (ONE / Qmag)
          C(3) = Ozcm * (ONE / Qmag)
!!
!! Compute angle of rotation and evaluate the trigonometric functions
!! SinPhi = Sin(Phi) and CosPm1 = Cos(Phi) - 1.
!!
          Phi = TIMSIM%DTnext * Qmag
          Ph2 = Phi*Phi
          CosPm1 = -Ph2*(C2F - Ph2*(C4F - C6F*Ph2))
          SinPhi =  Phi*(ONE - Ph2*(C3F - C5F*Ph2))
!!
!! Construct rotation operator RTX. R = (Cos(phi)-1.0)*(I-CC) + Sin(phi)*(C x .)
!! The operator RTX is an "incremental" proper-orthogonal rotation. It is
!! constructed to be an "incremental operator" to avoid round-off.
!!
          DO n = 1,3
            DO m = 1,3
              RTX(m,n) = 0.0
              PRX(m,n) = -C(m) * C(n)
            ENDDO
          ENDDO
          DO n = 1,3
            PRX(n,n) = ONE + PRX(n,n)
          ENDDO
          RTX(1,2) = -(SinPhi*C(3))
          RTX(1,3) =  (SinPhi*C(2))
          RTX(2,1) =  (SinPhi*C(3))
          RTX(2,3) = -(SinPhi*C(1))
          RTX(3,1) = -(SinPhi*C(2))
          RTX(3,2) =  (SinPhi*C(1))
          DO n = 1,3
            DO m = 1,3
              RTX(m,n) = CosPm1*PRX(m,n) + RTX(m,n)
            ENDDO
          ENDDO
!!
!! Rotate and translate each nodal point in the rigid body. The
!! incremental displacement of the "rigid link" with end point M is
!! computed by applying a rigid body rotation RTX to the vector from
!! the center of mass to the point M along with the translation of
!! the center of mass.
!!
          M = RIGID_BODY(i)%FirstNP
          DO WHILE (M .GT. 0)
!!
            DO n = 1,3
              dU(n) = RTX(n,1) * (MOTION(M)%Px + MOTION(M)%Ux - Xcm)           &
     &              + RTX(n,2) * (MOTION(M)%Py + MOTION(M)%Uy - Ycm)           &
     &              + RTX(n,3) * (MOTION(M)%Pz + MOTION(M)%Uz - Zcm)
            ENDDO
!!
!! Back calculate translational nodal accelerations at M.
!!
            MOTION(M)%Ax = ((dU(1) + Uxcm) * D1 - MOTION(M)%Vx) * D2
            MOTION(M)%Ay = ((dU(2) + Uycm) * D1 - MOTION(M)%Vy) * D2
            MOTION(M)%Az = ((dU(3) + Uzcm) * D1 - MOTION(M)%Vz) * D2
!!
!! Back calculate angular nodal accelerations at M.
!!
            IF (NODE(M)%IRT .GT. 0) THEN
              MOTION(NODE(M)%IRT)%Ax = (Oxcm-MOTION(NODE(M)%IRT)%Vx)*D2
              MOTION(NODE(M)%IRT)%Ay = (Oycm-MOTION(NODE(M)%IRT)%Vy)*D2
              MOTION(NODE(M)%IRT)%Az = (Ozcm-MOTION(NODE(M)%IRT)%Vz)*D2
            ENDIF
!!
            M = NODE(M)%IRB
          ENDDO
!!
!! Update components of inertia tensor from time n to n+1 to acccount for
!! rotational change of the configuration wrt the spatial coordinates.
!! The incremental rotation operator RTX needs to be augmented with the
!! identity operator:
!!
          RTX(1,1) = ONE + RTX(1,1)
          RTX(2,2) = ONE + RTX(2,2)
          RTX(3,3) = ONE + RTX(3,3)

          DO n = 1,3
            A(1,n) = RTX(1,1) * RIGID_BODY(i)%B(1,n)                           &
     &             + RTX(1,2) * RIGID_BODY(i)%B(2,n)                           &
     &             + RTX(1,3) * RIGID_BODY(i)%B(3,n)
            A(2,n) = RTX(2,1) * RIGID_BODY(i)%B(1,n)                           &
     &             + RTX(2,2) * RIGID_BODY(i)%B(2,n)                           &
     &             + RTX(2,3) * RIGID_BODY(i)%B(3,n)
            A(3,n) = RTX(3,1) * RIGID_BODY(i)%B(1,n)                           &
     &             + RTX(3,2) * RIGID_BODY(i)%B(2,n)                           &
     &             + RTX(3,3) * RIGID_BODY(i)%B(3,n)
          ENDDO
!!
          RIGID_BODY(i)%B(1,1) =                                               &
     &          A(1,1)*RTX(1,1) + A(1,2)*RTX(1,2) + A(1,3)*RTX(1,3)

          RIGID_BODY(i)%B(1,2) =                                               &
     &          A(1,1)*RTX(2,1) + A(1,2)*RTX(2,2) + A(1,3)*RTX(2,3)

          RIGID_BODY(i)%B(1,3) =                                               &
     &          A(1,1)*RTX(3,1) + A(1,2)*RTX(3,2) + A(1,3)*RTX(3,3)

          RIGID_BODY(i)%B(2,2) =                                               &
     &          A(2,1)*RTX(2,1) + A(2,2)*RTX(2,2) + A(2,3)*RTX(2,3)

          RIGID_BODY(i)%B(2,3) =                                               &
     &          A(2,1)*RTX(3,1) + A(2,2)*RTX(3,2) + A(2,3)*RTX(3,3)

          RIGID_BODY(i)%B(3,3) =                                               &
     &          A(3,1)*RTX(3,1) + A(3,2)*RTX(3,2) + A(3,3)*RTX(3,3)

          RIGID_BODY(i)%B(2,1) = RIGID_BODY(i)%B(1,2)
          RIGID_BODY(i)%B(3,1) = RIGID_BODY(i)%B(1,3)
          RIGID_BODY(i)%B(3,2) = RIGID_BODY(i)%B(2,3)

        ELSE
!!
!! Rotation is too small; perform translation only.
!!
          M = RIGID_BODY(i)%FirstNP
          DO WHILE (M .GT. 0)
!!
!! Back calculate translational nodal accelerations at M.
!!
            MOTION(M)%Ax = (Uxcm * D1 - MOTION(M)%Vx) * D2
            MOTION(M)%Ay = (Uycm * D1 - MOTION(M)%Vy) * D2
            MOTION(M)%Az = (Uzcm * D1 - MOTION(M)%Vz) * D2
!!
            M = NODE(M)%IRB
          ENDDO
        ENDIF
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE INIT_CONCENTRATED_MASSES
!!
!! Copyright (c) by KEY Associates, 4-APR-1992 18:34:40
!!
!! Purpose: Initialize nodal point masses with user "concentrated" nodal
!! point masses. (The inertia component Bkk/3.0 = (Ixx + Iyy + Izz)/3.0
!! is place-holder. The rotational acceleration, velocity and displacement
!! are calculated in a spearate module for nodal points with user-specified
!! concentrated inertia.)
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE force_
      USE nodal_point_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: RTID
!!
      DO i = 1,NUMCM
!!
!! Complete inertia tensor.
!!
        NODAL_POINT_MASS(i)%B(2,1) = NODAL_POINT_MASS(i)%B(1,2)
        NODAL_POINT_MASS(i)%B(3,1) = NODAL_POINT_MASS(i)%B(1,3)
        NODAL_POINT_MASS(i)%B(3,2) = NODAL_POINT_MASS(i)%B(2,3)
!!
!! Put concentrated mass and inertia into glodal mass "vector." The inertia
!! component Bkk/3.0 is used to provide a non-zero mass in case there is only
!! a torsional spring acting at the nodal point and used to give the rotational
!! displacement boundary conditions the chance to get the modified torques
!! about right.
!!
        NPID = NODAL_POINT_MASS(i)%NPID
        NODE(NPID)%Mass = NODAL_POINT_MASS(i)%Mass
        IF (NODE(NPID)%IRT .GT. 0) THEN
          NODE(NODE(NPID)%IRT)%Mass =                                          &
     &          (                                                              &
     &          NODAL_POINT_MASS(i)%B(1,1) +                                   &
     &          NODAL_POINT_MASS(i)%B(2,2) +                                   &
     &          NODAL_POINT_MASS(i)%B(3,3)                                     &
     &          ) / 3.0D+0
        ENDIF
!!
!! Initialize nodal-point-mass data structure with nodal point initial
!! conditions.
!!
        NODAL_POINT_MASS(i)%Pzero(1) = MOTION(NPID)%Px
        NODAL_POINT_MASS(i)%Pzero(2) = MOTION(NPID)%Py
        NODAL_POINT_MASS(i)%Pzero(3) = MOTION(NPID)%Pz
        NODAL_POINT_MASS(i)%Disp(1)  = MOTION(NPID)%Ux
        NODAL_POINT_MASS(i)%Disp(2)  = MOTION(NPID)%Uy
        NODAL_POINT_MASS(i)%Disp(3)  = MOTION(NPID)%Uz
        NODAL_POINT_MASS(i)%Vel(1)   = MOTION(NPID)%Vx
        NODAL_POINT_MASS(i)%Vel(2)   = MOTION(NPID)%Vy
        NODAL_POINT_MASS(i)%Vel(3)   = MOTION(NPID)%Vz
!!
        IF (NODE(NPID)%IRT .GT. 0) THEN
          RTID = NODE(NPID)%IRT
          NODAL_POINT_MASS(i)%Omega(1) = MOTION(RTID)%Vx
          NODAL_POINT_MASS(i)%Omega(2) = MOTION(RTID)%Vy
          NODAL_POINT_MASS(i)%Omega(3) = MOTION(RTID)%Vz
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE UPDATE1_CONCENTRATED_MASSES
!!
!! Copyright (c) by KEY Associates, 4-APR-1992 18:34:40
!!
!! Purpose: Compute proper angular acceleration based on the torque applied
!! to the concentrated inertia.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE force_
      USE nodal_point_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: RTID
!!
      DO i = 1,NUMCM
        NPID = NODAL_POINT_MASS(i)%NPID
        IF (NODE(NPID)%On) THEN
          DTlast = NODE(NPID)%DTlast
          DTnext = NODE(NPID)%DTnext
          DTaver = 0.5D+0 * (DTlast + DTnext)
!!
!! Compute rotational acceleration. First, check to see if this nodal point
!! has rotational degrees of freedom and a user-specified concentrated
!! inertia tensor.
!!
          Bxx = NODAL_POINT_MASS(i)%B(1,1)
          IF (NODE(NPID)%IRT .GT. 0 .AND. Bxx .GT. 0.0) THEN
            RTID = NODE(NPID)%IRT
!!
!! Retrieve angular velocity at time n-1/2 and torque at time n at nodal point.
!!
            Ox = NODAL_POINT_MASS(i)%Omega(1)
            Oy = NODAL_POINT_MASS(i)%Omega(2)
            Oz = NODAL_POINT_MASS(i)%Omega(3)
            Tx = FORCE(RTID)%Xext-FORCE(RTID)%Xint
            Ty = FORCE(RTID)%Yext-FORCE(RTID)%Yint
            Tz = FORCE(RTID)%Zext-FORCE(RTID)%Zint
!!
            Bxx = NODAL_POINT_MASS(i)%B(1,1)
            Byy = NODAL_POINT_MASS(i)%B(2,2)
            Bzz = NODAL_POINT_MASS(i)%B(3,3)
            Bxy = NODAL_POINT_MASS(i)%B(1,2)
            Bxz = NODAL_POINT_MASS(i)%B(1,3)
            Byz = NODAL_POINT_MASS(i)%B(2,3)
!!
!! Invert the inertia tensor B.
!!
            Det = (Bxx*Byy)*Bzz                                                &
     &          + (Bxy*Byz)*(Bxz+Bxz)                                          &
     &          - (Bxz*Bxz)*Byy                                                &
     &          - (Byz*Byz)*Bxx                                                &
     &          - (Bxy*Bxy)*Bzz
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
!!
!! Using Euler's equation, compute angular acceleration at time n about center
!! of mass; dO/dt = Binv (T - OxB*O)
!!
            Q1 = Bxx*Ox + Bxy*Oy + Bxz*Oz
            Q2 = Bxy*Ox + Byy*Oy + Byz*Oz
            Q3 = Bxz*Ox + Byz*Oy + Bzz*Oz
            Tx = Tx + Q3*Oy - Q2*Oz
            Ty = Ty + Q1*Oz - Q3*Ox
            Tz = Tz + Q2*Ox - Q1*Oy
!!
!! Compute new angular acceleration of the center of mass.
!!
            MOTION(RTID)%Ax = (Cxx*Tx + Cxy*Ty + Cxz*Tz)
            MOTION(RTID)%Ay = (Cxy*Tx + Cyy*Ty + Cyz*Tz)
            MOTION(RTID)%Az = (Cxz*Tx + Cyz*Ty + Czz*Tz)
!!
          ENDIF
!!
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE UPDATE2_CONCENTRATED_MASSES
!!
!! Copyright (c) by KEY Associates, 4-APR-1992 18:34:40
!!
!! Purpose: Compute acceleration, velocity and displacement for nodal points
!! with user-specified concentrated masses.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE force_
      USE nodal_point_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: RTID
!!
      DO i = 1,NUMCM
        NPID = NODAL_POINT_MASS(i)%NPID
        IF (NODE(NPID)%On) THEN
          DTlast = NODE(NPID)%DTlast
          DTnext = NODE(NPID)%DTnext
          DTaver = 0.5D+0 * (DTlast + DTnext)
!!
!! For the record, record translational motion.
!!
          NODAL_POINT_MASS(i)%Force(1) = FORCE(NPID)%Xext                      &
     &                                 - FORCE(NPID)%Xint
          NODAL_POINT_MASS(i)%Force(2) = FORCE(NPID)%Yext                      &
     &                                 - FORCE(NPID)%Yint
          NODAL_POINT_MASS(i)%Force(3) = FORCE(NPID)%Zext                      &
     &                                 - FORCE(NPID)%Zint
          NODAL_POINT_MASS(i)%Disp(1)  = MOTION(NPID)%Ux
          NODAL_POINT_MASS(i)%Disp(2)  = MOTION(NPID)%Uy
          NODAL_POINT_MASS(i)%Disp(3)  = MOTION(NPID)%Uz
          NODAL_POINT_MASS(i)%Vel(1)   = MOTION(NPID)%Vx
          NODAL_POINT_MASS(i)%Vel(2)   = MOTION(NPID)%Vy
          NODAL_POINT_MASS(i)%Vel(3)   = MOTION(NPID)%Vz
          NODAL_POINT_MASS(i)%Accel(1) = MOTION(NPID)%Ax
          NODAL_POINT_MASS(i)%Accel(2) = MOTION(NPID)%Ay
          NODAL_POINT_MASS(i)%Accel(3) = MOTION(NPID)%Az
!!
!! Compute rotational acceleration, velocity and displacement for nodal points
!! with user-specified concentrated inertia.
!!
          Bxx = NODAL_POINT_MASS(i)%B(1,1)
          IF (NODE(NPID)%IRT .GT. 0 .AND. Bxx .GT. 0.0) THEN
            RTID = NODE(NPID)%IRT
!!
!! Retrieve angular velocity at time n-1/2 and torque at time n at nodal point.
!! Note: the torque may have been modified by kinematic boundary conditions.
!!
            Ox = NODAL_POINT_MASS(i)%Omega(1)
            Oy = NODAL_POINT_MASS(i)%Omega(2)
            Oz = NODAL_POINT_MASS(i)%Omega(3)
            Tx = FORCE(RTID)%Xext-FORCE(RTID)%Xint
            Ty = FORCE(RTID)%Yext-FORCE(RTID)%Yint
            Tz = FORCE(RTID)%Zext-FORCE(RTID)%Zint
!!
!! Retrive components of the inertia tensor at time n.
!!
            Bxx = NODAL_POINT_MASS(i)%B(1,1)
            Byy = NODAL_POINT_MASS(i)%B(2,2)
            Bzz = NODAL_POINT_MASS(i)%B(3,3)
            Bxy = NODAL_POINT_MASS(i)%B(1,2)
            Bxz = NODAL_POINT_MASS(i)%B(1,3)
            Byz = NODAL_POINT_MASS(i)%B(2,3)
!!
!! Invert the inertia tensor B.
!!
            Det = (Bxx*Byy)*Bzz                                                &
     &          + (Bxy*Byz)*(Bxz+Bxz)                                          &
     &          - (Bxz*Bxz)*Byy                                                &
     &          - (Byz*Byz)*Bxx                                                &
     &          - (Bxy*Bxy)*Bzz
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
!!
!! Using Euler's equation, compute new angular velocity about center
!! of mass; dO/dt = Binv (T - OxB*O)
!!
            Q1 = Bxx*Ox + Bxy*Oy + Bxz*Oz
            Q2 = Bxy*Ox + Byy*Oy + Byz*Oz
            Q3 = Bxz*Ox + Byz*Oy + Bzz*Oz
            Tx = Tx + Q3*Oy - Q2*Oz
            Ty = Ty + Q1*Oz - Q3*Ox
            Tz = Tz + Q2*Ox - Q1*Oy
!!
!! Compute new angular velocity of the center of mass. This calculation will
!! give the correct motion to the concentrated inertia based on the torque.
!! The kinematic boundary conditions will be only approximately correct.
!! However, at each time step they will attempt to compensate for past
!! inaccuracies. The result should be OK.
!!
            Ox = Ox + DTaver * (Cxx*Tx + Cxy*Ty + Cxz*Tz)
            Oy = Oy + DTaver * (Cxy*Tx + Cyy*Ty + Cyz*Tz)
            Oz = Oz + DTaver * (Cxz*Tx + Cyz*Ty + Czz*Tz)
!!
!! Update angular displacement.
!!
            Ux = NODAL_POINT_MASS(i)%Theta(1) + DTnext * Ox
            Uy = NODAL_POINT_MASS(i)%Theta(2) + DTnext * Oy
            Uz = NODAL_POINT_MASS(i)%Theta(3) + DTnext * Oz
!!
!! Update angular velocity and displacement of the center of mass
!!
            MOTION(RTID)%Vx = Ox
            MOTION(RTID)%Vy = Oy
            MOTION(RTID)%Vz = Oz
            MOTION(RTID)%Ux = Ux
            MOTION(RTID)%Uy = Uy
            MOTION(RTID)%Uz = Uz
            NODAL_POINT_MASS(i)%Omega(1)  = Ox
            NODAL_POINT_MASS(i)%Omega(2)  = Oy
            NODAL_POINT_MASS(i)%Omega(3)  = Oz
            NODAL_POINT_MASS(i)%Theta(1)  = Ux
            NODAL_POINT_MASS(i)%Theta(2)  = Uy
            NODAL_POINT_MASS(i)%Theta(3)  = Uz
            NODAL_POINT_MASS(i)%Torque(1) = Tx
            NODAL_POINT_MASS(i)%Torque(2) = Ty
            NODAL_POINT_MASS(i)%Torque(3) = Tz
!!
!! Update components of inertia tensor from time n to n+1 to acccount for
!! rotational change of the congifuration wrt the spatial coordinates.
!!
            CALL ROTATE_INERTIA (Ox,Oy,Oz,DTnext,NODAL_POINT_MASS(i)%B)
!!
          ENDIF
!!
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE ROTATE_INERTIA (Ox,Oy,Oz,DT,B)
!!
!! Copyright (c) by KEY Associates, 4-APR-1992 21:00:17
!!
!! Purpose: Rotate the inertia tensor B through the angle Dt*|Omega|.
!! The components of Omega are (Ox,Oy,Oz).
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN)    :: Ox,Oy,Oz,DT
      REAL(KIND(0D0)), INTENT(INOUT) :: B(3,3)
!!
!! Local variables
      REAL(KIND(0D0)) :: A(3,3),C(3),PRX(3,3),RTX(3,3)
!!
      REAL(KIND(0D0)), PARAMETER :: C2F = (1.0D+0 /   2.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C3F = (1.0D+0 /   6.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C4F = (1.0D+0 /  24.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C5F = (1.0D+0 / 120.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C6F = (1.0D+0 / 720.0D+0)
!!
      Qmag = SQRT (Ox*Ox + Oy*Oy + Oz*Oz)
      IF (Qmag .GT. 1.0D-20) THEN
!!
!! Construct unit vector C aligned with angular velocity "vector."
!!
        C(1) = Ox * (ONE / Qmag)
        C(2) = Oy * (ONE / Qmag)
        C(3) = Oz * (ONE / Qmag)
!!
!! Compute angle of rotation and evaluate the trigonometric functions
!! SinPhi = Sin(Phi) and CosPm1 = Cos(Phi) - 1.
!!
        Phi = DT * Qmag
        Ph2 = Phi*Phi
        CosPm1 = -Ph2*(C2F - Ph2*(C4F - C6F*Ph2))
        SinPhi =  Phi*(ONE - Ph2*(C3F - C5F*Ph2))
!!
!! Construct rotation operator RTX. R = (Cos(phi)-1.0)*(I-CC) + Sin(phi)*(C x .)
!! The operator RTX is an "incremental" proper-orthogonal rotation. It is
!! constructed to be an "incremental operator" to avoid round-off.
!!
        DO n = 1,3
          DO m = 1,3
            RTX(m,n) = 0.0
            PRX(m,n) = -C(m) * C(n)
          ENDDO
        ENDDO
        DO n = 1,3
          PRX(n,n) = ONE + PRX(n,n)
        ENDDO
        RTX(1,2) = -(SinPhi*C(3))
        RTX(1,3) =  (SinPhi*C(2))
        RTX(2,1) =  (SinPhi*C(3))
        RTX(2,3) = -(SinPhi*C(1))
        RTX(3,1) = -(SinPhi*C(2))
        RTX(3,2) =  (SinPhi*C(1))
        DO n = 1,3
          DO m = 1,3
            RTX(m,n) = CosPm1*PRX(m,n) + RTX(m,n)
          ENDDO
        ENDDO
!!
!! Augment incremental rotation operator with identity operator.
!!
        RTX(1,1) = ONE + RTX(1,1)
        RTX(2,2) = ONE + RTX(2,2)
        RTX(3,3) = ONE + RTX(3,3)
!!
!! Rotate inertia tensor B.
!!
        DO n = 1,3
          A(1,n) = RTX(1,1)*B(1,n) + RTX(1,2)*B(2,n) + RTX(1,3)*B(3,n)
          A(2,n) = RTX(2,1)*B(1,n) + RTX(2,2)*B(2,n) + RTX(2,3)*B(3,n)
          A(3,n) = RTX(3,1)*B(1,n) + RTX(3,2)*B(2,n) + RTX(3,3)*B(3,n)
        ENDDO
!!
        B(1,1) = A(1,1)*RTX(1,1) + A(1,2)*RTX(1,2) + A(1,3)*RTX(1,3)
        B(1,2) = A(1,1)*RTX(2,1) + A(1,2)*RTX(2,2) + A(1,3)*RTX(2,3)
        B(1,3) = A(1,1)*RTX(3,1) + A(1,2)*RTX(3,2) + A(1,3)*RTX(3,3)
        B(2,2) = A(2,1)*RTX(2,1) + A(2,2)*RTX(2,2) + A(2,3)*RTX(2,3)
        B(2,3) = A(2,1)*RTX(3,1) + A(2,2)*RTX(3,2) + A(2,3)*RTX(3,3)
        B(3,3) = A(3,1)*RTX(3,1) + A(3,2)*RTX(3,2) + A(3,3)*RTX(3,3)

        B(2,1) = B(1,2)
        B(3,1) = B(1,3)
        B(3,2) = B(2,3)
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_DISPLACEMENT_BC
!!
!! Copyright (c) by KEY Associates, 27-MAR-1991 20:18:00
!!
!! Purpose: Apply displacement B.C. to individual nodal points.
!!
      USE shared_common_data
      USE displacement_bc_
      USE node_set_
      USE motion_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!!      I. Translational Degrees Of Freedom:
!!
!!                Code = 1, Prescribed movement along the X-axis,
!!                          free to move in the Y,Z-plane.
!!
!!                     = 2, prescribed movement along the Y-axis,
!!                          free to move in the X,Z-plane.
!!
!!                     = 3, Prescribed movement along the Z-axis,
!!                          free to move in the X,Y-plane.
!!
!!                     = 4, Prescribed movement along a line defined
!!                          by the unit vector  D with components
!!                          Dx, Dy, and Dz. Free to move in the plane
!!                          at right angles to the unit vector D.
!!
!!                     =10, Free to move along the X-axis,
!!                          No movement permitted in the Y,Z-plane.
!!
!!                     =20, Free to move along the Y-axis,
!!                          No movement permitted in the X,Z-plane.
!!
!!                     =30, Free to move along the Z-axis,
!!                          No movement permitted in the X,Y-plane.
!!
!!                     =40, Free to move along a line defined by the
!!                          unit vector  D with components Dx, Dy, and
!!                          Dz. No movement permitted in the plane at
!!                          right angles to the unit vector D.
!!
!!                     =11, Prescribed movement along the X-axis,
!!                          No movement permitted in the Y,Z-plane.
!!
!!                     =22, prescribed movement along the Y-axis,
!!                          No movement permitted in the X,Z-plane.
!!
!!                     =33, Prescribed movement along the Z-axis,
!!                          No movement permitted in the X,Y-plane.
!!
!!                     =44, Prescribed movement along a line defined
!!                          by the unit vector  D with components
!!                          Dx, Dy, and Dz. No movement permitted
!!                          in the plane at right angles to the unit
!!                          vector D.
!!
!!      II. Rotational Degrees Of Freedom:
!!
!!                Code = 5, Prescribed rotation about the X-axis, free
!!                          to rotate about the axes in the Y,Z-plane.
!!
!!                     = 6, prescribed rotation about the Y-axis, free
!!                          to rotate about the axes in the X,Z-plane.
!!
!!                     = 7, Prescribed rotation about the Z-axis, free
!!                          to rotate about the axes in the X,Y-plane.
!!
!!                     = 8, Prescribed rotation about a line defined
!!                          by the unit vector D with components Dx,
!!                          Dy, and Dz. Free to rotate about the axes
!!                          in the plane at right angles to the unit
!!                          vector D.
!!
!!                     =50, Free to rotate about the X-axis, No rotation
!!                          permitted about the axes in the Y,Z-plane.
!!
!!                     =60, Free to rotate about the Y-axis, No rotation
!!                          permitted about the axes in the X,Z-plane.
!!
!!                     =70, Free to rotate about the Z-axis, No rotation
!!                          permitted about the axes in the X,Y-plane.
!!
!!                     =80, Free to rotate about a line defined by the
!!                          unit vector  D with components Dx, Dy, and
!!                          Dz. No rotation permitted about the axes in
!!                          the plane at right angles to the unit vector
!!                           D.
!!
!!                     =55, Prescribed rotation about the X-axis, No
!!                          rotation permitted about the axes in the
!!                          Y,Z-plane.
!!
!!                     =66, prescribed rotation about the Y-axis, No
!!                          rotation permitted about the axes in the
!!                          X,Z-plane.
!!
!!                     =77, Prescribed rotation about the Z-axis, No
!!                          rotation permitted about the axes in the
!!                          X,Y-plane.
!!
!!                     =88, Prescribed rotation about a line defined by
!!                          the unit vector  D with components Dx, Dy,
!!                          and Dz. No rotation permitted about the axes
!!                          in the plane at right angles to the unit
!!                          vector D.
!!
!!                Kavd = 1, Prescribed acceleration condition.
!!
!!                     = 2, Prescribed velocity condition.
!!
!!                     = 3, Prescribed displacement condition.
!!
!!_
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          Code,                                                          &
     &          HstID,                                                         &
     &          Kavd
      REAL(KIND(0D0))                                                          &
     &          TABLE_LOOK_UP
      LOGICAL                                                                  &
     &          NEXT_NP_ID,                                                    &
     &          Rotation

      REAL(KIND(0D0)), PARAMETER :: DTR = 1.745329252D-2

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Check for rotational BC's on non-shell nodal points.
!!
        DO NDC = 1,NUMDC
          SetID = DISPLACEMENT_BC(NDC)%SetID
          Code  = DISPLACEMENT_BC(NDC)%Code
          Rotation = (Code .GE. 5 .AND. Code .LE. 8) .OR. Code .GE. 50
!!
!! If this is a kinematic boundary condition for a rotational degree of
!! freedom, examine all nodes in the set to insure they are nodes with
!! rotational degrees of freedom.
!!
          IF (Rotation) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (NODE(N)%IRT .EQ. 0) THEN
                WRITE (MSG1,'(I8)') DISPLACEMENT_BC(NDC)%DBCID
                WRITE (MSG2,'(I8)') NODE(N)%ID
                CALL USER_MESSAGE                                              &
     &      (                                                                  &
     &      MSGL//'FATAL'//                                                    &
     &      MSGL//'IMPOSE_DISPLACEMENT_BC.001.00'//                            &
     &      MSGL//'DISPBC Input Record ID:'//MSG1//                            &
     &      MSGL//'Specifies A Rotational BC On Nodal Point:'//MSG2//          &
     &      MSGL//'Which Does Not Have Rotational Degrees Of Freedom.'         &
     &      )
              ENDIF
            ENDDO
          ENDIF
        ENDDO
!!
!! Bring acceleration/velocity/displacement flag Kavd into range.
!!
        DO N = 1,NUMDC
          DISPLACEMENT_BC(N)%Kavd =                                            &
     &      MAX (1,MIN (3,DISPLACEMENT_BC(N)%Kavd))
        ENDDO
!!
!! Normalize direction vector A.
!!
        DO N = 1,NUMDC
          Amag = SQRT                                                          &
     &          (                                                              &
     &          DISPLACEMENT_BC(N)%Ax * DISPLACEMENT_BC(N)%Ax +                &
     &          DISPLACEMENT_BC(N)%Ay * DISPLACEMENT_BC(N)%Ay +                &
     &          DISPLACEMENT_BC(N)%Az * DISPLACEMENT_BC(N)%Az                  &
     &          )
          IF (Amag .EQ. 0.0) Amag = ONE
          DISPLACEMENT_BC(N)%Ax = DISPLACEMENT_BC(N)%Ax * (ONE / Amag)
          DISPLACEMENT_BC(N)%Ay = DISPLACEMENT_BC(N)%Ay * (ONE / Amag)
          DISPLACEMENT_BC(N)%Az = DISPLACEMENT_BC(N)%Az * (ONE / Amag)
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all kinematic boundary condition specifications.
!! DISPLACEMENT_BC,DBCID,SetID,Code,HstID,Kavd,Scale,Ax,Ay,Az
!!
      DO 200 NDC = 1,NUMDC
        SetID = DISPLACEMENT_BC(NDC)%SetID
        Code  = DISPLACEMENT_BC(NDC)%Code
        HstID = DISPLACEMENT_BC(NDC)%HstID
        Kavd  = DISPLACEMENT_BC(NDC)%Kavd
        Scale = DISPLACEMENT_BC(NDC)%Scale
        Dx    = DISPLACEMENT_BC(NDC)%Ax
        Dy    = DISPLACEMENT_BC(NDC)%Ay
        Dz    = DISPLACEMENT_BC(NDC)%Az
        Rotation = (Code .GE. 5 .AND. Code .LE. 8) .OR. Code .GE. 50
!!
!! Process all nodes in set.
!!
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Apply kinematic boundary condition only to those nodal points that are
!! being integrated.
!!
          IF (NODE(N)%On) THEN
!!
!! Retrieve time increments.
!!
            DTlast = NODE(N)%DTlast
            DTnext = NODE(N)%DTnext
            Dniv = ONE / DTnext
            Daiv = 2.0D+0 / (DTlast + DTnext)
!!
!! Evaluate history function.
!!
            IF (Kavd .EQ. 1) THEN
              Tx = TIMSIM%Total
            ELSE IF (Kavd .EQ. 2) THEN
              Tx = TIMSIM%Total + 0.5D+0 * DTnext
            ELSE IF (Kavd .EQ. 3) THEN
              Tx = TIMSIM%Total + DTnext
            ENDIF
            Ft = Scale * TABLE_LOOK_UP (HstID,Tx)
!!
!! Discriminate between translations and rotations. Convert degrees into
!! radians and locate rotational degree of freedom.
!!
            IF (Rotation) THEN
              Ft = DTR * Ft
              N = NODE(N)%IRT
            ENDIF
!!
!! Apply constraint.
!!
            IF (Code .EQ. 1 .OR. Code .EQ. 5) THEN
              IF (Kavd .EQ. 1) THEN
                Pa = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Pa = (Ft-MOTION(N)%Vx)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Pa = ((Ft-MOTION(N)%Ux)*Dniv-MOTION(N)%Vx)*Daiv
              ENDIF
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext-(MOTION(N)%Ax-Pa)*NODE(N)%Mass
              MOTION(N)%Ax = Pa
            ELSE IF (Code .EQ. 2 .OR. Code .EQ. 6) THEN
              IF (Kavd .EQ. 1) THEN
                Pa = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Pa = (Ft-MOTION(N)%Vy)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Pa = ((Ft-MOTION(N)%Uy)*Dniv-MOTION(N)%Vy)*Daiv
              ENDIF
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext-(MOTION(N)%Ay-Pa)*NODE(N)%Mass
              MOTION(N)%Ay = Pa
            ELSE IF (Code .EQ. 3 .OR. Code .EQ. 7) THEN
              IF (Kavd .EQ. 1) THEN
                Pa = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Pa = (Ft-MOTION(N)%Vz)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Pa = ((Ft-MOTION(N)%Uz)*Dniv-MOTION(N)%Vz)*Daiv
              ENDIF
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pa)*NODE(N)%Mass
              MOTION(N)%Az = Pa
            ELSE IF (Code .EQ. 4 .OR. Code .EQ. 8) THEN
              IF (Kavd .EQ. 1) THEN
                DA = Dx*MOTION(N)%Ax+Dy*MOTION(N)%Ay+Dz*MOTION(N)%Az
                Pa = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                DA = Dx*MOTION(N)%Ax+Dy*MOTION(N)%Ay+Dz*MOTION(N)%Az
                DV = Dx*MOTION(N)%Vx+Dy*MOTION(N)%Vy+Dz*MOTION(N)%Vz
                Pa = (Ft-DV)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                DA = Dx*MOTION(N)%Ax+Dy*MOTION(N)%Ay+Dz*MOTION(N)%Az
                DU = Dx*MOTION(N)%Ux+Dy*MOTION(N)%Uy+Dz*MOTION(N)%Uz
                DV = Dx*MOTION(N)%Vx+Dy*MOTION(N)%Vy+Dz*MOTION(N)%Vz
                Pa = ((Ft-DU)*Dniv-DV)*Daiv
              ENDIF
              FORCE(N)%Xext = FORCE(N)%Xext-(DA-Pa)*Dx*NODE(N)%Mass
              FORCE(N)%Yext = FORCE(N)%Yext-(DA-Pa)*Dy*NODE(N)%Mass
              FORCE(N)%Zext = FORCE(N)%Zext-(DA-Pa)*Dz*NODE(N)%Mass
              MOTION(N)%Ax = MOTION(N)%Ax-(DA-Pa)*Dx
              MOTION(N)%Ay = MOTION(N)%Ay-(DA-Pa)*Dy
              MOTION(N)%Az = MOTION(N)%Az-(DA-Pa)*Dz
            ELSE IF (Code .EQ. 10 .OR. Code .EQ. 50) THEN
              FORCE(N)%Yext = FORCE(N)%Yext-MOTION(N)%Ay*NODE(N)%Mass
              FORCE(N)%Zext = FORCE(N)%Zext-MOTION(N)%Az*NODE(N)%Mass
              MOTION(N)%Uy = 0.0
              MOTION(N)%Uz = 0.0
              MOTION(N)%Vy = 0.0
              MOTION(N)%Vz = 0.0
              MOTION(N)%Ay = 0.0
              MOTION(N)%Az = 0.0
            ELSE IF (Code .EQ. 20 .OR. Code .EQ. 60) THEN
              FORCE(N)%Xext = FORCE(N)%Xext-MOTION(N)%Ax*NODE(N)%Mass
              FORCE(N)%Zext = FORCE(N)%Zext-MOTION(N)%Az*NODE(N)%Mass
              MOTION(N)%Ux = 0.0
              MOTION(N)%Uz = 0.0
              MOTION(N)%Vx = 0.0
              MOTION(N)%Vz = 0.0
              MOTION(N)%Ax = 0.0
              MOTION(N)%Az = 0.0
            ELSE IF (Code .EQ. 30 .OR. Code .EQ. 70) THEN
              FORCE(N)%Xext = FORCE(N)%Xext-MOTION(N)%Ax*NODE(N)%Mass
              FORCE(N)%Yext = FORCE(N)%Yext-MOTION(N)%Ay*NODE(N)%Mass
              MOTION(N)%Ux = 0.0
              MOTION(N)%Uy = 0.0
              MOTION(N)%Vx = 0.0
              MOTION(N)%Vy = 0.0
              MOTION(N)%Ax = 0.0
              MOTION(N)%Ay = 0.0
            ELSE IF (Code .EQ. 40 .OR. Code .EQ. 80) THEN
              Pa = Dx*MOTION(N)%Ax + Dy*MOTION(N)%Ay + Dz*MOTION(N)%Az
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext-(MOTION(N)%Ax-Pa*Dx)*NODE(N)%Mass
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext-(MOTION(N)%Ay-Pa*Dy)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pa*Dz)*NODE(N)%Mass
              MOTION(N)%Ax = Pa*Dx
              MOTION(N)%Ay = Pa*Dy
              MOTION(N)%Az = Pa*Dz
            ELSE
              IF (Code .EQ. 11 .OR. Code .EQ. 55) THEN
                Fx = Ft
                Fy = 0.0
                Fz = 0.0
              ELSE IF (Code .EQ. 22 .OR. Code .EQ. 66) THEN
                Fx = 0.0
                Fy = Ft
                Fz = 0.0
              ELSE IF (Code .EQ. 33 .OR. Code .EQ. 77) THEN
                Fx = 0.0
                Fy = 0.0
                Fz = Ft
              ELSE IF (Code .EQ. 44 .OR. Code .EQ. 88) THEN
                Fx = Dx*Ft
                Fy = Dy*Ft
                Fz = Dz*Ft
              ELSE
                WRITE (MSG1,'(I8)') DISPLACEMENT_BC(NDC)%DBCID
                WRITE (MSG2,'(I8)') Code
                CALL USER_MESSAGE                                              &
     &            (                                                            &
     &            MSGL//'FATAL'//                                              &
     &            MSGL//'IMPOSE_DISPLACEMENT_BC.002.00'//                      &
     &            MSGL//'DISPBC Input Record ID:'//MSG1//                      &
     &            MSGL//'Contains Unknown BC Code:'//MSG2                      &
     &            )
              ENDIF
              IF (Kavd .EQ. 1) THEN
                Px = Fx
                Py = Fy
                Pz = Fz
              ELSE IF (Kavd .EQ. 2) THEN
                Px = (Fx - MOTION(N)%Vx)*Daiv
                Py = (Fy - MOTION(N)%Vy)*Daiv
                Pz = (Fz - MOTION(N)%Vz)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Px = ((Fx - MOTION(N)%Ux)*Dniv-MOTION(N)%Vx)*Daiv
                Py = ((Fy - MOTION(N)%Uy)*Dniv-MOTION(N)%Vy)*Daiv
                Pz = ((Fz - MOTION(N)%Uz)*Dniv-MOTION(N)%Vz)*Daiv
              ENDIF
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext - (MOTION(N)%Ax-Px)*NODE(N)%Mass
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext - (MOTION(N)%Ay-Py)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext - (MOTION(N)%Az-Pz)*NODE(N)%Mass
              MOTION(N)%Ax = Px
              MOTION(N)%Ay = Py
              MOTION(N)%Az = Pz
            ENDIF
          ENDIF
        ENDDO
 200  ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_TIED_BC
!!
!! Copyright (c) by KEY Associates, 27-MAR-1991 20:18:00
!!
!! Purpose: Apply tied B.C. to individual nodal point sets.
!!
      USE shared_common_data
      USE tied_bc_
      USE node_set_
      USE motion_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!!      I. Translational Degrees Of Freedom:
!!
!!                Code = 1, Tied movement along the X-axis,
!!                          free to move in the Y,Z-plane.
!!
!!                     = 2, Tied movement along the Y-axis,
!!                          free to move in the X,Z-plane.
!!
!!                     = 3, Tied movement along the Z-axis,
!!                          free to move in the X,Y-plane.
!!
!!                     = 4, Tied movement along a line defined
!!                          by the unit vector  D with components
!!                          Dx, Dy, and Dz. Free to move in the plane
!!                          at right angles to the unit vector D.
!!
!!                     =10, Free to move along the X-axis,
!!                          Tied movement in the Y,Z-plane.
!!
!!                     =20, Free to move along the Y-axis,
!!                          Tied movement in the X,Z-plane.
!!
!!                     =30, Free to move along the Z-axis,
!!                          Tied movement in the X,Y-plane.
!!
!!                     =40, Free to move along a line defined by the
!!                          unit vector  D with components Dx, Dy, and Dz.
!!                          Tied movement in the plane at right angles to
!!                          the unit vector D.
!!
!!                     =11, Tied movement along the X-axis,
!!                          Tied movement in the Y,Z-plane.
!!
!!                     =22, Tied movement along the Y-axis,
!!                          Tied movement in the X,Z-plane.
!!
!!                     =33, Tied movement along the Z-axis,
!!                          Tied movement in the X,Y-plane.
!!
!!                     =44, Tied movement along a line defined by the
!!                          unit vector  D with components Dx, Dy, and Dz.
!!                          Tied movement in the plane at right angles to
!!                          the unit vector D.
!!
!!      II. Rotational Degrees Of Freedom:
!!
!!                Code = 5, Tied rotation about the X-axis, free
!!                          to rotate about the axes in the Y,Z-plane.
!!
!!                     = 6, Tied rotation about the Y-axis, free
!!                          to rotate about the axes in the X,Z-plane.
!!
!!                     = 7, Tied rotation about the Z-axis, free
!!                          to rotate about the axes in the X,Y-plane.
!!
!!                     = 8, Tied rotation about a line defined
!!                          by the unit vector D with components Dx,
!!                          Dy, and Dz. Free to rotate about the axes
!!                          in the plane at right angles to the unit
!!                          vector D.
!!
!!                     =50, Free to rotate about the X-axis. Tied
!!                          rotation about the axes in the Y,Z-plane.
!!
!!                     =60, Free to rotate about the Y-axis, Tied
!!                          rotation about the axes in the X,Z-plane.
!!
!!                     =70, Free to rotate about the Z-axis, Tied
!!                          rotation about the axes in the X,Y-plane.
!!
!!                     =80, Free to rotate about a line defined by the
!!                          unit vector  D with components Dx, Dy, and
!!                          Dz. No rotation permitted about the axes in
!!                          the plane at right angles to the unit vector D.
!!
!!                     =55, Tied rotation about the X-axis,
!!                          Tied rotation about the axes in the Y,Z-plane.
!!
!!                     =66, Tied rotation about the Y-axis.
!!                          Tied rotation about the axes in the X,Z-plane.
!!
!!                     =77, Tied rotation about the Z-axis.
!!                          Tied rotation about the axes in the X,Y-plane.
!!
!!                     =88, Tied rotation about a line defined by
!!                          the unit vector  D with components Dx, Dy,
!!                          and Dz. Tied rotation about the axes
!!                          in the plane at right angles to the unit
!!                          vector D.
!!
!!_
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          Code
      LOGICAL                                                                  &
     &          NEXT_NP_ID,                                                    &
     &          Rotation,                                                      &
     &          Constrained

      REAL(KIND(0D0)), PARAMETER :: DTR = 1.745329252D-2

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Check for rotational BC's on non-shell nodal points.
!!
        DO NTC = 1,NUMTC
          SetID = TIED_BC(NTC)%SetID
          Code  = TIED_BC(NTC)%Code
          Rotation = (Code .GE. 5 .AND. Code .LE. 8) .OR. Code .GE. 50
!!
!! If this is a tied boundary condition for a rotational degree of
!! freedom, examine all nodes in the set to insure they are nodes with
!! rotational degrees of freedom.
!!
          IF (Rotation) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (NODE(N)%IRT .EQ. 0) THEN
                WRITE (MSG1,'(I8)') TIED_BC(NTC)%TBCID
                WRITE (MSG2,'(I8)') NODE(N)%ID
                CALL USER_MESSAGE                                              &
     &   (                                                                     &
     &   MSGL//'FATAL'//                                                       &
     &   MSGL//'IMPOSE_TIED_BC.001.00'//                                       &
     &   MSGL//'TIEDBC Input Record ID:'//MSG1//                               &
     &   MSGL//'Specifies A Rotational BC On Nodal Point:'//MSG2//             &
     &   MSGL//'Which Does Not Have Rotational Degrees Of Freedom.'            &
     &   )
              ENDIF
            ENDDO
          ENDIF
        ENDDO
!!
!! Normalize direction vector A.
!!
        DO NTC = 1,NUMTC
          Amag = SQRT                                                          &
     &          (                                                              &
     &          TIED_BC(NTC)%Ax * TIED_BC(NTC)%Ax +                            &
     &          TIED_BC(NTC)%Ay * TIED_BC(NTC)%Ay +                            &
     &          TIED_BC(NTC)%Az * TIED_BC(NTC)%Az                              &
     &          )
          IF (Amag .EQ. 0.0) Amag = ONE
          TIED_BC(NTC)%Ax = TIED_BC(NTC)%Ax * (ONE / Amag)
          TIED_BC(NTC)%Ay = TIED_BC(NTC)%Ay * (ONE / Amag)
          TIED_BC(NTC)%Az = TIED_BC(NTC)%Az * (ONE / Amag)
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all tied boundary condition specifications.
!! TIED_BC,TBCID,SetID,Code,Fmax,Ax,Ay,Az
!!
      DO 200 NTC = 1,NUMTC
        SetID = TIED_BC(NTC)%SetID
        Code  = TIED_BC(NTC)%Code
        Fmax  = TIED_BC(NTC)%Fmax
        Dx    = TIED_BC(NTC)%Ax
        Dy    = TIED_BC(NTC)%Ay
        Dz    = TIED_BC(NTC)%Az
        Rotation = (Code .GE. 5 .AND. Code .LE. 8) .OR. Code .GE. 50
!!
!! If "Fmax" is negative, the tied constraint has "broken" and is no longer
!! operating.
!!
        IF (Fmax .LT. 0.0) GO TO 200
!!
!! For all nodes in the set, determine total force(torque) and mass(inertia).
!!
        N = 0
        Fx = 0.0
        Fy = 0.0
        Fz = 0.0
        Amass = 0.0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Discriminate between translations and rotations.
!!
          IF (Rotation) N = NODE(N)%IRT
!!
!! Add up forces currently acting on nodes and sum nodal point masses.
!!
          Fx = Fx + FORCE(N)%Xext - FORCE(N)%Xint
          Fy = Fy + FORCE(N)%Yext - FORCE(N)%Yint
          Fz = Fz + FORCE(N)%Zext - FORCE(N)%Zint
          Amass = Amass + NODE(N)%Mass
!!
        ENDDO
!!
!! Apply constraint.
!!
        IF (Code .EQ. 1 .OR. Code .EQ. 5) THEN
          Pa = Fx / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              PaMass = ABS((MOTION(N)%Ax-Pa)*NODE(N)%Mass)
              Rmax = MAX (Rmax,PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext-(MOTION(N)%Ax-Pa)*NODE(N)%Mass
              MOTION(N)%Ax = Pa
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 2 .OR. Code .EQ. 6) THEN
          Pa = Fy / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              PaMass = ABS((MOTION(N)%Ay-Pa)*NODE(N)%Mass)
              Rmax = MAX (Rmax,PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext-(MOTION(N)%Ay-Pa)*NODE(N)%Mass
              MOTION(N)%Ay = Pa
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 3 .OR. Code .EQ. 7) THEN
          Pa = Fz / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              PaMass = ABS((MOTION(N)%Az-Pa)*NODE(N)%Mass)
              Rmax = MAX (Rmax,PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pa)*NODE(N)%Mass
              MOTION(N)%Az = Pa
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 4 .OR. Code .EQ. 8) THEN
          Pa = (Dx*Fx + Dy*Fy + Dz*Fz) / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              DA = Dx*MOTION(N)%Ax + Dy*MOTION(N)%Ay + Dz*MOTION(N)%Az
              PaMass = ABS((DA-Pa)*NODE(N)%Mass)
              Rmax = MAX (Rmax,PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              DA = Dx*MOTION(N)%Ax + Dy*MOTION(N)%Ay + Dz*MOTION(N)%Az
              FORCE(N)%Xext = FORCE(N)%Xext-(DA-Pa)*Dx*NODE(N)%Mass
              FORCE(N)%Yext = FORCE(N)%Yext-(DA-Pa)*Dy*NODE(N)%Mass
              FORCE(N)%Zext = FORCE(N)%Zext-(DA-Pa)*Dz*NODE(N)%Mass
              MOTION(N)%Ax = MOTION(N)%Ax-(DA-Pa)*Dx
              MOTION(N)%Ay = MOTION(N)%Ay-(DA-Pa)*Dy
              MOTION(N)%Az = MOTION(N)%Az-(DA-Pa)*Dz
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 10 .OR. Code .EQ. 50) THEN
          Py = Fy / Amass
          Pz = Fz / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              Pa = SQRT                                                        &
     &          (                                                              &
     &          (MOTION(N)%Ay-Py)*(MOTION(N)%Ay-Py) +                          &
     &          (MOTION(N)%Az-Pz)*(MOTION(N)%Az-Pz)                            &
     &          )
              PaMass = Pa*NODE(N)%Mass
              Rmax = MAX (Rmax, PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext-(MOTION(N)%Ay-Py)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pz)*NODE(N)%Mass
              MOTION(N)%Ay = Py
              MOTION(N)%Az = Pz
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 20 .OR. Code .EQ. 60) THEN
          Px = Fx / Amass
          Pz = Fz / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              Pa = SQRT                                                        &
     &          (                                                              &
     &          (MOTION(N)%Ax-Px)*(MOTION(N)%Ax-Px) +                          &
     &          (MOTION(N)%Az-Pz)*(MOTION(N)%Az-Pz)                            &
     &          )
              PaMass = Pa*NODE(N)%Mass
              Rmax = MAX (Rmax, PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext-(MOTION(N)%Ax-Px)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pz)*NODE(N)%Mass
              MOTION(N)%Ax = Px
              MOTION(N)%Az = Pz
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 30 .OR. Code .EQ. 70) THEN
          Px = Fx / Amass
          Py = Fy / Amass
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              Pa = SQRT                                                        &
     &          (                                                              &
     &          (MOTION(N)%Ax-Px)*(MOTION(N)%Ax-Px) +                          &
     &          (MOTION(N)%Ay-Py)*(MOTION(N)%Ay-Py)                            &
     &          )
              PaMass = Pa*NODE(N)%Mass
              Rmax = MAX (Rmax, PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Xext =                                                  &
     &        FORCE(N)%Xext-(MOTION(N)%Ax-Px)*NODE(N)%Mass
              FORCE(N)%Yext =                                                  &
     &        FORCE(N)%Yext-(MOTION(N)%Ay-Py)*NODE(N)%Mass
              MOTION(N)%Ax = Px
              MOTION(N)%Ay = Py
            ENDDO
          ENDIF
        ELSE IF (Code .EQ. 40 .OR. Code .EQ. 80) THEN
          Px = Fx / Amass
          Py = Fy / Amass
          Pz = Fz / Amass
          Pa = Dx*Px + Dy*Py + Dz*Pz
          Px = Px - Pa*Dx
          Py = Py - Pa*Dy
          Pz = Pz - Pa*Dz
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              Pa = Dx*MOTION(N)%Ax + Dy*MOTION(N)%Ay + Dz*MOTION(N)%Az
              Pxn = Px + Pa*Dx
              Pyn = Py + Pa*Dy
              Pzn = Pz + Pa*Dz
              Pa = SQRT                                                        &
     &          (                                                              &
     &          (MOTION(N)%Ax-Pxn)*(MOTION(N)%Ax-Pxn) +                        &
     &          (MOTION(N)%Ay-Pyn)*(MOTION(N)%Ay-Pyn) +                        &
     &          (MOTION(N)%Az-Pzn)*(MOTION(N)%Az-Pzn)                          &
     &          )
              PaMass = Pa*NODE(N)%Mass
              Rmax = MAX (Rmax, PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              Pa = Dx*MOTION(N)%Ax + Dy*MOTION(N)%Ay + Dz*MOTION(N)%Az
              Pxn = Px + Pa*Dx
              Pyn = Py + Pa*Dy
              Pzn = Pz + Pa*Dz
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext-(MOTION(N)%Ax-Pxn)*NODE(N)%Mass
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext-(MOTION(N)%Ay-Pyn)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext-(MOTION(N)%Az-Pzn)*NODE(N)%Mass
              MOTION(N)%Ax = Pxn
              MOTION(N)%Ay = Pyn
              MOTION(N)%Az = Pzn
            ENDDO
          ENDIF
        ELSE
          IF (Code .EQ. 11 .OR. Code .EQ. 55) THEN
            Px = Fx / Amass
            Py = Fy / Amass
            Pz = Fz / Amass
          ELSE IF (Code .EQ. 22 .OR. Code .EQ. 66) THEN
            Px = Fx / Amass
            Py = Fy / Amass
            Pz = Fz / Amass
          ELSE IF (Code .EQ. 33 .OR. Code .EQ. 77) THEN
            Px = Fx / Amass
            Py = Fy / Amass
            Pz = Fz / Amass
          ELSE IF (Code .EQ. 44 .OR. Code .EQ. 88) THEN
            Px = Fx / Amass
            Py = Fy / Amass
            Pz = Fz / Amass
          ELSE
            WRITE (MSG1,'(I8)') TIED_BC(NTC)%TBCID
            WRITE (MSG2,'(I8)') Code
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_TIED_BC.002.00'//                                &
     &          MSGL//'TIEDBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown BC Code:'//MSG2                        &
     &          )
          ENDIF
          IF (Fmax .GT. 0.0) THEN
            N = 0
            Rmax = 0.0
            DO WHILE (NEXT_NP_ID(SetID,N))
              Pa = SQRT                                                        &
     &          (                                                              &
     &          (MOTION(N)%Ax-Px)*(MOTION(N)%Ax-Px) +                          &
     &          (MOTION(N)%Ay-Py)*(MOTION(N)%Ay-Py) +                          &
     &          (MOTION(N)%Az-Pz)*(MOTION(N)%Az-Pz)                            &
     &          )
              PaMass = Pa*NODE(N)%Mass
              Rmax = MAX (Rmax, PaMass)
            ENDDO
            Constrained = Rmax .LT. Fmax
            IF (.NOT.Constrained) Fmax = -ONE
          ELSE
            Constrained = .TRUE.
          ENDIF
          IF (Constrained) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (Rotation)  N = NODE(N)%IRT
              FORCE(N)%Xext =                                                  &
     &          FORCE(N)%Xext - (MOTION(N)%Ax-Px)*NODE(N)%Mass
              FORCE(N)%Yext =                                                  &
     &          FORCE(N)%Yext - (MOTION(N)%Ay-Py)*NODE(N)%Mass
              FORCE(N)%Zext =                                                  &
     &          FORCE(N)%Zext - (MOTION(N)%Az-Pz)*NODE(N)%Mass
              MOTION(N)%Ax = Px
              MOTION(N)%Ay = Py
              MOTION(N)%Az = Pz
            ENDDO
          ENDIF
        ENDIF
 200    ENDDO
!!
      RETURN
      END
!!_
        SUBROUTINE BUILD_EFFECTIVE_MASS
!!
!! Copyright (c) by KEY Associates; 20-OCT-1995 19:37:50.00
!!
!! Purpose: Build effective masses for the nodal points between which con-
!! strained midside nodal points lie by distrubuting the mass of the midside
!! nodes to the end points.
!!
      USE shared_common_data
      USE node_
      USE constrained_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Use NODE(*)%Minv as a scratch location for accumulation; initialize it
!! with the mass at the nodes.
!!
      DO i = 1,NUMRT
        NODE(i)%Minv = NODE(i)%Mass
      ENDDO
!!
!! Add in the mass from the constrained midside nodes.
!!
      DO i = 1,NUMNC
        IF (CONSTRAINED_NODE(i)%CNID .GT. 0) THEN
          IMS = CONSTRAINED_NODE(i)%CNID
          ID1 = CONSTRAINED_NODE(i)%NPID(1)
          ID2 = CONSTRAINED_NODE(i)%NPID(2)
          Wg1 = CONSTRAINED_NODE(i)%Weight
          Wg2 = ONE - Wg1
!!
!! Transfer mass of midside nodal point to end points. (Note that
!! the midside nodal point is NOT assumed to be half way between
!! the end nodes.)
!!
          NODE(ID1)%Minv = NODE(ID1)%Minv + Wg1 * NODE(IMS)%Mass
          NODE(ID2)%Minv = NODE(ID2)%Minv + Wg2 * NODE(IMS)%Mass
!!
!! Check for a constrained node with rotational degrees of freedom.
!!
          IF (NODE(IMS)%IRT .GT. 0) THEN
            IMS = NODE(IMS)%IRT
            ID1 = NODE(ID1)%IRT
            ID2 = NODE(ID2)%IRT
!!
!! Transfer inertia components of midside nodal point to end points. (Note
!! that the midside nodal point is NOT assumed to be half way between the
!! end nodes.)
!!
            NODE(ID1)%Minv = NODE(ID1)%Minv + Wg1 * NODE(IMS)%Mass
            NODE(ID2)%Minv = NODE(ID2)%Minv + Wg2 * NODE(IMS)%Mass
          ENDIF
!!
        ENDIF
      ENDDO
!!
!! Save effective mass with the constrained node data structure.
!!
      DO i = 1,NUMNC
        IF (CONSTRAINED_NODE(i)%CNID .GT. 0) THEN
          IMS = CONSTRAINED_NODE(i)%CNID
          ID1 = CONSTRAINED_NODE(i)%NPID(1)
          ID2 = CONSTRAINED_NODE(i)%NPID(2)
          CONSTRAINED_NODE(i)%TMass(1) = NODE(ID1)%Minv
          CONSTRAINED_NODE(i)%TMass(2) = NODE(ID2)%Minv
          IF (NODE(IMS)%IRT .GT. 0) THEN
            ID1 = NODE(ID1)%IRT
            ID2 = NODE(ID2)%IRT
            CONSTRAINED_NODE(i)%RMass(1) = NODE(ID1)%Minv
            CONSTRAINED_NODE(i)%RMass(2) = NODE(ID2)%Minv
          ENDIF
        ENDIF
      ENDDO
!!
!! Restore inverse mass values.
!!
      DO i = 1,NUMRT
        IF (NODE(i)%Mass .GT. 0.0) NODE(i)%Minv = ONE / NODE(i)%Mass
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_MIDSIDE_CONSTRAINTS1
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 20:58:25.00
!!
!! Purpose: For those nodal points influenced by the presences of constraints
!! on midside nodal points, compute modifications in external forces and
!! accelerations.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE force_
      USE constrained_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! ### FIRST, MODIFY ALL BOUNDING NODES
!! Loop over all constraints.
!!
      DO i = 1,NUMNC
!!
        IF (CONSTRAINED_NODE(i)%CNID .GT. 0) THEN
          IMS = CONSTRAINED_NODE(i)%CNID
          ID1 = CONSTRAINED_NODE(i)%NPID(1)
          ID2 = CONSTRAINED_NODE(i)%NPID(2)
          Wg1 = CONSTRAINED_NODE(i)%Weight
          Wg2 = ONE - Wg1
!!
!! Apply kinematic boundary condition only to those nodal points that are
!! being integrated.
!!
          IF (NODE(IMS)%On) THEN
!!
            Fmsx = Wg1 * (FORCE(IMS)%Xext - FORCE(IMS)%Xint)
            Fmsy = Wg1 * (FORCE(IMS)%Yext - FORCE(IMS)%Yint)
            Fmsz = Wg1 * (FORCE(IMS)%Zext - FORCE(IMS)%Zint)
!!
!! Recompute acceleration of the first end point.
!!
            Qinv = ONE / CONSTRAINED_NODE(i)%TMass(1)
            MOTION(ID1)%Ax =                                                   &
     &        (FORCE(ID1)%Xext - FORCE(ID1)%Xint + Fmsx) * Qinv
            MOTION(ID1)%Ay =                                                   &
     &        (FORCE(ID1)%Yext - FORCE(ID1)%Yint + Fmsy) * Qinv
            MOTION(ID1)%Az =                                                   &
     &        (FORCE(ID1)%Zext - FORCE(ID1)%Zint + Fmsz) * Qinv
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
            FORCE(ID1)%Xext =                                                  &
     &        NODE(ID1)%Mass * MOTION(ID1)%Ax + FORCE(ID1)%Xint
            FORCE(ID1)%Yext =                                                  &
     &        NODE(ID1)%Mass * MOTION(ID1)%Ay + FORCE(ID1)%Yint
            FORCE(ID1)%Zext =                                                  &
     &        NODE(ID1)%Mass * MOTION(ID1)%Az + FORCE(ID1)%Zint
!!
            Fmsx = Wg2 * (FORCE(IMS)%Xext - FORCE(IMS)%Xint)
            Fmsy = Wg2 * (FORCE(IMS)%Yext - FORCE(IMS)%Yint)
            Fmsz = Wg2 * (FORCE(IMS)%Zext - FORCE(IMS)%Zint)
!!
!! Recompute acceleration of the second end point.
!!
            Qinv = ONE / CONSTRAINED_NODE(i)%TMass(2)
            MOTION(ID2)%Ax =                                                   &
     &        (FORCE(ID2)%Xext - FORCE(ID2)%Xint + Fmsx) * Qinv
            MOTION(ID2)%Ay =                                                   &
     &        (FORCE(ID2)%Yext - FORCE(ID2)%Yint + Fmsy) * Qinv
            MOTION(ID2)%Az =                                                   &
     &        (FORCE(ID2)%Zext - FORCE(ID2)%Zint + Fmsz) * Qinv
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
            FORCE(ID2)%Xext =                                                  &
     &        NODE(ID2)%Mass * MOTION(ID2)%Ax + FORCE(ID2)%Xint
            FORCE(ID2)%Yext =                                                  &
     &        NODE(ID2)%Mass * MOTION(ID2)%Ay + FORCE(ID2)%Yint
            FORCE(ID2)%Zext =                                                  &
     &        NODE(ID2)%Mass * MOTION(ID2)%Az + FORCE(ID2)%Zint
!!
!! Check to see if this is a nodal point with rotational degress of freedom.
!!
            IF (NODE(IMS)%IRT .GT. 0) THEN
              IMS = NODE(IMS)%IRT
              ID1 = NODE(ID1)%IRT
              ID2 = NODE(ID2)%IRT
!!
!! Recompute acceleration of the first end point.
!!
              Fmsx = Wg1 * (FORCE(IMS)%Xext - FORCE(IMS)%Xint)
              Fmsy = Wg1 * (FORCE(IMS)%Yext - FORCE(IMS)%Yint)
              Fmsz = Wg1 * (FORCE(IMS)%Zext - FORCE(IMS)%Zint)
!!
              Qinv = ONE / CONSTRAINED_NODE(i)%RMass(1)
              MOTION(ID1)%Ax =                                                 &
     &          (FORCE(ID1)%Xext - FORCE(ID1)%Xint + Fmsx) * Qinv
              MOTION(ID1)%Ay =                                                 &
     &          (FORCE(ID1)%Yext - FORCE(ID1)%Yint + Fmsy) * Qinv
              MOTION(ID1)%Az =                                                 &
     &          (FORCE(ID1)%Zext - FORCE(ID1)%Zint + Fmsz) * Qinv
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
              FORCE(ID1)%Xext =                                                &
     &          NODE(ID1)%Mass * MOTION(ID1)%Ax + FORCE(ID1)%Xint
              FORCE(ID1)%Yext =                                                &
     &          NODE(ID1)%Mass * MOTION(ID1)%Ay + FORCE(ID1)%Yint
              FORCE(ID1)%Zext =                                                &
     &          NODE(ID1)%Mass * MOTION(ID1)%Az + FORCE(ID1)%Zint
!!
!! Recompute acceleration of the second end point.
!!
              Fmsx = Wg2 * (FORCE(IMS)%Xext - FORCE(IMS)%Xint)
              Fmsy = Wg2 * (FORCE(IMS)%Yext - FORCE(IMS)%Yint)
              Fmsz = Wg2 * (FORCE(IMS)%Zext - FORCE(IMS)%Zint)
!!
              Qinv = ONE / CONSTRAINED_NODE(i)%RMass(2)
              MOTION(ID2)%Ax =                                                 &
     &          (FORCE(ID2)%Xext - FORCE(ID2)%Xint + Fmsx) * Qinv
              MOTION(ID2)%Ay =                                                 &
     &          (FORCE(ID2)%Yext - FORCE(ID2)%Yint + Fmsy) * Qinv
              MOTION(ID2)%Az =                                                 &
     &          (FORCE(ID2)%Zext - FORCE(ID2)%Zint + Fmsz) * Qinv
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
              FORCE(ID2)%Xext =                                                &
     &          NODE(ID2)%Mass * MOTION(ID2)%Ax + FORCE(ID2)%Xint
              FORCE(ID2)%Yext =                                                &
     &          NODE(ID2)%Mass * MOTION(ID2)%Ay + FORCE(ID2)%Yint
              FORCE(ID2)%Zext =                                                &
     &          NODE(ID2)%Mass * MOTION(ID2)%Az + FORCE(ID2)%Zint
!!
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! ### SECOND, MODIFY ALL MIDSIDE NODAL POINTS
!! Loop over all constraints.
!!
      DO i = 1,NUMNC
!!
        IF (CONSTRAINED_NODE(i)%CNID .GT. 0) THEN
          IMS = CONSTRAINED_NODE(i)%CNID
          ID1 = CONSTRAINED_NODE(i)%NPID(1)
          ID2 = CONSTRAINED_NODE(i)%NPID(2)
          Wg1 = CONSTRAINED_NODE(i)%Weight
          Wg2 = ONE - Wg1
!!
!! Apply kinematic boundary condition only to those nodal points that are
!! being integrated.
!!
          IF (NODE(IMS)%On) THEN
!!
!! Interpolate to get acceleration of midside node.
!!
            MOTION(IMS)%Ax = Wg1 * MOTION(ID1)%Ax + Wg2 * MOTION(ID2)%Ax
            MOTION(IMS)%Ay = Wg1 * MOTION(ID1)%Ay + Wg2 * MOTION(ID2)%Ay
            MOTION(IMS)%Az = Wg1 * MOTION(ID1)%Az + Wg2 * MOTION(ID2)%Az
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
            FORCE(IMS)%Xext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Ax + FORCE(IMS)%Xint
            FORCE(IMS)%Yext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Ay + FORCE(IMS)%Yint
            FORCE(IMS)%Zext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Az + FORCE(IMS)%Zint
!!
!! Check to see if this is a nodal point with rotational degress of freedom.
!!
            IF (NODE(IMS)%IRT .GT. 0) THEN
              IMS = NODE(IMS)%IRT
              ID1 = NODE(ID1)%IRT
              ID2 = NODE(ID2)%IRT
!!
!! Interpolate to get acceleration of midside node.
!!
              MOTION(IMS)%Ax =                                                 &
     &          Wg1 * MOTION(ID1)%Ax + Wg2 * MOTION(ID2)%Ax
              MOTION(IMS)%Ay =                                                 &
     &          Wg1 * MOTION(ID1)%Ay + Wg2 * MOTION(ID2)%Ay
              MOTION(IMS)%Az =                                                 &
     &          Wg1 * MOTION(ID1)%Az + Wg2 * MOTION(ID2)%Az
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
              FORCE(IMS)%Xext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Ax + FORCE(IMS)%Xint
              FORCE(IMS)%Yext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Ay + FORCE(IMS)%Yint
              FORCE(IMS)%Zext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Az + FORCE(IMS)%Zint
!!
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_MIDSIDE_CONSTRAINTS2
!!
!! Copyright (c) by KEY Associates; 19-OCT-1995 20:58:34.00
!!
!! Purpose: For those nodal points located on the midsides of adjacent elements
!! compute the constrained motion. This routine is a "last pass" that insures
!! the midside node acceleration is exactly half of the controlling end nodes.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE force_
      USE constrained_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Loop over all constraints.
!!
      DO i = 1,NUMNC
!!
        IF (CONSTRAINED_NODE(i)%CNID .GT. 0) THEN
          IMS = CONSTRAINED_NODE(i)%CNID
          ID1 = CONSTRAINED_NODE(i)%NPID(1)
          ID2 = CONSTRAINED_NODE(i)%NPID(2)
          Wg1 = CONSTRAINED_NODE(i)%Weight
          Wg2 = ONE - Wg1
!!
!! Apply kinematic boundary condition only to those nodal points that are
!! being integrated.
!!
          IF (NODE(IMS)%On) THEN
!!
!! Redefine acceleration of the midside nodal point.
!!
            MOTION(IMS)%Ax = Wg1*MOTION(ID1)%Ax + Wg2*MOTION(ID2)%Ax
            MOTION(IMS)%Ay = Wg1*MOTION(ID1)%Ay + Wg2*MOTION(ID2)%Ay
            MOTION(IMS)%Az = Wg1*MOTION(ID1)%Az + Wg2*MOTION(ID2)%Az
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
            FORCE(IMS)%Xext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Ax + FORCE(IMS)%Xint
            FORCE(IMS)%Yext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Ay + FORCE(IMS)%Yint
            FORCE(IMS)%Zext =                                                  &
     &        NODE(IMS)%Mass * MOTION(IMS)%Az + FORCE(IMS)%Zint
!!
!! Check to see if this is a nodal point with rotational degress of freedom.
!!
            IF (NODE(IMS)%IRT .GT. 0) THEN
              IMS = NODE(IMS)%IRT
              ID1 = NODE(ID1)%IRT
              ID2 = NODE(ID2)%IRT
!!
!! Redefine acceleration of the midside nodal point.
!!
              MOTION(IMS)%Ax =                                                 &
     &          Wg1 * MOTION(ID1)%Ax + Wg2 * MOTION(ID2)%Ax
              MOTION(IMS)%Ay =                                                 &
     &          Wg1 * MOTION(ID1)%Ay + Wg2 * MOTION(ID2)%Ay
              MOTION(IMS)%Az =                                                 &
     &          Wg1 * MOTION(ID1)%Az + Wg2 * MOTION(ID2)%Az
!!
!! Recompute the external forces to reflect the presence of the constraints.
!!
              FORCE(IMS)%Xext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Ax + FORCE(IMS)%Xint
              FORCE(IMS)%Yext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Ay + FORCE(IMS)%Yint
              FORCE(IMS)%Zext =                                                &
     &          NODE(IMS)%Mass * MOTION(IMS)%Az + FORCE(IMS)%Zint
!!
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_SPOT_WELD
!!
!! Copyright (c) by KEY Associates; 16-JUN-1994 20:47:52.00
!!
!! Purpose: Impose acceleration constraints that duplicate the attachments
!! between parts (triangular and quadrilateral shell elements) that a spot
!! weld provides.
!!
      USE shared_common_data
      USE spot_weld_
      USE platt_
      USE platq_
      USE motion_
      USE force_
      USE node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                         &
     &          NP(3),                & ! -/- Spot welded nodal points         &
     &          EL(3),                & ! -/- Spot welded elements             &
     &          IX(12)                  ! -/- Constrained nodal points
      REAL(KIND(0D0))                 &
     &          W(12),                & ! -/- Acceleration aligned w/ Tx,Ty,Tz &
     &          Ax(12),dAx(12),       & ! -/- X-componet of acceleration       &
     &          Ay(12),dAy(12),       & ! -/- Y-componet of acceleration       &
     &          Az(12),dAz(12),       & ! -/- Z-componet of acceleration       &
     &          C(12,4),              & ! -/- Constraint matrix                &
     &          Scale(3,4),           & ! -/- Constraint magnitudes**2         &
     &          Ith_Constraint_Force,                                          &
     &          Max_Constraint_Force
      LOGICAL                                                                  &
     &          FAILED

      LOGICAL, SAVE :: FIRST
!!
!! Loop over all spot weld specifications. Note that up to three adjacent
!! nodal points or shell elements can be spot welded together.
!! SPOT_WELD,swid,NODES,npid1,npid2,npid3,fail
!! SPOT_WELD,swid,POSITION,px,py,pz,fail
!!
      DO NSW = 1,NUMSW
        IF (.NOT.SPOT_WELD(NSW)%FAILED) THEN
          Fmax = SPOT_WELD(NSW)%Fmax

          IF (SPOT_WELD(NSW)%PLACE .EQ. 'NODES') THEN
            NP(1) = SPOT_WELD(NSW)%NPID(1)
            NP(2) = SPOT_WELD(NSW)%NPID(2)
            NP(3) = SPOT_WELD(NSW)%NPID(3)
!!
!! Find the total force F acting on the spot welded nodal points.
!!
            Fx = 0.0
            Fy = 0.0
            Fz = 0.0
            Qm = 0.0
            Gx = 0.0
            Gy = 0.0
            Gz = 0.0
            Rm = 0.0
            DO i = 1,3
              IF (NP(i) .GT. 0) THEN
                NTR = NP(i)
                Fx = Fx + FORCE(NTR)%Xext - FORCE(NTR)%Xint
                Fy = Fy + FORCE(NTR)%Yext - FORCE(NTR)%Yint
                Fz = Fz + FORCE(NTR)%Zext - FORCE(NTR)%Zint
                Qm = Qm + NODE(NTR)%Mass
!!
!! Check for rotational degrees of freedom.
!!
                IF (NODE(NP(i))%IRT .GT. 0) THEN
                  IRT = NODE(NP(i))%IRT
                  Gx = Gx + FORCE(IRT)%Xext - FORCE(IRT)%Xint
                  Gy = Gy + FORCE(IRT)%Yext - FORCE(IRT)%Yint
                  Gz = Gz + FORCE(IRT)%Zext - FORCE(IRT)%Zint
                  Rm = Rm + NODE(IRT)%Mass
                ENDIF
              ENDIF
            ENDDO
!!
!! Find the translational acceleration of the center of mass of the
!! spot welded nodal points.
!!
            Axspw = Fx / Qm
            Ayspw = Fy / Qm
            Azspw = Fz / Qm
!!
!! Find the rotational acceleration of the center of mass of the spot
!! welded nodal points.
!!
            IF (Rm .GT. 0.0) THEN
              Bxspw = Gx / Rm
              Byspw = Gy / Rm
              Bzspw = Gz / Rm
            ELSE
              Bxspw = 0.0
              Byspw = 0.0
              Bzspw = 0.0
            ENDIF
!!
!! Find the maximum "external" force generated by the spot weld constraint.
!! (Ignores constraint torques related to rotational degrees of freedom.)
!! (This failure criterion does not distinguish between shear and normal.)
!! (It also does not distinguish between tension or compression.)
!!
            Max_Constraint_Force = 0.0
            DO i = 1,3
              IF (NP(i) .GT. 0) THEN
                NTR = NP(i)
                Fx = NODE(NTR)%Mass * Axspw -                                  &
     &            (FORCE(NTR)%Xext - FORCE(NTR)%Xint)
                Fy = NODE(NTR)%Mass * Ayspw -                                  &
     &            (FORCE(NTR)%Yext - FORCE(NTR)%Yint)
                Fz = NODE(NTR)%Mass * Azspw -                                  &
     &            (FORCE(NTR)%Zext - FORCE(NTR)%Zint)
                Ith_Constraint_Force = SQRT (Fx*Fx +Fy*Fy +Fz*Fz)
                Max_Constraint_Force =                                         &
     &            MAX (Max_Constraint_Force, Ith_Constraint_Force)
              ENDIF
            ENDDO
            SPOT_WELD(NSW)%Force = Max_Constraint_Force
!!
!! Check for a tensile failure of the spot weld.
!!
            FAILED = (Fmax.GT.0.0) .AND. (Max_Constraint_Force.GE.Fmax)
            SPOT_WELD(NSW)%FAILED = FAILED

            IF (.NOT.SPOT_WELD(NSW)%FAILED) THEN
!!
!! Distinguish between coincident nodal points and non-coincident nodal points.
!! Non-coincident nodal points must be treated as a rigid body "dumb bell."
!!
              IF (SPOT_WELD(NSW)%Coincident) THEN
!!
!! Coincident nodal points. Treat the "spot welded" nodal points as one.
!! Make the motion of the spot welded nodal points identical.
!!
                DO i = 1,3
                  IF (NP(i) .GT. 0) THEN
                    NTR = NP(i)
                    MOTION(NTR)%Ax = Axspw
                    MOTION(NTR)%Ay = Ayspw
                    MOTION(NTR)%Az = Azspw
                    FORCE(NTR)%Xext =                                          &
     &                NODE(NTR)%Mass * Axspw + FORCE(NTR)%Xint
                    FORCE(NTR)%Yext =                                          &
     &                NODE(NTR)%Mass * Ayspw + FORCE(NTR)%Yint
                    FORCE(NTR)%Zext =                                          &
     &                NODE(NTR)%Mass * Azspw + FORCE(NTR)%Zint
                    IF (NODE(NP(i))%IRT .GT. 0) THEN
                      IRT = NODE(NP(i))%IRT
                      MOTION(IRT)%Ax = Bxspw
                      MOTION(IRT)%Ay = Byspw
                      MOTION(IRT)%Az = Bzspw
                      FORCE(IRT)%Xext =                                        &
     &                  NODE(IRT)%Mass * Bxspw + FORCE(IRT)%Xint
                      FORCE(IRT)%Yext =                                        &
     &                  NODE(IRT)%Mass * Byspw + FORCE(IRT)%Yint
                      FORCE(IRT)%Zext =                                        &
     &                  NODE(IRT)%Mass * Bzspw + FORCE(IRT)%Zint
                    ENDIF
                  ENDIF
                ENDDO

              ELSE IF (.NOT.SPOT_WELD(NSW)%Coincident) THEN
!!
!! Non-coincident nodal points. Treat the "spot welded" nodal points as
!! a rigid body dumb bell. (During initialization, the intermediate nodal
!! point NP(2) was placed on a line between between the two end-points.
!!
!! (Temporary fatal error message. Needs to be replaced with coding.)
!!
                WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
                CALL USER_MESSAGE                                              &
     &                  (                                                      &
     &                  MSGL//'FATAL'//                                        &
     &                  MSGL//'IMPOSE_SPOT_WELD.001.01'//                      &
     &                  MSGL//'SPOTWELD Input Record ID:'//MSG1//              &
     &                  MSGL//'References Non-Coincident Nodal Points.'        &
     &                  )
              ENDIF
            ENDIF

          ELSE IF (SPOT_WELD(NSW)%PLACE .EQ. 'POSITION') THEN
            EL(1) = SPOT_WELD(NSW)%EleID(1)
            EL(2) = SPOT_WELD(NSW)%EleID(2)
            EL(3) = SPOT_WELD(NSW)%EleID(3)
!!
!! Construct constraint matrix C. The constraint matrix has the
!! following structure --
!!  Cols          1             2             3             4
!!  Rows
!! 1  5  9   N1(Xi1,Xi2) N1,x(Xi1,Xi2) N1,y(Xi1,Xi2) N1,z(Xi1,Xi2)
!! 2  6 10   N2(Xi1,Xi2) N2,x(Xi1,Xi2) N2,y(Xi1,Xi2) N2,z(Xi1,Xi2)
!! 3  7 11   N3(Xi1,Xi2) N3,x(Xi1,Xi2) N3,y(Xi1,Xi2) N3,z(Xi1,Xi2)
!! 4  8 12   N4(Xi1,Xi2) N4,x(Xi1,Xi2) N4,y(Xi1,Xi2) N4,z(Xi1,Xi2)
!!
!! El(1): Rows 1 -  4
!! El(2): Rows 5 -  8
!! El(3): Rows 9 - 12
!!
!! Column 1 contains the respective element shape functions
!! evaluated at the spot weld position and are used to impose equal
!! translational accelerations. Example of use with Ax --
!!
!!          Ax = A1x*N1 + A2x*N2 + A3x*N3 + A4x*N4   (*)
!!
!! where Ax = SUM (MiAxi) / SUM (Mi) , i={all elements}.
!! That is, Ax is the mass averaged acceleration. An inner product
!! with Ni's is used to obtain the unconstrained acceleration in
!! each element. The nodal point accelerations are modified with
!! an acceleration increment, the profile of which is given by the
!! Ni's, so that constraint equation (*) is satisfied.
!!
!! Columns 2, 3 and 4 contain the x, y, and z derivatives of the
!! shape functions evaluated at the spot weld position and are used
!! to impose equal tilts and in-plane spin.
!!
            DO n = 1,12
              IX(n) = 0
            ENDDO
            DO n = 1,48
              C(n,1) = 0.0
            ENDDO
            FIRST = .TRUE.

            DO i = 1,3
              IF (EL(i) .GT. 0) THEN
!!
!! Evaluate shape functions.
!!
                Xi1 = SPOT_WELD(NSW)%Xi1(i)
                Xi2 = SPOT_WELD(NSW)%Xi2(i)
                IF (SPOT_WELD(NSW)%Type(i) .EQ. 0) THEN
                  A1 = ONE - Xi1 - Xi2
                  A2 = Xi1
                  A3 = Xi2
                  A4 = 0.0
                  A11 = -ONE
                  A21 = +ONE
                  A31 =  0.0
                  A41 =  0.0
                  A12 = -ONE
                  A22 =  0.0
                  A32 = +ONE
                  A42 =  0.0
                ELSE IF (SPOT_WELD(NSW)%Type(i) .EQ. 1) THEN
                  A1 = 0.25*(ONE - Xi1)*(ONE - Xi2)
                  A2 = 0.25*(ONE + Xi1)*(ONE - Xi2)
                  A3 = 0.25*(ONE + Xi1)*(ONE + Xi2)
                  A4 = 0.25*(ONE - Xi1)*(ONE + Xi2)
                  A11 = -0.25*(ONE - Xi2)
                  A21 = +0.25*(ONE - Xi2)
                  A31 = +0.25*(ONE + Xi2)
                  A41 = -0.25*(ONE + Xi2)
                  A12 = -0.25*(ONE - Xi1)
                  A22 = -0.25*(ONE + Xi1)
                  A32 = +0.25*(ONE + Xi1)
                  A42 = +0.25*(ONE - Xi1)
                ENDIF
!!
!! Retrieve element coordinates and construct [dXi_i/dX_j] matrix.
!!
                IF (SPOT_WELD(NSW)%Type(i) .EQ. 0) THEN
                  I1 = PLATT(EL(i))%PAR%IX(1)
                  I2 = PLATT(EL(i))%PAR%IX(2)
                  I3 = PLATT(EL(i))%PAR%IX(3)
                  I4 = 0

                  X1 = MOTION(I1)%Px + MOTION(I1)%Ux
                  X2 = MOTION(I2)%Px + MOTION(I2)%Ux
                  X3 = MOTION(I3)%Px + MOTION(I3)%Ux
                  X4 = 0.0

                  Y1 = MOTION(I1)%Py + MOTION(I1)%Uy
                  Y2 = MOTION(I2)%Py + MOTION(I2)%Uy
                  Y3 = MOTION(I3)%Py + MOTION(I3)%Uy
                  Y4 = 0.0

                  Z1 = MOTION(I1)%Pz + MOTION(I1)%Uz
                  Z2 = MOTION(I2)%Pz + MOTION(I2)%Uz
                  Z3 = MOTION(I3)%Pz + MOTION(I3)%Uz
                  Z4 = 0.0

                ELSE IF (SPOT_WELD(NSW)%Type(i) .EQ. 1) THEN
                  I1 = PLATQ(EL(i))%PAR%IX(1)
                  I2 = PLATQ(EL(i))%PAR%IX(2)
                  I3 = PLATQ(EL(i))%PAR%IX(3)
                  I4 = PLATQ(EL(i))%PAR%IX(4)

                  X1 = MOTION(I1)%Px + MOTION(I1)%Ux
                  X2 = MOTION(I2)%Px + MOTION(I2)%Ux
                  X3 = MOTION(I3)%Px + MOTION(I3)%Ux
                  X4 = MOTION(I4)%Px + MOTION(I4)%Ux

                  Y1 = MOTION(I1)%Py + MOTION(I1)%Uy
                  Y2 = MOTION(I2)%Py + MOTION(I2)%Uy
                  Y3 = MOTION(I3)%Py + MOTION(I3)%Uy
                  Y4 = MOTION(I4)%Py + MOTION(I4)%Uy

                  Z1 = MOTION(I1)%Pz + MOTION(I1)%Uz
                  Z2 = MOTION(I2)%Pz + MOTION(I2)%Uz
                  Z3 = MOTION(I3)%Pz + MOTION(I3)%Uz
                  Z4 = MOTION(I4)%Pz + MOTION(I4)%Uz
                ENDIF
!!
!! Construct derivative matrix [Dij] = [dX_i/dXi_j].
!!
                D11 = X1*A11 + X2*A21 + X3*A31 + X4*A41
                D21 = Y1*A11 + Y2*A21 + Y3*A31 + Y4*A41
                D31 = Z1*A11 + Z2*A21 + Z3*A31 + Z4*A41

                D12 = X1*A12 + X2*A22 + X3*A32 + X4*A42
                D22 = Y1*A12 + Y2*A22 + Y3*A32 + Y4*A42
                D32 = Z1*A12 + Z2*A22 + Z3*A32 + Z4*A42
!!
!! Define the unit vector T normal to the element by taking the
!! cross product of two in-plane tangent vectors at the spot
!! weld location. (Using T = Tx,Ty,Tz for D13,D23,D33 is tanta-
!! mont to saying that Xi3 is a true-length coordinate perpen-
!! dicular to the shell reference surface.)
!!
                IF (FIRST) THEN
                  Rx = D11
                  Ry = D21
                  Rz = D31
                  Sx = D12
                  Sy = D22
                  Sz = D32
                  Tx = Ry*Sz - Sy*Rz
                  Ty = Rz*Sx - Sz*Rx
                  Tz = Rx*Sy - Sx*Ry
                  Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)

                  IF (Tmag .EQ. 0.0) THEN
                    WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
                    CALL USER_MESSAGE                                          &
     &              (                                                          &
     &              MSGL//'FATAL'//                                            &
     &              MSGL//'IMPOSE_SPOT_WELD.001.02'//                          &
     &              MSGL//'SPOTWELD Input Record ID:'//MSG1//                  &
     &              MSGL//'Element Geometry Has An Undefined Normal.'          &
     &              )
                  ENDIF

                  Tx = Tx * (ONE / Tmag)
                  Ty = Ty * (ONE / Tmag)
                  Tz = Tz * (ONE / Tmag)
                  D13 = Tx
                  D23 = Ty
                  D33 = Tz
                  FIRST = .FALSE.
                ENDIF
!!
!! Invert [Dij] to get [Bij] = [dXi_i/dX_j].
!!
                DET = D11*D22*D33+D12*D23*D31+D13*D21*D32                      &
     &              - D31*D22*D13-D32*D23*D11-D33*D21*D12
                B11 = (D22*D33-D32*D23)*(ONE/DET)
                B12 = (D32*D13-D12*D33)*(ONE/DET)
                B13 = (D12*D23-D13*D22)*(ONE/DET)
                B21 = (D31*D23-D21*D33)*(ONE/DET)
                B22 = (D11*D33-D31*D13)*(ONE/DET)
                B23 = (D21*D13-D11*D23)*(ONE/DET)
!!                 B31 = (D21*D32-D31*D22)*(ONE/DET)
!!                 B32 = (D31*D12-D11*D32)*(ONE/DET)
!!                 B33 = (D22*D11-D12*D21)*(ONE/DET)
!!
!! Construct A1,x = A1,1*B11 + A1,2*B21 + A1,3*B31 et cetera.
!! (Please note that An,3 is always zero.)
!!
                A1x = A11*B11 + A12*B21
                A1y = A11*B12 + A12*B22
                A1z = A11*B13 + A12*B23

                A2x = A21*B11 + A22*B21
                A2y = A21*B12 + A22*B22
                A2z = A21*B13 + A22*B23

                A3x = A31*B11 + A32*B21
                A3y = A31*B12 + A32*B22
                A3z = A31*B13 + A32*B23

                A4x = A41*B11 + A42*B21
                A4y = A41*B12 + A42*B22
                A4z = A41*B13 + A42*B23
!!
!! Remove derivative normal to shell reference surfaces.  The normal
!! direction is defined by the unit vector T with components (Tx,Ty,Tz).
!!
                A1t = A1x*Tx + A1y*Ty + A1z*Tz
                A1x = A1x - A1t*Tx
                A1y = A1y - A1t*Ty
                A1z = A1z - A1t*Tz

                A2t = A2x*Tx + A2y*Ty + A2z*Tz
                A2x = A2x - A2t*Tx
                A2y = A2y - A2t*Ty
                A2z = A2z - A2t*Tz

                A3t = A3x*Tx + A3y*Ty + A3z*Tz
                A3x = A3x - A3t*Tx
                A3y = A3y - A3t*Ty
                A3z = A3z - A3t*Tz

                A4t = A4x*Tx + A4y*Ty + A4z*Tz
                A4x = A4x - A4t*Tx
                A4y = A4y - A4t*Ty
                A4z = A4z - A4t*Tz
!!
!! Store shape functions in constraint matrix C for later processing.
!!
                IF (i .EQ. 1) THEN
                  C(1,1) = A1
                  C(2,1) = A2
                  C(3,1) = A3
                  C(4,1) = A4
                  C(1,2) = A1x
                  C(2,2) = A2x
                  C(3,2) = A3x
                  C(4,2) = A4x
                  C(1,3) = A1y
                  C(2,3) = A2y
                  C(3,3) = A3y
                  C(4,3) = A4y
                  C(1,4) = A1z
                  C(2,4) = A2z
                  C(3,4) = A3z
                  C(4,4) = A4z
                  IX(1)  = I1
                  IX(2)  = I2
                  IX(3)  = I3
                  IX(4)  = I4
                ENDIF
                IF (i .EQ. 2) THEN
                  C(5,1) = A1
                  C(6,1) = A2
                  C(7,1) = A3
                  C(8,1) = A4
                  C(5,2) = A1x
                  C(6,2) = A2x
                  C(7,2) = A3x
                  C(8,2) = A4x
                  C(5,3) = A1y
                  C(6,3) = A2y
                  C(7,3) = A3y
                  C(8,3) = A4y
                  C(5,4) = A1z
                  C(6,4) = A2z
                  C(7,4) = A3z
                  C(8,4) = A4z
                  IX(5)  = I1
                  IX(6)  = I2
                  IX(7)  = I3
                  IX(8)  = I4
                ENDIF
                IF (i .EQ. 3) THEN
                  C( 9,1) = A1
                  C(10,1) = A2
                  C(11,1) = A3
                  C(12,1) = A4
                  C( 9,2) = A1x
                  C(10,2) = A2x
                  C(11,2) = A3x
                  C(12,2) = A4x
                  C( 9,3) = A1y
                  C(10,3) = A2y
                  C(11,3) = A3y
                  C(12,3) = A4y
                  C( 9,4) = A1z
                  C(10,4) = A2z
                  C(11,4) = A3z
                  C(12,4) = A4z
                  IX( 9)  = I1
                  IX(10)  = I2
                  IX(11)  = I3
                  IX(12)  = I4
                ENDIF
              ENDIF
            ENDDO
!!
!! Compute scale factors for constraint vectors (The magnitude
!! of the constraint vector squared is used to get the correct
!! magnitude for the incremental acceleration correction.)
!!
            DO n = 1,12
              Scale(n,1) = 0.0
            ENDDO
            DO n = 1,4
              Scale(1,1) = Scale(1,1) + C(n,1)*C(n,1)
              Scale(1,2) = Scale(1,2) + C(n,2)*C(n,2)
              Scale(1,3) = Scale(1,3) + C(n,3)*C(n,3)
              Scale(1,4) = Scale(1,4) + C(n,4)*C(n,4)
            ENDDO
            DO n = 5,8
              Scale(2,1) = Scale(2,1) + C(n,1)*C(n,1)
              Scale(2,2) = Scale(2,2) + C(n,2)*C(n,2)
              Scale(2,3) = Scale(2,3) + C(n,3)*C(n,3)
              Scale(2,4) = Scale(2,4) + C(n,4)*C(n,4)
            ENDDO
            DO n = 9,12
              Scale(3,1) = Scale(3,1) + C(n,1)*C(n,1)
              Scale(3,2) = Scale(3,2) + C(n,2)*C(n,2)
              Scale(3,3) = Scale(3,3) + C(n,3)*C(n,3)
              Scale(3,4) = Scale(3,4) + C(n,4)*C(n,4)
            ENDDO
            DO n = 1,12
              IF (Scale(n,1) .NE. 0.0) Scale(n,1) = ONE / Scale(n,1)
            ENDDO
!!
!! Retrieve nodal point accelerations, and compute mass averaging
!! coefficients for constraint values for the spot weld motion.
!!
            DO n = 1,12
              IF (IX(n) .EQ. 0) THEN
                Ax(n) = 0.0
                Ay(n) = 0.0
                Az(n) = 0.0
              ELSE
                Ax(n) = MOTION(IX(n))%Ax
                Ay(n) = MOTION(IX(n))%Ay
                Az(n) = MOTION(IX(n))%Az
              ENDIF
            ENDDO

            Amwt1 = 0.0
            Amwt2 = 0.0
            Amwt3 = 0.0
            DO n = 1,4
              IF (IX(n  ) .NE. 0) THEN
                Amwt1 = Amwt1 + NODE(IX(n  ))%Mass
              ENDIF
              IF (IX(n+4) .NE. 0) THEN
                Amwt2 = Amwt2 + NODE(IX(n+4))%Mass
              ENDIF
              IF (IX(n+8) .NE. 0) THEN
                Amwt3 = Amwt3 + NODE(IX(n+8))%Mass
              ENDIF
            ENDDO
            Amass = Amwt1 + Amwt2 + Amwt3
            Amwt1 = Amwt1 / Amass
            Amwt2 = Amwt2 / Amass
            Amwt3 = Amwt3 / Amass
!!
!! Impose constraints on translation. Constraints are imposed by finding
!! the component of the acceleration related to the spot weld points in
!! the respective elements taking different paths and removing it from the
!! nodal point accelerations. Restraint forces are added to the exterior
!! forces. The restraint forces represent the constraint of the spot_weld.
!! First, find unrestrained accelerations at spot weld.
!!
            C1dotAx = 0.0
            C1dotAy = 0.0
            C1dotAz = 0.0
            C2dotAx = 0.0
            C2dotAy = 0.0
            C2dotAz = 0.0
            C3dotAx = 0.0
            C3dotAy = 0.0
            C3dotAz = 0.0
            DO n = 1,4
              C1dotAx = C1dotAx + C(n  ,1)*Ax(n  )
              C1dotAy = C1dotAy + C(n  ,1)*Ay(n  )
              C1dotAz = C1dotAz + C(n  ,1)*Az(n  )
              C2dotAx = C2dotAx + C(n+4,1)*Ax(n+4)
              C2dotAy = C2dotAy + C(n+4,1)*Ay(n+4)
              C2dotAz = C2dotAz + C(n+4,1)*Az(n+4)
              C3dotAx = C3dotAx + C(n+8,1)*Ax(n+8)
              C3dotAy = C3dotAy + C(n+8,1)*Ay(n+8)
              C3dotAz = C3dotAz + C(n+8,1)*Az(n+8)
            ENDDO
!!
!! Compute weighted average of accelerations to obtain spot weld value.
!!
            Axspw = Amwt1 * C1dotAx + Amwt2 * C2dotAx + Amwt3 * C3dotAx
            Ayspw = Amwt1 * C1dotAy + Amwt2 * C2dotAy + Amwt3 * C3dotAy
            Azspw = Amwt1 * C1dotAz + Amwt2 * C2dotAz + Amwt3 * C3dotAz
!!
!! Compute incremental accelerations due to translational constraint.
!!
            DO n = 1,4
              dAx(n  ) = (Scale(1,1) * (Axspw - C1dotAx)) * C(n  ,1)
              dAy(n  ) = (Scale(1,1) * (Ayspw - C1dotAy)) * C(n  ,1)
              dAz(n  ) = (Scale(1,1) * (Azspw - C1dotAz)) * C(n  ,1)
              dAx(n+4) = (Scale(2,1) * (Axspw - C2dotAx)) * C(n+4,1)
              dAy(n+4) = (Scale(2,1) * (Ayspw - C2dotAy)) * C(n+4,1)
              dAz(n+4) = (Scale(2,1) * (Azspw - C2dotAz)) * C(n+4,1)
              dAx(n+8) = (Scale(3,1) * (Axspw - C3dotAx)) * C(n+8,1)
              dAy(n+8) = (Scale(3,1) * (Ayspw - C3dotAy)) * C(n+8,1)
              dAz(n+8) = (Scale(3,1) * (Azspw - C3dotAz)) * C(n+8,1)
            ENDDO
!!
!! Compute normal acceleration (acceleration aligned with T = (Tx,Ty,Tz)).
!!
            DO n = 1,12
              W(n) = Tx*Ax(n) + Ty*Ay(n) + Tz*Az(n)
            ENDDO
!!
!! Impose constraints on tilt. Constraints are imposed by finding the
!! component of the acceleration related to the spot weld points in the
!! respective elements taking different paths and removing it from the
!! nodal point accelerations. Restraint forces are added to the exterior
!! forces. The restraint forces represent the constraint of the spot_weld.
!! First, find unrestrained tilt accelerations at spot weld. (The gradient
!! operator has already had the ,t removed, that is, ,x ,y and ,z only con-
!! tain the in-plane derivatives ,r and ,s. Thus, W,x W,y and W,z describe
!! the tilt (W,r and W,s) at the spot weld location. The spot weld tilt
!! constraints mean that the elements remain "tangent" at the spot weld.)
!!
            C1dotWx = 0.0
            C1dotWy = 0.0
            C1dotWz = 0.0
            C2dotWx = 0.0
            C2dotWy = 0.0
            C2dotWz = 0.0
            C3dotWx = 0.0
            C3dotWy = 0.0
            C3dotWz = 0.0
            DO n = 1,4
              C1dotWx = C1dotWx + C(n  ,2)*W(n  )
              C1dotWy = C1dotWy + C(n  ,3)*W(n  )
              C1dotWz = C1dotWz + C(n  ,4)*W(n  )
              C2dotWx = C2dotWx + C(n+4,2)*W(n+4)
              C2dotWy = C2dotWy + C(n+4,3)*W(n+4)
              C2dotWz = C2dotWz + C(n+4,4)*W(n+4)
              C3dotWx = C3dotWx + C(n+8,2)*W(n+8)
              C3dotWy = C3dotWy + C(n+8,3)*W(n+8)
              C3dotWz = C3dotWz + C(n+8,4)*W(n+8)
            ENDDO
!!
!! Compute weighted average of accelerations to obtain spot weld value.
!!
            Wxspw = Amwt1 * C1dotWx + Amwt2 * C2dotWx + Amwt3 * C3dotWx
            Wyspw = Amwt1 * C1dotWy + Amwt2 * C2dotWy + Amwt3 * C3dotWy
            Wzspw = Amwt1 * C1dotWz + Amwt2 * C2dotWz + Amwt3 * C3dotWz
!!
!! Compute incremental accelerations due to tilt constraint.
!!
            DO n = 1,4
              dW1 = (Scale(1,2) * (Wxspw - C1dotWx)) * C(n  ,2)                &
     &            + (Scale(1,3) * (Wyspw - C1dotWy)) * C(n  ,3)                &
     &            + (Scale(1,4) * (Wzspw - C1dotWz)) * C(n  ,4)
              dW2 = (Scale(2,2) * (Wxspw - C2dotWx)) * C(n+4,2)                &
     &            + (Scale(2,3) * (Wyspw - C2dotWy)) * C(n+4,3)                &
     &            + (Scale(2,4) * (Wzspw - C2dotWz)) * C(n+4,4)
              dW3 = (Scale(3,2) * (Wxspw - C3dotWx)) * C(n+8,2)                &
     &            + (Scale(3,3) * (Wyspw - C3dotWy)) * C(n+8,3)                &
     &            + (Scale(3,4) * (Wzspw - C3dotWz)) * C(n+8,4)
              dAx(n  ) = dAx(n  ) + Tx * dW1
              dAy(n  ) = dAy(n  ) + Ty * dW1
              dAz(n  ) = dAz(n  ) + Tz * dW1
              dAx(n+4) = dAx(n+4) + Tx * dW2
              dAy(n+4) = dAy(n+4) + Ty * dW2
              dAz(n+4) = dAz(n+4) + Tz * dW2
              dAx(n+8) = dAx(n+8) + Tx * dW3
              dAy(n+8) = dAy(n+8) + Ty * dW3
              dAz(n+8) = dAz(n+8) + Tz * dW3
            ENDDO
!!
!! Impose constraints on spin. Constraints are imposed by finding the
!! component of the acceleration related to the spot weld points in the
!! respective elements taking different paths and removing it from the
!! nodal point accelerations. Restraint forces are added to the exterior
!! forces. The restraint forces represent the constraint of the spot_weld.
!! (Spin is the skew-symmetric part of the in-plane gradient of the in-
!! plane velocity.) First, remove the normal acceleration W from Ax,Ay,Az.
!!
            DO n = 1,12
              Ax(n) = Ax(n) - Tx * W(n)
              Ay(n) = Ay(n) - Ty * W(n)
              Az(n) = Az(n) - Tz * W(n)
            ENDDO
!!
!! Second compute acceleration gradients at the spot weld.
!!
            C1Axy = 0.0
            C1Axz = 0.0
            C1Ayx = 0.0
            C1Ayz = 0.0
            C1Azx = 0.0
            C1Azy = 0.0
            DO n = 1,4
              C1Axy = C1Axy + C(n,3)*Ax(n)
              C1Axz = C1Axz + C(n,4)*Ax(n)
              C1Ayx = C1Ayx + C(n,2)*Ay(n)
              C1Ayz = C1Ayz + C(n,4)*Ay(n)
              C1Azx = C1Azx + C(n,2)*Az(n)
              C1Azy = C1Azy + C(n,3)*Az(n)
            ENDDO

            C2Axy = 0.0
            C2Axz = 0.0
            C2Ayx = 0.0
            C2Ayz = 0.0
            C2Azx = 0.0
            C2Azy = 0.0
            DO n = 5,8
              C2Axy = C2Axy + C(n,3)*Ax(n)
              C2Axz = C2Axz + C(n,4)*Ax(n)
              C2Ayx = C2Ayx + C(n,2)*Ay(n)
              C2Ayz = C2Ayz + C(n,4)*Ay(n)
              C2Azx = C2Azx + C(n,2)*Az(n)
              C2Azy = C2Azy + C(n,3)*Az(n)
            ENDDO

            C3Axy = 0.0
            C3Axz = 0.0
            C3Ayx = 0.0
            C3Ayz = 0.0
            C3Azx = 0.0
            C3Azy = 0.0
            DO n = 9,12
              C3Axy = C3Axy + C(n,3)*Ax(n)
              C3Axz = C3Axz + C(n,4)*Ax(n)
              C3Ayx = C3Ayx + C(n,2)*Ay(n)
              C3Ayz = C3Ayz + C(n,4)*Ay(n)
              C3Azx = C3Azx + C(n,2)*Az(n)
              C3Azy = C3Azy + C(n,3)*Az(n)
            ENDDO
!!
!! Compute weighted average of spin acceleration to obtain spot weld value.
!! (Since both the normal acceleration has been removed and the normal deri-
!! vative has been removed the only remaining derivatives are in-plane deri-
!! vatives of in-plane acceleration -- in particular Ar,s and As,r)
!!
            C1Wxy = 0.5D+0 * (C1Axy-C1Ayx)
            C1Wyz = 0.5D+0 * (C1Ayz-C1Azy)
            C1Wxz = 0.5D+0 * (C1Axz-C1Azx)
            C2Wxy = 0.5D+0 * (C2Axy-C2Ayx)
            C2Wyz = 0.5D+0 * (C2Ayz-C2Azy)
            C2Wxz = 0.5D+0 * (C2Axz-C2Azx)
            C3Wxy = 0.5D+0 * (C3Axy-C3Ayx)
            C3Wyz = 0.5D+0 * (C3Ayz-C3Azy)
            C3Wxz = 0.5D+0 * (C3Axz-C3Azx)
            Wxyspw = Amwt1 * C1Wxy + Amwt2 * C2Wxy + Amwt3 * C3Wxy
            Wyzspw = Amwt1 * C1Wyz + Amwt2 * C2Wyz + Amwt3 * C3Wyz
            Wxzspw = Amwt1 * C1Wxz + Amwt2 * C2Wxz + Amwt3 * C3Wxz
!!
!! Compute incremental accelerations at nodal points requiring the spin
!! at the spot weld to be the same in all elements.
!!
            DO n = 1,4
              dAx(n  ) = dAx(n  )                                              &
     &                 + (Scale(1,3)*(Wxyspw - C1Wxy)) * C(n  ,3)              &
     &                 + (Scale(1,4)*(Wxzspw - C1Wxz)) * C(n  ,4)
              dAy(n  ) = dAy(n  )                                              &
     &                 - (Scale(1,2)*(Wxyspw - C1Wxy)) * C(n  ,2)              &
     &                 + (Scale(1,4)*(Wyzspw - C1Wyz)) * C(n  ,4)
              dAz(n  ) = dAz(n  )                                              &
     &                 - (Scale(1,3)*(Wyzspw - C1Wyz)) * C(n  ,3)              &
     &                 - (Scale(1,2)*(Wxzspw - C1Wxz)) * C(n  ,2)

              dAx(n+4) = dAx(n+4)                                              &
     &                 + (Scale(2,3)*(Wxyspw - C2Wxy)) * C(n+4,3)              &
     &                 + (Scale(2,4)*(Wxzspw - C2Wxz)) * C(n+4,4)
              dAy(n+4) = dAy(n+4)                                              &
     &                 - (Scale(2,2)*(Wxyspw - C2Wxy)) * C(n+4,2)              &
     &                 + (Scale(2,4)*(Wyzspw - C2Wyz)) * C(n+4,4)
              dAz(n+4) = dAz(n+4)                                              &
     &                 - (Scale(2,3)*(Wyzspw - C2Wyz)) * C(n+4,3)              &
     &                 - (Scale(2,2)*(Wxzspw - C2Wxz)) * C(n+4,2)

              dAx(n+8) = dAx(n+8)                                              &
     &                 + (Scale(3,3)*(Wxyspw - C3Wxy)) * C(n+8,3)              &
     &                 + (Scale(3,4)*(Wxzspw - C3Wxz)) * C(n+8,4)
              dAy(n+8) = dAy(n+8)                                              &
     &                 - (Scale(3,2)*(Wxyspw - C3Wxy)) * C(n+8,2)              &
     &                 + (Scale(3,4)*(Wyzspw - C3Wyz)) * C(n+8,4)
              dAz(n+8) = dAz(n+8)                                              &
     &                 - (Scale(3,3)*(Wyzspw - C3Wyz)) * C(n+8,3)              &
     &                 - (Scale(3,2)*(Wxzspw - C3Wxz)) * C(n+8,2)
            ENDDO
!!
!! Restore normal acceleration to the local vectors Ax, Ay and Az.
!!
            DO n = 1,12
              Ax(n) = Ax(n) + Tx * W(n)
              Ay(n) = Ay(n) + Ty * W(n)
              Az(n) = Az(n) + Tz * W(n)
            ENDDO
!!
!! Find the maximum "external" force generated by the spot weld constraint.
!! (This failure criterion does not distinguish between shear and normal.)
!! (It also does not distinguish between tension or compression.)
!!
            Fx = 0.0
            Fy = 0.0
            Fz = 0.0
            DO n = 1,4
              IF (IX(n) .GT. 0) THEN
                Fx = Fx + NODE(IX(n))%Mass * dAx(n)
                Fy = Fy + NODE(IX(n))%Mass * dAy(n)
                Fz = Fz + NODE(IX(n))%Mass * dAz(n)
              ENDIF
            ENDDO
            Ith_Constraint_Force = 0.0
            Fsquared = Fx*Fx +Fy*Fy +Fz*Fz
            IF (Fsquared .NE. 0.0) Ith_Constraint_Force = SQRT(Fsquared)
            Max_Constraint_Force = Ith_Constraint_Force
            Fx = 0.0
            Fy = 0.0
            Fz = 0.0
            DO n = 5,8
              IF (IX(n) .GT. 0) THEN
                Fx = Fx + NODE(IX(n))%Mass * dAx(n)
                Fy = Fy + NODE(IX(n))%Mass * dAy(n)
                Fz = Fz + NODE(IX(n))%Mass * dAz(n)
              ENDIF
            ENDDO
            Ith_Constraint_Force = 0.0
            Fsquared = Fx*Fx +Fy*Fy +Fz*Fz
            IF (Fsquared.NE.0.0) Ith_Constraint_Force = SQRT(Fsquared)
            Max_Constraint_Force =                                             &
     &        MAX (Max_Constraint_Force, Ith_Constraint_Force)
            Fx = 0.0
            Fy = 0.0
            Fz = 0.0
            DO n = 9,12
              IF (IX(n) .GT. 0) THEN
                Fx = Fx + NODE(IX(n))%Mass * dAx(n)
                Fy = Fy + NODE(IX(n))%Mass * dAy(n)
                Fz = Fz + NODE(IX(n))%Mass * dAz(n)
              ENDIF
            ENDDO
            Ith_Constraint_Force = 0.0
            Fsquared = Fx*Fx +Fy*Fy +Fz*Fz
            IF (Fsquared.NE.0.0) Ith_Constraint_Force = SQRT(Fsquared)
            Max_Constraint_Force =                                             &
     &        MAX (Max_Constraint_Force, Ith_Constraint_Force)

            SPOT_WELD(NSW)%Force = Max_Constraint_Force
!!
!! Check for a tensile failure of the spot weld.
!!
            FAILED = (Fmax.GT.0.0) .AND.                                       &
     &        (Max_Constraint_Force.GE.Fmax)
            SPOT_WELD(NSW)%FAILED = FAILED

            IF (.NOT.SPOT_WELD(NSW)%FAILED) THEN
!!
!! At long last, the incremental accelerations have been constructed which
!! will insure that the spot weld point in all elements will move in unison.
!! The incremental acceleration does not "remove" the individual element's
!! ability to strain.
!!
              DO n = 1,12
                IF (IX(n) .GT. 0) THEN
                  MOTION(IX(n))%Ax = (Ax(n) + dAx(n))
                  MOTION(IX(n))%Ay = (Ay(n) + dAy(n))
                  MOTION(IX(n))%Az = (Az(n) + dAz(n))
                  FORCE(IX(n))%Xext = NODE(IX(n))%Mass *                       &
     &              (Ax(n) + dAx(n)) + FORCE(IX(n))%Xint
                  FORCE(IX(n))%Yext = NODE(IX(n))%Mass *                       &
     &              (Ay(n) + dAy(n)) + FORCE(IX(n))%Yint
                  FORCE(IX(n))%Zext = NODE(IX(n))%Mass *                       &
     &              (Az(n) + dAz(n)) + FORCE(IX(n))%Zint
                ENDIF
              ENDDO
            ENDIF
!!
!! Spot weld position if-check.
!!
          ENDIF
!!
!! Spot weld not-failed if-check.
!!
        ENDIF
!!
!! Spot weld do-loop
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_RIGID_WALL_BC
!!
!! Copyright (c) by KEY Associates, 27-MAY-1991 13:42:41
!!
!! Purpose: This is a "look-ahead" calculation. If a nodal point in the next
!! time step crosses into the rigid wall exclusion domain, reaction forces
!! are added which will cause the nodal point in question to end up located on
!! the rigid wall. If the nodal point "leaves" the rigid wall, no action is
!! taken. If the rigid wall is inertial, that is, a heavy mass that moves as
!! a result of the impact, the wall reaction forces are used to accelerate the
!! wall center of mass. At the next time step the wall will actually be in a
!! different spot than the "look-ahead" calculation used for computing the
!! impact forces. As a result, an inertial rigid wall must have a mass that
!! is significantly larger than the impacting object to produce an accurate
!! result.
!!
      USE shared_common_data
      USE rigid_wall_bc_
      USE node_set_
      USE motion_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          Kode,                                                          &
     &          Code,                                                          &
     &          HstID,                                                         &
     &          Kavd
      REAL(KIND(0D0))                                                          &
     &          C(3),                                                          &
     &          Cn(3),                                                         &
     &          PRX(3,3),                                                      &
     &          RTX(3,3),                                                      &
     &          TABLE_LOOK_UP
      LOGICAL                                                                  &
     &          ON_WALL,                                                       &
     &          ON_WALL_W,                                                     &
     &          ON_WALL_H,                                                     &
     &          NEXT_NP_ID

      REAL(KIND(0D0)), PARAMETER :: DTR = 1.745329252D-2

      REAL(KIND(0D0)), PARAMETER :: C2F = (1.0D+0 /   2.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C3F = (1.0D+0 /   6.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C4F = (1.0D+0 /  24.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C5F = (1.0D+0 / 120.0D+0)
      REAL(KIND(0D0)), PARAMETER :: C6F = (1.0D+0 / 720.0D+0)

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Initialize rigid wall boundary condition data structure.
!!
      IF (FIRST .AND. CONTROL%RDSTAR .NE. 1) THEN
!!
!! Convert wall points P2 and P3 into unit vectors pointing from P1 to P2
!! and pointing from P1 to P3. Construct the wall normal vector Cn point-
!! ing "into" the wall; Cn = - (P2-P1) x (P3-P1).
!!
        DO N = 1,NUMWC
          Ax = RIGID_WALL_BC(N)%P2(1) - RIGID_WALL_BC(N)%P1(1)
          Ay = RIGID_WALL_BC(N)%P2(2) - RIGID_WALL_BC(N)%P1(2)
          Az = RIGID_WALL_BC(N)%P2(3) - RIGID_WALL_BC(N)%P1(3)
          Amag = SQRT (Ax*Ax + Ay*Ay + Az*Az)
          Ax = Ax * (ONE / Amag)
          Ay = Ay * (ONE / Amag)
          Az = Az * (ONE / Amag)
          RIGID_WALL_BC(N)%P2(1) = Ax
          RIGID_WALL_BC(N)%P2(2) = Ay
          RIGID_WALL_BC(N)%P2(3) = Az
          Bx = RIGID_WALL_BC(N)%P3(1) - RIGID_WALL_BC(N)%P1(1)
          By = RIGID_WALL_BC(N)%P3(2) - RIGID_WALL_BC(N)%P1(2)
          Bz = RIGID_WALL_BC(N)%P3(3) - RIGID_WALL_BC(N)%P1(3)
          Bmag = SQRT (Bx*Bx + By*By + Bz*Bz)
          Bx = Bx * (ONE / Bmag)
          By = By * (ONE / Bmag)
          Bz = Bz * (ONE / Bmag)
          RIGID_WALL_BC(N)%P3(1) = Bx
          RIGID_WALL_BC(N)%P3(2) = By
          RIGID_WALL_BC(N)%P3(3) = Bz
          Cx = Ay*Bz - By*Az
          Cy = Az*Bx - Bz*Ax
          Cz = Ax*By - Bx*Ay
          Cmag = SQRT (Cx*Cx + Cy*Cy + Cz*Cz)
          RIGID_WALL_BC(N)%Cn(1) = -Cx * (ONE / Cmag)
          RIGID_WALL_BC(N)%Cn(2) = -Cy * (ONE / Cmag)
          RIGID_WALL_BC(N)%Cn(3) = -Cz * (ONE / Cmag)
        ENDDO
!!
!! Bring acceleration/velocity/displacement flag Kavd into range.
!!
        DO N = 1,NUMWC
          RIGID_WALL_BC(N)%Kavd = MAX (1,MIN (3,RIGID_WALL_BC(N)%Kavd))
        ENDDO
!!
!! If the inertia tensor B(1,1) is non-zero, invert it.
!!
        DO N = 1,NUMWC
          IF (RIGID_WALL_BC(N)%B(1,1) .GT. 0.0) THEN
            Bxx = RIGID_WALL_BC(N)%B(1,1)
            Bxy = RIGID_WALL_BC(N)%B(2,1)
            Bxz = RIGID_WALL_BC(N)%B(3,1)
            Byy = RIGID_WALL_BC(N)%B(2,2)
            Byz = RIGID_WALL_BC(N)%B(3,2)
            Bzz = RIGID_WALL_BC(N)%B(3,3)
            Det = (Bxx*Byy)*Bzz                                                &
     &          + (Bxy*Byz)*(Bxz+Bxz)                                          &
     &          - (Bxz*Bxz)*Byy                                                &
     &          - (Byz*Byz)*Bxx                                                &
     &          - (Bxy*Bxy)*Bzz
            IF (ABS (Det) .LT. 1.0D-30) THEN
              WRITE (MSG1,'(I8)') RIGID_WALL_BC(N)%RWID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_RIGID_WALL_BC.001.00'//                          &
     &          MSGL//'Inertia Tensor Is Singular.'//                          &
     &          MSGL//'For Rigid Wall ID:'//MSG1                               &
     &          )
            ENDIF
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
            RIGID_WALL_BC(N)%C(1,1) = Cxx
            RIGID_WALL_BC(N)%C(2,1) = Cxy
            RIGID_WALL_BC(N)%C(3,1) = Cxz
            RIGID_WALL_BC(N)%C(1,2) = Cxy
            RIGID_WALL_BC(N)%C(2,2) = Cyy
            RIGID_WALL_BC(N)%C(3,2) = Cyz
            RIGID_WALL_BC(N)%C(1,3) = Cxz
            RIGID_WALL_BC(N)%C(2,3) = Cyz
            RIGID_WALL_BC(N)%C(3,3) = Czz
          ENDIF
        ENDDO
!!
!! Normalize direction vector A.
!!
        DO N = 1,NUMWC
          Amag = SQRT                                                          &
     &          (                                                              &
     &          RIGID_WALL_BC(N)%Ax * RIGID_WALL_BC(N)%Ax +                    &
     &          RIGID_WALL_BC(N)%Ay * RIGID_WALL_BC(N)%Ay +                    &
     &          RIGID_WALL_BC(N)%Az * RIGID_WALL_BC(N)%Az                      &
     &          )
          IF (Amag .EQ. 0.0) Amag = ONE
          RIGID_WALL_BC(N)%Ax = RIGID_WALL_BC(N)%Ax * (ONE / Amag)
          RIGID_WALL_BC(N)%Ay = RIGID_WALL_BC(N)%Ay * (ONE / Amag)
          RIGID_WALL_BC(N)%Az = RIGID_WALL_BC(N)%Az * (ONE / Amag)
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all rigid wall boundary condition specifications.
!! RIGID_WALL_BC,RWID,SetID,P1,P2,P3,Cn,Width,Height,CoF,Kode,Mass,B,
!! Disp,Vel,Accel,Omega,Force,Torque,Code,HstID,Kavd,Scale,Ax,Ay,Az
!!
      DTaver = 0.5D+0 * (TIMSIM%DTlast + TIMSIM%DTnext)
      Dniv   = ONE / TIMSIM%DTnext
      Daiv   = ONE / DTaver
!!
      DO NWC = 1,NUMWC
        SetID = RIGID_WALL_BC(NWC)%SetID
        Cof   = RIGID_WALL_BC(NWC)%CoF
        Kode  = RIGID_WALL_BC(NWC)%Kode
        DO i = 1,3
          Cn(i) = RIGID_WALL_BC(NWC)%Cn(i)
          RIGID_WALL_BC(NWC)%Force(i) = 0.0
          RIGID_WALL_BC(NWC)%Torque(i) = 0.0
        ENDDO
        Px = RIGID_WALL_BC(NWC)%P1(1) + RIGID_WALL_BC(NWC)%Disp(1)
        Py = RIGID_WALL_BC(NWC)%P1(2) + RIGID_WALL_BC(NWC)%Disp(2)
        Pz = RIGID_WALL_BC(NWC)%P1(3) + RIGID_WALL_BC(NWC)%Disp(3)
        Ax = RIGID_WALL_BC(NWC)%P2(1)
        Ay = RIGID_WALL_BC(NWC)%P2(2)
        Az = RIGID_WALL_BC(NWC)%P2(3)
        Bx = RIGID_WALL_BC(NWC)%P3(1)
        By = RIGID_WALL_BC(NWC)%P3(2)
        Bz = RIGID_WALL_BC(NWC)%P3(3)
!!
!! Note velocity of rigid wall at time n-1/2 and project position of rigid
!! wall at time n+1.
!!
        Vwx = RIGID_WALL_BC(NWC)%Vel(1)
        Vwy = RIGID_WALL_BC(NWC)%Vel(2)
        Vwz = RIGID_WALL_BC(NWC)%Vel(3)
        Pwx = Px + TIMSIM%DTnext * Vwx
        Pwy = Py + TIMSIM%DTnext * Vwy
        Pwz = Pz + TIMSIM%DTnext * Vwz
!!
!! Process all nodes in set.
!!
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Retrieve position of nodal point at time n.
!!
          Xm = MOTION(N)%Px + MOTION(N)%Ux
          Ym = MOTION(N)%Py + MOTION(N)%Uy
          Zm = MOTION(N)%Pz + MOTION(N)%Uz
!!
!! Compute velocity of nodal point at time n+1/2.
!!
          Vx = MOTION(N)%Vx + DTaver * MOTION(N)%Ax
          Vy = MOTION(N)%Vy + DTaver * MOTION(N)%Ay
          Vz = MOTION(N)%Vz + DTaver * MOTION(N)%Az
!!
!! Compute position of nodal point at time n+1.
!!
          Xn = Xm + TIMSIM%DTnext * Vx
          Yn = Ym + TIMSIM%DTnext * Vy
          Zn = Zm + TIMSIM%DTnext * Vz
!!
!! Compute depth of penetration by taking inner product between wall normal
!! vector Cn and position vector from the current location of wall reference
!! point, P1(1:3)+Disp(1:3), to nodal position (Xn,Yn,Zn). Note: Wall normal
!! vector points into the wall.
!!
          Depth = Cn(1)*(Xn-Pwx) + Cn(2)*(Yn-Pwy) + Cn(3)*(Zn-Pwz)
!!
!! If penetration has occured, compute force necessary at time n to place
!! node at the wall at time n+1. Recompute acceleration at time n.
!!
          IF (Depth .GT. 0.0) THEN
!!
!! Locate point of contact with rigid wall.
!!
            Alpha = Cn(1)*(Xn-Xm) + Cn(2)*(Yn-Ym) + Cn(3)*(Zn-Zm)
            Xw = Xm + (Alpha - Depth) * (Xn - Xm)
            Yw = Ym + (Alpha - Depth) * (Yn - Ym)
            Zw = Zm + (Alpha - Depth) * (Zn - Zm)
!!
!! Check wall width.
!!
            IF (RIGID_WALL_BC(NWC)%Width .GT. 0.0) THEN
              DA = Ax*(Xw-Pwx) + Ay*(Yw-Pwy) + Az*(Zw-Pwz)
              ON_WALL_W = (DA+DA) .LE. RIGID_WALL_BC(NWC)%Width
            ELSE
              ON_WALL_W = .TRUE.
            ENDIF
!!
!! Check wall height.
!!
            IF (RIGID_WALL_BC(NWC)%Height .GT. 0.0) THEN
              DB = Bx*(Xw-Pwx) + By*(Yw-Pwy) + Bz*(Zw-Pwz)
              ON_WALL_H = (DB+DB) .LE. RIGID_WALL_BC(NWC)%Height
            ELSE
              ON_WALL_H = .TRUE.
            ENDIF
!!
            ON_WALL = ON_WALL_W .AND. ON_WALL_H
!!
!! For those nodal points striking the wall alter their motion to reflect
!! their impact.
!!
           IF (ON_WALL) THEN
!!
!! Compute velocity of nodal point normal to the wall at time n+1/2.
!! !!! Must include wall rotation !!!
!!
              Vn = Cn(1)*(Vx-Vwx) + Cn(2)*(Vy-Vwy) + Cn(3)*(Vz-Vwz)
!!
!! Compute incremental acceleration dAm at time n needed to place nodal point
!! at the wall at time n+1.
!!
              dAm = -Depth*Dniv*Daiv
              Fm  = NODE(N)%Mass*dAm
!!
!! Check for non-zero friction coefficient.
!!
              IF (CoF .GT. 0.0) THEN
!!
!! Compute relative tangential velocity of nodal point at time n+1/2.
!!
                Vtx  = (Vx-Vwx) - Vn * Cn(1)
                Vty  = (Vy-Vwy) - Vn * Cn(1)
                Vtz  = (Vz-Vwz) - Vn * Cn(1)
                Vtan = SQRT (Vtx*Vtx + Vty*Vty + Vtz*Vtz)
                Ftan = Shear_Traction (CoF,ABS(Fm),Vtan)
                dAt  = Ftan * NODE(N)%Minv
              ELSE
                Vtx  = 0.0
                Vty  = 0.0
                Vtz  = 0.0
                Vtan = ONE
                Ftan = 0.0
                dAt  = 0.0
              ENDIF
!!
!! Add incremental acceleration dA to acceleration of nodal point at time n.
!!
              MOTION(N)%Ax = MOTION(N)%Ax + dAm*Cn(1)-(dAt/Vtan)*Vtx
              MOTION(N)%Ay = MOTION(N)%Ay + dAm*Cn(2)-(dAt/Vtan)*Vty
              MOTION(N)%Az = MOTION(N)%Az + dAm*Cn(3)-(dAt/Vtan)*Vtz
!!
!! Compute incremental forces acting on nodal point.
!!
              Fx = Fm * Cn(1) - (Ftan/Vtan) * Vtx
              Fy = Fm * Cn(2) - (Ftan/Vtan) * Vty
              Fz = Fm * Cn(3) - (Ftan/Vtan) * Vtz
!!
!! Modify external force acting on nodal point at time n to produce incre-
!! mental acceleration.
!!
              FORCE(N)%Xext = FORCE(N)%Xext + Fx
              FORCE(N)%Yext = FORCE(N)%Yext + Fy
              FORCE(N)%Zext = FORCE(N)%Zext + Fz
!!
!! Accumulate incremental force to obtain total reaction force on wall.
!!
              RIGID_WALL_BC(NWC)%Force(1) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(1) - Fx
              RIGID_WALL_BC(NWC)%Force(2) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(2) - Fy
              RIGID_WALL_BC(NWC)%Force(3) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(3) - Fz
!!
!! Compute position vector from P1 ("wall center point") to nodal point
!! position on the wall.
!!
              Dx = Xw - Pwx
              Dy = Yw - Pwy
              Dz = Zw - Pwz
!!
!! Compute torque at P1 needed to produce incremental external force acting
!! on nodal point at time n.
!!
              Tx = -(Dy*Fz - Fy*Dz)
              Ty = -(Dz*Fx - Fz*Dx)
              Tz = -(Dx*Fy - Fx*Dy)
!!
!! Accumulate incremental torque to obtain total reaction torque on wall.
!!
              RIGID_WALL_BC(NWC)%Torque(1) =                                   &
     &          RIGID_WALL_BC(NWC)%Torque(1) + Tx
              RIGID_WALL_BC(NWC)%Torque(2) =                                   &
     &          RIGID_WALL_BC(NWC)%Torque(2) + Ty
              RIGID_WALL_BC(NWC)%Torque(3) =                                   &
     &          RIGID_WALL_BC(NWC)%Torque(3) + Tz
            ENDIF
          ENDIF
        ENDDO
!!
!! Check for a rigid wall that is allowed to move. If rigid wall is inertial
!! or has a prescribed motion, reposition wall for next time step, time n+1.
!! (Kode = 0/1/2 = rigid/inertial/prescribed movement)
!!
        IF (Kode .GT. 0) THEN
!!
!! Compute translational acceleration for wall center-of-mass (located
!! at P1+Disp).
!!
          Ax = (ONE/RIGID_WALL_BC(NWC)%Mass)*RIGID_WALL_BC(NWC)%Force(1)
          Ay = (ONE/RIGID_WALL_BC(NWC)%Mass)*RIGID_WALL_BC(NWC)%Force(2)
          Az = (ONE/RIGID_WALL_BC(NWC)%Mass)*RIGID_WALL_BC(NWC)%Force(3)
!!
!! Retrieve velocity and displacement of wall center-of-mass.
!!
          Vx = RIGID_WALL_BC(NWC)%Vel(1)
          Vy = RIGID_WALL_BC(NWC)%Vel(2)
          Vz = RIGID_WALL_BC(NWC)%Vel(2)
          Ux = RIGID_WALL_BC(NWC)%Disp(1)
          Uy = RIGID_WALL_BC(NWC)%Disp(2)
          Uz = RIGID_WALL_BC(NWC)%Disp(2)
!!
!! Check for kinematic boundary conditions on motion of wall center-of-mass.
!!
          IF (Kode .EQ. 2) THEN
            Code  = RIGID_WALL_BC(NWC)%Code
            HstID = RIGID_WALL_BC(NWC)%HstID
            Kavd  = RIGID_WALL_BC(NWC)%Kavd
            Scale = RIGID_WALL_BC(NWC)%Scale
            Dx    = RIGID_WALL_BC(NWC)%Ax
            Dy    = RIGID_WALL_BC(NWC)%Ay
            Dz    = RIGID_WALL_BC(NWC)%Az
!!
!! Evaluate history function.
!!
            IF (Kavd .EQ. 1) THEN
              Tx = TIMSIM%Total
            ELSE IF (Kavd .EQ. 2) THEN
              Tx = TIMSIM%Total + 0.5D+0 * TIMSIM%DTnext
            ELSE IF (Kavd .EQ. 3) THEN
              Tx = TIMSIM%Total + TIMSIM%DTnext
            ENDIF
            Ft = Scale * TABLE_LOOK_UP (HstID,Tx)
!!
!! Apply constraint.
!!
            IF (Code .EQ. 1) THEN
              IF (Kavd .EQ. 1) THEN
                Ax = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Ax = (Ft - Vx)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Ax = ((Ft - Ux)*Dniv - Vx)*Daiv
              ENDIF
              RIGID_WALL_BC(NWC)%Force(1) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(1) + Ax*RIGID_WALL_BC(NWC)%Mass
            ELSE IF (Code .EQ. 2) THEN
              IF (Kavd .EQ. 1) THEN
                Ay = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Ay = (Ft - Vy)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Ay = ((Ft - Uy)*Dniv - Vy)*Daiv
              ENDIF
              RIGID_WALL_BC(NWC)%Force(2) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(2) + Ay*RIGID_WALL_BC(NWC)%Mass
            ELSE IF (Code .EQ. 3) THEN
              IF (Kavd .EQ. 1) THEN
                Az = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                Az = (Ft - Vz)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Az = ((Ft - Uz)*Dniv - Vz)*Daiv
              ENDIF
              RIGID_WALL_BC(NWC)%Force(3) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(3) + Az*RIGID_WALL_BC(NWC)%Mass
            ELSE IF (Code .EQ. 4) THEN
              IF (Kavd .EQ. 1) THEN
                DA = Dx*Ax + Dy*Ay + Dz*Az
                Pa = Ft
              ELSE IF (Kavd .EQ. 2) THEN
                DA = Dx*Ax + Dy*Ay + Dz*Az
                DV = Dx*Vx + Dy*Vy + Dz*Vz
                Pa = (Ft - DV)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                DA = Dx*Ax + Dy*Ay + Dz*Az
                DU = Dx*Ux + Dy*Uy + Dz*Uz
                DV = Dx*Vx + Dy*Vy + Dz*Vz
                Pa = ((Ft - DU)*Dniv - DV)*Daiv
              ENDIF
              RIGID_WALL_BC(NWC)%Force(1) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(1) +                                  &
     &          Pa*Dx * RIGID_WALL_BC(NWC)%Mass
              RIGID_WALL_BC(NWC)%Force(2) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(2) +                                  &
     &          Pa*Dy * RIGID_WALL_BC(NWC)%Mass
              RIGID_WALL_BC(NWC)%Force(3) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(3) +                                  &
     &          Pa*Dz * RIGID_WALL_BC(NWC)%Mass
              Ax = Ax - (DA - Pa)*Dx
              Ay = Ay - (DA - Pa)*Dy
              Az = Az - (DA - Pa)*Dz
            ELSE IF (Code .EQ. 10) THEN
              Ay = 0.0
              Az = 0.0
            ELSE IF (Code .EQ. 20) THEN
              Ax = 0.0
              Az = 0.0
            ELSE IF (Code .EQ. 30) THEN
              Ax = 0.0
              Ay = 0.0
            ELSE IF (Code .EQ. 40) THEN
              Pa = Dx*Ax + Dy*Ay + Dz*Az
              Ax = Pa*Dx
              Ay = Pa*Dy
              Az = Pa*Dz
            ELSE
              IF (Code .EQ. 11) THEN
                Fx = Ft
                Fy = 0.0
                Fz = 0.0
              ELSE IF (Code .EQ. 22) THEN
                Fx = 0.0
                Fy = Ft
                Fz = 0.0
              ELSE IF (Code .EQ. 33) THEN
                Fx = 0.0
                Fy = 0.0
                Fz = Ft
              ELSE IF (Code .EQ. 44) THEN
                Fx = Dx*Ft
                Fy = Dy*Ft
                Fz = Dz*Ft
              ELSE
                WRITE (MSG1,'(I8)') RIGID_WALL_BC(NWC)%RWID
                WRITE (MSG2,'(I8)') Code
                CALL USER_MESSAGE                                              &
     &            (                                                            &
     &            MSGL//'FATAL'//                                              &
     &            MSGL//'IMPOSE_RIGID_WALL_BC.002.00'//                        &
     &            MSGL//'WALLBC Input Record ID:'//MSG1//                      &
     &            MSGL//'Contains Unknown BC Code:'//MSG2                      &
     &            )
              ENDIF
              IF (Kavd .EQ. 1) THEN
                Px = Fx
                Py = Fy
                Pz = Fz
              ELSE IF (Kavd .EQ. 2) THEN
                Px = (Fx - Vx)*Daiv
                Py = (Fy - Vy)*Daiv
                Pz = (Fz - Vz)*Daiv
              ELSE IF (Kavd .EQ. 3) THEN
                Px = ((Fx - Ux)*Dniv - Vx)*Daiv
                Py = ((Fy - Uy)*Dniv - Vy)*Daiv
                Pz = ((Fz - Uz)*Dniv - Vz)*Daiv
              ENDIF
              RIGID_WALL_BC(NWC)%Force(1) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(1) +                                  &
     &          Px * RIGID_WALL_BC(NWC)%Mass
              RIGID_WALL_BC(NWC)%Force(2) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(2) +                                  &
     &          Py * RIGID_WALL_BC(NWC)%Mass
              RIGID_WALL_BC(NWC)%Force(3) =                                    &
     &          RIGID_WALL_BC(NWC)%Force(3) +                                  &
     &          Pz * RIGID_WALL_BC(NWC)%Mass
              Ax = Px
              Ay = Py
              Az = Pz
            ENDIF
          ENDIF
!!
!! Save accelerations and integrate center-of-mass equations of motion.
!!
          RIGID_WALL_BC(NWC)%Accel(1) = Ax
          RIGID_WALL_BC(NWC)%Accel(2) = Ay
          RIGID_WALL_BC(NWC)%Accel(3) = Az
          RIGID_WALL_BC(NWC)%Vel(1) = (Vx + DTaver * Ax)
          RIGID_WALL_BC(NWC)%Vel(2) = (Vy + DTaver * Ay)
          RIGID_WALL_BC(NWC)%Vel(3) = (Vz + DTaver * Az)
          RIGID_WALL_BC(NWC)%Disp(1) =                                         &
     &      Ux + TIMSIM%DTnext * (Vx + DTaver * Ax)
          RIGID_WALL_BC(NWC)%Disp(2) =                                         &
     &      Uy + TIMSIM%DTnext * (Vy + DTaver * Ay)
          RIGID_WALL_BC(NWC)%Disp(3) =                                         &
     &      Uz + TIMSIM%DTnext * (Vz + DTaver * Az)
!!
!! Compute angular accelerations about wall center-of-mass using Euler's
!! equation; dO/dt = Binv (T - OxB*O)
!!
          IF (RIGID_WALL_BC(NWC)%B(1,1) .GT. 0.0) THEN
            Bxx = RIGID_WALL_BC(NWC)%B(1,1)
            Bxy = RIGID_WALL_BC(NWC)%B(2,1)
            Bxz = RIGID_WALL_BC(NWC)%B(3,1)
            Byy = RIGID_WALL_BC(NWC)%B(2,2)
            Byz = RIGID_WALL_BC(NWC)%B(3,2)
            Bzz = RIGID_WALL_BC(NWC)%B(3,3)
            Ox = RIGID_WALL_BC(NWC)%Omega(1)
            Oy = RIGID_WALL_BC(NWC)%Omega(2)
            Oz = RIGID_WALL_BC(NWC)%Omega(3)
            Q1 = Bxx*Ox + Bxy*Oy + Bxz*Oz
            Q2 = Bxy*Ox + Byy*Oy + Byz*Oz
            Q3 = Bxz*Ox + Byz*Oy + Bzz*Oz
            Tx = RIGID_WALL_BC(NWC)%Torque(1)
            Ty = RIGID_WALL_BC(NWC)%Torque(2)
            Tz = RIGID_WALL_BC(NWC)%Torque(3)
            Tx = Tx + Q3*Oy - Q2*Oz
            Ty = Ty + Q1*Oz - Q3*Ox
            Tz = Tz + Q2*Ox - Q1*Oy
            Cxx = RIGID_WALL_BC(NWC)%C(1,1)
            Cxy = RIGID_WALL_BC(NWC)%C(2,1)
            Cxz = RIGID_WALL_BC(NWC)%C(3,1)
            Cyy = RIGID_WALL_BC(NWC)%C(2,2)
            Cyz = RIGID_WALL_BC(NWC)%C(3,2)
            Czz = RIGID_WALL_BC(NWC)%C(3,3)
            Oxdot = (Cxx*Tx + Cxy*Ty + Cxz*Tz)
            Oydot = (Cxy*Tx + Cyy*Ty + Cyz*Tz)
            Ozdot = (Cxz*Tx + Cyz*Ty + Czz*Tz)
          ELSE
            Oxdot = 0.0
            Oydot = 0.0
            Ozdot = 0.0
          ENDIF
!!
!! If angular acceleration is non-zero, make rotational updates.
!!
          IF ((ABS(Oxdot)+ABS(Oydot)+ABS(Ozdot)) .GT. 0.0) THEN
!!
!! Compute new angular velocity of the center-of-mass.
!!
            Ox = Ox + DTaver * Oxdot
            Oy = Oy + DTaver * Oydot
            Oz = Oz + DTaver * Ozdot
!!
!! Update angular velocity of the center-of-mass.
!!
            RIGID_WALL_BC(NWC)%Omega(1) = Ox
            RIGID_WALL_BC(NWC)%Omega(2) = Oy
            RIGID_WALL_BC(NWC)%Omega(3) = Oz
!!
!! Construct unit vector C aligned with angular velocity "vector."
!!
            Qmag = SQRT (Ox*Ox + Oy*Oy + Oz*Oz)
            C(1) = Ox * (ONE / Qmag)
            C(2) = Oy * (ONE / Qmag)
            C(3) = Oz * (ONE / Qmag)
!!
!! Compute angle of rotation and evaluate the trigonometric functions
!! SinPhi = Sin(Phi) and CosPm1 = Cos(Phi) - 1.
!!
            Phi = TIMSIM%DTnext * Qmag
            Ph2 = Phi*Phi
            CosPm1 = -Ph2*(C2F - Ph2*(C4F - C6F*Ph2))
            SinPhi =  Phi*(ONE - Ph2*(C3F - C5F*Ph2))
!!
!! Construct rotation operator RTX. R = (Cos(phi)-1.0)*(I-CC) + Sin(phi)*(C x .)
!! The operator RTX is an "incremental" proper-orthogonal rotation. It is
!! constructed to be an "incremental operator" to avoid round-off.
!!
            DO n = 1,3
              DO m = 1,3
                RTX(m,n) = 0.0
                PRX(m,n) = -C(m) * C(n)
              ENDDO
            ENDDO
            DO n = 1,3
              PRX(n,n) = ONE + PRX(n,n)
            ENDDO
            RTX(1,2) = -(SinPhi*C(3))
            RTX(1,3) =  (SinPhi*C(2))
            RTX(2,1) =  (SinPhi*C(3))
            RTX(2,3) = -(SinPhi*C(1))
            RTX(3,1) = -(SinPhi*C(2))
            RTX(3,2) =  (SinPhi*C(1))
            DO n = 1,3
              DO m = 1,3
                RTX(m,n) = CosPm1*PRX(m,n) + RTX(m,n)
              ENDDO
            ENDDO
!!
!! Update orientation of wall normal vector Cn.
!!
            C(1) = RIGID_WALL_BC(NWC)%Cn(1)
            C(2) = RIGID_WALL_BC(NWC)%Cn(2)
            C(3) = RIGID_WALL_BC(NWC)%Cn(3)
            DO n = 1,3
              dC = RTX(n,1)*C(1) + RTX(n,2)*C(2)+ RTX(n,3)*C(3)
              RIGID_WALL_BC(NWC)%Cn(n) = RIGID_WALL_BC(NWC)%Cn(n) + dC
            ENDDO
!!
!! Update orientation of wall vector A = (Ax,Ay,Az).
!!
            C(1) = RIGID_WALL_BC(NWC)%P2(1)
            C(2) = RIGID_WALL_BC(NWC)%P2(2)
            C(3) = RIGID_WALL_BC(NWC)%P2(3)
            DO n = 1,3
              dP = RTX(n,1)*C(1) + RTX(n,2)*C(2)+ RTX(n,3)*C(3)
              RIGID_WALL_BC(NWC)%P2(n) = RIGID_WALL_BC(NWC)%P2(n) + dP
            ENDDO
!!
!! Update orientation of wall vector B = (Bx,By,Bz).
!!
            C(1) = RIGID_WALL_BC(NWC)%P3(1)
            C(2) = RIGID_WALL_BC(NWC)%P3(2)
            C(3) = RIGID_WALL_BC(NWC)%P3(3)
            DO n = 1,3
              dP = RTX(n,1)*C(1) + RTX(n,2)*C(2)+ RTX(n,3)*C(3)
              RIGID_WALL_BC(NWC)%P3(n) = RIGID_WALL_BC(NWC)%P3(n) + dP
            ENDDO
!!
!! Update components of inertia tensor from time n to n+1 to acccount for
!! rotational change of the congifuration wrt the spatial coordinates.
!! First, augment incremental rotation operator with identity operator.
!!
            RTX(1,1) = ONE + RTX(1,1)
            RTX(2,2) = ONE + RTX(2,2)
            RTX(3,3) = ONE + RTX(3,3)

            Bxx = RIGID_WALL_BC(NWC)%B(1,1)
            Bxy = RIGID_WALL_BC(NWC)%B(2,1)
            Bxz = RIGID_WALL_BC(NWC)%B(3,1)
            Byy = RIGID_WALL_BC(NWC)%B(2,2)
            Byz = RIGID_WALL_BC(NWC)%B(3,2)
            Bzz = RIGID_WALL_BC(NWC)%B(3,3)
            Txx = RTX(1,1)*Bxx + RTX(1,2)*Bxy + RTX(1,3)*Bxz + Bxx
            Tyx = RTX(2,1)*Bxx + RTX(2,2)*Bxy + RTX(2,3)*Bxz + Bxy
            Tzx = RTX(3,1)*Bxx + RTX(3,2)*Bxy + RTX(3,3)*Bxz + Bxz
            Txy = RTX(1,1)*Bxy + RTX(1,2)*Byy + RTX(1,3)*Byz + Bxy
            Tyy = RTX(2,1)*Bxy + RTX(2,2)*Byy + RTX(2,3)*Byz + Byy
            Tzy = RTX(3,1)*Bxy + RTX(3,2)*Byy + RTX(3,3)*Byz + Byz
            Txz = RTX(1,1)*Bxz + RTX(1,2)*Byz + RTX(1,3)*Bzz + Bxz
            Tyz = RTX(2,1)*Bxz + RTX(2,2)*Byz + RTX(2,3)*Bzz + Byz
            Tzz = RTX(3,1)*Bxz + RTX(3,2)*Byz + RTX(3,3)*Bzz + Bzz
            Bxx = Txx*RTX(1,1) + Txy*RTX(1,2) + Txz*RTX(1,3) + Txx
            Bxy = Txx*RTX(2,1) + Txy*RTX(2,2) + Txz*RTX(2,3) + Txy
            Bxz = Txx*RTX(3,1) + Txy*RTX(3,2) + Txz*RTX(3,3) + Txz
            Byy = Tyx*RTX(2,1) + Tyy*RTX(2,2) + Tyz*RTX(2,3) + Tyy
            Byz = Tyx*RTX(3,1) + Tyy*RTX(3,2) + Tyz*RTX(3,3) + Tyz
            Bzz = Tzx*RTX(3,1) + Tzy*RTX(3,2) + Tzz*RTX(3,3) + Tzz
            RIGID_WALL_BC(NWC)%B(1,1) = Bxx
            RIGID_WALL_BC(NWC)%B(2,1) = Bxy
            RIGID_WALL_BC(NWC)%B(3,1) = Bxz
            RIGID_WALL_BC(NWC)%B(1,2) = Bxy
            RIGID_WALL_BC(NWC)%B(2,2) = Byy
            RIGID_WALL_BC(NWC)%B(3,2) = Byz
            RIGID_WALL_BC(NWC)%B(1,3) = Bxz
            RIGID_WALL_BC(NWC)%B(2,3) = Byz
            RIGID_WALL_BC(NWC)%B(3,3) = Bzz
!!
!! Invert inertia tensor.
!!
            Det = (Bxx*Byy)*Bzz                                                &
     &          + (Bxy*Byz)*(Bxz+Bxz)                                          &
     &          - (Bxz*Bxz)*Byy                                                &
     &          - (Byz*Byz)*Bxx                                                &
     &          - (Bxy*Bxy)*Bzz
            Cxx = ( Byy*Bzz  - (Byz*Byz)) * (ONE / Det)
            Cyy = ( Bxx*Bzz  - (Bxz*Bxz)) * (ONE / Det)
            Czz = ((Bxx*Byy) - (Bxy*Bxy)) * (ONE / Det)
            Cxy = ( Bxz*Byz  -  Bxy*Bzz ) * (ONE / Det)
            Cxz = ((Bxy*Byz) -  Bxz*Byy ) * (ONE / Det)
            Cyz = ( Bxy*Bxz  -  Byz*Bxx ) * (ONE / Det)
            RIGID_WALL_BC(NWC)%C(1,1) = Cxx
            RIGID_WALL_BC(NWC)%C(2,1) = Cxy
            RIGID_WALL_BC(NWC)%C(3,1) = Cxz
            RIGID_WALL_BC(NWC)%C(1,2) = Cxy
            RIGID_WALL_BC(NWC)%C(2,2) = Cyy
            RIGID_WALL_BC(NWC)%C(3,2) = Cyz
            RIGID_WALL_BC(NWC)%C(1,3) = Cxz
            RIGID_WALL_BC(NWC)%C(2,3) = Cyz
            RIGID_WALL_BC(NWC)%C(3,3) = Czz
          ENDIF
!!
!! End of inertial wall if-test.
!!
        ENDIF
!!
!! End of rigid wall BC do-loop
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_PERIODIC_BC
!!
!! Copyright (c) by KEY Associates; 19-APR-1994 19:52:13.24
!!
!! Purpose: Apply periodic B.C. to individual nodal points on
!! the bounding faces.
!!
      USE shared_common_data
      USE periodic_bc_
      USE node_set_
      USE motion_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID
      REAL(KIND(0D0))                                                          &
     &          RTX(3,3),                                                      &
     &          Fx1,Fy1,Fz1,                                                   &
     &          Fx2,Fy2,Fz2,                                                   &
     &          Gx1,Gy1,Gz1,                                                   &
     &          Gx2,Gy2,Gz2
!!
!! Loop over all periodic boundary condition specifications.
!! PERIODIC_BC,PerID,NPSET=S1ID,NPSET=S2ID,Type,Ax,Ay,Az,Px,Py,Pz,Theta,Advance
!!
      DO NCC = 1,NUMCC
!!
!! Linear Periodicity --
!!
        IF (PERIODIC_BC(NCC)%Type .EQ. 'LINEAR') THEN
!!
!! Process node pairs, one from each face.
!!
          SetID = PERIODIC_BC(NCC)%S1ID
          N1 = NODE_SET(SetID)%Istart - 1
          SetID = PERIODIC_BC(NCC)%S2ID
          N2 = NODE_SET(SetID)%Istart - 1
          DO i = NODE_SET(SetID)%Istart,NODE_SET(SetID)%Iend
            N1 = N1 + 1
            N2 = N2 + 1
            NP1 = NNPSETS(N1)
            NP2 = NNPSETS(N2)
!!
!! Compute common accelerations based on the loads applied on the
!! respective surfaces.
!!
            Fx1 = FORCE(NP1)%Xext - FORCE(NP1)%Xint
            Fy1 = FORCE(NP1)%Yext - FORCE(NP1)%Yint
            Fz1 = FORCE(NP1)%Zext - FORCE(NP1)%Zint
            Fx2 = FORCE(NP2)%Xext - FORCE(NP2)%Xint
            Fy2 = FORCE(NP2)%Yext - FORCE(NP2)%Yint
            Fz2 = FORCE(NP2)%Zext - FORCE(NP2)%Zint

            Ax = (Fx1 + Fx2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Ay = (Fy1 + Fy2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Az = (Fz1 + Fz2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

            MOTION(NP1)%Ax  = Ax
            MOTION(NP1)%Ay  = Ay
            MOTION(NP1)%Az  = Az
            FORCE(NP1)%Xext = Ax * NODE(NP1)%Mass + FORCE(NP1)%Xint
            FORCE(NP1)%Yext = Ay * NODE(NP1)%Mass + FORCE(NP1)%Yint
            FORCE(NP1)%Zext = Az * NODE(NP1)%Mass + FORCE(NP1)%Zint

            MOTION(NP2)%Ax  = Ax
            MOTION(NP2)%Ay  = Ay
            MOTION(NP2)%Az  = Az
            FORCE(NP2)%Xext = Ax * NODE(NP2)%Mass + FORCE(NP2)%Xint
            FORCE(NP2)%Yext = Ay * NODE(NP2)%Mass + FORCE(NP2)%Yint
            FORCE(NP2)%Zext = Az * NODE(NP2)%Mass + FORCE(NP2)%Zint
!!
!! Check to see if this nodal point has rotational degress of freedom.
!!
            IF (NODE(NP1)%IRT .GT. 0) THEN
              NP1 = NODE(NP1)%IRT
              NP2 = NODE(NP2)%IRT
!!
!! Compute common accelerations based on the torques applied on the
!! respective surfaces.
!!
              Fx1 = FORCE(NP1)%Xext - FORCE(NP1)%Xint
              Fy1 = FORCE(NP1)%Yext - FORCE(NP1)%Yint
              Fz1 = FORCE(NP1)%Zext - FORCE(NP1)%Zint
              Fx2 = FORCE(NP2)%Xext - FORCE(NP2)%Xint
              Fy2 = FORCE(NP2)%Yext - FORCE(NP2)%Yint
              Fz2 = FORCE(NP2)%Zext - FORCE(NP2)%Zint

              Ax = (Fx1 + Fx2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Ay = (Fy1 + Fy2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Az = (Fz1 + Fz2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

              MOTION(NP1)%Ax  = Ax
              MOTION(NP1)%Ay  = Ay
              MOTION(NP1)%Az  = Az
              FORCE(NP1)%Xext = Ax * NODE(NP1)%Mass + FORCE(NP1)%Xint
              FORCE(NP1)%Yext = Ay * NODE(NP1)%Mass + FORCE(NP1)%Yint
              FORCE(NP1)%Zext = Az * NODE(NP1)%Mass + FORCE(NP1)%Zint

              MOTION(NP2)%Ax  = Ax
              MOTION(NP2)%Ay  = Ay
              MOTION(NP2)%Az  = Az
              FORCE(NP2)%Xext = Ax * NODE(NP2)%Mass + FORCE(NP2)%Xint
              FORCE(NP2)%Yext = Ay * NODE(NP2)%Mass + FORCE(NP2)%Yint
              FORCE(NP2)%Zext = Az * NODE(NP2)%Mass + FORCE(NP2)%Zint
            ENDIF
          ENDDO
!!
!! Cyclic Periodicity --
!!
        ELSE IF (PERIODIC_BC(NCC)%Type .EQ. 'CYCLIC') THEN
!!
!! Retrieve orthogonal rotation constructed during initialization.
!!
          DO i = 1,9
            RTX(i,1) = PERIODIC_BC(NCC)%RTX(i,1)
          ENDDO
!!
!! Process node pairs, one from each face.
!!
          SetID = PERIODIC_BC(NCC)%S1ID
          N1 = NODE_SET(SetID)%Istart - 1
          SetID = PERIODIC_BC(NCC)%S2ID
          N2 = NODE_SET(SetID)%Istart - 1
          DO i = NODE_SET(SetID)%Istart,NODE_SET(SetID)%Iend
            N1 = N1 + 1
            N2 = N2 + 1
            NP1 = NNPSETS(N1)
            NP2 = NNPSETS(N2)
!!
!! Compute common accelerations based on the loads applied on the
!! respective surfaces.
!!
            Fx1 = FORCE(NP1)%Xext - FORCE(NP1)%Xint
            Fy1 = FORCE(NP1)%Yext - FORCE(NP1)%Yint
            Fz1 = FORCE(NP1)%Zext - FORCE(NP1)%Zint
            Fx2 = FORCE(NP2)%Xext - FORCE(NP2)%Xint
            Fy2 = FORCE(NP2)%Yext - FORCE(NP2)%Yint
            Fz2 = FORCE(NP2)%Zext - FORCE(NP2)%Zint
!!
!! Rotate force on face 1 to face 2.
!!
            Gx1 = Fx1 + RTX(1,1)*Fx1 + RTX(1,2)*Fy1 + RTX(1,3)*Fz1
            Gy1 = Fy1 + RTX(2,1)*Fx1 + RTX(2,2)*Fy1 + RTX(2,3)*Fz1
            Gz1 = Fz1 + RTX(3,1)*Fx1 + RTX(3,2)*Fy1 + RTX(3,3)*Fz1
!!
!! Update acceleration and external force on face 2.
!!
            Ax  = (Fx2 + Gx1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Ay  = (Fy2 + Gy1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Az  = (Fz2 + Gz1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

            MOTION(NP2)%Ax  = Ax
            MOTION(NP2)%Ay  = Ay
            MOTION(NP2)%Az  = Az
            FORCE(NP2)%Xext = Ax * NODE(NP2)%Mass + FORCE(NP2)%Xint
            FORCE(NP2)%Yext = Ay * NODE(NP2)%Mass + FORCE(NP2)%Yint
            FORCE(NP2)%Zext = Az * NODE(NP2)%Mass + FORCE(NP2)%Zint
!!
!! Rotate Force on face 2 to face 1.
!!
            Gx2 = Fx2 + RTX(1,1)*Fx2 + RTX(2,1)*Fy2 + RTX(3,1)*Fz2
            Gy2 = Fy2 + RTX(1,2)*Fx2 + RTX(2,2)*Fy2 + RTX(3,2)*Fz2
            Gz2 = Fz2 + RTX(1,3)*Fx2 + RTX(2,3)*Fy2 + RTX(3,3)*Fz2
!!
!! Update acceleration and external force on face 1.
!!
            Ax = (Fx1 + Gx2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Ay = (Fy1 + Gy2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
            Az = (Fz1 + Gz2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

            MOTION(NP1)%Ax  = Ax
            MOTION(NP1)%Ay  = Ay
            MOTION(NP1)%Az  = Az
            FORCE(NP1)%Xext = Ax * NODE(NP1)%Mass + FORCE(NP1)%Xint
            FORCE(NP1)%Yext = Ay * NODE(NP1)%Mass + FORCE(NP1)%Yint
            FORCE(NP1)%Zext = Az * NODE(NP1)%Mass + FORCE(NP1)%Zint
!!
!! Check to see if this nodal point has rotational degress of freedom.
!!
            IF (NODE(NP1)%IRT .GT. 0) THEN
              NP1 = NODE(NP1)%IRT
              NP2 = NODE(NP2)%IRT
!!
!! Compute common accelerations based on the torques applied on the
!! respective surfaces.
!!
              Fx1 = FORCE(NP1)%Xext - FORCE(NP1)%Xint
              Fy1 = FORCE(NP1)%Yext - FORCE(NP1)%Yint
              Fz1 = FORCE(NP1)%Zext - FORCE(NP1)%Zint
              Fx2 = FORCE(NP2)%Xext - FORCE(NP2)%Xint
              Fy2 = FORCE(NP2)%Yext - FORCE(NP2)%Yint
              Fz2 = FORCE(NP2)%Zext - FORCE(NP2)%Zint
!!
!! Rotate torque on face 1 to face 2.
!!
              Gx1 = Fx1 + RTX(1,1)*Fx1 + RTX(1,2)*Fy1 + RTX(1,3)*Fz1
              Gy1 = Fy1 + RTX(2,1)*Fx1 + RTX(2,2)*Fy1 + RTX(2,3)*Fz1
              Gz1 = Fz1 + RTX(3,1)*Fx1 + RTX(3,2)*Fy1 + RTX(3,3)*Fz1
!!
!! Update acceleration and external torque on face 2.
!!
              Ax  = (Fx2 + Gx1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Ay  = (Fy2 + Gy1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Az  = (Fz2 + Gz1) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

              MOTION(NP2)%Ax  = Ax
              MOTION(NP2)%Ay  = Ay
              MOTION(NP2)%Az  = Az
              FORCE(NP2)%Xext = Ax * NODE(NP2)%Mass + FORCE(NP2)%Xint
              FORCE(NP2)%Yext = Ay * NODE(NP2)%Mass + FORCE(NP2)%Yint
              FORCE(NP2)%Zext = Az * NODE(NP2)%Mass + FORCE(NP2)%Zint
!!
!! Rotate Torque on face 2 to face 1.
!!
              Gx2 = Fx2 + RTX(1,1)*Fx2 + RTX(2,1)*Fy2 + RTX(3,1)*Fz2
              Gy2 = Fy2 + RTX(1,2)*Fx2 + RTX(2,2)*Fy2 + RTX(3,2)*Fz2
              Gz2 = Fz2 + RTX(1,3)*Fx2 + RTX(2,3)*Fy2 + RTX(3,3)*Fz2
!!
!! Update acceleration and external torque on face 1.
!!
              Ax = (Fx1 + Gx2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Ay = (Fy1 + Gy2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)
              Az = (Fz1 + Gz2) / (NODE(NP1)%Mass + NODE(NP2)%Mass)

              MOTION(NP1)%Ax  = Ax
              MOTION(NP1)%Ay  = Ay
              MOTION(NP1)%Az  = Az
              FORCE(NP1)%Xext = Ax * NODE(NP1)%Mass + FORCE(NP1)%Xint
              FORCE(NP1)%Yext = Ay * NODE(NP1)%Mass + FORCE(NP1)%Yint
              FORCE(NP1)%Zext = Az * NODE(NP1)%Mass + FORCE(NP1)%Zint
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_BODY_FORCE
!!
!! Copyright (c) by KEY Associates, 24-DEC-1991 10:19:09
!!
!! Purpose: Apply body forces to individual nodal points.
!!
      USE shared_common_data
      USE body_force_
      USE node_set_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          TABLE_LOOK_UP
      LOGICAL                                                                  &
     &          NEXT_NP_ID

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Normalize direction vector A.
!!
        DO N = 1,NUMBF
          Amag = SQRT                                                          &
     &          (                                                              &
     &          BODY_FORCE(N)%Ax * BODY_FORCE(N)%Ax +                          &
     &          BODY_FORCE(N)%Ay * BODY_FORCE(N)%Ay +                          &
     &          BODY_FORCE(N)%Az * BODY_FORCE(N)%Az                            &
     &          )
          IF (Amag .EQ. 0.0) THEN
            WRITE (MSG1,'(I8)') BODY_FORCE(N)%BFID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'IMPOSE_BODY_FORCE.001.00'//                             &
     &          MSGL//'BODYFORCE Input Record ID:'//MSG1//                     &
     &          MSGL//'Has A Null Direction Vector,'//                         &
     &          MSGL//'No Body Forces Will Be Generated.'                      &
     &          )
            Amag = ONE
          ENDIF
          BODY_FORCE(N)%Ax = BODY_FORCE(N)%Ax * (ONE / Amag)
          BODY_FORCE(N)%Ay = BODY_FORCE(N)%Ay * (ONE / Amag)
          BODY_FORCE(N)%Az = BODY_FORCE(N)%Az * (ONE / Amag)
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all body force specifications.
!! BODY_FORCE,BFID,SetID,HstID,Scale,Gravity,Ax,Ay,Az,Delay
!!
      DO NBF = 1,NUMBF
        SetID   = BODY_FORCE(NBF)%SetID
        HstID   = BODY_FORCE(NBF)%HstID
        Scale   = BODY_FORCE(NBF)%Scale
        Gravity = BODY_FORCE(NBF)%Gravity
        Ax      = BODY_FORCE(NBF)%Ax
        Ay      = BODY_FORCE(NBF)%Ay
        Az      = BODY_FORCE(NBF)%Az
!!
!! Evaluate history function.
!!
        Tx = TIMSIM%Total - BODY_FORCE(NBF)%Delay
        Ft = Gravity * Scale * TABLE_LOOK_UP (HstID,Tx)
!!
!! Store components of nodal body force in global arrays. Note: A zero force
!! is a do-nothing and, thus, is skipped.
!!
        IF (Ft .NE. 0.0) THEN
!!
!! Process all nodes in set.
!!
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Add body force to existing external forces at the nodal point.
!!
            FORCE(N)%Xext = FORCE(N)%Xext + (Ax * Ft) * NODE(N)%Mass
            FORCE(N)%Yext = FORCE(N)%Yext + (Ay * Ft) * NODE(N)%Mass
            FORCE(N)%Zext = FORCE(N)%Zext + (Az * Ft) * NODE(N)%Mass
          ENDDO
        ENDIF
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_PRESSURE_BC
!!
!! Copyright (c) by KEY Associates, 27-APR-1991 15:45:23
!!
!! Purpose: Compute consistent nodal forces from pressure applied normal to
!! element surface.
!!
      USE shared_common_data
      USE pressure_bc_
      USE motion_
      USE force_
      USE segment_
      USE segment_set_
      USE enumerated_sets_, ONLY: NSGSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          HstID
      REAL(KIND(0D0))                                                          &
     &          TABLE_LOOK_UP,                                                 &
     &          Cross_Dot_Product,                                             &
     &          Xel(4),Yel(4),Zel(4),                                          &
     &          A(3),B(3),C(3),D(3),                                           &
     &          R(3),S(3),G(4),Qn(3)
      LOGICAL                                                                  &
     &          NEXT_SEG_ID,                                                   &
     &          QUADRILATERAL

      REAL(KIND(0D0)), PARAMETER :: ONE3RD = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: ONE6TH = (1.0D+0 / 6.0D+0)
      REAL(KIND(0D0)), PARAMETER :: ONE8TH = (1.0D+0 / 8.0D+0)
!!
!! Loop over all pressure boundary condition specifications.
!! PRESSURE_BC,PBCID,SetID,HstID,Scale,PI,PJ,PK,PL,Delay
!!

      DO NPC = 1,NUMPC
        SetID = PRESSURE_BC(NPC)%SetID
        HstID = PRESSURE_BC(NPC)%HstID
        Scale = PRESSURE_BC(NPC)%Scale
!!
!! Evaluate history function.
!!
        Tx = TIMSIM%Total - PRESSURE_BC(NPC)%Delay
        Ft = Scale * TABLE_LOOK_UP (HstID,Tx)
!!
!! A zero pressure is a do-nothing and, thus, is skipped.
!!
        IF (Ft .NE. 0.0) THEN
!!
!! Process all nodes in set.
!!
          N = 0
          DO WHILE (NEXT_SEG_ID(SetID,N))
!!
!! Distinguish between a quadrilateral segment and a triangular segment.
!!
            QUADRILATERAL = (SEGMENT(N)%PAR%IX(4) .NE. 0)
            IF (QUADRILATERAL) THEN
!!
!! Gather quadrilateral segment nodal coordinates.
!!
              DO i = 1,4
                NPID = SEGMENT(N)%PAR%IX(i)
                Xel(i) = MOTION(NPID)%Px + MOTION(NPID)%Ux
                Yel(i) = MOTION(NPID)%Py + MOTION(NPID)%Uy
                Zel(i) = MOTION(NPID)%Pz + MOTION(NPID)%Uz
              ENDDO
!!
!! Define vectors A through C along loaded segment edges.
!!
              A(1) = Xel(2)-Xel(1)
              A(2) = Yel(2)-Yel(1)
              A(3) = Zel(2)-Zel(1)
              B(1) = Xel(3)-Xel(2)
              B(2) = Yel(3)-Yel(2)
              B(3) = Zel(3)-Zel(2)
              C(1) = Xel(4)-Xel(3)
              C(2) = Yel(4)-Yel(3)
              C(3) = Zel(4)-Zel(3)
              D(1) = Xel(1)-Xel(4)
              D(2) = Yel(1)-Yel(4)
              D(3) = Zel(1)-Zel(4)
!!
!! Define vectors R and S across element diagonals.
!!
              R(1) = Xel(3)-Xel(1)
              R(2) = Yel(3)-Yel(1)
              R(3) = Zel(3)-Zel(1)
              S(1) = Xel(4)-Xel(2)
              S(2) = Yel(4)-Yel(2)
              S(3) = Zel(4)-Zel(2)
!!
!! Define Qn normal to loaded element face, Qn = R x S. The magnitude of Qn
!! is twice the area of the quadralateral defined by R and S. Thus, it con-
!! verts the stress to a total force.
!!
              Qn(1) = R(2)*S(3)-S(2)*R(3)
              Qn(2) = R(3)*S(1)-S(3)*R(1)
              Qn(3) = R(1)*S(2)-S(1)*R(2)
!!
!! Compute components of surface traction force.
!!
              Pcntr = PRESSURE_BC(NPC)%PI + PRESSURE_BC(NPC)%PJ                &
     &              + PRESSURE_BC(NPC)%PK + PRESSURE_BC(NPC)%PL
              Qx = (ONE8TH*Ft*Pcntr)*Qn(1)
              Qy = (ONE8TH*Ft*Pcntr)*Qn(2)
              Qz = (ONE8TH*Ft*Pcntr)*Qn(3)
!!
!! Compute nodal reactions for element face to be loaded which are in static
!! equilibrium with Pcntr and which satisfy the constraints:
!!
!!              Pi + Pj + Pk + Pl = Pcntr,  Equilibrium,
!!              Pi - Pj + Pk - Pl = 0,      Anti-hourglassing,
!!              Pi * Li + Pk * Lk = 0,      Zero moment,
!!              Pj * Lj + Pl * Ll = 0,      Zero moment.
!!
              AI = Cross_Dot_Product (D,S,Qn)
              AJ = Cross_Dot_Product (A,R,Qn)
              AK = Cross_Dot_Product (B,S,Qn)
              AL = Cross_Dot_Product (C,R,Qn)
              G(4) = -AJ*(0.5D+0/(AL-AJ))
              G(3) = -AI*(0.5D+0/(AK-AI))
              G(2) =  AL*(0.5D+0/(AL-AJ))
              G(1) =  AK*(0.5D+0/(AK-AI))
!!
!! Store local values of nodal loads in global arrays. Note: The pressure
!! acts opposite to the direction Qn, hence, the minus sign in the accumu-
!! lation.
!!
              DO i = 1,4
                NPID = SEGMENT(N)%PAR%IX(i)
                FORCE(NPID)%Xext = FORCE(NPID)%Xext - Qx*G(i)
                FORCE(NPID)%Yext = FORCE(NPID)%Yext - Qy*G(i)
                FORCE(NPID)%Zext = FORCE(NPID)%Zext - Qz*G(i)
              ENDDO
!!
            ELSE
!!
!! Gather triangular segment nodal coordinates.
!!
              DO i = 1,3
                NPID = SEGMENT(N)%PAR%IX(i)
                Xel(i) = MOTION(NPID)%Px + MOTION(NPID)%Ux
                Yel(i) = MOTION(NPID)%Py + MOTION(NPID)%Uy
                Zel(i) = MOTION(NPID)%Pz + MOTION(NPID)%Uz
              ENDDO
!!
!! Calculate element area and normal vector. Define element edge vectors.
!!
              Rx = Xel(2)-Xel(1)
              Ry = Yel(2)-Yel(1)
              Rz = Zel(2)-Zel(1)
              Sx = Xel(3)-Xel(1)
              Sy = Yel(3)-Yel(1)
              Sz = Zel(3)-Zel(1)
!!
!! Define the unit vector T normal to the element. The magnitude of T is
!! twice the area of the triangular segment.
!!
              Tx = Ry*Sz - Sy*Rz
              Ty = Rz*Sx - Sz*Rx
              Tz = Rx*Sy - Sx*Ry
!!
!! Compute nodal components of surface traction force.
!!
              Pmultiplier = ONE3RD *                                           &
     &          (                                                              &
     &          PRESSURE_BC(NPC)%PI +                                          &
     &          PRESSURE_BC(NPC)%PJ +                                          &
     &          PRESSURE_BC(NPC)%PK                                            &
     &          )
              Qx = (ONE6TH * Pmultiplier * Ft) * Tx
              Qy = (ONE6TH * Pmultiplier * Ft) * Ty
              Qz = (ONE6TH * Pmultiplier * Ft) * Tz
!!
!! Store local values of nodal loads in global arrays. Note: The pressure
!! acts opposite to the direction T, hence, the minus sign in the accumu-
!! lation.
!!
              DO i = 1,3
                NPID = SEGMENT(N)%PAR%IX(i)
                FORCE(NPID)%Xext = FORCE(NPID)%Xext - Qx
                FORCE(NPID)%Yext = FORCE(NPID)%Yext - Qy
                FORCE(NPID)%Zext = FORCE(NPID)%Zext - Qz
              ENDDO
!!
!! End of quadrilateral/triangle if-test.
!!
            ENDIF
!!
!! End of segment set do-while.
!!
          ENDDO
!!
!! End of zero pressure if-test.
!!
        ENDIF
!!
!! End of pressure BC do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_FORCE_BC
!!
!! Copyright (c) by KEY Associates, 6-APR-1991 16:06:40
!!
!! Purpose: Apply force B.C. to individual nodal points.
!!
      USE shared_common_data
      USE force_bc_
      USE node_set_
      USE motion_
      USE force_
      USE node_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          HstID,                                                         &
     &          Type,                                                          &
     &          Follow
      REAL(KIND(0D0))                                                          &
     &          TABLE_LOOK_UP
      LOGICAL                                                                  &
     &          NEXT_NP_ID,                                                    &
     &          Moment

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Bring force/moment flag Type into range.
!!
        DO N = 1,NUMFC
          FORCE_BC(N)%Type = MAX (0,MIN (1,FORCE_BC(N)%Type))
        ENDDO
!!
!! Check for concentrated moments on non-shell nodal points.
!!
        DO NFC = 1,NUMFC
          SetID = FORCE_BC(NFC)%SetID
          Moment = (FORCE_BC(NFC)%Type .EQ. 1)
!!
!! If this is a moment boundary condition for a rotational degree of
!! freedom, examine all nodes in the set to insure they are nodes with
!! rotational degrees of freedom.
!!
          IF (Moment) THEN
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              IF (NODE(N)%IRT .EQ. 0) THEN
                WRITE (MSG1,'(I8)') FORCE_BC(NFC)%CFID
                WRITE (MSG2,'(I8)') NODE(N)%ID
                CALL USER_MESSAGE                                              &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'IMPOSE_FORCE_BC.001.00'//                                     &
     &    MSGL//'FORCEBC Input Record ID:'//MSG1//                             &
     &    MSGL//'Specifies A Moment BC On Nodal Point:'//MSG2//                &
     &    MSGL//'Which Does Not Have Rotational Degrees Of Freedom.'           &
     &    )
              ENDIF
            ENDDO
          ENDIF
        ENDDO
!!
!! Normalize direction vector C.
!!
        DO N = 1,NUMFC
          Cmag = SQRT                                                          &
     &          (                                                              &
     &          FORCE_BC(N)%Cx * FORCE_BC(N)%Cx +                              &
     &          FORCE_BC(N)%Cy * FORCE_BC(N)%Cy +                              &
     &          FORCE_BC(N)%Cz * FORCE_BC(N)%Cz                                &
     &          )
          IF (Cmag .EQ. 0.0) Cmag = ONE
          FORCE_BC(N)%Cx = FORCE_BC(N)%Cx * (ONE / Cmag)
          FORCE_BC(N)%Cy = FORCE_BC(N)%Cy * (ONE / Cmag)
          FORCE_BC(N)%Cz = FORCE_BC(N)%Cz * (ONE / Cmag)
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all concentrated force boundary condition specifications.
!! FORCE_BC,CFID,SetID,HstID,Type,Follow,Force,Cx,Cy,Cz,Delay
!!
      DO NFC = 1,NUMFC
        SetID  = FORCE_BC(NFC)%SetID
        HstID  = FORCE_BC(NFC)%HstID
        Follow = FORCE_BC(NFC)%Follow
        Cx     = FORCE_BC(NFC)%Cx
        Cy     = FORCE_BC(NFC)%Cy
        Cz     = FORCE_BC(NFC)%Cz
        Moment = (FORCE_BC(NFC)%Type .EQ. 1)
!!
!! Evaluate history function.
!!
        Tx = TIMSIM%Total - FORCE_BC(NFC)%Delay
        Ft = FORCE_BC(NFC)%Force * TABLE_LOOK_UP (HstID,Tx)
!!
!! Store components of concentrated nodal force in global arrays.
!! Note: A zero force or moment is a do-nothing and, thus, is skipped.
!!
        IF (Ft .NE. 0.0) THEN
!!
!! Process all nodes in set.
!!
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Discriminate between a force and a moment. Locate rotational degree
!! of freedom.
!!
            IF (Moment) N = NODE(N)%IRT
!!
!! Add force/moment to existing external forces at the nodal point.
!!
            FORCE(N)%Xext = FORCE(N)%Xext + Cx * Ft
            FORCE(N)%Yext = FORCE(N)%Yext + Cy * Ft
            FORCE(N)%Zext = FORCE(N)%Zext + Cz * Ft
          ENDDO
        ENDIF
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_SPRING_BC
!!
!! Copyright (c) by KEY Associates, 27-JUL-1991 13:06:22
!!
!! Purpose: Apply spring B.C. to individual nodal points.
!!
      USE shared_common_data
      USE spring_bc_
      USE material_
      USE motion_
      USE force_
      USE node_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          MatID,                                                         &
     &          Type,                                                          &
     &          Follow
      REAL(KIND(0D0))                                                          &
     &          Length
      LOGICAL                                                                  &
     &          MOMENT,                                                        &
     &          NEXT_NP_ID

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Bring force/moment flag Type and follower flag Follow into range.
!!
        DO NSC = 1,NUMSC
          SPRING_BC(NSC)%Type   = MAX (0,MIN (1,SPRING_BC(NSC)%Type  ))
          SPRING_BC(NSC)%Follow = MAX (0,MIN (1,SPRING_BC(NSC)%Follow))
          SPRING_BC(NSC)%RES%Delta = 0.0
        ENDDO
!!
!! Check for concentrated moments on non-shell nodal points.
!!
        DO NSC = 1,NUMSC
          SetID = SPRING_BC(NSC)%SetID
          MOMENT = (SPRING_BC(NSC)%Type .EQ. 1)
!!
!! Process all nodes in set.
!!
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! If this is a moment boundary condition for a rotational degree of freedom,
!! examine node to insure it is a node with rotational degrees of freedom.
!!
            IF (MOMENT .AND. NODE(N)%IRT .EQ. 0) THEN
              WRITE (MSG1,'(I8)') SPRING_BC(NSC)%SprID
              WRITE (MSG2,'(I8)') NODE(N)%ID
              CALL USER_MESSAGE                                                &
     &   (                                                                     &
     &   MSGL//'FATAL'//                                                       &
     &   MSGL//'IMPOSE_SPRING_BC.001.01'//                                     &
     &   MSGL//'Spring BC (SPRINGBC) Input Record ID:'//MSG1//                 &
     &   MSGL//'Specifies A Moment BC On Nodal Point ID:'//MSG2//              &
     &   MSGL//'Which Does Not Have Rotational Degrees Of Freedom.'            &
     &   )
            ENDIF
          ENDDO
        ENDDO
!!
!! Check to make certain that NO node sets are used. (Hopefully this is
!! temporary.)
!!
        DO NSC = 1,NUMSC
          IF (SPRING_BC(NSC)%SetID .GT. 0) THEN
            WRITE (MSG1,'(I8)') SPRING_BC(NSC)%SprID
            WRITE (MSG2,'(I8)') SPRING_BC(NSC)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_SPRING_BC.001.02'//                              &
     &          MSGL//'Spring BC (SPRINGBC) Input Record ID:'//MSG1//          &
     &          MSGL//'References Nodal Set ID:'//MSG2//                       &
     &          MSGL//'Each Nodal Point Requires A Separate Spring BC.'        &
     &          )
          ENDIF
        ENDDO
!!
!! Normalize direction vector Axis(1:3).
!!
        DO NSC = 1,NUMSC
          Amag = SQRT                                                          &
     &          (                                                              &
     &          SPRING_BC(NSC)%Axis(1) * SPRING_BC(NSC)%Axis(1) +              &
     &          SPRING_BC(NSC)%Axis(2) * SPRING_BC(NSC)%Axis(2) +              &
     &          SPRING_BC(NSC)%Axis(3) * SPRING_BC(NSC)%Axis(3)                &
     &          )
          IF (Amag .EQ. 0.0) Amag = ONE
          SPRING_BC(NSC)%Axis(1) = SPRING_BC(NSC)%Axis(1) * (ONE / Amag)
          SPRING_BC(NSC)%Axis(2) = SPRING_BC(NSC)%Axis(2) * (ONE / Amag)
          SPRING_BC(NSC)%Axis(3) = SPRING_BC(NSC)%Axis(3) * (ONE / Amag)
        ENDDO
!!
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all spring boundary condition specifications.
!! SPRING_BC,SprID,SetID,MatID,Type,Follow,Axis(1:3),Force,Delta,Int_Eng
!!
      DO NSC = 1,NUMSC
        SetID  = SPRING_BC(NSC)%SetID
        MatID  = SPRING_BC(NSC)%MatID
        Type   = SPRING_BC(NSC)%Type
        Follow = SPRING_BC(NSC)%Follow
        Rx     = SPRING_BC(NSC)%Axis(1)
        Ry     = SPRING_BC(NSC)%Axis(2)
        Rz     = SPRING_BC(NSC)%Axis(3)
        MOMENT = (Type .EQ. 1)
!!
!! Process all nodes in set.
!!
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Discriminate between a force and a moment. Locate rotational degrees
!! of freedom.
!!
          IF (MOMENT) THEN
            NP = NODE(N)%IRT
          ELSE
            NP = N
          ENDIF
!!
!! The spring will act in one of two "directions" depending on the value
!! of *.Follow specified. Thus, if
!!
!!       Follow = 0, Spring force is aligned with the displacement.
!!              = 1, Spring force is aligned with the direction Axis(1:3).
!!
          IF (Follow .EQ. 0) THEN
            Rx = MOTION(NP)%Ux
            Ry = MOTION(NP)%Uy
            Rz = MOTION(NP)%Uz
            Length = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! Note: Axis(1:3) is used in the event the displacements are zero. Thus,
!! Axis(1:3) will be used on the very first time step (the initialization
!! pass at TIMSIM%Step = 0). For preloaded springs, that is, springs with force-
!! displacement relations which have a non-zero force at zero displacement,
!! it is important that Axis(1:3) be defined in order to impose the preload
!! during the initialization pass.
!!
            IF (Length .LT. 1.0D-25) THEN
              Rx = SPRING_BC(NSC)%Axis(1)
              Ry = SPRING_BC(NSC)%Axis(2)
              Rz = SPRING_BC(NSC)%Axis(3)
            ELSE
              Rx = Rx * (ONE / Length)
              Ry = Ry * (ONE / Length)
              Rz = Rz * (ONE / Length)
            ENDIF
          ELSE
            Rx = SPRING_BC(NSC)%Axis(1)
            Ry = SPRING_BC(NSC)%Axis(2)
            Rz = SPRING_BC(NSC)%Axis(3)
          ENDIF
!!
!! Relative axial velocity using Vx,Vy,Vz transformed to local R coordinate.
!!
          Drr = Rx*MOTION(NP)%Vx + Ry*MOTION(NP)%Vy + Rz*MOTION(NP)%Vz
!!
!! Update spring displacement.
!!
          dU = TIMSIM%DTnext*Drr
          SPRING_BC(NSC)%RES%Delta = SPRING_BC(NSC)%RES%Delta + dU
!!
!! Update force-displacement model to obtain new force.
!!
          IF (MATERIAL(MatID)%Type .EQ. 1) THEN
            CALL MATERIAL_S1                                                   &
     &          (                                                              &
     &          SPRING_BC(NSC)%RES%Force,                                      &
     &          STATE_VARIABLES(SPRING_BC(NSC)%Isv),                           &
     &          SPRING_BC(NSC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 2) THEN
            CALL MATERIAL_S2                                                   &
     &          (                                                              &
     &          SPRING_BC(NSC)%RES%Force,                                      &
     &          STATE_VARIABLES(SPRING_BC(NSC)%Isv),                           &
     &          SPRING_BC(NSC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 3) THEN
            CALL MATERIAL_S3                                                   &
     &          (                                                              &
     &          SPRING_BC(NSC)%RES%Force,                                      &
     &          STATE_VARIABLES(SPRING_BC(NSC)%Isv),                           &
     &          SPRING_BC(NSC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 5) THEN
            CALL MATERIAL_S5                                                   &
     &          (                                                              &
     &          SPRING_BC(NSC)%RES%Force,                                      &
     &          STATE_VARIABLES(SPRING_BC(NSC)%Isv),                           &
     &          SPRING_BC(NSC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE
            WRITE (MSG1,'(I8)') SPRING_BC(NSC)%SprID
            WRITE (MSG2,'(I8)') MATERIAL(MatID)%MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_SPRING_BC.001.03'//                              &
     &          MSGL//'SPRINGBC Input Record ID:'//MSG1//                      &
     &          MSGL//'References Material ID:'//MSG2//                        &
     &          MSGL//'With An Invalid Material Type.'                         &
     &          )
          ENDIF
!!
!! Add force/moment to existing external forces at the nodal point.
!!
          FORCE(NP)%Xext =                                                     &
     &    FORCE(NP)%Xext - Rx * SPRING_BC(NSC)%RES%Force
          FORCE(NP)%Yext =                                                     &
     &    FORCE(NP)%Yext - Ry * SPRING_BC(NSC)%RES%Force
          FORCE(NP)%Zext =                                                     &
     &    FORCE(NP)%Zext - Rz * SPRING_BC(NSC)%RES%Force
!!
!! Critical time step calculation. (This assumes that half the mass is
!! "associated" with the spring similar to what occurs with a truss element.
!! The formulae used are (1) dTcrit < 2/Omega and (2) Omega**2 = k/m.
!!
          Spring_Constant = MAX (ONE, SOUND_SPEED%RCL2)
          DTelt = SQRT ((NODE(NP)%Mass+NODE(NP)%Mass)/Spring_Constant)
!!
!! Compare spring time step with minimum and maximum spring element time step.
!!
          DTSBx = TIMSIM%DTSBx
          TIMSIM%DTSBx = MAX (DTSBx,DTelt)
          IF (DTelt .LT. TIMSIM%DTSBC(1)) THEN
            TIMSIM%DTSBC(1) = DTelt
            TIMSIM%SprBC(1) = NSC
          ENDIF
!!
!! End of node set do-while
!!
        ENDDO
!!
!! End of spring BC do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_DAMPER_BC
!!
!! Copyright (c) by KEY Associates, 27-JUL-1991 13:06:22
!!
!! Purpose: Apply damper B.C. to individual nodal points.
!!
      USE shared_common_data
      USE damper_bc_
      USE material_
      USE motion_
      USE force_
      USE node_
      USE node_set_
      USE state_variables_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      COMMON /SPRDMP/                                                          &
     &          Rx,Ry,Rz,Drr
!!
      INTEGER                                                                  &
     &          SetID,                                                         &
     &          MatID,                                                         &
     &          Type,                                                          &
     &          Follow
      REAL(KIND(0D0))                                                          &
     &          Length
      LOGICAL                                                                  &
     &          MOMENT,                                                        &
     &          NEXT_NP_ID

      LOGICAL, SAVE :: FIRST = .TRUE.
!!
      IF (FIRST) THEN
!!
!! Bring force/moment flag Type and follower flag Follow into range.
!!
        DO N = 1,NUMVC
          DAMPER_BC(N)%Type   = MAX (0,MIN (1,DAMPER_BC(N)%Type  ))
          DAMPER_BC(N)%Follow = MAX (0,MIN (2,DAMPER_BC(N)%Follow))
          DAMPER_BC(N)%RES%Delta = 0.0
        ENDDO
!!
!! Check for concentrated moments on non-shell nodal points.
!!
        DO NVC = 1,NUMVC
          SetID = DAMPER_BC(NVC)%SetID
          MOMENT = (DAMPER_BC(NVC)%Type .EQ. 1)
!!
!! Process all nodes in set.
!!
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! If this is a moment boundary condition for a rotational degree of freedom,
!! examine node to insure it is a node with rotational degrees of freedom.
!!
            IF (MOMENT .AND. NODE(N)%IRT .EQ. 0) THEN
              WRITE (MSG1,'(I8)') DAMPER_BC(NVC)%DprID
              WRITE (MSG2,'(I8)') NODE(N)%ID
              CALL USER_MESSAGE                                                &
     &   (                                                                     &
     &   MSGL//'FATAL'//                                                       &
     &   MSGL//'IMPOSE_DAMPER_BC.001.00'//                                     &
     &   MSGL//'Damper BC (DAMPERBC) Input Record ID:'//MSG1//                 &
     &   MSGL//'Specifies A Moment BC On Nodal Point ID:'//MSG2//              &
     &   MSGL//'Which Does Not Have Rotational Degrees Of Freedom.'            &
     &            )
            ENDIF
          ENDDO
        ENDDO
!!
!! Check to make certain damper BC's reference a valid matertial model.
!!
        DO N = 1,NUMVC
          IF                                                                   &
     &          (                                                              &
     &          MATERIAL(DAMPER_BC(N)%MatID)%Type .LT. 6                       &
     &          .OR.                                                           &
     &          MATERIAL(DAMPER_BC(N)%MatID)%Type .GT. 9                       &
     &          )                                                              &
     &      THEN
            WRITE (MSG1,'(I8)') DAMPER_BC(N)%DprID
            WRITE (MSG2,'(I8)') MATERIAL(DAMPER_BC(N)%MatID)%MatID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'IMPOSE_DAMPER_BC.001.00'//                              &
     &          MSGL//'Damper BC (DAMPERBC) Input Record ID:'//MSG1//          &
     &          MSGL//'References Material ID:'//MSG2//                        &
     &          MSGL//'With An Invalid Material Type.'                         &
     &          )
          ENDIF
        ENDDO
!!
!! Normalize direction vector Axis(1:3).
!!
        DO N = 1,NUMVC
          Amag = SQRT                                                          &
     &          (                                                              &
     &          DAMPER_BC(N)%Axis(1) * DAMPER_BC(N)%Axis(1) +                  &
     &          DAMPER_BC(N)%Axis(2) * DAMPER_BC(N)%Axis(2) +                  &
     &          DAMPER_BC(N)%Axis(3) * DAMPER_BC(N)%Axis(3)                    &
     &          )
          IF (Amag .EQ. 0.0) Amag = ONE
          DAMPER_BC(N)%Axis(1) = DAMPER_BC(N)%Axis(1) * (ONE / Amag)
          DAMPER_BC(N)%Axis(2) = DAMPER_BC(N)%Axis(2) * (ONE / Amag)
          DAMPER_BC(N)%Axis(3) = DAMPER_BC(N)%Axis(3) * (ONE / Amag)
        ENDDO
!!
        FIRST = .FALSE.
      ENDIF
!!
!! Loop over all damper boundary condition specifications.
!! DAMPER_BC,DprID,SetID,MatID,Type,Follow,Axis(3)
!!
      DO NVC = 1,NUMVC
        SetID  = DAMPER_BC(NVC)%SetID
        MatID  = DAMPER_BC(NVC)%MatID
        Type   = DAMPER_BC(NVC)%Type
        Follow = DAMPER_BC(NVC)%Follow
        Rx     = DAMPER_BC(NVC)%Axis(1)
        Ry     = DAMPER_BC(NVC)%Axis(2)
        Rz     = DAMPER_BC(NVC)%Axis(3)
        MOMENT = (DAMPER_BC(NVC)%Type .EQ. 1)
!!
!! Process all nodes in set. Note: This works because none of the damper
!! material models require any state variables to function.
!!
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!!
!! Discriminate between a force and a moment. Locate rotational degrees
!! of freedom.
!!
          IF (MOMENT) THEN
            NP = NODE(N)%IRT
          ELSE
            NP = N
          ENDIF
!!
!! The damper will act in one of three "directions" depending on the value
!! of *.Follow specified. Thus, if
!!
!!       Follow = 0, Damper force is aligned with the displacement.
!!              = 1, Damper force is aligned with the velocity.
!!              = 2, Damper force is aligned with the direction Axis(1:3).
!!
          IF (Follow .EQ. 0) THEN
            Rx = MOTION(NP)%Ux
            Ry = MOTION(NP)%Uy
            Rz = MOTION(NP)%Uz
            Length = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! Note: Axis(1:3) is used in the event the displacements are zero. Thus,
!! Axis(1:3) will used on the very first time step (the initialization pass at
!! TIMSIM%Step = 0). For "preloaded" dampers, that is, dampers attached to nodes
!! with a non-zero velocity IC (non-zero velocity at TIMSIM%Step = 0), it is
!! important that Axis(1:3) be defined in order to impose the preload during
!! the initialization pass.
!!
            IF (Length .LT. 1.0D-25) THEN
              Rx = DAMPER_BC(NVC)%Axis(1)
              Ry = DAMPER_BC(NVC)%Axis(2)
              Rz = DAMPER_BC(NVC)%Axis(3)
            ELSE
              Rx = Rx * (ONE / Length)
              Ry = Ry * (ONE / Length)
              Rz = Rz * (ONE / Length)
            ENDIF
          ELSE IF (Follow .EQ. 1) THEN
            Rx = MOTION(NP)%Vx
            Ry = MOTION(NP)%Vy
            Rz = MOTION(NP)%Vz
            Length = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! Note: Axis(1:3) is used in the event the velocities are zero.
!!
            IF (Length .LT. 1.0D-25) THEN
              Rx = DAMPER_BC(NVC)%Axis(1)
              Ry = DAMPER_BC(NVC)%Axis(2)
              Rz = DAMPER_BC(NVC)%Axis(3)
            ELSE
              Rx = Rx * (ONE / Length)
              Ry = Ry * (ONE / Length)
              Rz = Rz * (ONE / Length)
            ENDIF
          ELSE IF (Follow .EQ. 2) THEN
            Rx = DAMPER_BC(NVC)%Axis(1)
            Ry = DAMPER_BC(NVC)%Axis(2)
            Rz = DAMPER_BC(NVC)%Axis(3)
          ENDIF
!!
!! Relative axial velocity using Vx,Vy,Vz transformed to local R coordinate.
!!
          Drr = Rx*MOTION(NP)%Vx + Ry*MOTION(NP)%Vy + Rz*MOTION(NP)%Vz
!!
!! Update damper velocity.
!!
          DAMPER_BC(NVC)%RES%Delta = Drr
!!
!! Interogate force-velocity model to obtain new force.
!!
          IF (MATERIAL(MatID)%Type .EQ. 6) THEN
            CALL MATERIAL_D6                                                   &
     &          (                                                              &
     &          DAMPER_BC(NVC)%RES%Force,                                      &
     &          STATE_VARIABLES(DAMPER_BC(NVC)%Isv),                           &
     &          DAMPER_BC(NVC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 7) THEN
            CALL MATERIAL_D7                                                   &
     &          (                                                              &
     &          DAMPER_BC(NVC)%RES%Force,                                      &
     &          STATE_VARIABLES(DAMPER_BC(NVC)%Isv),                           &
     &          DAMPER_BC(NVC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 8) THEN
            CALL MATERIAL_D8                                                   &
     &          (                                                              &
     &          DAMPER_BC(NVC)%RES%Force,                                      &
     &          STATE_VARIABLES(DAMPER_BC(NVC)%Isv),                           &
     &          DAMPER_BC(NVC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ELSE IF (MATERIAL(MatID)%Type .EQ. 9) THEN
            CALL MATERIAL_D9                                                   &
     &          (                                                              &
     &          DAMPER_BC(NVC)%RES%Force,                                      &
     &          STATE_VARIABLES(DAMPER_BC(NVC)%Isv),                           &
     &          DAMPER_BC(NVC)%RES%Int_Eng,                                    &
     &          TIMSIM%DTnext,                                                 &
     &          MatID                                                          &
     &          )
          ENDIF
!!
!! Add force/moment to existing external forces at the nodal point.
!!
          FORCE(NP)%Xext = FORCE(NP)%Xext - Rx*DAMPER_BC(NVC)%RES%Force
          FORCE(NP)%Yext = FORCE(NP)%Yext - Ry*DAMPER_BC(NVC)%RES%Force
          FORCE(NP)%Zext = FORCE(NP)%Zext - Rz*DAMPER_BC(NVC)%RES%Force
!!
!! Check critical time step.
!!
!!
!! End of node set do-while
!!
        ENDDO
!!
!! End of damper BC do-loop.
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE IMPOSE_NONREFLECTING_BC
!!
!! Copyright (c) by KEY Associates, 9-DEC-1991 20:53:07
!!
!! Purpose: Generate viscous surface tractions in "equilibrium" with the motion
!! reaching the boundary:
!!
      USE shared_common_data
      USE nonreflecting_bc_
      USE nrbc_data_
      USE segment_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER         :: SegID
      REAL(KIND(0D0)) :: Xel(4),Yel(4),Zel(4)

      REAL(KIND(0D0)), PARAMETER :: ONE3RD = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: ONE4TH = (1.0D+0 / 4.0D+0)
!!
!! NOTE: The default value of the constant PARAMVALUE%NRBC_Q is used to make
!! the viscous surface tractions five percent greater.
!!
!! Loop over nonreflecting boundary condition data.
!!
      DO N = 1,NUMND
        SegID = NRBC_DATA(N)%SegID
!!
!! Distinguish between quadrilateral and triangular boundary segments.
!!
        IF (SEGMENT(SegID)%PAR%IX(4) .GT. 0) THEN
          KQT = 4
          Aweight = ONE4TH
        ELSE
          KQT = 3
          Aweight = ONE3RD
        ENDIF
!!
!! Gather boundary segment nodal point coordinates.
!!
        DO k = 1,KQT
          IX = SEGMENT(SegID)%PAR%IX(k)
          Xel(k) = MOTION(IX)%Px + MOTION(IX)%Ux
          Yel(k) = MOTION(IX)%Py + MOTION(IX)%Uy
          Zel(k) = MOTION(IX)%Pz + MOTION(IX)%Uz
        ENDDO
!!
        IF (KQT .EQ. 4) THEN
!!
!! Define Q normal to the loaded boundary segment by taking the cross product
!! of the vectors across the segment diagonals. The magnitude of Q is twice the
!! area of the quadrilateral.
!!
          Qx = (Yel(3)-Yel(1))*(Zel(4)-Zel(2))                                 &
     &         - (Yel(4)-Yel(2))*(Zel(3)-Zel(1))
          Qy = (Zel(3)-Zel(1))*(Xel(4)-Xel(2))                                 &
     &         - (Zel(4)-Zel(2))*(Xel(3)-Xel(1))
          Qz = (Xel(3)-Xel(1))*(Yel(4)-Yel(2))                                 &
     &         - (Xel(4)-Xel(2))*(Yel(3)-Yel(1))
!!
        ELSE
!!
!! Define Q normal to the loaded boundary segment by taking the cross product
!! of the vectors along the segment edges. The magnitude of Q is twice the
!! area of the triangle.
!!
          Qx = (Yel(2)-Yel(1))*(Zel(3)-Zel(1))                                 &
     &         - (Yel(3)-Yel(1))*(Zel(2)-Zel(1))
          Qy = (Zel(2)-Zel(1))*(Xel(3)-Xel(1))                                 &
     &         - (Zel(3)-Zel(1))*(Xel(2)-Xel(1))
          Qz = (Xel(2)-Xel(1))*(Yel(3)-Yel(1))                                 &
     &         - (Xel(3)-Xel(1))*(Yel(2)-Yel(1))
!!
        ENDIF
!!
        Area = 0.5D+0 * SQRT (Qx*Qx + Qy*Qy + Qz*Qz)
!!
!! Normalize the vector Q.
!!
        Qx = Qx * (ONE/(Area+Area))
        Qy = Qy * (ONE/(Area+Area))
        Qz = Qz * (ONE/(Area+Area))
!!
!! Compute the viscous surface tractions in "equilibrium" with the motion
!! reaching the boundary:
!!
!! Fnormal  = rho * longitudinal_wave_speed * Vnormal  * Area * Aweight
!! Ftangent = rho *  transverse_wave_speed  * Vtangent * Area * Aweight
!!
!! Store values of nodal loads in global arrays.
!!
        QLA = (PARAMVALUE%NRBC_Q * Aweight * Area) * NRBC_DATA(N)%RCL
        QTA = (PARAMVALUE%NRBC_Q * Aweight * Area) * NRBC_DATA(N)%RCS
        DO k = 1,KQT
          IX = SEGMENT(SegID)%PAR%IX(k)
          Vn = Qx*MOTION(IX)%Vx + Qy*MOTION(IX)%Vy + Qz*MOTION(IX)%Vz
          Fn = QLA * Vn
          Fx = QTA * (MOTION(IX)%Vx - Qx * Vn)
          Fy = QTA * (MOTION(IX)%Vy - Qy * Vn)
          Fz = QTA * (MOTION(IX)%Vz - Qz * Vn)
          FORCE(IX)%Xext = FORCE(IX)%Xext - (Fx + Qx*Fn)
          FORCE(IX)%Yext = FORCE(IX)%Yext - (Fy + Qy*Fn)
          FORCE(IX)%Zext = FORCE(IX)%Zext - (Fz + Qz*Fn)
        ENDDO
!!
!! End of do-loop over nonreflecting boundary condition data
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE SAVE_COMPLIANCE_FOR_NRBC (MEL,Itype)
!!
!! Copyright (c) by KEY Associates, 15-DEC-1991 15:37:33
!!
!! Purpose: Save element-specific density times the longitudinal and transverse
!! wave speeds for later use in computing nonreflecting boundary viscous loads.
!!
      USE shared_common_data
      USE nrbc_data_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: MEL    ! I/- Current element number (internal ID)
      INTEGER, INTENT(IN) :: Itype  ! I/- Element type 0/1/2=hexah/penta/tetra
!!
!! Local variables.
      LOGICAL :: FOUND   ! -/- Flag for match between element and NRBC segment
!!
!! Test current element for nonrelecting boundary condition.
!!
      N = 0
      FOUND = .FALSE.
      DO WHILE (.NOT.FOUND .AND. N .LT. NUMND)
        N = N + 1
        FOUND = (MEL.EQ.NRBC_DATA(N)%MEL .AND.                                 &
     &    NRBC_DATA(N)%Type.EQ.Itype)
      ENDDO
!!
      IF (FOUND) THEN
!!
!! Save the density times the longitudinal and transverse wave speeds for use
!! later in computing viscous surface tractions.
!!
!!        RCL = rho * longitudinal_wave_speed
!!        RCS = rho * transverse_wave_speed
!!
        NRBC_DATA(N)%RCL = SQRT (SOUND_SPEED%Density * SOUND_SPEED%RCL2)
        NRBC_DATA(N)%RCS = SQRT (SOUND_SPEED%Density * SOUND_SPEED%RCS2)
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION NEXT_NP_ID (SetID,N)
!!
!! Copyright (c) by KEY Associates; 28-MAR-1991 21:47:37
!! Copyright (c) by KEY Associates; 25-DEC-1993 16:36:11.20
!!
!! Purpose: Return .TRUE. and next nodal point ID, N. Terminate with .FALSE.
!! Note: This module must be entered with N = 0 the first time in order to
!! be initialized properly. The module cannot process two overlapped sets.
!! M is an internal counter that allows N to be modified by the calling
!! module without disrupting the "indexing" provided by this function and
!! the do-while-enddo operation.
!!
      USE shared_common_data
      USE node_set_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)    :: SetID
      INTEGER, INTENT(INOUT) :: N
!!
!! Local variables.
      INTEGER, SAVE :: M,MODE,Istart,Iend
!!
!! Use "MODE" to minimize IF-argument evaluations for large node sets.
!!  MODE = 0, Effectively a do-loop from 1 to NUMNP in increments of 1
!!       = 1, Effectively a do-loop stepping through the ID's in NNPSETS
!!       = 2, Converts a negative set ID into a single nodal point ID
!!
      IF (N .EQ. 0) THEN
        M = 0
        IF (SetID .LT. 0) THEN
          MODE = 2
        ELSE IF (SetID .EQ. 0) THEN
          MODE = 0
        ELSE IF (INDEX(NODE_SET(SetID)%Flag,'ALL') .NE. 0) THEN
          MODE = 0
        ELSE
          MODE = 1
          Istart = NODE_SET(SetID)%Istart
          Iend   = NODE_SET(SetID)%Iend
        ENDIF
      ENDIF
!!
!! Process "longest" mode first (MODE equal to 0 loops over all nodes.)
!!
      NEXT_NP_ID = .TRUE.
      IF (MODE .EQ. 0) THEN
        IF (M .LT. NUMNP) THEN
          M = M + 1
          N = M
        ELSE
          N = 0
          NEXT_NP_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 1) THEN
        IF (M .EQ. 0) THEN
          M = Istart
          N = NNPSETS(M)
        ELSE IF (M .LT. Iend) THEN
          M = M + 1
          N = NNPSETS(M)
        ELSE
          N = 0
          NEXT_NP_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 2) THEN
        IF (M .EQ. 0) THEN
          M = 1
          N = ABS(SetID)
        ELSE
          N = 0
          NEXT_NP_ID = .FALSE.
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION NEXT_EL_ID (SetID,N)
!!
!! Copyright (c) by KEY Associates; 28-MAR-1991 21:47:37
!! Copyright (c) by KEY Associates; 25-DEC-1993 16:47:39.58
!!
!! Purpose: Return .TRUE. and next element ID, N. Terminate with .FALSE.
!! Note: This module must be entered with N = 0 the first time in order to
!! be initialized properly. The module cannot process two overlapped sets.
!! M is an internal counter that allows N to be modified by the calling
!! module without disrupting the "indexing" provided by this function and
!! the do-while-enddo operation.
!!
      USE shared_common_data
      USE element_set_
      USE enumerated_sets_, ONLY: NELSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)    :: SetID
      INTEGER, INTENT(INOUT) :: N
!!
!! Local variables.
      INTEGER, SAVE :: M,MODE,Istart,Iend
!!
!! Use "MODE" to minimize IF-argument evaluations for large element sets.
!!  MODE = 0, Effectively a do-loop from 1 to NUMEL in increments of 1
!!       = 1, Effectively a do-loop stepping through the ID's in NELSETS
!!       = 2, Converts a negative set ID into a single nodal point ID
!!
      IF (N .EQ. 0) THEN
        M = 0
        IF (SetID .LT. 0) THEN
          MODE = 2
        ELSE IF (SetID .EQ. 0) THEN
          MODE = 0
        ELSE IF (INDEX(ELEMENT_SET(SetID)%Flag,'ALL') .NE. 0) THEN
          MODE = 0
        ELSE
          MODE = 1
          Istart = ELEMENT_SET(SetID)%Istart
          Iend   = ELEMENT_SET(SetID)%Iend
        ENDIF
      ENDIF
!!
!! Process "longest" mode first (MODE equal to 0 loops over all elements.)
!!
      NEXT_EL_ID = .TRUE.
      IF (MODE .EQ. 0) THEN
        IF (M .LT. NUMEL) THEN
          M = M + 1
          N = M
        ELSE
          N = 0
          NEXT_EL_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 1) THEN
        IF (M .EQ. 0) THEN
          M = Istart
          N = NELSETS(M)
        ELSE IF (M .LT. Iend) THEN
          M = M + 1
          N = NELSETS(M)
        ELSE
          N = 0
          NEXT_EL_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 2) THEN
        IF (M .EQ. 0) THEN
          M = 1
          N = ABS(SetID)
        ELSE
          N = 0
          NEXT_EL_ID = .FALSE.
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION NEXT_SEG_ID (SetID,N)
!!
!! Copyright (c) by KEY Associates; 28-MAR-1991 21:47:37
!! Copyright (c) by KEY Associates; 25-DEC-1993 16:47:16.67
!!
!! Purpose: Return .TRUE. and next segment ID, N. Terminate with .FALSE.
!! Note: This module must be entered with N = 0 the first time in order to
!! be initialized properly. The module cannot process two overlapped sets.
!! M is an internal counter that allows N to be modified by the calling
!! module without disrupting the "indexing" provided by this function and
!! the do-while-enddo operation.
!!
      USE shared_common_data
      USE segment_set_
      USE enumerated_sets_, ONLY: NSGSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)    :: SetID
      INTEGER, INTENT(INOUT) :: N
!!
!! Local variables.
      INTEGER, SAVE :: M,MODE,Istart,Iend
!!
!! Use "MODE" to minimize IF-argument evaluations for large segment sets.
!!  MODE = 0, Effectively a do-loop from 1 to NUMSG in increments of 1
!!       = 1, Effectively a do-loop stepping through the ID's in NSGSETS
!!       = 2, Converts a negative set ID into a single nodal point ID
!!
      IF (N .EQ. 0) THEN
        M = 0
        IF (SetID .LT. 0) THEN
          MODE = 2
        ELSE IF (SetID .EQ. 0) THEN
          MODE = 0
        ELSE IF (INDEX(SEGMENT_SET(SetID)%Flag,'ALL') .NE. 0) THEN
          MODE = 0
        ELSE
          MODE = 1
          Istart = SEGMENT_SET(SetID)%Istart
          Iend   = SEGMENT_SET(SetID)%Iend
        ENDIF
      ENDIF
!!
!! Process "longest" mode first (MODE equal to 0 loops over all segments.)
!!
      NEXT_SEG_ID = .TRUE.
      IF (MODE .EQ. 0) THEN
        IF (M .LT. NUMSG) THEN
          M = M + 1
          N = M
        ELSE
          N = 0
          NEXT_SEG_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 1) THEN
        IF (M .EQ. 0) THEN
          M = Istart
          N = NSGSETS(M)
        ELSE IF (M .LT. Iend) THEN
          M = M + 1
          N = NSGSETS(M)
        ELSE
          N = 0
          NEXT_SEG_ID = .FALSE.
        ENDIF
      ELSE IF (MODE .EQ. 2) THEN
        IF (M .EQ. 0) THEN
          M = 1
          N = ABS(SetID)
        ELSE
          N = 0
          NEXT_SEG_ID = .FALSE.
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION Ipts_PLATT (NEL)
!!
!! Copyright (c) by KEY Associates, 24-APR-1992 13:50:19.89
!!
!! Purpose: Return number of integration stations through the thickness of the
!! shell element.
!!
      USE platt_
      USE section_2d_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, INTENT(IN) :: NEL
!!
      IF (SECTION_2D(PLATT(NEL)%PAR%SecID)%Ipts .NE. 0) THEN
        Ipts = SECTION_2D(PLATT(NEL)%PAR%SecID)%Ipts
      ELSE
        Ipts = TABULATED_FUNCTION(SECTION_2D(                                  &
     &    PLATT(NEL)%PAR%SecID)%Irule)%Number_of_pairs
      ENDIF
!!
      Ipts_PLATT = Ipts
!!
      RETURN
      END
!!_
      INTEGER FUNCTION Ipts_PLATQ (NEL)
!!
!! Copyright (c) by KEY Associates, 24-APR-1992 13:50:19.89
!!
!! Purpose: Return number of integration stations through the thickness of the
!! shell element.
!!
      USE platq_
      USE section_2d_
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, INTENT(IN) :: NEL
!!
      IF (SECTION_2D(PLATQ(NEL)%PAR%SecID)%Ipts .NE. 0) THEN
        Ipts = SECTION_2D(PLATQ(NEL)%PAR%SecID)%Ipts
      ELSE
        Ipts = TABULATED_FUNCTION(SECTION_2D(                                  &
     &    PLATQ(NEL)%PAR%SecID)%Irule)%Number_of_pairs
      ENDIF
!!
      Ipts_PLATQ = Ipts
!!
      RETURN
      END
!!_
      SUBROUTINE PRINTER_OUTPUT
!!
!! Copyright (c) by KEY Associates, 6-MAR-1991 20:27:37
!!
!! Purpose: Provide 132-column, Fortan-carriage-control printer output.
!!
      USE shared_common_data
!!
!! The complete simulation data set.
!!
      USE indx_;          USE node_;          USE tabulated_function_;
      USE beam_;          USE coord_;         USE sliding_interface_;
      USE force_;         USE hexah_;         USE nonreflecting_bc_;
      USE penta_;         USE tetra_;         USE nodal_point_mass_;
      USE lsold_;         USE membt_;         USE constrained_node_;
      USE membq_;         USE truss_;         USE displacement_bc_;
      USE platt_;         USE platq_;         USE rigid_body_mass_;
      USE motion_;        USE stress_;        USE enumerated_sets_;
      USE spring_;        USE damper_;        USE contact_surface_;
      USE segment_;       USE tied_bc_;       USE state_variables_;
      USE results_;       USE gauge1d_;       USE rigid_wall_bc_;
      USE gauge2d_;       USE gauge3d_;       USE contact_node_;
      USE node_set_;      USE force_bc_;      USE sliding_node_;
      USE material_;      USE layering_;      USE segment_set_;
      USE massprop_;      USE spring_bc_;     USE element_set_;
      USE damper_bc_;     USE spot_weld_;     USE periodic_bc_;
      USE qa_record_;     USE nrbc_data_;     USE pressure_bc_;
      USE plate_pair_;    USE section_2d_;    USE section_1d_;
      USE rigid_body_;    USE body_force_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER ::                                                               &
     &          SetID,                                                         &
     &          PAGE_NUMBER = 0        ! Local counter for output pages
      LOGICAL                                                                  &
     &          NEXT_NP_ID,                                                    &
     &          NEXT_EL_ID,                                                    &
     &          ROTATIONAL_DATA

      LOGICAL, SAVE :: FIRST = .TRUE.

      CHARACTER(14), PARAMETER :: &
     &  ORIGIN(3) = (/'Center of Mass','Global Origin ','User Input Pt.'/)
!!
!! Open computed results output file.
!!
!SPEC_CPU2000     OPEN
!SPEC_CPU2000    &          (
!SPEC_CPU2000    &          UNIT   =  IO_UNIT%LCRO,
!SPEC_CPU2000    &          FILE   = 'fmacro',
!SPEC_CPU2000    &          STATUS = 'UNKNOWN',
!SPEC_CPU2000    &          FORM   = 'FORMATTED'
!SPEC_CPU2000    &          )
!!
!! Put header and summary of time step information at start of printed output.
!!
      IF (FIRST) THEN

        CALL INTERCALATION (IO_UNIT%LCRO)

!SPEC_CPU2000        CALL TOTAL_MASS_REPORT (IO_UNIT%LCRO)

!SPEC_omp2001        CALL TIME_STEP_SUMMARY (IO_UNIT%LCRO)

        FIRST = .FALSE.
      ENDIF
!!
!! New page with header.
!!
      CALL NEW_PAGE (IO_UNIT%LCRO, PAGE_NUMBER)
!!
!! Output crurrent time and step counter.
!!
      WRITE (IO_UNIT%LCRO,100) TIMSIM%Total,TIMSIM%Step
!!
!! Print nodal displacements, velocities, and accelerations.
!!
      IF (PRINT%NODES .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,110)
        DO N = 1,NUMNP
          WRITE (IO_UNIT%LCRO,115)                                             &
     &      NODE(N)%ID,                                                        &
     &      MOTION(N)%Ux,MOTION(N)%Uy,MOTION(N)%Uz,                            &
     &      MOTION(N)%Vx,MOTION(N)%Vy,MOTION(N)%Vz,                            &
     &      MOTION(N)%Ax,MOTION(N)%Ay,MOTION(N)%Az,                            &
     &      NODE(N)%ID
        ENDDO
      ELSE IF (PRINT%NODES .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,110)
        SetID = PRINT%NODES
        N = 0
        DO WHILE (NEXT_NP_ID(SetID,N))
!SPEC_CPU2000          WRITE (IO_UNIT%LCRO,115)                                             &
!SPEC_CPU2000     &      NODE(N)%ID,                                                        &
!SPEC_CPU2000     &      MOTION(N)%Ux,MOTION(N)%Uy,MOTION(N)%Uz,                            &
!SPEC_CPU2000     &      MOTION(N)%Vx,MOTION(N)%Vy,MOTION(N)%Vz,                            &
!SPEC_CPU2000     &      MOTION(N)%Ax,MOTION(N)%Ay,MOTION(N)%Az,                            &
!SPEC_CPU2000     &      NODE(N)%ID
          WRITE (IO_UNIT%LCRO,115)                                             &
     &      NODE(N)%ID,                                                        &
     &      MOTION(N)%Ux,MOTION(N)%Uy,MOTION(N)%Uz,                            &
     &      0.0,0.0,0.0,                                                       &
     &      0.0,0.0,0.0,                                                       &
     &      NODE(N)%ID
        ENDDO
      ENDIF
!!
!! Print 8-node hexahedron element stresses.
!!
      IF (NUMHX .GT. 0 .AND. PRINT%HEXAH .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,120)
        DO N = 1,NUMHX
          MatID = HEXAH(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      HEXAH(N)%PAR%EleID,                                                &
     &      (HEXAH(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      HEXAH(N)%RES%Int_Eng,                                              &
     &      HEXAH(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMHX .GT. 0 .AND. PRINT%HEXAH .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,120)
        SetID = PRINT%HEXAH
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = HEXAH(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      HEXAH(N)%PAR%EleID,                                                &
     &      (HEXAH(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      HEXAH(N)%RES%Int_Eng,                                              &
     &      HEXAH(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 6-node pentahedron element stresses.
!!
      IF (NUMPX .GT. 0 .AND. PRINT%PENTA .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,121)
        DO N = 1,NUMPX
          MatID = PENTA(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      PENTA(N)%PAR%EleID,                                                &
     &      (PENTA(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      PENTA(N)%RES%Int_Eng,                                              &
     &      PENTA(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMPX .GT. 0 .AND. PRINT%PENTA .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,121)
        SetID = PRINT%PENTA
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = PENTA(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      PENTA(N)%PAR%EleID,                                                &
     &      (PENTA(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      PENTA(N)%RES%Int_Eng,                                              &
     &      PENTA(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 4-node tetrahedron element stresses.
!!
      IF (NUMTX .GT. 0 .AND. PRINT%TETRA .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,122)
        DO N = 1,NUMTX
          MatID = TETRA(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      TETRA(N)%PAR%EleID,                                                &
     &      (TETRA(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      TETRA(N)%RES%Int_Eng,                                              &
     &      TETRA(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMTX .GT. 0 .AND. PRINT%TETRA .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,122)
        SetID = PRINT%TETRA
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = TETRA(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,125)                                             &
     &      TETRA(N)%PAR%EleID,                                                &
     &      (TETRA(N)%RES%Stress(i),i=1,6),                                    &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      TETRA(N)%RES%Int_Eng,                                              &
     &      TETRA(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 3-node membrane element stresses.
!!
      IF (NUMM3 .GT. 0 .AND. PRINT%MEMBT .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,130)
        DO N = 1,NUMM3
          MatID = MEMBT(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,135)                                             &
     &      MEMBT(N)%PAR%EleID,                                                &
     &      MEMBT(N)%RES%Stress(1),                                            &
     &      MEMBT(N)%RES%Stress(2),                                            &
     &      0.0,                                                               &
     &      MEMBT(N)%RES%Stress(3),                                            &
     &      0.0,                                                               &
     &      0.0,                                                               &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      MEMBT(N)%RES%Int_Eng,                                              &
     &      MEMBT(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMM3 .GT. 0 .AND. PRINT%MEMBT .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,130)
        SetID = PRINT%MEMBT
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = MEMBT(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,135)                                             &
     &      MEMBT(N)%PAR%EleID,                                                &
     &      MEMBT(N)%RES%Stress(1),                                            &
     &      MEMBT(N)%RES%Stress(2),                                            &
     &      0.0,                                                               &
     &      MEMBT(N)%RES%Stress(3),                                            &
     &      0.0,                                                               &
     &      0.0,                                                               &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      MEMBT(N)%RES%Int_Eng,                                              &
     &      MEMBT(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 4-node membrane element stresses.
!!
      IF (NUMM4 .GT. 0 .AND. PRINT%MEMBQ .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,130)
        DO N = 1,NUMM4
          MatID = MEMBQ(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,135)                                             &
     &      MEMBQ(N)%PAR%EleID,                                                &
     &      MEMBQ(N)%RES%Stress(1),                                            &
     &      MEMBQ(N)%RES%Stress(2),                                            &
     &      0.0,                                                               &
     &      MEMBQ(N)%RES%Stress(3),                                            &
     &      0.0,                                                               &
     &      0.0,                                                               &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      MEMBQ(N)%RES%Int_Eng,                                              &
     &      MEMBQ(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMM4 .GT. 0 .AND. PRINT%MEMBQ .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,130)
        SetID = PRINT%MEMBQ
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = MEMBQ(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,135)                                             &
     &      MEMBQ(N)%PAR%EleID,                                                &
     &      MEMBQ(N)%RES%Stress(1),                                            &
     &      MEMBQ(N)%RES%Stress(2),                                            &
     &      0.0,                                                               &
     &      MEMBQ(N)%RES%Stress(3),                                            &
     &      0.0,                                                               &
     &      0.0,                                                               &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      MEMBQ(N)%RES%Int_Eng,                                              &
     &      MEMBQ(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 2-node truss element stresses.
!!
      IF (NUMTR .GT. 0 .AND. PRINT%TRUSS .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,140)
        DO N = 1,NUMTR
          MatID = TRUSS(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,145)                                             &
     &      TRUSS(N)%PAR%EleID,                                                &
     &      TRUSS(N)%RES%Stress,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      TRUSS(N)%RES%Int_Eng,                                              &
     &      TRUSS(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMTR .GT. 0 .AND. PRINT%TRUSS .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,140)
        SetID = PRINT%TRUSS
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = TRUSS(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,145)                                             &
     &      TRUSS(N)%PAR%EleID,                                                &
     &      TRUSS(N)%RES%Stress,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      TRUSS(N)%RES%Int_Eng,                                              &
     &      TRUSS(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 4-node plate element stresses; bottom 1 to top N.
!!
      IF (NUMP4 .GT. 0 .AND. PRINT%PLATQ .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,150)
        DO N = 1,NUMP4
          NPL = N
          Ist = PLATQ(N)%PAR%Ist
          MatID = PLATQ(N)%PAR%MatID
          Ipts = Ipts_PLATQ(NPL)
          DO i = 1,Ipts
            m = Ist + (i-1)
            WRITE (IO_UNIT%LCRO,155)                                           &
     &        PLATQ(N)%PAR%EleID,                                              &
     &        i,                                                               &
     &        Stress(1,m),Stress(2,m),0.0,Stress(4,m),0.0,0.0,                 &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        PLATQ(N)%RES%Int_Eng,                                            &
     &        PLATQ(N)%PAR%EleID
          ENDDO
        ENDDO
      ELSE IF (NUMP4 .GT. 0 .AND. PRINT%PLATQ .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,150)
        SetID = PRINT%PLATQ
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          Ist = PLATQ(N)%PAR%Ist
          MatID = PLATQ(N)%PAR%MatID
          Ipts = Ipts_PLATQ(N)
          DO i = 1,Ipts
            m = Ist + (i-1)
            WRITE (IO_UNIT%LCRO,155)                                           &
     &        PLATQ(N)%PAR%EleID,                                              &
     &        i,                                                               &
     &        Stress(1,m),Stress(2,m),0.0,Stress(4,m),0.0,0.0,                 &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        PLATQ(N)%RES%Int_Eng,                                            &
     &        PLATQ(N)%PAR%EleID
          ENDDO
        ENDDO
      ENDIF
!!
!! Print 3-node plate element stresses; bottom 1 to top N.
!!
      IF (NUMP3 .GT. 0 .AND. PRINT%PLATT .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,150)
        DO N = 1,NUMP3
          NPL = N
          Ist = PLATT(N)%PAR%Ist
          MatID = PLATT(N)%PAR%MatID
          Ipts = Ipts_PLATT(NPL)
          DO i = 1,Ipts
            m = Ist + (i-1)
            WRITE (IO_UNIT%LCRO,155)                                           &
     &        PLATT(N)%PAR%EleID,                                              &
     &        i,                                                               &
     &        Stress(1,m),Stress(2,m),0.0,Stress(4,m),0.0,0.0,                 &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        PLATT(N)%RES%Int_Eng,                                            &
     &        PLATT(N)%PAR%EleID
          ENDDO
        ENDDO
      ELSE IF (NUMP3 .GT. 0 .AND. PRINT%PLATT .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,150)
        SetID = PRINT%PLATT
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          Ist = PLATT(N)%PAR%Ist
          MatID = PLATT(N)%PAR%MatID
          Ipts = Ipts_PLATT(N)
          DO i = 1,Ipts
            m = Ist + (i-1)
            WRITE (IO_UNIT%LCRO,155)                                           &
     &        PLATT(N)%PAR%EleID,                                              &
     &        i,                                                               &
     &        Stress(1,m),Stress(2,m),0.0,Stress(4,m),0.0,0.0,                 &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        PLATT(N)%RES%Int_Eng,                                            &
     &        PLATT(N)%PAR%EleID
          ENDDO
        ENDDO
      ENDIF
!!
!! Print 3-node plate element stresses; bottom 1 to top N.
!!
      IF (NUMBM .GT. 0 .AND. PRINT%BEAMS .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,160)
        DO N = 1,NUMBM
          MatID = BEAM(N)%PAR%MatID
          IF (MATERIAL(MatID)%Type .EQ. 50) THEN
            Ipts = 8
          ELSE
            Ipts = 16
          ENDIF
          DO i = 1,Ipts
            WRITE (IO_UNIT%LCRO,165)                                           &
     &        BEAM(N)%PAR%EleID,                                               &
     &        i,                                                               &
     &        BEAM(N)%RES%Axial(i),0.0,0.0,                                    &
     &        BEAM(N)%RES%Shear(i),BEAM(N)%RES%Trs,BEAM(N)%RES%Trt,            &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        BEAM(N)%RES%Int_Eng,                                             &
     &        BEAM(N)%PAR%EleID
          ENDDO
        ENDDO
      ELSE IF (NUMBM .GT. 0 .AND. PRINT%BEAMS .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,160)
        SetID = PRINT%BEAMS
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = BEAM(N)%PAR%MatID
          IF (MATERIAL(MatID)%Type .EQ. 50) THEN
            Ipts = 8
          ELSE
            Ipts = 16
          ENDIF
          DO i = 1,Ipts
            WRITE (IO_UNIT%LCRO,165)                                           &
     &        BEAM(N)%PAR%EleID,                                               &
     &        i,                                                               &
     &        BEAM(N)%RES%Axial(i),0.0,0.0,                                    &
     &        BEAM(N)%RES%Shear(i),BEAM(N)%RES%Trs,BEAM(N)%RES%Trt,            &
     &        MATERIAL(MatID)%MatID,                                           &
     &        MATERIAL(MatID)%Type,                                            &
     &        BEAM(N)%RES%Int_Eng,                                             &
     &        BEAM(N)%PAR%EleID
          ENDDO
        ENDDO
      ENDIF
!!
!! Print 2-node spring element stresses.
!!
      IF (NUMSP .GT. 0 .AND. PRINT%SPRING .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,170)
        DO N = 1,NUMSP
          MatID = SPRING(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,175)                                             &
     &      SPRING(N)%PAR%EleID,                                               &
     &      SPRING(N)%RES%Force,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      SPRING(N)%RES%Int_Eng,                                             &
     &      SPRING(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMSP .GT. 0 .AND. PRINT%SPRING .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,170)
        SetID = PRINT%SPRING
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = SPRING(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,175)                                             &
     &      SPRING(N)%PAR%EleID,                                               &
     &      SPRING(N)%RES%Force,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      SPRING(N)%RES%Int_Eng,                                             &
     &      SPRING(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print 2-node damper element stresses.
!!
      IF (NUMDM .GT. 0 .AND. PRINT%DAMPER .LT. 0) THEN
        WRITE (IO_UNIT%LCRO,180)
        DO N = 1,NUMDM
          MatID = DAMPER(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,185)                                             &
     &      DAMPER(N)%PAR%EleID,                                               &
     &      DAMPER(N)%RES%Force,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      DAMPER(N)%RES%Int_Eng,                                             &
     &      DAMPER(N)%PAR%EleID
        ENDDO
      ELSE IF (NUMDM .GT. 0 .AND. PRINT%DAMPER .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,180)
        SetID = PRINT%DAMPER
        N = 0
        DO WHILE (NEXT_EL_ID(SetID,N))
          MatID = DAMPER(N)%PAR%MatID
          WRITE (IO_UNIT%LCRO,185)                                             &
     &      DAMPER(N)%PAR%EleID,                                               &
     &      DAMPER(N)%RES%Force,                                               &
     &      (0.0,i=1,5),                                                       &
     &      MATERIAL(MatID)%MatID,                                             &
     &      MATERIAL(MatID)%Type,                                              &
     &      DAMPER(N)%RES%Int_Eng,                                             &
     &      DAMPER(N)%PAR%EleID
        ENDDO
      ENDIF
!!
!! Print mass properties.
!!
      IF (NUMMP .GT. 0) THEN
        WRITE (IO_UNIT%LCRO,190)
        DO N = 1,NUMMP
          WRITE (IO_UNIT%LCRO,195)                                             &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%Mass,                                              &
     &          MASSPROP(N)%Xcm,                                               &
     &          MASSPROP(N)%Ycm,                                               &
     &          MASSPROP(N)%Zcm,                                               &
     &          MASSPROP(N)%Vxcm,                                              &
     &          MASSPROP(N)%Vycm,                                              &
     &          MASSPROP(N)%Vzcm,                                              &
     &          MASSPROP(N)%MPID
        ENDDO
        WRITE (IO_UNIT%LCRO,200)
        DO N = 1,NUMMP
          WRITE (IO_UNIT%LCRO,205)                                             &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%Xmv,                                               &
     &          MASSPROP(N)%Ymv,                                               &
     &          MASSPROP(N)%Zmv,                                               &
     &          MASSPROP(N)%KE,                                                &
     &          MASSPROP(N)%MPID
        ENDDO
        ROTATIONAL_DATA = .FALSE.
        DO N = 1,NUMMP
          ROTATIONAL_DATA = ROTATIONAL_DATA .OR. (MASSPROP(N)%Irot.GT.0)
        ENDDO
        IF (ROTATIONAL_DATA) THEN
          WRITE (IO_UNIT%LCRO,210)
          DO N = 1,NUMMP
            WRITE (IO_UNIT%LCRO,215)                                           &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%Xnert,                                             &
     &          MASSPROP(N)%Ynert,                                             &
     &          MASSPROP(N)%Znert,                                             &
     &          MASSPROP(N)%Vxnert,                                            &
     &          MASSPROP(N)%Vynert,                                            &
     &          MASSPROP(N)%Vznert,                                            &
     &          ORIGIN(MASSPROP(N)%Irot),                                      &
     &          MASSPROP(N)%MPID
          ENDDO
          WRITE (IO_UNIT%LCRO,220)
          DO N = 1,NUMMP
            WRITE (IO_UNIT%LCRO,225)                                           &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%B(1),                                              &
     &          MASSPROP(N)%B(2),                                              &
     &          MASSPROP(N)%B(3),                                              &
     &          MASSPROP(N)%B(4),                                              &
     &          MASSPROP(N)%B(5),                                              &
     &          MASSPROP(N)%B(6),                                              &
     &          ORIGIN(MASSPROP(N)%Irot),                                      &
     &          MASSPROP(N)%MPID
          ENDDO
          WRITE (IO_UNIT%LCRO,230)
          DO N = 1,NUMMP
            WRITE (IO_UNIT%LCRO,235)                                           &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%Omega,                                             &
     &          MASSPROP(N)%Ax,                                                &
     &          MASSPROP(N)%Ay,                                                &
     &          MASSPROP(N)%Az,                                                &
     &          ORIGIN(MASSPROP(N)%Irot),                                      &
     &          MASSPROP(N)%MPID
          ENDDO
          WRITE (IO_UNIT%LCRO,240)
          DO N = 1,NUMMP
            WRITE (IO_UNIT%LCRO,245)                                           &
     &          MASSPROP(N)%MPID,                                              &
     &          MASSPROP(N)%Oxmv,                                              &
     &          MASSPROP(N)%Oymv,                                              &
     &          MASSPROP(N)%Ozmv,                                              &
     &          ORIGIN(MASSPROP(N)%Irot),                                      &
     &          MASSPROP(N)%MPID
          ENDDO
        ENDIF
      ENDIF
!!
!SPEC_CPU2000      CLOSE (UNIT=IO_UNIT%LCRO, STATUS='KEEP')
!!
      RETURN
!!
 100    FORMAT                                                                 &
     &    (                                                                    &
     &    '    Time T =',1PE12.5,'   Time Step =',I8                           &
     &    )
 110    FORMAT                                                                 &
     &    (                                                                    &
     &    '0','   NP    X-DISPLACE   Y-DISPLACE   Z-DISPLACE   ',              &
     &    'X-VELOCITY   Y-VELOCITY   Z-VELOCITY  X-ACCELERATE ',               &
     &    'Y-ACCELERATE Z-ACCELERATE     NP'                                   &
     &    )
 120    FORMAT                                                                 &
     &    (                                                                    &
     &    //5X,'Hex EL',5X,'Sig_X',9X,'Sig_Y',9X,'Sig_Z',8X,                   &
     &    'Sig_XY',8X,'Sig_XZ',8X,'Sig_YZ',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 121    FORMAT                                                                 &
     &    (                                                                    &
     &    //5X,'Pen EL',5X,'Sig_X',9X,'Sig_Y',9X,'Sig_Z',8X,                   &
     &    'Sig_XY',8X,'Sig_XZ',8X,'Sig_YZ',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 122    FORMAT                                                                 &
     &    (                                                                    &
     &    //5X,'Tet EL',5X,'Sig_X',9X,'Sig_Y',9X,'Sig_Z',8X,                   &
     &    'Sig_XY',8X,'Sig_XZ',8X,'Sig_YZ',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 130    FORMAT                                                                 &
     &    (                                                                    &
     &    //5X,'Mbr EL',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',8X,                   &
     &    'Sig_12',8X,'Sig_13',8X,'Sig_23',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 140    FORMAT                                                                 &
     &    (                                                                    &
     &    //3X,'Truss EL',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',8X,                 &
     &    'Sig_12',8X,'Sig_13',8X,'Sig_23',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 150    FORMAT                                                                 &
     &    (                                                                    &
     &    //3X,'Shl EL/j',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',8X,                 &
     &    'Sig_12',8X,'Sig_13',8X,'Sig_23',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 160    FORMAT                                                                 &
     &    (                                                                    &
     &    //4X,'Beam/ j',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',9X,                  &
     &    'Shear',8X,'Tau_XY',8X,'Tau_XZ',10X,'Mat Type  Eng-Den',             &
     &    5X,'EL'                                                              &
     &    )
 170    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Spring EL',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',8X,                &
     &    'Sig_12',8X,'Sig_13',8X,'Sig_23',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 180    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Damper EL',5X,'Sig_1',9X,'Sig_2',9X,'Sig_3',8X,                &
     &    'Sig_12',8X,'Sig_13',8X,'Sig_23',10X,'Mat Type  Eng-Den',            &
     &    5X,'EL'                                                              &
     &    )
 190    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',5X,'Mass',10X,'X_cm',10X,'Y_cm',10X,                &
     &    'Z_cm',9X,'Vx_cm',9X,'Vy_cm',9X,'Vz_cm',10X,'Mass Prop'              &
     &    )
 200    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',5X,'X_Mom',9X,'Y_Mom',9X,'Z_Mom',5X,                &
     &    'Kinetic Eng',4X,'Mass Prop'                                         &
     &    )
 210    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',5X,'X-Ref',9X,'Y-Ref',9X,'Z-Ref',8X,                &
     &    'Vx-Ref',8X,'Vy-Ref',8X,'Vz-Ref',5X,'Origin f/ Rot',                 &
     &    '     Mass Prop'                                                     &
     &    )
 220    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',6X,'Ixx',11X,'Iyy',11X,'Izz',11X,                   &
     &    'Ixy',11X,'Ixz',11X,'Iyz',5X,'Origin f/ Rot',                        &
     &    '      Mass Prop'                                                    &
     &    )
 230    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',4X,'Ang Vel',9X,'Ax ',11X,'Ay ',11X,                &
     &    'Az ',5X,'Origin f/ Rot','      Mass Prop'                           &
     &    )
 240    FORMAT                                                                 &
     &    (                                                                    &
     &    //2X,'Mass Prop',5X,'Ox-Mom',8X,'Oy-Mom',8X,'Oz-Mom',                &
     &    5X,'Origin f/ Rot','      Mass Prop'                                 &
     &    )
 115    FORMAT                                                                 &
     &    (                                                                    &
     &    1X,I8,1X,9(1PE13.4),3X,I8                                            &
     &    )
 125    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),6X,I4,I8,1PE11.2,I8                                    &
     &    )
 135    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),6X,I4,I8,1PE11.2,I8                                    &
     &    )
 145    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),6X,I4,I8,1PE11.2,I8                                    &
     &    )
 155    FORMAT                                                                 &
     &    (                                                                    &
     &    I8,'/',I1,6(1PE14.4),6X,I4,I8,1PE11.2,I8                             &
     &    )
 165    FORMAT                                                                 &
     &    (                                                                    &
     &    I8,'/',I2,6(1PE14.4),6X,I4,I8,1PE11.2,I8                             &
     &    )
 175    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),6X,I4,I8,1PE11.2,I8                                    &
     &    )
 185    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),6X,I4,I8,1PE11.2,I8                                    &
     &    )
 195    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,7(1PE14.4),6X,I8                                                  &
     &    )
 205    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,4(1PE14.4),6X,I8                                                  &
     &    )
 215    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),3X,A,4X,I8                                             &
     &    )
 225    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,6(1PE14.4),2X,A,4X,I8                                             &
     &    )
 235    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,4(1PE14.4),2X,A,4X,I8                                             &
     &    )
 245    FORMAT                                                                 &
     &    (                                                                    &
     &    I9,3(1PE14.4),4X,A,4X,I8                                             &
     &    )
!!
      END
