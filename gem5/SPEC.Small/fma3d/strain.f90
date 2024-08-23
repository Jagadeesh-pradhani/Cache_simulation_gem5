      SUBROUTINE STRAIN_GAUGE_INITIALIZATION
!!
!! Copyright (c) Copyright by KEY Associates, 3-DEC-1991 20:46:20
!!
      USE shared_common_data
      USE gauge1d_
      USE gauge2d_
      USE gauge3d_
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
      USE node_
      USE section_1d_
      USE section_2d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          EleID,                                                         &
     &          TB_Flag,                                                       &
     &          HPT_Flag,                                                      &
     &          MPQT_Flag
!!
!! Set error count to zero.
!!
      ERROR%COUNT = 0
!!
!! INITIALIZE 1-D STRAIN GAUGES
!! Loop over one-dimensional strain gauges.
!!
      DO N = 1,NUMG1
!!
!! Bring truss/beam flag into range.
!!
        GAUGE1D(N)%PAR%TBflg = MAX (0, MIN (2,GAUGE1D(N)%PAR%TBflg))
!!
!! Retrieve truss/beam flag.
!!
        TB_Flag = GAUGE1D(N)%PAR%TBflg
!!
!! If the strain gauge is defined by a truss element or beam element,
!! gather gauge nodal ID's from the appropriate element.
!!
        EleID = GAUGE1D(N)%PAR%EleID
        IF (TB_Flag .EQ. 0) THEN
          NIX = 0
          DO i = 1,4
            IF (GAUGE1D(N)%PAR%IX(i) .GT. 0) NIX = NIX + 1
          ENDDO
          GAUGE1D(N)%PAR%NUMIX = NIX
!!
!! Convert translational ID's to rotational ID's.
!!
          IF (NIX .EQ. 4) THEN
            IF (NODE(GAUGE1D(N)%PAR%IX(3))%IRT .GT. 0) THEN
              GAUGE1D(N)%PAR%IX(3) = NODE(GAUGE1D(N)%PAR%IX(3))%IRT
            ELSE
              ERROR%COUNT = ERROR%COUNT + 1
              WRITE (MSG1,'(I8)') GAUGE1D(N)%PAR%GauID
              WRITE (MSG2,'(I8)') NODE(GAUGE1D(N)%PAR%IX(3))%ID
              CALL USER_MESSAGE                                                &
     &  (                                                                      &
     &  MSGL//'WARN'//                                                         &
     &  MSGL//'STRAIN_GAUGE_INITIALIZATION.002.01'//                           &
     &  MSGL//'GAUGE1D Input Record ID:'//MSG1//                               &
     &  MSGL//'References NPT ID:'//MSG2//' Without Rotational DOF''s.'        &
     &  )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
            IF (NODE(GAUGE1D(N)%PAR%IX(4))%IRT .GT. 0) THEN
              GAUGE1D(N)%PAR%IX(4) = NODE(GAUGE1D(N)%PAR%IX(4))%IRT
            ELSE
              ERROR%COUNT = ERROR%COUNT + 1
              WRITE (MSG1,'(I8)') GAUGE1D(N)%PAR%GauID
              WRITE (MSG2,'(I8)') NODE(GAUGE1D(N)%PAR%IX(4))%ID
              CALL USER_MESSAGE                                                &
     &  (                                                                      &
     &  MSGL//'WARN'//                                                         &
     &  MSGL//'STRAIN_GAUGE_INITIALIZATION.002.02'//                           &
     &  MSGL//'GAUGE1D Input Record ID:'//MSG1//                               &
     &  MSGL//'References NPT ID:'//MSG2//' Without Rotational DOF''s.'        &
     &  )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ELSE IF (TB_Flag .EQ. 1) THEN
          DO i = 1,2
            GAUGE1D(N)%PAR%IX(i) = TRUSS(EleID)%PAR%IX(i)
          ENDDO
          GAUGE1D(N)%PAR%NUMIX = 2
        ELSE IF (TB_Flag .EQ. 2) THEN
          DO i = 1,4
            GAUGE1D(N)%PAR%IX(i) = BEAM(EleID)%PAR%IX(i)
          ENDDO
          GAUGE1D(N)%PAR%NUMIX = 4
        ENDIF
!!
!! Check to make certain a sufficient number of nodal points has been defined.
!!
        IF (GAUGE1D(N)%PAR%NUMIX .LT. 2) THEN
          WRITE (MSG1,'(I8)') GAUGE1D(N)%PAR%GauID
          CALL USER_MESSAGE                                                    &
     &    (                                                                    &
     &    MSGL//'WARN'//                                                       &
     &    MSGL//'STRAIN_GAUGE_INITIALIZATION.001.01'//                         &
     &    MSGL//'1-D Strain Gauge (GAUGE1D) Input Record ID:'//MSG1//          &
     &    MSGL//'Is Defined With Less Than 2 Nodal Points.'                    &
     &    )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Initialize strain to zero.
!!
        GAUGE1D(N)%RES%Strain = 0.0
        GAUGE1D(N)%RES%Torsion = 0.0
      ENDDO
!!
!! INITIALIZE 2-D STRAIN GAUGES
!! Loop over two-dimensional strain gauges.
!!
      DO N = 1,NUMG2
!!
!! Bring membrane/plate flag into range.
!!
        GAUGE2D(N)%PAR%MPQTflg = MAX (0, MIN (4,GAUGE2D(N)%PAR%MPQTflg))
!!
!! Retrieve membrane/plate flag.
!!
        MPQT_Flag = GAUGE2D(N)%PAR%MPQTflg
!!
!! If the strain gauge is defined by a membrane element or plate element,
!! gather gauge nodal ID's from the appropriate element.
!!
        EleID = GAUGE2D(N)%PAR%EleID
        IF (MPQT_Flag .EQ. 0) THEN
          NIX = 0
          DO i = 1,4
            IF (GAUGE2D(N)%PAR%IX(i) .GT. 0) NIX = NIX + 1
          ENDDO
          GAUGE2D(N)%PAR%NUMIX = NIX
        ELSE IF (MPQT_Flag .EQ. 1) THEN
          DO i = 1,3
            GAUGE2D(N)%PAR%IX(i) = MEMBT(EleID)%PAR%IX(i)
          ENDDO
          GAUGE2D(N)%PAR%NUMIX = 3
          GAUGE2D(N)%PAR%Thickness = 0.0
          GAUGE2D(N)%PAR%Refloc    = 0.0
        ELSE IF (MPQT_Flag .EQ. 2) THEN
          DO i = 1,4
            GAUGE2D(N)%PAR%IX(i) = MEMBQ(EleID)%PAR%IX(i)
          ENDDO
          GAUGE2D(N)%PAR%NUMIX = 4
          GAUGE2D(N)%PAR%Thickness = 0.0
          GAUGE2D(N)%PAR%Refloc    = 0.0
        ELSE IF (MPQT_Flag .EQ. 3) THEN
          DO i = 1,6
            GAUGE2D(N)%PAR%IX(i) = PLATT(EleID)%PAR%IX(i)
          ENDDO
          GAUGE2D(N)%PAR%NUMIX = 6
          GAUGE2D(N)%PAR%Thickness =                                           &
     &      SECTION_2D(PLATT(EleID)%PAR%SecID)%Thickness
          GAUGE2D(N)%PAR%Refloc    =                                           &
     &      SECTION_2D(PLATT(EleID)%PAR%SecID)%Refloc
        ELSE IF (MPQT_Flag .EQ. 4) THEN
          DO i = 1,8
            GAUGE2D(N)%PAR%IX(i) = PLATQ(EleID)%PAR%IX(i)
          ENDDO
          GAUGE2D(N)%PAR%NUMIX = 8
          GAUGE2D(N)%PAR%Thickness =                                           &
     &      SECTION_2D(PLATQ(EleID)%PAR%SecID)%Thickness
          GAUGE2D(N)%PAR%Refloc    =                                           &
     &      SECTION_2D(PLATQ(EleID)%PAR%SecID)%Refloc
        ENDIF
!!
!! Check to make certain a sufficient number of nodal points has been defined.
!!
        IF (GAUGE2D(N)%PAR%NUMIX .LT. 3) THEN
          WRITE (MSG1,'(I8)') GAUGE2D(N)%PAR%GauID
          CALL USER_MESSAGE                                                    &
     &      (                                                                  &
     &      MSGL//'WARN'//                                                     &
     &      MSGL//'STRAIN_GAUGE_INITIALIZATION.001.02'//                       &
     &      MSGL//'2-D Strain Gauge (GAUGE2D) Input Record ID:'//MSG1//        &
     &      MSGL//'Is Defined With Less Than 3 Nodal Points.'                  &
     &      )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Initialize strain to zero.
!!
        DO i = 1,3
          GAUGE2D(N)%RES%Strain(i) = 0.0
        ENDDO
      ENDDO
!!
!! INITIALIZE 3-D STRAIN GAUGES
!! Loop over three-dimensional strain gauges.
!!
      DO N = 1,NUMG3
!!
!! Bring hexahedron/pentahedron/tetrahedron flag into range.
!!
        GAUGE3D(N)%PAR%HPTflg = MAX (0, MIN (3,GAUGE3D(N)%PAR%HPTflg))
!!
!! Retrieve hexahedron/pentahedron/tetrahedron flag.
!!
        HPT_Flag = GAUGE3D(N)%PAR%HPTflg
!!
!! If the strain gauge is defined by a solid element gather gauge nodal ID's
!! from the appropriate element.
!!
        EleID = GAUGE3D(N)%PAR%EleID
        IF (HPT_Flag .EQ. 0) THEN
          NIX = 0
          DO i = 1,8
            IF (GAUGE3D(N)%PAR%IX(i) .GT. 0) NIX = NIX + 1
          ENDDO
          GAUGE3D(N)%PAR%NUMIX = NIX
        ELSE IF (HPT_Flag .EQ. 1) THEN
          DO i = 1,8
            GAUGE3D(N)%PAR%IX(i) = HEXAH(EleID)%PAR%IX(i)
          ENDDO
          GAUGE3D(N)%PAR%NUMIX = 8
        ELSE IF (HPT_Flag .EQ. 2) THEN
          DO i = 1,6
            GAUGE3D(N)%PAR%IX(i) = PENTA(EleID)%PAR%IX(i)
          ENDDO
          GAUGE3D(N)%PAR%NUMIX = 6
        ELSE IF (HPT_Flag .EQ. 3) THEN
          DO i = 1,4
            GAUGE3D(N)%PAR%IX(i) = TETRA(EleID)%PAR%IX(i)
          ENDDO
          GAUGE3D(N)%PAR%NUMIX = 4
        ENDIF
!!
!! Check to make certain a sufficient number of nodal points has been defined.
!!
        IF (GAUGE3D(N)%PAR%NUMIX .LT. 4) THEN
          WRITE (MSG1,'(I8)') GAUGE3D(N)%PAR%GauID
          CALL USER_MESSAGE                                                    &
     &      (                                                                  &
     &      MSGL//'WARN'//                                                     &
     &      MSGL//'STRAIN_GAUGE_INITIALIZATION.001.03'//                       &
     &      MSGL//'3-D Strain Gauge (GAUGE3D) Input Record ID:'//MSG1//        &
     &      MSGL//'Is Defined With Less Than 4 Nodal Points.'                  &
     &      )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Initialize strain to zero.
!!
        DO i = 1,6
          GAUGE3D(N)%RES%Strain(i) = 0.0
        ENDDO
      ENDDO
!!
!! Check to see if any errors were found.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'STRAIN_GAUGE_INITIALIZATION.001.04'//                   &
     &          MSGL//'Total Number Of Strain Gauge Errors:'//MSG1//           &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE STRAIN_GAUGE_INTEGRATION
!!
!! (c) Copyright by KEY Associates, 3-DEC-1991 20:46:20
!!
      USE shared_common_data
      USE gauge1d_
      USE gauge2d_
      USE gauge3d_
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
      USE node_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          TB_Flag,                                                       &
     &          HPT_Flag,                                                      &
     &          MPQT_Flag
!!
!! PROCESS 1-D STRAIN GAUGES
!!
      DO N = 1,NUMG1
!!
!! Gather gauge nodal motion.
!!
        DO i = 1,GAUGE1D(N)%PAR%NUMIX
          EMOTION(i) = MOTION(GAUGE1D(N)%PAR%IX(i))
        ENDDO
!!
!! For subcycling, scale nodal positions to current time.
!!
        IF (CONTROL%SUBCYC .NE. 0) THEN
          DO i = 1,GAUGE1D(N)%PAR%NUMIX
            QA = NODE(GAUGE1D(N)%PAR%IX(i))%Time - TIMSIM%Total
            EMOTION(i)%Ux =                                                    &
     &        EMOTION(i)%Ux - QA * EMOTION(i)%Vx
            EMOTION(i)%Uy =                                                    &
     &        EMOTION(i)%Uy - QA * EMOTION(i)%Vy
            EMOTION(i)%Uz =                                                    &
     &        EMOTION(i)%Uz - QA * EMOTION(i)%Vz
          ENDDO
        ENDIF
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Count gauge execution.
!!
        COUNTER%GAUGE1D = COUNTER%GAUGE1D + 1
!!
!! Update gauge strain.
!!
        CALL GAUGE1D_STRAIN_INTEGRATION ( NEL )
!!
!! End of one-dimensional gauge do-loop.
!!
      ENDDO
!!
!! Record time spent in integrating 1-D strain gauges.
!!
!SPEC_CPU2000      IF (NUMG1 .GT. 0) CALL TIMER (23)
!!
!! PROCESS 2-D STRAIN GAUGES
!!
      DO N = 1,NUMG2
!!
!! Gather gauge nodal motion.
!!
        DO i = 1,GAUGE2D(N)%PAR%NUMIX
          EMOTION(i) = MOTION(GAUGE2D(N)%PAR%IX(i))
        ENDDO
!!
!! For subcycling, scale nodal positions to current time.
!!
        IF (CONTROL%SUBCYC .NE. 0) THEN
          DO i = 1,GAUGE2D(N)%PAR%NUMIX
            QA = NODE(GAUGE2D(N)%PAR%IX(i))%Time - TIMSIM%Total
            EMOTION(i)%Ux =                                                    &
     &        EMOTION(i)%Ux - QA * EMOTION(i)%Vx
            EMOTION(i)%Uy =                                                    &
     &        EMOTION(i)%Uy - QA * EMOTION(i)%Vy
            EMOTION(i)%Uz =                                                    &
     &        EMOTION(i)%Uz - QA * EMOTION(i)%Vz
          ENDDO
        ENDIF
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Count gauge execution.
!!
        COUNTER%GAUGE2D = COUNTER%GAUGE2D + 1
!!
!! Update gauge strain.
!!
        CALL GAUGE2D_STRAIN_INTEGRATION ( NEL )
!!
!! End of two-dimensional gauge do-loop.
!!
      ENDDO
!!
!! Record time spent in integrating 2-D strain gauges.
!!
!SPEC_CPU2000      IF (NUMG2 .GT. 0) CALL TIMER (24)
!!
!! PROCESS 3-D STRAIN GAUGES
!!
      DO N = 1,NUMG3
!!
!! Gather gauge nodal motion.
!!
        DO i = 1,GAUGE3D(N)%PAR%NUMIX
          EMOTION(i) = MOTION(GAUGE3D(N)%PAR%IX(i))
        ENDDO
!!
!! For subcycling, scale nodal positions to current time.
!!
        IF (CONTROL%SUBCYC .NE. 0) THEN
          DO i = 1,GAUGE3D(N)%PAR%NUMIX
            QA = NODE(GAUGE3D(N)%PAR%IX(i))%Time - TIMSIM%Total
            EMOTION(i)%Ux =                                                    &
     &        EMOTION(i)%Ux - QA * EMOTION(i)%Vx
            EMOTION(i)%Uy =                                                    &
     &        EMOTION(i)%Uy - QA * EMOTION(i)%Vy
            EMOTION(i)%Uz =                                                    &
     &        EMOTION(i)%Uz - QA * EMOTION(i)%Vz
          ENDDO
        ENDIF
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Count gauge execution.
!!
        COUNTER%GAUGE3D = COUNTER%GAUGE3D + 1
!!
!! Update gauge strain.
!!
        CALL GAUGE3D_STRAIN_INTEGRATION ( NEL )
!!
!! End of three-dimensional gauge do-loop.
!!
      ENDDO
!!
!! Record time spent in integrating 3-D strain gauges.
!!
!SPEC_CPU2000      IF (NUMG3 .GT. 0) CALL TIMER (25)
!!
      RETURN
      END
!!_
      SUBROUTINE GAUGE1D_STRAIN_INTEGRATION ( NEL )
!!
!! (c) Copyright by KEY Associates, 3-DEC-1991 21:25:36
!!
!! Purpose: Calculate the gradient of the axial velocity and axial twist.
!!
      USE shared_common_data
      USE gauge1d_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          Vx(2),Vy(2),Vz(2),X(2),Y(2),Z(2),Wx(2),Wy(2),Wz(2)
!!
!! Mid-interval position of nodal points, and current translational velocities.
!!
      DT = 0.5 * TIMSIM%DTnext
      DO i = 1,2
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! Define a basis vector R aligned with the gauge. Save length calculation
!! until gradient operation is performed. A square root operation is saved.
!!
      Rx = X(2) - X(1)
      Ry = Y(2) - Y(1)
      Rz = Z(2) - Z(1)
!!
!! Axial stretching using velocity Vx,Vy,Vz transformed to local R coordinate.
!!
      dVr = Rx*(Vx(2)-Vx(1)) + Ry*(Vy(2)-Vy(1)) + Rz*(Vz(2)-Vz(1))
      Drr = dVr / (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! Update gauge axial strain.
!!
      GAUGE1D(NEL)%RES%Strain =                                                &
     &  GAUGE1D(NEL)%RES%Strain + TIMSIM%DTnext * Drr
!!
!! Check to see if torsional strain is required.
!!
      IF (GAUGE1D(NEL)%PAR%TBflg .EQ. 2) THEN
!!
!! Retrieve current angular velocity.
!!
        Wx(1) = EMOTION(3)%Vx
        Wy(1) = EMOTION(3)%Vy
        Wz(1) = EMOTION(3)%Vz
        Wx(2) = EMOTION(4)%Vx
        Wy(2) = EMOTION(4)%Vy
        Wz(2) = EMOTION(4)%Vz
!!
!! Axial twisting using angular velocity Wx,Wy,Wz transformed to local R
!! coordinate.
!!
        dWr = Rx*(Wx(2)-Wx(1)) + Ry*(Wy(2)-Wy(1)) + Rz*(Wz(2)-Wz(1))
        Wrr = dWr / (Rx*Rx + Ry*Ry + Rz*Rz)
!!
!! Update gauge torsional strain.
!!
        GAUGE1D(NEL)%RES%Torsion =                                             &
     &    GAUGE1D(NEL)%RES%Torsion + TIMSIM%DTnext * Wrr
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE GAUGE2D_STRAIN_INTEGRATION ( NEL )
!!
!! (c) Copyright by KEY Associates, 3-DEC-1991 21:25:36
!!
!! Purpose: Calculate the gradient of the axial velocity and axial twist.
!!
      USE shared_common_data
      USE gauge2d_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8),Vx(8),Vy(8),Vz(8),                              &
     &          R(8),S(8),Vr(4),Vs(4),Wr(4),Ws(4),                             &
     &          Br(4),Bs(4),Strain(3)
!!
!! Mid-interval position of nodal points, and current velocities.
!!
      DT = 0.5 * TIMSIM%DTnext
      DO i = 1,GAUGE2D(NEL)%PAR%NUMIX
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
      DO i = 1,3
        Strain(i) = GAUGE2D(NEL)%RES%Strain(i)
      ENDDO
!!
!! Select appropriate geometry for strain rate computation.
!! GAUGE2D(NEL)%PAR%NUMIX = 3, triangular    membrane
!! GAUGE2D(NEL)%PAR%NUMIX = 4, quadrilateral membrane
!! GAUGE2D(NEL)%PAR%NUMIX = 6, triangular    plate
!! GAUGE2D(NEL)%PAR%NUMIX = 8, quadrilateral plate
!!
!! Triangular membrane and plate case.
!!
      IF (GAUGE2D(NEL)%PAR%NUMIX .EQ. 3 .OR.                                   &
     &    GAUGE2D(NEL)%PAR%NUMIX .EQ. 6) THEN
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
          WRITE (MSG1,'(I8)') GAUGE2D(NEL)%PAR%GauID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'GAUGE2D_STRAIN_INTEGRATION.001.01'//                    &
     &          MSGL//'Gauge Geometry Has An Undefined Normal.'//              &
     &          MSGL//'Gauge ID:'//MSG1                                        &
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
          Wr(i) = Rx*Vx(i+3) + Ry*Vy(i+3) + Rz*Vz(i+3)
          Ws(i) = Sx*Vx(i+3) + Sy*Vy(i+3) + Sz*Vz(i+3)
        ENDDO
!!
!! Gradient operator.
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
        Ain = 2.0 / (R12*S13 - R13*S12)
!!
!! Construct reference surface velocity gradients.
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
!! Check for shell surface gauge.
!!
        IF (GAUGE2D(NEL)%PAR%NUMIX .EQ. 6) THEN
!!
!! Check surface designator and compute distance from reference surface
!! to gauge location.
!!
          IF (GAUGE2D(NEL)%PAR%GauLoc .LT. 0) THEN
            H = (GAUGE2D(NEL)%PAR%Refloc - 0.5) *                              &
     &        GAUGE2D(NEL)%PAR%Thickness
          ELSE IF (GAUGE2D(NEL)%PAR%GauLoc .GT. 0) THEN
            H = (0.5 - GAUGE2D(NEL)%PAR%Refloc) *                              &
     &        GAUGE2D(NEL)%PAR%Thickness
          ELSE
            H = 0.0
          ENDIF
!!
          IF (H .NE. 0.0) THEN
!!
!! Angular velocity gradients.
!!
            Wrr = 0.0
            Wsr = 0.0
            Wrs = 0.0
            Wss = 0.0
            DO i = 1,3
              Wrr = Wrr + Wr(i)*Br(i)
              Wsr = Wsr + Ws(i)*Br(i)
              Wrs = Wrs + Wr(i)*Bs(i)
              Wss = Wss + Ws(i)*Bs(i)
            ENDDO
            Wrr = Ain * Wrr
            Wsr = Ain * Wsr
            Wrs = Ain * Wrs
            Wss = Ain * Wss
!!
!! Surface velocity gradients.
!!
            Vrr = Vrr + H * Wsr
            Vsr = Vsr - H * Wrr
            Vrs = Vrs + H * Wss
            Vss = Vss - H * Wrs
!!
          ENDIF
        ENDIF
!!
!! Quadrilateral membrane/plate case.
!!
      ELSE
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
          WRITE (MSG1,'(I8)') GAUGE2D(NEL)%PAR%GauID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'GAUGE2D_STRAIN_INTEGRATION.001.02'//                    &
     &          MSGL//'Gauge Geometry Has An Undefined Normal.'//              &
     &          MSGL//'Gauge ID:'//MSG1                                        &
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
          Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
          Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
          Wr(i) = Rx*Vx(i+4) + Ry*Vy(i+4) + Rz*Vz(i+4)
          Ws(i) = Sx*Vx(i+4) + Sy*Vy(i+4) + Sz*Vz(i+4)
        ENDDO
!!
!! Membrane gradient operator.
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
        Ain = 2.0 / (R13*S24 - S13*R24)
!!
!! Construct velocity gradients.
!!
        Vrr = 0.0
        Vsr = 0.0
        Vrs = 0.0
        Vss = 0.0
        DO i = 1,4
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
!! Check for shell surface gauge.
!!
        IF (GAUGE2D(NEL)%PAR%NUMIX .EQ. 8) THEN
!!
!! Check surface designator and compute distance from reference surface
!! to gauge location.
!!
          IF (GAUGE2D(NEL)%PAR%GauLoc .LT. 0) THEN
            H = (GAUGE2D(NEL)%PAR%Refloc - 0.5) *                              &
     &        GAUGE2D(NEL)%PAR%Thickness
          ELSE IF (GAUGE2D(NEL)%PAR%GauLoc .GT. 0) THEN
            H = (0.5 - GAUGE2D(NEL)%PAR%Refloc) *                              &
     &        GAUGE2D(NEL)%PAR%Thickness
          ELSE
            H = 0.0
          ENDIF
!!
          IF (H .NE. 0.0) THEN
!!
!! Angular velocity gradients.
!!
            Wrr = 0.0
            Wsr = 0.0
            Wrs = 0.0
            Wss = 0.0
            DO i = 1,4
              Wrr = Wrr + Wr(i)*Br(i)
              Wsr = Wsr + Ws(i)*Br(i)
              Wrs = Wrs + Wr(i)*Bs(i)
              Wss = Wss + Ws(i)*Bs(i)
            ENDDO
            Wrr = Ain * Wrr
            Wsr = Ain * Wsr
            Wrs = Ain * Wrs
            Wss = Ain * Wss
!!
!! Surface velocity gradients.
!!
            Vrr = Vrr + H * Wsr
            Vsr = Vsr - H * Wrr
            Vrs = Vrs + H * Wss
            Vss = Vss - H * Wrs
!!
          ENDIF
        ENDIF
      ENDIF
!!
!! Stretching components in local coordinates.
!!
      Drr = Vrr
      Dss = Vss
      Drs = 0.5 * (Vrs + Vsr)
!!
!! POLAR DECOMPOSITION OF THE DEFORMATION GRADIENT.
!! Compute the polar decomposition of the deformation gradient beteen
!! time (n) and time (n+1). Construct components of the deformation
!! gradient; compute the sum of F11 and F22, and the difference of F12
!! and F21.  If r(R,S,t) and s(R,S,t), then QSF = r,R + s,S  and
!! QDF = r,S - s,R. Hin = 1.0 / (2.0*Area(0)) cancels from the quotient
!! in the expression for the inverse tangent function, ATan.
!!
      DTnext = TIMSIM%DTnext
      dBeta = ATAN (DTnext*(Vrs-Vsr)/(2.0 + DTnext*(Vrr+Vss)))
      GAUGE2D(NEL)%RES%Beta = GAUGE2D(NEL)%RES%Beta + dBeta
!!
!! "De-rotate" strain components to account for the fact that the
!! "co-rotational" axes have been redrawn from the last time step
!! along the same material fiber. The required rotational rate comes
!! from Vsr. (The rotation based on dBeta contains the finite rotations
!! of the element relative to the co-rotational coordinate axes, but
!! does not account for the fact that the co-rotational axes rotated
!! with the element.)
!!
!! Rotate the stress from time n to time n+1. This is a proper orthognal
!! rotation using the rotation obtained from the polar decomposition of
!! the deformation gradient.
!!
      Xi = dBeta + TIMSIM%DTnext * Vsr
      CB = COS (Xi)
      SB = SIN (Xi)
      S2B = SB*SB
      C2B = CB*CB
      SCB = SB*CB
      Q1 = (Strain(3)*SCB) + (Strain(3)*SCB)
      Err = Strain(1)*C2B + Strain(2)*S2B + Q1
      Ess = Strain(2)*C2B + Strain(1)*S2B - Q1
      Ers = Strain(3)*(C2B-S2B) + SCB*(Strain(2)-Strain(1))
!!
!! Return strains to global storage and accumulate strain increment.
!!
      GAUGE2D(NEL)%RES%Strain(1) = Err + TIMSIM%DTnext * Drr
      GAUGE2D(NEL)%RES%Strain(2) = Ess + TIMSIM%DTnext * Dss
      GAUGE2D(NEL)%RES%Strain(3) = Ers + TIMSIM%DTnext * Drs
!!
      RETURN
      END
!!_
      SUBROUTINE GAUGE3D_STRAIN_INTEGRATION ( NEL )
!!
!! (c) Copyright by KEY Associates, 3-DEC-1991 21:25:36
!!
!! Purpose: Calculate the gradient of the axial velocity and axial twist.
!!
      USE shared_common_data
      USE gauge3d_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          I2(6),I3(6),I4(6),I5(6),I6(6),I7(6),I8(6),                     &
     &          J2(8),J3(8),J4(8),J5(8),J6(8),J7(8),J8(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8),Vx(8),Vy(8),Vz(8)
!!
      DATA                                                                     &
     &          I2 /2,3,1,5,5,6/, I3 /3,1,1,7,5,5/, I4 /1,1,2,6,7,5/,          &
     &          I5 /5,6,7,1,2,3/, I6 /6,7,5,1,1,2/, I7 /7,5,5,3,1,1/,          &
     &          I8 /5,5,6,2,3,1/
      DATA                                                                     &
     &          J2 /2,3,4,1,8,5,6,7/, J3 /3,4,1,2,7,8,5,6/,                    &
     &          J4 /4,1,2,3,6,7,8,5/, J5 /5,6,7,8,1,2,3,4/,                    &
     &          J6 /6,7,8,5,4,1,2,3/, J7 /7,8,5,6,3,4,1,2/,                    &
     &          J8 /8,5,6,7,2,3,4,1/
!!
!! Mid-interval position of nodal points, and current translational velocities.
!!
      DT = 0.5 * TIMSIM%DTnext
      DO i = 1,GAUGE3D(NEL)%PAR%NUMIX
        X(i)  = EMOTION(i)%Px + (EMOTION(i)%Ux - DT * EMOTION(i)%Vx)
        Y(i)  = EMOTION(i)%Py + (EMOTION(i)%Uy - DT * EMOTION(i)%Vy)
        Z(i)  = EMOTION(i)%Pz + (EMOTION(i)%Uz - DT * EMOTION(i)%Vz)
        Vx(i) = EMOTION(i)%Vx
        Vy(i) = EMOTION(i)%Vy
        Vz(i) = EMOTION(i)%Vz
      ENDDO
!!
!! Select appropriate geometry for strain rate computation.
!! GAUGE3D(NEL)%PAR%NUMIX = 4, Tetrahedral solid
!! GAUGE3D(NEL)%PAR%NUMIX = 6, Pentahedral solid
!! GAUGE3D(NEL)%PAR%NUMIX = 8, Hexahedral  solid
!!
      IF (GAUGE3D(NEL)%PAR%NUMIX .EQ. 4) THEN
!!
!! Construct useful differences.
!!
        X12 = X(1) - X(2)
        X13 = X(1) - X(3)
        X14 = X(1) - X(4)
        X24 = X(2) - X(4)
        X34 = X(3) - X(4)
        Y12 = Y(1) - Y(2)
        Y13 = Y(1) - Y(3)
        Y14 = Y(1) - Y(4)
        Y24 = Y(2) - Y(4)
        Y34 = Y(3) - Y(4)
        Z12 = Z(1) - Z(2)
        Z13 = Z(1) - Z(3)
        Z14 = Z(1) - Z(4)
        Z24 = Z(2) - Z(4)
        Z34 = Z(3) - Z(4)
!!
!! Gradient operators.
!!
        Bx(1) = 0.16666666667 * (Y34*Z24 - Y24*Z34)
        Bx(2) = 0.16666666667 * (Y13*Z14 - Y14*Z13)
        Bx(3) = 0.16666666667 * (Y14*Z12 - Y12*Z14)
        Bx(4) = 0.16666666667 * (Y12*Z13 - Y13*Z12)
!!
        By(1) = 0.16666666667 * (Z34*X24 - Z24*X34)
        By(2) = 0.16666666667 * (Z13*X14 - Z14*X13)
        By(3) = 0.16666666667 * (Z14*X12 - Z12*X14)
        By(4) = 0.16666666667 * (Z12*X13 - Z13*X12)
!!
        Bz(1) = 0.16666666667 * (X34*Y24 - X24*Y34)
        Bz(2) = 0.16666666667 * (X13*Y14 - X14*Y13)
        Bz(3) = 0.16666666667 * (X14*Y12 - X12*Y14)
        Bz(4) = 0.16666666667 * (X12*Y13 - X13*Y12)
!!
      ELSE IF (GAUGE3D(NEL)%PAR%NUMIX .EQ. 6) THEN
!!
!! Pentahedron gradient operators.
!!
        DO i = 1,6
          Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))             &
     &               +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))           &
     &               +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))           &
     &               +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                             &
     &               +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                             &
     &               +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
        ENDDO
        DO i = 1,6
          By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))             &
     &               +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))           &
     &               +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))           &
     &               +Z(I3(i))*(X(I2(i))-X(I4(i)))                             &
     &               +Z(I6(i))*(X(I5(i))-X(I2(i)))                             &
     &               +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
        ENDDO
        DO i = 1,6
          Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))             &
     &               +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))           &
     &               +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))           &
     &               +X(I3(i))*(Y(I2(i))-Y(I4(i)))                             &
     &               +X(I6(i))*(Y(I5(i))-Y(I2(i)))                             &
     &               +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
        ENDDO
!!
      ELSE IF (GAUGE3D(NEL)%PAR%NUMIX .EQ. 8) THEN
!!
!! Hexahedron gradient operators.
!!
        DO i = 1,8
          Bx(i) = ( Y(J2(i))*(Z(J6(i))-Z(J3(i))+Z(J5(i))-Z(J4(i)))             &
     &               +Y(J4(i))*(Z(J3(i))-Z(J8(i))+Z(J2(i))-Z(J5(i)))           &
     &               +Y(J5(i))*(Z(J8(i))-Z(J6(i))+Z(J4(i))-Z(J2(i)))           &
     &               +Y(J3(i))*(Z(J2(i))-Z(J4(i)))                             &
     &               +Y(J6(i))*(Z(J5(i))-Z(J2(i)))                             &
     &               +Y(J8(i))*(Z(J4(i))-Z(J5(i))) ) * 0.08333333333
        ENDDO
        DO i = 1,8
          By(i) = ( Z(J2(i))*(X(J6(i))-X(J3(i))+X(J5(i))-X(J4(i)))             &
     &               +Z(J4(i))*(X(J3(i))-X(J8(i))+X(J2(i))-X(J5(i)))           &
     &               +Z(J5(i))*(X(J8(i))-X(J6(i))+X(J4(i))-X(J2(i)))           &
     &               +Z(J3(i))*(X(J2(i))-X(J4(i)))                             &
     &               +Z(J6(i))*(X(J5(i))-X(J2(i)))                             &
     &               +Z(J8(i))*(X(J4(i))-X(J5(i))) ) * 0.08333333333
        ENDDO
        DO i = 1,8
          Bz(i) = ( X(J2(i))*(Y(J6(i))-Y(J3(i))+Y(J5(i))-Y(J4(i)))             &
     &               +X(J4(i))*(Y(J3(i))-Y(J8(i))+Y(J2(i))-Y(J5(i)))           &
     &               +X(J5(i))*(Y(J8(i))-Y(J6(i))+Y(J4(i))-Y(J2(i)))           &
     &               +X(J3(i))*(Y(J2(i))-Y(J4(i)))                             &
     &               +X(J6(i))*(Y(J5(i))-Y(J2(i)))                             &
     &               +X(J8(i))*(Y(J4(i))-Y(J5(i))) ) * 0.08333333333
        ENDDO
!!
      ENDIF
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      Volume = 0.0
      DO i = 1,GAUGE3D(NEL)%PAR%NUMIX
        Volume = Volume + Bx(i)*X(i)
      ENDDO
      Vin = 1.0 / Volume
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
      DO i = 1,GAUGE3D(NEL)%PAR%NUMIX
        Vxx = Vxx + Vx(i) * Bx(i)
        Vyx = Vyx + Vy(i) * Bx(i)
        Vzx = Vzx + Vz(i) * Bx(i)
        Vxy = Vxy + Vx(i) * By(i)
        Vyy = Vyy + Vy(i) * By(i)
        Vzy = Vzy + Vz(i) * By(i)
        Vxz = Vxz + Vx(i) * Bz(i)
        Vyz = Vyz + Vy(i) * Bz(i)
        Vzz = Vzz + Vz(i) * Bz(i)
      ENDDO
      Vxx = Vin * Vxx
      Vyx = Vin * Vyx
      Vzx = Vin * Vzx
      Vxy = Vin * Vxy
      Vyy = Vin * Vyy
      Vzy = Vin * Vzy
      Vxz = Vin * Vxz
      Vyz = Vin * Vyz
      Vzz = Vin * Vzz

      Dxx = Vxx
      Dyy = Vyy
      Dzz = Vzz
      Dxy = 0.5 * (Vxy + Vyx)
      Dxz = 0.5 * (Vxz + Vzx)
      Dyz = 0.5 * (Vyz + Vzy)
      Wxy = 0.5 * (Vxy - Vyx)
      Wxz = 0.5 * (Vxz - Vzx)
      Wyz = 0.5 * (Vyz - Vzy)
!!
!! Generate the required rotation operator (and new stretching components
!! if Ipolard equals one or two).
!!
      Igenrot = 1
      Ipolard = CONTROL%POLARD
!!
      CALL SYMMETRIC_TENSOR_ROTATION                                           &
     &          (                                                              &
     &          GAUGE3D(NEL)%RES%Strain,Igenrot,Ipolard,TIMSIM%DTnext          &
     &          )
!!
!! Rotate strain from the configuration at time n to the configuration
!! at time n+1 (Ipolard=0).
!!
      IF (Ipolard .EQ. 0) CALL SYMMETRIC_TENSOR_ROTATION                       &
     &          (                                                              &
     &          GAUGE3D(NEL)%RES%Strain,Igenrot,Ipolard,TIMSIM%DTnext          &
     &          )
!!
!! Accumulate strain increment.
!!
      GAUGE3D(NEL)%RES%Strain(1) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(1) + TIMSIM%DTnext * Dxx
      GAUGE3D(NEL)%RES%Strain(2) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(2) + TIMSIM%DTnext * Dyy
      GAUGE3D(NEL)%RES%Strain(3) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(3) + TIMSIM%DTnext * Dzz
      GAUGE3D(NEL)%RES%Strain(4) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(4) + TIMSIM%DTnext * Dxy
      GAUGE3D(NEL)%RES%Strain(5) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(5) + TIMSIM%DTnext * Dxz
      GAUGE3D(NEL)%RES%Strain(6) =                                             &
     &  GAUGE3D(NEL)%RES%Strain(6) + TIMSIM%DTnext * Dyz
!!
!! If the end-of-the-interval Rashid approximate polar decomposition is
!! being used (Ipolard=1,2), rotate strain from the configuration at
!! time n to the configuration at time n+1. Note: the rotation operator R
!! has been saved from the generation entry (Igenrot=1).
!!
      IF (Ipolard .NE. 0) THEN
        CALL SYMMETRIC_TENSOR_ROTATION                                         &
     &    (GAUGE3D(NEL)%RES%Strain,Igenrot,Ipolard,TIMSIM%DTnext)
      ENDIF
!!
      RETURN
      END
