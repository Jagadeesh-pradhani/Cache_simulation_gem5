      SUBROUTINE LSOLD_INITIALIZATION
!!
!! Copyright (c) by KEY Associates;  4-MAR-1994 18:07:36.96
!!
      USE shared_common_data
      USE lsold_
      USE hexah_, ONLY:LSHEX, hexah_type
      USE membq_, ONLY:LSMBQ, membq_type
      USE layering_
      USE section_2d_
      USE material_
      USE node_
      USE motion_
      USE force_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Local record structures. Note: MXLSL is a PARAMETER in "lsold_"
!!
      TYPE (motion_type), DIMENSION(8)     :: ELEMENT_MOTION
      TYPE (force_type),  DIMENSION(8)     :: ELEMENT_FORCE
      TYPE (motion_type), DIMENSION(8)     :: LAYER_MOTION
      TYPE (node_type),   DIMENSION(8)     :: LAYER_NODE
      TYPE (hexah_type),  DIMENSION(MXLSL) :: LYHEX
      TYPE (membq_type),  DIMENSION(MXLSL) :: LYMBQ
!!
      INTEGER                                                                  &
     &          IX,                                                            &
     &          Isv,                                                           &
     &          SecID,                                                         &
     &          LupID,                                                         &
     &          MatID
      REAL(KIND(0D0))                                                          &
     &          Internal_Energy,                                               &
     &          Membrane_Factor,                                               &
     &          Element_Thickness,                                             &
     &          Membrane_Thickness,                                            &
     &          B1(4),B2(4),T1(4),T2(4),                                       &
     &          X(8),Y(8),Z(8),Ux(8),Uy(8),Uz(8),Vx(8),Vy(8),Vz(8),            &
     &          R(8),S(8),T(8),Ur(8),Us(8),Ut(8),Vr(8),Vs(8),Vt(8)
      LOGICAL                                                                  &
     &          FOUND
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Loop over all layered solid elements.
!!
      DO N = 1,NUMLS
!!
!! Gather element motion; clear element internal force record.
!!
        DO i = 1,8
          ELEMENT_MOTION(i) = MOTION(LSOLD(N)%PAR%IX(i))
          ELEMENT_FORCE(i) = force_type (0,0,0,0,0,0)
        ENDDO
!!
!! Initialize element clock and time step
!!
        LSOLD(N)%RES%Time   = 0.0
        LSOLD(N)%RES%DTnext = 0.0
!!
!! Access element do-loop index for use in subroutine calls.
!!
        NEL = N
!!
!! Obtain element volume and generalized element dimensions for stability
!! calculations.
!!
        CALL LSOLD_VOLUME_OPERATOR ( NEL,ELEMENT_MOTION,DSQ )
!!
!! Save initial element volume for later volume strain calculations.
!!
        LSOLD(N)%PAR%Volume = LSOLD(N)%RES%Volume
!!
!! Distinguish between a rigid body element and a deformable element.
!!
        IF (LSOLD(N)%PAR%ParID .LT. 0) THEN
!!
!! "Critical time step." Used to control integration when only rigid body
!! domains are present.
!!
          LSOLD(N)%RES%H(1,1) = SQRT (DSQ)
          MatID = LAYERING(LSOLD(N)%PAR%LupID)%MatID(1)
          Density = MATERIAL(MatID)%PVAL(1)
          Ymod = MATERIAL(MatID)%PVAL(6)
          LSOLD(N)%RES%DTelt = LSOLD(N)%RES%H(1,1) / SQRT (Ymod/Density)
!!
        ELSE
!!
!! Define the current position of the layered solid nodal points.
!!
          DO i = 1,8
            X(i)  = ELEMENT_MOTION(i)%Px + ELEMENT_MOTION(i)%Ux
            Y(i)  = ELEMENT_MOTION(i)%Py + ELEMENT_MOTION(i)%Uy
            Z(i)  = ELEMENT_MOTION(i)%Pz + ELEMENT_MOTION(i)%Uz
            Ux(i) = ELEMENT_MOTION(i)%Ux
            Uy(i) = ELEMENT_MOTION(i)%Uy
            Uz(i) = ELEMENT_MOTION(i)%Uz
            Vx(i) = ELEMENT_MOTION(i)%Vx
            Vy(i) = ELEMENT_MOTION(i)%Vy
            Vz(i) = ELEMENT_MOTION(i)%Vz
          ENDDO
!!
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element sides. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
          Rx = (X(3)-X(1)) + (X(2)-X(4)) + (X(7)-X(5)) + (X(6)-X(8))
          Ry = (Y(3)-Y(1)) + (Y(2)-Y(4)) + (Y(7)-Y(5)) + (Y(6)-Y(8))
          Rz = (Z(3)-Z(1)) + (Z(2)-Z(4)) + (Z(7)-Z(5)) + (Z(6)-Z(8))
          Sx = (X(3)-X(1)) - (X(2)-X(4)) + (X(7)-X(5)) - (X(6)-X(8))
          Sy = (Y(3)-Y(1)) - (Y(2)-Y(4)) + (Y(7)-Y(5)) - (Y(6)-Y(8))
          Sz = (Z(3)-Z(1)) - (Z(2)-Z(4)) + (Z(7)-Z(5)) - (Z(6)-Z(8))
          Rmag  = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
          Rx = Rx * (ONE / Rmag)
          Ry = Ry * (ONE / Rmag)
          Rz = Rz * (ONE / Rmag)
          Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
          Sx = Sx * (ONE / Smag)
          Sy = Sy * (ONE / Smag)
          Sz = Sz * (ONE / Smag)
!!
!! Define the unit vector T normal to the element.
!!
          Tx = Ry*Sz - Sy*Rz
          Ty = Rz*Sx - Sz*Rx
          Tz = Rx*Sy - Sx*Ry
          Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
          IF (Tmag .EQ. 0.0) THEN
            WRITE (MSG1,'(I8)') LSOLD(N)%PAR%EleID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LSOLD_INITIALIZATION.001.00'//                          &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
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
!! Transform position X,Y,Z, displacement Ux,Uy,Uz, and translational
!! velocity Vx,Vy,Vz to the local R,S,T-coordinate system.
!!
          DO i = 1,8
            R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
            S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
            T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
            Ur(i) = Rx*Ux(i) + Ry*Uy(i) + Rz*Uz(i)
            Us(i) = Sx*Ux(i) + Sy*Uy(i) + Sz*Uz(i)
            Ut(i) = Tx*Ux(i) + Ty*Uy(i) + Tz*Uz(i)
            Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
            Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
            Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
          ENDDO
!!
!! "Convert" global motion to local directions.
!!
          DO i = 1,8
            ELEMENT_MOTION(i)%Px = R(i) - Ur(i)
            ELEMENT_MOTION(i)%Py = S(i) - Us(i)
            ELEMENT_MOTION(i)%Pz = T(i) - Ut(i)
            ELEMENT_MOTION(i)%Ux = Ur(i)
            ELEMENT_MOTION(i)%Uy = Us(i)
            ELEMENT_MOTION(i)%Uz = Ut(i)
            ELEMENT_MOTION(i)%Vx = Vr(i)
            ELEMENT_MOTION(i)%Vy = Vs(i)
            ELEMENT_MOTION(i)%Vz = Vt(i)
          ENDDO
!!
!! Retrieve lay-up pointer LupID and number of layers NUMLY in lay-up.
!!
          LupID = LSOLD(N)%PAR%LupID
          NUMLY = LAYERING(LupID)%Number_of_Layers
!!
!! Normalize element corner thicknesses and store with element. Special
!! attention is paid to membrane thicknesses to preserve their percentage
!! in the proportional factors. The thicknesses of the hexahedral sub-
!! elements are scaled to "fit" the element dimension.
!!
          Membrane_Thickness = 0.0
          DO m = 1,NUMLY
            IF (LAYERING(LupID)%Ltype(m) .NE. 0) THEN
              Membrane_Thickness = Membrane_Thickness +                        &
     &          SECTION_2D(LAYERING(LupID)%Ltype(m))%Thickness
            ENDIF
          ENDDO
!!
!! Loop over all four corners of the layered solid.
!!
          DO i = 1,4
!!
!! Compute element corner thickness. This calculation assumes
!! that the membrane layers are, more or less, perpendicular to
!! the element "vertical" edges.
!!
            DX = (ELEMENT_MOTION(i)%Px - ELEMENT_MOTION(i+4)%Px)
            DY = (ELEMENT_MOTION(i)%Py - ELEMENT_MOTION(i+4)%Py)
            DZ = (ELEMENT_MOTION(i)%Pz - ELEMENT_MOTION(i+4)%Pz)
            Element_Thickness = SQRT (DX*DX + DY*DY + DZ*DZ)

            Membrane_Factor =                                                  &
     &        (Element_Thickness-Membrane_Thickness) / Element_Thickness

            Scale = 0.0
            DO m = 1,NUMLY
              IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
                Scale = Scale + LAYERING(LupID)%H(i,m)
              ENDIF
            ENDDO
            Scale = Scale / Membrane_Factor

            DO m = 1,NUMLY
              IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
                LSOLD(N)%RES%H(i,m) = (ONE/Scale)*LAYERING(LupID)%H(i,m)
              ELSE
                LSOLD(N)%RES%H(i,m) =                                          &
     &            SECTION_2D(LAYERING(LupID)%Ltype(m))%Thickness /             &
     &            Element_Thickness
              ENDIF
            ENDDO

          ENDDO
!!
!! Gather "sub-element" data and store in local LYHEX and LYMBQ records.
!!
          DO m = 1,NUMLY
            IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
              LYHEX(m) = LSHEX(LSOLD(N)%PAR%ID(m))
            ELSE
              LYMBQ(m) = LSMBQ(LSOLD(N)%PAR%ID(m))
            ENDIF
          ENDDO
!!
!! In preparation for looping over all of the layers in a layered solid
!! element, set-up the motion and weightings for transferring data between
!! the layer nodal points and the LSOLD nodal points.
!!
          DO i = 1,4
            LAYER_MOTION(i+4) = ELEMENT_MOTION(i)
            T1(i) = ONE
            T2(i) = 0.0
          ENDDO
          CSQmax = 0.0
!!
!! Loop over all layers in layered solid element processing each layer as
!! a hexahedron element or a membrane element.
!!
          DO M = 1,NUMLY
!!
!! Construct "sub-element." The do-loop starts at the layer next to the "1"
!! surface of the layered solid, IX(1:4), and ends at the layer next to the
!! "2" surface of the layered solid, IX(5:8). The weights T1 and T2 are for
!! use with data associted with the top of the layer, and B1 and B2 are for
!! use with data associated with the bottom of the layer. The "1" in B1 and
!! T1 refer to the "1" surface of the layered solid, IX(1:4). The "2" in B2
!! and T2 refer to the "2" surface of the layered solid, IX(5:8). The weights
!! are used to transfer results between the layer nodal points and the LSOLD
!! nodal points. Note that each corner of the LSOLD hexahedron can have dif-
!! ferent thicknesses associated with the layer, but the same number of layers
!! is expected at each corner. The interpolation at each corner is linear.
!!
            DO i = 1,4
              B1(i) = T1(i)
              B2(i) = T2(i)
              T1(i) = T1(i) - LSOLD(N)%RES%H(i,M)
              T2(i) = ONE - T1(i)
              LAYER_MOTION(i) = LAYER_MOTION(i+4)
              LAYER_MOTION(i+4)%Px = T1(i) * ELEMENT_MOTION(i  )%Px            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Px

              LAYER_MOTION(i+4)%Py = T1(i) * ELEMENT_MOTION(i  )%Py            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Py

              LAYER_MOTION(i+4)%Pz = T1(i) * ELEMENT_MOTION(i  )%Pz            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Pz

              LAYER_MOTION(i+4)%Ux = T1(i) * ELEMENT_MOTION(i  )%Ux            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Ux

              LAYER_MOTION(i+4)%Uy = T1(i) * ELEMENT_MOTION(i  )%Uy            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Uy

              LAYER_MOTION(i+4)%Uz = T1(i) * ELEMENT_MOTION(i  )%Uz            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Uz

              LAYER_MOTION(i+4)%Vx = T1(i) * ELEMENT_MOTION(i  )%Vx            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Vx

              LAYER_MOTION(i+4)%Vy = T1(i) * ELEMENT_MOTION(i  )%Vy            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Vy

              LAYER_MOTION(i+4)%Vz = T1(i) * ELEMENT_MOTION(i  )%Vz            &
     &                             + T2(i) * ELEMENT_MOTION(i+4)%Vz
            ENDDO
!!
!! Distinguish between hexahedron and quadrilateral membrane sub-elements.
!!
            IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
!!
!! Case (a): HEXAHERDAL SUBLAYER
!! Process layer as a hexahedron. Retrieve material ID and state variable
!! starting location.
!!
              MatID = LYHEX(M)%PAR%MatID
              Isv   = LYHEX(M)%PAR%Isv
!!
!! Initialize "sub-element" IX array. (HEXAH routines will be used to build
!! a layer mass matrix after which this routine will distribute the mass to
!! the LSOLD nodal points in the global array NODE.)
!!
              DO i = 1,8
                LYHEX(M)%PAR%IX(i) = i
              ENDDO
!!
!! Gradient operator, stretching, and rotation.
!!
              CALL LYHEX_GRADIENT_OPERATOR (LYHEX(M),LAYER_MOTION)
!!
!! Save initial element volume for later volume strain calculations.
!!
              LYHEX(M)%PAR%Volume = LYHEX(M)%RES%Volume
!!
!! Compute element mass matrix and store in global mass matrix.
!!
              DO i = 1,8
                LAYER_NODE(i)%Mass = 0.0
              ENDDO
              CALL LYHEX_MASS                                                  &
     &            (LYHEX(M),MatID,LAYER_NODE,LAYER_MOTION)
              DO i = 1,4
                IX = LSOLD(N)%PAR%IX(i)
                NODE(IX)%Mass = NODE(IX)%Mass +                                &
     &            B1(i)*LAYER_NODE(i)%Mass + T1(i)*LAYER_NODE(i+4)%Mass
                IX = LSOLD(N)%PAR%IX(i+4)
                NODE(IX)%Mass = NODE(IX)%Mass +                                &
     &            B2(i)*LAYER_NODE(i)%Mass + T2(i)*LAYER_NODE(i+4)%Mass
              ENDDO
!!
!! Find initial sound speeds and stress.
!!
              MTRL_TYPE = MATERIAL(MatID)%Type
              SELECT CASE (MTRL_TYPE)
              CASE (40)
                CALL MATERIAL_40                                               &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
              CASE (41)
                CALL MATERIAL_41                                               &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID,RCL2,RCS2                                        &
     &                  )
              CASE (45)
!!
!! Convert global material direction vectors to local element coordinates.
!!
                CALL GLOBAL_TO_LOCAL                                           &
     &            (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))

                CALL MATERIAL_45                                               &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
              CASE (47)
                CALL MATERIAL_47                                               &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
              CASE DEFAULT
                WRITE (MSG1,'(I8)') LSOLD(N)%PAR%EleID
                WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
                CALL USER_MESSAGE                                              &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LSOLD_INITIALIZATION.001.01'//                          &
     &          MSGL//'LSEL (8-Node Layered Solid) Element ID:'//MSG1//        &
     &          MSGL//'References An Unknown Material Type:'//MSG2//           &
     &          MSGL//'Via a Hexahedron LAYERING Specification.'               &
     &          )
              END SELECT
!!
!! Compute initial stress divergence and time step.
!!
              CALL LYHEX_STRESS_DIVERGENCE (LYHEX(M),MatID)
!!
!! Collect maximum layer tangent sound speed squared. (CSQ = C**2)
!!
              Density = MATERIAL(MatID)%PVAL(1)
              CSQtmp = SOUND_SPEED%RCL2/Density
              CSQmax = MAX (CSQmax,CSQtmp)

            ELSE
!!
!! Case (b): QUADRILATERAL MEMBRANE SUBLAYER
!! Process layer as a quadrilateral. Retrieve material ID and state variable
!! starting location.
!!
              MatID = LYMBQ(M)%PAR%MatID
              Isv   = LYMBQ(M)%PAR%Isv
!!
!! Initialize "sub-element" IX array. (MEMBQ routines will be used to build
!! a layer mass matrix after which this routine will distribute the mass to
!! the LSOLD nodal points in the global array NODE%)
!!
              DO i = 1,4
                LYMBQ(M)%PAR%IX(i) = i
              ENDDO
!!
!! Gradient operator, stretching, and rotation.
!!
              DO i = 1,4
                LAYER_MOTION(i)%Px = 0.5 * ( LAYER_MOTION(i  )%Px +            &
     &                                       LAYER_MOTION(i+4)%Px )
                LAYER_MOTION(i)%Py = 0.5 * ( LAYER_MOTION(i  )%Py +            &
     &                                       LAYER_MOTION(i+4)%Py )
                LAYER_MOTION(i)%Pz = 0.5 * ( LAYER_MOTION(i  )%Pz +            &
     &                                       LAYER_MOTION(i+4)%Pz )
                LAYER_MOTION(i)%Ux = 0.5 * ( LAYER_MOTION(i  )%Ux +            &
     &                                       LAYER_MOTION(i+4)%Ux )
                LAYER_MOTION(i)%Uy = 0.5 * ( LAYER_MOTION(i  )%Uy +            &
     &                                       LAYER_MOTION(i+4)%Uy )
                LAYER_MOTION(i)%Uz = 0.5 * ( LAYER_MOTION(i  )%Uz +            &
     &                                       LAYER_MOTION(i+4)%Uz )
                LAYER_MOTION(i)%Vx = 0.5 * ( LAYER_MOTION(i  )%Vx +            &
     &                                       LAYER_MOTION(i+4)%Vx )
                LAYER_MOTION(i)%Vy = 0.5 * ( LAYER_MOTION(i  )%Vy +            &
     &                                       LAYER_MOTION(i+4)%Vy )
                LAYER_MOTION(i)%Vz = 0.5 * ( LAYER_MOTION(i  )%Vz +            &
     &                                       LAYER_MOTION(i+4)%Vz )
              ENDDO
              CALL LYMBQ_GRADIENT_OPERATOR (LYMBQ(M),LAYER_MOTION)
!!
!! Save initial element volume for later volume strain calculations.
!!
              LYMBQ(M)%PAR%Area = LYMBQ(M)%RES%Area
!!
!! Compute element mass matrix and store in global mass matrix.
!!
              DO i = 1,4
                LAYER_NODE(i)%Mass = 0.0
              ENDDO
              SecID = LAYERING(LupID)%Ltype(M)
              CALL LYMBQ_MASS                                                  &
     &            (LYMBQ(M),SecID,MatID,LAYER_NODE,LAYER_MOTION)
              DO i = 1,4
                IX = LSOLD(N)%PAR%IX(i)
                NODE(IX)%Mass =                                                &
     &            NODE(IX)%Mass + 0.5*(B1(i)+T1(i))*LAYER_NODE(i)%Mass
                IX = LSOLD(N)%PAR%IX(i+4)
                NODE(IX)%Mass =                                                &
     &            NODE(IX)%Mass + 0.5*(B2(i)+T2(i))*LAYER_NODE(i)%Mass
              ENDDO
!!
!! Find initial sound speeds and stress.
!!
              MTRL_TYPE = MATERIAL(MatID)%Type
              SELECT CASE (MTRL_TYPE)
              CASE (20)
                CALL MATERIAL_20                                               &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
              CASE (21)
               CALL MATERIAL_21                                                &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
              CASE (22)
!!
!! Convert global material direction vectors to local element coordinates.
!!
                CALL GLOBAL_TO_LOCAL                                           &
     &            (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))

                CALL MATERIAL_22                                               &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
              CASE (25)
!!
!! Convert global material direction vectors to local element coordinates.
!!
                CALL GLOBAL_TO_LOCAL                                           &
     &            (Rx,Ry,Rz,Sx,Sy,Sz,STATE_VARIABLES(Isv))

                CALL MATERIAL_25                                               &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
              CASE (27)
                CALL MATERIAL_27                                               &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
              CASE DEFAULT
                WRITE (MSG1,'(I8)') LSOLD(N)%PAR%EleID
                WRITE (MSG2,'(I8)') MATERIAL(MatID)%Type
                CALL USER_MESSAGE                                              &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LSOLD_INITIALIZATION.001.02'//                          &
     &          MSGL//'LSEL (8-Node Layered Solid) Element ID:'//MSG1//        &
     &          MSGL//'References An Unknown Material Type:'//MSG2//           &
     &          MSGL//'Via a Hexahedron LAYERING Specification.'               &
     &          )
              END SELECT
!!
!! Compute initial stress divergence and time step.
!!
              SecID = LAYERING(LupID)%Ltype(M)
              CALL LYMBQ_STRESS_DIVERGENCE (LYMBQ(M),SecID,MatID)
!!
!! Collect maximum layer tangent sound speed squared. (CSQ = C**2)
!!
              Density = MATERIAL(MatID)%PVAL(1)
              CSQtmp = SOUND_SPEED%RCL2/Density
              CSQmax = MAX (CSQmax,CSQtmp)
!!
!! End of layer type if-test.
!!
            ENDIF
!!
!! End of layer do-loop.
!!
          ENDDO
!!
!! Update global "sub-element" storage.
!!
          DO M = 1,NUMLY
            IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
              LSHEX(LSOLD(N)%PAR%ID(M)) = LYHEX(M)
            ELSE
              LSMBQ(LSOLD(N)%PAR%ID(M)) = LYMBQ(M)
            ENDIF
          ENDDO
!!
!! In preparation for looping over all of the layers in a layered solid
!! element, set-up the weightings for transferring data between the layer
!! nodal points and the LSOLD nodal points.
!!
          T1(1) = 0.0
          T1(2) = 0.0
          T1(3) = 0.0
          T1(4) = 0.0
!!
!! Accumulate nodal forces. Distribute layer internal forces to LSOLD nodal poin
!! based on force and moment equilibrium.
!!
          DO i = 1,4
            DO M = 1,NUMLY
              B1(i) = T1(i)
              T1(i) = T1(i) + LSOLD(N)%RES%H(i,M)
!!
!! Distinguish between a hexahedron layer and membrane layer.
!!
              IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
!!
!! Forces at the bottom of the layer distributed to the LSOLD nodes.
!!
                Fx = B1(i) * LYHEX(M)%RES%Xint(i)
                Fy = B1(i) * LYHEX(M)%RES%Yint(i)
                Fz = B1(i) * LYHEX(M)%RES%Zint(i)
                ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                Fx = LYHEX(M)%RES%Xint(i) - Fx
                Fy = LYHEX(M)%RES%Yint(i) - Fy
                Fz = LYHEX(M)%RES%Zint(i) - Fz
                ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz
!!
!! Forces at the top of the layer distributed to the LSOLD nodes.
!!
                Fx = T1(i) * LYHEX(M)%RES%Xint(i+4)
                Fy = T1(i) * LYHEX(M)%RES%Yint(i+4)
                Fz = T1(i) * LYHEX(M)%RES%Zint(i+4)
                ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                Fx = LYHEX(M)%RES%Xint(i+4) - Fx
                Fy = LYHEX(M)%RES%Yint(i+4) - Fy
                Fz = LYHEX(M)%RES%Zint(i+4) - Fz
                ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz

              ELSE
!!
!! Forces in the membrane layer distributed to the LSOLD nodes.
!!
                Fx = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Xint(i)
                Fy = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Yint(i)
                Fz = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Zint(i)
                ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                Fx = LYMBQ(M)%RES%Xint(i) - Fx
                Fy = LYMBQ(M)%RES%Yint(i) - Fy
                Fz = LYMBQ(M)%RES%Zint(i) - Fz
                ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz

              ENDIF
            ENDDO
          ENDDO
!!
!! Transfer element nodal forces to element data structure and at the same
!! time transform the element nodal forces from the local r,s,t-coordinates
!! to the global x,y,z-coordinates.
!!
          DO i = 1,8
            LSOLD(N)%RES%Xint(i) = Rx*ELEMENT_FORCE(i)%Xint                    &
     &                           + Sx*ELEMENT_FORCE(i)%Yint                    &
     &                           + Tx*ELEMENT_FORCE(i)%Zint
            LSOLD(N)%RES%Yint(i) = Ry*ELEMENT_FORCE(i)%Xint                    &
     &                           + Sy*ELEMENT_FORCE(i)%Yint                    &
     &                           + Ty*ELEMENT_FORCE(i)%Zint
            LSOLD(N)%RES%Zint(i) = Rz*ELEMENT_FORCE(i)%Xint                    &
     &                           + Sz*ELEMENT_FORCE(i)%Yint                    &
     &                           + Tz*ELEMENT_FORCE(i)%Zint
          ENDDO
!!
!! Construct element critical time step. (In-plane the layer moduli act
!! in parallel; perpendicular the layer moduli act in series.)
!!
          LSOLD(N)%RES%DTelt = SQRT (DSQ/CSQmax)
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
          IF (NUMNR .NE. 0) THEN
            Itype = -1  !  Layered Solid Hexahedron
!! fix-up sound speed data!!
            CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
          ENDIF
!!
!! Compute internal energy density of layered solid from "sub-elements."
!!
          Internal_Energy = 0.0
          DO M = 1,NUMLY
            IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
              Internal_Energy = Internal_Energy +                              &
     &          LYHEX(M)%RES%Volume * LYHEX(M)%RES%Int_Eng
            ELSE
              Internal_Energy = Internal_Energy + LYMBQ(M)%PAR%Area *          &
     &          SECTION_2D(LAYERING(LupID)%Ltype(m))%Thickness *               &
     &          LYMBQ(M)%RES%Int_Eng
            ENDIF
          ENDDO
          LSOLD(N)%RES%Int_Eng = Internal_Energy / LSOLD(N)%RES%Volume
!!
!! End of rigid/deformable element if-then-else-endif.
!!
        ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!! Note: NPNDT is a PARAMETER in "shared_common_data.f90"
!!
        TIMSIM%DTLSx = MAX (TIMSIM%DTLSx,LSOLD(N)%RES%DTelt)
        i = 0
        FOUND = .FALSE.
        DO WHILE (.NOT.FOUND .AND. i.LT.NPNDT)
          i = i + 1
          FOUND = LSOLD(N)%RES%DTelt .LT. TIMSIM%DTLYS(i)
        ENDDO
        IF (FOUND) THEN
          IF (i .LT. NPNDT) THEN
            DO j = NPNDT-1,i,-1
              TIMSIM%DTLYS(j + 1) = TIMSIM%DTLYS(j)
              TIMSIM%LSold(j + 1) = TIMSIM%LSold(j)
            ENDDO
          ENDIF
          TIMSIM%DTLYS(i) = LSOLD(N)%RES%DTelt
          TIMSIM%LSold(i) = N
        ENDIF
!!
!! End of element do-loop
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE LSOLD_INTERNAL_FORCES
!!
!! Copyright (c) by KEY Associates;  4-MAR-1994 18:07:47.17
!!
      USE shared_common_data
      USE lsold_
      USE hexah_, ONLY: LSHEX, hexah_type
      USE membq_, ONLY: LSMBQ, membq_type
      USE layering_
      USE section_2d_
      USE material_
      USE node_
      USE motion_
      USE force_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Local record structures. Note: MXLSL is a PARAMETER in "lsold_"
!!
      TYPE (motion_type), DIMENSION(8)     :: ELEMENT_MOTION
      TYPE (force_type),  DIMENSION(8)     :: ELEMENT_FORCE
      TYPE (motion_type), DIMENSION(8)     :: LAYER_MOTION
      TYPE (hexah_type),  DIMENSION(MXLSL) :: LYHEX
      TYPE (membq_type),  DIMENSION(MXLSL) :: LYMBQ
!!
      INTEGER                                                                  &
     &          Isv,                                                           &
     &          SecID,                                                         &
     &          LupID,                                                         &
     &          MatID
      REAL(KIND(0D0))                                                          &
     &          Normal_Stress,                                                 &
     &          Internal_Energy,                                               &
     &          B1(4),B2(4),T1(4),T2(4),                                       &
     &          X(8),Y(8),Z(8),Ux(8),Uy(8),Uz(8),Vx(8),Vy(8),Vz(8),            &
     &          R(8),S(8),T(8),Ur(8),Us(8),Ut(8),Vr(8),Vs(8),Vt(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Loop over all layered solid elements.
!!
      DO N = 1,NUMLS
!!
!! Test element subcycling index for time-to-update.
!!
        IF (MOD (TIMSIM%Cycle,LSOLD(N)%RES%ISI) .EQ. 0) THEN
!!
!! Distinguish between a rigid body element and a deformable element.
!!
          IF (LSOLD(N)%PAR%ParID .GE. 0) THEN
!!
!! Count element execution
!!
            COUNTER%LSOLD = COUNTER%LSOLD + 1
!!
!! Gather element motion; clear element internal force record.
!!
            DO i = 1,8
              ELEMENT_MOTION(i) = MOTION(LSOLD(N)%PAR%IX(i))
              ELEMENT_FORCE(i)  = force_type (0,0,0,0,0,0)
            ENDDO
!!
!! Increment element clock.
!!
            LSOLD(N)%RES%Time = LSOLD(N)%RES%Time + LSOLD(N)%RES%DTnext
!!
!! Scale nodal positions to current element time.
!!
            IF (CONTROL%SUBCYC .GT. 0) THEN
              DO i = 1,8
                QA = NODE(LSOLD(N)%PAR%IX(i))%Time - LSOLD(N)%RES%Time
                ELEMENT_MOTION(i)%Ux =                                         &
     &          ELEMENT_MOTION(i)%Ux - QA * ELEMENT_MOTION(i)%Vx
                ELEMENT_MOTION(i)%Uy =                                         &
     &          ELEMENT_MOTION(i)%Uy - QA * ELEMENT_MOTION(i)%Vy
                ELEMENT_MOTION(i)%Uz =                                         &
     &          ELEMENT_MOTION(i)%Uz - QA * ELEMENT_MOTION(i)%Vz
              ENDDO
            ENDIF
!!
!! Access element do-loop index for use in subroutine calls.
!!
            NEL = N
!!
!! Obtain element volume and generalized element dimensions for stability
!! calculations.
!!
            CALL LSOLD_VOLUME_OPERATOR ( NEL,ELEMENT_MOTION,DSQ )
!!
!! Compute volume strain for use latter in constructing the normal stress
!! thickness constraint.
!!
            Volume_Strain =                                                    &
     &        (LSOLD(N)%RES%Volume - LSOLD(N)%PAR%Volume) /                    &
     &        LSOLD(N)%PAR%Volume
!!
!! Define the current position of the layered solid nodal points.
!!
            DO i = 1,8
              X(i)  = ELEMENT_MOTION(i)%Px + ELEMENT_MOTION(i)%Ux
              Y(i)  = ELEMENT_MOTION(i)%Py + ELEMENT_MOTION(i)%Uy
              Z(i)  = ELEMENT_MOTION(i)%Pz + ELEMENT_MOTION(i)%Uz
              Ux(i) = ELEMENT_MOTION(i)%Ux
              Uy(i) = ELEMENT_MOTION(i)%Uy
              Uz(i) = ELEMENT_MOTION(i)%Uz
              Vx(i) = ELEMENT_MOTION(i)%Vx
              Vy(i) = ELEMENT_MOTION(i)%Vy
              Vz(i) = ELEMENT_MOTION(i)%Vz
            ENDDO
!!
!! Define an orthonormal set of basis vectors with the vectors R and S in
!! the plane of the element and T perpendicular to the element. Initially,
!! the vectors R and S are defined along element bisectors. As a last step,
!! the vector S is redefined to be perpendicular to R and T.
!!
            Rx = (X(3)-X(1)) + (X(2)-X(4)) + (X(7)-X(5)) + (X(6)-X(8))
            Ry = (Y(3)-Y(1)) + (Y(2)-Y(4)) + (Y(7)-Y(5)) + (Y(6)-Y(8))
            Rz = (Z(3)-Z(1)) + (Z(2)-Z(4)) + (Z(7)-Z(5)) + (Z(6)-Z(8))
            Sx = (X(3)-X(1)) - (X(2)-X(4)) + (X(7)-X(5)) - (X(6)-X(8))
            Sy = (Y(3)-Y(1)) - (Y(2)-Y(4)) + (Y(7)-Y(5)) - (Y(6)-Y(8))
            Sz = (Z(3)-Z(1)) - (Z(2)-Z(4)) + (Z(7)-Z(5)) - (Z(6)-Z(8))
            Rmag  = SQRT (Rx*Rx + Ry*Ry + Rz*Rz)
            Rx = Rx * (ONE / Rmag)
            Ry = Ry * (ONE / Rmag)
            Rz = Rz * (ONE / Rmag)
            Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
            Sx = Sx * (ONE / Smag)
            Sy = Sy * (ONE / Smag)
            Sz = Sz * (ONE / Smag)
!!
!! Define the unit vector T normal to the element.
!!
            Tx = Ry*Sz - Sy*Rz
            Ty = Rz*Sx - Sz*Rx
            Tz = Rx*Sy - Sx*Ry
            Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
            IF (Tmag .EQ. 0.0) THEN
              WRITE (MSG1,'(I8)') LSOLD(N)%PAR%EleID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LSOLD_INTERNAL_FORCES.001.00'//                         &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
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
!! Transform position X,Y,Z, displacement Ux,Uy,Uz, and translational
!! velocity Vx,Vy,Vz to the local R,S,T-coordinate system.
!!
            DO i = 1,8
              R(i)  = Rx*X(i)  + Ry*Y(i)  + Rz*Z(i)
              S(i)  = Sx*X(i)  + Sy*Y(i)  + Sz*Z(i)
              T(i)  = Tx*X(i)  + Ty*Y(i)  + Tz*Z(i)
              Ur(i) = Rx*Ux(i) + Ry*Uy(i) + Rz*Uz(i)
              Us(i) = Sx*Ux(i) + Sy*Uy(i) + Sz*Uz(i)
              Ut(i) = Tx*Ux(i) + Ty*Uy(i) + Tz*Uz(i)
              Vr(i) = Rx*Vx(i) + Ry*Vy(i) + Rz*Vz(i)
              Vs(i) = Sx*Vx(i) + Sy*Vy(i) + Sz*Vz(i)
              Vt(i) = Tx*Vx(i) + Ty*Vy(i) + Tz*Vz(i)
            ENDDO
!!
!! "Convert" global motion to local directions.
!!
            DO i = 1,8
              ELEMENT_MOTION(i)%Px = R(i) - Ur(i)
              ELEMENT_MOTION(i)%Py = S(i) - Us(i)
              ELEMENT_MOTION(i)%Pz = T(i) - Ut(i)
              ELEMENT_MOTION(i)%Ux = Ur(i)
              ELEMENT_MOTION(i)%Uy = Us(i)
              ELEMENT_MOTION(i)%Uz = Ut(i)
              ELEMENT_MOTION(i)%Vx = Vr(i)
              ELEMENT_MOTION(i)%Vy = Vs(i)
              ELEMENT_MOTION(i)%Vz = Vt(i)
            ENDDO
!!
!! Retrieve lay-up pointer LupID and number of layers NUMLY in lay-up.
!!
            LupID = LSOLD(N)%PAR%LupID
            NUMLY = LAYERING(LupID)%Number_of_Layers
!!
!! Gather "sub-element" data and store in local LYHEX and LYMBQ records.
!!
            DO m = 1,NUMLY
              IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
                LYHEX(m) = LSHEX(LSOLD(N)%PAR%ID(m))
              ELSE
                LYMBQ(m) = LSMBQ(LSOLD(N)%PAR%ID(m))
              ENDIF
            ENDDO
!!
!! In preparation for looping over all layers in layered solid element, set-up
!! the motion and weightings for transfering data between the layer nodal
!! points and the LSOLD nodal points.
!!
            DO i = 1,4
              LAYER_MOTION(i+4) = ELEMENT_MOTION(i)
              T1(i) = ONE
              T2(i) = 0.0
            ENDDO
            CSQmax = 0.0
!!
!! Loop over all layers in layered solid element processing each layer as
!! either a hexahedron or a membrane element.
!!
            DO M = 1,NUMLY
!!
!! Construct "sub-element." The do-loop starts at the layer next to the "1"
!! surface of the layered solid, IX(1:4), and ends at the layer next to the
!! "2" surface of the layered solid, IX(5:8). The weights T1 and T2 are for
!! use with data associted with the top of the layer, and B1 and B2 are for
!! use with data associated with the bottom of the layer. The "1" in B1 and
!! T1 refer to the "1" surface of the layered solid, IX(1:4). The "2" in B2
!! and T2 refer to the "2" surface of the layered solid, IX(5:8). The weights
!! are used to transfer results between the layer nodal points and the LSOLD
!! nodal points. Note that each corner of the LSOLD hexahedron can have dif-
!! ferent thicknesses associated with the layer, but the same number of layers
!! is expected at each corner. The interpolation at each corner is linear.
!!
              DO i = 1,4
                B1(i) = T1(i)
                B2(i) = T2(i)
                T1(i) = T1(i) - LSOLD(N)%RES%H(i,M)
                T2(i) = ONE - T1(i)
                LAYER_MOTION(i) = LAYER_MOTION(i+4)
                LAYER_MOTION(i+4)%Px = T1(i) * ELEMENT_MOTION(i  )%Px          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Px

                LAYER_MOTION(i+4)%Py = T1(i) * ELEMENT_MOTION(i  )%Py          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Py

                LAYER_MOTION(i+4)%Pz = T1(i) * ELEMENT_MOTION(i  )%Pz          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Pz

                LAYER_MOTION(i+4)%Ux = T1(i) * ELEMENT_MOTION(i  )%Ux          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Ux

                LAYER_MOTION(i+4)%Uy = T1(i) * ELEMENT_MOTION(i  )%Uy          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Uy

                LAYER_MOTION(i+4)%Uz = T1(i) * ELEMENT_MOTION(i  )%Uz          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Uz

                LAYER_MOTION(i+4)%Vx = T1(i) * ELEMENT_MOTION(i  )%Vx          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Vx

                LAYER_MOTION(i+4)%Vy = T1(i) * ELEMENT_MOTION(i  )%Vy          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Vy

                LAYER_MOTION(i+4)%Vz = T1(i) * ELEMENT_MOTION(i  )%Vz          &
     &                               + T2(i) * ELEMENT_MOTION(i+4)%Vz
              ENDDO
!!
!! Distinguish between hexahedron and quadrilateral membrane sub-elements.
!!
              IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
!!
!! Case (a): HEXAHERDAL SUBLAYER
!! Process layer as hexahedron. Retrieve material ID and state variable
!! starting location.
!!
                MatID = LYHEX(M)%PAR%MatID
                Isv   = LYHEX(M)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation.
!!
                CALL LYHEX_GRADIENT_OPERATOR (LYHEX(M),LAYER_MOTION)
!!
!! Update sound speeds and stress.
!!
                MTRL_TYPE = MATERIAL(MatID)%Type
                SELECT CASE (MTRL_TYPE)
                CASE (40)
                  CALL MATERIAL_40                                             &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
                CASE (41)
                  CALL MATERIAL_41                                             &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID,RCL2,RCS2                                        &
     &                  )
                CASE (45)
                  CALL MATERIAL_45                                             &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
                CASE (47)
                  CALL MATERIAL_47                                             &
     &                  (                                                      &
     &                  LYHEX(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  Dxx,Dyy,Dxy,Dxz,Dyz,Wxy,                               &
     &                  MatID                                                  &
     &                  )
                END SELECT
!!
!! Compute Normal_Stress to add to plane stress material model results.
!! The Normal_Stress is a "thickness constraint" designed to maintain a
!! constant element volume.
!!
                Normal_Stress = SOUND_SPEED%RCL2 * Volume_Strain
!!
!! Impose thickness constraint (actually a constant volume constraint).
!!
                LYHEX(M)%RES%STRESS(3) = Normal_Stress
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
                IF (MATERIAL(MatID)%PVAL(5) .GT. 0.0) THEN
                  CALL LYHEX_HOURGLASS_FORCES (LYHEX(M),MatID)
                ENDIF
!!
!! Divergence operator evaluated at time n+1.
!!
!!!               IF (CONTROL%MIDINT .NE. 0) THEN
!!!                 CALL LYHEX_DIVERGENCE_OPERATOR (LYHEX(M),LAYER_MOTION)
!!!               ENDIF
!!
!! Update stress divergence, viscosity stress and time step.
!!
                CALL LYHEX_STRESS_DIVERGENCE (LYHEX(M),MatID)
!!
!! Collect maximum layer tangent sound speed squared. (CSQ = C**2)
!!
                Density = MATERIAL(MatID)%PVAL(1)
                CSQtmp = SOUND_SPEED%RCL2/Density
                CSQmax = MAX (CSQmax,CSQtmp)

              ELSE
!!
!! Case (b): QUADRILATERAL MEMBRANE SUBLAYER
!! Process layer as a quadrilateral. Retrieve material ID and state variable
!! starting location.
!!
                MatID = LYMBQ(M)%PAR%MatID
                Isv   = LYMBQ(M)%PAR%Isv
!!
!! Gradient operator, stretching, and rotation.
!!
                DO i = 1,4
                  LAYER_MOTION(i)%Px =                                         &
     &              0.5*(LAYER_MOTION(i)%Px + LAYER_MOTION(i+4)%Px)
                  LAYER_MOTION(i)%Py =                                         &
     &              0.5*(LAYER_MOTION(i)%Py + LAYER_MOTION(i+4)%Py)
                  LAYER_MOTION(i)%Pz =                                         &
     &              0.5*(LAYER_MOTION(i)%Pz + LAYER_MOTION(i+4)%Pz)
                  LAYER_MOTION(i)%Ux =                                         &
     &              0.5*(LAYER_MOTION(i)%Ux + LAYER_MOTION(i+4)%Ux)
                  LAYER_MOTION(i)%Uy =                                         &
     &              0.5*(LAYER_MOTION(i)%Uy + LAYER_MOTION(i+4)%Uy)
                  LAYER_MOTION(i)%Uz =                                         &
     &              0.5*(LAYER_MOTION(i)%Uz + LAYER_MOTION(i+4)%Uz)
                  LAYER_MOTION(i)%Vx =                                         &
     &              0.5*(LAYER_MOTION(i)%Vx + LAYER_MOTION(i+4)%Vx)
                  LAYER_MOTION(i)%Vy =                                         &
     &              0.5*(LAYER_MOTION(i)%Vy + LAYER_MOTION(i+4)%Vy)
                  LAYER_MOTION(i)%Vz =                                         &
     &              0.5*(LAYER_MOTION(i)%Vz + LAYER_MOTION(i+4)%Vz)
                ENDDO
                CALL LYMBQ_GRADIENT_OPERATOR (LYMBQ(M),LAYER_MOTION)
!!
!! Update sound speeds and stress.
!!
                MTRL_TYPE = MATERIAL(MatID)%Type
                SELECT CASE (MTRL_TYPE)
                CASE (20)
                  CALL MATERIAL_20                                             &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
                CASE (21)
                  CALL MATERIAL_21                                             &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
                 CASE (22)
                   CALL MATERIAL_22                                            &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
                 CASE (25)
                   CALL MATERIAL_25                                            &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
                 CASE (27)
                   CALL MATERIAL_27                                            &
     &                  (                                                      &
     &                  LYMBQ(M)%RES%STRESS,                                   &
     &                  STATE_VARIABLES(Isv),                                  &
     &                  LYMBQ(M)%RES%Int_Eng,                                  &
     &                  LSOLD(N)%RES%DTnext,                                   &
     &                  MatID                                                  &
     &                  )
                 END SELECT
!!
!! Update divergence operator.
!!
                CALL LYMBQ_DIVERGENCE_OPERATOR (LYMBQ(M),LAYER_MOTION)
!!
!! Update hourglass control forces provided hourglass stiffness is non-zero.
!!
                IF (MATERIAL(MatID)%PVAL(5) .NE. 0.0) THEN
                  CALL LYMBQ_HOURGLASS_FORCES (LYMBQ(M),MatID)
                ENDIF
!!
!! Update stress divergence, viscosity stress and time step.
!!
                SecID = LAYERING(LupID)%Ltype(M)
                CALL LYMBQ_STRESS_DIVERGENCE (LYMBQ(M),SecID,MatID)
!!
!! Collect maximum layer tangent sound speed squared. (CSQ = C**2)
!!
                Density = MATERIAL(MatID)%PVAL(1)
                CSQtmp = SOUND_SPEED%RCL2/Density
                CSQmax = MAX (CSQmax,CSQtmp)
!!
!! End of layer type if-test.
!!
              ENDIF
!!
!! End of layer do-loop.
!!
            ENDDO
!!
!! Update global "sub-element" storage.
!!
            DO M = 1,NUMLY
              IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
                LSHEX(LSOLD(N)%PAR%ID(M)) = LYHEX(M)
              ELSE
                LSMBQ(LSOLD(N)%PAR%ID(M)) = LYMBQ(M)
              ENDIF
            ENDDO
!!
!! Transfer force data from the layer nodal points to the LSOLD nodal points.
!! In preparation for looping over all of the layers in a layered solid
!! element, set-up the weightings for transferring data between the layer
!! nodal points and the LSOLD nodal points.
!!
            T1(1) = 0.0
            T1(2) = 0.0
            T1(3) = 0.0
            T1(4) = 0.0
!!
!! Accumulate nodal forces. Distribute layer internal forces to LSOLD nodal poin
!! based on force and moment equilibrium.
!!
            DO i = 1,4
              DO M = 1,NUMLY
                B1(i) = T1(i)
                T1(i) = T1(i) + LSOLD(N)%RES%H(i,M)
!!
!! Distinguish between a hexahedron layer and membrane layer.
!!
                IF (LAYERING(LupID)%Ltype(m) .EQ. 0) THEN
!!
!! Forces at the bottom of the layer distributed to the LSOLD nodes.
!!
                  Fx = B1(i) * LYHEX(M)%RES%Xint(i)
                  Fy = B1(i) * LYHEX(M)%RES%Yint(i)
                  Fz = B1(i) * LYHEX(M)%RES%Zint(i)
                  ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                  ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                  ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                  Fx = LYHEX(M)%RES%Xint(i) - Fx
                  Fy = LYHEX(M)%RES%Yint(i) - Fy
                  Fz = LYHEX(M)%RES%Zint(i) - Fz
                  ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                  ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                  ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz
!!
!! Forces at the top of the layer distributed to the LSOLD nodes.
!!
                  Fx = T1(i) * LYHEX(M)%RES%Xint(i+4)
                  Fy = T1(i) * LYHEX(M)%RES%Yint(i+4)
                  Fz = T1(i) * LYHEX(M)%RES%Zint(i+4)
                  ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                  ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                  ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                  Fx = LYHEX(M)%RES%Xint(i+4) - Fx
                  Fy = LYHEX(M)%RES%Yint(i+4) - Fy
                  Fz = LYHEX(M)%RES%Zint(i+4) - Fz
                  ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                  ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                  ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz

                ELSE
!!
!! Forces in the membrane layer distributed to the LSOLD nodes.
!!
                  Fx = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Xint(i)
                  Fy = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Yint(i)
                  Fz = 0.5 * (B1(i) + T1(i)) * LYMBQ(M)%RES%Zint(i)
                  ELEMENT_FORCE(i+4)%Xint = ELEMENT_FORCE(i+4)%Xint + Fx
                  ELEMENT_FORCE(i+4)%Yint = ELEMENT_FORCE(i+4)%Yint + Fy
                  ELEMENT_FORCE(i+4)%Zint = ELEMENT_FORCE(i+4)%Zint + Fz
                  Fx = LYMBQ(M)%RES%Xint(i) - Fx
                  Fy = LYMBQ(M)%RES%Yint(i) - Fy
                  Fz = LYMBQ(M)%RES%Zint(i) - Fz
                  ELEMENT_FORCE(i)%Xint = ELEMENT_FORCE(i)%Xint + Fx
                  ELEMENT_FORCE(i)%Yint = ELEMENT_FORCE(i)%Yint + Fy
                  ELEMENT_FORCE(i)%Zint = ELEMENT_FORCE(i)%Zint + Fz

                ENDIF
              ENDDO
            ENDDO
!!
!! Transfer element nodal forces to element data structure and at the same
!! time transform the element nodal forces from the local r,s,t-coordinates
!! to the global x,y,z-coordinates.
!!
            DO i = 1,8
              LSOLD(N)%RES%Xint(i) = Rx*ELEMENT_FORCE(i)%Xint                  &
     &                             + Sx*ELEMENT_FORCE(i)%Yint                  &
     &                             + Tx*ELEMENT_FORCE(i)%Zint
              LSOLD(N)%RES%Yint(i) = Ry*ELEMENT_FORCE(i)%Xint                  &
     &                             + Sy*ELEMENT_FORCE(i)%Yint                  &
     &                             + Ty*ELEMENT_FORCE(i)%Zint
              LSOLD(N)%RES%Zint(i) = Rz*ELEMENT_FORCE(i)%Xint                  &
     &                             + Sz*ELEMENT_FORCE(i)%Yint                  &
     &                             + Tz*ELEMENT_FORCE(i)%Zint
            ENDDO
!!
!! Construct element critical time step. (In-plane the layer moduli act
!! in parallel; perpendicular the layer moduli act in series.)
!!
            LSOLD(N)%RES%DTelt = SQRT (DSQ/CSQmax)
!!
!! Save sound speed data for nonreflecting boundary condition.
!!
            IF (NUMNR .NE. 0) THEN
              Itype = -1  !  Layered Solid Hexahedron
!! fix-up sound_speed data!!
              CALL SAVE_COMPLIANCE_FOR_NRBC ( NEL,Itype )
            ENDIF
!!
!! Compute internal energy density of layered solid from sub-elements.
!!
            Internal_Energy = 0.0
            DO M = 1,NUMLY
              IF (LAYERING(LupID)%Ltype(M) .EQ. 0) THEN
                Internal_Energy = Internal_Energy +                            &
     &            LYHEX(M)%RES%Volume * LYHEX(M)%RES%Int_Eng
              ELSE
                Internal_Energy = Internal_Energy +                            &
     &            LYMBQ(M)%PAR%Area *                                          &
     &            SECTION_2D(LAYERING(LupID)%Ltype(m))%Thickness *             &
     &            LYMBQ(M)%RES%Int_Eng
              ENDIF
            ENDDO
            LSOLD(N)%RES%Int_Eng = Internal_Energy / LSOLD(N)%RES%Volume
!!
!! End of skip-rigid-body-element if-test.
!!
          ENDIF
!!
!! Compare element time step with minimum and maximum time step.
!!
          TIMSIM%DTLSx = MAX (TIMSIM%DTLSx,LSOLD(N)%RES%DTelt)
          IF (LSOLD(N)%RES%DTelt .LT. TIMSIM%DTLYS(1)) THEN
            TIMSIM%DTLYS(1) = LSOLD(N)%RES%DTelt
            TIMSIM%LSold(1) = N
          ENDIF
!!
!! End of time-to-subcycle if-test.
!!
        ENDIF
!!
!! End of element do-loop
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE LSOLD_VOLUME_OPERATOR ( NEL,ELEMENT_MOTION,DSQ )
!!
!! Copyright (c) by KEY Associates;  4-MAR-1994 18:07:57.09
!!
      USE shared_common_data
      USE lsold_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (motion_type), DIMENSION(8) :: ELEMENT_MOTION
!!
      INTEGER                                                                  &
     &          I2(8),I3(8),I4(8),I5(8),I6(8),I8(8)
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8),Bx(8),By(8),Bz(8)
!!
      DATA                                                                     &
     &          I2 /2,3,4,1,8,5,6,7/, I3 /3,4,1,2,7,8,5,6/,                    &
     &          I4 /4,1,2,3,6,7,8,5/, I5 /5,6,7,8,1,2,3,4/,                    &
     &          I6 /6,7,8,5,4,1,2,3/, I8 /8,5,6,7,2,3,4,1/
!!
!! Compute position of nodal points at time n+1.
!!
      DO i = 1,8
        X(i) = ELEMENT_MOTION(i)%Px + ELEMENT_MOTION(i)%Ux
        Y(i) = ELEMENT_MOTION(i)%Py + ELEMENT_MOTION(i)%Uy
        Z(i) = ELEMENT_MOTION(i)%Pz + ELEMENT_MOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,8
        Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))               &
     &           +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))               &
     &           +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))               &
     &           +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                                 &
     &           +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                                 &
     &           +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))               &
     &           +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))               &
     &           +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))               &
     &           +Z(I3(i))*(X(I2(i))-X(I4(i)))                                 &
     &           +Z(I6(i))*(X(I5(i))-X(I2(i)))                                 &
     &           +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))               &
     &           +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))               &
     &           +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))               &
     &           +X(I3(i))*(Y(I2(i))-Y(I4(i)))                                 &
     &           +X(I6(i))*(Y(I5(i))-Y(I2(i)))                                 &
     &           +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume.
!!
      LSOLD(NEL)%RES%Volume = 0.0
      DO i = 1,8
        LSOLD(NEL)%RES%Volume = LSOLD(NEL)%RES%Volume + Bx(i)*X(i)
      ENDDO
!!
!! Calculate generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,8
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      DSQ = (LSOLD(NEL)%RES%Volume * LSOLD(NEL)%RES%Volume) / (Dx+Dx)
!!
      RETURN
      END
!!_
      SUBROUTINE LYHEX_GRADIENT_OPERATOR (LYHEX,LAYER_MOTION)
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE hexah_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (hexah_type) :: LYHEX
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
!!
      INTEGER                                                                  &
     &          I2(8),I3(8),I4(8),I5(8),I6(8),I8(8)
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,4,1,8,5,6,7/, I3 /3,4,1,2,7,8,5,6/,                    &
     &          I4 /4,1,2,3,6,7,8,5/, I5 /5,6,7,8,1,2,3,4/,                    &
     &          I6 /6,7,8,5,4,1,2,3/, I8 /8,5,6,7,2,3,4,1/
!!
!! Compute position of nodal points at time n+1 (whole-interval evaluation).
!!
      DO i = 1,8
        X(i) = LAYER_MOTION(i)%Px + LAYER_MOTION(i)%Ux
        Y(i) = LAYER_MOTION(i)%Py + LAYER_MOTION(i)%Uy
        Z(i) = LAYER_MOTION(i)%Pz + LAYER_MOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,8
        Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))               &
     &           +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))               &
     &           +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))               &
     &           +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                                 &
     &           +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                                 &
     &           +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))               &
     &           +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))               &
     &           +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))               &
     &           +Z(I3(i))*(X(I2(i))-X(I4(i)))                                 &
     &           +Z(I6(i))*(X(I5(i))-X(I2(i)))                                 &
     &           +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))               &
     &           +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))               &
     &           +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))               &
     &           +X(I3(i))*(Y(I2(i))-Y(I4(i)))                                 &
     &           +X(I6(i))*(Y(I5(i))-Y(I2(i)))                                 &
     &           +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      LYHEX%RES%Volume = 0.0
      DO i = 1,8
        LYHEX%RES%Volume = LYHEX%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = ONE / LYHEX%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,8
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
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
      DO i = 1,8
        Vxx = Vxx + LAYER_MOTION(i)%Vx * Bx(i)
        Vyx = Vyx + LAYER_MOTION(i)%Vy * Bx(i)
        Vzx = Vzx + LAYER_MOTION(i)%Vz * Bx(i)
        Vxy = Vxy + LAYER_MOTION(i)%Vx * By(i)
        Vyy = Vyy + LAYER_MOTION(i)%Vy * By(i)
        Vzy = Vzy + LAYER_MOTION(i)%Vz * By(i)
        Vxz = Vxz + LAYER_MOTION(i)%Vx * Bz(i)
        Vyz = Vyz + LAYER_MOTION(i)%Vy * Bz(i)
        Vzz = Vzz + LAYER_MOTION(i)%Vz * Bz(i)
      ENDDO
      Vxx = Vxx * Vin
      Vyx = Vyx * Vin
      Vzx = Vzx * Vin
      Vxy = Vxy * Vin
      Vyy = Vyy * Vin
      Vzy = Vzy * Vin
      Vxz = Vxz * Vin
      Vyz = Vyz * Vin
      Vzz = Vzz * Vin
!!
      Dxx = Vxx
      Dyy = Vyy
      Dzz = Vzz
      Dxy = 0.5 * (Vxy + Vyx)
      Dxz = 0.5 * (Vxz + Vzx)
      Dyz = 0.5 * (Vyz + Vzy)
      Wxy = Vxy - Dxy
      Wxz = Vxz - Dxz
      Wyz = Vyz - Dyz
!!
!! Anti-hourglass gradients.
!!
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(4)) -(X(5)+X(6)) +(X(7)+X(8))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(4)) -(X(5)-X(6)) +(X(7)-X(8))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(4)) -(Y(5)+Y(6)) +(Y(7)+Y(8))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(4)) -(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(4)) -(Z(5)+Z(6)) +(Z(7)+Z(8))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(4)) -(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
!!
      DO i = 1,4
        Gx(i) = -(Vxx*Hx(i) + Vxy*Hy(i) + Vxz*Hz(i)) * Delta
        Gy(i) = -(Vyx*Hx(i) + Vyy*Hy(i) + Vyz*Hz(i)) * Delta
        Gz(i) = -(Vzx*Hx(i) + Vzy*Hy(i) + Vzz*Hz(i)) * Delta
      ENDDO
      Gx(1) = Gx(1) +                                                          &
     &          ( (LAYER_MOTION(1)%Vx + LAYER_MOTION(2)%Vx)                    &
     &          - (LAYER_MOTION(3)%Vx + LAYER_MOTION(4)%Vx)                    &
     &          - (LAYER_MOTION(5)%Vx + LAYER_MOTION(6)%Vx)                    &
     &          + (LAYER_MOTION(7)%Vx + LAYER_MOTION(8)%Vx) )*Delta
      Gx(2) = Gx(2) +                                                          &
     &          ( (LAYER_MOTION(1)%Vx - LAYER_MOTION(2)%Vx)                    &
     &          - (LAYER_MOTION(3)%Vx - LAYER_MOTION(4)%Vx)                    &
     &          - (LAYER_MOTION(5)%Vx - LAYER_MOTION(6)%Vx)                    &
     &          + (LAYER_MOTION(7)%Vx - LAYER_MOTION(8)%Vx) )*Delta
      Gx(3) = Gx(3) +                                                          &
     &          ( (LAYER_MOTION(1)%Vx - LAYER_MOTION(2)%Vx)                    &
     &          + (LAYER_MOTION(3)%Vx - LAYER_MOTION(4)%Vx)                    &
     &          + (LAYER_MOTION(5)%Vx - LAYER_MOTION(6)%Vx)                    &
     &          + (LAYER_MOTION(7)%Vx - LAYER_MOTION(8)%Vx) )*Delta
      Gx(4) = Gx(4) -                                                          &
     &          ( (LAYER_MOTION(1)%Vx - LAYER_MOTION(2)%Vx)                    &
     &          + (LAYER_MOTION(3)%Vx - LAYER_MOTION(4)%Vx)                    &
     &          - (LAYER_MOTION(5)%Vx - LAYER_MOTION(6)%Vx)                    &
     &          - (LAYER_MOTION(7)%Vx - LAYER_MOTION(8)%Vx) )*Delta
!!
      Gy(1) = Gy(1) +                                                          &
     &          ( (LAYER_MOTION(1)%Vy + LAYER_MOTION(2)%Vy)                    &
     &          - (LAYER_MOTION(3)%Vy + LAYER_MOTION(4)%Vy)                    &
     &          - (LAYER_MOTION(5)%Vy + LAYER_MOTION(6)%Vy)                    &
     &          + (LAYER_MOTION(7)%Vy + LAYER_MOTION(8)%Vy) )*Delta
      Gy(2) = Gy(2) +                                                          &
     &          ( (LAYER_MOTION(1)%Vy - LAYER_MOTION(2)%Vy)                    &
     &          - (LAYER_MOTION(3)%Vy - LAYER_MOTION(4)%Vy)                    &
     &          - (LAYER_MOTION(5)%Vy - LAYER_MOTION(6)%Vy)                    &
     &          + (LAYER_MOTION(7)%Vy - LAYER_MOTION(8)%Vy) )*Delta
      Gy(3) = Gy(3) +                                                          &
     &          ( (LAYER_MOTION(1)%Vy - LAYER_MOTION(2)%Vy)                    &
     &          + (LAYER_MOTION(3)%Vy - LAYER_MOTION(4)%Vy)                    &
     &          + (LAYER_MOTION(5)%Vy - LAYER_MOTION(6)%Vy)                    &
     &          + (LAYER_MOTION(7)%Vy - LAYER_MOTION(8)%Vy) )*Delta
      Gy(4) = Gy(4) -                                                          &
     &          ( (LAYER_MOTION(1)%Vy - LAYER_MOTION(2)%Vy)                    &
     &          + (LAYER_MOTION(3)%Vy - LAYER_MOTION(4)%Vy)                    &
     &          - (LAYER_MOTION(5)%Vy - LAYER_MOTION(6)%Vy)                    &
     &          - (LAYER_MOTION(7)%Vy - LAYER_MOTION(8)%Vy) )*Delta
!!
      Gz(1) = Gz(1) +                                                          &
     &          ( (LAYER_MOTION(1)%Vz + LAYER_MOTION(2)%Vz)                    &
     &          - (LAYER_MOTION(3)%Vz + LAYER_MOTION(4)%Vz)                    &
     &          - (LAYER_MOTION(5)%Vz + LAYER_MOTION(6)%Vz)                    &
     &          + (LAYER_MOTION(7)%Vz + LAYER_MOTION(8)%Vz) )*Delta
      Gz(2) = Gz(2) +                                                          &
     &          ( (LAYER_MOTION(1)%Vz - LAYER_MOTION(2)%Vz)                    &
     &          - (LAYER_MOTION(3)%Vz - LAYER_MOTION(4)%Vz)                    &
     &          - (LAYER_MOTION(5)%Vz - LAYER_MOTION(6)%Vz)                    &
     &          + (LAYER_MOTION(7)%Vz - LAYER_MOTION(8)%Vz) )*Delta
      Gz(3) = Gz(3) +                                                          &
     &          ( (LAYER_MOTION(1)%Vz - LAYER_MOTION(2)%Vz)                    &
     &          + (LAYER_MOTION(3)%Vz - LAYER_MOTION(4)%Vz)                    &
     &          + (LAYER_MOTION(5)%Vz - LAYER_MOTION(6)%Vz)                    &
     &          + (LAYER_MOTION(7)%Vz - LAYER_MOTION(8)%Vz) )*Delta
      Gz(4) = Gz(4) -                                                          &
     &          ( (LAYER_MOTION(1)%Vz - LAYER_MOTION(2)%Vz)                    &
     &          + (LAYER_MOTION(3)%Vz - LAYER_MOTION(4)%Vz)                    &
     &          - (LAYER_MOTION(5)%Vz - LAYER_MOTION(6)%Vz)                    &
     &          - (LAYER_MOTION(7)%Vz - LAYER_MOTION(8)%Vz) )*Delta
!!
      RETURN
      END
!!_
      SUBROUTINE LYHEX_DIVERGENCE_OPERATOR (LYHEX,LAYER_MOTION)
!!
!! Copyright (c) by KEY Associates, 16-FEB-1991 20:31:21
!!
      USE shared_common_data
      USE hexah_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (hexah_type) :: LYHEX
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
!!
      INTEGER                                                                  &
     &          I2(8),I3(8),I4(8),I5(8),I6(8),I8(8)
      REAL(KIND(0D0))                                                          &
     &          X(8),Y(8),Z(8)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
      DATA                                                                     &
     &          I2 /2,3,4,1,8,5,6,7/, I3 /3,4,1,2,7,8,5,6/,                    &
     &          I4 /4,1,2,3,6,7,8,5/, I5 /5,6,7,8,1,2,3,4/,                    &
     &          I6 /6,7,8,5,4,1,2,3/, I8 /8,5,6,7,2,3,4,1/
!!
!! Current position of nodal points, time n+1.
!!
      DO i = 1,8
        X(i) = LAYER_MOTION(i)%Px + LAYER_MOTION(i)%Ux
        Y(i) = LAYER_MOTION(i)%Py + LAYER_MOTION(i)%Uy
        Z(i) = LAYER_MOTION(i)%Pz + LAYER_MOTION(i)%Uz
      ENDDO
!!
!! Gradient operators.
!!
      DO i = 1,8
        Bx(i) = ( Y(I2(i))*(Z(I6(i))-Z(I3(i))+Z(I5(i))-Z(I4(i)))               &
     &           +Y(I4(i))*(Z(I3(i))-Z(I8(i))+Z(I2(i))-Z(I5(i)))               &
     &           +Y(I5(i))*(Z(I8(i))-Z(I6(i))+Z(I4(i))-Z(I2(i)))               &
     &           +Y(I3(i))*(Z(I2(i))-Z(I4(i)))                                 &
     &           +Y(I6(i))*(Z(I5(i))-Z(I2(i)))                                 &
     &           +Y(I8(i))*(Z(I4(i))-Z(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        By(i) = ( Z(I2(i))*(X(I6(i))-X(I3(i))+X(I5(i))-X(I4(i)))               &
     &           +Z(I4(i))*(X(I3(i))-X(I8(i))+X(I2(i))-X(I5(i)))               &
     &           +Z(I5(i))*(X(I8(i))-X(I6(i))+X(I4(i))-X(I2(i)))               &
     &           +Z(I3(i))*(X(I2(i))-X(I4(i)))                                 &
     &           +Z(I6(i))*(X(I5(i))-X(I2(i)))                                 &
     &           +Z(I8(i))*(X(I4(i))-X(I5(i))) ) * 0.08333333333
      ENDDO
      DO i = 1,8
        Bz(i) = ( X(I2(i))*(Y(I6(i))-Y(I3(i))+Y(I5(i))-Y(I4(i)))               &
     &           +X(I4(i))*(Y(I3(i))-Y(I8(i))+Y(I2(i))-Y(I5(i)))               &
     &           +X(I5(i))*(Y(I8(i))-Y(I6(i))+Y(I4(i))-Y(I2(i)))               &
     &           +X(I3(i))*(Y(I2(i))-Y(I4(i)))                                 &
     &           +X(I6(i))*(Y(I5(i))-Y(I2(i)))                                 &
     &           +X(I8(i))*(Y(I4(i))-Y(I5(i))) ) * 0.08333333333
      ENDDO
!!
!! Calculate current element volume, Vin = 1.0/Volume
!!
      LYHEX%RES%Volume = 0.0
      DO i = 1,8
        LYHEX%RES%Volume = LYHEX%RES%Volume + Bx(i)*X(i)
      ENDDO
      Vin = ONE / LYHEX%RES%Volume
!!
!! Calculate inverse of generalized element dimension.
!!
      Dx = 0.0
      DO i = 1,8
        Dx = Dx + Bx(i)*Bx(i) + By(i)*By(i) + Bz(i)*Bz(i)
      ENDDO
      Delta = Vin * SQRT(Dx+Dx)
!!
!! Anti-hourglass divergence operator.
!!
      Hx(1) =  (X(1)+X(2)) -(X(3)+X(4)) -(X(5)+X(6)) +(X(7)+X(8))
      Hx(2) =  (X(1)-X(2)) -(X(3)-X(4)) -(X(5)-X(6)) +(X(7)-X(8))
      Hx(3) =  (X(1)-X(2)) +(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
      Hx(4) = -(X(1)-X(2)) -(X(3)-X(4)) +(X(5)-X(6)) +(X(7)-X(8))
!!
      Hy(1) =  (Y(1)+Y(2)) -(Y(3)+Y(4)) -(Y(5)+Y(6)) +(Y(7)+Y(8))
      Hy(2) =  (Y(1)-Y(2)) -(Y(3)-Y(4)) -(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(3) =  (Y(1)-Y(2)) +(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
      Hy(4) = -(Y(1)-Y(2)) -(Y(3)-Y(4)) +(Y(5)-Y(6)) +(Y(7)-Y(8))
!!
      Hz(1) =  (Z(1)+Z(2)) -(Z(3)+Z(4)) -(Z(5)+Z(6)) +(Z(7)+Z(8))
      Hz(2) =  (Z(1)-Z(2)) -(Z(3)-Z(4)) -(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(3) =  (Z(1)-Z(2)) +(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
      Hz(4) = -(Z(1)-Z(2)) -(Z(3)-Z(4)) +(Z(5)-Z(6)) +(Z(7)-Z(8))
!!
      RETURN
      END
!!_
      SUBROUTINE LYHEX_MASS (LYHEX,MatID,LAYER_NODE,LAYER_MOTION)
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 19:47:17
!!
!! Purpose: Compute element mass matrix. The simplest of mass lumpings
!! is used - one eighth at each node.
!!
      USE shared_common_data
      USE hexah_
      USE material_
      USE node_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (hexah_type) :: LYHEX
      TYPE (node_type),   DIMENSION(8) :: LAYER_NODE
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Compute one eighth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = 0.125 * Density * LYHEX%RES%Volume
!!
!! Accumulate mass at each nodal point and accumulate mass properties
!! for each material domain.
!!
      DO i = 1,8
        LAYER_NODE(LYHEX%PAR%IX(i))%Mass =                                     &
     &    LAYER_NODE(LYHEX%PAR%IX(i))%Mass + QMass
        MATERIAL(MatID)%Mass = MATERIAL(MatID)%Mass + QMass
        Px = LAYER_MOTION(i)%Px
        Py = LAYER_MOTION(i)%Py
        Pz = LAYER_MOTION(i)%Pz
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
      SUBROUTINE LYHEX_STRESS_DIVERGENCE (LYHEX,MatID )
!!
!! Copyright (c) by KEY Associates, 19-FEB-1991 20:03:56
!!
      USE shared_common_data
      USE hexah_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (hexah_type) :: LYHEX
!!
      REAL(KIND(0D0))                                                          &
     &          Fx(8),Fy(8),Fz(8),PGx(4),PGy(4),PGz(4)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
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
      SOUND_SPEED%Density =                                                    &
     &          Density * LYHEX%PAR%Volume / LYHEX%RES%Volume
!!
      CSQ = MAX(SOUND_SPEED%RCL2,SOUND_SPEED%RCS2) / SOUND_SPEED%Density
      IF (CSQ .GT. 0.0) THEN
        CX = SQRT(CSQ)
      ELSE
        WRITE (MSG1,'(I8)') LYHEX%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'LYHEX_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'HXEL (8-Node Hexahedron) Element ID:'//MSG1//           &
     &          MSGL//'Sound Speed Imaginary (C**2 Is Zero Or Negative)'       &
     &          )
        CX = TIMSIM%DTlast * 1.0D-6
      ENDIF
!!
!! Calculate generalized element dimension.
!!
      Dx = ONE / Delta
!!
!! Artificial bulk viscosity pressure.
!!
      Dkk = Dxx + Dyy + Dzz
      IF (Dkk .LT. 0.0) THEN
        QC = Bulk_Qd * Bulk_Qd * Dx * ABS(Dkk) + Bulk_Ln * CX
        QP = SOUND_SPEED%Density * Dx * Dkk * QC
      ELSE
        QC = 0.0
        QP = 0.0
      ENDIF
!!
!! Hourglass stiffness forces combined with hourglass viscosity forces.
!! Note: The inverse of the generalized element dimension used here is
!! based on the geometry at time n+1.
!!
      Qhg = SOUND_SPEED%Density * Dx * HG_Visc * CX
      DO i = 1,4
        PGx(i) = (LYHEX%RES%Px(i) + Qhg*Gx(i)) * Delta
        PGy(i) = (LYHEX%RES%Py(i) + Qhg*Gy(i)) * Delta
        PGz(i) = (LYHEX%RES%Pz(i) + Qhg*Gz(i)) * Delta
      ENDDO
!!
!! Critical time step calculation.
!!
      LYHEX%RES%DTelt = Dx / (QC + SQRT (QC*QC + CX*CX))
!!
!! Divergence of the hourglass forces.
!!
      Fx(1) = ( (PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * LYHEX%RES%Volume
      Fx(2) = ( (PGx(1)-PGx(2))-(PGx(3)-PGx(4))) * LYHEX%RES%Volume
      Fx(3) = (-(PGx(1)+PGx(2))+(PGx(3)-PGx(4))) * LYHEX%RES%Volume
      Fx(4) = (-(PGx(1)-PGx(2))-(PGx(3)-PGx(4))) * LYHEX%RES%Volume
      Fx(5) = (-(PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * LYHEX%RES%Volume
      Fx(6) = (-(PGx(1)-PGx(2))-(PGx(3)+PGx(4))) * LYHEX%RES%Volume
      Fx(7) = ( (PGx(1)+PGx(2))+(PGx(3)+PGx(4))) * LYHEX%RES%Volume
      Fx(8) = ( (PGx(1)-PGx(2))-(PGx(3)+PGx(4))) * LYHEX%RES%Volume
!!
      Fy(1) = ( (PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * LYHEX%RES%Volume
      Fy(2) = ( (PGy(1)-PGy(2))-(PGy(3)-PGy(4))) * LYHEX%RES%Volume
      Fy(3) = (-(PGy(1)+PGy(2))+(PGy(3)-PGy(4))) * LYHEX%RES%Volume
      Fy(4) = (-(PGy(1)-PGy(2))-(PGy(3)-PGy(4))) * LYHEX%RES%Volume
      Fy(5) = (-(PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * LYHEX%RES%Volume
      Fy(6) = (-(PGy(1)-PGy(2))-(PGy(3)+PGy(4))) * LYHEX%RES%Volume
      Fy(7) = ( (PGy(1)+PGy(2))+(PGy(3)+PGy(4))) * LYHEX%RES%Volume
      Fy(8) = ( (PGy(1)-PGy(2))-(PGy(3)+PGy(4))) * LYHEX%RES%Volume
!!
      Fz(1) = ( (PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * LYHEX%RES%Volume
      Fz(2) = ( (PGz(1)-PGz(2))-(PGz(3)-PGz(4))) * LYHEX%RES%Volume
      Fz(3) = (-(PGz(1)+PGz(2))+(PGz(3)-PGz(4))) * LYHEX%RES%Volume
      Fz(4) = (-(PGz(1)-PGz(2))-(PGz(3)-PGz(4))) * LYHEX%RES%Volume
      Fz(5) = (-(PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * LYHEX%RES%Volume
      Fz(6) = (-(PGz(1)-PGz(2))-(PGz(3)+PGz(4))) * LYHEX%RES%Volume
      Fz(7) = ( (PGz(1)+PGz(2))+(PGz(3)+PGz(4))) * LYHEX%RES%Volume
      Fz(8) = ( (PGz(1)-PGz(2))-(PGz(3)+PGz(4))) * LYHEX%RES%Volume
!!
!! Divergence of the (mean) stress plus the orthogonalized correction Hi(1:4).
!!
      Qxx = LYHEX%RES%Stress(1) + QP
      Qyx = LYHEX%RES%Stress(4)
      Qzx = LYHEX%RES%Stress(5)
      DO j = 1,4
        Qxx = Qxx - Hx(j)*PGx(j)
        Qyx = Qyx - Hy(j)*PGx(j)
        Qzx = Qzx - Hz(j)*PGx(j)
      ENDDO
!!
      Qxy = LYHEX%RES%Stress(4)
      Qyy = LYHEX%RES%Stress(2) + QP
      Qzy = LYHEX%RES%Stress(6)
      DO j = 1,4
        Qxy = Qxy - Hx(j)*PGy(j)
        Qyy = Qyy - Hy(j)*PGy(j)
        Qzy = Qzy - Hz(j)*PGy(j)
      ENDDO
!!
      Qxz = LYHEX%RES%Stress(5)
      Qyz = LYHEX%RES%Stress(6)
      Qzz = LYHEX%RES%Stress(3) + QP
      DO j = 1,4
        Qxz = Qxz - Hx(j)*PGz(j)
        Qyz = Qyz - Hy(j)*PGz(j)
        Qzz = Qzz - Hz(j)*PGz(j)
      ENDDO
!!
!! Accumulate element divergence results.
!!
      DO i = 1,8
        LYHEX%RES%Xint(i) = Fx(i) + Bx(i)*Qxx + By(i)*Qyx + Bz(i)*Qzx
        LYHEX%RES%Yint(i) = Fy(i) + Bx(i)*Qxy + By(i)*Qyy + Bz(i)*Qzy
        LYHEX%RES%Zint(i) = Fz(i) + Bx(i)*Qxz + By(i)*Qyz + Bz(i)*Qzz
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE LYHEX_HOURGLASS_FORCES (LYHEX,MatID )
!!
!! Copyright (c) by KEY Associates, 26-MAY-1991 13:19:09
!!
!! Purpose: Increment stiffness based hourglass control forces.
!!
      USE shared_common_data
      USE hexah_
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (hexah_type) :: LYHEX
!!
      REAL(KIND(0D0))                                                          &
     &          Q1(4),Q2(4),Q3(4)
!!
      COMMON /SOLID/                                                           &
     &          Bx(8),By(8),Bz(8),Hx(4),Hy(4),Hz(4),Gx(4),Gy(4),Gz(4),         &
     &          Vxx,Vxy,Vxz,Vyx,Vyy,Vyz,Vzx,Vzy,Vzz,Delta,Dxx,Dyy,Dzz,         &
     &          Dxy,Dxz,Dyz,Wxy,Wxz,Wyz,Rotation(3,3)
!!
!! Rotation of hour glass forces based on spin, matches Jaumann stress
!! flux.
!!
      IF (CONTROL%POLARD .EQ. 0) THEN
        dWxy = LYHEX%RES%DTnext * Wxy
        dWxz = LYHEX%RES%DTnext * Wxz
        dWyz = LYHEX%RES%DTnext * Wyz
        DO i = 1,4
          Q1(i) =  dWxy * LYHEX%RES%Py(i) + dWxz * LYHEX%RES%Pz(i)
          Q2(i) = -dWxy * LYHEX%RES%Px(i) + dWyz * LYHEX%RES%Pz(i)
          Q3(i) = -dWxz * LYHEX%RES%Px(i) - dWyz * LYHEX%RES%Py(i)
        ENDDO
      ELSE
        DO i = 1,4
          Q1(i) = 0.0
          Q2(i) = 0.0
          Q3(i) = 0.0
        ENDDO
      ENDIF
!!
!! Increment elastic anti-hourglassing forces.
!!
        HG_Stiff = MATERIAL(MatID)%PVAL(5)
        Qe = LYHEX%RES%DTnext * HG_Stiff * SOUND_SPEED%RCS2
        DO i = 1,4
          LYHEX%RES%Px(i) = LYHEX%RES%Px(i) + Qe * Gx(i) + Q1(i)
          LYHEX%RES%Py(i) = LYHEX%RES%Py(i) + Qe * Gy(i) + Q2(i)
          LYHEX%RES%Pz(i) = LYHEX%RES%Pz(i) + Qe * Gz(i) + Q3(i)
        ENDDO
!!
!! Rotation of hour glass forces based on rotation from polar decomposi-
!! tion of the deformation gradient, matches Green-McInnes stress flux.
!! This is a rotation forward from time n to the end of the interval at
!! time n+1 and thus, uses the inverse of the rotation R  generated in
!! SYMMETRIC_TENSOR_ROTATION That is, Rt * P is used in place of R * P,
!! the transformation we would have used if we had been able to construct
!! the "forward" deformation gradient in SYMMETRIC_TENSOR_ROTATION.
!!
      IF (CONTROL%POLARD .GT. 0) THEN
        DO i = 1,4
          Q1(i) = Rotation(1,1) * LYHEX%RES%Px(i)                              &
     &          + Rotation(2,1) * LYHEX%RES%Py(i)                              &
     &          + Rotation(3,1) * LYHEX%RES%Pz(i)

          Q2(i) = Rotation(1,2) * LYHEX%RES%Px(i)                              &
     &          + Rotation(2,2) * LYHEX%RES%Py(i)                              &
     &          + Rotation(3,2) * LYHEX%RES%Pz(i)

          Q3(i) = Rotation(1,3) * LYHEX%RES%Px(i)                              &
     &          + Rotation(2,3) * LYHEX%RES%Py(i)                              &
     &          + Rotation(3,3) * LYHEX%RES%Pz(i)
        ENDDO
        DO i = 1,4
          LYHEX%RES%Px(i) = Q1(i)
          LYHEX%RES%Py(i) = Q2(i)
          LYHEX%RES%Pz(i) = Q3(i)
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE LYMBQ_GRADIENT_OPERATOR ( LYMBQ,LAYER_MOTION )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
      USE shared_common_data
      USE membq_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (membq_type) :: LYMBQ
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
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
      DT = 0.5 * LYMBQ%RES%DTnext
      DO i = 1,4
        X(i)  = LAYER_MOTION(i)%Px +                                           &
     &    (LAYER_MOTION(i)%Ux - DT * LAYER_MOTION(i)%Vx)
        Y(i)  = LAYER_MOTION(i)%Py +                                           &
     &    (LAYER_MOTION(i)%Uy - DT * LAYER_MOTION(i)%Vy)
        Z(i)  = LAYER_MOTION(i)%Pz +                                           &
     &    (LAYER_MOTION(i)%Uz - DT * LAYER_MOTION(i)%Vz)
        Vx(i) = LAYER_MOTION(i)%Vx
        Vy(i) = LAYER_MOTION(i)%Vy
        Vz(i) = LAYER_MOTION(i)%Vz
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
      Rx = Rx * (ONE / Rmag)
      Ry = Ry * (ONE / Rmag)
      Rz = Rz * (ONE / Rmag)
      Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (ONE / Smag)
      Sy = Sy * (ONE / Smag)
      Sz = Sz * (ONE / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') LYMBQ%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LYMBQ_GRADIENT_OPERATOR.001.00'//                       &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
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
      LYMBQ%RES%Area = 0.5 * (R13*S24 - S13*R24)
      Ain = ONE / LYMBQ%RES%Area
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
      DTnext = LYMBQ%RES%DTnext
      dBeta = ATAN (DTnext*(Vrs-Vsr)/(2.0 + DTnext*(Vrr+Vss)))
      LYMBQ%RES%Beta = LYMBQ%RES%Beta + dBeta
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
      SUBROUTINE LYMBQ_DIVERGENCE_OPERATOR ( LYMBQ,LAYER_MOTION )
!!
!! Migrated by: S W Key, 20-APR-1991 15:44:50
!!
      USE shared_common_data
      USE membq_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TYPE (membq_type) :: LYMBQ
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
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
        X(i)  = LAYER_MOTION(i)%Px + LAYER_MOTION(i)%Ux
        Y(i)  = LAYER_MOTION(i)%Py + LAYER_MOTION(i)%Uy
        Z(i)  = LAYER_MOTION(i)%Pz + LAYER_MOTION(i)%Uz
        Vx(i) = LAYER_MOTION(i)%Vx
        Vy(i) = LAYER_MOTION(i)%Vy
        Vz(i) = LAYER_MOTION(i)%Vz
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
      Rx = Rx * (ONE / Rmag)
      Ry = Ry * (ONE / Rmag)
      Rz = Rz * (ONE / Rmag)
      Smag  = SQRT (Sx*Sx + Sy*Sy + Sz*Sz)
      Sx = Sx * (ONE / Smag)
      Sy = Sy * (ONE / Smag)
      Sz = Sz * (ONE / Smag)
!!
!! Define the unit vector T normal to the element.
!!
      Tx = Ry*Sz - Sy*Rz
      Ty = Rz*Sx - Sz*Rx
      Tz = Rx*Sy - Sx*Ry
      Tmag = SQRT (Tx*Tx + Ty*Ty + Tz*Tz)
!!
      IF (Tmag .EQ. 0.0) THEN
        WRITE (MSG1,'(I8)') LYMBQ%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'LYMBQ_DIVERGENCE_OPERATOR.001.00'//                     &
     &          MSGL//'Element Geometry Has An Undefined Normal.'//            &
     &          MSGL//'Element ID:'//MSG1                                      &
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
      LYMBQ%RES%Area = 0.5 * (R13*S24 - S13*R24)
      Ain = ONE / LYMBQ%RES%Area
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
      SUBROUTINE LYMBQ_MASS (LYMBQ,SecID,MatID,LAYER_NODE,LAYER_MOTION)
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
      TYPE (membq_type) :: LYMBQ
      TYPE (node_type),   DIMENSION(8) :: LAYER_NODE
      TYPE (motion_type), DIMENSION(8) :: LAYER_MOTION
!!
      INTEGER                 &
     &          SecID,        & ! Current element section index                
     &          MatID           ! Current element material index
!!
!! Compute one fourth total mass.
!!
      Density = MATERIAL(MatID)%PVAL(1)
      QMass = 0.25*Density*LYMBQ%PAR%Area*SECTION_2D(SecID)%Thickness
!!
!! Accumulate mass at each nodal point.
!!
      DO i = 1,4
        LAYER_NODE(LYMBQ%PAR%IX(i))%Mass =                                     &
     &    LAYER_NODE(LYMBQ%PAR%IX(i))%Mass + QMass
        MATERIAL(MatID)%Mass = MATERIAL(MatID)%Mass + QMass
        Px = LAYER_MOTION(i)%Px
        Py = LAYER_MOTION(i)%Py
        Pz = LAYER_MOTION(i)%Pz
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
      SUBROUTINE LYMBQ_STRESS_DIVERGENCE (LYMBQ,SecID,MatID)
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
      TYPE (membq_type) :: LYMBQ
!!
      INTEGER                 &
     &          SecID,        & ! Current element section index                
     &          MatID           ! Current element material index
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
        WRITE (MSG1,'(I8)') LYMBQ%PAR%EleID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'LYMBQ_STRESS_DIVERGENCE.001.00'//                       &
     &          MSGL//'Sound Speed Imaginary, C**2 Negative.'//                &
     &          MSGL//'M4EL (4-Node Membrane) Element ID:'//MSG1               &
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
      LYMBQ%RES%DTelt = Dx / (QC + SQRT (QC*QC + Cv*Cv))
!!
!! Artificial hourglass viscosity.
!!
      Qhg = Density * Dx * HG_Visc * Cv
      PGr = (LYMBQ%RES%Pr + Qhg*Gr) * Delta
      PGs = (LYMBQ%RES%Ps + Qhg*Gs) * Delta
      PGt = (LYMBQ%RES%Pt + Qhg*Gt) * Delta
!!
!! Compute current element thickness based on constant volume.
!!
      Thickness = SECTION_2D(SecID)%Thickness *                                &
     &          LYMBQ%PAR%Area / LYMBQ%RES%Area
!!
!! Divergence of the membrane stress resultants.
!!
      Vrr = Thickness * (LYMBQ%RES%Stress(1) + QP - Hr*PGr)
      Vsr = Thickness * (LYMBQ%RES%Stress(3)      - Hs*PGr)
      Vrs = Thickness * (LYMBQ%RES%Stress(3)      - Hr*PGs)
      Vss = Thickness * (LYMBQ%RES%Stress(2) + QP - Hs*PGs)
      Vrt = Thickness * (                         - Hr*PGt)
      Vst = Thickness * (                         - Hs*PGt)
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
        LYMBQ%RES%Xint(i) = Rx*Fr(i) + Sx*Fs(i)+ Tx*Ft(i)
        LYMBQ%RES%Yint(i) = Ry*Fr(i) + Sy*Fs(i)+ Ty*Ft(i)
        LYMBQ%RES%Zint(i) = Rz*Fr(i) + Sz*Fs(i)+ Tz*Ft(i)
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE LYMBQ_HOURGLASS_FORCES (LYMBQ,MatID)
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
      TYPE (membq_type) :: LYMBQ
!!
      INTEGER                                 &
     &          MatID           ! Current element material index
!!
      COMMON /MEMBX/                          &
     &          Br(4),Bs(4),                  & ! Gradient Operators           
     &          Drr,Dss,Drs,Wrs,              & ! In-plane stretching          
     &          Delta,                        & ! Generalized element size     
     &          dBeta,                        & ! Incremental rotation         
     &          Hr,Hs,Ht,Gr,Gs,Gt,            & ! Anti-hg gradients (4-node)   
     &          Rx,Ry,Rz,Sx,Sy,Sz,Tx,Ty,Tz      ! Element basis vectors
!!
      dWrs = LYMBQ%RES%DTnext * Wrs
      HG_Stiff = MATERIAL(MatID)%PVAL(5)
      Qmod = LYMBQ%RES%DTnext * HG_Stiff * SOUND_SPEED%RCS2
!!
      Qr = Qmod * Gr + dWrs * LYMBQ%RES%Ps
      Qs = Qmod * Gs - dWrs * LYMBQ%RES%Pr
      Qt = Qmod * Gt
!!
      LYMBQ%RES%Pr = LYMBQ%RES%Pr + Qr
      LYMBQ%RES%Ps = LYMBQ%RES%Ps + Qs
      LYMBQ%RES%Pt = LYMBQ%RES%Pt + Qt
!!
      RETURN
      END
