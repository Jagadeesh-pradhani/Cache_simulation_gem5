      SUBROUTINE SLIDING_INTERFACE_CONTACT
!!
!! Copyright (c) by KEY Associates; 20-DEC-1991 12:52:44
!!
!! Purpose: Monitor the position of independently meshed domains and when
!! overlap is detected introduce interaction forces reflecting the fact
!! that contact has occurred and acting to "eliminate" the overlap.
!!
!! This is the executive module that determines and processes "active"
!! sliding interfaces. It distinguishes between regular and single-surface
!! sliding interfaces, and between "master-slave" and symmetric processing.
!!
!! The algorithm has three parts: (1) proximity checking, (2) contact
!! monitoring, and (3) computing interaction forces. The primary data
!! structures are SLIDING_NODE(1:NUMSN), CONTACT_SURFACE(1:NUMCE), and
!! CONTACT_NODE(1:NUMCN)% The algorithm works by examining the sliding
!! nodal points from the list SLIDING_NODE to see if they have penetrated
!! the element facets listed in CONTACT_SURFACE which are currently op-
!! posing them. The list CONTACT_NODE contains a list of all nodes parti-
!! cipating in a sliding interface whether as a sliding node or as a node
!! defining a contact element of a contact surface. The technique is
!! based on "looking ahead" to see if contact will occur. If contact is
!! detected, restoring forces are introduced which will eliminate the
!! overlap.
!!
!! SLIDING_NODE(i)%Nsn = Nodal point Nsn, a node on a sliding surface.
!! SLIDING_NODE(i)%Mce = Contact-element Mce, the element with which
!!                      the sliding nodal point is currentlly in
!!                      "contact." If Mce is negative, the sliding
!!                      nodal point is "off" the edge of the sliding
!!                      surface "adjacent" to the contact-element Mce.
!!                      Mce points to an element 4-tuple in the list
!!                      CONTACT_SURFACE(Mce)%IX(1:4)
!!
!! CONTACT_SURFACE(Mce)%IX(1:4) is the collection of elements (hexahedral
!! surface facets, shell elements, and membrane elements) defined in terms of
!! 4-tuples from which one or more sliding interfaces is defined. Either
!! symmetric or "master-slave" sliding interfaces are supported with these
!! data structures.
!!
!! CONTACT_SURFACE(Mce)%IX(1) = I-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(2) = J-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(3) = K-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(4) = L-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%NX(1) = Surface element neighbor across side IJ.
!! CONTACT_SURFACE(Mce)%NX(2) = Surface element neighbor across side JK.
!! CONTACT_SURFACE(Mce)%NX(3) = Surface element neighbor across side KL.
!! CONTACT_SURFACE(Mce)%NX(4) = Surface element neighbor across side LI.
!!
      USE shared_common_data
      USE sliding_interface_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
      USE coord_
      USE indx_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Symmetric,                                                     &
     &          Slave_Master,                                                  &
     &          Master_Slave,                                                  &
     &          Contact_Test,                                                  &
     &          Interface_Cycle,                                               &
     &          TwoSurface_OneSided,  & ! Two Solid Bodies In Contact          &
     &          OneSurface_OneSided,  & ! One Solid Body Self Contact          &
     &          OneSurface_TwoSided     ! Shell Buckling Two-Sided Contact
      LOGICAL                                                                  &
     &          Active
      DATA                                                                     &
     &          Symmetric           /0/                                        &
     &          Slave_Master        /1/                                        &
     &          Master_Slave        /2/                                        &
     &          Interface_Cycle     /0/                                        &
     &          TwoSurface_OneSided /0/                                        &
     &          OneSurface_OneSided /1/                                        &
     &          OneSurface_TwoSided /2/
!!
!! Check to see if any sliding interfaces are active.
!!
      Active = .FALSE.

      DO Nsi = 1,NUMSI
!!
!! Check to see if this is an active interface.
!!
        Active = Active .OR.                                                   &
     &          (                                                              &
     &          TIMSIM%Total .GE. SLIDING_INTERFACE(Nsi)%Begin                 &
     &          .AND.                                                          &
     &          TIMSIM%Total .LE. SLIDING_INTERFACE(Nsi)%End                   &
     &          )

      ENDDO

      IF (.NOT.Active) RETURN
!!
!! Index sliding interface execution counter.
!!
      Interface_Cycle = Interface_Cycle + 1
!!
!! Define constants used in contact force calculations.
!!
      DTaver = 0.5D+0 * (TIMSIM%DTlast + TIMSIM%DTnext)
      Eta = ONE / (DTaver*TIMSIM%DTnext)
!!
!! CONTACT FLAGS FOR POSTPROCESSING.
!! Clear contact flags used for postprocessing.
!!
      DO N = 1,NUMNP
        NODE(N)%ICF = 0
      ENDDO
!!
!! PROJECTED VELOCITIES AND CONFIGURATION.
!! Compute projected velocites at time t(n+1/2).
!! Compute projected positions at time t(n+1).
!!
      DO Ncn = 1,NUMCN
        ID = CONTACT_NODE(Ncn)%NPID
        Vx = MOTION(ID)%Vx + DTaver * MOTION(ID)%Ax
        Vy = MOTION(ID)%Vy + DTaver * MOTION(ID)%Ay
        Vz = MOTION(ID)%Vz + DTaver * MOTION(ID)%Az
        CONTACT_NODE(Ncn)%Vx = Vx
        CONTACT_NODE(Ncn)%Vy = Vy
        CONTACT_NODE(Ncn)%Vz = Vz
        Px = MOTION(ID)%Px + MOTION(ID)%Ux + TIMSIM%DTnext * Vx
        Py = MOTION(ID)%Py + MOTION(ID)%Uy + TIMSIM%DTnext * Vy
        Pz = MOTION(ID)%Pz + MOTION(ID)%Uz + TIMSIM%DTnext * Vz
        CONTACT_NODE(Ncn)%Px = Px
        CONTACT_NODE(Ncn)%Py = Py
        CONTACT_NODE(Ncn)%Pz = Pz
      ENDDO
!!
!! CONTACT ELEMENT NORMAL VECTORS.
!! Define vectors normal to contact elements in the projected configuration.
!!
      DO Mce = 1,NUMCE
!!
!! Define vectors R and S across the contact-element diagonals.
!! Locate each contact element's current center point.
!!
        I1 = CONTACT_SURFACE(Mce)%IX(1)
        I2 = CONTACT_SURFACE(Mce)%IX(2)
        I3 = CONTACT_SURFACE(Mce)%IX(3)
        I4 = CONTACT_SURFACE(Mce)%IX(4)

        IF (I4 .EQ. 0) THEN
          R1 = CONTACT_NODE(I2)%Px - CONTACT_NODE(I1)%Px
          R2 = CONTACT_NODE(I2)%Py - CONTACT_NODE(I1)%Py
          R3 = CONTACT_NODE(I2)%Pz - CONTACT_NODE(I1)%Pz
          S1 = CONTACT_NODE(I3)%Px - CONTACT_NODE(I1)%Px
          S2 = CONTACT_NODE(I3)%Py - CONTACT_NODE(I1)%Py
          S3 = CONTACT_NODE(I3)%Pz - CONTACT_NODE(I1)%Pz

          CONTACT_SURFACE(Mce)%Xave = 0.3333333333D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px + CONTACT_NODE(I2)%Px +                    &
     &          CONTACT_NODE(I3)%Px                                            &
     &          )
          CONTACT_SURFACE(Mce)%Yave = 0.3333333333D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py + CONTACT_NODE(I2)%Py +                    &
     &          CONTACT_NODE(I3)%Py                                            &
     &          )
          CONTACT_SURFACE(Mce)%Zave = 0.3333333333D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz + CONTACT_NODE(I2)%Pz +                    &
     &          CONTACT_NODE(I3)%Pz                                            &
     &          )
        ELSE
          R1 = CONTACT_NODE(I3)%Px - CONTACT_NODE(I1)%Px
          R2 = CONTACT_NODE(I3)%Py - CONTACT_NODE(I1)%Py
          R3 = CONTACT_NODE(I3)%Pz - CONTACT_NODE(I1)%Pz
          S1 = CONTACT_NODE(I4)%Px - CONTACT_NODE(I2)%Px
          S2 = CONTACT_NODE(I4)%Py - CONTACT_NODE(I2)%Py
          S3 = CONTACT_NODE(I4)%Pz - CONTACT_NODE(I2)%Pz

          CONTACT_SURFACE(Mce)%Xave = 0.2500000000D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px + CONTACT_NODE(I2)%Px +                    &
     &          CONTACT_NODE(I3)%Px + CONTACT_NODE(I4)%Px                      &
     &          )
          CONTACT_SURFACE(Mce)%Yave = 0.2500000000D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py + CONTACT_NODE(I2)%Py +                    &
     &          CONTACT_NODE(I3)%Py + CONTACT_NODE(I4)%Py                      &
     &          )
          CONTACT_SURFACE(Mce)%Zave = 0.2500000000D+0 *                        &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz + CONTACT_NODE(I2)%Pz +                    &
     &          CONTACT_NODE(I3)%Pz + CONTACT_NODE(I4)%Pz                      &
     &          )
        ENDIF
!!
!! Define An(1:3) normal to the contact-element Mce from R x S. It is assumed
!! that the right-hand rule for the cross product will generate a vector
!! pointing out of the body bounded by the n-tuple. The magnitude of R x S
!! is twice the area of the triangle/quadrilateral defined by R and S.
!!
        Ax = R2*S3 - S2*R3
        Ay = R3*S1 - S3*R1
        Az = R1*S2 - S1*R2
        Area = 0.5D+0 * SQRT (Ax*Ax + Ay*Ay + Az*Az)

        CONTACT_SURFACE(Mce)%Area  = Area
        CONTACT_SURFACE(Mce)%An(1) = Ax * (0.5D+0 / Area)
        CONTACT_SURFACE(Mce)%An(2) = Ay * (0.5D+0 / Area)
        CONTACT_SURFACE(Mce)%An(3) = Az * (0.5D+0 / Area)

      ENDDO
!!
!! INDIVIDUAL SLIDING INTERFACE PROCESSING.
!! Scan all sliding interfaces and check for sliding nodes that will, in the
!! next time step, penetrate an opposing contact element.
!!
      DO Nsi = 1,NUMSI
!!
!! Check to see if this is an active interface.
!!
        Active = TIMSIM%Total .GE. SLIDING_INTERFACE(Nsi)%Begin                &
     &       .AND. TIMSIM%Total .LE. SLIDING_INTERFACE(Nsi)%End
!!
        IF (Active) THEN
!!
          CoF     = SLIDING_INTERFACE(Nsi)%CoF
          Capture = SLIDING_INTERFACE(Nsi)%Capture
          Factor  = SLIDING_INTERFACE(Nsi)%Factor
          Border  = SLIDING_INTERFACE(Nsi)%Border
          Isym    = SLIDING_INTERFACE(Nsi)%Isym
          Ityp    = SLIDING_INTERFACE(Nsi)%Type
          Contact_Test = 2  !  (touch test)
!!
!! Select Scale on the basis of a symmetric (Scale = 0.5) or a master-slave
!! (Scale = 1.0) interface. (0=Symmetric; 1=Slave-Master; 2=Master-Slave)
!!
          IF (Isym .EQ. 0) THEN
            Scale = 0.5D+0
          ELSE
            Scale = ONE
          ENDIF
!!
!! Check for Type of sliding interface contact.
!!
          IF (Ityp .EQ. TwoSurface_OneSided) THEN
!!
!! Loop over sliding nodal points on Side 1 of sliding interface Nsi.
!!
            IF (Isym.EQ.Symmetric .OR. Isym.EQ.Master_Slave) THEN
!!
              ISN1bgn = SLIDING_INTERFACE(Nsi)%ISN1bgn
              ISN1end = SLIDING_INTERFACE(Nsi)%ISN1end
              ICE2bgn = SLIDING_INTERFACE(Nsi)%ICE2bgn
              ICE2end = SLIDING_INTERFACE(Nsi)%ICE2end
!!
!! Scan all sliding interfaces and check for sliding nodes that have moved
!! to an adjacent contact element. Update pointers if movement has occurred.
!!
              CALL TWO_SURFACE_CONTACT_SEARCH                                  &
     &        (                                                                &
     &        INDX(1,1),INDX(1,2),INDX(1,3),INDX(1,4),INDX(1,5),               &
     &        INDX(1,6),INDX(1,7),INDX(1,8),COORD(1,1),COORD(1,2),             &
     &        COORD(1,3),Capture,Contact_Test,Border,ISN1bgn,ISN1end,          &
     &        ICE2bgn,ICE2end                                                  &
     &        )

              CALL TWO_SURFACE_INTERACTION                                     &
     &        (                                                                &
     &        Eta,Scale,Capture,Contact_Test,Factor,CoF,ISN1bgn,ISN1end        &
     &        )
            ENDIF
!!
!! Loop over sliding nodal points on Side 2 of sliding interface Nsi.
!!
            IF (Isym.EQ.Symmetric .OR. Isym.EQ.Slave_Master) THEN
!!
              ISN2bgn = SLIDING_INTERFACE(Nsi)%ISN2bgn
              ISN2end = SLIDING_INTERFACE(Nsi)%ISN2end
              ICE1bgn = SLIDING_INTERFACE(Nsi)%ICE1bgn
              ICE1end = SLIDING_INTERFACE(Nsi)%ICE1end
!!
!! Scan all sliding interfaces and check for sliding nodes that have moved
!! to an adjacent contact element. Update pointers if movement has occurred.
!!
              CALL TWO_SURFACE_CONTACT_SEARCH                                  &
     &        (                                                                &
     &        INDX(1,1),INDX(1,2),INDX(1,3),INDX(1,4),INDX(1,5),               &
     &        INDX(1,6),INDX(1,7),INDX(1,8),COORD(1,1),COORD(1,2),             &
     &        COORD(1,3),Capture,Contact_Test,Border,ISN2bgn,ISN2end,          &
     &        ICE1bgn,ICE1end                                                  &
     &        )

              CALL TWO_SURFACE_INTERACTION                                     &
     &        (                                                                &
     &        Eta,Scale,Capture,Contact_Test,Factor,CoF,ISN2bgn,ISN2end        &
     &        )
            ENDIF
!!
          ELSE IF (Ityp .EQ. OneSurface_OneSided) THEN
!!
!! Loop over sliding nodal points on Side 1 of sliding interface Nsi.
!!
            ISN1bgn = SLIDING_INTERFACE(Nsi)%ISN1bgn
            ISN1end = SLIDING_INTERFACE(Nsi)%ISN1end
            ICE2bgn = SLIDING_INTERFACE(Nsi)%ICE2bgn
            ICE2end = SLIDING_INTERFACE(Nsi)%ICE2end
!!
!! Scan sliding interface and check for sliding nodes that have moved to
!! an adjacent contact element. Update pointers if movement has occurred.
!!
            CALL SINGLE_SURFACE_CONTACT_SEARCH                                 &
     &      (                                                                  &
     &      INDX(1,1),INDX(1,2),INDX(1,3),INDX(1,4),INDX(1,5),                 &
     &      INDX(1,6),INDX(1,7),INDX(1,8),COORD(1,1),COORD(1,2),               &
     &      COORD(1,3),Capture,Border,ISN1bgn,ISN1end,ICE2bgn,ICE2end          &
     &      )

            CALL SINGLE_SURFACE_INTERACTION                                    &
     &      (                                                                  &
     &      Eta,Scale,Capture,Factor,CoF,ISN1bgn,ISN1end,ICE2bgn,ICE2end       &
     &      )
!!
          ELSE IF (Ityp .EQ. OneSurface_TwoSided) THEN

          ENDIF
!!
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE TWO_SURFACE_CONTACT_SEARCH                                    &
     &  (                                                                      &
     &  LISN,LMCE,INDX_X,INDX_Y,INDX_Z,LNDX_X,LNDX_Y,LNDX_Z,                   &
     &  XORD,YORD,ZORD,Capture,Contact_Test,Border,ISNbgn,ISNend,              &
     &  ICEbgn,ICEend                                                          &
     &  )
!!
!! Copyright (c) by KEY Associates;  5-DEC-1993 15:01:35.37
!!
!! Purpose: Search for the nodal points that are currently closest
!! to each contact element.
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(OUT) :: LISN(*)
      INTEGER,         INTENT(OUT) :: LMCE(*)
      INTEGER,         INTENT(OUT) :: INDX_X(*),INDX_Y(*),INDX_Z(*)
      INTEGER,         INTENT(OUT) :: LNDX_X(*),LNDX_Y(*),LNDX_Z(*)
      REAL(KIND(0D0)), INTENT(OUT) ::   XORD(*),  YORD(*),  ZORD(*)
      REAL(KIND(0D0)), INTENT(IN)  :: Capture
      INTEGER,         INTENT(IN)  :: Contact_Test
      REAL(KIND(0D0)), INTENT(IN)  :: Border
      INTEGER,         INTENT(IN)  :: ISNbgn,ISNend,ICEbgn,ICEend
!!
!! Local variables.      
      INTEGER :: IORDLOC
!!
!! Clear sliding nodes of past contact element ID's.
!!
      DO ISN = ISNbgn,ISNend
        SLIDING_NODE(ISN)%Mce = 0
      ENDDO
!!
!! Sort each sliding node coordinate direction from minimum
!! to maximum coordinate value.
!!
      L = 0
      DO ISN = ISNbgn,ISNend
        L = L + 1
        INDX_X(L) = L
        INDX_Y(L) = L
        INDX_Z(L) = L
        XORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Px
        YORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Py
        ZORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Pz
      ENDDO
      LENGTH = L

      CALL RVSORT (XORD,INDX_X,LENGTH,IO_UNIT%LELO)
      CALL RVSORT (YORD,INDX_Y,LENGTH,IO_UNIT%LELO)
      CALL RVSORT (ZORD,INDX_Z,LENGTH,IO_UNIT%LELO)
!!
!! Construct reverse maps for INDX_X, INDX_Y and INDX_Z.
!!
      DO L = 1,LENGTH
        LNDX_X(INDX_X(L)) = L
        LNDX_Y(INDX_Y(L)) = L
        LNDX_Z(INDX_Z(L)) = L
      ENDDO
!!
!! Loop over all contact elements. For each contact element,
!! find all sliding nodes in its "box."
!!
      L = 0
      DO 200 Mce = ICEbgn,ICEend

        I1 = CONTACT_SURFACE(Mce)%IX(1)
        I2 = CONTACT_SURFACE(Mce)%IX(2)
        I3 = CONTACT_SURFACE(Mce)%IX(3)
        I4 = CONTACT_SURFACE(Mce)%IX(4)

        IF (I4 .EQ. 0) THEN
          Xmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px                                            &
     &          )
          Xmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px                                            &
     &          )

          Ymin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py                                            &
     &          )
          Ymax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py                                            &
     &          )

          Zmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz                                            &
     &          )
          Zmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz                                            &
     &          )
        ELSE
          Xmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px, CONTACT_NODE(I4)%Px                       &
     &          )
          Xmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px, CONTACT_NODE(I4)%Px                       &
     &          )

          Ymin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py, CONTACT_NODE(I4)%Py                       &
     &          )
          Ymax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py, CONTACT_NODE(I4)%Py                       &
     &          )

          Zmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz, CONTACT_NODE(I4)%Pz                       &
     &          )
          Zmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz, CONTACT_NODE(I4)%Pz                       &
     &          )
        ENDIF
!!
!! Expand box containing the contact element by an increment based
!! on the capture distance.
!!
        Delta = Capture * MAX (Xmax-Xmin, Ymax-Ymin, Zmax-Zmin)

        Xmin = Xmin - Delta
        Ymin = Ymin - Delta
        Zmin = Zmin - Delta

        Xmax = Xmax + Delta
        Ymax = Ymax + Delta
        Zmax = Zmax + Delta

        Imin = IORDLOC(XORD,LENGTH,Xmin,'MIN')
        Jmin = IORDLOC(YORD,LENGTH,Ymin,'MIN')
        Kmin = IORDLOC(ZORD,LENGTH,Zmin,'MIN')

        Imax = IORDLOC(XORD,LENGTH,Xmax,'MAX')
        Jmax = IORDLOC(YORD,LENGTH,Ymax,'MAX')
        Kmax = IORDLOC(ZORD,LENGTH,Zmax,'MAX')
!!
!! Take shortest list and build common set, that is, the set of indicies
!! common to all three intervals.
!!
        Lmin = MIN ((Imax-Imin), (Jmax-Jmin), (Kmax-Kmin))
        IF (Lmin .EQ. (Imax-Imin)) THEN
          DO I = Imin,Imax
            IX = INDX_X(I)
            IF (Jmin .LE. LNDX_Y(IX) .AND. LNDX_Y(IX) .LE. Jmax) THEN
              IF (Kmin .LE. LNDX_Z(IX) .AND. LNDX_Z(IX) .LE. Kmax) THEN
                ISN = IX + (ISNbgn - 1)
!!                Nsn = SLIDING_NODE(ISN)%Nsn
                L = L + 1
                LISN(L) = ISN
                LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                IF (L .EQ. NUMCX) THEN
                  CALL TWO_SURFACE_CONTACT_TESTING                             &
     &                  (Contact_Test,Border,LISN,LMCE,NUMCX)
                  L = 0
                ENDIF

              ENDIF
            ENDIF
          ENDDO
        ELSE IF (Lmin .EQ. (Jmax-Jmin)) THEN
          DO J = Jmin,Jmax
            JX = INDX_Y(J)
            IF (Imin .LE. LNDX_X(JX) .AND. LNDX_X(JX) .LE. Imax) THEN
              IF (Kmin .LE. LNDX_Z(JX) .AND. LNDX_Z(JX) .LE. Kmax) THEN
                ISN = JX + (ISNbgn - 1)
!!                Nsn = SLIDING_NODE(ISN)%Nsn
                L = L + 1
                LISN(L) = ISN
                LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                IF (L .EQ. NUMCX) THEN
                  CALL TWO_SURFACE_CONTACT_TESTING                             &
     &                  (Contact_Test,Border,LISN,LMCE,NUMCX)
                  L = 0
                ENDIF

              ENDIF
            ENDIF
          ENDDO
        ELSE
          DO K = Kmin,Kmax
            KX = INDX_Z(K)
            IF (Imin .LE. LNDX_X(KX) .AND. LNDX_X(KX) .LE. Imax) THEN
              IF (Jmin .LE. LNDX_Y(KX) .AND. LNDX_Y(KX) .LE. Jmax) THEN
                ISN = KX + (ISNbgn - 1)
!!                Nsn = SLIDING_NODE(ISN)%Nsn
                L = L + 1
                LISN(L) = ISN
                LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                IF (L .EQ. NUMCX) THEN
                  CALL TWO_SURFACE_CONTACT_TESTING                             &
     &                  (Contact_Test,Border,LISN,LMCE,NUMCX)
                  L = 0
                ENDIF

              ENDIF
            ENDIF
          ENDDO
        ENDIF

 200    ENDDO
!!
!! Process the remaining "short stack" of candidate contacts.
!!
      IF (L .GT. 0) THEN
        LMAX = L
        CALL TWO_SURFACE_CONTACT_TESTING                                       &
     &          (Contact_Test,Border,LISN,LMCE,LMAX)
      ENDIF

      RETURN
      END
!!_
      INTEGER FUNCTION IORDLOC (VECTOR,N,X,MINMAX)
!!
!! Copyright (c) by KEY Associates;  4-DEC-1993 10:39:04.15
!!
!! Purpose: Locate X in the array VECTOR(1:N), that is, when MINMAX
!! equals 'MIN' find the first location N in the array VECTOR for
!! which VECTOR(N) .GE. X, and when MINMAX equals 'MAX' find the last
!! location N in the array VECTOR for which VECTOR(N) .LE. X. The
!! values in array VECTOR are assumed to be in ascending order. A
!! binary search is conducted.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: VECTOR(N) ! Ascending values to be searched
      INTEGER,         INTENT(IN) :: N         ! Dimension of the array VECTOR
      REAL(KIND(0D0)), INTENT(IN) :: X         ! Value to be located in VECTOR
      CHARACTER(3),    INTENT(IN) :: MINMAX    ! 'MIN' look for VECTOR(i) .GE. Xmin
                                               ! 'MAX' look for VECTOR(i) .LE. Xmax
!!
      IF (X .LE. VECTOR(1)) THEN
        IORDLOC = 1
        RETURN
      ELSE IF (X .GE. VECTOR(N)) THEN
        IORDLOC = N
        RETURN
      ENDIF

      Ione = 1
      Itwo = N

 100    Imid = (Ione + Itwo) / 2

      IF (Ione+1 .LT. Itwo) THEN
        IF (X .GT. VECTOR(Imid)) THEN
          Ione = Imid
          GO TO 100
        ELSE IF (X .LT. VECTOR(Imid)) THEN
          Itwo = Imid
          GO TO 100
        ELSE
          IORDLOC = Imid
          RETURN
        ENDIF
      ENDIF

      IF (MINMAX .EQ. 'MIN') THEN
        IORDLOC = Imid + 1
      ELSE ! IF (MINMAX .EQ. 'MAX') THEN
        IORDLOC = Imid
      ENDIF
!!
!! A less efficient linear scan, but can be used for checking.
!!
!!!     IORDLOC = 0
!!!     IF (MINMAX .EQ. 'MIN') THEN
!!!       DO i = 1,N,+1
!!!         IF (VECTOR(i) .GE. X) THEN
!!!           IORDLOC = i
!!!           RETURN
!!!         ENDIF
!!!       ENDDO
!!!       IORDLOC = N
!!!     ELSE ! IF (MINMAX .EQ. 'MAX') THEN
!!!       DO i = N,1,-1
!!!         IF (VECTOR(i) .LE. X) THEN
!!!           IORDLOC = i
!!!           RETURN
!!!         ENDIF
!!!       ENDDO
!!!       IORDLOC = 1
!!!     ENDIF

      RETURN
      END
!!_
      SUBROUTINE TWO_SURFACE_CONTACT_TESTING                                   &
     &          (Contact_Test,Border,LISN,LMCE,LMAX)
!!
!! Copyright (c) by KEY Associates; 12-JUN-1995 19:28:10.00
!!
!! Purpose: Retain old contact element or select new contact element
!! for each sliding node in the list LISN(1:Lmax).
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(IN) :: Contact_Test
      REAL(KIND(0D0)), INTENT(IN) :: Border
      INTEGER,         INTENT(IN) :: LISN(*)
      INTEGER,         INTENT(IN) :: LMCE(*)
      INTEGER,         INTENT(IN) :: LMAX
!!
!! Local variables.
      INTEGER, PARAMETER :: Tube_Test  = 1
      INTEGER, PARAMETER :: Touch_Test = 2
!!
!! Distinguish between the following two methods for deciding if
!! a sliding nodal point is opposite the candidate contact element
!! or off the edge of the candidate contact element:
!!
!!      (1) Tube Test. Using the average element normals of the
!!      candidate contact element and the adjacent contact surface
!!      elements, construct a "tube" and test for "in" or "out."
!!      (Does NOT use Border to extend element surface.)
!!
!!      (2) Touch Test. Using a detailed calculation of the
!!      time and place that the sliding nodal point and "plane"
!!      of the candidate contact element touch, determine if the
!!      contact is "on" or "off" the candidate contact element.
!!      (Does use Border to extend element surface.)
!!
      IF (Contact_Test .EQ. Tube_Test) THEN

        CALL TWO_SURFACE_TUBE_TEST (Border,LISN,LMCE,LMAX)

      ELSE IF (Contact_Test .EQ. Touch_Test) THEN

        CALL TWO_SURFACE_TOUCH_TEST (Border,LISN,LMCE,LMAX)

      ENDIF

      RETURN
      END
!!_
      SUBROUTINE TWO_SURFACE_TUBE_TEST (Border,LISN,LMCE,LMAX)
!!
!! Copyright (c) by KEY Associates; 12-JUN-1995 19:28:10.00
!!
!! Purpose: Retain old contact element or select new contact element
!! for each sliding node in the list LISN(1:Lmax).
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: Border
      INTEGER,         INTENT(IN) :: LISN(*)
      INTEGER,         INTENT(IN) :: LMCE(*)
      INTEGER,         INTENT(IN) :: LMAX
!!
!! Local variables.
      REAL(KIND(0D0)) :: A(3),B(3),C(3)
!!
!! Loop over sliding nodal points contained in the contact element's box.
!!
      DO L = 1,LMAX
        ISN = LISN(L)
        Mce = LMCE(L)
!!
!! Check for location of sliding nodal point wrt candidate contact
!! element. This check is based on checking if the sliding nodal
!! point is inside the contact elements "tube." The checks are
!! carried out one edge at a time.
!!
        Px = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Px
        Py = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Py
        Pz = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Pz

        I1 = CONTACT_SURFACE(Mce)%IX(1)
        I2 = CONTACT_SURFACE(Mce)%IX(2)
        I3 = CONTACT_SURFACE(Mce)%IX(3)
        I4 = CONTACT_SURFACE(Mce)%IX(4)

        A1 = CONTACT_SURFACE(Mce)%An(1)
        A2 = CONTACT_SURFACE(Mce)%An(2)
        A3 = CONTACT_SURFACE(Mce)%An(3)
!!
!! Distinguish between a triangular contact element (I4=0), and
!! a quadrilateral contact element (I4>0).
!!
        IF (I4 .EQ. 0) THEN
!!
!! Case (a): Triangular contact element.
!!
          Neighbor = CONTACT_SURFACE(Mce)%NX(3)
          IF (Neighbor .GT. 0) THEN
            A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
            A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
            A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
          ELSE
            A(1) = A1
            A(2) = A2
            A(3) = A3
          ENDIF

          B(1) = CONTACT_NODE(I1)%Px - CONTACT_NODE(I3)%Px
          B(2) = CONTACT_NODE(I1)%Py - CONTACT_NODE(I3)%Py
          B(3) = CONTACT_NODE(I1)%Pz - CONTACT_NODE(I3)%Pz

          C(1)  = Px - CONTACT_NODE(I3)%Px
          C(2)  = Py - CONTACT_NODE(I3)%Py
          C(3)  = Pz - CONTACT_NODE(I3)%Pz

          IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

            Neighbor = CONTACT_SURFACE(Mce)%NX(1)
            IF (Neighbor .GT. 0) THEN
              A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
              A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
              A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
            ELSE
              A(1) = A1
              A(2) = A2
              A(3) = A3
            ENDIF

            B(1) = CONTACT_NODE(I2)%Px - CONTACT_NODE(I1)%Px
            B(2) = CONTACT_NODE(I2)%Py - CONTACT_NODE(I1)%Py
            B(3) = CONTACT_NODE(I2)%Pz - CONTACT_NODE(I1)%Pz

            C(1)  = Px - CONTACT_NODE(I1)%Px
            C(2)  = Py - CONTACT_NODE(I1)%Py
            C(3)  = Pz - CONTACT_NODE(I1)%Pz

            IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

              Neighbor = CONTACT_SURFACE(Mce)%NX(2)
              IF (Neighbor .GT. 0) THEN
                A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
                A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
                A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
              ELSE
                A(1) = A1
                A(2) = A2
                A(3) = A3
              ENDIF

              B(1) = CONTACT_NODE(I3)%Px - CONTACT_NODE(I2)%Px
              B(2) = CONTACT_NODE(I3)%Py - CONTACT_NODE(I2)%Py
              B(3) = CONTACT_NODE(I3)%Pz - CONTACT_NODE(I2)%Pz

              C(1)  = Px - CONTACT_NODE(I2)%Px
              C(2)  = Py - CONTACT_NODE(I2)%Py
              C(3)  = Pz - CONTACT_NODE(I2)%Pz

              IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

                SLIDING_NODE(ISN)%Mce = Mce

              ENDIF
            ENDIF
          ENDIF

        ELSE
!!
!! Case (b): Quadrilateral contact element.
!!
          Neighbor = CONTACT_SURFACE(Mce)%NX(4)
          IF (Neighbor .GT. 0) THEN
            A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
            A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
            A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
          ELSE
            A(1) = A1
            A(2) = A2
            A(3) = A3
          ENDIF

          B(1) = CONTACT_NODE(I1)%Px - CONTACT_NODE(I4)%Px
          B(2) = CONTACT_NODE(I1)%Py - CONTACT_NODE(I4)%Py
          B(3) = CONTACT_NODE(I1)%Pz - CONTACT_NODE(I4)%Pz

          C(1)  = Px - CONTACT_NODE(I4)%Px
          C(2)  = Py - CONTACT_NODE(I4)%Py
          C(3)  = Pz - CONTACT_NODE(I4)%Pz

          IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

            Neighbor = CONTACT_SURFACE(Mce)%NX(1)
            IF (Neighbor .GT. 0) THEN
              A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
              A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
              A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
            ELSE
              A(1) = A1
              A(2) = A2
              A(3) = A3
            ENDIF

            B(1) = CONTACT_NODE(I2)%Px - CONTACT_NODE(I1)%Px
            B(2) = CONTACT_NODE(I2)%Py - CONTACT_NODE(I1)%Py
            B(3) = CONTACT_NODE(I2)%Pz - CONTACT_NODE(I1)%Pz

            C(1)  = Px - CONTACT_NODE(I1)%Px
            C(2)  = Py - CONTACT_NODE(I1)%Py
            C(3)  = Pz - CONTACT_NODE(I1)%Pz

            IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

              Neighbor = CONTACT_SURFACE(Mce)%NX(2)
              IF (Neighbor .GT. 0) THEN
                A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
                A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
                A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
              ELSE
                A(1) = A1
                A(2) = A2
                A(3) = A3
              ENDIF

              B(1) = CONTACT_NODE(I3)%Px - CONTACT_NODE(I2)%Px
              B(2) = CONTACT_NODE(I3)%Py - CONTACT_NODE(I2)%Py
              B(3) = CONTACT_NODE(I3)%Pz - CONTACT_NODE(I2)%Pz

              C(1)  = Px - CONTACT_NODE(I2)%Px
              C(2)  = Py - CONTACT_NODE(I2)%Py
              C(3)  = Pz - CONTACT_NODE(I2)%Pz

              IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

                Neighbor = CONTACT_SURFACE(Mce)%NX(3)
                IF (Neighbor .GT. 0) THEN
                  A(1) = A1 + CONTACT_SURFACE(Neighbor)%An(1)
                  A(2) = A2 + CONTACT_SURFACE(Neighbor)%An(2)
                  A(3) = A3 + CONTACT_SURFACE(Neighbor)%An(3)
                ELSE
                  A(1) = A1
                  A(2) = A2
                  A(3) = A3
                ENDIF

                B(1) = CONTACT_NODE(I4)%Px - CONTACT_NODE(I3)%Px
                B(2) = CONTACT_NODE(I4)%Py - CONTACT_NODE(I3)%Py
                B(3) = CONTACT_NODE(I4)%Pz - CONTACT_NODE(I3)%Pz

                C(1)  = Px - CONTACT_NODE(I3)%Px
                C(2)  = Py - CONTACT_NODE(I3)%Py
                C(3)  = Pz - CONTACT_NODE(I3)%Pz

                IF (Cross_Dot_Product(A,B,C) .GE. 0.0) THEN

                  SLIDING_NODE(ISN)%Mce = Mce

                ENDIF
              ENDIF
            ENDIF
          ENDIF

        ENDIF
      ENDDO

      RETURN
      END
!!_
      SUBROUTINE TWO_SURFACE_TOUCH_TEST (Border,LISN,LMCE,LMAX)
!!
!! Copyright (c) by KEY Associates; 15-DEC-1993 21:08:55.36
!!
!! Purpose: Retain old contact element or select new contact element
!! for each sliding node in the list LISN(1:Lmax).
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: Border
      INTEGER,         INTENT(IN) :: LISN(*)
      INTEGER,         INTENT(IN) :: LMCE(*)
      INTEGER,         INTENT(IN) :: LMAX
!!
!! Local variables.
      INTEGER            :: INTEXT
      INTEGER, PARAMETER :: INTERIOR = 1
      LOGICAL            :: INTERCEPT
!!
!! Loop over sliding nodal points contained in the contact element's box
!! and see if this contact element should replace the "old" contact
!! element opposing the sliding nodal point.
!!
      DO L = 1,LMAX
        ISN = LISN(L)
        Mce = LMCE(L)
!!
!! Check for location and time of intercept between sliding nodal point
!! and candidate contact element.
!!
        IF(INTERCEPT(ISN,Mce,Xcontact,Ycontact,Zcontact,                       &
     &      Tcontact,Border,Xdepth,Ydepth,Zdepth,INTEXT)) THEN
!!
!! Check to see if this is the first contact detected. If it is,
!! save contact state data.
!!
          IF (SLIDING_NODE(ISN)%Mce .EQ. 0) THEN

            SLIDING_NODE(ISN)%Mce = Mce
            SLIDING_NODE(ISN)%Tcontact = Tcontact
            SLIDING_NODE(ISN)%Xcontact = Xcontact
            SLIDING_NODE(ISN)%Ycontact = Ycontact
            SLIDING_NODE(ISN)%Zcontact = Zcontact
            SLIDING_NODE(ISN)%Xdepth   = Xdepth
            SLIDING_NODE(ISN)%Ydepth   = Ydepth
            SLIDING_NODE(ISN)%Zdepth   = Zdepth
            SLIDING_NODE(ISN)%INTEXT   = INTEXT
          ELSE
!!
!! Compare current candidate contact with previously detected contact.
!!
            IF (INTEXT .EQ. SLIDING_NODE(ISN)%INTEXT) THEN
!!
!! Case (a) Both contacts are exterior.
!! Case (b) Both contacts are interior.
!!
!! ...Take largest penetration depth.
!!
             DSQ_old = SLIDING_NODE(ISN)%Xdepth*SLIDING_NODE(ISN)%Xdepth       &
     &               + SLIDING_NODE(ISN)%Ydepth*SLIDING_NODE(ISN)%Ydepth       &
     &               + SLIDING_NODE(ISN)%Zdepth*SLIDING_NODE(ISN)%Zdepth

             DSQ_new = Xdepth*Xdepth                                           &
     &               + Ydepth*Ydepth                                           &
     &               + Zdepth*Zdepth

             IF (DSQ_new .GT. DSQ_old) THEN

               SLIDING_NODE(ISN)%Mce = Mce
               SLIDING_NODE(ISN)%Tcontact = Tcontact
               SLIDING_NODE(ISN)%Xcontact = Xcontact
               SLIDING_NODE(ISN)%Ycontact = Ycontact
               SLIDING_NODE(ISN)%Zcontact = Zcontact
               SLIDING_NODE(ISN)%Xdepth   = Xdepth
               SLIDING_NODE(ISN)%Ydepth   = Ydepth
               SLIDING_NODE(ISN)%Zdepth   = Zdepth
               SLIDING_NODE(ISN)%INTEXT   = INTEXT
             ENDIF
           ELSE
!!
!! Case (c) One contact is interior, one contact is exterior.
!!
!! ...Take interior contact.
!!
              IF (INTEXT .EQ. INTERIOR) THEN

                SLIDING_NODE(ISN)%Mce = Mce
                SLIDING_NODE(ISN)%Tcontact = Tcontact
                SLIDING_NODE(ISN)%Xcontact = Xcontact
                SLIDING_NODE(ISN)%Ycontact = Ycontact
                SLIDING_NODE(ISN)%Zcontact = Zcontact
                SLIDING_NODE(ISN)%Xdepth   = Xdepth
                SLIDING_NODE(ISN)%Ydepth   = Ydepth
                SLIDING_NODE(ISN)%Zdepth   = Zdepth
                SLIDING_NODE(ISN)%INTEXT   = INTEXT
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!!_
      SUBROUTINE TWO_SURFACE_INTERACTION                                       &
     &  (Eta,Scale,Capture,Contact_Test,Factor,CoF,ISNbgn,ISNend)
!!
!! Copyright (c) by KEY Associates; 19-DEC-1991 22:03:31
!!
!! Purpose: When overlap is detected introduce interaction forces reflect-
!! ing the fact that contact has occurred and acting to "eliminate" the
!! overlap.
!!
!! SLIDING_NODE(i)%Nsn = Nodal point Nsn, a node on a sliding surface.
!! SLIDING_NODE(i)%Mce = Contact-element Mce, the element with which
!!                      the sliding nodal point is currentlly in
!!                      "contact." If Mce is negative, the sliding
!!                      nodal point is "off" the edge of the sliding
!!                      surface "adjacent" to the contact-element Mce.
!!                      Mce points to an element 4-tuple in the list
!!                      CONTACT_SURFACE(Mce)%IX(1:4)
!!
!! CONTACT_SURFACE(Mce)%IX(1:4) is the collection of elements (hexahedral
!! surface facets, shell elements, and membrane elements) defined in terms of
!! 4-tuples from which one or more sliding interfaces is defined. Either
!! symmetric or "master-slave" sliding interfaces are supported with these
!! data structures.
!!
!! CONTACT_SURFACE(Mce)%IX(1) = I-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(2) = J-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(3) = K-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%IX(4) = L-th nodal point of surface element M.
!! CONTACT_SURFACE(Mce)%NX(1) = Surface element neighbor across side IJ.
!! CONTACT_SURFACE(Mce)%NX(2) = Surface element neighbor across side JK.
!! CONTACT_SURFACE(Mce)%NX(3) = Surface element neighbor across side KL.
!! CONTACT_SURFACE(Mce)%NX(4) = Surface element neighbor across side LI.
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: Eta
      REAL(KIND(0D0)), INTENT(IN) :: Scale
      REAL(KIND(0D0)), INTENT(IN) :: Capture
      INTEGER,         INTENT(IN) :: Contact_Test
      REAL(KIND(0D0)), INTENT(IN) :: Factor
      REAL(KIND(0D0)), INTENT(IN) :: CoF
      INTEGER,         INTENT(IN) :: ISNbgn,ISNend
!!
!! Local variables.
      INTEGER, PARAMETER :: Tube_Test  = 1
      INTEGER, PARAMETER :: Touch_Test = 2
      REAL(KIND(0D0))    :: Qn(3),Qm,Qmass,Qminv,Fsn,Fx,Fy,Fz
      LOGICAL            :: TRIANGLE

      REAL(KIND(0D0)), PARAMETER :: OneUpon3 = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: OneUpon4 = (1.0D+0 / 4.0D+0)
!!
      DTaver = 0.5D+0 * (TIMSIM%DTlast + TIMSIM%DTnext)
!!
!! Initialize mass of contact nodes and zero out force/acceleration values.
!!
      DO Ncn = 1,NUMCN
        CONTACT_NODE(Ncn)%Mass = NODE(CONTACT_NODE(Ncn)%NPID)%Mass
        CONTACT_NODE(Ncn)%Ax = 0.0
        CONTACT_NODE(Ncn)%Ay = 0.0
        CONTACT_NODE(Ncn)%Az = 0.0
      ENDDO
      DO Nsn = 1,NUMSN
        SLIDING_NODE(Nsn)%Force = 0.0
      ENDDO
!!
!! PASS 1. COMPUTE THE FORCE ACTING ON EACH SLIDING NODE TO PLACE IT
!! BACK ON THE CONTACT ELEMENT AS IF THE CONTACT ELEMENT IS UNAFFECTED
!! BY THE IMPACT.
!!
!! Loop over sliding nodal points on Side i of sliding interface Nsi. (See
!! calling sequence for i equal to 1 or 2.)
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes "off and way off the surface.")
!!
        IF (Mce .GT. 0) THEN
!!
!! Retrieve contact element normal vector.
!!
          Qn(1) = CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = CONTACT_SURFACE(Mce)%An(3)
!!
!! Compute depth of penetration based on data available from contact test.
!!
          IF (Contact_Test .EQ. Tube_Test) THEN
!!
!! Gather sliding nodal point Nsn coordinates.
!!
            Xsn = CONTACT_NODE(Nsn)%Px
            Ysn = CONTACT_NODE(Nsn)%Py
            Zsn = CONTACT_NODE(Nsn)%Pz
!!
!! Retrieve contact element average coordinates.
!!
            Xcn = CONTACT_SURFACE(Mce)%Xave
            Ycn = CONTACT_SURFACE(Mce)%Yave
            Zcn = CONTACT_SURFACE(Mce)%Zave
!!
!! Take dot-product between contact-element normal vector and
!! vector connecting "contact point" and position of sliding
!! nodal point at time t(n+1).
!!
            Delta = Qn(1)*(Xsn-Xcn) + Qn(2)*(Ysn-Ycn) + Qn(3)*(Zsn-Zcn)

          ELSE IF (Contact_Test .EQ. Touch_Test) THEN
!!
!! Retrieve sliding node "depth vector."
!!
            Xdp = SLIDING_NODE(i)%Xdepth
            Ydp = SLIDING_NODE(i)%Ydepth
            Zdp = SLIDING_NODE(i)%Zdepth
!!
!! Take dot-product between contact-element normal vector and
!! vector connecting "contact point" and position of sliding
!! nodal point at time t(n+1).
!!
            Delta = Qn(1)*Xdp + Qn(2)*Ydp + Qn(3)*Zdp

          ENDIF
!!
!! If no penetration has occurred, set contact-element index
!! Mce to zero.
!!
          IF (Delta .GE. 0.0) THEN
            SLIDING_NODE(i)%Mce = 0   ! No-penetration flag
!!
!! Compute interaction "force." (The negative sign cancels the
!! negative penetration distance to make Fsn a positive number.
!! Fsn acts in the direction of the contact element's outward
!! normal Qn.) (Factor is between 0.0 and 1.0; Factor less than
!! 1.0 means that it will take several time steps to correct
!! the penetration.)
!!
          ELSE
            Qmass = CONTACT_NODE(Nsn)%Mass
            Fsn = -(Scale * Eta) * Delta * Qmass
            SLIDING_NODE(i)%Force = Factor * Fsn
!!
!! End of if-test for penetration (negative penetration depth).
!!
          ENDIF
!!
!! End of if-test for sliding nodal points on the surface.
!!
        ENDIF
!!
!! End of do-loop over sliding nodal points on Side i of sliding interface Nsi.
!!
      ENDDO
!!
!! PASS 2. ON THE BASIS THAT (1) THE "IMPACTS" WILL MOVE THE
!! CONTACT ELEMENTS, AND (2) THE SLIDING NODES AND THE CONTACT
!! SURFACE WILL MOVE TOGETHER (BE TOGETHER AT THE END OF THE
!! TIME STEP), DISTRIBUTE THE SLIDING NODAL MASSES AND FORCES
!! TO THE NODES OF THE CONTACT SURFACE.
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes that have not penetrated a contact element.)
!!
        IF (Mce .GT. 0) THEN
!!
!! Retrieve contact element normal vector.
!!
          Qn(1) = CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = CONTACT_SURFACE(Mce)%An(3)
!!
!! Distinguish between a quadrilateral and triangular contact
!! element.
!!
          TRIANGLE = (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0)
          IF (TRIANGLE) THEN
            kmax = 3
            Cavr = OneUpon3
          ELSE
            kmax = 4
            Cavr = OneUpon4
          ENDIF
!!
!! Construct nodal point weighting to distribute impact forces
!! and sliding nodal point mass to contact-element nodes. Mass
!! weighting is used in order to produce the same acceleration
!! change at each node of the contact element.
!!
          Wght = 0.0
          DO k = 1,kmax
            NP = CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID
            Wght = Wght + NODE(NP)%Mass
          ENDDO
          DO k = 1,kmax
            NP = CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID
            SLIDING_NODE(i)%P(k) = NODE(NP)%Mass / Wght
          ENDDO
!!
!! Mark nodes (for plotting) that are participating in
!! the contact.
!!
          NODE(CONTACT_NODE(Nsn)%NPID)%ICF = 1
          DO k = 1,kmax
            NODE(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%ICF = 2
          ENDDO
!!
!! Distribute sliding node force and mass to contact element
!! nodes. Note that a negative sign in the force calculation
!! "reverses" the direction of Fsn to produce forces acting
!! on the contact element.
!!
          Fsn = SLIDING_NODE(i)%Force
          Qmass = NODE(CONTACT_NODE(Nsn)%NPID)%Mass
          DO k = 1,kmax
            Qm = Qmass * SLIDING_NODE(i)%P(k)
            Fx = -(Fsn * Qn(1)) * SLIDING_NODE(i)%P(k)
            Fy = -(Fsn * Qn(2)) * SLIDING_NODE(i)%P(k)
            Fz = -(Fsn * Qn(3)) * SLIDING_NODE(i)%P(k)
            Ncn = CONTACT_SURFACE(Mce)%IX(k)
            CONTACT_NODE(Ncn)%Mass = CONTACT_NODE(Ncn)%Mass + Qm
            CONTACT_NODE(Ncn)%Ax   = CONTACT_NODE(Ncn)%Ax   + Fx
            CONTACT_NODE(Ncn)%Ay   = CONTACT_NODE(Ncn)%Ay   + Fy
            CONTACT_NODE(Ncn)%Az   = CONTACT_NODE(Ncn)%Az   + Fz
          ENDDO
!!
!! End of if-test for sliding nodal points on the surface.
!!
        ENDIF
!!
!! End of do-loop over sliding nodal points on Side i of sliding interface Nsi.
!!
      ENDDO
!!
!! PASS 3. COMPUTE INCREMENTAL ACCELERATIONS OF THE CONTACT-SURFACE
!! NODES BASED ON ALL OF THE IMPACTS PRODUCED BY THE SLIDING NODES.
!!
!! Compute acceleration increments on contact-surface nodes. (Note:
!! There is no easy way to process just those contact nodes modified
!! by impact.)
!!
      DO Ncn = 1,NUMCN
        Qminv = ONE / CONTACT_NODE(Ncn)%Mass
        CONTACT_NODE(Ncn)%Ax = Qminv * CONTACT_NODE(Ncn)%Ax
        CONTACT_NODE(Ncn)%Ay = Qminv * CONTACT_NODE(Ncn)%Ay
        CONTACT_NODE(Ncn)%Az = Qminv * CONTACT_NODE(Ncn)%Az
!!
!! Modify forces acting on contact-surface nodes to reflect sliding
!! interface contact. (None of the data for the sliding nodes in
!! CONTACT_NODE(*) has been modified. That is, the acceleration
!! increments are zero.)
!!
        NP = CONTACT_NODE(Ncn)%NPID
        Qmass = NODE(NP)%Mass
        FORCE(NP)%Xext = FORCE(NP)%Xext + Qmass * CONTACT_NODE(Ncn)%Ax
        FORCE(NP)%Yext = FORCE(NP)%Yext + Qmass * CONTACT_NODE(Ncn)%Ay
        FORCE(NP)%Zext = FORCE(NP)%Zext + Qmass * CONTACT_NODE(Ncn)%Az
      ENDDO
!!
!! PASS 4. MODIFY THE FORCES ACTING ON THE SLIDING NODES TO REFLECT
!! THE MOVEMENT OF THE CONTACT SURFACE.
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes "off and way off the surface.")
!!
        IF (Mce .GT. 0) THEN

          TRIANGLE = (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0)
          IF (TRIANGLE) THEN
            kmax = 3
            Cavr = OneUpon3
          ELSE
            kmax = 4
            Cavr = OneUpon4
          ENDIF
!!
!! Collect contact-element Mce normal vector.
!!
          Qn(1) = CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = CONTACT_SURFACE(Mce)%An(3)
!!
!! Distinguish between velocity compatibility and position
!! compatibility at time n+1.
!!
!!!         IF (POSITION_COMPATIBILITY) THEN
!!
!! Compute the incremental acceleration from the contact surface
!! at the sliding nodal point. Modify the force acting on the
!! sliding nodal point to reflect the movement of the contact
!! surface and, thus, the movement of the sliding nodal point
!! so that it will be co-located with the contact surface at the
!! end of the time step.
!!
            Ax = 0.0
            Ay = 0.0
            Az = 0.0
            DO k = 1,kmax
              Ncn = CONTACT_SURFACE(Mce)%IX(k)
              Ax = Ax + CONTACT_NODE(Ncn)%Ax
              Ay = Ay + CONTACT_NODE(Ncn)%Ay
              Az = Az + CONTACT_NODE(Ncn)%Az
            ENDDO
            NP = CONTACT_NODE(Nsn)%NPID
            Fsn = SLIDING_NODE(i)%Force + Cavr *                               &
     &          NODE(NP)%Mass * (Qn(1)*Ax + Qn(2)*Ay + Qn(3)*Az)
!!
!!!         ELSE IF (VELOCITY_COMPATIBILITY) THEN
!!!
!!!         ENDIF
!!
!! Evaluate shear traction due to friction. The shear traction acts opposite
!! to the relative tangential velocity dVelx, dVely, dVelz. The magnitude of
!! the relative tangential velocity is dVtan. A nonzero coefficient of fric-
!! tion CoF is required to obtain a nonzero shear traction.
!!
          IF (CoF .GT. 0.0) THEN
!!
!! Compute components of the relative velocity normal and tangent to the
!! contact-element Mce.
!!
            Vxce = 0.0
            Vyce = 0.0
            Vzce = 0.0
            DO k = 1,kmax
              Ncn = CONTACT_SURFACE(Mce)%IX(k)
              Vxce = Vxce + CONTACT_NODE(Ncn)%Vx
              Vyce = Vyce + CONTACT_NODE(Ncn)%Vy
              Vzce = Vzce + CONTACT_NODE(Ncn)%Vz
            ENDDO
            dVelx = Cavr * Vxce - CONTACT_NODE(Nsn)%Vx
            dVely = Cavr * Vyce - CONTACT_NODE(Nsn)%Vy
            dVelz = Cavr * Vzce - CONTACT_NODE(Nsn)%Vz
            dVnrm = Qn(1)*dVelx + Qn(2)*dVely + Qn(3)*dVelz
            dVelx = dVelx - dVnrm * Qn(1)
            dVely = dVely - dVnrm * Qn(2)
            dVelz = dVelz - dVnrm * Qn(3)
            dVtan = SQRT (dVelx*dVelx + dVely*dVely + dVelz*dVelz)
            Area  = CONTACT_SURFACE(Mce)%Area
            Sigma = ABS (Fsn / Area)
            Ssn   = Shear_Traction (CoF,Sigma,dVtan) * Area
            IF (dVtan .LE. 1.0D-20) dVtan = ONE
!!
!! Put shear forces due to friction on contact-element Mce nodes.
!!
            DO k = 1,kmax
              Fx  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVelx)
              Fy  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVely)
              Fz  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVelz)
              Ncn = CONTACT_SURFACE(Mce)%IX(k)
              NP  = CONTACT_NODE(Ncn)%NPID
              FORCE(NP)%Xext = FORCE(NP)%Xext - Fx
              FORCE(NP)%Yext = FORCE(NP)%Yext - Fy
              FORCE(NP)%Zext = FORCE(NP)%Zext - Fz
            ENDDO
          ELSE
            Ssn   = 0.0
            dVtan = ONE
          ENDIF
!!
!! Compute components of interaction force for sliding node.
!!
          Fx = Fsn * Qn(1) + (Ssn / dVtan) * dVelx
          Fy = Fsn * Qn(2) + (Ssn / dVtan) * dVely
          Fz = Fsn * Qn(3) + (Ssn / dVtan) * dVelz
!!
          NP = CONTACT_NODE(Nsn)%NPID
          FORCE(NP)%Xext = FORCE(NP)%Xext + Fx
          FORCE(NP)%Yext = FORCE(NP)%Yext + Fy
          FORCE(NP)%Zext = FORCE(NP)%Zext + Fz
!!
        ENDIF
!!
      ENDDO
!!
      RETURN
      END
!!_
      REAL(KIND(0D0)) FUNCTION Cross_Dot_Product (A,B,C)
!!
!! Copyright (c) by KEY Associates; 27-APR-1991 22:17:28
!!
!! Purpose: Compute "A Cross B Dot C" from the vectors A, B, and C.
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: A(3),B(3),C(3)
!!
      Cross_Dot_Product = C(1)*(A(2)*B(3)-B(2)*A(3))                           &
     &                  + C(2)*(A(3)*B(1)-B(3)*A(1))                           &
     &                  + C(3)*(A(1)*B(2)-B(1)*A(2))
!!
      RETURN
      END
!!_
      REAL(KIND(0D0)) FUNCTION Shear_Traction (CoF,Sigma,dVel)
!!
!! Copyright (c) by KEY Associates; 17-JUN-1991 19:42:09
!!
!! Purpose: Implements the sliding interface friction model.
!! Current model: Mohr-Coulomb, velocity independent.
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: CoF   ! Coefficient of friction.
      REAL(KIND(0D0)), INTENT(IN) :: Sigma ! Normal stress, compression positive.
      REAL(KIND(0D0)), INTENT(IN) :: dVel  ! Relative tangential velocity magnitude.
!!
      Shear_Traction = CoF * MAX (Sigma,0.0D+0)
!!
      RETURN
      END
!!_
      SUBROUTINE SINGLE_SURFACE_CONTACT_SEARCH                                 &
     &  (                                                                      &
     &  LISN,LMCE,INDX_X,INDX_Y,INDX_Z,LNDX_X,LNDX_Y,LNDX_Z,                   &
     &  XORD,YORD,ZORD,Capture,Border,ISNbgn,ISNend,ICEbgn,ICEend              &
     &  )
!!
!! Copyright (c) by KEY Associates;  5-DEC-1993 15:01:35.37
!!
!! Purpose: Search for the nodal points that are currently closest
!! to each contact element.
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER,         INTENT(OUT) :: LISN(*)
      INTEGER,         INTENT(OUT) :: LMCE(*)
      INTEGER,         INTENT(OUT) :: INDX_X(*),INDX_Y(*),INDX_Z(*)
      INTEGER,         INTENT(OUT) :: LNDX_X(*),LNDX_Y(*),LNDX_Z(*)
      REAL(KIND(0D0)), INTENT(OUT) ::   XORD(*),  YORD(*),  ZORD(*)
      REAL(KIND(0D0)), INTENT(IN)  :: Capture
      REAL(KIND(0D0)), INTENT(IN)  :: Border
      INTEGER,         INTENT(IN)  :: ISNbgn,ISNend,ICEbgn,ICEend
!!
!! Local variables.      
      INTEGER :: IORDLOC
      LOGICAL :: Skip
!!
!! Clear sliding nodes of past contact element ID's.
!!
      DO ISN = ISNbgn,ISNend
        SLIDING_NODE(ISN)%Mce = 0
      ENDDO
!!
!! Sort each sliding node coordinate direction from minimum
!! to maximum coordinate value.
!!
      L = 0
      DO ISN = ISNbgn,ISNend
        L = L + 1
        INDX_X(L) = L
        INDX_Y(L) = L
        INDX_Z(L) = L
        XORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Px
        YORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Py
        ZORD(L) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Pz
      ENDDO
      LENGTH = L

      CALL RVSORT (XORD,INDX_X,LENGTH,IO_UNIT%LELO)
      CALL RVSORT (YORD,INDX_Y,LENGTH,IO_UNIT%LELO)
      CALL RVSORT (ZORD,INDX_Z,LENGTH,IO_UNIT%LELO)
!!
!! Construct reverse maps for INDX_X, INDX_Y and INDX_Z.
!!
      DO L = 1,LENGTH
        LNDX_X(INDX_X(L)) = L
        LNDX_Y(INDX_Y(L)) = L
        LNDX_Z(INDX_Z(L)) = L
      ENDDO
!!
!! Loop over all contact elements. For each contact element find
!! all sliding nodes in its "box." (Because this is a single surface
!! the nodal points defining the contact element will also be
!! found in the "box." They are excluded from consideration as a
!! sliding nodal point as they are found.)
!!
      L = 0
      DO 200 Mce = ICEbgn,ICEend

        I1 = CONTACT_SURFACE(Mce)%IX(1)
        I2 = CONTACT_SURFACE(Mce)%IX(2)
        I3 = CONTACT_SURFACE(Mce)%IX(3)
        I4 = CONTACT_SURFACE(Mce)%IX(4)

        IF (I4 .EQ. 0) THEN
          Xmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px                                            &
     &          )
          Xmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px                                            &
     &          )

          Ymin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py                                            &
     &          )
          Ymax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py                                            &
     &          )

          Zmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz                                            &
     &          )
          Zmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz                                            &
     &          )
        ELSE
          Xmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px, CONTACT_NODE(I4)%Px                       &
     &          )
          Xmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Px, CONTACT_NODE(I2)%Px,                      &
     &          CONTACT_NODE(I3)%Px, CONTACT_NODE(I4)%Px                       &
     &          )

          Ymin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py, CONTACT_NODE(I4)%Py                       &
     &          )
          Ymax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Py, CONTACT_NODE(I2)%Py,                      &
     &          CONTACT_NODE(I3)%Py, CONTACT_NODE(I4)%Py                       &
     &          )

          Zmin = MIN                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz, CONTACT_NODE(I4)%Pz                       &
     &          )
          Zmax = MAX                                                           &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Pz, CONTACT_NODE(I2)%Pz,                      &
     &          CONTACT_NODE(I3)%Pz, CONTACT_NODE(I4)%Pz                       &
     &          )
        ENDIF
!!
!! Expand box containing the contact element by an increment based
!! on the capture distance.
!!
        Delta = Capture * MAX (Xmax-Xmin, Ymax-Ymin, Zmax-Zmin)

        Xmin = Xmin - Delta
        Ymin = Ymin - Delta
        Zmin = Zmin - Delta

        Xmax = Xmax + Delta
        Ymax = Ymax + Delta
        Zmax = Zmax + Delta

        Imin = IORDLOC(XORD,LENGTH,Xmin,'MIN')
        Jmin = IORDLOC(YORD,LENGTH,Ymin,'MIN')
        Kmin = IORDLOC(ZORD,LENGTH,Zmin,'MIN')

        Imax = IORDLOC(XORD,LENGTH,Xmax,'MAX')
        Jmax = IORDLOC(YORD,LENGTH,Ymax,'MAX')
        Kmax = IORDLOC(ZORD,LENGTH,Zmax,'MAX')
!!
!! Take shortest list and build common set, that is, the set of indicies
!! common to all three intervals. Since this is a single-surface sliding
!! interface, skip those sliding nodal points that define the contact
!! element.
!!
        Lmin = MIN ((Imax-Imin), (Jmax-Jmin), (Kmax-Kmin))
        IF (Lmin .EQ. (Imax-Imin)) THEN
          DO I = Imin,Imax
            IX = INDX_X(I)
            IF (Jmin .LE. LNDX_Y(IX) .AND. LNDX_Y(IX) .LE. Jmax) THEN
              IF (Kmin .LE. LNDX_Z(IX) .AND. LNDX_Z(IX) .LE. Kmax) THEN
                ISN = IX + (ISNbgn - 1)
                Nsn = SLIDING_NODE(ISN)%Nsn
                Skip = Nsn.EQ.I1.OR.Nsn.EQ.I2.OR.Nsn.EQ.I3.OR.Nsn.EQ.I4
                IF (.NOT.Skip) THEN
                  L = L + 1
                  LISN(L) = ISN
                  LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                  IF (L .EQ. NUMCX) THEN
                    CALL SINGLE_SURFACE_CONTACT_TESTING                        &
     &                (Border,LISN,LMCE,NUMCX)
                    L = 0
                  ENDIF

                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ELSE IF (Lmin .EQ. (Jmax-Jmin)) THEN
          DO J = Jmin,Jmax
            JX = INDX_Y(J)
            IF (Imin .LE. LNDX_X(JX) .AND. LNDX_X(JX) .LE. Imax) THEN
              IF (Kmin .LE. LNDX_Z(JX) .AND. LNDX_Z(JX) .LE. Kmax) THEN
                ISN = JX + (ISNbgn - 1)
                Nsn = SLIDING_NODE(ISN)%Nsn
                Skip = Nsn.EQ.I1.OR.Nsn.EQ.I2.OR.Nsn.EQ.I3.OR.Nsn.EQ.I4
                IF (.NOT.Skip) THEN
                  L = L + 1
                  LISN(L) = ISN
                  LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                  IF (L .EQ. NUMCX) THEN
                    CALL SINGLE_SURFACE_CONTACT_TESTING                        &
     &                (Border,LISN,LMCE,NUMCX)
                    L = 0
                  ENDIF

                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ELSE
          DO K = Kmin,Kmax
            KX = INDX_Z(K)
            IF (Imin .LE. LNDX_X(KX) .AND. LNDX_X(KX) .LE. Imax) THEN
              IF (Jmin .LE. LNDX_Y(KX) .AND. LNDX_Y(KX) .LE. Jmax) THEN
                ISN = KX + (ISNbgn - 1)
                Nsn = SLIDING_NODE(ISN)%Nsn
                Skip = Nsn.EQ.I1.OR.Nsn.EQ.I2.OR.Nsn.EQ.I3.OR.Nsn.EQ.I4
                IF (.NOT.Skip) THEN
                  L = L + 1
                  LISN(L) = ISN
                  LMCE(L) = Mce
!!
!! Process list of candidate contacts when LISN and LMCE are full.
!!
                  IF (L .EQ. NUMCX) THEN
                    CALL SINGLE_SURFACE_CONTACT_TESTING                        &
     &                (Border,LISN,LMCE,NUMCX)
                    L = 0
                  ENDIF

                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF

 200    ENDDO
!!
!! Process the remaining "short stack" of candidate contacts.
!!
      IF (L .GT. 0) THEN
        LMAX = L
        CALL SINGLE_SURFACE_CONTACT_TESTING (Border,LISN,LMCE,LMAX)
      ENDIF

      RETURN
      END
!!_
      SUBROUTINE SINGLE_SURFACE_CONTACT_TESTING (Border,LISN,LMCE,LMAX)
!!
!! Copyright (c) by KEY Associates; 15-DEC-1993 21:08:55.36
!!
!! Purpose: Retain old contact element or select new contact element
!! for each sliding node in the list LISN(1:Lmax).
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: Border
      INTEGER,         INTENT(IN) :: LISN(*)
      INTEGER,         INTENT(IN) :: LMCE(*)
      INTEGER,         INTENT(IN) :: LMAX
!!
!! Local variables.
      INTEGER            :: INTEXT
      INTEGER, PARAMETER :: INTERIOR = 1
      LOGICAL            :: INTERCEPT
      LOGICAL            :: NEW_SUFACE
!!
!! Loop over sliding nodal points contained in the contact element's box
!! and see if this contact element should replace the "old" contact
!! element opposing the sliding nodal point.
!!
      DO L = 1,LMAX
        ISN = LISN(L)
        Mce = LMCE(L)
!!
!! Check for location and time of intercept between sliding nodal point
!! and candidate contact element.
!!
        IF (INTERCEPT(ISN,Mce,Xcontact,Ycontact,Zcontact,                      &
     &      Tcontact,Border,Xdepth,Ydepth,Zdepth,INTEXT)) THEN
!!
!! Check to see if this is the first contact detected. If it is,
!! save contact state data.
!!
          IF (SLIDING_NODE(ISN)%Mce .EQ. 0) THEN
!!
!! This is the first time this sliding nodal point has encounterd a
!! contact element.  Note from which side of the contact surface the
!! sliding node Nsn is approaching the contact element. (Use inner
!! product between surface outward normal and sliding nodal point
!! position wrt point of contact.) (The "direction of approach" is
!! introduced to allow for single-surface sliding interfaces.)
!!
!! Retrieve projected sliding nodal point coordinates.
!!
          Xsn  = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Px
          Ysn  = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Py
          Zsn  = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)%Pz
!!
!! Compute length of normal vector from contact element to sliding nodal point.
!!
          Pn = (Xsn-Xcontact) * CONTACT_SURFACE(Mce)%An(1) +                   &
     &         (Ysn-Ycontact) * CONTACT_SURFACE(Mce)%An(2) +                   &
     &         (Zsn-Zcontact) * CONTACT_SURFACE(Mce)%An(3)
!!
!! Since we are using projected geometry at time (n+1), the sliding nodal
!! point has "already" penetrated and Pn will be "negative" for an outward
!! pointing normal. Since SIGN is used to make An(1:3) an "outward" point-
!! ing normal, we must change the sign on Pn to get SLIDING_NODE(ISN)%Sign
!! correct.
!!
            SLIDING_NODE(ISN)%Mce      = Mce
            SLIDING_NODE(ISN)%Sign     = SIGN (ONE,-Pn)
            SLIDING_NODE(ISN)%Tcontact = Tcontact
            SLIDING_NODE(ISN)%Xcontact = Xcontact
            SLIDING_NODE(ISN)%Ycontact = Ycontact
            SLIDING_NODE(ISN)%Zcontact = Zcontact
            SLIDING_NODE(ISN)%Xdepth   = Xdepth
            SLIDING_NODE(ISN)%Ydepth   = Ydepth
            SLIDING_NODE(ISN)%Zdepth   = Zdepth
            SLIDING_NODE(ISN)%INTEXT   = INTEXT
          ELSE
!!
!! Compare current candidate contact with previously detected contact.
!!
            IF (INTEXT .EQ. SLIDING_NODE(ISN)%INTEXT) THEN
!!
!! Case (a) Both contacts are exterior.
!! Case (b) Both contacts are interior.
!!
!! ...Take largest penetration depth.
!!
             DSQ_old = SLIDING_NODE(ISN)%Xdepth*SLIDING_NODE(ISN)%Xdepth       &
     &               + SLIDING_NODE(ISN)%Ydepth*SLIDING_NODE(ISN)%Ydepth       &
     &               + SLIDING_NODE(ISN)%Zdepth*SLIDING_NODE(ISN)%Zdepth

             DSQ_new = Xdepth*Xdepth                                           &
     &               + Ydepth*Ydepth                                           &
     &               + Zdepth*Zdepth

             IF (DSQ_new .GT. DSQ_old) THEN

               SLIDING_NODE(ISN)%Mce = Mce
               SLIDING_NODE(ISN)%Tcontact = Tcontact
               SLIDING_NODE(ISN)%Xcontact = Xcontact
               SLIDING_NODE(ISN)%Ycontact = Ycontact
               SLIDING_NODE(ISN)%Zcontact = Zcontact
               SLIDING_NODE(ISN)%Xdepth   = Xdepth
               SLIDING_NODE(ISN)%Ydepth   = Ydepth
               SLIDING_NODE(ISN)%Zdepth   = Zdepth
               SLIDING_NODE(ISN)%INTEXT   = INTEXT
             ENDIF
           ELSE
!!
!! Case (c) One contact is interior, one contact is exterior.
!!
!! ...Take interior contact.
!!
              IF (INTEXT .EQ. INTERIOR) THEN

                SLIDING_NODE(ISN)%Mce = Mce
                SLIDING_NODE(ISN)%Tcontact = Tcontact
                SLIDING_NODE(ISN)%Xcontact = Xcontact
                SLIDING_NODE(ISN)%Ycontact = Ycontact
                SLIDING_NODE(ISN)%Zcontact = Zcontact
                SLIDING_NODE(ISN)%Xdepth   = Xdepth
                SLIDING_NODE(ISN)%Ydepth   = Ydepth
                SLIDING_NODE(ISN)%Zdepth   = Zdepth
                SLIDING_NODE(ISN)%INTEXT   = INTEXT
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
!!_
      SUBROUTINE SINGLE_SURFACE_INTERACTION                                    &
     &  (Eta,Scale,Capture,Factor,CoF,ISNbgn,ISNend,ICEbgn,ICEend)
!!
!! Copyright (c) by KEY Associates; 19-DEC-1991 22:03:31
!!
!! Purpose: When overlap is detected introduce interaction forces reflect-
!! ing the fact that contact has occurred and acting to "eliminate" the
!! overlap.
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_
      USE motion_
      USE force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      REAL(KIND(0D0)), INTENT(IN) :: Eta
      REAL(KIND(0D0)), INTENT(IN) :: Scale
      REAL(KIND(0D0)), INTENT(IN) :: Capture
      REAL(KIND(0D0)), INTENT(IN) :: Factor
      REAL(KIND(0D0)), INTENT(IN) :: CoF
      INTEGER,         INTENT(IN) :: ISNbgn,ISNend,ICEbgn,ICEend
!!
!! Local variables.
      REAL(KIND(0D0))    :: Qn(3),Qm,Qmass,Qminv,Fsn,Fx,Fy,Fz,Sgn
      LOGICAL            :: TRIANGLE

      REAL(KIND(0D0)), PARAMETER :: OneUpon3 = (1.0D+0 / 3.0D+0)
      REAL(KIND(0D0)), PARAMETER :: OneUpon4 = (1.0D+0 / 4.0D+0)
!!
!! Initialize mass of contact nodes and zero out force/acceleration values.
!!
      DO Ncn = 1,NUMCN
        CONTACT_NODE(Ncn)%Mass = NODE(CONTACT_NODE(Ncn)%NPID)%Mass
        CONTACT_NODE(Ncn)%Ax = 0.0
        CONTACT_NODE(Ncn)%Ay = 0.0
        CONTACT_NODE(Ncn)%Az = 0.0
      ENDDO
      DO i = 1,NUMSN
        SLIDING_NODE(i)%Force = 0.0
      ENDDO
!!
!! PASS 1. COMPUTE THE FORCE ACTING ON EACH SLIDING NODE TO PLACE IT
!! BACK ON THE CONTACT ELEMENT AS IF THE CONTACT ELEMENT IS UNAFFECTED
!! BY THE IMPACT.
!!
!! Loop over sliding nodal points on Side i of sliding interface Nsi. (See
!! calling sequence for i equal to 1 or 2.)
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
        Sgn = SLIDING_NODE(i)%Sign
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes "off and way off the surface.")
!!
        IF (Mce .GT. 0) THEN
!!
!! Retrieve contact element normal vector. The value of Sgn (+1 or -1)
!! assures that the normal vector to the contact element points
!! "outward" relative to the direction of contact.
!!
          Qn(1) = Sgn * CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = Sgn * CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = Sgn * CONTACT_SURFACE(Mce)%An(3)
!!
!! Retrieve sliding node "depth vector."
!!
          Xdp = SLIDING_NODE(i)%Xdepth
          Ydp = SLIDING_NODE(i)%Ydepth
          Zdp = SLIDING_NODE(i)%Zdepth
!!
!! Take dot-product between contact-element normal vector and
!! vector connecting "contact point" and position of sliding
!! nodal point at time t(n+1).
!!
          Delta = Qn(1)*Xdp + Qn(2)*Ydp + Qn(3)*Zdp
!!
!! If no penetration has occurred, set contact-element index
!! Mce to zero.
!!
          IF (Delta .GE. 0.0) THEN
            SLIDING_NODE(i)%Mce = 0   ! No-penetration flag
!!
!! Compute interaction "force." (The negative sign cancels the
!! negative penetration distance to make Fsn a positive number.
!! Fsn acts in the direction of the contact element's outward
!! normal Qn.) (Factor is between 0.0 and 1.0; Factor less than
!! 1.0 means that it will take several time steps to correct
!! the penetration.)
!!
          ELSE

!!!           write (*,'(a14,3i6,1pe15.7)')
!!!     &     ' Nsn,NPID,Mce,Delta:',Nsn,CONTACT_NODE(Nsn)%NPID,Mce,delta

            Qmass = CONTACT_NODE(Nsn)%Mass
            Fsn = -(Scale * Eta) * Delta * Qmass
            SLIDING_NODE(i)%Force = Factor * Fsn
!!
!! End of if-test for penetration (negative penetration depth).
!!
          ENDIF
!!
!! End of if-test for sliding nodal points on the surface.
!!
        ENDIF
!!
!! End of do-loop over sliding nodal points on Side i of sliding interface Nsi.
!!
      ENDDO
!!
!! PASS 2. ON THE BASIS THAT (1) THE "IMPACTS" WILL MOVE THE
!! CONTACT ELEMENTS, AND (2) THE SLIDING NODES AND THE CONTACT
!! SURFACE WILL MOVE TOGETHER (BE TOGETHER AT THE END OF THE
!! TIME STEP), DISTRIBUTE THE SLIDING NODAL MASSES AND FORCES
!! TO THE NODES OF THE CONTACT SURFACE.
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
        Sgn = SLIDING_NODE(i)%Sign
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes that have not penetrated a contact element.)
!!
        IF (Mce .GT. 0) THEN
!!
!! Retrieve contact element normal vector.
!!
          Qn(1) = Sgn * CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = Sgn * CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = Sgn * CONTACT_SURFACE(Mce)%An(3)
!!
!! Distinguish between a quadrilateral and triangular contact element.
!!
          TRIANGLE = (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0)
          IF (TRIANGLE) THEN
            kmax = 3
            Cavr = OneUpon3
          ELSE
            kmax = 4
            Cavr = OneUpon4
          ENDIF
!!
!! Construct nodal point weighting to distribute impact forces
!! and sliding nodal point mass to contact-element nodes. Mass
!! weighting is used in order to produce the same acceleration
!! change at each node of the contact element.
!!
          Wght = 0.0
          DO k = 1,kmax
            NP = CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID
            Wght = Wght + NODE(NP)%Mass
          ENDDO
          DO k = 1,kmax
            NP = CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID
            SLIDING_NODE(i)%P(k) = NODE(NP)%Mass / Wght
          ENDDO
!!
!! Mark nodes (for plotting) that are participating in
!! the contact.
!!
          NODE(CONTACT_NODE(Nsn)%NPID)%ICF = 1
          DO k = 1,kmax
            NODE(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%ICF = 2
          ENDDO
!!
!! Distribute sliding node force and mass to contact element
!! nodes. Note that a negative sign in the force calculation
!! "reverses" the direction of Fsn to produce forces acting
!! on the contact element.
!!
          Fsn = SLIDING_NODE(i)%Force
          Qmass = NODE(CONTACT_NODE(Nsn)%NPID)%Mass
          DO k = 1,kmax
            Qm = Qmass * SLIDING_NODE(i)%P(k)
            Fx = -(Fsn * Qn(1)) * SLIDING_NODE(i)%P(k)
            Fy = -(Fsn * Qn(2)) * SLIDING_NODE(i)%P(k)
            Fz = -(Fsn * Qn(3)) * SLIDING_NODE(i)%P(k)
            Ncn = CONTACT_SURFACE(Mce)%IX(k)
            CONTACT_NODE(Ncn)%Mass = CONTACT_NODE(Ncn)%Mass + Qm
            CONTACT_NODE(Ncn)%Ax   = CONTACT_NODE(Ncn)%Ax   + Fx
            CONTACT_NODE(Ncn)%Ay   = CONTACT_NODE(Ncn)%Ay   + Fy
            CONTACT_NODE(Ncn)%Az   = CONTACT_NODE(Ncn)%Az   + Fz
          ENDDO
!!
!! End of if-test for sliding nodal points on the surface.
!!
        ENDIF
!!
!! End of do-loop over sliding nodal points on Side i of sliding interface Nsi.
!!
      ENDDO
!!
!! PASS 3. COMPUTE INCREMENTAL ACCELERATIONS OF THE CONTACT-SURFACE
!! NODES BASED ON ALL OF THE IMPACTS PRODUCED BY THE SLIDING NODES.
!!
!! Compute acceleration increments on contact-surface nodes. (Note:
!! There is no easy way to process just those contact nodes modified
!! by impact.)
!!
      DO Ncn = 1,NUMCN
        Qminv = ONE / CONTACT_NODE(Ncn)%Mass
        CONTACT_NODE(Ncn)%Ax = Qminv * CONTACT_NODE(Ncn)%Ax
        CONTACT_NODE(Ncn)%Ay = Qminv * CONTACT_NODE(Ncn)%Ay
        CONTACT_NODE(Ncn)%Az = Qminv * CONTACT_NODE(Ncn)%Az
!!
!! Modify forces acting on contact-surface nodes to reflect sliding
!! interface contact. (None of the data for the sliding nodes in
!! CONTACT_NODE(*) has been modified. That is, the acceleration
!! increments are zero.)
!!
        NP = CONTACT_NODE(Ncn)%NPID
        Qmass = NODE(NP)%Mass
        FORCE(NP)%Xext = FORCE(NP)%Xext + Qmass * CONTACT_NODE(Ncn)%Ax
        FORCE(NP)%Yext = FORCE(NP)%Yext + Qmass * CONTACT_NODE(Ncn)%Ay
        FORCE(NP)%Zext = FORCE(NP)%Zext + Qmass * CONTACT_NODE(Ncn)%Az
      ENDDO
!!
!! PASS 4. MODIFY THE FORCES ACTING ON THE SLIDING NODES TO REFLECT
!! THE MOVEMENT OF THE CONTACT SURFACE.
!!
      DO i = ISNbgn,ISNend
        Nsn = SLIDING_NODE(i)%Nsn
        Mce = SLIDING_NODE(i)%Mce
        Sgn = SLIDING_NODE(i)%Sign
!!
!! Process only those sliding nodes on the contact surface. (Skip sliding
!! nodes "off and way off the surface.")
!!
        IF (Mce .GT. 0) THEN
!!
!! Collect contact-element Mce normal vector.
!!
          Qn(1) = Sgn * CONTACT_SURFACE(Mce)%An(1)
          Qn(2) = Sgn * CONTACT_SURFACE(Mce)%An(2)
          Qn(3) = Sgn * CONTACT_SURFACE(Mce)%An(3)
!!
!! Distinguish between a quadrilateral and triangular contact element.
!!
          TRIANGLE = (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0)
          IF (TRIANGLE) THEN
            kmax = 3
            Cavr = OneUpon3
          ELSE
            kmax = 4
            Cavr = OneUpon4
          ENDIF
!!
!! Compute the incremental acceleration from the contact surface
!! at the sliding nodal point. Modify the force acting on the
!! sliding nodal point to reflect the movement of the contact
!! surface and, thus, the movement of the sliding nodal point
!! so that it will be co-located with the contact surface at the
!! end of the time step.
!!
          Ax = 0.0
          Ay = 0.0
          Az = 0.0
          DO k = 1,kmax
            Ncn = CONTACT_SURFACE(Mce)%IX(k)
            Ax = Ax + CONTACT_NODE(Ncn)%Ax
            Ay = Ay + CONTACT_NODE(Ncn)%Ay
            Az = Az + CONTACT_NODE(Ncn)%Az
          ENDDO
          NP = CONTACT_NODE(Nsn)%NPID
          Fsn = SLIDING_NODE(i)%Force + Cavr *                                 &
     &          NODE(NP)%Mass * (Qn(1)*Ax + Qn(2)*Ay + Qn(3)*Az)
!!
!! Evaluate shear traction due to friction. The shear traction acts opposite
!! to the relative tangential velocity dVelx, dVely, dVelz. The magnitude of
!! the relative tangential velocity is dVtan. A nonzero coefficient of fric-
!! tion CoF is required to obtain a nonzero shear traction.
!!
          IF (CoF .GT. 0.0) THEN
!!
!! Compute components of the relative velocity normal and tangent to the
!! contact-element Mce.
!!
            Vxce = 0.0
            Vyce = 0.0
            Vzce = 0.0
            DO k = 1,kmax
              Ncn = CONTACT_SURFACE(Mce)%IX(k)
              Vxce = Vxce + CONTACT_NODE(Ncn)%Vx
              Vyce = Vyce + CONTACT_NODE(Ncn)%Vy
              Vzce = Vzce + CONTACT_NODE(Ncn)%Vz
            ENDDO
            dVelx = Cavr * Vxce - CONTACT_NODE(Nsn)%Vx
            dVely = Cavr * Vyce - CONTACT_NODE(Nsn)%Vy
            dVelz = Cavr * Vzce - CONTACT_NODE(Nsn)%Vz
            dVnrm = Qn(1)*dVelx + Qn(2)*dVely + Qn(3)*dVelz
            dVelx = dVelx - dVnrm * Qn(1)
            dVely = dVely - dVnrm * Qn(2)
            dVelz = dVelz - dVnrm * Qn(3)
            dVtan = SQRT (dVelx*dVelx + dVely*dVely + dVelz*dVelz)
            Area  = CONTACT_SURFACE(Mce)%Area
            Sigma = ABS (Fsn / Area)
            Ssn   = Shear_Traction (CoF,Sigma,dVtan) * Area
            IF (dVtan .LE. 1.0D-20) dVtan = ONE
!!
!! Put shear forces due to friction on contact-element Mce nodes.
!!
            DO k = 1,kmax
              Fx  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVelx)
              Fy  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVely)
              Fz  = SLIDING_NODE(i)%P(k) * ((Ssn / dVtan) * dVelz)
              Ncn = CONTACT_SURFACE(Mce)%IX(k)
              NP  = CONTACT_NODE(Ncn)%NPID
              FORCE(NP)%Xext = FORCE(NP)%Xext - Fx
              FORCE(NP)%Yext = FORCE(NP)%Yext - Fy
              FORCE(NP)%Zext = FORCE(NP)%Zext - Fz
            ENDDO
          ELSE
            Ssn   = 0.0
            dVtan = ONE
          ENDIF
!!
!! Compute components of interaction force for sliding node.
!!
          Fx = Fsn * Qn(1) + (Ssn / dVtan) * dVelx
          Fy = Fsn * Qn(2) + (Ssn / dVtan) * dVely
          Fz = Fsn * Qn(3) + (Ssn / dVtan) * dVelz
!!
          NP = CONTACT_NODE(Nsn)%NPID
          FORCE(NP)%Xext = FORCE(NP)%Xext + Fx
          FORCE(NP)%Yext = FORCE(NP)%Yext + Fy
          FORCE(NP)%Zext = FORCE(NP)%Zext + Fz
!!
        ENDIF
!!
      ENDDO
!!
      RETURN
      END
!!_
      LOGICAL FUNCTION INTERCEPT                                               &
     &          (ISN,Mce,Xcontact,Ycontact,Zcontact,Tcontact,Border,           &
     &          Xdepth,Ydepth,Zdepth,INTEXT)
!!
!! Copyright (c) by KEY Associates; 25-FEB-1995 21:16:02.00
!!
!! Purpose: Compute intersection time and location for the proposed
!! sliding-nodal-point/contact-element pair. (To process triangles
!! and to obtain a fast, easy intercept calculation, quadrilateral
!! surface facets are subdivided into four triangles.)
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! For reference:
!!      TYPE (contact_surface_type)             :: TRIANGLE_SURFACE
!!      TYPE (sliding_node_type)                :: TRIANGLE_SLINODE
!!      TYPE (contact_node_type), DIMENSION(4)  :: TRIANGLE_NODE
!!
      INTERCEPT = .FALSE.
!!
!! Distinguish between a triangular and a quadrilateral contact element.
!!
      I1 = CONTACT_SURFACE(Mce)%IX(1)
      I2 = CONTACT_SURFACE(Mce)%IX(2)
      I3 = CONTACT_SURFACE(Mce)%IX(3)
      I4 = CONTACT_SURFACE(Mce)%IX(4)
!!
!! If this is a triangular contact element (I4 = 0); look for intercept
!! directly.
!!
      IF (I4 .EQ. 0) THEN

        TRIANGLE_SLINODE = SLIDING_NODE(ISN)
        TRIANGLE_SLINODE%Nsn = 4
        TRIANGLE_SURFACE = CONTACT_SURFACE(Mce)
        TRIANGLE_SURFACE%IX(1) = 1
        TRIANGLE_SURFACE%IX(2) = 2
        TRIANGLE_SURFACE%IX(3) = 3
        TRIANGLE_NODE(1) = CONTACT_NODE(I1)
        TRIANGLE_NODE(2) = CONTACT_NODE(I2)
        TRIANGLE_NODE(3) = CONTACT_NODE(I3)
        TRIANGLE_NODE(4) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)

        CALL TRIANGULAR_SURFACE_INTERCEPT                                      &
     &          (Xi1,Xi2,Xcontact,Ycontact,Zcontact,Tcontact,                  &
     &          Xdepth,Ydepth,Zdepth)
!!
!! If an intercept has been found that occurs within the next time step
!! and occurs within the element, return .TRUE. and coordinate values.
!!
        IF (Tcontact .LT. TIMSIM%Total+TIMSIM%DTnext) THEN
          IF (Xi1+(Border-ONE).GE.Xi2 .AND. Xi1.LE.+Border                     &
     &        .AND. Xi2.GE.-Border) THEN

            INTERCEPT = .TRUE.
            IF (Xi1.GE.Xi2 .AND. Xi1.LE.+ONE .AND. Xi2.GE.-ONE) THEN
              INTEXT = 1
            ELSE
              INTEXT = 0
            ENDIF

          ENDIF
        ENDIF
!!
!! This is a quadrilateral contact element; subdivide into four
!! triangles and process triangles individually.
!!
      ELSE
        Vx = 0.25D+0 *                                                         &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Vx + CONTACT_NODE(I2)%Vx +                    &
     &          CONTACT_NODE(I3)%Vx + CONTACT_NODE(I4)%Vx                      &
     &          )
        Vy = 0.25D+0 *                                                         &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Vy + CONTACT_NODE(I2)%Vy +                    &
     &          CONTACT_NODE(I3)%Vy + CONTACT_NODE(I4)%Vy                      &
     &          )
        Vz = 0.25D+0 *                                                         &
     &          (                                                              &
     &          CONTACT_NODE(I1)%Vz + CONTACT_NODE(I2)%Vz +                    &
     &          CONTACT_NODE(I3)%Vz + CONTACT_NODE(I4)%Vz                      &
     &          )

        TRIANGLE_SLINODE = SLIDING_NODE(ISN)
        TRIANGLE_SLINODE%Nsn = 4
        TRIANGLE_SURFACE = CONTACT_SURFACE(Mce)
        TRIANGLE_SURFACE%IX(1) = 1
        TRIANGLE_SURFACE%IX(2) = 2
        TRIANGLE_SURFACE%IX(3) = 3
        TRIANGLE_NODE(3)%Vx = Vx
        TRIANGLE_NODE(3)%Vy = Vy
        TRIANGLE_NODE(3)%Vz = Vz
        TRIANGLE_NODE(3)%Px = CONTACT_SURFACE(Mce)%Xave
        TRIANGLE_NODE(3)%Py = CONTACT_SURFACE(Mce)%Yave
        TRIANGLE_NODE(3)%Pz = CONTACT_SURFACE(Mce)%Zave
        TRIANGLE_NODE(4) = CONTACT_NODE(SLIDING_NODE(ISN)%Nsn)

        INTEXT = -1
        K2 = I4
        DO k = 1,4
          K1 = K2
          K2 = CONTACT_SURFACE(Mce)%IX(k)
          TRIANGLE_NODE(1) = CONTACT_NODE(K1)
          TRIANGLE_NODE(2) = CONTACT_NODE(K2)

!!      write (IO_UNIT%LELO,*) ' Triangle EL #', k
!!      write (IO_UNIT%LELO,*) ' NODE 1:',TRIANGLE_NODE(1)%NPID
!!      write (IO_UNIT%LELO,*) ' NODE 2:',TRIANGLE_NODE(2)%NPID
!!      write (IO_UNIT%LELO,*) ' NODE 3:',TRIANGLE_NODE(3)%NPID
!!      write (IO_UNIT%LELO,*) ' SLIDE#:',TRIANGLE_NODE(4)%NPID

          CALL TRIANGULAR_SURFACE_INTERCEPT                                    &
     &          (Xi1,Xi2,Xcont,Ycont,Zcont,Tcontact,Xdep,Ydep,Zdep)
!!
!! If an intercept has been found that occurs within the next time step
!! and occurs within the element, return .TRUE. and the coordinate values.
!! (Note: The use of "Border" in the test for Xi1 only is based on the
!! usage in TRIANGULAR_SURFACE_INTERCEPT where the side "Xi1 = +1" is
!! the quadrilateral exterior edge.)
!!
          IF (Tcontact .LT. TIMSIM%Total+TIMSIM%DTnext) THEN
            IF (Xi1+(Border-ONE).GE.Xi2 .AND. Xi1.LE.+Border                   &
     &          .AND. Xi2.GE.-Border) THEN

              INTERCEPT = .TRUE.
              DSQ_new = Xdep*Xdep + Ydep*Ydep + Zdep*Zdep
              IF (Xi1 .GE. Xi2 .AND. Xi1.LE.+ONE .AND. Xi2.GE.-ONE) THEN
                INEX = 1
              ELSE
                INEX = 0
              ENDIF

              IF (INTEXT .EQ. -1) THEN
                Xcontact = Xcont
                Ycontact = Ycont
                Zcontact = Zcont
                Xdepth = Xdep
                Ydepth = Ydep
                Zdepth = Zdep
                INTEXT = INEX
                DSQ_old = DSQ_new
              ELSE
                IF (INEX .EQ. INTEXT) THEN
                  IF (DSQ_new .GT. DSQ_old) THEN
                    Xcontact = Xcont
                    Ycontact = Ycont
                    Zcontact = Zcont
                    Xdepth = Xdep
                    Ydepth = Ydep
                    Zdepth = Zdep
                    INTEXT = INEX
                    DSQ_old = DSQ_new
                  ENDIF
                ELSE IF (INEX .EQ. 1) THEN
                  Xcontact = Xcont
                  Ycontact = Ycont
                  Zcontact = Zcont
                  Xdepth = Xdep
                  Ydepth = Ydep
                  Zdepth = Zdep
                  INTEXT = INEX
                  DSQ_old = DSQ_new
                ENDIF
              ENDIF

            ENDIF
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END
!!_
      SUBROUTINE TRIANGULAR_SURFACE_INTERCEPT                                  &
     &          (Xi1,Xi2,Xcontact,Ycontact,Zcontact,Tcontact,                  &
     &          Xdepth,Ydepth,Zdepth)
!!
!! Copyright (c) by KEY Associates; 29-JAN-1994 12:48:30.92
!!
!! Purpose: Compute intersection time and location for the proposed
!! sliding nodal point-contact element pair.
!!
      USE shared_common_data
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! For reference:
!!      TYPE (contact_surface_type)             :: TRIANGLE_SURFACE
!!      TYPE (sliding_node_type)                :: TRIANGLE_SLINODE
!!      TYPE (contact_node_type), DIMENSION(4)  :: TRIANGLE_NODE
!!
!! Retrieve contact-element corner coordinates and sliding
!! nodal point coordinates at time (n+1), i.e., projected.
!! (Note: The use of "Border" in the test for Xi1 only in
!! INTERCEPT is based on the usage here where the side
!! "Xi1 = +1" is the quadrilateral exterior edge. That is,
!! if the definitions of I1,I2,I3 are changed, then the
!! tests using "Border" in the calling module will have
!! to be updated.)
!!
      I1 = TRIANGLE_SURFACE%IX(3)  !  For a quadrilateral, the center node.
      I2 = TRIANGLE_SURFACE%IX(1)  !  For a quadrilateral, an exterior node.
      I3 = TRIANGLE_SURFACE%IX(2)  !  For a quadrilateral, an exterior node.

!!      write (IO_UNIT%LELO,*) ' I1:',I1
!!      write (IO_UNIT%LELO,*) ' I2:',I2
!!      write (IO_UNIT%LELO,*) ' I3:',I3

      PxI = TRIANGLE_NODE(I1)%Px
      PyI = TRIANGLE_NODE(I1)%Py
      PzI = TRIANGLE_NODE(I1)%Pz
      VxI = TRIANGLE_NODE(I1)%Vx
      VyI = TRIANGLE_NODE(I1)%Vy
      VzI = TRIANGLE_NODE(I1)%Vz

      PxJ = TRIANGLE_NODE(I2)%Px
      PyJ = TRIANGLE_NODE(I2)%Py
      PzJ = TRIANGLE_NODE(I2)%Pz
      VxJ = TRIANGLE_NODE(I2)%Vx
      VyJ = TRIANGLE_NODE(I2)%Vy
      VzJ = TRIANGLE_NODE(I2)%Vz

      PxK = TRIANGLE_NODE(I3)%Px
      PyK = TRIANGLE_NODE(I3)%Py
      PzK = TRIANGLE_NODE(I3)%Pz
      VxK = TRIANGLE_NODE(I3)%Vx
      VyK = TRIANGLE_NODE(I3)%Vy
      VzK = TRIANGLE_NODE(I3)%Vz

      Nsn = TRIANGLE_SLINODE%Nsn

      PxS = TRIANGLE_NODE(Nsn)%Px
      PyS = TRIANGLE_NODE(Nsn)%Py
      PzS = TRIANGLE_NODE(Nsn)%Pz
      VxS = TRIANGLE_NODE(Nsn)%Vx
      VyS = TRIANGLE_NODE(Nsn)%Vy
      VzS = TRIANGLE_NODE(Nsn)%Vz
!!
!! The triangular contact surface is "augmented" with a
!! fourth node L located in the plane defined by I,J,K.
!! The corresponding skewed coordinated system (constant
!! skewness) allows us some easy algebra below. (We don't
!! have to construct anything for L or use L explicitly,
!! just use it conceptually, like "constant skewness.")
!!
!!         +1^Xi2
!!      L...|...K
!!       .   |  /|
!!       .   | / |
!!     -1.   |/  |+1
!!     --|---/---|->Xi1
!!       .  /| *=Sliding Nodal Point
!!       . / |   |
!!       ./  |   |=Contact Element
!!       I---|---J
!!         -1|
!!
!! Find point in I,J,K-plane closest to PxS,PyS,PzS The intercept
!! coordinates Xi1 and Xi2 are computed on the basis of a skewed
!! coordiate system (constant skewness).
!!
      Ax = 0.5D+00 * (PxI + PxK)
      Ay = 0.5D+00 * (PyI + PyK)
      Az = 0.5D+00 * (PzI + PzK)

      Bx = 0.5D+00 * (PxJ - PxI)
      By = 0.5D+00 * (PyJ - PyI)
      Bz = 0.5D+00 * (PzJ - PzI)

      Cx = 0.5D+00 * (PxK - PxJ)
      Cy = 0.5D+00 * (PyK - PyJ)
      Cz = 0.5D+00 * (PzK - PzJ)

      PB = (PxS-Ax)*Bx + (PyS-Ay)*By + (PzS-Az)*Bz
      PC = (PxS-Ax)*Cx + (PyS-Ay)*Cy + (PzS-Az)*Cz
      BB = Bx*Bx + By*By + Bz*Bz
      CC = Cx*Cx + Cy*Cy + Cz*Cz
      BC = Bx*Cx + By*Cy + Bz*Cz

      DET = BB*CC - BC*BC
!!
!! This is the point in the I,J,K-plane that is closest to PxS,PyS,PzS
!! PxP = Xi1*Bx + Xi2*Cx; PyP = Xi1*By + Xi2*Cy; PzP = Xi1*Bz + Xi2*Cz
!!
      IF (ABS(DET) .LT. 1.0D-20) THEN
        WRITE (MSG1,'(I8)') TRIANGLE_NODE(Nsn)%NPID  !  TRIANGLE_NODE(Nsn)%USID
        WRITE (MSG2,'(I8)') TIMSIM%Step
        WRITE (MSGF,'(1PE12.4)') TIMSIM%Total
        CALL USER_MESSAGE                                                      &
     &  (                                                                      &
     &  MSGL//'FATAL'//                                                        &
     &  MSGL//'TRIANGULAR_SURFACE_INTERCEPT.001.01'//                          &
     &  MSGL//'Contact Element With Singular Local Coordinate System.'//       &
     &  MSGL//'Element Opposite Sliding Node ID:'//MSG1//                      &
     &  MSGL//'Time Step:'//MSG2//'  Event Time:'//MSGF                        &
     &  )
      ELSE
        Xi1 = (CC*PB - BC*PC)/DET
        Xi2 = (BB*PC - BC*PB)/DET
      ENDIF

!!      if (contacted) then
!!        write(IO_UNIT%LELO,*) 'Xi1     :',Xi1
!!        write(IO_UNIT%LELO,*) 'Xi2     :',Xi2
!!      endif
!!
!! Compute contact point.
!!
      Xcontact = Ax + (Xi1*Bx + Xi2*Cx)
      Ycontact = Ay + (Xi1*By + Xi2*Cy)
      Zcontact = Az + (Xi1*Bz + Xi2*Cz)
!!
!! Estimate contact time by taking distance from the projected point in
!! the triangular contact facet to the actual location of the sliding
!! nodal point, all at time (n+1). (Delta_Time = Distance / Velocity)
!!
      Dx = PxS - Xcontact
      Dy = PyS - Ycontact
      Dz = PzS - Zcontact
!!
!! Cheaper...
!!      Vx = VxS - 0.333333333333333333D+0 * (VxI + VxJ + VxK)
!!      Vy = VyS - 0.333333333333333333D+0 * (VyI + VyJ + VyK)
!!      Vz = VzS - 0.333333333333333333D+0 * (VzI + VzJ + VzK)
!!
!! Better...
      Vx = VxS - 0.5D+0 * (VxI+VxK + Xi1*(VxJ-VxI) + Xi2*(VxK-VxI))
      Vy = VyS - 0.5D+0 * (VyI+VyK + Xi1*(VyJ-VyI) + Xi2*(VyK-VyI))
      Vz = VzS - 0.5D+0 * (VzI+VzK + Xi1*(VzJ-VzI) + Xi2*(VzK-VzI))

      DdotD = (Dx*Dx + Dy*Dy + Dz*Dz)
      DdotV = (Dx*Vx + Dy*Vy + Dz*Vz)

      IF (DdotD.LE.(1.0D-20) .OR. ABS(DdotV).LE.(1.0D-20)) THEN
        Tcontact = TIMSIM%Total + TIMSIM%DTnext

      ELSE
        Delta_t = DdotD / DdotV
        Tcontact = TIMSIM%Total + TIMSIM%DTnext - Delta_t

      ENDIF
!!
!! Record depth vector.
!!
      Xdepth = Dx
      Ydepth = Dy
      Zdepth = Dz
!!
      RETURN
      END
