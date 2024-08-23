      SUBROUTINE FORMATTED_PRINT_OF_INPUT
!!
!! Copyright (c) by KEY Associates; 13-DEC-1991 21:07:41
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
      USE property_
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: PAGE_NUMBER = 0  ! Local counter to number output pages
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered FORMATTED_PRINT_OF_INPUT.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open simulation data output (...sdo) file.
!!
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDO,                                        &
     &          FILE   = 'fmasdo',                                             &
     &          STATUS = 'UNKNOWN',                                            &
     &          FORM   = 'FORMATTED'                                           &
     &          )
!!
      CALL INTERCALATION (IO_UNIT%LSDO)
!!
!! List values for all counters.
!!
      CALL PRINT_COUNTERS (IO_UNIT%LSDO,'Number of ',NUMIF,PAGE_NUMBER)
      CALL PRINT_COUNTERS (IO_UNIT%LSDO,'Mesh with ',MSHIF,PAGE_NUMBER)
      CALL PRINT_COUNTERS (IO_UNIT%LSDO,'Restart w/',NRSIF,PAGE_NUMBER)
!!
!! Included files.
!!
      IF (NUMIF .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'INCLUDE_FILE (INCLUDE):'//                                          &
     &    '  File N, Line N,================================Input '//          &
     &    'File Name================================'

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,3X,A))')                               &
     &    (                                                                    &
     &    'INCLUDE_FILE (INCLUDE):',                                           &
     &    INCLUDE_FILE(n)%File_Number,                                         &
     &    INCLUDE_FILE(n)%Line_Number,                                         &
     &    INCLUDE_FILE(n)%Full_Name,                                           &
     &    n = 1,NUMIF                                                          &
     &    )

      ENDIF
!!
!! QA records.
!!
      IF (NUMQA .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'QA_RECORD (QAREC):'//                                               &
     &    '=============================User''s Input QA Text======='//        &
     &    '======================'

        WRITE (IO_UNIT%LSDO,'((1X,A,A))')                                      &
     &    (                                                                    &
     &    'QA_RECORD (QAREC):',                                                &
     &    QA_RECORD(n)%LINE(1:QA_RECORD(n)%LSIZE),                             &
     &    n = 1,NUMQA                                                          &
     &    )

      ENDIF
!!
!! Mass property calculation requests.
!!
      IF (NUMMP .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'MASSPROP (MASSPROP):'//                                             &
     &    '   MP ID, Set ID,Rot Typ,=====================X,Y,Z & Vx,'//        &
     &    'Vy,Vz for Rotational Properties====================='

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,6(1PE14.4)))')                         &
     &    (                                                                    &
     &    'MASSPROP (MASSPROP):',                                              &
     &    MASSPROP(n)%MPID,     & ! Mass property request ID                   &
     &    MASSPROP(n)%SetID,    & ! Nodal point set ID                         &
     &    MASSPROP(n)%Xnert,    & ! X-coordinate for Irot = 3                  &
     &    MASSPROP(n)%Ynert,    & ! Y-coordinate for Irot = 3                  &
     &    MASSPROP(n)%Znert,    & ! Z-coordinate for Irot = 3                  &
     &    MASSPROP(n)%Vxnert,   & ! X-velocity   for Irot = 3                  &
     &    MASSPROP(n)%Vynert,   & ! Y-velocity   for Irot = 3                  &
     &    MASSPROP(n)%Vznert,   & ! Z-velocity   for Irot = 3                  &
     &    n = 1,NUMMP                                                          &
     &    )

      ENDIF
!!
!! Material model input.
!!
      IF (NUMMT .GT. 0) THEN
        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)
        DO n = 1,NUMMT
          WRITE (IO_UNIT%LSDO,'(/1X,A)')                                       &
     &      'MATERIAL (MATERIAL):'//                                           &
     &      '   MatID,   Type,  Label...........................,  S'//        &
     &      'tate Var''s'
!!
          WRITE (IO_UNIT%LSDO,'(1X,A,2(I8),3X,A,I8/)')                         &
     &      n,                                                                 &
     &      MATERIAL(n)%MatID,                                                 &
     &      MATERIAL(n)%Type,                                                  &
     &      MATERIAL(n)%Label,                                                 &
     &      MATERIAL(n)%NSV
!!
          DO m = 1,5
            WRITE (IO_UNIT%LSDO,'((1X,2A,1PE14.4))')                           &
     &        'MATERIAL (MATERIAL): ',                                         &
     &        PROPERTY%NAME(m)//' =',                                          &
     &        MATERIAL(n)%PVAL(PROPERTY%Location(m))
          ENDDO
          DO m = 6,PROPERTY%Number_of_Entries
            IF (MATERIAL(n)%Set .EQ. PROPERTY%Set(m)) THEN
              WRITE (IO_UNIT%LSDO,'((1X,2A,1PE14.4))')                         &
     &          'MATERIAL (MATERIAL): ',                                       &
     &          PROPERTY%NAME(m)//' =',                                        &
     &          MATERIAL(n)%PVAL(PROPERTY%Location(m))
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!!
!! 1-D strain gauge.
!!
      IF (NUMG1 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'GAUGE1D (GAUGE1D):'//                                               &
     &    '   GauID,TB Flag,  El ID,Gau Loc,=============IX(1:4)===='//        &
     &    '========,IX Count'

        WRITE (IO_UNIT%LSDO,'((1X,A,9I8))')                                    &
     &    (                                                                    &
     &    'GAUGE1D (GAUGE1D):',                                                &
     &    GAUGE1D(n)%PAR%GauID, & ! Gauge ID (User defined value retained)       &
     &    GAUGE1D(n)%PAR%TBflg, & ! Truss/Beam Flag (0/1/2=no/truss/beam)        &
     &    GAUGE1D(n)%PAR%EleID, & ! Element ID (Reset: points to TRUSS/BEAMS)    &
     &    GAUGE1D(n)%PAR%GauLoc,& ! Cross section location (fiber strain)        &
     &    GAUGE1D(n)%PAR%IX(1), & ! Nodal points defining gauge section          &
     &    GAUGE1D(n)%PAR%IX(2), & ! Nodal points defining gauge section          &
     &    GAUGE1D(n)%PAR%IX(3), & ! Nodal points defining gauge section          &
     &    GAUGE1D(n)%PAR%IX(4), & ! Nodal points defining gauge section          &
     &    GAUGE1D(n)%PAR%NUMIX, & ! Actual number of nodes specified             &
     &    n = 1,NUMG1                                                          &
     &    )

      ENDIF
!!
!! 2-D strain gauge.
!!
      IF (NUMG2 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'GAUGE2D (GAUGE2D):'//                                               &
     &    '   GauID,MPQTFlg,  El ID,Gau Loc,========================'//        &
     &    '====IX(1:8)============================,IX Count'

        WRITE (IO_UNIT%LSDO,'((1X,A,13I8))')                                   &
     &    (                                                                    &
     &    'GAUGE2D (GAUGE2D):',                                                &
     &    GAUGE2D(n)%PAR%GauID,    & ! Gauge ID (User defined value retained)    &
     &    GAUGE2D(n)%PAR%MPQTflg,  & ! Membrane/Plate/Quad/Triangle flag         &
     &    GAUGE2D(n)%PAR%EleID,    & ! Element ID (Reset: points to MEMB*,PLAT*) &
     &    GAUGE2D(n)%PAR%GauLoc,   & ! Plate surface (-1,0,+1; fiber strain)     &
     &    GAUGE2D(n)%PAR%IX(1),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(2),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(3),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(4),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(5),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(6),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(7),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%IX(8),    & ! Nodal DOF's defining strain gauge.        &
     &    GAUGE2D(n)%PAR%NUMIX,    & ! Actual number of nodes specified
     &    n = 1,NUMG2                                                          &
     &    )

      ENDIF
!!
!! 3-D Strain gauge.
!!
      IF (NUMG3 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'GAUGE3D (GAUGE3D):'//                                               &
     &    '   GauID,HPT Flg,  El ID,============================IX(1'//        &
     &    ':8)============================,IX Count'

        WRITE (IO_UNIT%LSDO,'((1X,A,12I8))')                                   &
     &    (                                                                    &
     &    'GAUGE3D (GAUGE3D):',                                                &
     &    GAUGE3D(n)%PAR%GauID, & ! Gauge ID (User defined value retained)       &
     &    GAUGE3D(n)%PAR%HPTflg,& ! Hexa/Penta/Tetra flag (0/1/2/3=no/H/P/T)     &
     &    GAUGE3D(n)%PAR%EleID, & ! Element ID (Reset: pts. to HEXA,PENTA,TETRA) &
     &    GAUGE3D(n)%PAR%IX(1), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(2), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(3), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(4), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(5), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(6), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(7), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%IX(8), & ! Nodal points defining strain.                &
     &    GAUGE3D(n)%PAR%NUMIX, & ! Actual number of nodes specified             &
     &    n = 1,NUMG3                                                          &
     &    )

      ENDIF
!!
!! Solid element layering.
!!
      IF (NUMLU .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        DO n = 1,NUMLU
          WRITE (IO_UNIT%LSDO,'(/1X,A,I8,A,I8,3X,73(''=''))')                  &
     &      'LAY-UP (Lay-Up ID):',                                             &
     &      n,                                                                 &
     &      '(Coordinate System):',                                            &
     &      LAYERING(n)%Isys     ! Coordinate sytem, (0/n=global/special)

          WRITE (IO_UNIT%LSDO,'(/1X,A/)')                                      &
     &      ' Lay ID,LayTyp,Mat ID,=============Layer Thicknesses==='//        &
     &    '==========,=================Orthotropic Material Directio'//        &
     &    'ns================='

          WRITE (IO_UNIT%LSDO,'((1X,3I7,10(1PE11.3)))')                        &
     &      (                                                                  &
     &      n,                                                                 &
     &      LAYERING(n)%Ltype(i),& ! Layer type (0/n=solid/membr. section ID)    &
     &      LAYERING(n)%MatID(i),& ! Material ID (must match materials avail.)   &
     &      LAYERING(n)%H(1,i),  & ! Corner thickness, defines sub-elements      &
     &      LAYERING(n)%H(2,i),  & ! Corner thickness, defines sub-elements      &
     &      LAYERING(n)%H(3,i),  & ! Corner thickness, defines sub-elements      &
     &      LAYERING(n)%H(4,i),  & ! Corner thickness, defines sub-elements      &
     &      LAYERING(n)%Ax(i),   & ! Fiber orientation vector A, x-component     &
     &      LAYERING(n)%Ay(i),   & ! Fiber orientation vector A, y-component     &
     &      LAYERING(n)%Az(i),   & ! Fiber orientation vector A, z-component     &
     &      LAYERING(n)%Bx(i),   & ! Fiber orientation vector B, x-component     &
     &      LAYERING(n)%By(i),   & ! Fiber orientation vector B, y-component     &
     &      LAYERING(n)%Bz(i),   & ! Fiber orientation vector B, z-component     &
     &      i = 1,LAYERING(n)%Number_of_Layers                                 &
     &      )
        ENDDO

      ENDIF
!!
!! Plate cross section data.
!!
      IF (NUMS2 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SECTION_2D (PSECTION):'//                                           &
     &    ' SecID,IPts/Rule,Thickness,RefLoc,C-sys,================='//        &
     &    'Orthotropic Material Directions================='

        WRITE (IO_UNIT%LSDO,                                                   &
     &  '((1X,A,I6,2I5,1PE10.2,0PF7.3,I6,6(1PE11.3)))')                        &
     &    (                                                                    &
     &    'SECTION_2D (PSECTION):',                                            &
     &    n,                                                                   &
     &    SECTION_2D(n)%Ipts,     & ! Number of thickness integration points     &
     &    SECTION_2D(n)%Irule,    & ! Ftn: zeta coord's & weights f/ integration &
     &    SECTION_2D(n)%Thickness,& ! Thickness                                  &
     &    SECTION_2D(n)%RefLoc,   & ! Reference surface location (middle = 0.0)  &
     &    SECTION_2D(n)%Isys,     & ! Coordinate sytem, (0/n=global/special)     &
     &    SECTION_2D(n)%Ax,       & ! Fiber orientation vector A, x-component    &
     &    SECTION_2D(n)%Ay,       & ! Fiber orientation vector A, y-component    &
     &    SECTION_2D(n)%Az,       & ! Fiber orientation vector A, z-component    &
     &    SECTION_2D(n)%Bx,       & ! Fiber orientation vector B, x-component    &
     &    SECTION_2D(n)%By,       & ! Fiber orientation vector B, y-component    &
     &    SECTION_2D(n)%Bz,       & ! Fiber orientation vector B, z-component    &
     &    n = 1,NUMS2                                                          &
     &    )

      ENDIF
!!
!! Beam and truss cross section data.
!!
      IF (NUMS1 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SECTION_1D (BSECTION):'//                                           &
     &    '  Sec ID,Section,PropTyp,Ref Loc,====Width====,===Height='//        &
     &    '===,===T-Wall====,==T-Flange==='

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8,4(1PE14.4)))')                         &
     &    (                                                                    &
     &    'SECTION_1D (BSECTION):',                                            &
     &    n,                                                                   &
     &    SECTION_1D(n)%Section, & ! Beam cross section number (1,2,3,...)     &
     &    SECTION_1D(n)%Iprop,   & ! Geometric properties (0/1=dimension/inertia)&
     &    SECTION_1D(n)%NPLoc,   & ! Nodal point location (0,1,2,3,4)          &
     &    SECTION_1D(n)%Width,   & ! Cross section width                       &
     &    SECTION_1D(n)%Height,  & ! Cross section height                      &
     &    SECTION_1D(n)%Twall,   & ! Wall or web thickness                     &
     &    SECTION_1D(n)%Tflange, & ! Flange thickness (I-beam)                 &
     &    n = 1,NUMS1                                                          &
     &    )


        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SECTION_1D (BSECTION):'//                                           &
     &    '  Sec ID,Section,PropTyp,Ref Loc,====Area=====,==I(r**2)d'//        &
     &    'a==,==I(y**2)da==,==I(z**2)da=='

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8,4(1PE14.4)))')                         &
     &    (                                                                    &
     &    'SECTION_1D (BSECTION):',                                            &
     &    n,                                                                   &
     &    SECTION_1D(n)%Section, & ! Beam cross section number (1,2,3,...)     &
     &    SECTION_1D(n)%Iprop,   & ! Geometric properties (0/1=dimension/inertia)&
     &    SECTION_1D(n)%NPLoc,   & ! Nodal point location (0,1,2,3,4)          &
     &    SECTION_1D(n)%Area,    & ! Cross section area                        &
     &    SECTION_1D(n)%Br,      & ! I(r**2)da inertia about x-axis            &
     &    SECTION_1D(n)%By,      & ! I(y**2)da inertia about z-axis            &
     &    SECTION_1D(n)%Bz,      & ! I(z**2)da inertia about y-axis            &
     &    n = 1,NUMS1                                                          &
     &    )

      ENDIF
!!
!! Rigid body specification.
!!
      IF (NUMRB .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'RIGID_BODY (RBODY):   RB ID,Part ID,MassFlg,RBMASS ID'

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8))')                                    &
     &    (                                                                    &
     &    'RIGID_BODY (RBODY):',                                               &
     &    RIGID_BODY(n)%RBID,     & ! Rigid body ID (User defined value retained)&
     &    RIGID_BODY(n)%ParID,    & ! Part ID, Defines rigid body domain         &
     &    RIGID_BODY(n)%Prop,     & ! Computed mass properties not_used/used (0/1&
     &    RIGID_BODY(n)%CMID,     & ! Pts. to user spec'd added/subst'd masses
     &    n = 1,NUMRB                                                          &
     &    )

      ENDIF
!!
!! Rigid body mass prescription.
!!
      IF (NUMRM .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A,106(''=''))')                              &
     &    'RIGID_BODY_MASS (RBMASS):'

        WRITE (IO_UNIT%LSDO,'(/1X,A//(1X,I8,7(1PE14.4)))')                     &
     &    'RIGID_BODY_MASS (RBMASS):'//                                        &
     &    'RMass ID,====Mass=====,====I-xx=====,====I-yy=====,====I-'//        &
     &    'zz=====,====I-xy=====,====I-xz=====,====I-yz=====',                 &
     &    (                                                                    &
     &    'RIGID_BODY_MASS (RBMASS):',                                         &
     &    RIGID_BODY_MASS(n)%RMID,       & ! Rigid mass ID (User defined value ret
     &    RIGID_BODY_MASS(n)%Mass,       & ! Concentrated mass                   &
     &    RIGID_BODY_MASS(n)%B(1,1),     & ! Concentrated inertia                &
     &    RIGID_BODY_MASS(n)%B(2,2),     & ! Concentrated inertia                &
     &    RIGID_BODY_MASS(n)%B(3,3),     & ! Concentrated inertia                &
     &    RIGID_BODY_MASS(n)%B(1,2),     & ! Concentrated inertia                &
     &    RIGID_BODY_MASS(n)%B(1,3),     & ! Concentrated inertia                &
     &    RIGID_BODY_MASS(n)%B(2,3),     & ! Concentrated inertia
     &    n = 1,NUMRM                                                          &
     &    )


        WRITE (IO_UNIT%LSDO,'(/1X,A//(1X,2I7,9(1PE13.4)))')                    &
     &    'RIGID_BODY_MASS (RBMASS):'//                                        &
     &    'RMassID, NP ID,=====Px=====,=====Py=====,=====Pz=====,==='//        &
     &    '==Vx=====,=====Vy=====,=====Vz=====,=====Ox=====,=====Oy='//        &
     &    '====,=====Oz=====',                                                 &
     &    (                                                                    &
     &    'RIGID_BODY_MASS (RBMASS):',                                         &
     &    RIGID_BODY_MASS(n)%RMID,       & ! Rigid mass ID (User defined value re&
     &    NODE(RIGID_BODY_MASS(n)%NPID)%ID,      & ! Nodal point at which mass is&
     &    RIGID_BODY_MASS(n)%Pzero(1),   & ! Initial X,Y,Z position of C.M.      &
     &    RIGID_BODY_MASS(n)%Pzero(2),   & ! Initial X,Y,Z position of C.M.      &
     &    RIGID_BODY_MASS(n)%Pzero(3),   & ! Initial X,Y,Z position of C.M.
     &    RIGID_BODY_MASS(n)%Vel(1),     & ! Velocity of C.M.                    &
     &    RIGID_BODY_MASS(n)%Vel(2),     & ! Velocity of C.M.                    &
     &    RIGID_BODY_MASS(n)%Vel(3),     & ! Velocity of C.M.
     &    RIGID_BODY_MASS(n)%Omega(1),   & ! Angular velocity of C.M.            &
     &    RIGID_BODY_MASS(n)%Omega(2),   & ! Angular velocity of C.M.            &
     &    RIGID_BODY_MASS(n)%Omega(3),   & ! Angular velocity of C.M.
     &    n = 1,NUMRM                                                          &
     &    )

      ENDIF
!!
!! Nodal point concentrated masses.
!!
      IF (NUMCM .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NODAL_POINT_MASS (NPMASS):NMass ID,  NP ID,=====Mass===='

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,1(1PE14.4)))')                         &
     &    (                                                                    &
     &    'NODAL_POINT_MASS (NPMASS):',                                        &
     &    NODAL_POINT_MASS(n)%NMID,      & ! Nodal mass ID (User defined value re&
     &    NODE(NODAL_POINT_MASS(n)%NPID)%ID,     & ! Nodal point at which mass is
     &    NODAL_POINT_MASS(n)%Mass,      & ! Concentrated mass
     &    n = 1,NUMCM                                                          &
     &    )


        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NODAL_POINT_MASS (NPMASS):'//                                       &
     &    'NMass ID,  NP ID,====I-xx=====,====I-yy=====,====I-zz===='//        &
     &    '=,====I-xy=====,====I-xz=====,====I-yz====='

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,6(1PE14.4)))')                         &
     &    (                                                                    &
     &    'NODAL_POINT_MASS (NPMASS):',                                        &
     &    NODAL_POINT_MASS(n)%NMID,      & ! Nodal mass ID (User defined value re&
     &    NODAL_POINT_MASS(n)%NPID,      & ! Nodal point at which mass is concentr
     &    NODAL_POINT_MASS(n)%B(1,1),    & ! Concentrated inertia                &
     &    NODAL_POINT_MASS(n)%B(2,2),    & ! Concentrated inertia                &
     &    NODAL_POINT_MASS(n)%B(3,3),    & ! Concentrated inertia                &
     &    NODAL_POINT_MASS(n)%B(1,2),    & ! Concentrated inertia                &
     &    NODAL_POINT_MASS(n)%B(1,3),    & ! Concentrated inertia                &
     &    NODAL_POINT_MASS(n)%B(2,3),    & ! Concentrated inertia
     &    n = 1,NUMCM                                                          &
     &    )

      ENDIF
!!
!! Displacement boundary conditions.
!!
      IF (NUMDC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'DISPLACEMENT_BC (DISPBC):'//                                        &
     &    '   BC ID, Set ID,   Code,Hist ID,   Kavd,====Scale====,=='//        &
     &    '===Ax======,=====Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,5I8,4(1PE14.4)))')                         &
     &    (                                                                    &
     &    'DISPLACEMENT_BC (DISPBC):',                                         &
     &    DISPLACEMENT_BC(n)%DBCID,     & ! Displacement BC ID                   &
     &    DISPLACEMENT_BC(n)%SetID,     & ! Node set ID (-n/0/n=-NPID/all/Set ID)&
     &    DISPLACEMENT_BC(n)%Code,      & ! Constraint code                      &
     &    DISPLACEMENT_BC(n)%HstID,     & ! History ID (tabulated function ID)   &
     &    DISPLACEMENT_BC(n)%Kavd,      & ! Constrained kinematic variable (1/2/3&
     &    DISPLACEMENT_BC(n)%Scale,     & ! History function scale factor        &
     &    DISPLACEMENT_BC(n)%Ax,        & ! BC direction, x-component (Code=4,40,&
     &    DISPLACEMENT_BC(n)%Ay,        & ! BC direction, y-component  44,8,80,88&
     &    DISPLACEMENT_BC(n)%Az,        & ! BC direction, z-component            &
     &    n = 1,NUMDC                                                          &
     &    )

      ENDIF
!!
!! Tied boundary conditions.
!!
      IF (NUMTC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'TIED_BC (TIEDBC):'//                                                &
     &    '   BC ID, Set ID,   Code,====Fmax=====,=====Ax======,===='//        &
     &    '=Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,4(1PE14.4)))')                         &
     &    (                                                                    &
     &    'TIED_BC (TIEDBC):',                                                 &
     &    TIED_BC(n)%TBCID,     & ! Tied BC ID                                 &
     &    TIED_BC(n)%SetID,     & ! Node set ID (-n/0/n=-NPID/all/Set ID)      &
     &    TIED_BC(n)%Code,      & ! Constraint code                            &
     &    TIED_BC(n)%Fmax,      & ! Maximum constraint force b/ breaking       &
     &    TIED_BC(n)%Ax,        & ! BC direction, x-component (Code=4,40,      &
     &    TIED_BC(n)%Ay,        & ! BC direction, y-component  44,8,80,88)     &
     &    TIED_BC(n)%Az,        & ! BC direction, z-component                  &
     &    n = 1,NUMTC                                                          &
     &    )

      ENDIF
!!
!! Spot weld boundary conditions.
!!
      IF (NUMSW .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SPOT_WELD (SPOTWELD):'//                                            &
     &    '   SW ID,     PLACE     ,=====NP1=====,=====NP2=====,===='//        &
     &    '=NP3=====,====Fmax====='

        DO n = 1,NUMSW
          IF (SPOT_WELD(n)%PLACE .EQ. 'NODES') THEN
            WRITE (IO_UNIT%LSDO,'(1X,A,I8,4X,A,4X,3(I8,6X),1PE14.4)')          &
     &    'SPOT_WELD (SPOTWELD):',                                             &
     &    SPOT_WELD(n)%SWID,       & ! Spot weld ID
     &    SPOT_WELD(n)%PLACE,      & ! Keyword "NODE" or "POSITION"
     &    SPOT_WELD(n)%NPID(1),    & ! Node point ID 1
     &    SPOT_WELD(n)%NPID(2),    & ! Node point ID 2
     &    SPOT_WELD(n)%NPID(3),    & ! Node point ID 3
     &    SPOT_WELD(n)%Fmax          ! Fail constraint force, 0 => no-fail
          ENDIF
        ENDDO

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SPOT_WELD (SPOTWELD):'//                                            &
     &    '   SW ID,     PLACE     ,=====Px======,=====Py======,===='//        &
     &    '=Pz======,====Fmax====='

        DO n = 1,NUMSW
          IF (SPOT_WELD(n)%PLACE .EQ. 'POSITION') THEN
            WRITE (IO_UNIT%LSDO,'(1X,A,I8,4X,A,4X,4(1PE14.4))')                &
     &    'SPOT_WELD (SPOTWELD):',                                             &
     &    SPOT_WELD(n)%SWID,       & ! Spot weld ID
     &    SPOT_WELD(n)%PLACE,      & ! Keyword "NODE" or "POSITION"
     &    SPOT_WELD(n)%Xcm,        & ! X-coordinate of center of mass
     &    SPOT_WELD(n)%Ycm,        & ! Y-coordinate of center of mass
     &    SPOT_WELD(n)%Zcm,        & ! Z-coordinate of center of mass
     &    SPOT_WELD(n)%Fmax          ! Fail constraint force, 0 => no-fail
          ENDIF
        ENDDO

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SPOT_WELD (SPOTWELD):'//                                            &
     &    '   SW ID,     PLACE     ,====EleID====,0=tria/1=quad,===='//        &
     &    '=Xi1=====,====Xi2======'

        DO n = 1,NUMSW
          IF (SPOT_WELD(n)%PLACE .EQ. 'POSITION') THEN
            WRITE (IO_UNIT%LSDO,                                               &
     &      '(1X,A,I8,4X,A,4X,3(2(I8,6X),2(1PE14.4)/46X))')                    &
     &    'SPOT_WELD (SPOTWELD):',                                             &
     &    SPOT_WELD(n)%SWID,       & ! Spot weld ID
     &    SPOT_WELD(n)%PLACE,      & ! Keyword "NODE" or "POSITION"
     &    (                                                                    &
     &    SPOT_WELD(n)%EleID(i),   & ! Element ID's
     &    SPOT_WELD(n)%Type(i),    & ! 0/1 = triangle/quadrilateral
     &    SPOT_WELD(n)%Xi1(i),     & ! Isoparametric element location
     &    SPOT_WELD(n)%Xi2(i),     & ! Isoparametric element location
     &    i = 1,3                                                              &
     &    )
          ENDIF
        ENDDO

      ENDIF
!!
!! Wall boundary conditions.
!!
      IF (NUMWC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'RIGID_WALL_BC (WALLBC):'//                                          &
     &    '   BC ID, Set ID,===Px====,===Py====,===Pz====,===Wx====,'//        &
     &    '===Wy====,===Wz====,===Hx====,===Hy====,===Hz===='

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,9(1PE10.2)))')                         &
     &    (                                                                    &
     &    'RIGID_WALL_BC (WALLBC):',                                           &
     &    RIGID_WALL_BC(n)%RWID,     & ! Rigid wall ID                           &
     &    RIGID_WALL_BC(n)%SetID,    & ! Node set ID (-n/0/n=-NPID/all/Set ID)   &
     &    RIGID_WALL_BC(n)%P1(1),    & ! X,y,z of point on wall (center of mass) &
     &    RIGID_WALL_BC(n)%P1(2),    & ! X,y,z of point on wall (center of mass) &
     &    RIGID_WALL_BC(n)%P1(3),    & ! X,y,z of point on wall (center of mass) &
     &    RIGID_WALL_BC(n)%P2(1),    & ! X,y,z of point on wall (width direction)&
     &    RIGID_WALL_BC(n)%P2(2),    & ! X,y,z of point on wall (width direction)&
     &    RIGID_WALL_BC(n)%P2(3),    & ! X,y,z of point on wall (width direction)&
     &    RIGID_WALL_BC(n)%P3(1),    & ! X,y,z of point on wall (height direction&
     &    RIGID_WALL_BC(n)%P3(2),    & ! X,y,z of point on wall (height direction&
     &    RIGID_WALL_BC(n)%P3(3),    & ! X,y,z of point on wall (height direction)
     &    n = 1,NUMWC                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'RIGID_WALL_BC (WALLBC):'//                                          &
     &    '   BC ID,====Width====,===Height====,==Coef Fric=='

        WRITE (IO_UNIT%LSDO,'((1X,A,I8,3(1PE14.4)))')                          &
     &    (                                                                    &
     &    'RIGID_WALL_BC (WALLBC):',                                           &
     &    RIGID_WALL_BC(n)%RWID,                                               &
     &    RIGID_WALL_BC(n)%Width,                                              &
     &    RIGID_WALL_BC(n)%Height,                                             &
     &    RIGID_WALL_BC(n)%CoF,                                                &
     &    n = 1,NUMWC                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'RIGID_WALL_BC (WALLBC):'//                                          &
     &    '   BC ID,   Kode,====Mass=====,====I-xx=====,====I-yy===='//        &
     &    '=,====I-zz=====,====I-xy=====,====I-xz=====,====I-yz====='

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8,7(1PE14.4)))')                         &
     &    (                                                                    &
     &    'RIGID_WALL_BC (WALLBC):',                                           &
     &    RIGID_WALL_BC(n)%RWID,     & ! Rigid wall ID                         &
     &    RIGID_WALL_BC(n)%Kode,     & ! Constraint code (0/1/2=rigid/inertial/mo&
     &    RIGID_WALL_BC(n)%Mass,     & ! Concentrated mass at P1               &
     &    RIGID_WALL_BC(n)%B(1,1),   & ! Concentrated inertia tensor at P1     &
     &    RIGID_WALL_BC(n)%B(2,2),   & ! Concentrated inertia tensor at P1     &
     &    RIGID_WALL_BC(n)%B(3,3),   & ! Concentrated inertia tensor at P1     &
     &    RIGID_WALL_BC(n)%B(1,2),   & ! Concentrated inertia tensor at P1     &
     &    RIGID_WALL_BC(n)%B(1,3),   & ! Concentrated inertia tensor at P1     &
     &    RIGID_WALL_BC(n)%B(2,3),   & ! Concentrated inertia tensor at P1     &
     &    n = 1,NUMWC                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'RIGID_WALL_BC (WALLBC):'//                                          &
     &    '   BC ID,   Code,Hist ID,   Kavd,====Scale====,=====Ax==='//        &
     &    '===,=====Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8,4(1PE14.4)))')                         &
     &    (                                                                    &
     &    'RIGID_WALL_BC (WALLBC):',                                           &
     &    RIGID_WALL_BC(n)%RWID,     & ! Rigid wall ID                         &
     &    RIGID_WALL_BC(n)%Code,     & ! Constraint code on P1 (see Displacement &
     &    RIGID_WALL_BC(n)%HstID,    & ! History ID (tabulated function ID)      &
     &    RIGID_WALL_BC(n)%Kavd,     & ! Constrained kinematic variable (1/2/3=a/&
     &    RIGID_WALL_BC(n)%Scale,    & ! History function scale factor           &
     &    RIGID_WALL_BC(n)%Ax,       & ! BC direction, x-component (Code=4,40,   &
     &    RIGID_WALL_BC(n)%Ay,       & ! BC direction, y-component  44,8,80,88)  &
     &    RIGID_WALL_BC(n)%Az,       & ! BC direction, z-component               &
     &    n = 1,NUMWC                                                          &
     &    )

      ENDIF
!!
!! Midside nodal constraint data.
!!
      IF (NUMNC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'CONSTRAINED_NODE (NPCON1):'//                                       &
     &    '  Con ID,Ctrl NP,Left NP,Rght Np,=Tran''l Mass=,=Tran''l '//        &
     &    'Mass=,==Rot''l Mass=,==Rot''l Mass=,====Weight==='

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8,5(1PE14.4)))')                         &
     &    (                                                                    &
     &    'CONSTRAINED_NODE (NPCON1):',                                        &
     &    CONSTRAINED_NODE(n)%ID,       & ! Internal ID/pointer to next location &
     &    CONSTRAINED_NODE(n)%CNID,     & ! Constrained nodal point, ID          &
     &    CONSTRAINED_NODE(n)%NPID(1),  & ! Controlling nodal points, ID's       &
     &    CONSTRAINED_NODE(n)%NPID(2),  & ! Controlling nodal points, ID's       &
     &    CONSTRAINED_NODE(n)%TMass(1), & ! Effective mass at controlling nodes  &
     &    CONSTRAINED_NODE(n)%TMass(2), & ! Effective mass at controlling nodes  &
     &    CONSTRAINED_NODE(n)%RMass(1), & ! Effective mass at controlling nodes  &
     &    CONSTRAINED_NODE(n)%RMass(2), & ! Effective mass at controlling nodes  &
     &    CONSTRAINED_NODE(n)%Weight,   & ! Vcnid = V1*Weight + (1-Weight)*V2    &
     &    n = 1,NUMNC                                                          &
     &    )

      ENDIF
!!
!! Body forces.
!!
      IF (NUMBF .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'BODY_FORCE (BODYFORCE):'//                                          &
     &    '   BF ID, Set ID,Hist ID,====Scale====,===Gravity===,===='//        &
     &    '=Ax======,=====Ay======,=====Az======,====Delay===='

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,6(1PE14.4)))')                         &
     &    (                                                                    &
     &    'BODY_FORCE (BODYFORCE):',                                           &
     &    BODY_FORCE(n)%BFID,   & ! Pressure BC ID                               &
     &    BODY_FORCE(n)%SetID,  & ! Nodal point set ID (-n/n=-NPID/Set ID)       &
     &    BODY_FORCE(n)%HstID,  & ! History ID (tabulated function ID)           &
     &    BODY_FORCE(n)%Scale,  & ! History function scale factor                &
     &    BODY_FORCE(n)%Gravity,& ! Gravitational constant                       &
     &    BODY_FORCE(n)%Ax,     & ! Direction in which gravity acts, x-component &
     &    BODY_FORCE(n)%Ay,     & ! Direction in which gravity acts, y-component &
     &    BODY_FORCE(n)%Az,     & ! Direction in which gravity acts, z-component &
     &    BODY_FORCE(n)%Delay,  & ! Body force arrival time                      &
     &    n = 1,NUMBF                                                          &
     &    )

      ENDIF
!!
!! Pressure boundary conditions.
!!
      IF (NUMPC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PRESSURE_BC (PRESSBC):'//                                           &
     &    '   BC ID, Set ID,Hist ID,====Scale====,=====Pi======,===='//        &
     &    '=Pj======,=====Pk======,=====Pl======,====Delay===='

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,6(1PE14.4)))')                         &
     &    (                                                                    &
     &    'PRESSURE_BC (PRESSBC):',                                            &
     &    PRESSURE_BC(n)%PBCID, & ! Pressure BC ID                               &
     &    PRESSURE_BC(n)%SetID, & ! Segment set ID (-n/0/n=-SGID/all/Set ID)     &
     &    PRESSURE_BC(n)%HstID, & ! History ID (tabulated function ID)           &
     &    PRESSURE_BC(n)%Scale, & ! History function scale factor                &
     &    PRESSURE_BC(n)%PI,    & ! Pressure multiplier at node I/ on Seg Set    &
     &    PRESSURE_BC(n)%PJ,    & ! Pressure multiplier at node J                &
     &    PRESSURE_BC(n)%PK,    & ! Pressure multiplier at node K                &
     &    PRESSURE_BC(n)%PL,    & ! Pressure multiplier at node L                &
     &    PRESSURE_BC(n)%Delay, & ! Pressure arrival time                        &
     &    n = 1,NUMPC                                                          &
     &    )

      ENDIF
!!
!! Concentrated force boundary conditions.
!!
      IF (NUMFC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'FORCE_BC (FORCEBC):'//                                              &
     &    '   BC ID, Set ID,Hist ID,   Type, Follow,====Scale====,=='//        &
     &    '===Ax======,=====Ay======,=====Az======,====Delay===='

        WRITE (IO_UNIT%LSDO,'((1X,A,5I8,5(1PE14.4)))')                         &
     &    (                                                                    &
     &    'FORCE_BC (FORCEBC):',                                               &
     &    FORCE_BC(n)%CFID,     & ! Concentrated force ID                        &
     &    FORCE_BC(n)%SetID,    & ! Node set ID (-n/0/n=-NPID/all/Set ID)        &
     &    FORCE_BC(n)%HstID,    & ! History ID (force/torque = mag*f(t))         &
     &    FORCE_BC(n)%Type,     & ! Force or torque (0/1=force/torque)           &
     &    FORCE_BC(n)%Follow,   & ! Flag for follower force (0/1=no/yes)         &
     &    FORCE_BC(n)%Force,    & ! Concentrated force/torque magnitude          &
     &    FORCE_BC(n)%Cx,       & ! X-component of unit direction vector         &
     &    FORCE_BC(n)%Cy,       & ! Y-component of unit direction vector         &
     &    FORCE_BC(n)%Cz,       & ! Z-component of unit direction vector         &
     &    FORCE_BC(n)%Delay,    & ! Delay time to use with tabulated function    &
     &    n = 1,NUMFC                                                          &
     &    )

      ENDIF
!!
!! Spring boundary conditions.
!!
      IF (NUMSC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SPRING_BC (SPRINGBC):'//                                            &
     &    '   BC ID, Set ID, Mat ID,   Type, Follow,=====Ax======,=='//        &
     &    '===Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,5I8,3(1PE14.4)))')                         &
     &    (                                                                    &
     &    'SPRING_BC (SPRINGBC):',                                             &
     &    SPRING_BC(n)%SprID,   & ! Spring ID (User defined value retained)      &
     &    SPRING_BC(n)%SetID,   & ! Node set ID (-n/0/n=-NPID/all/Set ID)        &
     &    SPRING_BC(n)%MatID,   & ! Material ID  (Reset: points to MATERIAL)
     &    SPRING_BC(n)%Type,    & ! Axial or torsional (0/1=axial/torsional)     &
     &    SPRING_BC(n)%Follow,  & ! Flag for follower spring (0/1=U/Axis)        &
     &    SPRING_BC(n)%Axis(1), & ! Direction/hinge-axis orientation             &
     &    SPRING_BC(n)%Axis(2), & ! Direction/hinge-axis orientation             &
     &    SPRING_BC(n)%Axis(3), & ! Direction/hinge-axis orientation             &
     &    n = 1,NUMSC                                                          &
     &    )

      ENDIF
!!
!! Damper boundary conditions.
!!
      IF (NUMVC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'DAMPER_BC (DAMPERBC):'//                                            &
     &    '   BC ID, Set ID, Mat ID,  Type, Follow,=====Ax======,==='//        &
     &    '==Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,5I8,3(1PE14.4)))')                         &
     &    (                                                                    &
     &    'DAMPER_BC (DAMPERBC):',                                             &
     &    DAMPER_BC(n)%DprID,   & ! Damper ID (User defined value retained)      &
     &    DAMPER_BC(n)%SetID,   & ! Node set ID (-n/0/n=-NPID/all/Set ID)        &
     &    DAMPER_BC(n)%MatID,   & ! Material ID  (Reset: points to MATERIAL)
     &    DAMPER_BC(n)%Type,    & ! Axial or torsional (0/1=axial/torsional)     &
     &    DAMPER_BC(n)%Follow,  & ! Flag for follower Damper (0/1/2=U/V/Axis)    &
     &    DAMPER_BC(n)%Axis(1), & ! Direction/hinge-axis orientation             &
     &    DAMPER_BC(n)%Axis(2), & ! Direction/hinge-axis orientation             &
     &    DAMPER_BC(n)%Axis(3), & ! Direction/hinge-axis orientation             &
     &    n = 1,NUMVC                                                          &
     &    )

      ENDIF
!!
!! Periodic boundary conditions.
!!
      IF (NUMCC .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PERIODIC_BC (PERIODBC):'//                                          &
     &    '   BC ID,NPSetID,NPSetID,  Type,=====Ax======,=====Ay===='//        &
     &    '==,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,A8,3(1PE14.4)))')                      &
     &    (                                                                    &
     &    'PERIODIC_BC (PERIODBC):',                                           &
     &    PERIODIC_BC(n)%PerID,    & ! Periodic ID (User defined value retained)
     &    PERIODIC_BC(n)%S1ID,     & ! Node set ID (-n/n=-NPID/Set ID)
     &    PERIODIC_BC(n)%S2ID,     & ! Node set ID (-n/n=-NPID/Set ID)           &
     &    PERIODIC_BC(n)%Type,     & ! LINEAR or CYCLIC repeat                   &
     &    PERIODIC_BC(n)%Axis(1),  & ! Rotation axis orientation                 &
     &    PERIODIC_BC(n)%Axis(2),  & ! Rotation axis orientation                 &
     &    PERIODIC_BC(n)%Axis(3),  & ! Rotation axis orientation
     &    n = 1,NUMCC                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PERIODIC_BC (PERIODBC):'//                                          &
     &    '   BC ID,NPSetID,NPSetID,  Type,=====Px======,=====Py===='//        &
     &    '==,=====Pz======'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,A8,3(1PE14.4)))')                      &
     &    (                                                                    &
     &    'PERIODIC_BC (PERIODBC):',                                           &
     &    PERIODIC_BC(n)%PerID,    & ! Periodic ID (User defined value retained)
     &    PERIODIC_BC(n)%S1ID,     & ! Node set ID (-n/n=-NPID/Set ID)
     &    PERIODIC_BC(n)%S2ID,     & ! Node set ID (-n/n=-NPID/Set ID)           &
     &    PERIODIC_BC(n)%Type,     & ! LINEAR or CYCLIC repeat
     &    PERIODIC_BC(n)%Origin(1),& ! Rotation axis location                    &
     &    PERIODIC_BC(n)%Origin(2),& ! Rotation axis location                    &
     &    PERIODIC_BC(n)%Origin(3),& ! Rotation axis location
     &    n = 1,NUMCC                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PERIODIC_BC (PERIODBC):'//                                          &
     &    '   BC ID,NPSetID,NPSetID,  Type,====Theta====,===Advance==='

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,A8,2(1PE14.4)))')                      &
     &    (                                                                    &
     &    'PERIODIC_BC (PERIODBC):',                                           &
     &    PERIODIC_BC(n)%PerID,    & ! Periodic ID (User defined value retained)
     &    PERIODIC_BC(n)%S1ID,     & ! Node set ID (-n/n=-NPID/Set ID)
     &    PERIODIC_BC(n)%S2ID,     & ! Node set ID (-n/n=-NPID/Set ID)           &
     &    PERIODIC_BC(n)%Type,     & ! LINEAR or CYCLIC repeat
     &    PERIODIC_BC(n)%Theta,    & ! Rotation angle, degrees                   &
     &    PERIODIC_BC(n)%Advance,  & ! Helical advance per degree
     &    n = 1,NUMCC                                                          &
     &    )

      ENDIF
!!
!! Tabulated functions.
!!
      IF (NUMTF .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &   'TABULATED FUNCTION (TABFTN):       i       X(i)          Y(i)'

        DO n = 1,NUMTF
          WRITE (IO_UNIT%LSDO,'(/1X,A,I8,3X,25(''=''))')                       &
     &      'TABULATED FUNCTION (Ftn ID):',                                    &
     &      n
!!
          WRITE (IO_UNIT%LSDO,'((1X,A,I8,2(1PE14.4)))')                        &
     &     (                                                                   &
     &     'TABULATED FUNCTION (TABFTN):',                                     &
     &     i,                                                                  &
     &     TABULATED_FUNCTION(n)%X(i),                                         &
     &     TABULATED_FUNCTION(n)%Y(i),                                         &
     &     i = 1,TABULATED_FUNCTION(n)%Number_of_Pairs                         &
     &     )
        ENDDO

      ENDIF
!!
!! Nonreflecting boundary conditions.
!!
      IF (NUMNR .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NONREFLECTING_BC (NRBC):   BC ID, Set ID,  N-one,  N-two'

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8))')                                    &
     &    (                                                                    &
     &    'NONREFLECTING_BC (NRBC):',                                          &
     &    NONREFLECTING_BC(n)%NRID,     & ! Nonreflecting boundary condition ID  &
     &    NONREFLECTING_BC(n)%SetID,    & ! Segment set ID (-n/n=-SegID/SetID)   &
     &    NONREFLECTING_BC(n)%INRbgn,   & ! Starting location in NRBC_DATA       &
     &    NONREFLECTING_BC(n)%INRend,   & ! Ending location in NRBC_DATA         &
     &    n = 1,NUMNR                                                          &
     &    )

      ENDIF
!!
!! Nonreflecting boundary condition data.
!!
      IF (NUMND .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NRBC_DATA (...): N-1/N-2, Seg ID, HPT ID,HPT Typ'

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8))')                                    &
     &    (                                                                    &
     &    'NRBC_DATA (...):',                                                  &
     &    n,                                                                   &
     &    NRBC_DATA(n)%SegID,   & ! Boundary segment ID                        &
     &    NRBC_DATA(n)%Mel,     & ! Solid element w/w segment is associated    &
     &    NRBC_DATA(n)%Type,    & ! Element type, 0/1/2=hexa/penta/tetra
     &    n = 1,NUMND                                                          &
     &    )

      ENDIF
!!
!! Sliding interfaces.
!!
      IF (NUMSI .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SLIDING_INTERFACE (SLIDE):'//                                       &
     &    '   SI ID,Typ1, Set ID,Typ2, Set ID,Type,Sym,==Coef Fric=,'//        &
     &    '==T Begin===,===T End====,===Factor===,==Capture===,===Bo'//        &
     &    'rder==='

        WRITE (IO_UNIT%LSDO,'((1X,A,3(I8,I5),I4,6(1PE13.4)))')                 &
     &    (                                                                    &
     &    'SLIDING_INTERFACE (SLIDE):',                                        &
     &    SLIDING_INTERFACE(n)%SIID,     & ! Sliding interface ID                &
     &    SLIDING_INTERFACE(n)%Typ1,     & ! Type (0/1=segment-set/node-set)     &
     &    SLIDING_INTERFACE(n)%S1ID,     & ! Side 1 ID                           &
     &    SLIDING_INTERFACE(n)%Typ2,     & ! Type (0/1=segment-set/node-set)     &
     &    SLIDING_INTERFACE(n)%S2ID,     & ! Side 2 ID (0=single-surface interfac&
     &    SLIDING_INTERFACE(n)%Type,     & ! Type (n=0,1,2,3,...) [Type 0 only 7/&
     &    SLIDING_INTERFACE(n)%Isym,     & ! Flag, 0/1/2=symmetric/master-slave/v&
     &    SLIDING_INTERFACE(n)%CoF,      & ! Coefficient of friction             &
     &    SLIDING_INTERFACE(n)%Begin,    & ! Begin calculations                  &
     &    SLIDING_INTERFACE(n)%End,      & ! End calculations                    &
     &    SLIDING_INTERFACE(n)%Factor,   & ! Impact force relaxation factor      &
     &    SLIDING_INTERFACE(n)%Capture,  & ! Capture distance for penetration tes&
     &    SLIDING_INTERFACE(n)%Border,   & ! El exterior border for penetration te
     &    n = 1,NUMSI                                                          &
     &    )

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SLIDING_INTERFACE (SLIDE):'//                                       &
     &    '   SI ID,Typ1, Set ID,Typ2, Set ID,Type,Sym,ISN1bgn,ISN1e'//        &
     &    'nd,ISN2bgn,ISN2end,ICE1bgn,ICE1end,ICE2bgn,ICE2end'

        WRITE (IO_UNIT%LSDO,'((1X,A,3(I8,I5),I4,8(I8)))')                      &
     &    (                                                                    &
     &    'SLIDING_INTERFACE (SLIDE):',                                        &
     &    SLIDING_INTERFACE(n)%SIID,     & ! Sliding interface ID                &
     &    SLIDING_INTERFACE(n)%Typ1,     & ! Type (0/1=segment-set/node-set)     &
     &    SLIDING_INTERFACE(n)%S1ID,     & ! Side 1 ID                           &
     &    SLIDING_INTERFACE(n)%Typ2,     & ! Type (0/1=segment-set/node-set)     &
     &    SLIDING_INTERFACE(n)%S2ID,     & ! Side 2 ID (0=single-surface interfac&
     &    SLIDING_INTERFACE(n)%Type,     & ! Type (n=0,1,2,3,...) [Type 0 only 7/&
     &    SLIDING_INTERFACE(n)%Isym,     & ! Flag, 0/1/2=symmetric/master-slave/vi
     &    SLIDING_INTERFACE(n)%ISN1bgn,  & ! Starting location in SLIDING_NODE, s&
     &    SLIDING_INTERFACE(n)%ISN1end,  & ! Ending   location in SLIDING_NODE, s&
     &    SLIDING_INTERFACE(n)%ISN2bgn,  & ! Starting location in SLIDING_NODE, s&
     &    SLIDING_INTERFACE(n)%ISN2end,  & ! Ending   location in SLIDING_NODE, s&
     &    SLIDING_INTERFACE(n)%ICE1bgn,  & ! Starting location in CONTACT_SURFACE&
     &    SLIDING_INTERFACE(n)%ICE1end,  & ! Ending   location in CONTACT_SURFACE&
     &    SLIDING_INTERFACE(n)%ICE2bgn,  & ! Starting location in CONTACT_SURFACE&
     &    SLIDING_INTERFACE(n)%ICE2end,  & ! Ending   location in CONTACT_SURFACE&
     &    n = 1,NUMSI                                                          &
     &    )

      ENDIF
!!
!! Nodal point sets.
!!
      IF (NUMNS .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NODE_SET (NPSET):  Set ID,  N-one,  N-two, Flag'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,3X,A))')                               &
     &    (                                                                    &
     &    'NODE_SET (NPSET):',                                                 &
     &    n,                                                                   &
     &    NODE_SET(n)%Istart,   & ! Starting location                          &
     &    NODE_SET(n)%Iend,     & ! Ending location                            &
     &    NODE_SET(n)%Flag,     & ! ALL = All nodes, elements or segments.     &
     &    n = 1,NUMNS                                                          &
     &    )

      ENDIF
!!
!! Element sets.
!!
      IF (NUMES .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'ELEMENT_SET (ELSET):  Set ID,  N-one,  N-two, Flag'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,3X,A))')                               &
     &    (                                                                    &
     &    'ELEMENT_SET (ELSET):',                                              &
     &    n,                                                                   &
     &    ELEMENT_SET(n)%Istart,& ! Starting location                          &
     &    ELEMENT_SET(n)%Iend,  & ! Ending location                            &
     &    ELEMENT_SET(n)%Flag,  & ! ALL = All nodes, elements or segments.     &
     &    n = 1,NUMES                                                          &
     &    )

      ENDIF
!!
!! Segment sets.
!!
      IF (NUMSS .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SEGMENT_SET (SEGSET):  Set ID,  N-one,  N-two, Flag'

        WRITE (IO_UNIT%LSDO,'((1X,A,3I8,3X,A))')                               &
     &    (                                                                    &
     &    'SEGMENT_SET (SEGSET):',                                             &
     &    n,                                                                   &
     &    SEGMENT_SET(n)%Istart,& ! Starting location                          &
     &    SEGMENT_SET(n)%Iend,  & ! Ending location                            &
     &    SEGMENT_SET(n)%Flag,  & ! ALL = All nodes, elements or segments.     &
     &    n = 1,NUMSS                                                          &
     &    )

      ENDIF
!!
!! Concatenated nodal point sets.
!!
      IF (NUMNE .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A,80(''=''))')                               &
     &    'NNPSETS (...):'

        WRITE (IO_UNIT%LSDO,'((1X,A,10I8))')                                   &
     &    (                                                                    &
     &    'NNPSETS (...):',                                                    &
     &    (NNPSETS(i), i = n,MIN(NUMNE,n+9)),                                  &
     &    n = 1,NUMNE,10                                                       &
     &    )

      ENDIF
!!
!! Concatenated element sets.
!!
      IF (NUMEE .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A,80(''=''))')                               &
     &    'NELSETS (...):'

        WRITE (IO_UNIT%LSDO,'((1X,A,10I8))')                                   &
     &    (                                                                    &
     &    'NELSETS (...):',                                                    &
     &    (NELSETS(i), i = n,MIN(NUMEE,n+9)),                                  &
     &    n = 1,NUMEE,10                                                       &
     &    )

      ENDIF
!!
!! Concatenated segment sets.
!!
      IF (NUMSE .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A,80(''=''))')                               &
     &    'NSGSETS (...):'

        WRITE (IO_UNIT%LSDO,'((1X,A,10I8))')                                   &
     &    (                                                                    &
     &    'NSGSETS (...):',                                                    &
     &    (NSGSETS(i), i = n,MIN(NUMSE,n+9)),                                  &
     &    n = 1,NUMSE,10                                                       &
     &    )

      ENDIF
!!
!! Sliding node arrays.
!!
      IF (NUMSN .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SLIDING_NODE:     Nsn,    Mce'

        WRITE (IO_UNIT%LSDO,'((1X,A,2I8))')                                    &
     &    (                                                                    &
     &    'SLIDING_NODE:',                                                     &
     &    SLIDING_NODE(n)%Nsn,                                                 &
     &    SLIDING_NODE(n)%Mce,                                                 &
     &    n = 1,NUMSN                                                          &
     &    )

      ENDIF
!!
!! Sliding interface contact surface array.
!!
      IF (NUMCE .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'CONTACT_SURFACE: Surf ID,============IX(1:4)============,'//        &
     &    '============NX(1:4)============'

        WRITE (IO_UNIT%LSDO,'((1X,A,9I8))')                                    &
     &    (                                                                    &
     &    'CONTACT_SURFACE:',                                                  &
     &    n,                                                                   &
     &    CONTACT_SURFACE(n)%IX(1),                                            &
     &    CONTACT_SURFACE(n)%IX(2),                                            &
     &    CONTACT_SURFACE(n)%IX(3),                                            &
     &    CONTACT_SURFACE(n)%IX(4),                                            &
     &    CONTACT_SURFACE(n)%NX(1),                                            &
     &    CONTACT_SURFACE(n)%NX(2),                                            &
     &    CONTACT_SURFACE(n)%NX(3),                                            &
     &    CONTACT_SURFACE(n)%NX(4),                                            &
     &    n = 1,NUMCE                                                          &
     &    )

      ENDIF
!!
!! Sliding interface contact node array.
!!
      IF (NUMSN .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'CONTACT_NODE:   NP ID'

        WRITE (IO_UNIT%LSDO,'((1X,A,I8))')                                     &
     &    (                                                                    &
     &    'CONTACT_NODE:',                                                     &
     &    NODE(CONTACT_NODE(n)%NPID)%ID,                                       &
     &    n = 1,NUMCN                                                          &
     &    )

      ENDIF
!!
!! Segments.
!!
      IF (NUMSG .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SEGMENT (SEGMENT):  Seg ID,Part ID,=============IX(1:4)=='//        &
     &    '========='

        WRITE (IO_UNIT%LSDO,'((1X,A,6I8))')                                    &
     &    (                                                                    &
     &    'SEGMENT (SEGMENT):',                                                &
     &    n,                                                                   &
     &    SEGMENT(n)%PAR%ParID,  & ! Defining ID  (Reset: points to internal ID) &
     &    SEGMENT(n)%PAR%IX(1),  & ! NP Indicies  (Reset: internal node nums)  &
     &    SEGMENT(n)%PAR%IX(2),  & ! NP Indicies  (Reset: internal node nums)  &
     &    SEGMENT(n)%PAR%IX(3),  & ! NP Indicies  (Reset: internal node nums)  &
     &    SEGMENT(n)%PAR%IX(4),  & ! NP Indicies  (Reset: internal node nums)  &
     &    n = 1,NUMSG                                                          &
     &    )

      ENDIF
!!
!! 8-Node hexahedral solid elements.
!!
      IF (NUMHX .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'HEXAH (HXEL):   EleID,  ParID,  MatID,  SecID,==========='//        &
     &    '=================IX(1:8)============================='

        WRITE (IO_UNIT%LSDO,'((1X,A,12I8))')                                   &
     &    (                                                                    &
     &    'HEXAH (HXEL):',                                                     &
     &    n,                                                                   &
     &    HEXAH(n)%PAR%ParID,                                                  &
     &    HEXAH(n)%PAR%MatID,                                                  &
     &    HEXAH(n)%PAR%LupID,                                                  &
     &    HEXAH(n)%PAR%IX(1),                                                  &
     &    HEXAH(n)%PAR%IX(2),                                                  &
     &    HEXAH(n)%PAR%IX(3),                                                  &
     &    HEXAH(n)%PAR%IX(4),                                                  &
     &    HEXAH(n)%PAR%IX(5),                                                  &
     &    HEXAH(n)%PAR%IX(6),                                                  &
     &    HEXAH(n)%PAR%IX(7),                                                  &
     &    HEXAH(n)%PAR%IX(8),                                                  &
     &    n = 1,NUMHX                                                          &
     &    )

      ENDIF
!!
!! 6-Node pentahedral solid wedge elements.
!!
      IF (NUMPX .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PENTA (PXEL):   EleID,  ParID,  MatID,  SecID,==========='//        &
     &    '=========IX(1:6)====================='

        WRITE (IO_UNIT%LSDO,'((1X,A,10I8))')                                   &
     &    (                                                                    &
     &    'PENTA (PXEL):',                                                     &
     &    n,                                                                   &
     &    PENTA(n)%PAR%ParID,                                                  &
     &    PENTA(n)%PAR%MatID,                                                  &
     &    PENTA(n)%PAR%LupID,                                                  &
     &    PENTA(n)%PAR%IX(1),                                                  &
     &    PENTA(n)%PAR%IX(2),                                                  &
     &    PENTA(n)%PAR%IX(3),                                                  &
     &    PENTA(n)%PAR%IX(4),                                                  &
     &    PENTA(n)%PAR%IX(5),                                                  &
     &    PENTA(n)%PAR%IX(6),                                                  &
     &    n = 1,NUMPX                                                          &
     &    )

      ENDIF
!!
!! 4-Node tetrahedral solid elements.
!!
      IF (NUMTX .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'TETRA (TXEL):   EleID,  ParID,  MatID,  SecID,==========='//        &
     &    '=IX(1:4)============'

        WRITE (IO_UNIT%LSDO,'((1X,A,8I8))')                                    &
     &    (                                                                    &
     &    'TETRA (TXEL):',                                                     &
     &    n,                                                                   &
     &    TETRA(n)%PAR%ParID,                                                  &
     &    TETRA(n)%PAR%MatID,                                                  &
     &    TETRA(n)%PAR%LupID,                                                  &
     &    TETRA(n)%PAR%IX(1),                                                  &
     &    TETRA(n)%PAR%IX(2),                                                  &
     &    TETRA(n)%PAR%IX(3),                                                  &
     &    TETRA(n)%PAR%IX(4),                                                  &
     &    n = 1,NUMTX                                                          &
     &    )

      ENDIF
!!
!! 3-Node triangular membrane elements.
!!
      IF (NUMM3 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'MEMBT (M3EL):   EleID,  ParID,  MatID,  SecID,========IX('//        &
     &    '1:3)========'

        WRITE (IO_UNIT%LSDO,'((1X,A,7I8))')                                    &
     &    (                                                                    &
     &    'MEMBT (M3EL):',                                                     &
     &    n,                                                                   &
     &    MEMBT(n)%PAR%ParID,                                                  &
     &    MEMBT(n)%PAR%MatID,                                                  &
     &    MEMBT(n)%PAR%SecID,                                                  &
     &    MEMBT(n)%PAR%IX(1),                                                  &
     &    MEMBT(n)%PAR%IX(2),                                                  &
     &    MEMBT(n)%PAR%IX(3),                                                  &
     &    n = 1,NUMM3                                                          &
     &    )

      ENDIF
!!
!! 4-Node quadrilateral membrane elements.
!!
      IF (NUMM4 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'MEMBQ (M4EL):   EleID,  ParID,  MatID,  SecID,==========='//        &
     &    '=IX(1:4)============'

        WRITE (IO_UNIT%LSDO,'((1X,A,8I8))')                                    &
     &    (                                                                    &
     &    'MEMBQ (M4EL):',                                                     &
     &    n,                                                                   &
     &    MEMBQ(n)%PAR%ParID,                                                  &
     &    MEMBQ(n)%PAR%MatID,                                                  &
     &    MEMBQ(n)%PAR%SecID,                                                  &
     &    MEMBQ(n)%PAR%IX(1),                                                  &
     &    MEMBQ(n)%PAR%IX(2),                                                  &
     &    MEMBQ(n)%PAR%IX(3),                                                  &
     &    MEMBQ(n)%PAR%IX(4),                                                  &
     &    n = 1,NUMM4                                                          &
     &    )

      ENDIF
!!
!! 2-Node axial force truss elements.
!!
      IF (NUMTR .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'TRUSS (TRUSS):   EleID,  ParID,  MatID,  SecID,====IX(1:2'//        &
     &    ')===='

        WRITE (IO_UNIT%LSDO,'((1X,A,6I8))')                                    &
     &    (                                                                    &
     &    'TRUSS (TRUSS):',                                                    &
     &    n,                                                                   &
     &    TRUSS(n)%PAR%ParID,                                                  &
     &    TRUSS(n)%PAR%MatID,                                                  &
     &    TRUSS(n)%PAR%SecID,                                                  &
     &    TRUSS(n)%PAR%IX(1),                                                  &
     &    TRUSS(n)%PAR%IX(2),                                                  &
     &    n = 1,NUMTR                                                          &
     &    )

      ENDIF
!!
!! 3-Node triangular plate bending elements,
!!
      IF (NUMP3 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PLATT (P3EL):   EleID,  ParID,  MatID,  SecID,========IX('//        &
     &    '1:3)========='

        WRITE (IO_UNIT%LSDO,'((1X,A,7I8))')                                    &
     &    (                                                                    &
     &    'PLATT (P3EL):',                                                     &
     &    n,                                                                   &
     &    PLATT(n)%PAR%ParID,                                                  &
     &    PLATT(n)%PAR%MatID,                                                  &
     &    PLATT(n)%PAR%SecID,                                                  &
     &    PLATT(n)%PAR%IX(1),                                                  &
     &    PLATT(n)%PAR%IX(2),                                                  &
     &    PLATT(n)%PAR%IX(3),                                                  &
     &    n = 1,NUMP3                                                          &
     &    )

      ENDIF
!!
!! 4-Node quadrilateral plate bending elements.
!!
      IF (NUMP4 .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PLATQ (P4EL):   EleID,  ParID,  MatID,  SecID,==========='//        &
     &    '=IX(1:4)============'

        WRITE (IO_UNIT%LSDO,'((1X,A,8I8))')                                    &
     &    (                                                                    &
     &    'PLATQ (P4EL):',                                                     &
     &    n,                                                                   &
     &    PLATQ(n)%PAR%ParID,                                                  &
     &    PLATQ(n)%PAR%MatID,                                                  &
     &    PLATQ(n)%PAR%SecID,                                                  &
     &    PLATQ(n)%PAR%IX(1),                                                  &
     &    PLATQ(n)%PAR%IX(2),                                                  &
     &    PLATQ(n)%PAR%IX(3),                                                  &
     &    PLATQ(n)%PAR%IX(4),                                                  &
     &    n = 1,NUMP4                                                          &
     &    )

      ENDIF
!!
!! 2-Node axial force, torsional and bending beam.
!!
      IF (NUMBM .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'BEAMS (BEAM):   EleID,  ParID,  MatID,  SecID,====IX(1:2)'//        &
     &    '====,=====Ax======,=====Ay======,=====Az======'

        WRITE (IO_UNIT%LSDO,'((1X,A,6I8,3(1PE14.4)))')                         &
     &    (                                                                    &
     &    'BEAM  (BEAM):',                                                     &
     &    n,                                                                   &
     &    BEAM(n)%PAR%ParID,                                                   &
     &    BEAM(n)%PAR%MatID,                                                   &
     &    BEAM(n)%PAR%SecID,                                                   &
     &    BEAM(n)%PAR%IX(1),                                                   &
     &    BEAM(n)%PAR%IX(2),                                                   &
     &    BEAM(N)%RES%Zaxis(1),                                                &
     &    BEAM(N)%RES%Zaxis(2),                                                &
     &    BEAM(N)%RES%Zaxis(3),                                                &
     &    n = 1,NUMBM                                                          &
     &    )

      ENDIF
!!
!! 2-Node axial force spring element (massless).
!!
      IF (NUMSP .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'SPRING(SPRING):   EleID,  ParID,  MatID,   Type,====IX(1:'//        &
     &    '2)===='

        WRITE (IO_UNIT%LSDO,'((1X,A,6I8))')                                    &
     &    (                                                                    &
     &    'SPRING(SPRING):',                                                   &
     &    n,                                                                   &
     &    SPRING(n)%PAR%ParID,                                                 &
     &    SPRING(n)%PAR%MatID,                                                 &
     &    SPRING(n)%PAR%Type,                                                  &
     &    SPRING(n)%PAR%IX(1),                                                 &
     &    SPRING(n)%PAR%IX(2),                                                 &
     &    n = 1,NUMSP                                                          &
     &    )

      ENDIF
!!
!! 2-Node axial force damper element (massless).
!!
      IF (NUMDM .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'DAMPER(DAMPER):   EleID,  ParID,  MatID,   Type,====IX(1:'//        &
     &    '2)===='

        WRITE (IO_UNIT%LSDO,'((1X,A,6I8))')                                    &
     &    (                                                                    &
     &    'DAMPER(DAMPER):',                                                   &
     &    n,                                                                   &
     &    DAMPER(n)%PAR%ParID,                                                 &
     &    DAMPER(n)%PAR%MatID,                                                 &
     &    DAMPER(n)%PAR%Type,                                                  &
     &    DAMPER(n)%PAR%IX(1),                                                 &
     &    DAMPER(n)%PAR%IX(2),                                                 &
     &    n = 1,NUMDM                                                          &
     &    )

      ENDIF
!!
!! Paired plate data: plate ID's, and angles.
!!
      IF (NUMPP .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'PLATE_PAIR (...):'//                                                &
     &    ' Pair ID, 1st El,El Side, 2nd El,El Side,=Init CosPhi=,=C'//        &
     &    'urr CosPhi=,==CosPhi Dot='

        WRITE (IO_UNIT%LSDO,'((1X,A,5I8,3(1PE14.4)))')                         &
     &    (                                                                    &
     &    'PLATE_PAIR (...):',                                                 &
     &    PLATE_PAIR(n)%ID,             & ! Internal ID/pointer to next location &
     &    PLATE_PAIR(n)%IDS1%ID,        & ! First element ID                     &
     &    PLATE_PAIR(n)%IDS1%IS,        & ! Element side                         &
     &    PLATE_PAIR(n)%IDS2%ID,        & ! Second element ID                    &
     &    PLATE_PAIR(n)%IDS2%IS,        & ! Elements side                        &
     &    PLATE_PAIR(n)%CSPINI,         & ! Cosine of initial angle between plate&
     &    PLATE_PAIR(n)%CosPhi,         & ! Cosine of angle between plates       &
     &    PLATE_PAIR(n)%CosDot,         & ! Rate at which Cos(Phi) is changing   &
     &    n = 1,NUMPP                                                          &
     &    )

      ENDIF
!!
!! Nodal point data: ID's, rigid body linked lists, motion, force et cetera.
!!
      IF (NUMNP .GT. 0) THEN

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NODE (NPT):   NP ID,SubCycI,RB Link,Rot Ptr'

        WRITE (IO_UNIT%LSDO,'((1X,A,4I8))')                                    &
     &    (                                                                    &
     &    'NODE (NPT):',                                                       &
     &    n,                                                                   &
     &    NODE(n)%ISI,                                                         &
     &    NODE(n)%IRB,                                                         &
     &    NODE(n)%IRT,                                                         &
     &    n = 1,NUMNP                                                          &
     &    )

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'NODE (NPT):   NP ID,=====Time====,====DTlast===,====DTnex'//        &
     &    't===,=====Mass====,=====Minv===='

        WRITE (IO_UNIT%LSDO,'((1X,A,I8,5(1PE14.4)))')                          &
     &    (                                                                    &
     &    'NODE (NPT):',                                                       &
     &    n,                                                                   &
     &    NODE(n)%Time,                                                        &
     &    NODE(n)%DTlast,                                                      &
     &    NODE(n)%DTnext,                                                      &
     &    NODE(n)%Mass,                                                        &
     &    NODE(n)%Minv,                                                        &
     &    n = 1,NUMNP                                                          &
     &    )

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'MOTION (NPT):===ID===,=================Px,Py,Pz=========='//        &
     &    '======,=================Ux,Uy,Uz================'

        WRITE (IO_UNIT%LSDO,'((1X,A,I8,6(1PE14.4)))')                          &
     &    (                                                                    &
     &    'MOTION (NPT):',                                                     &
     &    n,                                                                   &
     &    MOTION(n)%Px,         & ! Initial x-position                         &
     &    MOTION(n)%Py,         & ! Initial y-position                         &
     &    MOTION(n)%Pz,         & ! Initial z-position                         &
     &    MOTION(n)%Ux,         & ! X displacement                             &
     &    MOTION(n)%Uy,         & ! Y displacement                             &
     &    MOTION(n)%Uz,         & ! Z displacement                             &
     &    n = 1,NUMNP                                                          &
     &    )

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'MOTION (NPT):===ID===,=================Vx,Vy,Vz=========='//        &
     &    '======,=================Ax,Ay,Az================'

        WRITE (IO_UNIT%LSDO,'((1X,A,I8,6(1PE14.4)))')                          &
     &    (                                                                    &
     &    'MOTION (NPT):',                                                     &
     &    n,                                                                   &
     &    MOTION(n)%Vx,         & ! X velocity                                 &
     &    MOTION(n)%Vy,         & ! Y velocity                                 &
     &    MOTION(n)%Vz,         & ! Z velocity                                 &
     &    MOTION(n)%Ax,         & ! X acceleration                             &
     &    MOTION(n)%Ay,         & ! Y acceleration                             &
     &    MOTION(n)%Az,         & ! Z acceleration                             &
     &    n = 1,NUMNP                                                          &
     &    )

        CALL NEW_PAGE (IO_UNIT%LSDO, PAGE_NUMBER)

        WRITE (IO_UNIT%LSDO,'(/1X,A)')                                         &
     &    'FORCE (...):===============Xext,Yext,Zext=============,=='//        &
     &    '============Xint,Yint,Zint============='

        WRITE (IO_UNIT%LSDO,'((1X,A,6(1PE14.4)))')                             &
     &    (                                                                    &
     &    'FORCE (...):',                                                      &
     &    FORCE(n)%Xext,                & ! X direction external force         &
     &    FORCE(n)%Yext,                & ! Y direction external force         &
     &    FORCE(n)%Zext,                & ! Z direction external force         &
     &    FORCE(n)%Xint,                & ! X direction internal force         &
     &    FORCE(n)%Yint,                & ! Y direction internal force         &
     &    FORCE(n)%Zint,                & ! Z direction internal force         &
     &    n = 1,NUMNP                                                          &
     &    )

      ENDIF
!!
      CLOSE (UNIT=IO_UNIT%LSDO, STATUS='KEEP')
!!
      RETURN
      END
!!_
      SUBROUTINE NEW_PAGE (IOUNIT,PAGE_NUMBER)
!!
!! Copyright (c) by KEY Associates; 18-APR-1992 11:41:00
!!
!! Purpose: Generate new listing page with heading. This module is used
!! with each of the three ASCII output files:
!!
!!      (1) fmaelo, Execution_Log_Output,
!!      (2) fmasdo, Simulation_Data_Output,
!!      (3) fmacro, Computed_Results_Output.
!!
!! Each of these files is required to keep its own page number register.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER, SAVE::                                                        &
     &           HEADLINE*132
      INTEGER                                                                  &
     &          IOUNIT,                                                        &
     &          PAGE_NUMBER
      LOGICAL, SAVE:: FIRST = .TRUE.
!!
      IF (FIRST) THEN
        LDA = LEN (JOB_ID_RECORD%CURRENT%DATEE)
        LTI = LEN (JOB_ID_RECORD%CURRENT%TIMEE)
        LNB = LAST_NONBLANK (JOB_ID_RECORD%CURRENT%TITLE)
        LST = 58-LNB/2
        LPR = LEN (JOB_ID_RECORD%PROGRAM)
        LVE = LEN (JOB_ID_RECORD%VERSION)
        NXT = 1
        HEADLINE(NXT:) = '1'//JOB_ID_RECORD%PROGRAM
        NXT = NXT + LPR + 2
        HEADLINE(NXT:) = JOB_ID_RECORD%VERSION
        NXT = LST
        HEADLINE(NXT:) = JOB_ID_RECORD%CURRENT%TITLE(1:LNB)
        NXT = 119 - LDA - LTI
!SPEC_CPU2000        HEADLINE(NXT:) = JOB_ID_RECORD%CURRENT%DATEE
        NXT = NXT + LDA + 1
!SPEC_CPU2000        HEADLINE(NXT:) = JOB_ID_RECORD%CURRENT%TIMEE
        NXT = NXT + LTI + 1
        HEADLINE(NXT:) = 'Page No: '
        FIRST = .FALSE.
      ENDIF
!!
      PAGE_NUMBER = PAGE_NUMBER + 1
      WRITE(HEADLINE(130:132),'(I3)') PAGE_NUMBER
!!
      WRITE (IOUNIT,'(A/)') HEADLINE
!!
      RETURN
      END
!!_
      SUBROUTINE PRINT_COUNTERS (IOUNIT,INTRODUCTION,ICOUNT,PAGE_NUMBER)
!!
!! Copyright (c) by KEY Associates;  2-APR-1993 19:35:28.50
!!
!! Purpose: Print table of counter values using "introducer" to form label.
!! Note: This module expects the counter array ICOUNT to be in a predifined
!! order such that the following "hard coded" counter descriptions make sense.
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          DARG*32,                                                       &
     &          LABEL*43,                                                      &
     &          INTRODUCER*10,                                                 &
     &          INTRODUCTION*(*)
      INTEGER                                                                  &
     &          IOUNIT,                                                        &
     &          ICOUNT(*),                                                     &
     &          PAGE_NUMBER
      CHARACTER*33 :: DOTS = ' . . . . . . . . . . . . . . . . '
!!
!! Function statements for building labels.
!!
!!BODO      LABEL(DARG) = INTRODUCER//DARG(:LEN(DARG))//DOTS(LEN(DARG)+1:)
!!
      CALL NEW_PAGE (IOUNIT, PAGE_NUMBER)
!!
!! Form introducer from introduction string.
!!
      INTRODUCER(1:) = INTRODUCTION(1:)
!!
!! List values for all counters.
!!
      WRITE (IOUNIT,'(4X,A//(2(8X,A,I8)))')                                    &
     &  INTRODUCER//' DATA ITEMS:',                                            &
     &  'Input Files                     ', ICOUNT(1),                         &
     &  'QA Records                      ', ICOUNT(2),                         &
     &  'Nodal Points                    ', ICOUNT(3),                         &
     &  'Nodal Points W/ Rotations       ', ICOUNT(61)-ICOUNT(3),              &
     &  'Elements, Total                 ', ICOUNT(4),                         &
     &  '8-Node Hexahedrons              ', ICOUNT(5),                         &
     &  '6-Node Pentahedrons             ', ICOUNT(6),                         &
     &  '4-Node Tetrahedrons             ', ICOUNT(7),                         &
     &  '8-Node Layered Solids           ', ICOUNT(8),                         &
     &  'Layered-Solid Hexahedrons       ', ICOUNT(9),                         &
     &  'Layered-Solid Membranes         ', ICOUNT(10),                        &
     &  '4-Node Membranes                ', ICOUNT(11),                        &
     &  '3-Node Membranes                ', ICOUNT(12),                        &
     &  '2-Node Trusses                  ', ICOUNT(13),                        &
     &  '4-Node Plates                   ', ICOUNT(14),                        &
     &  '3-Node Plates                   ', ICOUNT(15),                        &
     &  '2-Node Beams                    ', ICOUNT(16),                        &
     &  '2-Node Springs                  ', ICOUNT(17),                        &
     &  '2-Node Dampers                  ', ICOUNT(18),                        &
     &  '4-Node Segments                 ', ICOUNT(19),                        &
     &  'Displacement BC''s              ', ICOUNT(20),                        &
     &  'Tied BC''s                      ', ICOUNT(21),                        &
     &  'Spot Welds                      ', ICOUNT(22),                        &
     &  'Wall BC''s                      ', ICOUNT(23),                        &
     &  'Body Force Fields               ', ICOUNT(24),                        &
     &  'Pressure BC''s                  ', ICOUNT(25),                        &
     &  'Concentrated Force BC''s        ', ICOUNT(26),                        &
     &  'Spring BC''s                    ', ICOUNT(27),                        &
     &  'Damper BC''s                    ', ICOUNT(28),                        &
     &  'Periodic BC''s                  ', ICOUNT(29),                        &
     &  'Nonreflecting BC''s             ', ICOUNT(30),                        &
     &  'NRBC Data Entries               ', ICOUNT(31),                        &
     &  'Fastners                        ', ICOUNT(32),                        &
     &  'Interface Times                 ', ICOUNT(33),                        &
     &  'Sliding Interfaces              ', ICOUNT(34),                        &
     &  'Sliding Nodes                   ', ICOUNT(35),                        &
     &  'Contact Elements                ', ICOUNT(36),                        &
     &  'Contact Nodes                   ', ICOUNT(37),                        &
     &  'Max Contact Element Set         ', ICOUNT(38),                        &
     &  'Node Sets                       ', ICOUNT(39),                        &
     &  'Node Entries                    ', ICOUNT(40),                        &
     &  'Element Sets                    ', ICOUNT(41),                        &
     &  'Element Entries                 ', ICOUNT(42),                        &
     &  'Segment Sets                    ', ICOUNT(43),                        &
     &  'Segment Entries                 ', ICOUNT(44),                        &
     &  'Materials                       ', ICOUNT(45),                        &
     &  'Lay-Ups                         ', ICOUNT(46),                        &
     &  'Tabulated Functions             ', ICOUNT(47),                        &
     &  'Function Pairs                  ', ICOUNT(48),                        &
     &  'Results Files                   ', ICOUNT(49),                        &
     &  'Parameter Values                ', ICOUNT(50),                        &
     &  'Beam Sections                   ', ICOUNT(51),                        &
     &  'Plate Sections                  ', ICOUNT(52),                        &
     &  '1-D Strain Gauges               ', ICOUNT(53),                        &
     &  '2-D Strain Gauges               ', ICOUNT(54),                        &
     &  '3-D Strain Gauges               ', ICOUNT(55),                        &
     &  'Mass Properties                 ', ICOUNT(56),                        &
     &  'Substitute/Added RB Masses      ', ICOUNT(57),                        &
     &  'Concentrated Masses             ', ICOUNT(58),                        &
     &  'Initial Conditions              ', ICOUNT(59),                        &
     &  'Rigid Bodies                    ', ICOUNT(60),                        &
     &  'Shell Stresses                  ', ICOUNT(62),                        &
     &  'State Variables                 ', ICOUNT(63),                        &
     &  'Plate Pairs                     ', ICOUNT(64),                        &
     &  'Constrained Nodes               ', ICOUNT(65)
!!
      RETURN
      END
!!_
      SUBROUTINE TOTAL_MASS_REPORT (IOUNIT)
!!
!! Copyright (c) by KEY Associates;  4-FEB-1996 15:59:11.00
!!
!! Purpose: Report mass properties of each material wrt the center
!! of mass for each material. During element initialization, the
!! mass properties were "accumulated" by material. (For bodies based
!! on a small number of elements, these results will not match the
!! properties generated by continuum formulae.)
!!
      USE shared_common_data
      USE material_
      USE node_
      USE motion_
      USE massprop_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: IOUNIT  ! Logical file number for output
!!
      LOGICAL, SAVE:: FIRST = .TRUE.
!!
!! Finish calculations for each material domain.
!!
      IF (FIRST) THEN
        Total_Mass = 0.0
        DO N = 1,NUMMT
          IF (MATERIAL(N)%Mass .GT. 0.0) THEN
            QMass = MATERIAL(N)%Mass
            MATERIAL(N)%Xcm = MATERIAL(N)%Xcm/QMass
            MATERIAL(N)%Ycm = MATERIAL(N)%Ycm/QMass
            MATERIAL(N)%Zcm = MATERIAL(N)%Zcm/QMass
!!
!! Translate inertia tensor from the origin to the center of mass.
!!
            Xcm = MATERIAL(N)%Xcm
            Ycm = MATERIAL(N)%Ycm
            Zcm = MATERIAL(N)%Zcm
            RSQ = Xcm*Xcm + Ycm*Ycm + Zcm*Zcm
            MATERIAL(N)%Bxx = MATERIAL(N)%Bxx - (RSQ - Xcm*Xcm)*QMass
            MATERIAL(N)%Byy = MATERIAL(N)%Byy - (RSQ - Ycm*Ycm)*QMass
            MATERIAL(N)%Bzz = MATERIAL(N)%Bzz - (RSQ - Zcm*Zcm)*QMass
            MATERIAL(N)%Bxy = MATERIAL(N)%Bxy + (Xcm*Ycm)*QMass
            MATERIAL(N)%Bxz = MATERIAL(N)%Bxz + (Xcm*Zcm)*QMass
            MATERIAL(N)%Byz = MATERIAL(N)%Byz + (Ycm*Zcm)*QMass
            Total_Mass = Total_Mass + QMass
          ENDIF
        ENDDO
        Total_Mass_Xcm = 0.0
        Total_Mass_Ycm = 0.0
        Total_Mass_Zcm = 0.0
        IF (Total_Mass .GT. 0.0) THEN
          DO N = 1,NUMMT
            Total_Mass_Xcm = Total_Mass_Xcm +                                  &
     &        MATERIAL(N)%Mass * MATERIAL(N)%Xcm
            Total_Mass_Ycm = Total_Mass_Ycm +                                  &
     &        MATERIAL(N)%Mass * MATERIAL(N)%Ycm
            Total_Mass_Zcm = Total_Mass_Zcm +                                  &
     &        MATERIAL(N)%Mass * MATERIAL(N)%Zcm
          ENDDO
          Total_Mass_Xcm = Total_Mass_Xcm / Total_Mass
          Total_Mass_Ycm = Total_Mass_Ycm / Total_Mass
          Total_Mass_Zcm = Total_Mass_Zcm / Total_Mass
        ENDIF
        Bxx = 0.0
        Byy = 0.0
        Bzz = 0.0
        Bxy = 0.0
        Bxz = 0.0
        Byz = 0.0
        DO N = 1,NUMMT
          Xcm = MATERIAL(N)%Xcm
          Ycm = MATERIAL(N)%Ycm
          Zcm = MATERIAL(N)%Zcm
          QMass = MATERIAL(N)%Mass
          P2x = (Ycm - Total_Mass_Ycm)**2 + (Zcm - Total_Mass_Zcm)**2
          P2y = (Zcm - Total_Mass_Zcm)**2 + (Xcm - Total_Mass_Xcm)**2
          P2z = (Xcm - Total_Mass_Xcm)**2 + (Ycm - Total_Mass_Ycm)**2
          Bxx = Bxx + MATERIAL(N)%Bxx + P2x * QMass
          Byy = Byy + MATERIAL(N)%Byy + P2y * QMass
          Bzz = Bzz + MATERIAL(N)%Bzz + P2z * QMass
          Bxy = Bxy + MATERIAL(N)%Bxy -                                        &
     &      (Xcm-Total_Mass_Xcm)*(Ycm-Total_Mass_Ycm)*QMass
          Bxz = Bxz + MATERIAL(N)%Bxz -                                        &
     &      (Xcm-Total_Mass_Xcm)*(Zcm-Total_Mass_Zcm)*QMass
          Byz = Byz + MATERIAL(N)%Byz -                                        &
     &      (Ycm-Total_Mass_Ycm)*(Zcm-Total_Mass_Zcm)*QMass
        ENDDO
!!
!! Compute whole-body results based on nodal point masses. (We
!! do not expect NODE_SET and NNPSETS to be referenced in the
!! following call.)
!!
        MASSPROP(0)%MPID  = 1
        MASSPROP(0)%SetID = 0
        MASSPROP(0)%Irot  = 1
        NSAVE = NUMMP
        NUMMP = 0
        CALL MASS_PROPERTY_CALCULATION
        NUMMP = NSAVE
        FIRST = .FALSE.
      ENDIF
!!
!! Output material mass totals.
!!
      WRITE (IOUNIT,'(//1X,A)')                                                &
     &          '*** INFORMATION *** MASS TOTALS And INERTIA TENSORS'//        &
     &    ' By MATERIAL:'
!!
      WRITE (IOUNIT,'(/1X,A/(1X,I8,3X,4(1PE14.4)))')                           &
     &          'MATERIAL ID,==Total Mass=,=====Xcm=====,=====Ycm==='//        &
     &    '==,=====Zcm=====',                                                  &
     &          (                                                              &
     &          MATERIAL(n)%MatID,                                             &
     &          MATERIAL(n)%Mass,                                              &
     &          MATERIAL(n)%Xcm,                                               &
     &          MATERIAL(n)%Ycm,                                               &
     &          MATERIAL(n)%Zcm,                                               &
     &          n=1,NUMMT                                                      &
     &          )
!!
      WRITE (IOUNIT,'(/1X,A/(1X,I8,3X,6(1PE14.4)))')                           &
     &          'MATERIAL ID,=====Bxx=====,=====Byy=====,=====Bzz==='//        &
     &    '==,=====Bxy=====,=====Bxz=====,=====Byz=====',                      &
     &          (                                                              &
     &          MATERIAL(n)%MatID,                                             &
     &          MATERIAL(n)%Bxx,                                               &
     &          MATERIAL(n)%Byy,                                               &
     &          MATERIAL(n)%Bzz,                                               &
     &          MATERIAL(n)%Bxy,                                               &
     &          MATERIAL(n)%Bxz,                                               &
     &          MATERIAL(n)%Byz,                                               &
     &          n=1,NUMMT                                                      &
     &          )
      WRITE (IOUNIT,'(/(7X,A,1PE14.4))')                                       &
     &          'Total Mass For All Materials:',Total_Mass,                    &
     &          'Center Of Mass, X-Coordinate:',Total_Mass_Xcm,                &
     &          'Center Of Mass, Y-Coordinate:',Total_Mass_Ycm,                &
     &          'Center Of Mass, Z-Coordinate:',Total_Mass_Zcm,                &
     &          'Inertia (Center Of Mass) Bxx:',Bxx,                           &
     &          'Inertia (Center Of Mass) Byy:',Byy,                           &
     &          'Inertia (Center Of Mass) Bzz:',Bzz,                           &
     &          'Inertia (Center Of Mass) Bxy:',Bxy,                           &
     &          'Inertia (Center Of Mass) Bxz:',Bxz,                           &
     &          'Inertia (Center Of Mass) Byz:',Byz
      WRITE (IOUNIT,'(//A/A/(A))')                                             &
     &          'NOTES: (1) Inertia Tensor Components Are With Respe'//        &
     &    'ct To The Center Of Mass',                                          &
     &          '       (2) Does Not Include Concentrated Mass From '//        &
     &    'NPMASS Or RBMASS.',                                                 &
     &          '       (3) For Bodies Based On A Small Number Of El'//        &
     &    'ements, These Results',                                             &
     &          '           Will Not Match The Properties Generated '//        &
     &    'By Continuum Formulae.'

!!
!! Output results of mass property calculations.
!!
      WRITE (IOUNIT,'(//1X,A)')                                                &
     &          '*** INFORMATION *** MASS PROPERTIES DERIVED FROM AL'//        &
     &    'L NODAL POINTS:'
!!
      WRITE (IOUNIT,'(/(10X,A,1PE14.4))')                                      &
     &    'TOTAL MASS (ALL NODES):       ',MASSPROP(0)%Mass,                   &
     &    'Momentum, X-Direction:        ',MASSPROP(0)%Xmv,                    &
     &    'Momentum, Y-Direction:        ',MASSPROP(0)%Ymv,                    &
     &    'Momentum, Z-Direction:        ',MASSPROP(0)%Zmv,                    &
     &    'Center Of Mass, X-Coord:      ',MASSPROP(0)%Xcm,                    &
     &    'Center Of Mass, Y-Coord:      ',MASSPROP(0)%Ycm,                    &
     &    'Center Of Mass, Z-Coord:      ',MASSPROP(0)%Zcm,                    &
     &    'Center Of Mass, X-Velocity:   ',MASSPROP(0)%Vxcm,                   &
     &    'Center Of Mass, Y-Velocity:   ',MASSPROP(0)%Vycm,                   &
     &    'Center Of Mass, Z-Velocity:   ',MASSPROP(0)%Vzcm,                   &
     &    'Translational KE:             ',MASSPROP(0)%KE,                     &
     &    'Inertia (Center Of Mass) Bxx: ',MASSPROP(0)%B(1),                   &
     &    'Inertia (Center Of Mass) Byy: ',MASSPROP(0)%B(2),                   &
     &    'Inertia (Center Of Mass) Bzz: ',MASSPROP(0)%B(3),                   &
     &    'Inertia (Center Of Mass) Bxy: ',MASSPROP(0)%B(4),                   &
     &    'Inertia (Center Of Mass) Bxz: ',MASSPROP(0)%B(5),                   &
     &    'Inertia (Center Of Mass) Byz: ',MASSPROP(0)%B(6),                   &
     &    'Angular Momentum, X-Axis:     ',MASSPROP(0)%Oxmv,                   &
     &    'Angular Momentum, Y-Axis:     ',MASSPROP(0)%Oymv,                   &
     &    'Angular Momentum, Z-Axis:     ',MASSPROP(0)%Ozmv,                   &
     &    'Angular Velocity (rads/sec):  ',MASSPROP(0)%Omega,                  &
     &    'Axis Of Rotation, X-Component:',MASSPROP(0)%Ax,                     &
     &    'Axis Of Rotation, Y-Component:',MASSPROP(0)%Ay,                     &
     &    'Axis Of Rotation, Z-Component:',MASSPROP(0)%Az
!!
      RETURN
      END
!!_
      SUBROUTINE TIME_STEP_SUMMARY (IOUNIT)
!!
!! Copyright (c) by KEY Associates;  8-AUG-1993 15:43:49.13
!!
!! Purpose: Print summary of the smallest element critical time steps in each
!! element group.
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
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          IOUNIT
!!
      WRITE (IOUNIT,'(//1X,A)')                                                &
     &          '*** INFORMATION *** SUMMARY OF ELEMENT CRITICAL TIM'//        &
     &    'E STEPS:'
!!
!! Note: NPNDT is a PARAMETER in "common.inc"
!!
      IF (NUMHX .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'HXEL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'HXEL   ID:',                                                  &
     &          HEXAH(TIMSIM%Hexah(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTHex(i),                                               &
     &          i=1,MIN(NPNDT,NUMHX)                                           &
     &          )
      ENDIF
      IF (NUMPX .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'PXEL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'PXEL   ID:',                                                  &
     &          PENTA(TIMSIM%Penta(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTPen(i),                                               &
     &          i=1,MIN(NPNDT,NUMPX)                                           &
     &          )
      ENDIF
      IF (NUMTX .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'TXEL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'TXEL   ID:',                                                  &
     &          TETRA(TIMSIM%Tetra(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTTet(i),                                               &
     &          i=1,MIN(NPNDT,NUMTX)                                           &
     &          )
      ENDIF

      IF (NUMLS .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'LSEL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'LSEL   ID:',                                                  &
     &          LSOLD(TIMSIM%LSold(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTLYS(i),                                               &
     &          i=1,MIN(NPNDT,NUMLS)                                           &
     &          )
      ENDIF

      IF (NUMM3 .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'M3EL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'M3EL   ID:',                                                  &
     &          MEMBT(TIMSIM%Memb3(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTMb3(i),                                               &
     &          i=1,MIN(NPNDT,NUMM3)                                           &
     &          )
      ENDIF
      IF (NUMM4 .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'M4EL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'M4EL   ID:',                                                  &
     &          MEMBQ(TIMSIM%Memb4(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTMb4(i),                                               &
     &          i=1,MIN(NPNDT,NUMM4)                                           &
     &          )
      ENDIF

      IF (NUMTR .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'TRUSS Elements With The Smallest Critical Time Steps:',       &
     &          (                                                              &
     &          'TRUSS  ID:',                                                  &
     &          TRUSS(TIMSIM%Truss(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTTru(i),                                               &
     &          i=1,MIN(NPNDT,NUMTR)                                           &
     &          )
      ENDIF

      IF (NUMP3 .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'P3EL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'P3EL   ID:',                                                  &
     &          PLATT(TIMSIM%Plat3(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTPl3(i),                                               &
     &          i=1,MIN(NPNDT,NUMP3)                                           &
     &          )
      ENDIF
      IF (NUMP4 .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'P4EL Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'P4EL   ID:',                                                  &
     &          PLATQ(TIMSIM%Plat4(i))%PAR%EleID,                              &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTPl4(i),                                               &
     &          i=1,MIN(NPNDT,NUMP4)                                           &
     &          )
      ENDIF

      IF (NUMBM .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'BEAM Elements With The Smallest Critical Time Steps:',        &
     &          (                                                              &
     &          'BEAM   ID:',                                                  &
     &          BEAM(TIMSIM%BEAMS(i))%PAR%EleID,                               &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTBms(i),                                               &
     &          i=1,MIN(NPNDT,NUMBM)                                           &
     &          )
      ENDIF

      IF (NUMSP .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &         'SPRING Elements With The Smallest Critical Time Steps:',       &
     &          (                                                              &
     &          'SPRING ID:',                                                  &
     &          SPRING(TIMSIM%Spring(i))%PAR%EleID,                            &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTSpr(i),                                               &
     &          i=1,MIN(NPNDT,NUMSP)                                           &
     &          )
      ENDIF
      IF (NUMDM .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &         'DAMPER Elements With The Smallest Critical Time Steps:',       &
     &          (                                                              &
     &          'DAMPER ID:',                                                  &
     &          DAMPER(TIMSIM%Damper(i))%PAR%EleID,                            &
     &          'Element Critical Time Step:',                                 &
     &          TIMSIM%DTDmp(i),                                               &
     &          i=1,MIN(NPNDT,NUMDM)                                           &
     &          )
      ENDIF

      IF (NUMSC .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'SPRING BC''s With The Smallest Critical Time Steps:',         &
     &          (                                                              &
     &          'SPRING BC ID:',                                               &
     &          SPRING_BC(TIMSIM%SprBC(i))%SprID,                              &
     &          'BC      Critical Time Step:',                                 &
     &          TIMSIM%DTSBC(i),                                               &
     &          i=1,MIN(NPNDT,NUMSC)                                           &
     &          )
      ENDIF
      IF (NUMVC .GT. 0) THEN
        WRITE (IOUNIT,'(/5X,A//(1X,A,I8,5X,A,1PE14.4))')                       &
     &          'DAMPER BC''s With The Smallest Critical Time Steps:',         &
     &          (                                                              &
     &          'DAMPER BC ID:',                                               &
     &          DAMPER_BC(TIMSIM%DmpBC(i))%DprID,                              &
     &          'BC      Critical Time Step:',                                 &
     &          TIMSIM%DTDBC(i),                                               &
     &          i=1,MIN(NPNDT,NUMVC)                                           &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE LEFT_JUSTIFY (STRING,MAX)
!!
!! Left justify string.
!!
      INTEGER                                                                  &
     &          MAX
      CHARACTER                                                                &
     &          STRING*(*)
!!
      I = 1
      DO WHILE (STRING(I:I) .EQ. ' ' .AND. I .LT. MAX)
        I = I + 1
      ENDDO
      STRING(1:MAX) = STRING(I:MAX)
!!
      RETURN
      END
!!_
      SUBROUTINE CUPPER (STRING)
!!
      CHARACTER                                                                &
     &          CHAR*1,                                                        &
     &          STRING*(*)
!!vms   INTEGER*4
!!vms   2       STATUS,
!!vms   2       STR$UPCASE
!!vms   EXTERNAL
!!vms   2       STR$UPCASE
!!
!!vms   STATUS = STR$UPCASE (STRING,STRING)
!!vms   IF (.NOT.STATUS) CALL LIB$SIGNAL (%VAL(STATUS))
!!
      LN = LEN (STRING)
      DO L = 1,LN
        ICX = ICHAR (STRING(L:L))
        IF (ICX .GE. 97 .AND. ICX .LE. 122) THEN
          ICX = ICX - 32
          STRING(L:L) = CHAR(ICX)
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION LAST_NONBLANK (STRING)
!!
!! Copyright (c) by KEY Associates; 25-JAN-1992 09:58:11
!!
!! Purpose: Locate last non-blank character in STRING.
!!
      CHARACTER                                                                &
     &          CHAR*1,                                                        &
     &          NULL*1,                                                        &
     &          STRING*(*)
!!
      NB = 1
      NULL = CHAR(0)
      LN = LEN (STRING)
      DO I = 1,LN
        IF (STRING(I:I) .NE. NULL) THEN
          IF (STRING(I:I) .NE. ' ') THEN
            NB = I
          ENDIF
        ENDIF
      ENDDO
      LAST_NONBLANK = NB
!!
      RETURN
      END
!!_
      SUBROUTINE DIAGNOSTIC_PRINT (LAST_MODULE_CALLED)
!!
!! Copyright (c) by KEY Associates; 25-APR-1993 10:14:20.44
!!
!! Purpose: Print the current value of the counters in common block /BLOCK_03/.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: PAGE_NUMBER = 100
!!
      CHARACTER                                                                &
     &          LAST_MODULE_CALLED*(*)
!!
!! Label output.
!!
      L = LEN (LAST_MODULE_CALLED)
      WRITE (IO_UNIT%LELO,'(/1X,A)')                                           &
     &          'Last Module Called: '//LAST_MODULE_CALLED(1:L)
!!
!! List values for all counters.
!!
      CALL PRINT_COUNTERS                                                      &
     &  (IO_UNIT%LELO,'Number of ', NUMIF, PAGE_NUMBER)
      CALL PRINT_COUNTERS                                                      &
     &  (IO_UNIT%LELO,'Mesh with ', MSHIF, PAGE_NUMBER)
      CALL PRINT_COUNTERS                                                      &
     &  (IO_UNIT%LELO,'Restart w/', NRSIF, PAGE_NUMBER)
!!
      RETURN
      END
!!_
      SUBROUTINE USER_MESSAGE (MESSAGE)
!!
!! Copyright (c) by KEY Associates; 3-MAR-1991 12:01:42
!!
!! Purpose: Output MESSAGE to the print file and continue or terminate
!! depending on severity.
!!
      USE shared_common_data
!!
!! Argument.
      CHARACTER, INTENT(IN) :: MESSAGE*(*) ! I/- Message to be "parsed" and printed
!!
!! Local variables.
      CHARACTER(21)  :: LEADER(4)          ! -/- Severity leader for message.
      CHARACTER(132) :: MSG_TEXT*132       ! -/- Local buffer to left justify output.
      LOGICAL, SAVE  :: FIRST              ! -/- Flag for first line of message
!!
      DATA                                                                     &
     &          LEADER(1) /' *** WARNING MSG *** '/,                           &
     &          LEADER(2) /' *** FATAL ERROR *** '/,                           &
     &          LEADER(3) /' *** INFORMATION *** '/,                           &
     &          LEADER(4) /' *** SEVERITY ?? *** '/
!!
!! Put a blank line (record) in the printed output file.
!!
      WRITE (IO_UNIT%LELO,'()')
!!
!! First, check for valid message (A valid message must have at least one
!! Start-of-New-Line character.)
!!
      IF (INDEX(MESSAGE,MSGL) .EQ. 0) THEN
        LENGTH = MIN (100, LEN(MESSAGE))
        WRITE (IO_UNIT%LELO,200) MESSAGE(1:LENGTH)
        RETURN
      ENDIF
!!
!! Second, parse message into lines. (The minimum message is "MSGL" which
!! will produce nothing.)
!!
      FIRST = .TRUE.
      LM = LEN (MESSAGE)
      IL = INDEX(MESSAGE(1:LM),MSGL)
 100  IR = MIN (LM, IL+INDEX(MESSAGE(MIN(LM,IL+1):LM)//MSGL,MSGL))
      IF (FIRST) THEN
        IF (INDEX(MESSAGE(IL:IR),'WARN') .NE. 0) THEN
          K = 1
        ELSE IF (INDEX(MESSAGE(IL:IR),'FATAL') .NE. 0) THEN
          K = 2
        ELSE IF (INDEX(MESSAGE(IL:IR),'INFORM') .NE. 0) THEN
          K = 3
        ELSE
          K = 4
        ENDIF
        FIRST = .FALSE.
        IF (IR .LT. LM) THEN
          IL = IR
          GO TO 100
        ENDIF
      ELSE
        IF (MESSAGE(IR:IR) .EQ. MSGL) THEN
          MSG_TEXT = LEADER(K)//MESSAGE(IL+1:IR-1)
          WRITE (IO_UNIT%LELO,'(A)') TRIM (MSG_TEXT)
        ELSE
          MSG_TEXT = LEADER(K)//MESSAGE(IL+1:IR)
          WRITE (IO_UNIT%LELO,'(A)') TRIM (MSG_TEXT)
        ENDIF
        IF (IR .LT. LM) THEN
          IL = IR
          GO TO 100
        ENDIF
      ENDIF
!!
!! For fatal error messages, report CPU usage and terminate calculation.
!!
      IF (K .EQ. 2) THEN
!SPEC_CPU2000        CALL TIMER (100)
        PRINT *, ' FMA-3D> Fatal Error Reported. Job Terminated.'
        STOP
      ENDIF
!!
  200   FORMAT                                                                 &
     &    (                                                                    &
     &    ' *** WARNING MSG *** USER_MESSAGE.001.00'/                          &
     &    ' *** WARNING MSG *** Message Impossible To Parse.'/                 &
     &    ' *** WARNING MSG *** Message Has No '//                             &
     &    'Start-of-New-Line Character.'/                                      &
     &    ' *** WARNING MSG *** Message Follows:',A/                           &
     &    ' *** WARNING MSG *** Program Execution Will Continue.'              &
     &    )
!!
      RETURN
      END
!!_
      SUBROUTINE ERRMSG (ROUTINE,MESSAGE,SEVERITY)
!!
!! Copyright (c) by KEY Associates; 31-JUL-1990 18:47:42
!!
!! Purpose: Print error message. Continue or terminate depending on
!! severity. Note: This is an "old" module that can be deleted when all
!! of the Execution_Log_Output messages have converted to USER_MESSAGE.
!!
      USE shared_common_data
!!
      CHARACTER                                                                &
     &          ROUTINE*(*),                                                   &
     &          MESSAGE*(*),                                                   &
     &          SEVERITY*(*)
!!


      IF (SEVERITY(1:3) .EQ. 'WAR') THEN
        WRITE (IO_UNIT%LELO,100) ROUTINE,MESSAGE
        RETURN
      ENDIF
        IF (SEVERITY(1:3) .EQ. 'FAT') THEN
        WRITE (IO_UNIT%LELO,200) ROUTINE,MESSAGE
        PRINT *, 'Fatal error stop from error handler "ERRMSG."'
        STOP
      ELSE
        WRITE (IO_UNIT%LELO,300) ROUTINE,MESSAGE
        PRINT *, 'Fatal error stop from error handler "ERRMSG."'
        STOP
      ENDIF
!!
  100   FORMAT ('0',3X,'*** WARNING *** message from routine: ',               &
     &    A//4X,'<< ',A,' >>'//)
  200   FORMAT ('0',3X,'*** FATAL ERROR *** message from routine: ',           &
     &    A//4X,'<< ',A,' >>'//)
  300   FORMAT ('0',3X,'*** FATAL(?) ERROR *** message from routine: ',        &
     &    A//4X,'<< ',A,' >>'//)
!!
      END
