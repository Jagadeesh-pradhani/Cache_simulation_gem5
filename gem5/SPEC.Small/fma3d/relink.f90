      SUBROUTINE RELINK_USER_IDS
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:28:34
!! Copyright (c) by KEY Associates; 10-OCT-1997 08:39:19.00
!!
      USE shared_common_data
      USE results_
      USE qa_record_
      USE massprop_
      USE gauge1d_
      USE gauge2d_
      USE gauge3d_
      USE material_
      USE layering_
      USE section_2d_
      USE section_1d_
      USE rigid_body_
      USE rigid_body_mass_
      USE nodal_point_mass_
      USE displacement_bc_
      USE tied_bc_
      USE spot_weld_
      USE rigid_wall_bc_
      USE body_force_
      USE pressure_bc_
      USE force_bc_
      USE spring_bc_
      USE damper_bc_
      USE periodic_bc_
      USE nonreflecting_bc_
      USE sliding_interface_
      USE tabulated_function_
      USE node_set_
      USE element_set_
      USE segment_set_
      USE segment_
      USE node_
      USE motion_
      USE force_
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
      USE constrained_node_
      USE plate_pair_
      USE velocity_ic_
!!
      USE enumerated_sets_
      USE output_
      USE relink_scratch_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          NODAL_ID,                                                      &
     &          ID_TABFTN,                                                     &
     &          INTERNAL_ID,                                                   &
     &          RM_INTERNAL_ID,                                                &
     &          MP_INTERNAL_ID,                                                &
     &          NS_INTERNAL_ID,                                                &
     &          ES_INTERNAL_ID,                                                &
     &          SS_INTERNAL_ID,                                                &
     &          MT_INTERNAL_ID,                                                &
     &          TF_INTERNAL_ID,                                                &
     &          LY_INTERNAL_ID,                                                &
     &          S1_INTERNAL_ID,                                                &
     &          S2_INTERNAL_ID,                                                &
     &          SI_INTERNAL_ID
      EXTERNAL                                                                 &
     &          RM_INTERNAL_ID,                                                &
     &          MP_INTERNAL_ID,                                                &
     &          NS_INTERNAL_ID,                                                &
     &          ES_INTERNAL_ID,                                                &
     &          SS_INTERNAL_ID,                                                &
     &          MT_INTERNAL_ID,                                                &
     &          TF_INTERNAL_ID,                                                &
     &          LY_INTERNAL_ID,                                                &
     &          S1_INTERNAL_ID,                                                &
     &          S2_INTERNAL_ID,                                                &
     &          SI_INTERNAL_ID
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered RELINK_USER_IDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Sort external element ID's to facilitate relinking.
!!
      IF (NUMEL .GT. 0) THEN
        n = 0
        DO i = 1,NUMHX
          n = n + 1
          IELV(n) = HEXAH(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMPX
          n = n + 1
          IELV(n) = PENTA(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMTX
          n = n + 1
          IELV(n) = TETRA(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMLS
          n = n + 1
          IELV(n) = LSOLD(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMLX
          n = n + 1
          IELV(n) = LSHEX(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMLM
          n = n + 1
          IELV(n) = LSMBQ(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMM3
          n = n + 1
          IELV(n) = MEMBT(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMM4
          n = n + 1
          IELV(n) = MEMBQ(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMTR
          n = n + 1
          IELV(n) = TRUSS(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMP3
          n = n + 1
          IELV(n) = PLATT(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMP4
          n = n + 1
          IELV(n) = PLATQ(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMBM
          n = n + 1
          IELV(n) = BEAM(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMSP
          n = n + 1
          IELV(n) = SPRING(i)%PAR%EleID
          JELV(n) = i
        ENDDO
        DO i = 1,NUMDM
          n = n + 1
          IELV(n) = DAMPER(i)%PAR%EleID
          JELV(n) = i
        ENDDO
!!
!! Diagnostic write of element ID vector before sorting.
!!
!!
        CALL IVSORT (IELV,JELV,NUMEL,IO_UNIT%LELO)
!!
!! Diagnostic write of element ID vector after sorting.
!!
!!
!! Check to make sure element ID's are unique.
!!
        DO i = 1,NUMEL-1
          IF (IELV(i) .EQ. IELV(i+1))THEN
            WRITE (MSG1,'(I8)') IELV(i+1)
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.001.01'//                               &
     &          MSGL//'Repeated Element ID''s Found:'//MSG1                    &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDDO
      ENDIF
!!
!! Sort external node ID's to facilitate relinking.
!!
      IF (NUMNP .GT. 0) THEN
        DO n = 1,NUMNP
          INPV(n) = NODE(n)%ID
          JNPV(n) = n
        ENDDO
!!
!! Diagnostic write of nodal point ID vector before sorting.
!!
!!
        CALL IVSORT (INPV,JNPV,NUMNP,IO_UNIT%LELO)
!!
!! Diagnostic write of nodal point ID vector after sorting.
!!
!!
!! Check to make sure nodal point ID's are unique.
!!
        DO i = 1,NUMNP-1
          IF (INPV(i) .EQ. INPV(i+1))THEN
            WRITE (MSG1,'(I8)') INPV(i+1)
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.001.02'//                               &
     &          MSGL//'Repeated Nodal Point ID''s Found:'//MSG1                &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDDO
      ENDIF
!!
!! Sort external segment ID's to facilitate relinking.
!!
      IF (NUMSG .GT. 0) THEN
        DO n = 1,NUMSG
          ISGV(n) = SEGMENT(n)%PAR%SegID
          JSGV(n) = n
        ENDDO
!!
!! Diagnostic write of segment ID vector before sorting.
!!
!!
        CALL IVSORT (ISGV,JSGV,NUMSG,IO_UNIT%LELO)
!!
!! Diagnostic write of segment ID vector after sorting.
!!
!!
!! Check to make sure segment ID's are unique.
!!
        DO i = 1,NUMSG-1
          IF (ISGV(i) .EQ. ISGV(i+1))THEN
            WRITE (MSG1,'(I8)') ISGV(i+1)
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.001.03'//                               &
     &          MSGL//'Repeated Segment ID''s Found:'//MSG1                    &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDDO
      ENDIF
!!
!! Record time spent in scanning input records.
!!
!SPEC_CPU2000      CALL TIMER (1)
!!
!! Report total number of repeated node, element, and segment ID's and quit.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'RELINK_USER_IDS.001.04'//                               &
     &          MSGL//'Total Number of Repeated ID''s Found:'//MSG1            &
     &          )
      ENDIF
!!
!! Relink pointers in CONTROL.
!!
!!
      IF (CONTROL%DTTABF .NE. 0) THEN
        INTERNAL_ID = TF_INTERNAL_ID (CONTROL%DTTABF)
        IF (INTERNAL_ID .NE. 0) THEN
          CONTROL%DTTABF = INTERNAL_ID
        ELSE
!!
!! Matching tabulated function ID not found.
!!
          WRITE (MSG1,'(I8)') CONTROL%DTTABF
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.01'//                               &
     &          MSGL//'CONTROL Input Record, Entry: DT_TAB_FTN'//              &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG1                      &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDIF
!!
!! Relink pointers in PRINT.
!!
!!
      IF (PRINT%NODES .GT. 0) THEN
        INTERNAL_ID = NS_INTERNAL_ID (PRINT%NODES)
        IF (INTERNAL_ID .NE. 0) THEN
          PRINT%NODES = INTERNAL_ID
        ELSE
!!
!! Matching node set ID not found.
!!
          WRITE (MSG1,'(I8)') PRINT%NODES
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.02'//                               &
     &          MSGL//'PRINT Input Record Entry NPT'//                         &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG1                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDIF

      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword HXEL,',  0,PRINT%HEXAH )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword PXEL,',  0,PRINT%PENTA )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword TXEL,',  0,PRINT%TETRA )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword M3EL,',  0,PRINT%MEMBT )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword M4EL,',  0,PRINT%MEMBQ )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword TRUSS,', 0,PRINT%TRUSS )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword P3EL,',  0,PRINT%PLATT )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword P4EL,',  0,PRINT%PLATQ )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword BEAM,',  0,PRINT%BEAMS )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword SPRING,',0,PRINT%SPRING)
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword DAMPER,',0,PRINT%DAMPER)
!!
!! Relink pointers in RESULTS.
!!
!!
      DO i = 1,NUMRF
        IF (RESULTS(i)%NODES .GT. 0) THEN
          INTERNAL_ID = NS_INTERNAL_ID (RESULTS(i)%NODES)
          IF (INTERNAL_ID .NE. 0) THEN
            RESULTS(i)%NODES = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') RESULTS(i)%ResID
            WRITE (MSG2,'(I8)') RESULTS(i)%NODES
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.21'//                               &
     &          MSGL//'RESULTS Input Record ID:'//MSG1//                       &
     &          MSGL//'Input Record Parameter: '//OUTPUT%NAME(7)//             &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMRF
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword HXEL,',  RESULTS(i)%ResID,RESULTS(i)%HEXAH )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword PXEL,',  RESULTS(i)%ResID,RESULTS(i)%PENTA )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword TXEL,',  RESULTS(i)%ResID,RESULTS(i)%TETRA )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword LSEL,',  RESULTS(i)%ResID,RESULTS(i)%LSOLD )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword M3EL,',  RESULTS(i)%ResID,RESULTS(i)%MEMBT )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword M4EL,',  RESULTS(i)%ResID,RESULTS(i)%MEMBQ )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword TRUSS,', RESULTS(i)%ResID,RESULTS(i)%TRUSS )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword P3EL,',  RESULTS(i)%ResID,RESULTS(i)%PLATT )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword P4EL,',  RESULTS(i)%ResID,RESULTS(i)%PLATQ )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword BEAM,',  RESULTS(i)%ResID,RESULTS(i)%BEAMS )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword SPRING,',RESULTS(i)%ResID,RESULTS(i)%SPRING)
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword DAMPER,',RESULTS(i)%ResID,RESULTS(i)%DAMPER)
      ENDDO
!!
!! Relink pointers in MASSPROP.
!!
!!
      DO i = 1,NUMMP
        INTERNAL_ID = NS_INTERNAL_ID (MASSPROP(i)%SetID)
        IF (INTERNAL_ID .NE. 0) THEN
          MASSPROP(i)%SetID = INTERNAL_ID
        ELSE
!!
!! Matching node set ID not found.
!!
          WRITE (MSG1,'(I8)') MASSPROP(i)%MPID
          WRITE (MSG2,'(I8)') MASSPROP(i)%SetID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.04'//                               &
     &          MSGL//'MASSPROP Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
!! Relink pointers in GAUGExD.
!!
!!
      DO n = 1,NUMG1
        IF (GAUGE1D(n)%PAR%TBflg .GT. 0) THEN
          CALL RELINK_ELID                                                     &
     &(IELV,JELV,'GAUGE1D',GAUGE1D(n)%PAR%GauID,GAUGE1D(n)%PAR%EleID)
        ELSE
          i = 1
          DO WHILE (GAUGE1D(n)%PAR%IX(i) .GT. 0 .AND. i .LE. 4)
            CALL RELINK_NPID                                                   &
     &(INPV,JNPV,'GAUGE1D',GAUGE1D(n)%PAR%GauID,GAUGE1D(n)%PAR%IX(i))
            i = i + 1
          ENDDO
        ENDIF
      ENDDO
!!
      DO n = 1,NUMG2
        IF (GAUGE2D(n)%PAR%MPQTflg .GT. 0) THEN
          CALL RELINK_ELID                                                     &
     &(IELV,JELV,'GAUGE2D',GAUGE2D(n)%PAR%GauID,GAUGE2D(n)%PAR%EleID)
        ELSE
          i = 1
          DO WHILE ( i .LE. 4 .AND.GAUGE2D(n)%PAR%IX(i+1) .GT. 0)
            i = i + 1
            CALL RELINK_NPID                                                   &
     &(INPV,JNPV,'GAUGE2D',GAUGE2D(n)%PAR%GauID,GAUGE2D(n)%PAR%IX(i))
          ENDDO
        ENDIF
      ENDDO
!!
      DO n = 1,NUMG3
        IF (GAUGE3D(n)%PAR%HPTflg .GT. 0) THEN
          CALL RELINK_ELID                                                     &
     &(IELV,JELV,'GAUGE3D',GAUGE3D(n)%PAR%GauID,GAUGE3D(n)%PAR%EleID)
        ELSE
          i = 1
          DO WHILE (GAUGE3D(n)%PAR%IX(i) .GT. 0 .AND. i .LE. 8)
            CALL RELINK_NPID                                                   &
     &(INPV,JNPV,'GAUGE3D',GAUGE3D(n)%PAR%GauID,GAUGE3D(n)%PAR%IX(i))
            i = i + 1
          ENDDO
        ENDIF
      ENDDO
!!
!! Relink pointers in MATERIAL
!!
!!
      DO n = 1,NUMMT
        IF (MATERIAL(n)%Type .EQ. 3) THEN        ! Multi-linear spring
          MVAL = 6
        ELSE IF (MATERIAL(n)%Type .EQ. 5) THEN   ! User spring, UK1
          MVAL = 6
        ELSE IF (MATERIAL(n)%Type .EQ. 8) THEN   ! Multi-linear damper
          MVAL = 6
        ELSE IF (MATERIAL(n)%Type .EQ. 9) THEN   ! User damper, UD1
          MVAL = 6
        ELSE IF (MATERIAL(n)%Type .EQ. 17) THEN  ! Rubber shear modulus
          MVAL = 9
        ELSE IF (MATERIAL(n)%Type .EQ. 27) THEN  ! Rubber shear modulus
          MVAL = 9
        ELSE IF (MATERIAL(n)%Type .EQ. 37) THEN  ! Rubber shear modulus
          MVAL = 9
        ELSE IF (MATERIAL(n)%Type .EQ. 38) THEN  ! Soil & Crushable Foam
          MVAL = 11
        ELSE IF (MATERIAL(n)%Type .EQ. 47) THEN  ! Rubber shear modulus
          MVAL = 9
        ELSE IF (MATERIAL(n)%Type .EQ. 57) THEN  ! Rubber shear modulus
          MVAL = 9
        ELSE
          MVAL = 0
        ENDIF
        IF (MVAL .NE. 0) THEN
          ID_TABFTN = NINT (MATERIAL(n)%PVAL(MVAL))
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (ID_TABFTN)
          IF (INTERNAL_ID .NE. 0) THEN
            MATERIAL(n)%PVAL(MVAL) = DBLE (INTERNAL_ID)
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') MATERIAL(n)%MatID
            WRITE (MSG2,'(I8)') ID_TABFTN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.33'//                               &
     &          MSGL//'MATERIAL Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in LAYERING.
!!
!!
      DO n = 1,NUMLU
        DO i = 1,LAYERING(n)%Number_of_Layers
          CALL RELINK_MatID                                                    &
     &      ('LAYUP',LAYERING(n)%LupID,LAYERING(n)%MatID(i))
!!
!! Any references to plate/membrane section properties.
!!
          IF (LAYERING(n)%Ltype(i) .GT. 0) THEN
            INTERNAL_ID = S2_INTERNAL_ID (LAYERING(n)%Ltype(i))
            IF (INTERNAL_ID .NE. 0) THEN
              LAYERING(n)%Ltype(i) = INTERNAL_ID
            ELSE
!!
!! Matching section property ID not found.
!!
              WRITE (MSG1,'(I8)') LAYERING(n)%LupID
              WRITE (MSG2,'(I8)') LAYERING(n)%Ltype(i)
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.34'//                               &
     &          MSGL//'LAYUP Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!!
!! Relink pointers in SECTION_1D (Does not yet offer truss/beam
!! user integration rules.
!!
!!!     DO i = 1,NUMS1
!!!       IF (SECTION_1D(i)%Irule .GT. 0) THEN
!!!         INTERNAL_ID =
!!!     2   TF_INTERNAL_ID (SECTION_1D(i)%Irule)
!!!         IF (INTERNAL_ID .NE. 0) THEN
!!!           SECTION_1D(i)%Irule = INTERNAL_ID
!!!         ELSE
!!
!! Matching tabulated function ID not found.
!!
!!!           WRITE (MSG1,'(I8)') SECTION_1D(i)%SecID
!!!           WRITE (MSG2,'(I8)') SECTION_1D(i)%Irule
!!!           CALL USER_MESSAGE
!!!     2       (
!!!     2       MSGL//'WARNING'//
!!!     2       MSGL//'RELINK_USER_IDS.002.05'//
!!!     2       MSGL//'BSECTION Input Record ID:'//MSG1//
!!!     2       MSGL//'Contains Unknown TABFTN ID (Irule):'//MSG2
!!!     2       )
!!!           ERROR%COUNT = ERROR%COUNT + 1
!!!         ENDIF
!!!       ENDIF
!!!     ENDDO
!!
!! Relink pointers in SECTION_2D
!!
!!
      DO i = 1,NUMS2
        IF (SECTION_2D(i)%Irule .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (SECTION_2D(i)%Irule)
          IF (INTERNAL_ID .NE. 0) THEN
            SECTION_2D(i)%Irule = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') SECTION_2D(i)%SecID
            WRITE (MSG2,'(I8)') SECTION_2D(i)%Irule
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.06'//                               &
     &          MSGL//'PSECTION Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown TABFTN ID (Irule):'//MSG2              &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in RIGID_BODY.
!!
!!
      DO i = 1,NUMRB
        IF (RIGID_BODY(i)%Prop .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      RM_INTERNAL_ID (RIGID_BODY(i)%CMID)
          IF (INTERNAL_ID .NE. 0) THEN
            RIGID_BODY(i)%CMID = INTERNAL_ID
          ELSE
!!
!! Matching rigid body mass ID not found.
!!
            WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
            WRITE (MSG2,'(I8)') RIGID_BODY(i)%CMID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.05a'//                              &
     &          MSGL//'RBODY Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown RBMASS ID (mid):'//MSG2                &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in RIGID_BODY_MASS.
!!
      DO i = 1,NUMRM
        CALL RELINK_NPID (INPV,JNPV,'RBMASS',                                  &
     &    RIGID_BODY_MASS(i)%RMID,RIGID_BODY_MASS(i)%NPID)
      ENDDO
!!
!! Relink pointers in NODAL_POINT_MASS.
!!
      DO i = 1,NUMCM
        CALL RELINK_NPID (INPV,JNPV,'NPMASS',                                  &
     &    NODAL_POINT_MASS(i)%NMID,NODAL_POINT_MASS(i)%NPID)
      ENDDO
!!
!! Relink pointers in DISPLACEMENT_BC. (A SetID of zero means all nodes.)
!!
      DO i = 1,NUMDC
        IF (DISPLACEMENT_BC(i)%SetID .LT. 0) THEN
          NODAL_ID = -DISPLACEMENT_BC(i)%SetID
          CALL RELINK_NPID (INPV,JNPV,'DISPBC',                                &
     &      DISPLACEMENT_BC(i)%DBCID,NODAL_ID)
          DISPLACEMENT_BC(i)%SetID = -NODAL_ID
        ELSE IF (DISPLACEMENT_BC(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (DISPLACEMENT_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            DISPLACEMENT_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') DISPLACEMENT_BC(i)%DBCID
            WRITE (MSG2,'(I8)') DISPLACEMENT_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.07'//                               &
     &          MSGL//'DISPBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMDC
        IF (DISPLACEMENT_BC(i)%HstID .NE. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (DISPLACEMENT_BC(i)%HstID)
          IF (INTERNAL_ID .NE. 0) THEN
            DISPLACEMENT_BC(i)%HstID = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') DISPLACEMENT_BC(i)%DBCID
            WRITE (MSG2,'(I8)') DISPLACEMENT_BC(i)%HstID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.08'//                               &
     &          MSGL//'DISPBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in TIED_BC.
!!
      DO i = 1,NUMTC
        INTERNAL_ID = NS_INTERNAL_ID (TIED_BC(i)%SetID)
        IF (INTERNAL_ID .NE. 0) THEN
          TIED_BC(i)%SetID = INTERNAL_ID
        ELSE
!!
!! Matching node set ID not found.
!!
          WRITE (MSG1,'(I8)') TIED_BC(i)%TBCID
          WRITE (MSG2,'(I8)') TIED_BC(i)%SetID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.07'//                               &
     &          MSGL//'TIEDBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
!! Relink pointers in SPOT_WELD.
!!
      DO i = 1,NUMSW
        IF (SPOT_WELD(i)%PLACE .EQ. 'NODES') THEN
          NODAL_ID = SPOT_WELD(i)%NPID(1)
          IF (NODAL_ID .GT. 0) THEN
              CALL RELINK_NPID                                                 &
     &          (INPV,JNPV,'SPOTWELD',SPOT_WELD(i)%SWID,NODAL_ID)
          ENDIF
          SPOT_WELD(i)%NPID(1) = NODAL_ID

          NODAL_ID = SPOT_WELD(i)%NPID(2)
          IF (NODAL_ID .GT. 0) THEN
              CALL RELINK_NPID                                                 &
     &          (INPV,JNPV,'SPOTWELD',SPOT_WELD(i)%SWID,NODAL_ID)
          ENDIF
          SPOT_WELD(i)%NPID(2) = NODAL_ID

          NODAL_ID = SPOT_WELD(i)%NPID(3)
          IF (NODAL_ID .GT. 0) THEN
              CALL RELINK_NPID                                                 &
     &          (INPV,JNPV,'SPOTWELD',SPOT_WELD(i)%SWID,NODAL_ID)
          ENDIF
          SPOT_WELD(i)%NPID(3) = NODAL_ID
        ENDIF
      ENDDO
!!
!! Relink pointers in RIGID_WALL_BC. (A SetID of zero means all nodes.)
!!
      DO i = 1,NUMWC
        IF (RIGID_WALL_BC(i)%SetID .LT. 0) THEN
          NODAL_ID = -RIGID_WALL_BC(i)%SetID
          CALL RELINK_NPID (INPV,JNPV,'WALLBC',                                &
     &      RIGID_WALL_BC(i)%RWID,NODAL_ID)
          RIGID_WALL_BC(i)%SetID = -NODAL_ID
        ELSE IF (RIGID_WALL_BC(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (RIGID_WALL_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            RIGID_WALL_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') RIGID_WALL_BC(i)%RWID
            WRITE (MSG2,'(I8)') RIGID_WALL_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.09'//                               &
     &          MSGL//'WALLBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMWC
        IF (RIGID_WALL_BC(i)%HstID .NE. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (RIGID_WALL_BC(i)%HstID)
          IF (INTERNAL_ID .NE. 0) THEN
            RIGID_WALL_BC(i)%HstID = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') RIGID_WALL_BC(i)%RWID
            WRITE (MSG2,'(I8)') RIGID_WALL_BC(i)%HstID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.10'//                               &
     &          MSGL//'WALLBC Input Record ID:'//MSG1//                        &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in BODY_FORCE.
!!
      DO i = 1,NUMBF
        IF (BODY_FORCE(i)%SetID .LT. 0) THEN
          INTERNAL_ID = -BODY_FORCE(i)%SetID
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'BODYFORCE',BODY_FORCE(i)%BFID,INTERNAL_ID)
          BODY_FORCE(i)%SetID = -INTERNAL_ID
        ELSE IF (BODY_FORCE(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (BODY_FORCE(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            BODY_FORCE(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching nodal point set ID not found.
!!
            WRITE (MSG1,'(I8)') BODY_FORCE(i)%BFID
            WRITE (MSG2,'(I8)') BODY_FORCE(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.11'//                               &
     &          MSGL//'BODYFORCE Input Record ID:'//MSG1//                     &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMBF
        IF (BODY_FORCE(i)%HstID .NE. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (BODY_FORCE(i)%HstID)
          IF (INTERNAL_ID .NE. 0) THEN
            BODY_FORCE(i)%HstID = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') BODY_FORCE(i)%BFID
            WRITE (MSG2,'(I8)') BODY_FORCE(i)%HstID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.12'//                               &
     &          MSGL//'BODYFORCE Input Record ID:'//MSG1//                     &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in PRESSURE_BC.
!!
      DO i = 1,NUMPC
        IF (PRESSURE_BC(i)%SetID .LT. 0) THEN
          INTERNAL_ID = -PRESSURE_BC(i)%SetID
          CALL RELINK_SGID                                                     &
     &      (ISGV,JSGV,'PRESSBC',PRESSURE_BC(i)%PBCID,INTERNAL_ID)
          PRESSURE_BC(i)%SetID = -INTERNAL_ID
        ELSE IF (PRESSURE_BC(i)%SetID .GE. 0) THEN
          INTERNAL_ID =                                                        &
     &      SS_INTERNAL_ID (PRESSURE_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            PRESSURE_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching segment set ID not found.
!!
            WRITE (MSG1,'(I8)') PRESSURE_BC(i)%PBCID
            WRITE (MSG2,'(I8)') PRESSURE_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.13'//                               &
     &          MSGL//'PRESSBC Input Record ID:'//MSG1//                       &
     &          MSGL//'Contains Unknown SEGMENT SET ID:'//MSG2                 &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMPC
        IF (PRESSURE_BC(i)%HstID .NE. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (PRESSURE_BC(i)%HstID)
          IF (INTERNAL_ID .NE. 0) THEN
            PRESSURE_BC(i)%HstID = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') PRESSURE_BC(i)%PBCID
            WRITE (MSG2,'(I8)') PRESSURE_BC(i)%HstID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.14'//                               &
     &          MSGL//'PRESSBC Input Record ID:'//MSG1//                       &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in FORCE_BC. (A set ID of zero means all nodes.)
!!
      DO i = 1,NUMFC
        IF (FORCE_BC(i)%SetID .LT. 0) THEN
          NODAL_ID = -FORCE_BC(i)%SetID
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'FORCEBC',FORCE_BC(i)%CFID,NODAL_ID)
          FORCE_BC(i)%SetID = -NODAL_ID
        ELSE IF (FORCE_BC(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (FORCE_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            FORCE_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') FORCE_BC(i)%CFID
            WRITE (MSG2,'(I8)') FORCE_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.15'//                               &
     &          MSGL//'FORCEBC Input Record ID:'//MSG1//                       &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMFC
        IF (FORCE_BC(i)%HstID .NE. 0) THEN
          INTERNAL_ID =                                                        &
     &      TF_INTERNAL_ID (FORCE_BC(i)%HstID)
          IF (INTERNAL_ID .NE. 0) THEN
            FORCE_BC(i)%HstID = INTERNAL_ID
          ELSE
!!
!! Matching tabulated function ID not found.
!!
            WRITE (MSG1,'(I8)') FORCE_BC(i)%CFID
            WRITE (MSG2,'(I8)') FORCE_BC(i)%HstID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.16'//                               &
     &          MSGL//'FORCEBC Input Record ID:'//MSG1//                       &
     &          MSGL//'Contains Unknown TABFTN ID:'//MSG2                      &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in SPRING_BC. (A set ID of zero means all nodes.)
!!
      DO i = 1,NUMSC
        IF (SPRING_BC(i)%SetID .LT. 0) THEN
          NODAL_ID = -SPRING_BC(i)%SetID
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'SPRINGBC',SPRING_BC(i)%SprID,NODAL_ID)
          SPRING_BC(i)%SetID = -NODAL_ID
        ELSE IF (SPRING_BC(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (SPRING_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            SPRING_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') SPRING_BC(i)%SprID
            WRITE (MSG2,'(I8)') SPRING_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.17'//                               &
     &          MSGL//'SPRINGBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO n = 1,NUMSC
        CALL RELINK_MatID                                                      &
     &    ('SPRINGBC',SPRING_BC(n)%SprID,SPRING_BC(n)%MatID)
      ENDDO
!!
!! Relink pointers in DAMPER_BC. (A set ID of zero means all nodes.)
!!
      DO i = 1,NUMVC
        IF (DAMPER_BC(i)%SetID .LT. 0) THEN
          NODAL_ID = -DAMPER_BC(i)%SetID
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'DAMPERBC',DAMPER_BC(i)%DprID,NODAL_ID)
          DAMPER_BC(i)%SetID = -NODAL_ID
        ELSE IF (DAMPER_BC(i)%SetID .GT. 0) THEN
          INTERNAL_ID =                                                        &
     &      NS_INTERNAL_ID (DAMPER_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            DAMPER_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') DAMPER_BC(i)%DprID
            WRITE (MSG2,'(I8)') DAMPER_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.18'//                               &
     &          MSGL//'DAMPERBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO n = 1,NUMVC
        CALL RELINK_MatID                                                      &
     &    ('DAMPERBC',DAMPER_BC(n)%DprID,DAMPER_BC(n)%MatID)
      ENDDO
!!
!! Relink pointers in PERIODIC_BC.
!!
      DO i = 1,NUMCC
!!
!! Side 1 of the periodic BC.
!!
        IF (PERIODIC_BC(i)%Typ1 .EQ. 0) THEN
          IF (PERIODIC_BC(i)%S1ID .LT. 0) THEN
            INTERNAL_ID = -PERIODIC_BC(i)%S1ID
            CALL RELINK_SGID                                                   &
     &        (ISGV,JSGV,'PERIODBC',PERIODIC_BC(i)%PerID,INTERNAL_ID)
            PERIODIC_BC(i)%S1ID = -INTERNAL_ID
          ELSE IF (PERIODIC_BC(i)%S1ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        SS_INTERNAL_ID (PERIODIC_BC(i)%S1ID)
            IF (INTERNAL_ID .NE. 0) THEN
              PERIODIC_BC(i)%S1ID = INTERNAL_ID
            ELSE
!!
!! Matching segment set ID not found.
!!
              WRITE (MSG1,'(I8)') PERIODIC_BC(i)%PerID
              WRITE (MSG2,'(I8)') PERIODIC_BC(i)%S1ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.20'//                               &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown Side 1 SEGSET ID:'//MSG2               &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ELSE IF (PERIODIC_BC(i)%Typ1 .EQ. 1) THEN
          IF (PERIODIC_BC(i)%S1ID .LT. 0) THEN
            INTERNAL_ID = -PERIODIC_BC(i)%S1ID
            CALL RELINK_NPID                                                   &
     &        (INPV,JNPV,'PERIODBC',PERIODIC_BC(i)%PerID,INTERNAL_ID)
            PERIODIC_BC(i)%S1ID = -INTERNAL_ID
          ELSE IF (PERIODIC_BC(i)%S1ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        NS_INTERNAL_ID (PERIODIC_BC(i)%S1ID)
            IF (INTERNAL_ID .NE. 0) THEN
              PERIODIC_BC(i)%S1ID = INTERNAL_ID
            ELSE
!!
!! Matching nodal point set ID not found.
!!
              WRITE (MSG1,'(I8)') PERIODIC_BC(i)%PerID
              WRITE (MSG2,'(I8)') PERIODIC_BC(i)%S1ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.21'//                               &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown Side 1 NPSET ID:'//MSG2                &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ENDIF
!!
!! Side 2 of the periodic BC.
!!
        IF (PERIODIC_BC(i)%Typ2 .EQ. 0) THEN
          IF (PERIODIC_BC(i)%S2ID .LT. 0) THEN
            INTERNAL_ID = -PERIODIC_BC(i)%S2ID
            CALL RELINK_SGID                                                   &
     &        (ISGV,JSGV,'PERIODBC',PERIODIC_BC(i)%PerID,INTERNAL_ID)
            PERIODIC_BC(i)%S2ID = -INTERNAL_ID
          ELSE IF (PERIODIC_BC(i)%S2ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        SS_INTERNAL_ID (PERIODIC_BC(i)%S2ID)
            IF (INTERNAL_ID .NE. 0) THEN
              PERIODIC_BC(i)%S2ID = INTERNAL_ID
            ELSE
!!
!! Matching segment set ID not found.
!!
              WRITE (MSG1,'(I8)') PERIODIC_BC(i)%PerID
              WRITE (MSG2,'(I8)') PERIODIC_BC(i)%S2ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.22'//                               &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown Side 2 SEGSET ID:'//MSG2               &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ELSE IF (PERIODIC_BC(i)%Typ2 .EQ. 1) THEN
          IF (PERIODIC_BC(i)%S2ID .LT. 0) THEN
            INTERNAL_ID = -PERIODIC_BC(i)%S2ID
            CALL RELINK_NPID                                                   &
     &        (INPV,JNPV,'PERIODBC',PERIODIC_BC(i)%PerID,INTERNAL_ID)
            PERIODIC_BC(i)%S2ID = -INTERNAL_ID
          ELSE IF (PERIODIC_BC(i)%S2ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        NS_INTERNAL_ID (PERIODIC_BC(i)%S2ID)
            IF (INTERNAL_ID .NE. 0) THEN
              PERIODIC_BC(i)%S2ID = INTERNAL_ID
            ELSE
!!
!! Matching nodal point set ID not found.
!!
              WRITE (MSG1,'(I8)') PERIODIC_BC(i)%PerID
              WRITE (MSG2,'(I8)') PERIODIC_BC(i)%S2ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.23'//                               &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown Side 2 NPSET ID:'//MSG2                &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in NONREFLECTING_BC.
!!
      DO i = 1,NUMNR
        IF (NONREFLECTING_BC(i)%SetID .LT. 0) THEN
          INTERNAL_ID = -NONREFLECTING_BC(i)%SetID
          CALL RELINK_SGID                                                     &
     &      (ISGV,JSGV,'NRBC',NONREFLECTING_BC(i)%NRID,INTERNAL_ID)
          NONREFLECTING_BC(i)%SetID = -INTERNAL_ID
        ELSE
          INTERNAL_ID =                                                        &
     &      SS_INTERNAL_ID (NONREFLECTING_BC(i)%SetID)
          IF (INTERNAL_ID .NE. 0) THEN
            NONREFLECTING_BC(i)%SetID = INTERNAL_ID
          ELSE
!!
!! Matching segment set ID not found.
!!
            WRITE (MSG1,'(I8)') NONREFLECTING_BC(i)%NRID
            WRITE (MSG2,'(I8)') NONREFLECTING_BC(i)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.19'//                               &
     &          MSGL//'NRBC Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown SEGMENT SET ID:'//MSG2                 &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in INTERFACE_TIME.
!!
      DO n = 1,NUMIT
        INTERNAL_ID = SI_INTERNAL_ID (INTERFACE_TIME%SI(n)%ID)
        IF (INTERNAL_ID .NE. 0) THEN
          INTERFACE_TIME%SI(n)%ID = INTERNAL_ID
        ELSE
!!
!! Matching sliding interface ID not found.
!!
          WRITE (MSG2,'(I8)') INTERFACE_TIME%SI(n)%ID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.0002.12'//                              &
     &          MSGL//'An INTERFACE Input Record References'//                 &
     &          MSGL//'Unknown SLIDING_INTERFACE ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
!! Relink pointers in SLIDING_INTERFACE.
!!
      DO i = 1,NUMSI
!!
!! Side 1 of the sliding interface.
!!
        IF (SLIDING_INTERFACE(i)%Typ1 .EQ. 0) THEN
          IF (SLIDING_INTERFACE(i)%S1ID .LT. 0) THEN
            INTERNAL_ID = -SLIDING_INTERFACE(i)%S1ID
            CALL RELINK_SGID                                                   &
     &        (ISGV,JSGV,'SLIDE',SLIDING_INTERFACE(i)%SIID,INTERNAL_ID)
            SLIDING_INTERFACE(i)%S1ID = -INTERNAL_ID
          ELSE IF (SLIDING_INTERFACE(i)%S1ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        SS_INTERNAL_ID (SLIDING_INTERFACE(i)%S1ID)
            IF (INTERNAL_ID .NE. 0) THEN
              SLIDING_INTERFACE(i)%S1ID = INTERNAL_ID
            ELSE
!!
!! Matching segment set ID not found.
!!
              WRITE (MSG1,'(I8)') SLIDING_INTERFACE(i)%SIID
              WRITE (MSG2,'(I8)') SLIDING_INTERFACE(i)%S1ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.20'//                               &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown Side 1 SEGSET ID:'//MSG2               &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ELSE IF (SLIDING_INTERFACE(i)%Typ1 .EQ. 1) THEN
          IF (SLIDING_INTERFACE(i)%S1ID .LT. 0) THEN
            INTERNAL_ID = -SLIDING_INTERFACE(i)%S1ID
            CALL RELINK_NPID                                                   &
     &        (INPV,JNPV,'SLIDE',SLIDING_INTERFACE(i)%SIID,INTERNAL_ID)
            SLIDING_INTERFACE(i)%S1ID = -INTERNAL_ID
          ELSE IF (SLIDING_INTERFACE(i)%S1ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        NS_INTERNAL_ID (SLIDING_INTERFACE(i)%S1ID)
            IF (INTERNAL_ID .NE. 0) THEN
              SLIDING_INTERFACE(i)%S1ID = INTERNAL_ID
            ELSE
!!
!! Matching nodal point set ID not found.
!!
              WRITE (MSG1,'(I8)') SLIDING_INTERFACE(i)%SIID
              WRITE (MSG2,'(I8)') SLIDING_INTERFACE(i)%S1ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.21'//                               &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown Side 1 NPSET ID:'//MSG2                &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ENDIF
!!
!! Side 2 of the sliding interface.
!!
        IF (SLIDING_INTERFACE(i)%Typ2 .EQ. 0) THEN
          IF (SLIDING_INTERFACE(i)%S2ID .LT. 0) THEN
            INTERNAL_ID = -SLIDING_INTERFACE(i)%S2ID
            CALL RELINK_SGID                                                   &
     &        (ISGV,JSGV,'SLIDE',SLIDING_INTERFACE(i)%SIID,INTERNAL_ID)
            SLIDING_INTERFACE(i)%S2ID = -INTERNAL_ID
          ELSE IF (SLIDING_INTERFACE(i)%S2ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        SS_INTERNAL_ID (SLIDING_INTERFACE(i)%S2ID)
            IF (INTERNAL_ID .NE. 0) THEN
              SLIDING_INTERFACE(i)%S2ID = INTERNAL_ID
            ELSE
!!
!! Matching segment set ID not found.
!!
              WRITE (MSG1,'(I8)') SLIDING_INTERFACE(i)%SIID
              WRITE (MSG2,'(I8)') SLIDING_INTERFACE(i)%S2ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.22'//                               &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown Side 2 SEGSET ID:'//MSG2               &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ELSE IF (SLIDING_INTERFACE(i)%Typ2 .EQ. 1) THEN
          IF (SLIDING_INTERFACE(i)%S2ID .LT. 0) THEN
            INTERNAL_ID = -SLIDING_INTERFACE(i)%S2ID
            CALL RELINK_NPID                                                   &
     &        (INPV,JNPV,'SLIDE',SLIDING_INTERFACE(i)%SIID,INTERNAL_ID)
            SLIDING_INTERFACE(i)%S2ID = -INTERNAL_ID
          ELSE IF (SLIDING_INTERFACE(i)%S2ID .GT. 0) THEN
            INTERNAL_ID =                                                      &
     &        NS_INTERNAL_ID (SLIDING_INTERFACE(i)%S2ID)
            IF (INTERNAL_ID .NE. 0) THEN
              SLIDING_INTERFACE(i)%S2ID = INTERNAL_ID
            ELSE
!!
!! Matching nodal point set ID not found.
!!
              WRITE (MSG1,'(I8)') SLIDING_INTERFACE(i)%SIID
              WRITE (MSG2,'(I8)') SLIDING_INTERFACE(i)%S2ID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.23'//                               &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown Side 2 NPSET ID:'//MSG2                &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! Relink pointers in HEXAH.
!!
      DO n = 1,NUMHX
        CALL RELINK_MatID                                                      &
     &    ('HXEL',HEXAH(n)%PAR%EleID,HEXAH(n)%PAR%MatID)
!!
        IF (HEXAH(n)%PAR%LupID .GE. 0) THEN
          CALL RELINK_LupID                                                    &
     &      ('HXEL',HEXAH(n)%PAR%EleID,HEXAH(n)%PAR%LupID)
        ENDIF
!!
        DO i = 1,8
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'HXEL',HEXAH(n)%PAR%EleID,HEXAH(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in PENTA.
!!
      DO n = 1,NUMPX
        CALL RELINK_MatID                                                      &
     &    ('PXEL',PENTA(n)%PAR%EleID,PENTA(n)%PAR%MatID)
!!
        IF (PENTA(n)%PAR%LupID .GE. 0) THEN
          CALL RELINK_LupID                                                    &
     &      ('PXEL',PENTA(n)%PAR%EleID,PENTA(n)%PAR%LupID)
        ENDIF
!!
        DO i = 1,6
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'PXEL',PENTA(n)%PAR%EleID,PENTA(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in TETRA.
!!
      DO n = 1,NUMTX
        CALL RELINK_MatID                                                      &
     &    ('TXEL',TETRA(n)%PAR%EleID,TETRA(n)%PAR%MatID)
!!
        IF (TETRA(n)%PAR%LupID .GE. 0) THEN
          CALL RELINK_LupID                                                    &
     &      ('TXEL',TETRA(n)%PAR%EleID,TETRA(n)%PAR%LupID)
        ENDIF
!!
        DO i = 1,4
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'TXEL',TETRA(n)%PAR%EleID,TETRA(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in LSOLD.
!!
      DO n = 1,NUMLS
        CALL RELINK_LupID                                                      &
     &    ('LSEL',LSOLD(n)%PAR%EleID,LSOLD(n)%PAR%LupID)
!!
        DO i = 1,8
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'LSEL',LSOLD(n)%PAR%EleID,LSOLD(n)%PAR%IX(i))
        ENDDO
!!
        DO i = 1,LAYERING(LSOLD(n)%PAR%LupID)%Number_of_Layers
          CALL RELINK_ELID                                                     &
     &      (IELV,JELV,'LSEL',LSOLD(n)%PAR%EleID,LSOLD(n)%PAR%ID(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in LSHEX.
!!
      DO n = 1,NUMLX
        CALL RELINK_MatID                                                      &
     &    ('LSEL Sub Hex',LSHEX(n)%PAR%EleID,LSHEX(n)%PAR%MatID)
      ENDDO
!!
!! Relink pointers in LSMBQ.
!!
      DO n = 1,NUMLM
        CALL RELINK_MatID                                                      &
     &    ('LSEL Sub MBQ',LSMBQ(n)%PAR%EleID,LSMBQ(n)%PAR%MatID)
!!
        INTERNAL_ID = S2_INTERNAL_ID (LSMBQ(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          LSMBQ(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') LSMBQ(n)%PAR%EleID
          WRITE (MSG2,'(I8)') LSMBQ(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.25'//                               &
     &          MSGL//'LSEL Sub-Membrane ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
!! Relink pointers in MEMBQ.
!!
      DO n = 1,NUMM4
        CALL RELINK_MatID                                                      &
     &    ('M4EL',MEMBQ(n)%PAR%EleID,MEMBQ(n)%PAR%MatID)
      ENDDO
!!
      DO n = 1,NUMM4
        INTERNAL_ID = S2_INTERNAL_ID (MEMBQ(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          MEMBQ(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') MEMBQ(n)%PAR%EleID
          WRITE (MSG2,'(I8)') MEMBQ(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.26'//                               &
     &          MSGL//'M4EL Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
      DO n = 1,NUMM4
        DO i = 1,4
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'M4EL',MEMBQ(n)%PAR%EleID,MEMBQ(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in MEMBT.
!!
      DO n = 1,NUMM3
        CALL RELINK_MatID                                                      &
     &    ('M3EL',MEMBT(n)%PAR%EleID,MEMBT(n)%PAR%MatID)
!!
        INTERNAL_ID = S2_INTERNAL_ID (MEMBT(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          MEMBT(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') MEMBT(n)%PAR%EleID
          WRITE (MSG2,'(I8)') MEMBT(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.27'//                               &
     &          MSGL//'M3EL Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        DO i = 1,3
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'M3EL',MEMBT(n)%PAR%EleID,MEMBT(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in TRUSS.
!!
      DO n = 1,NUMTR
        CALL RELINK_MatID                                                      &
     &    ('TRUSS',TRUSS(n)%PAR%EleID,TRUSS(n)%PAR%MatID)
!!
        INTERNAL_ID = S1_INTERNAL_ID (TRUSS(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          TRUSS(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') TRUSS(n)%PAR%EleID
          WRITE (MSG2,'(I8)') TRUSS(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.28'//                               &
     &          MSGL//'TRUSS Input Record ID:'//MSG1//                         &
     &          MSGL//'Contains Unknown BSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        DO i = 1,2
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'TRUSS',TRUSS(n)%PAR%EleID,TRUSS(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in PLATQ.
!!
      DO n = 1,NUMP4
        CALL RELINK_MatID                                                      &
     &    ('P4EL',PLATQ(n)%PAR%EleID,PLATQ(n)%PAR%MatID)
!!
        INTERNAL_ID = S2_INTERNAL_ID (PLATQ(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          PLATQ(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') PLATQ(n)%PAR%EleID
          WRITE (MSG2,'(I8)') PLATQ(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.29'//                               &
     &          MSGL//'P4EL Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        DO i = 1,4
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'P4EL',PLATQ(n)%PAR%EleID,PLATQ(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in PLATT.
!!
      DO n = 1,NUMP3
        CALL RELINK_MatID                                                      &
     &    ('P3EL',PLATT(n)%PAR%EleID,PLATT(n)%PAR%MatID)
!!
        INTERNAL_ID = S2_INTERNAL_ID (PLATT(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          PLATT(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') PLATT(n)%PAR%EleID
          WRITE (MSG2,'(I8)') PLATT(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.30'//                               &
     &          MSGL//'P3EL Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown PSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        DO i = 1,3
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'P3EL',PLATT(n)%PAR%EleID,PLATT(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in BEAM.
!!
      DO n = 1,NUMBM
        CALL RELINK_MatID                                                      &
     &    ('BEAM',BEAM(n)%PAR%EleID,BEAM(n)%PAR%MatID)
!!
        INTERNAL_ID = S1_INTERNAL_ID (BEAM(n)%PAR%SecID)
        IF (INTERNAL_ID .NE. 0) THEN
          BEAM(n)%PAR%SecID = INTERNAL_ID
        ELSE
!!
!! Matching section property ID not found.
!!
          WRITE (MSG1,'(I8)') BEAM(n)%PAR%EleID
          WRITE (MSG2,'(I8)') BEAM(n)%PAR%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.31'//                               &
     &          MSGL//'BEAM Input Record ID:'//MSG1//                          &
     &          MSGL//'Contains Unknown BSECTION ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        DO i = 1,2
          CALL RELINK_NPID                                                     &
     &      (INPV,JNPV,'BEAM',BEAM(n)%PAR%EleID,BEAM(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in SPRING.
!!
      DO n = 1,NUMSP
        CALL RELINK_MatID                                                      &
     &    ('SPRING',SPRING(n)%PAR%EleID,SPRING(n)%PAR%MatID)
!!
        DO i = 1,2
          CALL RELINK_NPID                                                     &
     &(INPV,JNPV,'SPRING',SPRING(n)%PAR%EleID,SPRING(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in DAMPER.
!!
      DO n = 1,NUMDM
        CALL RELINK_MatID                                                      &
     &    ('DAMPER',DAMPER(n)%PAR%EleID,DAMPER(n)%PAR%MatID)
!!
        DO i = 1,2
          CALL RELINK_NPID                                                     &
     &(INPV,JNPV,'DAMPER',DAMPER(n)%PAR%EleID,DAMPER(n)%PAR%IX(i))
        ENDDO
      ENDDO
!!
!! Relink pointers in SEGMENT.
!!
      DO n = 1,NUMSG
        DO i = 1,3
          CALL RELINK_NPID                                                     &
     &(INPV,JNPV,'SEGMENT',SEGMENT(n)%PAR%SegID,SEGMENT(n)%PAR%IX(i))
        ENDDO
        IF (SEGMENT(n)%PAR%IX(4) .NE. 0) THEN
          CALL RELINK_NPID                                                     &
     &(INPV,JNPV,'SEGMENT',SEGMENT(n)%PAR%SegID,SEGMENT(n)%PAR%IX(4))
        ENDIF
      ENDDO
!!
!! Relink pointers in NODE_SET.
!!
      DO n = 1,NUMNS
        IF (NODE_SET(n)%Flag .NE. 'ALL') THEN
          DO i = NODE_SET(n)%Istart,NODE_SET(n)%Iend
            CALL RELINK_NPID                                                   &
     &        (INPV,JNPV,'NPSET',NODE_SET(n)%SetID,NNPSETS(i))
          ENDDO
        ENDIF
      ENDDO
!!
!! Relink pointers in ELEMENT_SET.
!!
      DO n = 1,NUMES
        IF (ELEMENT_SET(n)%Flag .NE. 'ALL') THEN
          DO i = ELEMENT_SET(n)%Istart,ELEMENT_SET(n)%Iend
            CALL RELINK_ELID                                                   &
     &        (IELV,JELV,'ELSET',ELEMENT_SET(n)%SetID,NELSETS(i))
          ENDDO
        ENDIF
      ENDDO
!!
!! Relink pointers in SEGMENT_SET.
!!
      DO n = 1,NUMSS
        IF (SEGMENT_SET(n)%Flag .NE. 'ALL') THEN
          DO i = SEGMENT_SET(n)%Istart,SEGMENT_SET(n)%Iend
            CALL RELINK_SGID                                                   &
     &        (ISGV,JSGV,'SEGSET',SEGMENT_SET(n)%SetID,NSGSETS(i))
          ENDDO
        ENDIF
      ENDDO
!!
!! Relink pointers in CONSTRAINED_NODE.
!!
      DO n = 1,NUMNC
        CALL RELINK_NPID                                                       &
     &    (INPV,JNPV,'NPCON1',CONSTRAINED_NODE(n)%ID,                          &
     &    CONSTRAINED_NODE(n)%CNID)
        CALL RELINK_NPID                                                       &
     &    (INPV,JNPV,'NPCON1',CONSTRAINED_NODE(n)%ID,                          &
     &    CONSTRAINED_NODE(n)%NPID(1))
        CALL RELINK_NPID                                                       &
     &    (INPV,JNPV,'NPCON1',CONSTRAINED_NODE(n)%ID,                          &
     &    CONSTRAINED_NODE(n)%NPID(2))
      ENDDO
!!
!! Relink pointers in VELOCITY_IC.
!!
      DO n = 1,NUMIC
        CALL RELINK_NPID                                                       &
     &    (INPV,JNPV,'VELIC',VELOCITY_IC(n)%ICID,VELOCITY_IC(n)%NPID)
      ENDDO
!!
!! Add time spent in relinking input records to time spent in
!! scanning, reading, and processing input records.
!!
!SPEC_CPU2000      CALL TIMER (1)
!!
      IF (ERROR%COUNT .GT. 0) THEN
!!
!! Print total number of relinking errors detected and quit.
!!
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_USER_IDS.002.32'//                               &
     &          MSGL//'Total Number of Relinking Errors:'//MSG1//              &
     &          MSGL//'Please Check ID''s Very Carefully.'                     &
     &          )
!SPEC_CPU2000        CALL TIMER (100)
      STOP ' FMA-3D/RELINK> Abnormal Termination; Relinking Errors.'
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_ELID                                                   &
     &          (IELV,JELV,RECORD_KW,RECORD_ID,ELEMENT_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1990 22:24:37
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          ELID,           & ! -/- Convert 2-byte integer to 4-byte 
     &          IELV(1),        & ! I/- Sorted list of element ID's      
     &          JELV(1),        & ! I/- Corresponding internal elemt ID's
     &          INTFEXT           ! I/- External integer function
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record key word
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch 
     &          ELEMENT_ID        ! I/O ELEMENT ID, relinked
!!
!! Look for internal ID. (The internal ID is the value in JELV(1:NUMEL)
!! "opposite" the external element ID in the sequence IELV(1:NUMEL).
!! The user specified ELEMENT_ID is located in IELV(1:NUMEL) by INTFEXT.
!!
      ELID = ELEMENT_ID
      INTERNAL_ID = INTFEXT (IELV,JELV,NUMEL,ELID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        ELEMENT_ID = INTERNAL_ID
      ELSE
!!
!! Input element ID not found among the ID's in IELV(1:NUMEL).
!!
        L = LEN (RECORD_KW)
        WRITE (MSG1,'(I8)') RECORD_ID
        WRITE (MSG2,'(I8)') ELEMENT_ID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_ELID.001.00'//                                   &
     &          MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//              &
     &          MSGL//'Contains Unknown Element ID:'//MSG2                     &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_NPID                                                   &
     &          (INPV,JNPV,RECORD_KW,RECORD_ID,NODAL_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1990 22:24:37
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          NPID,           & ! -/- Convert 2-byte integer to 4-byte 
     &          INPV(1),        & ! I/- Sorted list of nodal point ID's  
     &          JNPV(1),        & ! I/- Corresponding internal nodal ID's
     &          INTFEXT           ! I/- External integer function
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record key word
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch 
     &          NODAL_ID          ! I/O NODAL ID, relinked
!!
!! Look for internal ID. (The internal ID is the value in JNPV(1:NUMNP)
!! "opposite" the external nodal point ID in the sequence INPV(1:NUMNP).
!! The user specified NODAL_ID is located in INPV(1:NUMNP) by INTFEXT.
!!
      NPID = NODAL_ID
      INTERNAL_ID = INTFEXT (INPV,JNPV,NUMNP,NPID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        NODAL_ID = INTERNAL_ID
      ELSE
!!
!! Input nodal point ID not found among the ID's in INPV(1:NUMNP).
!!
        L = LEN (RECORD_KW)
        WRITE (MSG1,'(I8)') RECORD_ID
        WRITE (MSG2,'(I8)') NODAL_ID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_NPID.001.00'//                                   &
     &          MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//              &
     &          MSGL//'Contains Unknown Nodal Point ID:'//MSG2                 &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_SGID                                                   &
     &          (ISGV,JSGV,RECORD_KW,RECORD_ID,SEGMENT_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1990 22:24:37
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          SGID,           & ! -/- Convert 2-byte integer to 4-byte 
     &          ISGV(1),        & ! I/- Sorted list of segment ID's      
     &          JSGV(1),        & ! I/- Corresponding internal sgmt ID's 
     &          INTFEXT           ! I/- External integer function
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record key word
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch 
     &          SEGMENT_ID        ! I/O SEGMENT ID, relinked
!!
!! Look for internal ID. (The internal ID is the value in JSGV(1:NUMSG)
!! "opposite" the external segment ID in the sequence ISGV(1:NUMSG).
!! The user specified SEGMENT_ID is located in ISGV(1:NUMSG) by INTFEXT.
!!
      SGID = SEGMENT_ID
      INTERNAL_ID = INTFEXT (ISGV,JSGV,NUMSG,SGID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        SEGMENT_ID = INTERNAL_ID
      ELSE
!!
!! Input segment ID not found among the ID's in ISGV(1:NUMSG).
!!
        L = LEN (RECORD_KW)
        WRITE (MSG1,'(I8)') RECORD_ID
        WRITE (MSG2,'(I8)') SEGMENT_ID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_SGID.001.00'//                                   &
     &          MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//              &
     &          MSGL//'Contains Unknown SEGMENT ID:'//MSG2                     &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION RM_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:28:51
!!
      USE shared_common_data
      USE rigid_body_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMRM-1
          ID = RIGID_BODY_MASS(i)%RMID
          DO n = i+1,NUMRM
            IF (ID .EQ. RIGID_BODY_MASS(n)%RMID) THEN
              CALL NON_UNIQUE_ID ('RIGID_BODY_MASS',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      RM_INTERNAL_ID = 0
!!
      DO i = 1,NUMRM
        IF (ExtID .EQ. RIGID_BODY_MASS(i)%RMID) THEN
           RM_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION MP_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:28:51
!!
      USE shared_common_data
      USE massprop_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMMP-1
          ID = MASSPROP(i)%MPID
          DO n = i+1,NUMMP
            IF (ID .EQ. MASSPROP(n)%MPID) THEN
              CALL NON_UNIQUE_ID ('MASSPROP',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      MP_INTERNAL_ID = 0
!!
      DO i = 1,NUMMP
        IF (ExtID .EQ. MASSPROP(i)%MPID) THEN
           MP_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION NS_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:28:58
!!
      USE shared_common_data
      USE node_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMNS-1
          ID = NODE_SET(i)%SetID
          DO n = i+1,NUMNS
            IF (ID .EQ. NODE_SET(n)%SetID) THEN
              CALL NON_UNIQUE_ID ('NPSET',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      NS_INTERNAL_ID = 0
!!
      DO i = 1,NUMNS
        IF (ExtID .EQ. NODE_SET(i)%SetID) THEN
           NS_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_ELEMENT_SET_ID (RECORD_KW,RECORD_ID,ELEMSET_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1997 22:24:37
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record key word ID
     &          ELEMSET_ID,     & ! I/O Element set ID, relinked
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch    
     &          ES_INTERNAL_ID    ! I/- External integer function
!!
!! Look for internal ID. (The internal ID is the index n of the element
!! set record in the sequence ELEMENT_SET(1:NUMES) which has the user
!! specified ID, ELEMSET_ID.)
!!
      INTERNAL_ID = ES_INTERNAL_ID (ELEMSET_ID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        ELEMSET_ID = INTERNAL_ID
      ELSE
!!
!! Input element set ID not found among the element set records
!! ELEMENT_SET(1:NUMES). OK, only called by PRINT & RESULTS.
!!
!!        L = LEN (RECORD_KW)
!!        WRITE (MSG1,'(I8)') RECORD_ID
!!        WRITE (MSG2,'(I8)') ELEMSET_ID
!!        CALL USER_MESSAGE
!!      2       (
!!      2       MSGL//'WARNING'//
!!      2       MSGL//'RELINK_ELEMENT_SET_ID.001.00'//
!!      2       MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//
!!      2       MSGL//'Contains Unknown ELSET ID:'//MSG2
!!      2       )
!!        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION ES_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:08
!!
      USE shared_common_data
      USE element_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMES-1
          ID = ELEMENT_SET(i)%SetID
          DO n = i+1,NUMES
            IF (ID .EQ. ELEMENT_SET(n)%SetID) THEN
              CALL NON_UNIQUE_ID ('ELSET',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      ES_INTERNAL_ID = 0
!!
      DO i = 1,NUMES
        IF (ExtID .EQ. ELEMENT_SET(i)%SetID) THEN
           ES_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION SS_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:14
!!
      USE shared_common_data
      USE segment_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMSS-1
          ID = SEGMENT_SET(i)%SetID
          DO n = i+1,NUMSS
            IF (ID .EQ. SEGMENT_SET(n)%SetID) THEN
              CALL NON_UNIQUE_ID ('SEGSET',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      SS_INTERNAL_ID = 0
!!
      DO i = 1,NUMSS
        IF (ExtID .EQ. SEGMENT_SET(i)%SetID) THEN
           SS_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_MatID (RECORD_KW,RECORD_ID,MATERIAL_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1990 22:24:37
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record key word ID
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch    
     &          MATERIAL_ID,    & ! I/O Material ID, relinked   
     &          MT_INTERNAL_ID    ! I/- External integer function
!!
!! Look for internal ID. (The internal ID is the index n of the material
!! record in the sequence MATERIAL(1:NUMMT) which has the user specified
!! ID, MATERIAL_ID.)
!!
      INTERNAL_ID = MT_INTERNAL_ID (MATERIAL_ID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        MATERIAL_ID = INTERNAL_ID
      ELSE
!!
!! Input material ID not found among the material records MATERIAL(1:NUMMT)
!!
        L = LEN (RECORD_KW)
        WRITE (MSG1,'(I8)') RECORD_ID
        WRITE (MSG2,'(I8)') MATERIAL_ID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_MatID.001.00'//                                  &
     &          MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//              &
     &          MSGL//'Contains Unknown Material ID:'//MSG2                    &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION MT_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:23
!!
      USE shared_common_data
      USE material_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMMT-1
          ID = MATERIAL(i)%MatID
          DO n = i+1,NUMMT
            IF (ID .EQ. MATERIAL(n)%MatID) THEN
              CALL NON_UNIQUE_ID ('MATERIAL',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      MT_INTERNAL_ID = 0
!!
      DO i = 1,NUMMT
        IF (ExtID .EQ. MATERIAL(i)%MatID) THEN
           MT_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION TF_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:31
!!
      USE shared_common_data
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMTF-1
          ID = TABULATED_FUNCTION(i)%TFID
          DO n = i+1,NUMTF
            IF (ID .EQ. TABULATED_FUNCTION(n)%TFID) THEN
              CALL NON_UNIQUE_ID ('TABFTN',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      TF_INTERNAL_ID = 0
!!
      DO i = 1,NUMTF
        IF (ExtID .EQ. TABULATED_FUNCTION(i)%TFID) THEN
           TF_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE RELINK_LupID (RECORD_KW,RECORD_ID,LAYUP_ID)
!!
!! Copyright (c) by KEY Associates, 6-SEP-1990 22:24:37
!!
      USE shared_common_data
      USE layering_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                 &
     &          RECORD_KW*(*)     ! I/- Input record key word
      INTEGER                   &
     &          RECORD_ID,      & ! I/- Input record ID      
     &          INTERNAL_ID,    & ! -/- Internal ID, scratch 
     &          LAYUP_ID,       & ! I/O Layering ID, relinked
     &          LY_INTERNAL_ID    ! I/- External integer function
!!
!! Look for internal ID. (The internal ID is the index n of the layering
!! record in the sequence LAYERING(1:NUMLU) which has the user specified
!! ID, LAYUP_ID.)
!!
      INTERNAL_ID = LY_INTERNAL_ID (LAYUP_ID)
      IF (INTERNAL_ID .NE. 0) THEN
!!
!! Relink. (Change external ID into internal ID.)
!!
        LAYUP_ID = INTERNAL_ID
      ELSE
!!
!! Input Lay-up ID not found among the Lay-up records LAYERING(1:NUMLU)%
!!
        L = LEN (RECORD_KW)
        WRITE (MSG1,'(I8)') RECORD_ID
        WRITE (MSG2,'(I8)') LAYUP_ID
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RELINK_LupID.001.00'//                                  &
     &          MSGL//RECORD_KW(1:L)//' Input Record ID:'//MSG1//              &
     &          MSGL//'Contains Unknown LAYUP ID:'//MSG2                       &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ENDIF
!!
      RETURN
      END
!!_
      INTEGER FUNCTION LY_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:37
!!
      USE shared_common_data
      USE layering_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMLU-1
          ID = LAYERING(i)%LupID
          DO n = i+1,NUMLU
            IF (ID .EQ. LAYERING(n)%LupID) THEN
              CALL NON_UNIQUE_ID ('LAYUP',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      LY_INTERNAL_ID = 0
!!
      DO i = 1,NUMLU
        IF (ExtID .EQ. LAYERING(i)%LupID) THEN
           LY_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION S1_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:43
!!
      USE shared_common_data
      USE section_1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMS1-1
          ID = SECTION_1D(i)%SecID
          DO n = i+1,NUMS1
            IF (ID .EQ. SECTION_1D(n)%SecID) THEN
              CALL NON_UNIQUE_ID ('BSECTION',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      S1_INTERNAL_ID = 0
!!
      DO i = 1,NUMS1
        IF (ExtID .EQ. SECTION_1D(i)%SecID) THEN
           S1_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION S2_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:49
!!
      USE shared_common_data
      USE section_2d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMS2-1
          ID = SECTION_2D(i)%SecID
          DO n = i+1,NUMS2
            IF (ID .EQ. SECTION_2D(n)%SecID) THEN
              CALL NON_UNIQUE_ID ('PSECTION',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      S2_INTERNAL_ID = 0
!!
      DO i = 1,NUMS2
        IF (ExtID .EQ. SECTION_2D(i)%SecID) THEN
           S2_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION SI_INTERNAL_ID (ExtID)
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:29:58
!!
      USE shared_common_data
      USE sliding_interface_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: ExtID ! User assigned ID for which internal ID is sought.
!!
!! Local variables.
      INTEGER       :: ID          ! -/- External ID value
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! On first entry check for non-unique identifier.
!!
      IF (FIRST) THEN
        DO i = 1,NUMSI-1
          ID = SLIDING_INTERFACE(i)%SIID
          DO n = i+1,NUMSI
            IF (ID .EQ. SLIDING_INTERFACE(n)%SIID) THEN
              CALL NON_UNIQUE_ID ('SLIDE',ID)
            ENDIF
          ENDDO
        ENDDO
        FIRST = .FALSE.
      ENDIF
!!
!! Set return value to zero; flag for "not found."
!!
      SI_INTERNAL_ID = 0
!!
      DO i = 1,NUMSI
        IF (ExtID .EQ. SLIDING_INTERFACE(i)%SIID) THEN
           SI_INTERNAL_ID = i
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE NON_UNIQUE_ID (KEY_WORD,ID)
!!
!! Copyright (c) by KEY Associates, 4-APR-1992 08:56:44
!!
!! Purpose: Single call to initiate a user message for non-unique identifier.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          KEY_WORD*(*)    ! I/- Input record key word
      INTEGER                                                                  &
     &          ID              ! I/- Repeated identifier
!!
      L = LEN (KEY_WORD)
      WRITE (MSG1,'(I8)') ID
      CALL USER_MESSAGE                                                        &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'NON_UNIQUE_ID.001.00'//                                 &
     &          MSGL//'Repeated '//KEY_WORD(1:L)//' ID''s Found:'//MSG1        &
     &          )
      ERROR%COUNT = ERROR%COUNT + 1
!!
      RETURN
      END
!!_
      SUBROUTINE RESTART_RELINK_USER_IDS
!!
!! Copyright (c) by KEY Associates;  2-MAY-1992 08:25:26.82
!!
!! Purpose: Relink restart input using internal numbers in place of User ID's.
!!
      USE shared_common_data
      USE results_
      USE qa_record_
      USE massprop_
      USE gauge1d_
      USE gauge2d_
      USE gauge3d_
      USE material_
      USE layering_
      USE section_2d_
      USE section_1d_
      USE rigid_body_
      USE rigid_body_mass_
      USE nodal_point_mass_
      USE displacement_bc_
      USE tied_bc_
      USE spot_weld_
      USE rigid_wall_bc_
      USE body_force_
      USE pressure_bc_
      USE force_bc_
      USE spring_bc_
      USE damper_bc_
      USE periodic_bc_
      USE nonreflecting_bc_
      USE sliding_interface_
      USE tabulated_function_
      USE node_set_
      USE element_set_
      USE segment_set_
      USE segment_
      USE node_
      USE motion_
      USE force_
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
      USE constrained_node_
      USE plate_pair_
      USE enumerated_sets_
      USE output_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          NODAL_ID,                                                      &
     &          INTERNAL_ID,                                                   &
     &          NS_INTERNAL_ID,                                                &
     &          ES_INTERNAL_ID,                                                &
     &          SI_INTERNAL_ID
      EXTERNAL                                                                 &
     &          NS_INTERNAL_ID,                                                &
     &          ES_INTERNAL_ID,                                                &
     &          SI_INTERNAL_ID
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered RESTART_RELINK_USER_IDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Relink pointers in INTERFACE_TIME.
!!
      DO n = 1,NUMIT
        INTERNAL_ID = SI_INTERNAL_ID (INTERFACE_TIME%SI(n)%ID)
        IF (INTERNAL_ID .NE. 0) THEN
          INTERFACE_TIME%SI(n)%ID = INTERNAL_ID
        ELSE
!!
!! Matching sliding interface ID not found.
!!
          WRITE (MSG2,'(I8)') INTERFACE_TIME%SI(n)%ID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RESTART_RELINK_USER_IDS.001.01'//                       &
     &          MSGL//'An INTERFACE Input Record References'//                 &
     &          MSGL//'Unknown SLIDING_INTERFACE ID:'//MSG2                    &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
!! Relink pointers in PRINT.
!!
      IF (PRINT%NODES .GT. 0) THEN
        INTERNAL_ID = NS_INTERNAL_ID (PRINT%NODES)
        IF (INTERNAL_ID .NE. 0) THEN
          PRINT%NODES = INTERNAL_ID
        ELSE
!!
!! Matching node set ID not found.
!!
          WRITE (MSG1,'(I8)') PRINT%NODES
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RESTART_RELINK_USER_IDS.001.02'//                       &
     &          MSGL//'PRINT Input Record Entry NPT'//                         &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG1                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDIF

      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword HXEL,',  0,PRINT%HEXAH )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword PXEL,',  0,PRINT%PENTA )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword TXEL,',  0,PRINT%TETRA )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword M3EL,',  0,PRINT%MEMBT )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword M4EL,',  0,PRINT%MEMBQ )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword TRUSS,', 0,PRINT%TRUSS )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword P3EL,',  0,PRINT%PLATT )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword P4EL,',  0,PRINT%PLATQ )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword BEAM,',  0,PRINT%BEAMS )
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword SPRING,',0,PRINT%SPRING)
      CALL RELINK_ELEMENT_SET_ID                                               &
     &  ('PRINT Keyword DAMPER,',0,PRINT%DAMPER)
!!
!! Relink pointers in RESULTS.
!!
      DO i = 1,NUMRF
        IF (RESULTS(i)%NODES .GT. 0) THEN
          INTERNAL_ID = NS_INTERNAL_ID (RESULTS(i)%NODES)
          IF (INTERNAL_ID .NE. 0) THEN
            RESULTS(i)%NODES = INTERNAL_ID
          ELSE
!!
!! Matching node set ID not found.
!!
            WRITE (MSG1,'(I8)') RESULTS(i)%ResID
            WRITE (MSG2,'(I8)') RESULTS(i)%NODES
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RESTART_RELINK_USER_IDS.001.04'//                       &
     &          MSGL//'RESULTS Input Record ID:'//MSG1//                       &
     &          MSGL//'Input Record Parameter: '//OUTPUT%NAME(7)//             &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDIF
      ENDDO
!!
      DO i = 1,NUMRF
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword HXEL,', RESULTS(i)%ResID,RESULTS(i)%HEXAH )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword PXEL,', RESULTS(i)%ResID,RESULTS(i)%PENTA )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword TXEL,', RESULTS(i)%ResID,RESULTS(i)%TETRA )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword LSEL,', RESULTS(i)%ResID,RESULTS(i)%LSOLD )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword M3EL,', RESULTS(i)%ResID,RESULTS(i)%MEMBT )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword M4EL,', RESULTS(i)%ResID,RESULTS(i)%MEMBQ )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword TRUSS,',RESULTS(i)%ResID,RESULTS(i)%TRUSS )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword P3EL,', RESULTS(i)%ResID,RESULTS(i)%PLATT )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword P4EL,', RESULTS(i)%ResID,RESULTS(i)%PLATQ )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword BEAM,', RESULTS(i)%ResID,RESULTS(i)%BEAMS )
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword SPRING,',RESULTS(i)%ResID,RESULTS(i)%SPRING)
        CALL RELINK_ELEMENT_SET_ID                                             &
     &    ('RESULTS Keyword DAMPER,',RESULTS(i)%ResID,RESULTS(i)%DAMPER)
      ENDDO
!!
!! Relink pointers in MASSPROP.
!!
      DO i = 1,NUMMP
        INTERNAL_ID = NS_INTERNAL_ID (MASSPROP(i)%SetID)
        IF (INTERNAL_ID .NE. 0) THEN
          MASSPROP(i)%SetID = INTERNAL_ID
        ELSE
!!
!! Matching node set ID not found.
!!
          WRITE (MSG1,'(I8)') MASSPROP(i)%MPID
          WRITE (MSG2,'(I8)') MASSPROP(i)%SetID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARNING'//                                              &
     &          MSGL//'RESTART_RELINK_USER_IDS.001.06'//                       &
     &          MSGL//'MASSPROP Input Record ID:'//MSG1//                      &
     &          MSGL//'Contains Unknown NPSET ID:'//MSG2                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ENDDO
!!
      RETURN
      END
