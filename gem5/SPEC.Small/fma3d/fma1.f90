
      SUBROUTINE FMA_3D (MAXIF)
!!
!!                            FMA-3D  Version 10.01
!!
!!            Copyright (c) 1994 by Samuel W. Key dba KEY Associates
!!            Copyright (c) 1995 by Samuel W. Key dba KEY Associates
!!            Copyright (c) 1997 by Samuel W. Key dba KEY Associates
!!                        1851 Tramway Terrace Loop NE
!!                       Albuquerque, New Mexico  87122
!!                     Telephone and FAX:  (505) 856-1488
!!                        August  1994 (Version  8.01)
!!                        October 1995 (Version  9.01)
!!                        October 1997 (Version 10.01)
!!
!! All rights reserved. No part of this program may be reproduced, stored in
!! a retrieval system, or transmitted, in any form or by any means, elec-
!! tronic, mechanical, photocopying, recording or otherwise,  without the
!! prior written permission of KEY Associates. Printed in the United States
!! of America.
!!
!! FMA-3D is designed to compute the time-dependent displacements and stresses
!! within elastic or inelastic, three dimensional bodies of arbitrary shapes
!! and materials.  The Galerkin form of the finite element method is used to
!! generate the spatial discretization. A mesh composed of 2-node truss,
!! 2-node beam, 4-node and 3-node membrane, 4-node and 3-node plate, and
!! 8-node, 6-node and 4-node solid elements is used.  The displacements are
!! assumed to vary linearly, bilinearly, and trilinearly over the respective
!! elements using isoparametric coordinates. A mean stress quadrature and
!! orthogonal hourglass control due to Dennis Flanagan are used. The resulting
!! simultaneous equations are integrated in time using a central difference
!! method integrator. Courant subcycling is provided.
!!
!!* * * * * * * * * * * * * * *  DISCLAIMER * * * * * * * * * * * * * * * * *
!!* KEY Associates makes no representation, express or implied, with        *
!!* respect to this software or the documentation that describes it,        *
!!* including without limitation, any implied warranties of merchantabil-   *
!!* ity or fitness for a particular purpose, all of which are expressly     *
!!* disclaimed.  KEY Associates, its distributors and dealers shall in no   *
!!* event be liable for any indirect, incidental or consequential damages.  *
!!* The exclusion of implied warranties is not permitted by some statutes.  *
!!* The above exclusion therefore may not apply to you. This warranty       *
!!* provides you with specific legal rights. There may be other rights      *
!!* that you have which may vary from state to state.                       *
!!*                                                                         *
!!* Key Associates neither makes any warranty, express or implied, nor      *
!!* assumes any legal liability or responsibility for the accuracy, com-    *
!!* pleteness or usefulness of any information, apparatus, product, or      *
!!* process disclosed, or represents that its use would not infringe        *
!!* privately owned rights.                                                 *
!!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
      USE shared_common_data;
      USE indx_;           USE node_;           USE tabulated_function_;
      USE beam_;           USE coord_;          USE sliding_interface_;
      USE value_;          USE force_;          USE constrained_node_;
      USE hexah_;          USE penta_;          USE nonreflecting_bc_;
      USE tetra_;          USE lsold_;          USE nodal_point_mass_;
      USE membq_;          USE membt_;          USE rigid_body_mass_;
      USE truss_;          USE platq_;          USE state_variables_;
      USE platt_;          USE motion_;         USE enumerated_sets_;
      USE spring_;         USE damper_;         USE displacement_bc_;
      USE stress_;         USE segment_;        USE contact_surface_;
      USE tied_bc_;        USE results_;        USE relink_scratch_;
      USE gauge1d_;        USE gauge2d_;        USE rigid_wall_bc_;
      USE gauge3d_;        USE massprop_;       USE include_file_;
      USE material_;       USE layering_;       USE sliding_node_;
      USE force_bc_;       USE node_set_;       USE contact_node_;
      USE nrbc_data_;      USE spring_bc_;      USE periodic_bc_;
      USE damper_bc_;      USE spot_weld_;      USE pressure_bc_;
      USE qa_record_;      USE plate_pair_;     USE segment_set_;
      USE body_force_;     USE section_2d_;     USE element_set_;
      USE section_1d_;     USE rigid_body_;     USE plate_pair_;
      USE velocity_ic_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER :: MAXIF, IALLOC_FLAG
!!
!! Initilize derived data types in the USE file: shared_common_data
!!
      CALL INITIALIZE_SHARED_COMMON_DATA
!!
!! Identify the program and this version of the program.
!!
      CALL PROGRAM_ID                                                          &
     &          (                                                              &
     &          JOB_ID_RECORD%PROGRAM,                                         &
     &          JOB_ID_RECORD%VERSION,                                         &
     &          JOB_ID_RECORD%RELEASE                                          &
     &          )
!!
!! Initialize CPU timer; record current date and time.
!!
!SPEC_CPU2000      CALL TIMER (0)
!SPEC_CPU2000      CALL ROLEX
!SPEC_CPU2000     &          (
!SPEC_CPU2000     &          JOB_ID_RECORD%CURRENT%DATE,
!SPEC_CPU2000     &          JOB_ID_RECORD%CURRENT%TIME
!SPEC_CPU2000     &          )
!!
!! Open execution log output file.
!!
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LELO,                                        &
     &          FILE   = 'fmaelo',                                             &
     &          STATUS = 'UNKNOWN',                                            &
     &          FORM   = 'FORMATTED'                                           &
     &          )
!!
      CALL INTERCALATION (IO_UNIT%LELO)
!!
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!! PART I. SCAN INPUT FILES AND ALLOCATE THE REQUISITE STORAGE.
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
!! Allocate storage needed for MAXIF included files.
!!
      ALLOCATE ( INCLUDE_FILE(1:MAXIF), STAT=IALLOC_FLAG )
      CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'INCLUDE_FILE',MAXIF)
!$OMP PARALLEL DO
      DO i = 1,MAXIF
        INCLUDE_FILE(i) = include_file_type (0,0,0,0,0,0,0,0,' ')
      ENDDO
!!
!! Allocate storage needed for free-field reader to count input items.
!!
      ALLOCATE ( VALUE(1:2), STAT=IALLOC_FLAG )
      CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'VALUE', 2)
      DO i = 1,2
        VALUE(i) = value_type (0,' ',' ',0,0,0)
      ENDDO
!!
!! Initialize INCLUDE_FILE(1) with user's root input file.
!!
      CALL INIT_INCLUDE_FILE
!!
!! Scan included files to get data-record entry counts and array sizes.
!!
      CALL SCAN_INCLUDED_FILES (MAXIF)
!!
!! Allocate storage needed for free-field reader to read ALL input items.
!!
      DEALLOCATE ( VALUE )
      ALLOCATE ( VALUE(1:MAXRE), STAT=IALLOC_FLAG )
      CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'VALUE', MAXRE)
!$OMP PARALLEL DO
      DO i = 1,MAXRE
        VALUE(i) = value_type (0,' ',' ',0,0,0)
      ENDDO
!!
!! Pre-read layered solids and lay-up descriptions to count the number of
!! 8-node solid LSHEX and 4-node membrane LSMBQ sub-elements to reserve.
!!
      IF (NUMLS .GT. 0) THEN

        ALLOCATE ( LAYERING(1:NUMLU), STAT=IALLOC_FLAG )
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LAYERING',NUMLU)
!$OMP PARALLEL DO
        DO i = 1,NUMLU
          LAYERING(i)%LupID = 0
          LAYERING(i)%Number_of_Layers = MXNLY
          LAYERING(i)%Isys = 0
          DO j = 1,MXNLY
            LAYERING(i)%LayID(j) = 0
            LAYERING(i)%Ltype(j) = 0
            LAYERING(i)%MatID(j) = 0
            LAYERING(i)%H(1:4,j) = (/0,0,0,0/)
            LAYERING(i)%Ax(j)    = 0
            LAYERING(i)%Ay(j)    = 0
            LAYERING(i)%Az(j)    = 0
            LAYERING(i)%Bx(j)    = 0
            LAYERING(i)%By(j)    = 0
            LAYERING(i)%Bz(j)    = 0
          ENDDO
        ENDDO

        ALLOCATE ( LSOLD(1:NUMLS), STAT=IALLOC_FLAG )
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LSOLD',NUMLS)
!$OMP PARALLEL DO
        DO i = 1,NUMLS
          LSOLD(i)%PAR = PAR_lsold(0,0,0,(/(0,j=1,8)/),                        &
     &      (/(0,j=1,MXLSL)/),0,0)
          LSOLD(i)%RES = RES_lsold(0,0,RESHAPE((/(0,j=1,4*MXLSL)/),            &
     &      (/4,MXLSL/)),(/(0,k=1,8)/),(/(0,k=1,8)/),(/(0,k=1,8)/),            &
     &      0,0,0,0)
        ENDDO

        CALL PRE_READ_LAYUP_RECORDS
        CALL PRE_READ_LSOLD_RECORDS
!!
!! Include the layered solid sub-elements in total number of elements.
!!
        NUMEL = NUMEL + NUMLX + NUMLM
!!
      ENDIF
!!
!! Pre-read execution control records to see if there is a switch set to read
!! a mesh data input file and to see if any other control switches are set.
!!
      CALL PRE_READ_CONTROL_RECORDS
!!
!! If this is a restart from a rezoned problem, CONTROL%RZSTAR should be
!! non-zero. A restart from a rezoned problem differs from an unchanged
!! restart only on the "resume" address used in SOLVE. In all other
!! respects the two restarts should be the same. Thus, if CONTROL%RDSTAR
!! is zero, CONTROL%RZSTAR is used to set CONTROL%RDSTAR to 1
!!
      IF (CONTROL%RZSTAR .NE. 0) THEN
        CONTROL%RZSTAR = 1
        CONTROL%RDSTAR = 1
      ENDIF
!!
!! Scan mesh data file for counters (read parts 1 and 2 of mesh data file)
!! to increment storage requirements.
!!
      IF (CONTROL%RDMESH.GT.0 .AND. CONTROL%RDSTAR.NE.1) THEN
!!
        CALL SCAN_MESH_DATA
!!
      ELSE IF (CONTROL%RDSTAR .EQ. 1) THEN
!!
!! Scan restart data file for counters for an "unchanged" restart. (QA records
!! are additive, old RESULTS records are replaced with new ones, and old inter-
!! face times are replaced with new ones.)
!!
        CALL SCAN_RESTART_DATA
!!
      ENDIF
!!
!! Allot storage for derived data types. Initialize defaulted values,
!! set initial conditions on velocity to zero. Clear displacements,
!! accelerations, and forces.
!!
      CALL ALLOCATE_STORAGE
!!
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!! PART II. READ AND PROCESS INPUT DATA.
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
!! Read and process included files (ASCII input records).
!!
      CALL READ_INCLUDED_FILES
!!
!! If this is an "unchanged" restart, the reading and building of set member
!! lists, and the construction of derived data is skipped. In its place the
!! data for all of the arrays will be taken directly from the restart file.
!!
      IF (CONTROL%RDSTAR .NE. 1) THEN
!!
!! Read mesh data file to complete input data (read part 3 of mesh data input).
!!
        IF (CONTROL%RDMESH .GT. 0) THEN
          CALL READ_MESH_DATA
        ENDIF
!!
!! READ NODE SETS AND BUILT CONCATENATED LIST OF NODE SET MEMBERS.
!! Step 1: Count locations needed to store concatenated list of nodal
!! point members for all node sets.
!!
        CALL READ_AND_BUILD_NODE_SETS ('COUNT')
!!
!! Step 2: Allocate storage for concatenated list of node set members.
!!
        IF (NUMNE .GT. 0) THEN
          ALLOCATE ( NNPSETS(1:NUMNE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NNPSETS',NUMNE)
!$OMP PARALLEL DO
          DO i = 1,NUMNE
            NNPSETS(i) = 0
          ENDDO
!!
!! Step 3: Populate concatenated list of nodal point members for all node sets.
!!
          CALL READ_AND_BUILD_NODE_SETS ('BUILD')
        ENDIF
!!
!! READ ELEMENT SETS AND BUILT CONCATENATED LIST OF ELEMENT SET MEMBERS.
!! Step 1: Count locations needed to store concatenated list of element
!! members for all element sets.
!!
        CALL READ_AND_BUILD_ELEMENT_SETS ('COUNT')
!!
!! Step 2: Allocate storage for concatenated list of element set members
!!
        IF (NUMEE .GT. 0) THEN
          ALLOCATE ( NELSETS(1:NUMEE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NELSETS',NUMEE)
!$OMP PARALLEL DO
          DO i = 1,NUMEE
            NELSETS(i) = 0
          ENDDO
!!
!! Step 3: Populate concatenated list of element members for all element sets.
!!
          CALL READ_AND_BUILD_ELEMENT_SETS ('BUILD')
        ENDIF
!!
!! READ SEGMENT SETS AND BUILT CONCATENATED LIST OF SEGMENT SET MEMBERS.
!! Step 1: Count locations needed to store concatenated list of segment
!! members for all segment sets.
!!
        CALL READ_AND_BUILD_SEGMENT_SETS ('COUNT')
!!
!! Step 2: Allocate storage for concatenated list of segment set members.
!!
        IF (NUMSE .GT. 0) THEN
          ALLOCATE ( NSGSETS(1:NUMSE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NSGSETS',NUMSE)
!$OMP PARALLEL DO
          DO i = 1,NUMSE
            NSGSETS(i) = 0
          ENDDO
!!
!! Step 3: Populate concatenated list of segment members for all segment sets.
!!
          CALL READ_AND_BUILD_SEGMENT_SETS ('BUILD')
        ENDIF
!!
!! Close open mesh data file.
!!
        IF (CONTROL%RDMESH .GT. 0) THEN
          CLOSE (UNIT=IO_UNIT%LMDI,STATUS='KEEP')
        ENDIF
!!
!! Write mesh data.
!!
        IF (CONTROL%WRMESH .GT. 0) THEN
          CALL WRITE_MESH_DATA
          RETURN
        ENDIF
!!
!! Connect up layered solids to auxillary hexahedrons and membranes.
!!
        IF (NUMLS .GT. 0) THEN
          CALL CONNECT_LAYERED_SOLIDS
        ENDIF
!!
!! Allocate temporary storage for relinking.
!!
        IF (NUMNP .GT. 0) THEN
          ALLOCATE ( INPV(1:NUMNP), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'INPV',NUMNP)
          ALLOCATE ( JNPV(1:NUMNP), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'JNPV',NUMNP)
!$OMP PARALLEL DO
          DO i = 1,NUMNP
            INPV(i) = 0
            JNPV(i) = 0
          ENDDO
        ENDIF

        IF (NUMEL .GT. 0) THEN
          ALLOCATE ( IELV(1:NUMEL), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'IELV',NUMEL)
          ALLOCATE ( JELV(1:NUMEL), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'JELV',NUMEL)
!$OMP PARALLEL DO
          DO i = 1,NUMEL
            IELV(i) = 0
            JELV(i) = 0
          ENDDO
        ENDIF

        IF (NUMSG .GT. 0) THEN
          ALLOCATE ( ISGV(1:NUMSG), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'ISGV',NUMSG)
          ALLOCATE ( JSGV(1:NUMSG), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'JSGV',NUMSG)
!$OMP PARALLEL DO
          DO i = 1,NUMSG
            ISGV(i) = 0
            JSGV(i) = 0
          ENDDO
        ENDIF
!!
!! Relink using internal numbers in place of user assigned ID's.
!!
        CALL RELINK_USER_IDS
!!
!! Release temporary storage.
!!
        IF (NUMNP .GT. 0) DEALLOCATE ( INPV, JNPV )
        IF (NUMEL .GT. 0) DEALLOCATE ( IELV, JELV )
        IF (NUMSG .GT. 0) DEALLOCATE ( ISGV, JSGV )
        IF (MAXRE .GT. 0) DEALLOCATE ( VALUE )
!!
!! Compute the amount of storage needed for shell stresses and material
!! state variables. Each stress record has six components: STRESS(1:6,*)=
!! (Sigma_xx,Sigma_yy,Sigma_zz,Sigma_xy,Sigma_xz,Sigma_yz). State varia-
!! bles are stored sequentially in STATE_VARIABLES due to the variable
!! number of items needed by each material model. Compute the amount of
!! storage needed for the data arrays used to implement the sliding
!! interface calculations. Compute the amount of storage needed for
!! nonreflecting boundary conditions.
!!
        CALL GET_AUXILIARY_STORAGE_LENGTH
!!
!! Allocate storage for the shell stresses STRESS(1:6,1:NUMST), and for
!! the material state variable storage array STATE_VARIABLES(1:NUMAX).
!! Allocate storage for the the sliding interface data arrays SLIDING_NODE,
!! CONTACT_SURFACE, and CONTACT_NODE. Allocate storage for nonreflecting
!! boundary condition data NRBC_DATA(1:NUMND).
!!
        IF (NUMST .GT. 0) THEN
          ALLOCATE ( STRESS(1:6,1:NUMST), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'STRESS',NUMST)
!$OMP PARALLEL DO
          DO i = 1,NUMST
            STRESS(1:6,i) = (/0,0,0,0,0,0/)
          ENDDO
        ENDIF

        IF (NUMAX .GE. 0) THEN
          MXAX = MAX (1,NUMAX)
          ALLOCATE ( STATE_VARIABLES(1:MXAX), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'STATE_VARIABLES',NUMAX)
!$OMP PARALLEL DO
          DO i = 1,MXAX
            STATE_VARIABLES(i) = 0
          ENDDO
        ENDIF

        IF (NUMND .GT. 0) THEN
          ALLOCATE ( NRBC_DATA(1:NUMND), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NRBC_DATA',NUMND)
!$OMP PARALLEL DO
          DO i = 1,NUMND
            NRBC_DATA(i) = nrbc_data_type (0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMSN .GT. 0) THEN
          ALLOCATE ( SLIDING_NODE(1:NUMSN), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'SLIDING_NODE',NUMSN)
!$OMP PARALLEL DO
          DO i = 1,NUMSN
            SLIDING_NODE(i) = sliding_node_type (0,0,(/0,0,0,0/),              &
     &          0,0,0,0,0,0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCE .GT. 0) THEN
          ALLOCATE ( CONTACT_SURFACE(1:NUMCE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'CONTACT_ELEMENT',NUMCE)
!$OMP PARALLEL DO
          DO i = 1,NUMCE
            CONTACT_SURFACE(i) = contact_surface_type ((/0,0,0,0/),            &
     &          (/0,0,0,0/),(/0,0,0/),0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCN .GT. 0) THEN
          ALLOCATE ( CONTACT_NODE(1:NUMCN), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'CONTACT_NODE',NUMCN)
!$OMP PARALLEL DO
          DO i = 1,NUMCN
            CONTACT_NODE(i) = contact_node_type (0,0,0,0,0,0,0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCX .GT. 0) THEN
          ALLOCATE ( COORD(1:NUMCX,1:3), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'COORD(1:NUMCX,1:3)',NUMCX)
!$OMP PARALLEL DO
          DO i = 1,NUMCX
            COORD(i,1) = 0
            COORD(i,2) = 0
            COORD(i,3) = 0
          ENDDO
        ENDIF

        IF (NUMCX .GT. 0) THEN
          ALLOCATE ( INDX(1:NUMCX,1:8), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'INDX(1:NUMCX,1:8)',NUMCX)
!$OMP PARALLEL DO
          DO i = 1,NUMCX
            INDX(i,1) = 0
            INDX(i,2) = 0
            INDX(i,3) = 0
            INDX(i,4) = 0
            INDX(i,5) = 0
            INDX(i,6) = 0
            INDX(i,7) = 0
            INDX(i,8) = 0
          ENDDO
        ENDIF
!!
!! Build plate-pair records and allocate storage for plate pairs.
!!
        IF (CONTROL%REZONE .NE. 0) THEN

          CALL BUILD_PLATE_PAIRS ('COUNT')

          IF (NUMPP .GT. 0) THEN
            ALLOCATE ( PLATE_PAIR(1:NUMPP), STAT=IALLOC_FLAG )
            CALL REPORT_ALLOCATION_EVENT                                       &
     &        (IALLOC_FLAG,'PLATE_PAIR',NUMPP)
!$OMP PARALLEL DO
            DO i = 1,NUMPP
              PLATE_PAIR(i) = plate_pair_type (0,IBIDIS(0,0),                  &
     &          IBIDIS(0,0),0,0,0)
            ENDDO
            CALL BUILD_PLATE_PAIRS ('BUILD')
          ENDIF

        ENDIF
!!
!! Construct working arrays from data contained in the input arrays.
!!
        CALL INITIALIZE_WORKING_STORAGE
!!
!! Release temporary storage.
!!
        IF (NUMIC .GT. 0) DEALLOCATE ( VELOCITY_IC )
!!
!! As an alternative to generating the preceding data an "unchanged" restart
!! will take the data from the restart file.
!!
      ELSE
!!
!! Allocate the balance of storage needed to resume calculation.
!!
        IF (NUMNE .GT. 0) THEN
          ALLOCATE ( NNPSETS(1:NUMNE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NNPSETS',NUMNE)
!$OMP PARALLEL DO
          DO i = 1,NUMNE
            NNPSETS(i) = 0
          ENDDO
        ENDIF

        IF (NUMEE .GT. 0) THEN
          ALLOCATE ( NELSETS(1:NUMEE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NELSETS',NUMEE)
!$OMP PARALLEL DO
          DO i = 1,NUMEE
            NELSETS(i) = 0
          ENDDO
        ENDIF

        IF (NUMSE .GT. 0) THEN
          ALLOCATE ( NSGSETS(1:NUMSE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NSGSETS',NUMSE)
!$OMP PARALLEL DO
          DO i = 1,NUMSE
            NSGSETS(i) = 0
          ENDDO
        ENDIF

        IF (NUMST .GT. 0) THEN
          ALLOCATE ( STRESS(1:6,1:NUMST), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'STRESS',NUMST)
!$OMP PARALLEL DO
          DO i = 1,NUMST
            STRESS(1:6,i) = (/0,0,0,0,0,0/)
          ENDDO
        ENDIF

        IF (NUMAX .GE. 0) THEN
          MXAX = MAX (1,NUMAX)
          ALLOCATE ( STATE_VARIABLES(1:MXAX), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'STATE_VARIABLES',NUMAX)
!$OMP PARALLEL DO
          DO i = 1,MXAX
            STATE_VARIABLES(i) = 0
          ENDDO
        ENDIF

        IF (NUMND .GT. 0) THEN
          ALLOCATE ( NRBC_DATA(1:NUMND), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NRBC_DATA',NUMND)
!$OMP PARALLEL DO
          DO i = 1,NUMND
            NRBC_DATA(i) = nrbc_data_type (0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMSN .GT. 0) THEN
          ALLOCATE ( SLIDING_NODE(1:NUMSN), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'SLIDING_NODE',NUMSN)
!$OMP PARALLEL DO
          DO i = 1,NUMSN
            SLIDING_NODE(i) = sliding_node_type (0,0,(/0,0,0,0/),              &
     &          0,0,0,0,0,0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCE .GT. 0) THEN
          ALLOCATE ( CONTACT_SURFACE(1:NUMCE), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'CONTACT_ELEMENT',NUMCE)
!$OMP PARALLEL DO
          DO i = 1,NUMCE
            CONTACT_SURFACE(i) = contact_surface_type ((/0,0,0,0/),            &
     &          (/0,0,0,0/),(/0,0,0/),0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCN .GT. 0) THEN
          ALLOCATE ( CONTACT_NODE(1:NUMCN), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'CONTACT_NODE',NUMCN)
!$OMP PARALLEL DO
          DO i = 1,NUMCN
            CONTACT_NODE(i) = contact_node_type (0,0,0,0,0,0,0,0,0,0,0)
          ENDDO
        ENDIF

        IF (NUMCX .GT. 0) THEN
          ALLOCATE ( COORD(1:NUMCX,1:3), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'COORD(1:NUMCX,1:3)',NUMCX)
!$OMP PARALLEL DO
          DO i = 1,NUMCX
            COORD(i,1) = 0
            COORD(i,2) = 0
            COORD(i,3) = 0
          ENDDO
        ENDIF

        IF (NUMCX .GT. 0) THEN
          ALLOCATE ( INDX(1:NUMCX,1:8), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'INDX(1:NUMCX,1:8)',NUMCX)
!$OMP PARALLEL DO
          DO i = 1,NUMCX
            INDX(i,1) = 0
            INDX(i,2) = 0
            INDX(i,3) = 0
            INDX(i,4) = 0
            INDX(i,5) = 0
            INDX(i,6) = 0
            INDX(i,7) = 0
            INDX(i,8) = 0
          ENDDO
        ENDIF

        IF (NUMPP .GT. 0) THEN
          ALLOCATE ( PLATE_PAIR(1:NUMPP), STAT=IALLOC_FLAG )
          CALL REPORT_ALLOCATION_EVENT                                         &
     &      (IALLOC_FLAG,'PLATE_PAIR',NUMPP)
!$OMP PARALLEL DO
          DO i = 1,NUMPP
            PLATE_PAIR(i) = plate_pair_type (0,IBIDIS(0,0),                    &
     &        IBIDIS(0,0),0,0,0)
          ENDDO
        ENDIF

        CALL READ_RESTART_DATA
!!
!! Relink restart input using internal numbers in place of User ID's.
!!
        CALL RESTART_RELINK_USER_IDS
!!
!! Perform consistency checks on restart data. (Needed primarily by
!! restart data generated by a rezone.)
!!
        CALL RESTART_CONSISTENCY_CHECK
!!
!! Impose any changes requested in sliding interface beginning and ending times.
!!
        CALL RESTART_INTERFACE_TIMES

      ENDIF
!!
!! When requested, provide a formatted print of the data defining the problem.
!!
      IF (CONTROL%INECHO .NE. 0) CALL FORMATTED_PRINT_OF_INPUT
!!
!! Integrate forward in time to obtain transient dynamic response.
!!
      CALL SOLVE
!!
!! Report timing results.
!!
!SPEC_CPU2000      CALL TIMER (100)
!!
      END
!!_
      SUBROUTINE INITIALIZE_SHARED_COMMON_DATA
!!
!! Copyright (c) by KEY Associates;  5-NOV-1997 20:07:06.00
!!
!! Purpose: Initialize shared common data contained in the USE module:
!!
!!          ./fma3d/v*/use/shared_common_data.f90
!!
!! These common blocks and structures are used throughout the program.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Set all NUMxx counters to zero.
!!
      NUMIF = 0;  NUMQA = 0;  NUMNP = 0;  NUMEL = 0;  NUMHX = 0;
      NUMPX = 0;  NUMTX = 0;  NUMLS = 0;  NUMLX = 0;  NUMLM = 0;
      NUMM4 = 0;  NUMM3 = 0;  NUMTR = 0;  NUMP4 = 0;  NUMP3 = 0;
      NUMBM = 0;  NUMSP = 0;  NUMDM = 0;  NUMSG = 0;  NUMDC = 0;
      NUMTC = 0;  NUMSW = 0;  NUMWC = 0;  NUMBF = 0;  NUMPC = 0;
      NUMFC = 0;  NUMSC = 0;  NUMVC = 0;  NUMCC = 0;  NUMNR = 0;
      NUMND = 0;  NUMFS = 0;  NUMIT = 0;  NUMSI = 0;  NUMSN = 0;
      NUMCE = 0;  NUMCN = 0;  NUMCX = 0;  NUMNS = 0;  NUMNE = 0;
      NUMES = 0;  NUMEE = 0;  NUMSS = 0;  NUMSE = 0;  NUMMT = 0;
      NUMLU = 0;  NUMTF = 0;  NUMFP = 0;  NUMRF = 0;  NUMPV = 0;
      NUMS1 = 0;  NUMS2 = 0;  NUMG1 = 0;  NUMG2 = 0;  NUMG3 = 0;
      NUMMP = 0;  NUMRM = 0;  NUMCM = 0;  NUMIC = 0;  NUMRB = 0;
      NUMRT = 0;  NUMST = 0;  NUMAX = 0;  NUMPP = 0;  NUMNC = 0;
      MAXRE = 0;  MXEID = 0;
!!
!! MAUX(0:60) is an array that contains the number of state variable locations
!! needed at a stress point for each material model (constitutive equation).
!!
!!      Entries          Stress Behavior           Elements
!!
!!    0 through  9   Force-deflection behavior   Spring and damper elements
!!   10 through 19   Uniaxial stress behavior    Truss elements
!!   20 through 29   Plane stress behavior       Membrane elements
!!   30 through 39   Triaxial stress behavior    Solid elements
!!   40 through 49   Plane stress behavior       Shell elements
!!   50 through 59   Uniaxial stress behavior    Beam elements
!!
      MAUX( 0: 9) = (/ 0, 0, 2, 1, 0,15, 0, 0, 0,15 /)
      MAUX(10:19) = (/ 0, 2, 0, 0, 0, 0, 0, 1, 0, 0 /)
      MAUX(20:29) = (/ 0, 4, 8, 0, 0, 6, 0, 6, 0, 0 /)
      MAUX(30:39) = (/ 0, 7, 8, 3, 0, 9,12, 8, 3, 0 /)
      MAUX(40:49) = (/ 0, 4, 0, 0, 0, 6, 0, 8, 0, 0 /)
      MAUX(50:59) = (/ 0, 1, 0, 0, 0, 0, 0, 4, 0, 0 /)
!!
!! NOTE: The restart and rezone modules contain PARAMETER statements defining
!! the maximum number of state variable locations used by any material model.
!! The parameter is MAXSV = 320, where 320 equals 20 x 16 beam locations.
!!
!! MSET(0:60) is an array that contains the number of the property set for
!! each material. The array is only used for printer output.
!!
      MSET( 0: 9) = (/ 0, 1, 2, 3, 0, 5, 6, 7, 8, 9 /)
      MSET(10:19) = (/10,11, 0, 0, 0, 0, 0,17, 0, 0 /)
      MSET(20:29) = (/10,11,22, 0, 0,25, 0,17, 0, 0 /)
      MSET(30:39) = (/10,11,32,33, 0,25,36,17,23, 0 /)
      MSET(40:49) = (/10,11, 0, 0, 0,25, 0,17, 0, 0 /)
      MSET(50:59) = (/10,11, 0, 0, 0,25, 0,17, 0, 0 /)
!!
!! Set all MSHxx counters to zero.
!!
      MSHIF = 0;  MSHQA = 0;  MSHNP = 0;  MSHEL = 0;  MSHHX = 0;
      MSHPX = 0;  MSHTX = 0;  MSHLS = 0;  MSHLX = 0;  MSHLM = 0;
      MSHM4 = 0;  MSHM3 = 0;  MSHTR = 0;  MSHP4 = 0;  MSHP3 = 0;
      MSHBM = 0;  MSHSP = 0;  MSHDM = 0;  MSHSG = 0;  MSHDC = 0;
      MSHTC = 0;  MSHSW = 0;  MSHWC = 0;  MSHBF = 0;  MSHPC = 0;
      MSHFC = 0;  MSHSC = 0;  MSHVC = 0;  MSHCC = 0;  MSHNR = 0;
      MSHND = 0;  MSHFS = 0;  MSHIT = 0;  MSHSI = 0;  MSHSN = 0;
      MSHCE = 0;  MSHCN = 0;  MSHCX = 0;  MSHNS = 0;  MSHNE = 0;
      MSHES = 0;  MSHEE = 0;  MSHSS = 0;  MSHSE = 0;  MSHMT = 0;
      MSHLU = 0;  MSHTF = 0;  MSHFP = 0;  MSHRF = 0;  MSHPV = 0;
      MSHS1 = 0;  MSHS2 = 0;  MSHG1 = 0;  MSHG2 = 0;  MSHG3 = 0;
      MSHMP = 0;  MSHRM = 0;  MSHCM = 0;  MSHIC = 0;  MSHRB = 0;
      MSHRT = 0;  MSHST = 0;  MSHAX = 0;  MSHPP = 0;  MSHNC = 0;
      MSHRE = 0;  MSHID = 0;
!!
!! Set all NRSxx counters to zero.
!!
      NRSIF = 0;  NRSQA = 0;  NRSNP = 0;  NRSEL = 0;  NRSHX = 0;
      NRSPX = 0;  NRSTX = 0;  NRSLS = 0;  NRSLX = 0;  NRSLM = 0;
      NRSM4 = 0;  NRSM3 = 0;  NRSTR = 0;  NRSP4 = 0;  NRSP3 = 0;
      NRSBM = 0;  NRSSP = 0;  NRSDM = 0;  NRSSG = 0;  NRSDC = 0;
      NRSTC = 0;  NRSSW = 0;  NRSWC = 0;  NRSBF = 0;  NRSPC = 0;
      NRSFC = 0;  NRSSC = 0;  NRSVC = 0;  NRSCC = 0;  NRSNR = 0;
      NRSND = 0;  NRSFS = 0;  NRSIT = 0;  NRSSI = 0;  NRSSN = 0;
      NRSCE = 0;  NRSCN = 0;  NRSCX = 0;  NRSNS = 0;  NRSNE = 0;
      NRSES = 0;  NRSEE = 0;  NRSSS = 0;  NRSSE = 0;  NRSMT = 0;
      NRSLU = 0;  NRSTF = 0;  NRSFP = 0;  NRSRF = 0;  NRSPV = 0;
      NRSS1 = 0;  NRSS2 = 0;  NRSG1 = 0;  NRSG2 = 0;  NRSG3 = 0;
      NRSMP = 0;  NRSRM = 0;  NRSCM = 0;  NRSIC = 0;  NRSRB = 0;
      NRSRT = 0;  NRSST = 0;  NRSAX = 0;  NRSPP = 0;  NRSNC = 0;
      NRSRE = 0;  NRSID = 0;
!!
!! Set all IDXxx counters to zero.
!!
      IDXIF = 0;  IDXQA = 0;  IDXNP = 0;  IDXEL = 0;  IDXHX = 0;
      IDXPX = 0;  IDXTX = 0;  IDXLS = 0;  IDXLX = 0;  IDXLM = 0;
      IDXM4 = 0;  IDXM3 = 0;  IDXTR = 0;  IDXP4 = 0;  IDXP3 = 0;
      IDXBM = 0;  IDXSP = 0;  IDXDM = 0;  IDXSG = 0;  IDXDC = 0;
      IDXTC = 0;  IDXSW = 0;  IDXWC = 0;  IDXBF = 0;  IDXPC = 0;
      IDXFC = 0;  IDXSC = 0;  IDXVC = 0;  IDXCC = 0;  IDXNR = 0;
      IDXND = 0;  IDXFS = 0;  IDXIT = 0;  IDXSI = 0;  IDXSN = 0;
      IDXCE = 0;  IDXCN = 0;  IDXCX = 0;  IDXNS = 0;  IDXNE = 0;
      IDXES = 0;  IDXEE = 0;  IDXSS = 0;  IDXSE = 0;  IDXMT = 0;
      IDXLU = 0;  IDXTF = 0;  IDXFP = 0;  IDXRF = 0;  IDXPV = 0;
      IDXS1 = 0;  IDXS2 = 0;  IDXG1 = 0;  IDXG2 = 0;  IDXG3 = 0;
      IDXMP = 0;  IDXRM = 0;  IDXCM = 0;  IDXIC = 0;  IDXRB = 0;
      IDXRT = 0;  IDXST = 0;  IDXAX = 0;  IDXPP = 0;  IDXNC = 0;
      IDXRE = 0;  IDXID = 0;
!!
!! 1. USER's JOB_ID_RECORD: Current and restart file title, date, and time.
!! Set restart sequence counter to "zero." (The "1" becomes a "." later)
!!
      JOB_ID_RECORD%Title_Length            = NPNTL
      JOB_ID_RECORD%PROGRAM                 = ' '
      JOB_ID_RECORD%VERSION                 = ' '
      JOB_ID_RECORD%RELEASE                 = ' '
      JOB_ID_RECORD%CURRENT%TITLE           = ' '
      JOB_ID_RECORD%CURRENT%DATEE           = ' '
      JOB_ID_RECORD%CURRENT%TIMEE           = ' '
      JOB_ID_RECORD%CURRENT%SEQUENCE_NUMBER = 1000
      JOB_ID_RECORD%MESH%TITLE              = ' '
      JOB_ID_RECORD%MESH%DATEE              = ' '
      JOB_ID_RECORD%MESH%TIMEE              = ' '
      JOB_ID_RECORD%RESTART%TITLE           = ' '
      JOB_ID_RECORD%RESTART%DATEE           = ' '
      JOB_ID_RECORD%RESTART%TIMEE           = ' '
      JOB_ID_RECORD%RESTART%SEQUENCE_NUMBER = 1000
!!
!! 2. ANALYSIS CONTROL VARIABLES: Flags controling execution.
!!
      CONTROL%Number_of_Entries = NPNCV
      CONTROL%NAME( 1) = 'SPRINT      '
      CONTROL%NAME( 2) = 'INPUT_ECHO  '
      CONTROL%NAME( 3) = 'WR_RESTART  '
      CONTROL%NAME( 4) = 'RD_RESTART  '
      CONTROL%NAME( 5) = 'WR_MESH     '
      CONTROL%NAME( 6) = 'RD_MESH     '
      CONTROL%NAME( 7) = 'DATA_CHECK  '
      CONTROL%NAME( 8) = 'DT_TAB_FTN  '
      CONTROL%NAME( 9) = 'SUBCYCLING  '
      CONTROL%NAME(10) = 'POLAR_DECOMP'
      CONTROL%NAME(11) = 'MID_INTERVAL'
      CONTROL%NAME(12) = 'BT_PLATE    '
      CONTROL%NAME(13) = 'REZONING    '
      CONTROL%NAME(14) = 'RZ_RESTART  '
!!
!! Set default values for analysis control variables.
!!
      CONTROL%SPRINT  = 0 ! Skip unessential calculations: no
      CONTROL%INECHO  = 0 ! Write input echo: no
      CONTROL%WRSTAR  = 0 ! Write Restart: no
      CONTROL%RDSTAR  = 0 ! Read  Restart: no
      CONTROL%WRMESH  = 0 ! Write Mesh: no
      CONTROL%RDMESH  = 0 ! Read  Mesh: no
      CONTROL%DCHECK  = 0 ! Data check flag: no
      CONTROL%DTTABF  = 0 ! Time step tabulated ftn: none
      CONTROL%SUBCYC  = 0 ! Time integration Subcycling: no
      CONTROL%POLARD  = 0 ! Polar decomposition rotation: no
      CONTROL%MIDINT  = 0 ! Mid-interval gradient evaluation: no
      CONTROL%BTPLTQ  = 0 ! Belytscko-Tsay 4-node plate: no
      CONTROL%REZONE  = 0 ! Rezoning/remeshing processor: no
      CONTROL%RZSTAR  = 0 ! Rezoning/remeshing restart: no
!!
!! 3. SLIDING INTERFACE START-STOP CONTROL: Data to override data on
!! regular sliding interface input record.
!!
      INTERFACE_TIME% Number_of_Entries    = NPNIV
      INTERFACE_TIME% Number_of_Interfaces = NPNIT
      INTERFACE_TIME%NAME(1) = 'INTERFACE_ID'
      INTERFACE_TIME%NAME(2) = 'BEGIN       '
      INTERFACE_TIME%NAME(3) = 'END         '
!$OMP PARALLEL DO
      DO i = 1,NPNIT
        INTERFACE_TIME%SI(i) = SI_type (0,0,0)
      ENDDO
!!
!! 4. NUMERICAL PROCEDURE CONSTANTS: User adjustable coefficients.
!!
      PARAMVALUE%Number_of_Entries = NPNPA
      PARAMVALUE%NAME( 1) = 'LINEAR_Q    '
      PARAMVALUE%NAME( 2) = 'QUADRATIC_Q '
      PARAMVALUE%NAME( 3) = 'HG_VISCOSITY'
      PARAMVALUE%NAME( 4) = 'HG_STIFFNESS'
      PARAMVALUE%NAME( 5) = 'DT_RATIO    '
      PARAMVALUE%NAME( 6) = 'DT_SCALE    '
      PARAMVALUE%NAME( 7) = 'DT_GROW     '
      PARAMVALUE%NAME( 8) = 'STATUS      '
      PARAMVALUE%NAME( 9) = 'ENG_BAL     '
      PARAMVALUE%NAME(10) = 'SI_FACTOR   '
      PARAMVALUE%NAME(11) = 'SI_CAPTURE  '
      PARAMVALUE%NAME(12) = 'SI_BORDER   '
      PARAMVALUE%NAME(13) = 'SI_SORT_FREQ'
      PARAMVALUE%NAME(14) = 'CYCLE_RATIO '
      PARAMVALUE%NAME(15) = 'CYCLE_SPREAD'
      PARAMVALUE%NAME(16) = 'CYCLE_FREQ  '
      PARAMVALUE%NAME(17) = 'MAT32_STEPS '
      PARAMVALUE%NAME(18) = 'NRBC_SCALE  '
      PARAMVALUE%NAME(19) = 'LSOLD_EQ_TOL'
      PARAMVALUE%NAME(20) = 'MAT22_ROTATE'
      PARAMVALUE%NAME(21) = 'BCS_MAX_ITER'
      PARAMVALUE%NAME(22) = 'BCS_GROW_LIM'
      PARAMVALUE%NAME(23) = 'BCS_RELAX   '
      PARAMVALUE%NAME(24) = 'RZ_MAX_CPDOT'
      PARAMVALUE%NAME(25) = 'RZ_MIN_CPDOT'
      PARAMVALUE%NAME(26) = 'RZ_MAX_CPCHG'
      PARAMVALUE%NAME(27) = 'RZ_MIN_CPCHG'
!!
!! Set default user-adjustable parameter values.
!!
      PARAMVALUE%Blk1      = 0.06D+0     ! Linear bulk viscosity coefficeint
      PARAMVALUE%Blk2      = 1.20D+0     ! Quad.  bulk viscosity coefficeint
      PARAMVALUE%HGV       = 1.0D-3   ! Hourglass viscosity
      PARAMVALUE%HGK       = 1.0D-3   ! Hourglass stiffness
      PARAMVALUE%DTratio   = 1.0D-5   ! Limit ratio, DTmin(t)/Dtmin(0)
      PARAMVALUE%DTscale   = 1.0D+0      ! Scaling factor for critical Dt
      PARAMVALUE%DTgrow    = 1.1D+0      ! Step to step growth allowed in DTnext
      PARAMVALUE%STATUS    = 1.0D+5   ! Status request checking interval
      PARAMVALUE%ENG_BAL   = 1000.0D+0   ! Energy balance reporting interval
      PARAMVALUE%Factor    = 1.0D+0      ! Sliding Interface relaxation factor
      PARAMVALUE%Capture   = 0.10D+0     ! Sliding Interface capture distance
      PARAMVALUE%Border    = 1.0001D+0   ! Contact element "exterior" border
      PARAMVALUE%Sort_Freq = 100.0D+0    ! Sliding Interface sorting frequency
      PARAMVALUE%Cycle_R   = 2.0D+0      ! Subcycling group ratio
      PARAMVALUE%Cycle_S   = 1.0D+0      ! Subcycling interface spreading
      PARAMVALUE%Cycle_F   = 100.0D+0    ! Subcycling repartitioning frequency
      PARAMVALUE%Max_Steps = 20.0D+0     ! Material 32 substep upper limit
      PARAMVALUE%NRBC_Q    = 1.05D+0     ! Nonreflecting BC viscous scaling
      PARAMVALUE%EQUI_TOL  = 1.0D-4   ! Layered solid equilibrium tolerance
      PARAMVALUE%FIBROT    = 0.0D+0      ! Material 22 fiber rotation form
      PARAMVALUE%BCSMXIT   = 100.0D+0    ! Beam cross section maximum iterations
      PARAMVALUE%BCSGRLM   = 0.1D+0      ! Beam cross section iterate growth lim
      PARAMVALUE%BCSRLAX   = 0.25D+0     ! Beam cross section under relaxtion
      PARAMVALUE%CPDMAX    = 1.0D-4   ! Maximum rate of change for Cos(Phi)
      PARAMVALUE%CPDMIN    = 1.0D-4   ! Minimum rate of change for Cos(Phi)
      PARAMVALUE%CSPMAX    = 1.0D-4   ! Maximum change for Cos(Phi)
      PARAMVALUE%CSPMIN    = 1.0D-4   ! Minimum change for Cos(Phi)
!!
!! 5. SOUND SPEED: Element longitudinal and transverse sound speeds
!!
      SOUND_SPEED%Density = 0
      SOUND_SPEED%RCL2    = 0
      SOUND_SPEED%RCS2    = 0
!!
!! 6. SIMULATION TIME: Simulation time and associated critical time step info
!!
      TIMSIM%Stop    = 0
      TIMSIM%Total   = 0
      TIMSIM%Step    = 0
      TIMSIM%Cycle   = 0
      TIMSIM%Cymax   = 0
      TIMSIM%DTcntrl = 0
      TIMSIM%DTnext  = 0
      TIMSIM%DTlast  = 0
      TIMSIM%DTmin   = 0
      TIMSIM%DTmax   = 0
      TIMSIM%DTlim   = 0
      TIMSIM%DTsav   = 0
      TIMSIM%SOURCE  = ' '
      TIMSIM%DTHxx   = 0
      TIMSIM%DTPnx   = 0
      TIMSIM%DTTtx   = 0
      TIMSIM%DTLSx   = 0
      TIMSIM%DTM3x   = 0
      TIMSIM%DTM4x   = 0
      TIMSIM%DTTrx   = 0
      TIMSIM%DTP3x   = 0
      TIMSIM%DTP4x   = 0
      TIMSIM%DTBmx   = 0
      TIMSIM%DTSpx   = 0
      TIMSIM%DTDmx   = 0
      TIMSIM%DTSBx   = 0
      TIMSIM%DTDBx   = 0
!$OMP PARALLEL DO
      DO i = 1,NPNDT
        TIMSIM%DTHex(i) = 0
        TIMSIM%DTPen(i) = 0
        TIMSIM%DTTet(i) = 0
        TIMSIM%DTLYS(i) = 0
        TIMSIM%DTMb3(i) = 0
        TIMSIM%DTMb4(i) = 0
        TIMSIM%DTTru(i) = 0
        TIMSIM%DTPl3(i) = 0
        TIMSIM%DTPl4(i) = 0
        TIMSIM%DTBms(i) = 0
        TIMSIM%DTSpr(i) = 0
        TIMSIM%DTDmp(i) = 0
        TIMSIM%DTSBC(i) = 0
        TIMSIM%DTDBC(i) = 0
      ENDDO
      TIMSIM%Ncrit = 0
!$OMP PARALLEL DO
      DO i = 1,NPNDT
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
!! 7. ENERGY BALANCE: Energy balance terms.
!!
      ENERGY%External  = 0
      ENERGY%Internal  = 0
      ENERGY%Sliding   = 0
      ENERGY%Contact   = 0
      ENERGY%Hourglass = 0
      ENERGY%Bulk_Vis  = 0
      ENERGY%Kinetic   = 0
      ENERGY%Eng_Sum   = 0
      ENERGY%Balance   = 0
!!
!! 8. IO_UNIT NUMBERS: Integers assigned to I/O units. The use of terse,
!! lower-case, incomprehensible file names such as "fmasdi" is for
!! compatability with the mono-syllable argot of unix-speak and the
!! unix file system directory structure and file naming conventions.
!!
      IO_UNIT%LSDI =  5 ! fmasdi, simulation data input       (ASCII)
      IO_UNIT%LSRI = 12 ! fmasri, status request input        (ASCII)
      IO_UNIT%LMDI = 13 ! fmamdi, mesh data input             (binary)
      IO_UNIT%LRDI = 14 ! fmardi, restart data input          (binary)

      IO_UNIT%LELO = 21 ! fmaelo, execution log output        (ASCII)
      IO_UNIT%LSDO = 22 ! fmasdo, simulation data output      (ASCII)
      IO_UNIT%LCRO =  6 ! fmacro, computed results output     (ASCII)
      IO_UNIT%LSRO = 24 ! fmasro, status request output       (ASCII)
      IO_UNIT%LMDO = 25 ! fmamdo, mesh data output            (binary)
      IO_UNIT%LRDO = 26 ! fmardo, restart data output         (binary)
      IO_UNIT%LPDB = 27 ! fmapdb, plotting database           (binary)

      IO_UNIT%LRZD = 31 ! FMARZD, rezone restart data, out/in (binary)
!!
!! 9. PRINTED OUTPUT REQUESTS:
!!
      PRINT%Number_of_Entries = NPNPF
      PRINT%NAME( 1) = 'BEGIN '
      PRINT%NAME( 2) = 'END   '
      PRINT%NAME( 3) = 'DELTA '
      PRINT%NAME( 4) = 'TIME  '
      PRINT%NAME( 5) = 'NPT   '
      PRINT%NAME( 6) = 'HXEL  '
      PRINT%NAME( 7) = 'PXEL  '
      PRINT%NAME( 8) = 'TXEL  '
      PRINT%NAME( 9) = 'M3EL  '
      PRINT%NAME(10) = 'M4EL  '
      PRINT%NAME(11) = 'TRUSS '
      PRINT%NAME(12) = 'P3EL  '
      PRINT%NAME(13) = 'P4EL  '
      PRINT%NAME(14) = 'BEAM  '
      PRINT%NAME(15) = 'SPRING'
      PRINT%NAME(16) = 'DAMPER'
!!
      PRINT%Begin  = 0
      PRINT%End    = 0
      PRINT%Delta  = 0
      PRINT%Time   = 0
      PRINT%NODES  = 0
      PRINT%HEXAH  = 0
      PRINT%PENTA  = 0
      PRINT%TETRA  = 0
      PRINT%MEMBT  = 0
      PRINT%MEMBQ  = 0
      PRINT%TRUSS  = 0
      PRINT%PLATT  = 0
      PRINT%PLATQ  = 0
      PRINT%BEAMS  = 0
      PRINT%SPRING = 0
      PRINT%DAMPER = 0
!!
!! 10. OPERATION COUNTERS: Counters for number of times each element is
!! executed.
!!
      COUNTER%HEXAH   = 0
      COUNTER%PENTA   = 0
      COUNTER%TETRA   = 0
      COUNTER%LSOLD   = 0
      COUNTER%MEMBT   = 0
      COUNTER%MEMBQ   = 0
      COUNTER%TRUSS   = 0
      COUNTER%PLATT   = 0
      COUNTER%PLATQ   = 0
      COUNTER%BEAM    = 0
      COUNTER%SPRING  = 0
      COUNTER%DAMPER  = 0
      COUNTER%GAUGE1D = 0
      COUNTER%GAUGE2D = 0
      COUNTER%GAUGE3D = 0
!!
!! 11. INPUT RECORD EORRORS: Counter for number of errors on input records.
!!
      ERROR%COUNT = 0
!!
!! 12. REZONING PROCEDURE CONSTANTS: User adjustable coefficients.
!!
      REZONE%Number_of_Entries = NPNRZ
      REZONE%NAME( 1) = 'CHK_INTERVAL'
      REZONE%NAME( 2) = '<null>      '
      REZONE%NAME( 3) = '<null>      '
      REZONE%NAME( 4) = '<null>      '
      REZONE%NAME( 5) = 'MAX_FLAGS   '
      REZONE%NAME( 6) = 'MAX_FLAG_CYC'
      REZONE%NAME( 7) = 'MAX_LEVELS  '
      REZONE%NAME( 8) = 'WR_MESH     '
      REZONE%NAME( 9) = '<null>      '
      REZONE%NAME(10) = '<null>      '
      REZONE%NAME(11) = '<null>      '
      REZONE%NAME(12) = '<null>      '
      REZONE%NAME(13) = '<null>      '
!!
!! Set default values for user-adjustable rezone coefficients.
!!
      REZONE%INTERVAL    =  50     ! Master time step check interval
      REZONE%TMSTEP      =   0     ! Time step at which first flag set
      REZONE%FCOUNT      =   0     ! Flag counter
      REZONE%FLGCYC      =   0     ! Flag-cycles counter
      REZONE%MAXCNT      =  50     ! Maximum flags before rezone
      REZONE%MAXFCY      = 100     ! Maximum flag-cycles
      REZONE%MAXLEV      =   5     ! Maximum rezone levels
      REZONE%WRMESH      =   0     ! Write "rezone" mesh
      REZONE%FLAG_IS_SET = .FALSE. ! Rezoning control flag
      REZONE%FIRST       = .TRUE.  ! Flag for starting new statistics
      REZONE%REPARTITION = .FALSE. ! Flag to force a repartition
      REZONE%WRDATA      = .FALSE. ! Flag for writing rezone data
      REZONE%FAILED      = .FALSE. ! Flag to indicate array overflow
!!
!! Sub-process time accumulation array.
!!
      DO i =1,100
        SUBPROCESS_TIME(i) = 0
      ENDDO
!!
!! The character variable MSGL must be used to start each line of a message
!! written by USER_MESSAGE. The character variables MSG1, MSG2 and MSGF are
!! internal buffers for formatting data for use in constructing messages.
!!
!! Define Start-of-New-Line character MSGL for constructing messages.
!!
      MSGL = ACHAR(0)
      MSG1 = ' '
      MSG2 = ' '
      MSGF = ' '
!!
      RETURN
      END
!!_
      SUBROUTINE PROGRAM_ID (PROGRAM,VERSION,RELEASE)
!!
!! Copyright (c) by KEY Associates, 4-MAR-1991 21:17:11
!!
!! Purpose: Provide user a record of the current version of the program
!! performing his/her calculations.
!!
      CHARACTER, INTENT(OUT) :: PROGRAM*(*)  ! -/O Program name
      CHARACTER, INTENT(OUT) :: VERSION*(*)  ! -/O Current version of program
      CHARACTER, INTENT(OUT) :: RELEASE*(*)  ! -/O Release date for current version
!!
      PROGRAM = ' FMA-3D'
!!
      VERSION = ' V07.32 '
      RELEASE = 'April 15, 1992'
!!
!! Overhaul input and output file naming and usage to allow more orderly
!! input and output, and to allow an easier port to ubiquitous unix work
!! stations and the limitations of unix file systems.
!!
      VERSION = ' V07.33 '
      RELEASE = 'April 16, 1993'
!!
!! Expand truss/beam cross section property specification to allow either
!! section dimensions or section moments. In either case, the complimentary
!! properties are calculated. Introduce the NASTRAN beam geometry descriptors
!! ROD, TUBE, BAR, BOX and I.
!!
      VERSION = ' V07.34 '
      RELEASE = 'July 25, 1993'
!!
!! Add orthotropic and rubber elasticity material models to beam element.
!! Material models 55 and 17
!!
      VERSION = ' V07.35 '
      RELEASE = 'August 3, 1993'
!!
!! Add arbitrary reference axis location to beam element. (*.NPloc = 5 to
!! read and use *.Yrefloc and *.Zrefloc)
!!
      VERSION = ' V07.36 '
      RELEASE = 'August 21, 1993'
!!
!! Add "directionality" to SPRING and DAMPER elements. (*.Idir=0/1/2/3/4 to
!! get i,j/x/y/z/a orientations) This was done primarily for compatibility
!! with NASTRAN. In finite deformations, *.Idir=1/2/3/4 may not make sense.
!!
      VERSION = ' V07.37 '
      RELEASE = 'August 22, 1993'
!!
!! Rearrange beam calculations so that cross section unit shear vectors are
!! normalized upon first entry rather than every time they are used.
!!
      VERSION = ' V07.38 '
      RELEASE = 'October 3, 1993'
!!
!! Restructure the TABLE_LOOK_UP function so that interval slopes are pre-
!! calculated once, and so that every time a table is used the interval
!! used in the function evaluation is saved for the next time the table is
!! used.
!!
      VERSION = ' V07.39 '
      RELEASE = 'October 15, 1993'
!!
!! Add SEQUENCE_NUMBER to JOB_ID_RECORD to allow "serialization" of restart
!! writes. An "unchanged" restart will initialize its sequence counter with
!! the value found in the restart data.
!!
      VERSION = ' V07.40 '
      RELEASE = 'November 13, 1993'
!!
!! Add mass property calculated results to "line printer" output.
!!
      VERSION = ' V07.41 '
      RELEASE = 'December 19, 1993'
!!
!! Modify NEXT_NP_ID, NEXT_EL_ID and NEXT_SEG_ID to have internal
!! counters thereby allowing the returned index to be modified with-
!! out causing errors in the do-while-enddo indexing performed by
!! these three logical functions.
!!
      VERSION = ' V07.42 '
      RELEASE = 'December 25, 1993'
!!
!! Modify restart file to include subprocess time accumulations;
!! required changes to COMMON.INC and module TIMER.
!!
      VERSION = ' V07.43 '
      RELEASE = 'December 26, 1993'
!!
!! Added periodic boundary condition. Required the creation and
!! management of a new input keyword/record, PERIODBC, a new
!! structure PERIODIC_BC, and a new kinematic routine
!! IMPOSE_PERIODIC_BC. A new read routine for the input record.
!! New initialization code in INITIALIZE_WORKING_STORAGE that
!! checks and re-arranges input data, and generates a rotation
!! matrix (for "CYCLIC" periodic BC's).
!!
      VERSION = ' V07.44 '
      RELEASE = 'April 20, 1994'
!!
!! Changed logo from vertical #-characters to slanted _/-characters.
!!
      VERSION = ' V07.45 '
      RELEASE = 'June 14, 1994'
!!
!! Add spot weld constraint for shell elements. Required the creation
!! and management of a new input keyword/record, SPOTWELD, a new
!! structure SPOT_WELD, and a new kinematic routine IMPOSE_SPOT_WELD.
!! A new read routine for the input record. A new initialization call
!! in INITIALIZE_WORKING_STORAGE to a new module SPOT_WELD_INITILIZATION
!! that converts global location coordinates to element ID's and
!! associated isoparametric coordinate locations.
!!
      VERSION = ' V07.46 '
      RELEASE = 'July 23, 1994'
!!
!! The layered solid had membranes added to it and tested for isotropic
!! and orthotropic materials.
!!
      VERSION = ' V08.01 '
      RELEASE = 'August 16, 1994'
!!
!! Correct double usage of "Wrs" in platq.f and platt.f; change "Wrs"
!! to "Qrs" where it is used for mean membrane in-plane spin. "Wrs" is
!! now used exclusively for the s-derivative of omega_r. The change only
!! affected the Belytschko-Tsay 4-node bending quadrilateral.
!!
      VERSION = ' V08.02 '
      RELEASE = 'September 5, 1994'
!!
!! Fix error in printed output format for SEGMENT_SET header; remove
!! spurious "3X"
!!
      VERSION = ' V08.03 '
      RELEASE = 'September 12, 1994'
!!
!! Replace missing line in numerous calling lists in file "fma2.f"
!!
      VERSION = ' V08.04 '
      RELEASE = 'September 12, 1994'
!!
!! Improve information message in SCAN_RESTART_DATA to look more
!! like file naming conventions: "...Sequence Number #: .rdo.005"
!!
      VERSION = ' V08.05 '
      RELEASE = 'October 1, 1994'
!!
!! Correct simulation data printed output error for SEGMENT echo.
!! Add missing comma to get "output.f" dummy argument for
!! proper linkage to included-file argument list item. Improve
!! node, element and segment set table output.
!!
      VERSION = ' V08.06 '
      RELEASE = 'December 16, 1994'
!!
!! Fix sliding interface initialization error. Set *.Isym to two
!! for single surface contact (*.Type = 1).
!!
      VERSION = ' V08.07 '
      RELEASE = 'December 18, 1994'
!!
!! Add triangles to sliding interface and fix intercept calculation
!! to always use triangles.
!!
      VERSION = ' V08.08 '
      RELEASE = 'February 26, 1995'
!!
!! Streamline SPRING and DAMPER element calculations for orientation.
!!
      VERSION = ' V08.09 '
      RELEASE = 'April 2, 1995'
!!
!! Add Algebraic Input Function Evaluation (AIFE) capability.
!!
      VERSION = ' V08.10 '
      RELEASE = 'April 23, 1995'
!!
!! Correct fiber rotation treatment to add change in coordinate directions
!! from time t(n) to time t(n+1): change (Vn,t) to (Vn,t-Vs,r); fabric
!! material model 22.
!!
      VERSION = ' V08.11 '
      RELEASE = 'June 6, 1995'
!!
!! Correct initialization of PRINT%Time and RESULTS(*)%Time to start with
!! *.Begin before updating the next-time after TIMSIM%Total in module SOLVE.
!!
      VERSION = ' V08.12 '
      RELEASE = 'June 16, 1995'
!!
!! Correct number of record skips in READ_RESTART_DATA when reading state
!! variable data written element-by-element. (Omitted NUMSC and NUMVC in
!! count of records to skip.)
!!
      VERSION = ' V08.13 '
      RELEASE = 'July 30, 1995'
!!
!! Install Hondo III Soil & Crushable Foam material model as material_38.
!!
      VERSION = ' V08.14 '
      RELEASE = 'September 10, 1995'
!!
!! Correct rotation treatment in materials 20/40, 21/41, 25/45 and 27/47
!! to add change in coordinate directions from time t(n) to time t(n+1):
!! change (Wrs) to (Wrs+Vs,r) and change (dBeta) to (dBeta+Vy,x)
!!
      VERSION = ' V08.15 '
      RELEASE = 'October 9, 1995'
!!
!! Add MSHNS, MSHES, MSHSS offsets to do-loops in read-and-build modules
!! for node, element and segment sets to correct error that occurs when
!! mesh input file contains node, element and/or segment sets.
!!
      VERSION = ' V08.16 '
      RELEASE = 'February 20, 1996'
!!
!! Add skip for lupid = -1 in RELINK for PXEL and TXEL input data.
!!
      VERSION = ' V08.17 '
      RELEASE = 'May 20, 1996'
!!
!! Recode pentahedral element (PXEL) based on a proper derivation.
!!
      VERSION = ' V08.18 '
      RELEASE = 'June 28, 1996'
!!
!! Start Version 9: Add midside nodal constraints for adaptivity. Change
!! all INTEGER*2 to INTEGER. Add adaptive rezoning.
!!
      VERSION = ' V09.01 '
      RELEASE = 'October 15, 1995'
!!
!! Add element "border" to sliding interface algorithm. Add external
!! nodal point forces to Plotting_Database output.
!!
      VERSION = ' V09.02 '
      RELEASE = 'February 4, 1996'
!!
!! Change PENTA(*) to TETRA(*) in section of code processing tetrahedrons
!! in module INITIALIZE_WORKING_STORAGE.
!!
      VERSION = ' V09.03 '
      RELEASE = 'September 2, 1996'
!!
!! Add 32 byte character labels to material and set records.
!! Overhaul the reading and building of node sets and segment sets.
!!
      VERSION = ' V09.04 '
      RELEASE = 'November 11, 1996'
!!
!! Modify WRITE_TO_PLOTTING_DATABASE to also write element type with
!! connectivity n-tuple.
!!
      VERSION = ' V09.05 '
      RELEASE = 'December 26, 1996'
!!
!! Correct output initialization of PRINT and RESULTS records during
!! a restart.
!!
      VERSION = ' V09.06 '
      RELEASE = 'January 5, 1997'
!!
!! Expand error checks to distinguish beteen zero and negative nodal
!! point masses, and add checks for zero/negative element volumes.
!!
      VERSION = ' V09.07 '
      RELEASE = 'March 28, 1997'
!!
!! Expand integer output format from I5 to I7 in module TIMER.
!!
      VERSION = ' V09.08 '
      RELEASE = 'April 5, 1997'
!!
!! Correct errors in strain gauge input echo formats; error in
!! REZONE record entry count and initialization; error in plotting
!! results output for membranes (wrong index in accessing the
!! stress array *.STRESS(1:3)).
!!
      VERSION = ' V09.09 '
      RELEASE = 'May 31, 1997'
!!
!! Overhaul 2-D strain gauge integration. Insert omitted angular
!! velocity initialization for quadrilaterals.
!!
      VERSION = ' V09.10 '
      RELEASE = 'June 21, 1997'
!!
!! Placed PARAMETER statements outside of structure definitions, renamed
!! all PARAMETER values to NPxxx, removed a stray "," in a calling list.
!!
      VERSION = ' V09.11 '
      RELEASE = 'September 7, 1997'
!!
!! Converted to Fortran90.
!!
!SPEC_CPU2000      VERSION = ' V10.01 '
      VERSION = 'SPEC CPU'
      RELEASE = 'September 12, 1997'
!!
      RETURN
      END
!!_
      SUBROUTINE INTERCALATION (IOUNIT)
!!
!! Copyright (c) by KEY Associates, 18-APR-1992 09:29:25.85
!! Copyright (c) by KEY Associates; 14-JUN-1994 19:51:36.00
!!
!! Purpose: Identify printed output with program logo, execution time stamp,
!! and program version number.
!!
      USE shared_common_data

      CHARACTER                                                                &
     &          LOGO(10)*100
      INTEGER                                                                  &
     &          IOUNIT

      DATA                                                                     &
     &LOGO( 1)( 1: 50)                                                         &
     & /'         _/_/_/_/_/_/  _/_/          _/_/       _/'/                  &
     &LOGO( 1)(51:100)                                                         &
     & /'_/_/                   _/_/_/_/      _/_/_/_/_/_/ '/                  &
     &LOGO( 2)( 1: 50)                                                         &
     & /'        _/_/_/_/_/_/  _/_/_/      _/_/_/    _/_/_/'/                  &
     &LOGO( 2)(51:100)                                                         &
     & /'_/_/_/             _/_/_/_/_/_/_/   _/_/_/_/_/_/-/'/                  &
     &LOGO( 3)( 1: 50)                                                         &
     & /'       _/_/          _/_/_/_/  _/_/_/_/   _/_/    '/                  &
     &LOGO( 3)(51:100)                                                         &
     & /'  _/_/            _/_/       _/_/  _/_/       _/_/'/                  &
     &LOGO( 4)( 1: 50)                                                         &
     & /'      _/_/          _/_/ _/_/_/_/ _/_/  _/_/      '/                  &
     &LOGO( 4)(51:100)                                                         &
     & /'  _/_/                      _/_/  _/_/        _/_/'/                  &
     &LOGO( 5)( 1: 50)                                                         &
     & /'     _/_/_/_/_/    _/_/  _/_/_/  _/_/  _/_/       '/                  &
     &LOGO( 5)(51:100)                                                         &
     & /' _/_/  _/_/_/_/        _/_/_/_/  _/_/        _/_/ '/
      DATA                                                                     &
     &LOGO( 6)( 1: 50)                                                         &
     & /'    _/_/_/_/_/    _/_/   _/_/   _/_/  _/_/_/_/_/_/'/                  &
     &LOGO( 6)(51:100)                                                         &
     & /'_/_/  _/_/_/_/        _/_/_/_/  _/_/        _/_/  '/                  &
     &LOGO( 7)( 1: 50)                                                         &
     & /'   _/_/          _/_/          _/_/  _/_/_/_/_/_/_'/                  &
     &LOGO( 7)(51:100)                                                         &
     & /'/_/                      _/_/  _/_/        _/_/   '/                  &
     &LOGO( 8)( 1: 50)                                                         &
     & /'  _/_/          _/_/          _/_/  _/_/        _/'/                  &
     &LOGO( 8)(51:100)                                                         &
     & /'_/           _/_/       _/_/  _/_/       _/_/     '/                  &
     &LOGO( 9)( 1: 50)                                                         &
     & /' _/_/          _/_/          _/_/  _/_/        _/_'/                  &
     &LOGO( 9)(51:100)                                                         &
     & /'/           _/_/_/_/_/_/_/   _/_/_/_/_/_/_/       '/                  &
     &LOGO( 10)( 1: 50)                                                        &
     & /'_/_/          _/_/          _/_/  _/_/        _/_/'/                  &
     &LOGO( 10)(51:100)                                                        &
     & /'              _/_/_/_/      _/_/_/_/_/_/          '/

!SPEC_CPU2000      WRITE (IOUNIT,100)
!SPEC_CPU2000     &          (LOGO(i),i=1,10),
!SPEC_CPU2000     &          JOB_ID_RECORD%VERSION,
!SPEC_CPU2000     &          JOB_ID_RECORD%RELEASE,
!SPEC_CPU2000     &          JOB_ID_RECORD%CURRENT%DATE,
!SPEC_CPU2000     &          JOB_ID_RECORD%CURRENT%TIME

      WRITE (IOUNIT,100)                                                       &
     &          (LOGO(i),i=1,10),                                              &
     &          'SPEC CPU',                                                    &
     &          ''

!SPEC_CPU2000 100    FORMAT
!SPEC_CPU2000     &    (
!SPEC_CPU2000     &    //10(3x,A/)
!SPEC_CPU2000     &    /5x,'Version: ',A
!SPEC_CPU2000     &    /5x,'Released: ',A
!SPEC_CPU2000     &    /5x,'Double Precision, Fortran 90'
!SPEC_CPU2000     &    /5x,'This Problem Was Executed On: ',A,'   At ',A
!SPEC_CPU2000     &    )

 100    FORMAT                                                                 &
     &    (                                                                    &
     &    //10(3x,A/)                                                          &
     &    /5x,'Version: ',A                                                    &
     &    /5x,'Released: ',A                                                   &
     &    )

      RETURN
      END
!!_
      SUBROUTINE ALLOCATE_STORAGE
!!
!! Copyright (c) by KEY Associates, 30-AUG-1990 19:58:34
!!
!! Purpose: Zero-out all allocate arrays. Just to insure clean storage,
!! all allocated arrays are set zero. In a couple of cases max available
!! numbers are stored. In several cases, setting an array to zero is not
!! needed, but rather than wonder later if an array was set to zero, ALL
!! allocated storage is set to zero here
!!
      USE shared_common_data;
      USE indx_;           USE node_;           USE tabulated_function_;
      USE beam_;           USE coord_;          USE sliding_interface_;
      USE value_;          USE force_;          USE constrained_node_;
      USE hexah_;          USE penta_;          USE nonreflecting_bc_;
      USE tetra_;          USE lsold_;          USE nodal_point_mass_;
      USE membq_;          USE membt_;          USE rigid_body_mass_;
      USE truss_;          USE platq_;          USE state_variables_;
      USE platt_;          USE motion_;         USE enumerated_sets_;
      USE spring_;         USE damper_;         USE displacement_bc_;
      USE stress_;         USE segment_;        USE contact_surface_;
      USE tied_bc_;        USE results_;        USE relink_scratch_;
      USE gauge1d_;        USE gauge2d_;        USE rigid_wall_bc_;
      USE gauge3d_;        USE massprop_;       USE include_file_;
      USE material_;       USE layering_;       USE sliding_node_;
      USE force_bc_;       USE node_set_;       USE contact_node_;
      USE nrbc_data_;      USE spring_bc_;      USE periodic_bc_;
      USE damper_bc_;      USE spot_weld_;      USE pressure_bc_;
      USE qa_record_;      USE plate_pair_;     USE segment_set_;
      USE body_force_;     USE section_2d_;     USE element_set_;
      USE section_1d_;     USE rigid_body_;     USE plate_pair_;
      USE velocity_ic_;    USE location_;       USE mean_stress_;
      USE output_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered ALLOCATE_STORAGE.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! 1. INCLUDE FILE DATA: File names.
!!
!! At this juncture, INCLUDE_FILE(1:NUMIF) contains the user's input
!! file data and starting & ending addresses of selected input records.
!! The derived data type INCLUDE_FILE should NOT be zero'ed here.
!! This array is set to zero just after it is allocated in FMA_3D.
!!
!! 2. INPUT RECORD PARASING DATA.
!!
!! At this juncture, VALUE(1:MAXRE) has already been allocated and zero'ed.
!!
!! 3. QA RECORD: Captured text lines from users input file.
!!
      IF (NUMQA .GT. 0) THEN
        ALLOCATE (QA_RECORD(1:NUMQA), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'QA_RECORD',NUMQA)
!$OMP PARALLEL DO
        DO i = 1,NUMQA
          QA_RECORD(i) = qa_record_type (MXNQA,' ')
        ENDDO
      ENDIF
!!
!! 4. RESULTS OUTPUT (PLOTTING_DATABASE): Results file control structure.
!!
      IF (NUMRF .GT. 0) THEN
        ALLOCATE (RESULTS(1:NUMRF), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'RESULTS',NUMRF)
!$OMP PARALLEL DO
        DO i = 1,NUMRF
          RESULTS(i) = results_type (' ',0,0,MXNRS,0,0,0,0,0,0,0,0,0,0,        &
     &      0,0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 5. OUTPUT DATA NAMES: Used with input/output operations.
!!
      OUTPUT%Number_of_Entries = MXNRN
!$OMP PARALLEL DO
      DO i = 1,MXNRN
        OUTPUT%NAME(i) = ' '
      ENDDO
!!
!! 6. RESULTS OUTPUT/MASS PROPERTY: Mass property calculation request
!! structure. Note: There is always one mass property; it is used to
!! compute mass properties once at the start of the simulation for the
!! entire model.
!!
      IF (NUMMP .GE. 0) THEN
        ALLOCATE (MASSPROP(0:NUMMP), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'MASSPROP',NUMMP+1)
!$OMP PARALLEL DO
        DO i = 0,NUMMP
          MASSPROP(i) = massprop_type (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,        &
     &      0,0,0,0,(/0,0,0,0,0,0/),0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 7. 1-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      IF (NUMG1 .GT. 0) THEN
        ALLOCATE (GAUGE1D(1:NUMG1), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'GAUGE1D',NUMG1)
!$OMP PARALLEL DO
        DO i = 1,NUMG1
          GAUGE1D(i)%PAR = PAR_gauge1d (0,0,0,0,(/0,0,0,0/),0)
          GAUGE1D(i)%RES = RES_gauge1d (0,0)
        ENDDO
      ENDIF
!!
!! 8. 2-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      IF (NUMG2 .GT. 0) THEN
        ALLOCATE (GAUGE2D(1:NUMG2), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'GAUGE2D',NUMG2)
!$OMP PARALLEL DO
        DO i = 1,NUMG2
          GAUGE2D(i)%PAR = PAR_gauge2d (0,0,0,0,(/0,0,0,0,0,0,0,0/),0,         &
     &      0,0)
          GAUGE2D(i)%RES = RES_gauge2d (0,(/0,0,0/))
        ENDDO
      ENDIF
!!
!! 9. 3-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      IF (NUMG3 .GT. 0) THEN
        ALLOCATE (GAUGE3D(1:NUMG3), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'GAUGE3D',NUMG3)
!$OMP PARALLEL DO
        DO i = 1,NUMG3
          GAUGE3D(i)%PAR = PAR_gauge3d (0,0,0,(/0,0,0,0,0,0,0,0/),0)
          GAUGE3D(i)%RES = RES_gauge3d ((/0,0,0,0,0,0/))
        ENDDO
      ENDIF
!!
!! 10. MATERIAL PROPERTY DATA: Element material property structure.
!! Initialize material bulk and hourglass viscosity and stiffness values.
!!
      IF (NUMMT .GT. 0) THEN
        ALLOCATE (MATERIAL(1:NUMMT), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'MATERIAL',NUMMT)
!$OMP PARALLEL DO
        DO i = 1,NUMMT
          MATERIAL(i) = material_type (0,0,' ',0,0,0,0,0,0,0,0,0,0,0,0,        &
     &      (/(0,j=1,NPNMP)/))
          MATERIAL(i)%PVAL(2) = PARAMVALUE%Blk1
          MATERIAL(i)%PVAL(3) = PARAMVALUE%Blk2
          MATERIAL(i)%PVAL(4) = PARAMVALUE%HGV
          MATERIAL(i)%PVAL(5) = PARAMVALUE%HGK
        ENDDO
      ENDIF
!!
!! 11. MATERIAL PROPERTY DATA NAMES: Used with input/output operations.
!!
!! The derived data type PROPERTY is initialized in INITIALIZE_PROPERTY_NAMES,
!! called once from READ_MATERIAL_PROPERTIES upon first entry.
!!
!! 12. SECTION PROPERTY DATA: Layered solid section structures.
!!
!! The derived data type LAYERING is allocated and initialized in FMA_3D
!! for layered solids (LSOLD elements) and here for solid elements (HEXAH,
!! PENTA and TETRA elements) in the event there not any layered solids.
!!
      IF (.NOT.ALLOCATED(LAYERING) .AND. NUMLU.GT.0) THEN
        ALLOCATE ( LAYERING(1:NUMLU), STAT=IALLOC_FLAG )
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LAYERING',NUMLU )
!$OMP PARALLEL DO
        DO i = 1,NUMLU
          LAYERING(i)%LupID = 0
          LAYERING(i)%Number_of_Layers = MXNLY
          LAYERING(i)%Isys = 0
          DO j = 1,MXNLY
            LAYERING(i)%LayID(j) = 0
            LAYERING(i)%Ltype(j) = 0
            LAYERING(i)%MatID(j) = 0
            LAYERING(i)%H(1:4,j) = (/0,0,0,0/)
            LAYERING(i)%Ax(j)    = 0
            LAYERING(i)%Ay(j)    = 0
            LAYERING(i)%Az(j)    = 0
            LAYERING(i)%Bx(j)    = 0
            LAYERING(i)%By(j)    = 0
            LAYERING(i)%Bz(j)    = 0
          ENDDO
        ENDDO
      ENDIF
!!
!! 13. SECTION PROPERTY DATA: Plate and membrane section structures.
!!
      IF (NUMS2 .GT. 0) THEN
        ALLOCATE (SECTION_2D(1:NUMS2), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SECTION_2D',NUMS2)
!$OMP PARALLEL DO
        DO i = 1,NUMS2
          SECTION_2D(i) = section_2d_type (0,0,0,0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 14. SECTION PROPERTY DATA: Beam and truss section structures.
!!
      IF (NUMS1 .GT. 0) THEN
        ALLOCATE (SECTION_1D(1:NUMS1), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SECTION_1D',NUMS1)
!$OMP PARALLEL DO
        DO i = 1,NUMS1
          SECTION_1D(i) = section_1d_type (0,0,0,0,0,0,0,0,0,0,0,0,0,0,        &
     &      0)
        ENDDO
      ENDIF
!!
!! 15. RIGID BODY: Parameter and results structures.
!!
      IF (NUMRB .GT. 0) THEN
        ALLOCATE (RIGID_BODY(1:NUMRB), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'RIGID_BODY',NUMRB)
!$OMP PARALLEL DO
        DO i = 1,NUMRB
          RIGID_BODY(i) = rigid_body_type (0,0,0,0,0,(/0,0/),(/0,0/),          &
     &      (/0,0/),(/0,0/),(/0,0/),(/0,0/),(/0,0/),0,0,0,0,0,0,0,0,0,         &
     &      0,0,0,0,0,0,0,RESHAPE((/(0,j=1,9)/),(/3,3/)),(/0,0,0/),            &
     &      (/0,0,0/))
        ENDDO
      ENDIF
!!
!! 16. RIGID BODY MASS: Parameter and results structures.
!!
      IF (NUMRM .GT. 0) THEN
        ALLOCATE (RIGID_BODY_MASS(1:NUMRM), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'RIGID_BODY_MASS',NUMRM)
!$OMP PARALLEL DO
        DO i = 1,NUMRM
          RIGID_BODY_MASS(i) = rigid_body_mass_type (0,0,(/0,0,0/),            &
     &      (/0,0,0/),(/0,0,0/),(/0,0,0/),(/0,0,0/),0,RESHAPE(                 &
     &      (/(0,j=1,9)/),(/3,3/)),(/0,0,0/),(/0,0,0/))
        ENDDO
      ENDIF
!!
!! 17. CONCENTRATED NODAL MASS: Parameter and results structures.
!!
      IF (NUMCM .GT. 0) THEN
        ALLOCATE (NODAL_POINT_MASS(1:NUMCM), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'NODAL_POINT_MASS',NUMCM)
!$OMP PARALLEL DO
        DO i = 1,NUMCM
          NODAL_POINT_MASS(i) = nodal_point_mass_type (0,0,0,(/0,0,0/),        &
     &      (/0,0,0/),(/0,0,0/),(/0,0,0/),(/0,0,0/),(/0,0,0/),0,RESHAPE(       &
     &      (/(0,j=1,9)/),(/3,3/)),(/0,0,0/),(/0,0,0/))
        ENDDO
      ENDIF
!!
!! 18. DISPLACEMENT BC
!!
      IF (NUMDC .GT. 0) THEN
        ALLOCATE (DISPLACEMENT_BC(1:NUMDC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'DISPLACEMENT_BC',NUMDC)
!$OMP PARALLEL DO
        DO i = 1,NUMDC
          DISPLACEMENT_BC(i) = displacement_bc_type (0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 19. TIED BC
!!
      IF (NUMTC .GT. 0) THEN
        ALLOCATE (TIED_BC(1:NUMTC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'TIED_BC',NUMTC)
!$OMP PARALLEL DO
        DO i = 1,NUMTC
          TIED_BC(i) = tied_bc_type (0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 20. SPOT WELD
!!
      IF (NUMSW .GT. 0) THEN
        ALLOCATE (SPOT_WELD(1:NUMSW), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SPOT_WELD',NUMSW)
!$OMP PARALLEL DO
        DO i = 1,NUMSW
          SPOT_WELD(i) = spot_weld_type (0,' ',(/0,0,0/),.FALSE.,0,0,0,        &
     &      0,0,0,0,0,0,.FALSE.,(/0,0,0/),(/0,0,0/),(/0,0,0/),(/0,0,0/))
        ENDDO
      ENDIF
!!
!! 21. RIGID WALL BC
!!
      IF (NUMWC .GT. 0) THEN
        ALLOCATE (RIGID_WALL_BC(1:NUMWC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'RIGID_WALL_BC',NUMWC)
!$OMP PARALLEL DO
        DO i = 1,NUMWC
         RIGID_WALL_BC(i)= rigid_wall_bc_type (0,0,(/0,0,0/),(/0,0,0/),        &
     &    (/0,0,0/),(/0,0,0/),0,0,0,0,0,RESHAPE((/(0,j=1,9)/),(/3,3/)),        &
     &    RESHAPE((/(0,j=1,9)/),(/3,3/)),(/0,0,0/),(/0,0,0/),(/0,0,0/),        &
     &    (/0,0,0/),(/0,0,0/),(/0,0,0/),0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 22. BODY FORCE
!!
      IF (NUMBF .GT. 0) THEN
        ALLOCATE (BODY_FORCE(1:NUMBF), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'BODY_FORCE',NUMBF)
!$OMP PARALLEL DO
        DO i = 1,NUMBF
          BODY_FORCE(i) = body_force_type (0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 23. PRESSURE BC
!!
      IF (NUMPC .GT. 0) THEN
        ALLOCATE (PRESSURE_BC(1:NUMPC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'PRESSURE_BC',NUMPC)
!$OMP PARALLEL DO
        DO i = 1,NUMPC
          PRESSURE_BC(i) = pressure_bc_type (0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 24. CONCENTRATED FORCE BC:
!!
      IF (NUMFC .GT. 0) THEN
        ALLOCATE (FORCE_BC(1:NUMFC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'FORCE_BC',NUMFC)
!$OMP PARALLEL DO
        DO i = 1,NUMFC
          FORCE_BC(i) = force_bc_type (0,0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 25. SPRING BC: Nodal point axial/torsional spring restraint.
!!
      IF (NUMSC .GT. 0) THEN
        ALLOCATE (SPRING_BC(1:NUMSC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SPRING_BC',NUMSC)
!$OMP PARALLEL DO
        DO i = 1,NUMSC
          SPRING_BC(i) = spring_bc_type (0,0,0,0,0,0,(/0,0,0/),                &
     &      RES_spring_bc (0,0,0,0))
        ENDDO
      ENDIF
!!
!! 26. DAMPER BC: Nodal point axial/torsional viscous restraint.
!!
      IF (NUMVC .GT. 0) THEN
        ALLOCATE (DAMPER_BC(1:NUMVC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'DAMPER_BC',NUMVC)
!$OMP PARALLEL DO
        DO i = 1,NUMVC
          DAMPER_BC(i) = damper_bc_type (0,0,0,0,0,0,(/0,0,0/),                &
     &      RES_damper_bc (0,0,0,0))
        ENDDO
      ENDIF
!!
!! 27. PERIODIC BC: Periodic boundary condition, "repeated surfaces"
!!
      IF (NUMCC .GT. 0) THEN
        ALLOCATE (PERIODIC_BC(1:NUMCC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'PERIODIC_BC',NUMCC)
!$OMP PARALLEL DO
        DO i = 1,NUMCC
          PERIODIC_BC(i) = periodic_bc_type (0,0,0,0,0,' ',(/0,0,0/),          &
     &      (/0,0,0/),0,0,RESHAPE((/(0,j=1,9)/),(/3,3/)))
        ENDDO
      ENDIF
!!
!! 28. NONREFLECTING BC
!!
      IF (NUMNR .GT. 0) THEN
        ALLOCATE (NONREFLECTING_BC(1:NUMNR), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'NONREFLECTING_BC',NUMNR)
!$OMP PARALLEL DO
        DO i = 1,NUMNR
          NONREFLECTING_BC(i) = nonreflecting_bc_type (0,0,0,0)
        ENDDO
      ENDIF
!!
!! 29. NONREFLECTING BC DATA
!!
!!      DO i = 1,NUMND
!!        NRBC_DATA(i) = nrbc_data_type (0,0,0,0,0)
!!      ENDDO
!!
!! 30. SLIDING INTERFACE
!!
      IF (NUMSI .GT. 0) THEN
        ALLOCATE (SLIDING_INTERFACE(1:NUMSI), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'SLIDING_INTERFACE',NUMSI)
!$OMP PARALLEL DO
        DO i = 1,NUMSI
          SLIDING_INTERFACE(i) = sliding_interface_type (0,0,0,0,0,0,          &
     &      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 31. SLIDING INTERFACE: contact nodes.
!!
!!      DO i = 1,NUMCN
!!        CONTACT_NODE(i) = contact_node_type (0,0,0,0,0,0,0,0,0,0,0)
!!      ENDDO
!!
!! 32. SLIDING INTERFACE: sliding nodes.
!!
!!      DO i = 1,NUMSN
!!        SLIDING_NODE(i) = sliding_node_type (0,0,(/0,0,0,0/),0,0,0,0,
!!      2   0,0,0,0,0,0,
!!      ENDDO
!!
!! 33. SLIDING INTERFACE: contect elements.
!!
!!      DO i = 1,NUMCE
!!        CONTACT_SURFACE(i) = contact_surface_type ((/0,0,0,0/),(/0,0,0,0/),
!!      2   0,0,0,0,0)
!!      ENDDO
!!
!! 34. TABULATED FUNCTION
!!
      IF (NUMTF .GT. 0) THEN
        ALLOCATE (TABULATED_FUNCTION(1:NUMTF), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'TABULATED_FUNCTION',NUMTF)
!$OMP PARALLEL DO
        DO i = 1,NUMTF
          TABULATED_FUNCTION(i)%TFID = 0
          TABULATED_FUNCTION(i)%Number_of_Pairs = MXNTF
          DO j = 1,MXNTF
            TABULATED_FUNCTION(i)%X(j) = 0
            TABULATED_FUNCTION(i)%Y(j) = 0
            TABULATED_FUNCTION(i)%SLOPE(j) = 0
          ENDDO
        ENDDO
      ENDIF
!!
!! 35. NODE SET: set structure.
!!
      IF (NUMNS .GT. 0) THEN
        ALLOCATE (NODE_SET(1:NUMNS), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NODE_SET',NUMNS)
!$OMP PARALLEL DO
        DO i = 1,NUMNS
          NODE_SET(i) = node_set_type (0,0,0,0,0,' ',' ')
        ENDDO
      ENDIF
!!
!! 36. ELEMENT SET: set structure.
!!
      IF (NUMES .GT. 0) THEN
        ALLOCATE (ELEMENT_SET(1:NUMES), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'ELEMENT_SET',NUMES)
!$OMP PARALLEL DO
        DO i = 1,NUMES
          ELEMENT_SET(i) = element_set_type (0,0,0,0,0,' ',' ')
        ENDDO
      ENDIF
!!
!! 37. SEGMENT SET: set structure.
!!
      IF (NUMSS .GT. 0) THEN
        ALLOCATE (SEGMENT_SET(1:NUMSS), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SEGMENT_SET',NUMSS)
!$OMP PARALLEL DO
        DO i = 1,NUMSS
          SEGMENT_SET(i) = segment_set_type (0,0,0,0,0,' ',' ')
        ENDDO
      ENDIF
!!
!! 38. LOCATION: Initialize to "1" the LOCATION pointers for building
!! concatenated node, element and segment ID lists.
!!
      LOCATION = location_type (1,1,1)
!!
!! 39. NODAL POINT DATA/COORDINATES: Nodal point data structures.
!! 40. NODAL POINT DATA/MOTION: Nodal point data structures.
!! 41. NODAL POINT DATA/FORCE: Nodal point data struct
!! Initial conditions on velocity; default is zero. Clear displacement
!! and acceleration to zero.
!!
      IF (NUMRT .GT. 0) THEN
        ALLOCATE (NODE(1:NUMRT), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'NODE',NUMRT)
        ALLOCATE (FORCE(1:NUMRT), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'FORCE',NUMRT)
        ALLOCATE (MOTION(1:NUMRT), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'MOTION',NUMRT)
! Distribute arrays
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        DO i = 1,NUMRT
          NODE(i) = node_type (0,0,0,0,0,0,0,0,0,.FALSE.,0)
          FORCE(i) = force_type (0,0,0,0,0,0)
          MOTION(i) = motion_type (0,0,0,0,0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 42. HEXAHEDRON: Element parameter and results structures.
!!
      IF (NUMHX .GT. 0) THEN
        ALLOCATE (HEXAH(1:NUMHX), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'HEXAH',NUMHX)
!$OMP PARALLEL DO
        DO i = 1,NUMHX
          HEXAH(i) = hexah_type (PAR_hexah(0,0,0,0,(/(0,j=1,8)/),0,0,0),       &
     &      RES_hexah(0,0,(/0,0,0,0,0,0/),(/0,0,0,0/),(/0,0,0,0/),             &
     &      (/0,0,0,0/),(/(0,k=1,8)/),(/(0,k=1,8)/),(/(0,k=1,8)/),             &
     &      0,0,0,0))
        ENDDO
      ENDIF

      IF (NUMLX .GT. 0) THEN
        ALLOCATE (LSHEX(1:NUMLX), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LSHEX',NUMLX)
!$OMP PARALLEL DO
        DO i = 1,NUMLX
          LSHEX(i) = hexah_type (PAR_hexah(0,0,0,0,(/(0,j=1,8)/),0,0,0),       &
     &      RES_hexah(0,0,(/0,0,0,0,0,0/),(/0,0,0,0/),(/0,0,0,0/),             &
     &      (/0,0,0,0/),(/(0,k=1,8)/),(/(0,k=1,8)/),(/(0,k=1,8)/),             &
     &      0,0,0,0))
        ENDDO
      ENDIF
!!
!! 43. PENTAHEDRON: Element parameter and results structures.
!!
      IF (NUMPX .GT. 0) THEN
        ALLOCATE (PENTA(1:NUMPX), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'PENTA',NUMPX)
!$OMP PARALLEL DO
        DO i = 1,NUMPX
          PENTA(i) = penta_type (PAR_penta(0,0,0,0,(/(0,j=1,6)/),0,0,0),       &
     &      RES_penta(0,0,(/0,0,0,0,0,0/),(/0,0,0,0/),(/0,0,0,0/),             &
     &      (/0,0,0,0/),(/(0,k=1,6)/),(/(0,k=1,6)/),(/(0,k=1,6)/),             &
     &      0,0,0,0))
        ENDDO
      ENDIF
!!
!! 44. TETRAHEDRON: Element parameter and results structures.
!!
      IF (NUMTX .GT. 0) THEN
        ALLOCATE (TETRA(1:NUMTX), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'TETRA',NUMTX)
!$OMP PARALLEL DO
        DO i = 1,NUMTX
          TETRA(i) = tetra_type (PAR_tetra(0,0,0,0,(/(0,j=1,4)/),0,0,0),       &
     &      RES_tetra(0,0,(/0,0,0,0,0,0/),(/0,0,0,0/),(/0,0,0,0/),             &
     &      (/0,0,0,0/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 45. LAYERED SOLID: element parameter and results structures.
!! (Already allocated in main program during pre-read of LSOLD.)
!!
!!    IF (NUMLS .GT. 0) THEN
!!      ALLOCATE (LSOLD(1:NUMLS), STAT=IALLOC_FLAG)
!!      CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LSOLD',NUMLS)
!!      DO i = 1,NUMLS
!!        LSOLD(i) = lsold_type (PAR_lsold(0,0,0,(/(0,j=1,8)/),
!!   &      (/(0,j=1,MXLSL)/),0,0),RES_lsold(0,0,RESHAPE(
!!   &      (/(0,j=1,4*MXLSL)/),(/4,MXLSL/)),(/(0,k=1,8)/),
!!   &      (/(0,k=1,8)/),(/(0,k=1,8)/),0,0,0,0))
!!      ENDDO
!!    ENDIF
!!
!! 47. QUADRILATERAL MEMBRANE: Element parameter and results structures.
!!
      IF (NUMM4 .GT. 0) THEN
        ALLOCATE (MEMBQ(1:NUMM4), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'MEMBQ',NUMM4)
!$OMP PARALLEL DO
        DO i = 1,NUMM4
          MEMBQ(i) = membq_type (PAR_membq(0,0,0,0,(/(0,j=1,4)/),0,0,0),       &
     &      RES_membq(0,0,0,(/0,0,0/),0,0,0,(/(0,k=1,4)/),(/(0,k=1,4)/),       &
     &      (/(0,k=1,4)/),0,0,0,0))
        ENDDO
      ENDIF

      IF (NUMLM .GT. 0) THEN
        ALLOCATE (LSMBQ(1:NUMLM), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'LSMBQ',NUMLM)
!$OMP PARALLEL DO
        DO i = 1,NUMLM
          LSMBQ(i) = membq_type (PAR_membq(0,0,0,0,(/(0,j=1,4)/),0,0,0),       &
     &      RES_membq(0,0,0,(/0,0,0/),0,0,0,(/(0,k=1,4)/),(/(0,k=1,4)/),       &
     &      (/(0,k=1,4)/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 49. TRIANGULAR MEMBRANE: Element parameter and results structures.
!!
      IF (NUMM3 .GT. 0) THEN
        ALLOCATE (MEMBT(1:NUMM3), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'MEMBT',NUMM3)
!$OMP PARALLEL DO
        DO i = 1,NUMM3
          MEMBT(i) = membt_type (PAR_membt(0,0,0,0,(/(0,j=1,3)/),0,0,0),       &
     &      RES_membt(0,0,0,(/(0,k=1,3)/),(/(0,k=1,3)/),(/(0,k=1,3)/),         &
     &      (/(0,k=1,3)/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 50. AXIAL FORCE TRUSS: Element parameter and results structures.
!!
      IF (NUMTR .GT. 0) THEN
        ALLOCATE (TRUSS(1:NUMTR), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'TRUSS',NUMTR)
!$OMP PARALLEL DO
        DO i = 1,NUMTR
          TRUSS(i) = truss_type (PAR_truss(0,0,0,0,(/0,0/),0,0,0),             &
     &      RES_truss(0,0,0,(/0,0/),(/0,0/),(/0,0/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 51. QUADRILATERAL PLATE: Element parameter and results structures.
!!
      IF (NUMP4 .GT. 0) THEN
        ALLOCATE (PLATQ(1:NUMP4), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'PLATQ',NUMP4)
! Distribute Arrays
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        DO i = 1,NUMP4
          PLATQ(i) = platq_type (PAR_platq(0,0,0,0,(/(0,j=1,8)/),0,0,0,0),     &
     &      RES_platq(0,0,(/0,0/),(/0,0/),(/0,0/),(/0,0/),(/(0,k=1,8)/),       &
     &      (/(0,k=1,8)/),(/(0,k=1,8)/),0,0,0,0),ADP_platq(.FALSE.,0,0,        &
     &      0,0))
        ENDDO
      ENDIF
!!
!! 52. TRIANGULAR PLATE: Element parameter and results structures.
!!
      IF (NUMP3 .GT. 0) THEN
        ALLOCATE (PLATT(1:NUMP3), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'PLATT',NUMP3)
!$OMP PARALLEL DO
        DO i = 1,NUMP3
          PLATT(i) = platt_type (PAR_platt(0,0,0,0,(/(0,j=1,6)/),0,0,0,0),     &
     &      RES_platt(0,0,(/(0,k=1,6)/),(/(0,k=1,6)/),(/(0,k=1,6)/),           &
     &      0,0,0,0))
        ENDDO
      ENDIF
!!
!! 53. BEAM: Element parameters and results.
!!
      IF (NUMBM .GT. 0) THEN
        ALLOCATE (BEAM(1:NUMBM), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'BEAM',NUMBM)
!$OMP PARALLEL DO
        DO i = 1,NUMBM
          BEAM(i) = beam_type (PAR_beam(0,0,0,0,(/(0,j=1,4)/),0,0,0),          &
     &      RES_beam(0,0,(/0,0,0/),(/(0,k=1,16)/),(/(0,k=1,16)/),              &
     &      0,0,(/(0,k=1,4)/),(/(0,k=1,4)/),(/(0,k=1,4)/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 54. Mean_Stresses (= Stress_Resultants / Area)
!!
      MEAN = mean_stress_type (0,0,0,0,0,0)
!!
!! 55. AXIAL/TORSIONAL SPRING: Element parameters and results.
!!
      IF (NUMSP .GT. 0) THEN
        ALLOCATE (SPRING(1:NUMSP), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SPRING',NUMSP)
!$OMP PARALLEL DO
        DO i = 1,NUMSP
          SPRING(i) = spring_type (PAR_spring(0,0,0,0,(/0,0/),0,0,0,           &
     &      (/0,0,0/)),RES_spring(0,0,0,(/0,0,0/),(/0,0/),(/0,0/),             &
     &      (/0,0/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 56. AXIAL/TORSIONAL DAMPER: Element parameters and results.
!!
      IF (NUMDM .GT. 0) THEN
        ALLOCATE (DAMPER(1:NUMDM), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'DAMPER',NUMDM)
!$OMP PARALLEL DO
        DO i = 1,NUMDM
          DAMPER(i) = damper_type (PAR_damper(0,0,0,0,(/0,0/),0,0,0,           &
     &      (/0,0,0/)),RES_damper(0,0,0,(/0,0,0/),(/0,0/),(/0,0/),             &
     &      (/0,0/),0,0,0,0))
        ENDDO
      ENDIF
!!
!! 57. SEGMENT: Segment parameters
!!
      IF (NUMSG .GT. 0) THEN
        ALLOCATE (SEGMENT(1:NUMSG), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'SEGMENT',NUMSG)
!$OMP PARALLEL DO
        DO i = 1,NUMSG
          SEGMENT(i) = segment_type (PAR_segment(0,0,(/0,0,0,0/)))
        ENDDO
      ENDIF
!!
!! 58. VELOCITY IC: Initial condition data structure.
!!
      IF (NUMIC .GT. 0) THEN
        ALLOCATE (VELOCITY_IC(1:NUMIC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT (IALLOC_FLAG,'VELOCITY_IC',NUMIC)
!$OMP PARALLEL DO
        DO i = 1,NUMIC
          VELOCITY_IC(i) = velocity_ic_type (0,0,0,0,0,0,0,0)
        ENDDO
      ENDIF
!!
!! 59. REZONING: Plate pairs that share a common boundary.
!! This array is set to zero just after it is allocated in FMA_3D.
!!
!!      DO i = 1,NUMPP
!!        PLATE_PAIR(i) = plate_pair_type (0,IBIDIS(0,0),IBIDIS(0,0)0,0,0)
!!      ENDDO
!!
!! 60. REZONING: Constraints for midside nodes.
!!
      IF (NUMNC .GT. 0) THEN
        ALLOCATE (CONSTRAINED_NODE(1:NUMNC), STAT=IALLOC_FLAG)
        CALL REPORT_ALLOCATION_EVENT                                           &
     &    (IALLOC_FLAG,'CONSTRAINED_NODE',NUMNC)
!$OMP PARALLEL DO
        DO i = 1,NUMNC
          CONSTRAINED_NODE(i) = constrained_node_type (0,0,(/0,0/),            &
     &      (/0,0/),(/0,0/),0)
        ENDDO
      ENDIF
!!
!! 61. Concatenated ID lists for enumerated node, element and segment sets.
!! These arrays are set to zero just after they are allocated.
!!
!!      DO i = 1,NUMNE
!!        NNPSETS(i) = 0
!!      ENDDO
!!      DO i = 1,NUMEE
!!        NELSETS(i) = 0
!!      ENDDO
!!      DO i = 1,NUMSE
!!        NSGSETS(i) = 0
!!      ENDDO
!!
!! 62. Concatenated shell element thru-the-thickness stress components.
!!
!!      DO i = 1,NUMST
!!        STRESS(1:6,i) = (/0,0,0,0,0,0/)
!!      ENDDO
!!
!! 63. Concatenated material state variable storage.
!!
!!      DO i = 1,NUMAX
!!        STATE_VARIABLES(i) = 0
!!      ENDDO
!!
!! 64. Coordinate sorting arrays for sliding interface contact searchs.
!!
!!      DO i = 1,NUMCX
!!        COORD(i,1) = 0
!!        COORD(i,2) = 0
!!        COORD(i,3) = 0
!!      ENDDO
!!
!! 65. ID sorting arrays for sliding interface contact searchs.
!!
!!      DO i = 1,NUMCX
!!        INDX(i,1) = 0
!!        INDX(i,2) = 0
!!        INDX(i,3) = 0
!!        INDX(i,4) = 0
!!        INDX(i,5) = 0
!!        INDX(i,6) = 0
!!        INDX(i,7) = 0
!!        INDX(i,8) = 0
!!      ENDDO
!!
!! 66. Relinking scratch arrays for node, element and segment ID sorting.
!!
!!      DO i = 1,NUMNP
!!        INPV(i) = 0
!!        JNPV(i) = 0
!!      ENDDO
!!      DO i = 1,NUMEL
!!        IELV(i) = 0
!!        JELV(i) = 0
!!      ENDDO
!!      DO i = 1,NUMSG
!!        ISGV(i) = 0
!!        JSGV(i) = 0
!!      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE REPORT_ALLOCATION_EVENT (IALLOC_FLAG,DATA_TYPE,NUMBER)
!!
!! Copyright (c) by KEY Associates;  5-OCT-1997 15:11:02.00
!!
      USE shared_common_data
!!
      INTEGER         :: IALLOC_FLAG ! I/- ALLOCATE function status flag.
      CHARACTER(*)    :: DATA_TYPE   ! I/- Derived data type that failed.
      INTEGER         :: NUMBER      ! I/- Number of items allocated
!!
      L = LEN (DATA_TYPE)
      IF ( IALLOC_FLAG .EQ. 0) THEN
        WRITE (MSG1,'(I8)') NUMBER
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'REPORT_ALLOCATION_EVENT.001.00'//                       &
     &          MSGL//'ALLOCATE Function Succeeded For Derived Type: '         &
     &              //DATA_TYPE(1:L)//                                         &
     &          MSGL//'Number Of Items Allocated: '//MSG1                      &
     &          )
      ELSE
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'REPORT_ALLOCATION_EVENT.002.00'//                       &
     &          MSGL//'ALLOCATE Function Failed For Derived Type: '            &
     &              //DATA_TYPE(1:L)//                                         &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE INIT_INCLUDE_FILE
!!
!! Copyright (c) by KEY Associates, 20-JAN-1992 20:41:53
!!
!! Purpose: Initialize first entry in INCLUDE_FILE with user's root input
!! file.
!!
      USE shared_common_data
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          FORMATTED_READ*8,                                              &
     &          FULL_NAME*1024
      LOGICAL                                                                  &
     &          EXSTAT,                                                        &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered INIT_INCLUDE_FILE.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Inquire about file's existence and obtain full file name.
!!
      ERROR%COUNT = 0
      IOERROR = .TRUE.
      INQUIRE                                                                  &
     &          (                                                              &
     &          FILE        = 'fma3d.in',                                      &
     &          EXIST       = EXSTAT,                                          &
     &          FORMATTED   = FORMATTED_READ,                                  &
     &          NAME        = FULL_NAME,                                       &
     &          ERR         = 100                                              &
     &          )
      IOERROR = .FALSE.
 100    CONTINUE
      FULL_NAME = 'fma3d.in'
!!
!! Warning for failed INQUIRE operation.
!!
      IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'INIT_INCLUDE_FILE.001.01'//                             &
     &          MSGL//'Error Executing INQUIRE On: '//'fma3d.in'               &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
!!
!! Warning for non-existant file.
!!
      ELSE IF (.NOT.EXSTAT) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'INIT_INCLUDE_FILE.001.02'//                             &
     &          MSGL//'INQUIRE Unable To Find: '//'fma3d.in'                   &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
!!
!! Warning for file that isn't structured for formatted reads.
!!
!!!   ELSE IF (INDEX(FORMATTED_READ,'NO') .NE. 0) THEN
!!!     CALL USER_MESSAGE
!!!  &       (
!!!  &       MSGL//'WARN'//
!!!  &       MSGL//'INIT_INCLUDE_FILE.001.03'//
!!!  &       MSGL//'Unable To Perform Formatted Read On: '//TRIM (FULL_NAME)
!!!  &       )
!!!     ERROR%COUNT = ERROR%COUNT + 1
!!
      ENDIF
!!
!! Print total number of processing errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INIT_INCLUDE_FILE.002.00'//                             &
     &          MSGL//'Total Number Of File Name Errors:'//MSG1//              &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      INCLUDE_FILE(1)%Full_Name     = TRIM (FULL_NAME)
      INCLUDE_FILE(1)%File_Number   = 1
      INCLUDE_FILE(1)%Line_Number   = 0
      INCLUDE_FILE(1)%LAYUP_Begin   = 0
      INCLUDE_FILE(1)%LAYUP_End     = 0
      INCLUDE_FILE(1)%LSOLD_Begin   = 0
      INCLUDE_FILE(1)%LSOLD_End     = 0
      INCLUDE_FILE(1)%CONTROL_Begin = 0
      INCLUDE_FILE(1)%CONTROL_End   = 0
!!
      RETURN
      END
!!_
      SUBROUTINE SCAN_INCLUDED_FILES (MAXIF)
!!
!! Copyright (c) by KEY Associates, 19-JAN-1992 10:29:03
!!
!! Purpose: Open and scan the included files counting the storage needed to
!! conduct the calculation. Verify any new included files found.
!!
      USE shared_common_data
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: MAXIF
!!
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SCAN_INCLUDED_FILES.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Loop over included files. Each new INCLUDE record will extent the do-while.
!!
      NIF = 1
      NUMIF = 1
      ERROR%COUNT = 0
      DO WHILE (NIF .LE. NUMIF)
!!
!! Open included file and count storage needed to conduct the calculation.
!!
        IOERROR = .TRUE.
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(NIF)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
        IOERROR = .FALSE.
 200      CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
        IF (IOERROR) THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SCAN_INCLUDED_FILE.001.01'//                            &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of scan of included input record file.
!!
        ELSE
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'SCAN_INCLUDED_FILES.001.02'//                           &
     &          MSGL//'INCLUDE File Scanned: '                                 &
     &              //'fma3d.in'                                               &
     &          )
!!
!! Pre-read included file for the presence of algebraic input functions
!! that will require evaluation before the program's regular record
!! read takes place. (Algebraic input functions are set-off with curly
!! brackets: {...}. If any algebraic functions are found, a scratch
!! file will be created (with the algebraic functions evaluated) and used
!! in place of the original file. If no curly brackets are found, the
!! call to AIFE_FACILITY is a do-nothing.)
!!
!!          CALL AIFE_FACILITY (NIF)
!!
!! Scan included file input records. Save count of INCLUDE records in file
!! for use in comparing to count of INCLUDE records after scanning of the
!! included file.
!!
          NUMIF_OLD = NUMIF
          CALL SCAN_INPUT_RECORDS (NIF)
!!
!! Verify any new included files found.
!!
          IF (NUMIF .GT. NUMIF_OLD) THEN
            IF (NUMIF .LE. MAXIF) THEN
              DO i = NUMIF_OLD+1,NUMIF
                CALL VERIFY_INCLUDED_FILE (i)
                INCLUDE_FILE(i)%LAYUP_Begin   = 0
                INCLUDE_FILE(i)%LAYUP_End     = 0
                INCLUDE_FILE(i)%LSOLD_Begin   = 0
                INCLUDE_FILE(i)%LSOLD_End     = 0
                INCLUDE_FILE(i)%CONTROL_Begin = 0
                INCLUDE_FILE(i)%CONTROL_End   = 0
              ENDDO
!!
!! Error message for either too many included files or include-file infinite
!! loop.
!!
            ELSE
              WRITE (MSG1,'(I8)') MAXIF
              CALL USER_MESSAGE                                                &
     & (                                                                       &
     & MSGL//'FATAL'//                                                         &
     & MSGL//'FMA_3D.001.00'//                                                 &
     & MSGL//'The Maximum Number of Included Files Exceeded:'//MSG1//          &
     & MSGL//'Either Correct "Infinite Include Loop," Or'//                    &
     & MSGL//'Ask Program Adminstrator To Increase MAXIF.'                     &
     & )
!!
!! End of was-NUMIF-below-MAXIF if-test.
!!
            ENDIF
!!
!! End of was-NUMIF-greater-than-NUMIF_OLD if-test.
!!
          ENDIF
!!
!! Disconnect included file.
!!
          CLOSE (UNIT=IO_UNIT%LSDI, STATUS='KEEP')
!!
!! End of was-open-sucessful if-test.
!!
        ENDIF
!!
!! End of do-while loop on included files.
!!
        NIF = NIF + 1
      ENDDO
!!
!! Print total number of INCLUDE record processing errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SCAN_INCLUDED_FILES.004.00'//                           &
     &          MSGL//'Total Number Of INCLUDE Record Errors:'//MSG1//         &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE SCAN_INPUT_RECORDS (NIF)
!!
!! Copyright (c) by KEY Associates, 11-JUN-1989 11:30:44
!!
!! Purpose: Makes a pass through the input record file counting the storage
!! needed to conduct the calculation. A number of values from the input
!! records are obtained at this time.
!!
      USE shared_common_data
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NIF
!!
!! Local variables.
      INTEGER   :: ILCOUNT   ! S/S Input file line counter
      INTEGER   :: ILINDEX   ! S/S Start-of-Record line in input file
      CHARACTER :: TEXT*80,C_VALUE*32,KEY_WORD*12
      LOGICAL   :: EOF,RERROR

      LOGICAL, SAVE :: FIRST = .TRUE.

      EXTERNAL                                                                 &
     &          C_VALUE
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SCAN_INPUT_RECORDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Initialize counters to zero.
!!
      IF (FIRST) THEN
        NUMQA = 0     ! QA records
        NUMNP = 0     ! Nodal points
        NUMEL = 0     ! Elements, all
        NUMHX = 0     ! 8-Node hexahedron solids
        NUMPX = 0     ! 6-Node pentahedron solids
        NUMTX = 0     ! 4-Node tetrahedron solids
        NUMLS = 0     ! Layered solids
        NUMLX = 0     ! 8-Node hexahedrons for layered solids
        NUMLM = 0     ! 4-Node membranes for layered solids
        NUMM4 = 0     ! 4-Node quadrilateral membranes
        NUMM3 = 0     ! 3-Node triangular membranes
        NUMTR = 0     ! 2-Node trusses, neither bending nor twisting
        NUMP4 = 0     ! 4-Node quadrilateral plates
        NUMP3 = 0     ! 3-Node triangular plates
        NUMBM = 0     ! 2-Node beams
        NUMSP = 0     ! 2-Node springs
        NUMDM = 0     ! 2-Node dampers
        NUMSG = 0     ! Segments
        NUMDC = 0     ! Displacement BC's
        NUMTC = 0     ! Tied BC's
        NUMSW = 0     ! Spot welds
        NUMWC = 0     ! Rigid wall BC's
        NUMBF = 0     ! Body force fields
        NUMPC = 0     ! Pressure BC's
        NUMFC = 0     ! Concentrated force BC's
        NUMSC = 0     ! "Spring" BC's
        NUMVC = 0     ! "Damper" BC's
        NUMCC = 0     ! Periodic BC's
        NUMNR = 0     ! Non-reflecting BC's
        NUMND = 0     ! Non-reflecting BC data entries
        NUMFS = 0     ! Fastners
        NUMIT = 0     ! Interface times
        NUMSI = 0     ! Sliding interfaces
        NUMNS = 0     ! Node sets
        NUMNE = 0     ! Node set entries
        NUMES = 0     ! Element sets
        NUMEE = 0     ! Element set entries
        NUMSS = 0     ! Segment sets
        NUMSE = 0     ! Segment set entries
        NUMMT = 0     ! Materials, all (constitutive) models
        NUMLU = 0     ! Lay-up cross sections
        NUMTF = 0     ! Tabulated functions (history functions, etc.)
        NUMFP = 1     ! Maximum function x,y-pairs
        NUMRF = 0     ! Results files
        NUMPV = 0     ! Parameter value changes
        NUMS1 = 0     ! Beam/truss section properties
        NUMS2 = 0     ! Plate section properties
        NUMG1 = 0     ! 1-D strain gauge
        NUMG2 = 0     ! 2-D strain gauge
        NUMG3 = 0     ! 3-D strain gauge
        NUMMP = 0     ! Mass property requests
        NUMRM = 0     ! Substitute/added rigid body mass
        NUMCM = 0     ! Concentrated nodal point masses
        NUMIC = 0     ! Initial conditions on velocity
        NUMRB = 0     ! Rigid body domains
        NUMRT = 0     ! Nodal points + nodes with rotational dof's
        NUMIF = 1     ! Include files (The user's root file is 1.)
        NUMPP = 0     ! Plate pairs (element-element interfaces)
        NUMNC = 0     ! Nodal constraints for midside nodes
        MXEID = 0     ! Maximum element ID used.
        MAXRE = 0     ! Maximum number of record entries found.
        FIRST = .FALSE.
      ENDIF
!!
!! Read input records and look for keywords. Develop counts on the basis
!! of keywords, and node and element numbering. (The module GETIRC only
!! retrieves the key word and counts the total number of entries in the
!! record including the keyword.)
!!
      REWIND (IO_UNIT%LSDI)
      IERR = 0
      ILCOUNT = 0
      EOF = .FALSE.
      DO WHILE (.NOT.EOF)
        CALL GETIRC                                                            &
     & (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
        MAXRE = MAX (MAXRE,NVALS)
!!
!! Compare key word read with known key words. Note: "SPRINGBC" must
!! occur before "SPRING" and DAMPERBC" must occur before "DAMPER" in
!! the sequence of checks for the count to be correct.
!!
        KEY_WORD = C_VALUE(1)
!!
        IF (INDEX(KEY_WORD,'NPT'           ) .NE. 0) THEN
          NUMNP = NUMNP + 1
        ELSE IF (INDEX(KEY_WORD,'TIEDBC'   ) .NE. 0) THEN
          NUMTC = NUMTC + 1
        ELSE IF (INDEX(KEY_WORD,'SPOTWELD' ) .NE. 0) THEN
          NUMSW = NUMSW + 1
        ELSE IF (INDEX(KEY_WORD,'WALLBC'   ) .NE. 0) THEN
          NUMWC = NUMWC + 1
        ELSE IF (INDEX(KEY_WORD,'BODYFORCE') .NE. 0) THEN
          NUMBF = NUMBF + 1
        ELSE IF (INDEX(KEY_WORD,'PRESSBC'  ) .NE. 0) THEN
          NUMPC = NUMPC + 1
        ELSE IF (INDEX(KEY_WORD,'FORCEBC'  ) .NE. 0) THEN
          NUMFC = NUMFC + 1
        ELSE IF (INDEX(KEY_WORD,'SPRINGBC' ) .NE. 0) THEN
          NUMSC = NUMSC + 1
        ELSE IF (INDEX(KEY_WORD,'DAMPERBC' ) .NE. 0) THEN
          NUMVC = NUMVC + 1
        ELSE IF (INDEX(KEY_WORD,'PERIODBC' ) .NE. 0) THEN
          NUMCC = NUMCC + 1
        ELSE IF (INDEX(KEY_WORD,'NRBC'     ) .NE. 0) THEN
          NUMNR = NUMNR + 1
        ELSE IF (INDEX(KEY_WORD,'INTERFACE') .NE. 0) THEN
          NUMIT = NUMIT + 1
        ELSE IF (INDEX(KEY_WORD,'SLIDE'    ) .NE. 0) THEN
          NUMSI = NUMSI + 1
        ELSE IF (INDEX(KEY_WORD,'HXEL'     ) .NE. 0) THEN
          NUMHX = NUMHX + 1
        ELSE IF (INDEX(KEY_WORD,'PXEL'     ) .NE. 0) THEN
          NUMPX = NUMPX + 1
        ELSE IF (INDEX(KEY_WORD,'TXEL'     ) .NE. 0) THEN
          NUMTX = NUMTX + 1
        ELSE IF (INDEX(KEY_WORD,'LSEL'     ) .NE. 0) THEN
          NUMLS = NUMLS + 1
          IF (INCLUDE_FILE(NIF)%LSOLD_Begin .EQ. 0) THEN
            INCLUDE_FILE(NIF)%LSOLD_Begin = ILINDEX
          ENDIF
          INCLUDE_FILE(NIF)%LSOLD_End = ILINDEX
        ELSE IF (INDEX(KEY_WORD,'M4EL'     ) .NE. 0) THEN
          NUMM4 = NUMM4 + 1
        ELSE IF (INDEX(KEY_WORD,'M3EL'     ) .NE. 0) THEN
          NUMM3 = NUMM3 + 1
        ELSE IF (INDEX(KEY_WORD,'TRUSS'    ) .NE. 0) THEN
          NUMTR = NUMTR + 1
        ELSE IF (INDEX(KEY_WORD,'P4EL'     ) .NE. 0) THEN
          NUMP4 = NUMP4 + 1
        ELSE IF (INDEX(KEY_WORD,'P3EL'     ) .NE. 0) THEN
          NUMP3 = NUMP3 + 1
        ELSE IF (INDEX(KEY_WORD,'BEAM'     ) .NE. 0) THEN
          NUMBM = NUMBM + 1
        ELSE IF (INDEX(KEY_WORD,'SPRING'   ) .NE. 0) THEN
          NUMSP = NUMSP + 1
        ELSE IF (INDEX(KEY_WORD,'DAMPER'   ) .NE. 0) THEN
          NUMDM = NUMDM + 1
        ELSE IF (INDEX(KEY_WORD,'SEGMENT'  ) .NE. 0) THEN
          NUMSG = NUMSG + 1
        ELSE IF (INDEX(KEY_WORD,'DISPBC'   ) .NE. 0) THEN
          NUMDC = NUMDC + 1
        ELSE IF (INDEX(KEY_WORD,'NPSET'    ) .NE. 0) THEN
          NUMNS = NUMNS + 1
        ELSE IF (INDEX(KEY_WORD,'ELSET'    ) .NE. 0) THEN
          NUMES = NUMES + 1
        ELSE IF (INDEX(KEY_WORD,'SEGSET'   ) .NE. 0) THEN
          NUMSS = NUMSS + 1
        ELSE IF (INDEX(KEY_WORD,'MATERIAL' ) .NE. 0) THEN
          NUMMT = NUMMT + 1
        ELSE IF (INDEX(KEY_WORD,'LAYUP'    ) .NE. 0) THEN
          NUMLU = NUMLU + 1
          IF (INCLUDE_FILE(NIF)%LAYUP_Begin .EQ. 0) THEN
            INCLUDE_FILE(NIF)%LAYUP_Begin = ILINDEX
          ENDIF
          INCLUDE_FILE(NIF)%LAYUP_End = ILINDEX
        ELSE IF (INDEX(KEY_WORD,'TABFTN'   ) .NE. 0) THEN
          NUMTF = NUMTF + 1
          NUMFP = MAX (NUMFP,(NVALS-2)/2)
        ELSE IF (INDEX(KEY_WORD,'RESULTS'  ) .NE. 0) THEN
          NUMRF = NUMRF + 1
        ELSE IF (INDEX(KEY_WORD,'PARAM'    ) .NE. 0) THEN
          NUMPV = NUMPV + 1
        ELSE IF (INDEX(KEY_WORD,'CONTROL'  ) .NE. 0) THEN
          IF (INCLUDE_FILE(NIF)%CONTROL_Begin .EQ. 0) THEN
            INCLUDE_FILE(NIF)%CONTROL_Begin = ILINDEX
          ENDIF
          INCLUDE_FILE(NIF)%CONTROL_End = ILINDEX
        ELSE IF (INDEX(KEY_WORD,'BSECTION' ) .NE. 0) THEN
          NUMS1 = NUMS1 + 1
        ELSE IF (INDEX(KEY_WORD,'PSECTION' ) .NE. 0) THEN
          NUMS2 = NUMS2 + 1
        ELSE IF (INDEX(KEY_WORD,'GAUGE1D'  ) .NE. 0) THEN
          NUMG1 = NUMG1 + 1
        ELSE IF (INDEX(KEY_WORD,'GAUGE2D'  ) .NE. 0) THEN
          NUMG2 = NUMG2 + 1
        ELSE IF (INDEX(KEY_WORD,'GAUGE3D'  ) .NE. 0) THEN
          NUMG3 = NUMG3 + 1
        ELSE IF (INDEX(KEY_WORD,'MASSPROP' ) .NE. 0) THEN
          NUMMP = NUMMP + 1
        ELSE IF (INDEX(KEY_WORD,'VELIC'    ) .NE. 0) THEN
          NUMIC = NUMIC + 1
        ELSE IF (INDEX(KEY_WORD,'RBODY'    ) .NE. 0) THEN
          NUMRB = NUMRB + 1
        ELSE IF (INDEX(KEY_WORD,'RBMASS'   ) .NE. 0) THEN
          NUMRM = NUMRM + 1
        ELSE IF (INDEX(KEY_WORD,'NPMASS'   ) .NE. 0) THEN
          NUMCM = NUMCM + 1
        ELSE IF (INDEX(KEY_WORD,'QAREC'    ) .NE. 0) THEN
          NUMQA = NUMQA + 1
        ELSE IF (INDEX(KEY_WORD,'NPCON1'   ) .NE. 0) THEN
          NUMNC = NUMNC + 1
        ELSE IF (INDEX(KEY_WORD,'INCLUDE'  ) .NE. 0) THEN
          NUMIF = NUMIF + 1
          INCLUDE_FILE(NUMIF)%File_Number = NIF
          INCLUDE_FILE(NUMIF)%Line_Number = ILINDEX
        ELSE IF (INDEX(KEY_WORD,'ENDDATA'  ) .NE. 0) THEN
          EOF = .TRUE.
        ENDIF
      ENDDO
      ERROR%COUNT = IERR
!!
!! Add up the number of elements in the model. (The extra hexa's NUMLX
!! and 4-node membranes NUMLM needed to process the layered solids will
!! be counted and added to the total number of elements later.)
!!
      NUMEL = NUMHX + NUMPX + NUMTX   & ! Solid elements           
     &      + NUMLS                   & ! Layered solid elements   
     &      + NUMM4 + NUMM3           & ! Membrane elements        
     &      + NUMTR                   & ! Truss elements (axial force only)
     &      + NUMP4 + NUMP3           & ! Plate elements           
     &      + NUMBM                   & ! Beam elements            
     &      + NUMSP + NUMDM             ! Spring and Damper elements
!!
!! Determine if rotational degrees of freedom are present. This is an
!! esitimate of the number of locations needed since the connectivity
!! is not yet known.
!!
      IF (NUMP4+NUMP3+NUMBM+NUMSP+NUMDM .GT. 0) THEN
        NUMRT = MIN                                                            &
     &          (                                                              &
     &          (NUMNP + 4*NUMP4 + 3*NUMP3 + 2*NUMBM + NUMSP + NUMDM),         &
     &          (NUMNP + NUMNP)                                                &
     &          )
      ELSE
        NUMRT = NUMNP
      ENDIF
!!
!! Record time spent in scanning input records.
!!
!SPEC_CPU2000      CALL TIMER (1)
!!
!! Print total number of input record format errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SCAN_INPUT_RECORDS.001.00'//                            &
     &          MSGL//'Total Number Of Format Errors In Input:'//MSG1//        &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE VERIFY_INCLUDED_FILE (NIF)
!!
!! Copyright (c) by KEY Associates, 19-JAN-1992 10:29:03
!!
!! Purpose: Ascertain that included file exists and is "proper."
!!
      USE shared_common_data
      USE value_
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      INTEGER, INTENT(IN) :: NIF ! Included file index
!!
      CHARACTER                  &
     &          FNEC*1,          &  ! -/- File name enclosing character
     &          TEXT*80,         &  ! S/- First line of input record
     &          C_VALUE*32,      &  ! F/- Function that returns character value
     &          KEY_WORD*12,     &  ! -/- Key word from record
     &          FULL_NAME*1024,  &
     &          FORMATTED_READ*8

      INTEGER :: ILCOUNT   ! S/S Input file line counter
      INTEGER :: ILINDEX   ! S/S Start-of-Record line in input file

      EXTERNAL                                                                 &
     &          C_VALUE
      LOGICAL                                                                  &
     &          EOF,                                                           &
     &          EXSTAT,                                                        &
     &          IOERROR
      DATA                                                                     &
     &          FNEC /'"'/
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered VERIFY_INCLUDED_FILE.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Read and verify an INCLUDE record. The INCLUDE record has the form:
!! INCLUDE<delimiter>"file_name" (yes, inbetween double quotes)
!!
      IERR = 0
      ILCOUNT = 0
      EOF = .FALSE.
!!
!! Position input file to INCLUDE record.
!!
      REWIND (IO_UNIT%LSDI)
      ILINDEX = INCLUDE_FILE(NIF)%Line_Number
      DO i = 1,ILINDEX - 1
        ILCOUNT = ILCOUNT + 1
        READ (IO_UNIT%LSDI,100)
 100    FORMAT ()
      ENDDO
!!
!! Read INCLUDE record.
!!
      NVALS = 2
      CALL GETIRV                                                              &
     &(NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
      KEY_WORD = C_VALUE(1)
!!
      IF (INDEX(KEY_WORD,'INCLUDE') .NE. 0) THEN
!!
!! Locate included file name from data on the INCLUDE record. The INCLUDE
!! record has the following form: INCLUDE<delimiter>"file_name" Thus, at
!! a minimum VALUE(2)%VTYP should equal 'U' for unknown record entry type
!! due to the presence of apostrophies.
!!
        IF (VALUE(2)%VTYP .EQ. 'N') THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'VERIFY_INCLUDED_FILE.001.01'//                          &
     &          MSGL//'INCLUDE Record Has Null File Name.'                     &
     &          )
        ELSE
          LBGN = INDEX(TEXT,FNEC) + 1
          LEND = LBGN + INDEX(TEXT(LBGN:)//FNEC,FNEC) - 2
          IF (LEND .LE. LBGN) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'VERIFY_INCLUDED_FILE.001.02'//                          &
     &          MSGL//'INCLUDE Record Has Null File Name.'                     &
     &          )
          ENDIF
        ENDIF
!!
!! Inquire about file's existence and obtain full file name.
!!
        IOERROR = .TRUE.
        INQUIRE                                                                &
     &          (                                                              &
     &          FILE        = TEXT(LBGN:LEND),                                 &
     &          EXIST       = EXSTAT,                                          &
     &          FORMATTED   = FORMATTED_READ,                                  &
     &          NAME        = FULL_NAME,                                       &
     &          ERR         = 200                                              &
     &          )
        IOERROR = .FALSE.
 200      INCLUDE_FILE(NIF)%Full_Name = TRIM(FULL_NAME)
!!
!! Warning for failed INQUIRE operation.
!!
        IF (IOERROR) THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'VERIFY_INCLUDED_FILE.002.01'//                          &
     &          MSGL//'Error Executing INQUIRE On: '//TEXT(LBGN:LEND)          &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
!!
!! Warning for non-existent file.
!!
        ELSE IF (.NOT.EXSTAT) THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'VERIFY_INCLUDED_FILE.002.02'//                          &
     &          MSGL//'INQUIRE Unable To Find: '//TEXT(LBGN:LEND)              &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
!!
!! Warning for file that isn't structured for formatted reads.
!!
        ELSE IF (INDEX(FORMATTED_READ,'NO') .NE. 0) THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'VERIFY_INCLUDED_FILE.002.03'//                          &
     &          MSGL//'Unable To Perform Formatted Read On: '                  &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
      ELSE
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'VERIFY_INCLUDED_FILE.003.00'//                          &
     &          MSGL//'INCLUDE Input Record Expected.'//                       &
     &          MSGL//'Input Record Key-Word Found: '//KEY_WORD//              &
     &          MSGL//'Logic Error. Call KEY Associates.'                      &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE PRE_READ_LAYUP_RECORDS
!!
!! Copyright (c) by KEY Associates, 25-JAN-1992 14:07:28
!!
!! Purpose: Pre-read the layering records for use with the layered solid,
!! LSOLD, data structures to count the extra hexahedron and 4-node membrane
!! elements needed to implement the layered solid elements.
!!
      USE shared_common_data
      USE value_
      USE layering_
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER :: TEXT*80, C_VALUE*32, KEY_WORD*12
      INTEGER   :: ILCOUNT   ! S/S Input file line counter
      INTEGER   :: ILINDEX   ! S/S Start-of-Record line in input file
      LOGICAL   :: EOF
      LOGICAL   :: IOERROR
      LOGICAL   :: FIRST_CALL = .TRUE.

      EXTERNAL  C_VALUE
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered PRE_READ_LAYUP_RECORDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Read those included files that contain LAYUP records.
!!
      ERROR%COUNT = 0
      DO NIF = 1,NUMIF
        IF (INCLUDE_FILE(NIF)%LAYUP_Begin .GT. 0) THEN
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(NIF)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200      CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PRE_READ_LAYUP_RECORDS.001.01'//                        &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of pre-read of LAYUP records.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'PRE_READ_LAYUP_RECORDS.001.02'//                        &
     &          MSGL//'INCLUDE File Pre-Read For LAYUP Records: '              &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
          ENDIF
!!
!! Position input file to first LAYUP record.
!!
          IERR = 0
          ILCOUNT = 0
          EOF = .FALSE.
          REWIND (IO_UNIT%LSDI)
          ILINDEX = INCLUDE_FILE(NIF)%LAYUP_Begin
          DO i = 1,ILINDEX - 1
            ILCOUNT = ILCOUNT + 1
            READ (IO_UNIT%LSDI,100)
 100        FORMAT ()
          ENDDO
!!
!! Read input records and only process the LAYUP records.
!!
          DO WHILE (.NOT.EOF.AND.ILCOUNT.LE.INCLUDE_FILE(NIF)%LAYUP_End)
            NVALS = MAXRE
            CALL GETIRV                                                        &
     & (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
            KEY_WORD = C_VALUE(1)
!!
            IF (INDEX(KEY_WORD,'LAYUP') .NE. 0) THEN
!!
!! Get lay-up specifications for layered solid.
!!
              CALL READ_LAYERED_SOLID_LAYUP (FIRST_CALL,NVALS)
              FIRST_CALL = .FALSE.
!!
            ELSE IF (INDEX(KEY_WORD,'ENDDATA') .NE. 0) THEN
              EOF = .TRUE.
            ENDIF
          ENDDO
!!
!! Close included file.
!!
          CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
        ENDIF
      ENDDO
!!
!! Print total number of open errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'PRE_READ_LAYUP_RECORDS.004.00'//                        &
     &          MSGL//'Total Number Of Open Errors In Input:'//MSG1//          &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE PRE_READ_LSOLD_RECORDS
!!
!! Copyright (c) by KEY Associates, 25-JAN-1992 14:07:28
!!
!! Purpose: Pre-read the layered solid, LSOLD, data structures to count the
!! extra hexahedron and 4-node membrane elements needed to implement the
!! layered solid elements.
!!
      USE shared_common_data
      USE value_
      USE lsold_
      USE layering_
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER :: TEXT*80, C_VALUE*32, KEY_WORD*12
      INTEGER   :: INTERNAL_ID, LY_INTERNAL_ID
      INTEGER   :: ILCOUNT   ! S/S Input file line counter
      INTEGER   :: ILINDEX   ! S/S Start-of-Record line in input file
      LOGICAL   :: EOF
      LOGICAL   :: IOERROR
      LOGICAL   :: FIRST_CALL = .TRUE.

      EXTERNAL C_VALUE, LY_INTERNAL_ID
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered PRE_READ_LSOLD_RECORDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Read those included files that contain LSEL records.
!!
      ERROR%COUNT = 0
      DO NIF = 1,NUMIF
        IF (INCLUDE_FILE(NIF)%LSOLD_Begin .GT. 0) THEN
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(NIF)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200      CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PRE_READ_LSOLD_RECORDS.001.01'//                        &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of pre-read of LSEL records.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'PRE_READ_LSOLD_RECORDS.001.02'//                        &
     &          MSGL//'INCLUDE File Pre-Read For LSEL Records: '               &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
          ENDIF
!!
!! Position input file to first LSEL record.
!!
          IERR = 0
          ILCOUNT = 0
          EOF = .FALSE.
          REWIND (IO_UNIT%LSDI)
          ILINDEX = INCLUDE_FILE(NIF)%LSOLD_Begin
          DO i = 1,ILINDEX - 1
            ILCOUNT = ILCOUNT + 1
            READ (IO_UNIT%LSDI,100)
 100        FORMAT ()
          ENDDO
!!
!! Read input records and only process the LSEL records.
!!
          DO WHILE (.NOT.EOF.AND.ILCOUNT.LE.INCLUDE_FILE(NIF)%LSOLD_End)
            NVALS = MAXRE
            CALL GETIRV                                                        &
     & (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
            KEY_WORD = C_VALUE(1)
!!
            IF (INDEX(KEY_WORD,'LSEL') .NE. 0) THEN
!!
!! Get layered solid element definition. Here, FIRST_CALL remains .TRUE.
!! As a consequence all pre-reads are read into LSOLD(1). We are doing
!! this in order to count the total number of hexah layers NUMLX and
!! membrane layers NUMLM that are needed for all of the layered sold
!! elements. During the regular read of the input data, the individual
!! LSOLD records will be read and saved as one would expect.
!!
              CALL READ_LSEL_DEFINITION (FIRST_CALL,NVALS)
!!
              INTERNAL_ID = LY_INTERNAL_ID (LSOLD(1)%PAR%LupID)
              IF (INTERNAL_ID .NE. 0) THEN
                DO i = 1,LAYERING(INTERNAL_ID)%Number_of_Layers
                  IF (LAYERING(INTERNAL_ID)%Ltype(i) .EQ. 0) THEN
                    NUMLX = NUMLX + 1
                  ELSE
                    NUMLM = NUMLM + 1
                  ENDIF
                ENDDO
              ELSE
!!
!! Matching lay-up ID not found.
!!
                WRITE (MSG1,'(I8)') LSOLD(1)%PAR%EleID
                WRITE (MSG2,'(I8)') LSOLD(1)%PAR%LupID
                CALL USER_MESSAGE                                              &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'PRE_READ_LSOLD_RECORDS.001.00'//                      &
     &            MSGL//'LSEL (Layered Solid) Input Record ID:'//MSG1//        &
     &            MSGL//'Contains Unknown LAYUP ID:'//MSG2                     &
     &            )
                ERROR%COUNT = ERROR%COUNT + 1
              ENDIF
!!
            ELSE IF (INDEX(KEY_WORD,'ENDDATA') .NE. 0) THEN
              EOF = .TRUE.
            ENDIF
          ENDDO
!!
!! Close included file.
!!
          CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
        ENDIF
      ENDDO
!!
!! Print total number of open & lay-up ID errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &  (                                                                      &
     &  MSGL//'FATAL'//                                                        &
     &  MSGL//'PRE_READ_LSOLD_RECORDS.004.00'//                                &
     &  MSGL//'Total Number Of Open Or ID Errors In Input:'//MSG1//            &
     &  MSGL//'Execution Terminated By Program.'                               &
     &  )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE PRE_READ_CONTROL_RECORDS
!!
!! Copyright (c) by KEY Associates, 14-FEB-1992 09:46:40
!!
!! Purpose: Pre-read the execution control records to locate any RD_MESH
!! switches set to "on."
!!
      USE shared_common_data
      USE value_
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          TEXT*80,                                                       &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      INTEGER   :: ILCOUNT   ! S/S Input file line counter
      INTEGER   :: ILINDEX   ! S/S Start-of-Record line in input file

      EXTERNAL                                                                 &
     &          C_VALUE
      LOGICAL                                                                  &
     &          EOF,                                                           &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered PRE_READ_CONTROL_RECORDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Read only those included files that contain CONTROL records.
!!
      ERROR%COUNT = 0
      DO NIF = 1,NUMIF
        IF (INCLUDE_FILE(NIF)%CONTROL_Begin .GT. 0) THEN
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(NIF)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200      CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'PRE_READ_CONTROL_RECORDS.001.01'//                      &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of pre-read of CONTROL records.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'PRE_READ_CONTROL_RECORDS.001.02'//                      &
     &          MSGL//'INCLUDE File Pre-Read For CONTROL Records: '            &
     &              //'fma3d.in'                                               &
     &          )
          ENDIF
!!
!! Position input file to first CONTROL record.
!!
          IERR = 0
          ILCOUNT = 0
          EOF = .FALSE.
          REWIND (IO_UNIT%LSDI)
          ILINDEX = INCLUDE_FILE(NIF)%CONTROL_Begin
          DO i = 1,ILINDEX - 1
            ILCOUNT = ILCOUNT + 1
            READ (IO_UNIT%LSDI,100)
 100        FORMAT ()
          ENDDO
!!
!! Read input records and only process the CONTROL records.
!!
          DO WHILE (.NOT.EOF .AND.                                             &
     &      ILCOUNT.LE.INCLUDE_FILE(NIF)%CONTROL_End)
            NVALS = MAXRE
            CALL GETIRV                                                        &
     & (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
            KEY_WORD = C_VALUE(1)
!!
            IF (INDEX(KEY_WORD,'CONTROL') .NE. 0) THEN
!!
!! Read control record (to see if there is a switch set to read a mesh
!! data input file).
!!
              CALL READ_CONTROL_VALUES (NVALS)
!!
            ELSE IF (INDEX(KEY_WORD,'ENDDATA') .NE. 0) THEN
              EOF = .TRUE.
            ENDIF
          ENDDO
!!
!! Close included file.
!!
          CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
        ENDIF
      ENDDO
!!
!! Print total number of open errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'PRE_READ_CONTROL_RECORDS.004.00'//                            &
     &    MSGL//'Total Number Of Open Errors In Input:'//MSG1//                &
     &    MSGL//'Execution Terminated By Program.'                             &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE SCAN_MESH_DATA
!!
!! Copyright (c) by KEY Associates, 12-FEB-1992 20:06:03
!!
!! Purpose: Read mesh data counters to add to counter values obtained
!! from scanning the ASCII input file(s).
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SCAN_MESH_DATA.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open mesh data input file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LMDI,                                        &
     &          FILE   = 'fmamdi',                                             &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                        &
     &          ERR    =  100                                                  &
     &          )
      IOERROR = .FALSE.
!!
!! Warning exit for failed OPEN operation.
!!
 100    IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SCAN_MESH_DATA.001.01'//                                &
     &          MSGL//'Unable To Execute OPEN On: '//'fmamdi'                  &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ELSE
!!
!! Part 1: Read title and program identification information.
!!
        READ (IO_UNIT%LMDI,*)                                                    &
     &          JOB_ID_RECORD%MESH%TITLE,                                      &
     &          JOB_ID_RECORD%MESH%DATEE,                                      &
     &          JOB_ID_RECORD%MESH%TIMEE
!!
!! Inform user of scan of mesh data file.
!!
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'INFORM'//                                                     &
     &    MSGL//'SCAN_MESH_DATA.001.02'//                                      &
     &    MSGL//'Mesh Data Title: '//JOB_ID_RECORD%MESH%TITLE//                &
     &    MSGL//'Mesh Data  Date: '//JOB_ID_RECORD%MESH%DATEE//                &
     &    MSGL//'Mesh Data  Time: '//JOB_ID_RECORD%MESH%TIMEE                  &
     &    )
!!
!! Part 2: Read counters from the mesh data file.
!!
        READ (IO_UNIT%LMDI,*)                                                    &
     &          MSHIF,MSHQA,MSHNP,MSHEL,MSHHX,MSHPX,MSHTX,MSHLS,MSHLX,         &
     &          MSHLM,MSHM4,MSHM3,MSHTR,MSHP4,MSHP3,MSHBM,MSHSP,MSHDM,         &
     &          MSHSG,MSHDC,MSHTC,MSHSW,MSHWC,MSHBF,MSHPC,MSHFC,MSHSC,         &
     &          MSHVC,MSHCC,MSHNR,MSHND,MSHFS,MSHIT,MSHSI,MSHSN,MSHCE,         &
     &          MSHCN,MSHCX,MSHNS,MSHNE,MSHES,MSHEE,MSHSS,MSHSE,MSHMT,         &
     &          MSHLU,MSHTF,MSHFP,MSHRF,MSHPV,MSHS1,MSHS2,MSHG1,MSHG2,         &
     &          MSHG3,MSHMP,MSHRM,MSHCM,MSHIC,MSHRB,MSHRT,MSHST,MSHAX,         &
     &          MSHPP,MSHNC,MSHRE,MSHID
!!
!! Add mesh data counters to counters obtained from ASCII input file. All
!! counters for which there is no mesh data are listed with "redundant"
!! equality statements.
!!
        NUMIF = NUMIF
        NUMQA = NUMQA + MSHQA
        NUMNP = NUMNP + MSHNP
        NUMEL = NUMEL + MSHEL
        NUMHX = NUMHX + MSHHX
        NUMPX = NUMPX + MSHPX
        NUMTX = NUMTX + MSHTX
        NUMLS = NUMLS + MSHLS
        NUMLX = NUMLX + MSHLX
        NUMLM = NUMLM + MSHLM
        NUMM4 = NUMM4 + MSHM4
        NUMM3 = NUMM3 + MSHM3
        NUMTR = NUMTR + MSHTR
        NUMP4 = NUMP4 + MSHP4
        NUMP3 = NUMP3 + MSHP3
        NUMBM = NUMBM + MSHBM
        NUMSP = NUMSP + MSHSP
        NUMDM = NUMDM + MSHDM
        NUMSG = NUMSG + MSHSG
        NUMDC = NUMDC + MSHDC
        NUMTC = NUMTC + MSHTC
        NUMSW = NUMSW + MSHSW
        NUMWC = NUMWC + MSHWC
        NUMBF = NUMBF + MSHBF
        NUMPC = NUMPC + MSHPC
        NUMFC = NUMFC + MSHFC
        NUMSC = NUMSC + MSHSC
        NUMVC = NUMVC + MSHVC
        NUMCC = NUMCC + MSHCC
        NUMNR = NUMNR + MSHNR
        NUMND = NUMND
        NUMFS = NUMFS + MSHFS
        NUMIT = NUMIT
        NUMSI = NUMSI + MSHSI
        NUMSN = NUMSN
        NUMCE = NUMCE
        NUMCN = NUMCN
        NUMCX = NUMCX
        NUMNS = NUMNS + MSHNS
        NUMNE = NUMNE + MSHNE
        NUMES = NUMES + MSHES
        NUMEE = NUMEE + MSHEE
        NUMSS = NUMSS + MSHSS
        NUMSE = NUMSE + MSHSE
        NUMMT = NUMMT + MSHMT
        NUMLU = NUMLU + MSHLU
        NUMTF = NUMTF + MSHTF
        NUMFP = MAX (NUMFP,MSHFP)
        NUMRF = NUMRF
        NUMPV = NUMPV
        NUMS1 = NUMS1 + MSHS1
        NUMS2 = NUMS2 + MSHS2
        NUMG1 = NUMG1 + MSHG1
        NUMG2 = NUMG2 + MSHG2
        NUMG3 = NUMG3 + MSHG3
        NUMMP = NUMMP
        NUMRM = NUMRM + MSHRM
        NUMCM = NUMCM + MSHCM
        NUMIC = NUMIC + MSHIC
        NUMRB = NUMRB + MSHRB
        NUMST = NUMST
        NUMAX = NUMAX
        NUMRT = MIN (2*NUMNP,NUMRT+MSHRT)
        NUMPP = NUMPP
        NUMNC = NUMNC + MSHNC
        MAXRE = MAXRE
        MXEID = MAX (MXEID,MSHID)
!!
        CLOSE (UNIT=IO_UNIT%LMDI,STATUS='KEEP')
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_MESH_DATA
!!
!! Copyright (c) by KEY Associates, 12-FEB-1992 19:57:22
!!
!! Purpose: Read a mesh data file that will be combined with the material
!! from the standard ASCII input file.
!!
      USE shared_common_data;
      USE indx_;           USE node_;           USE tabulated_function_;
      USE beam_;           USE coord_;          USE sliding_interface_;
      USE value_;          USE force_;          USE constrained_node_;
      USE hexah_;          USE penta_;          USE nonreflecting_bc_;
      USE tetra_;          USE lsold_;          USE nodal_point_mass_;
      USE membq_;          USE membt_;          USE rigid_body_mass_;
      USE truss_;          USE platq_;          USE state_variables_;
      USE platt_;          USE motion_;         USE enumerated_sets_;
      USE spring_;         USE damper_;         USE displacement_bc_;
      USE stress_;         USE segment_;        USE contact_surface_;
      USE tied_bc_;        USE results_;        USE relink_scratch_;
      USE gauge1d_;        USE gauge2d_;        USE rigid_wall_bc_;
      USE gauge3d_;        USE massprop_;       USE include_file_;
      USE material_;       USE layering_;       USE sliding_node_;
      USE force_bc_;       USE node_set_;       USE contact_node_;
      USE nrbc_data_;      USE spring_bc_;      USE periodic_bc_;
      USE damper_bc_;      USE spot_weld_;      USE pressure_bc_;
      USE qa_record_;      USE plate_pair_;     USE segment_set_;
      USE body_force_;     USE section_2d_;     USE element_set_;
      USE section_1d_;     USE rigid_body_;     USE plate_pair_;
      USE velocity_ic_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered READ_MESH_DATA'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open mesh data input file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LMDI,                                        &
     &          FILE   = 'fmamdi',                                             &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                        &
     &          ERR    =  100                                                  &
     &          )
      IOERROR = .FALSE.
!!
!! Warning exit for failed OPEN operation.
!!
 100    IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_MESH_DATA.001.01'//                                &
     &          MSGL//'Unable To Execute OPEN On: '//'fmamdi'                  &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ELSE
!!
!! Inform user of read of mesh data file.
!!
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_MESH_DATA.001.02'//                                &
     &          MSGL//'Mesh Data Title: '//JOB_ID_RECORD%MESH%TITLE//          &
     &          MSGL//'Mesh Data  Date: '//JOB_ID_RECORD%MESH%DATEE//          &
     &          MSGL//'Mesh Data  Time: '//JOB_ID_RECORD%MESH%TIMEE            &
     &          )
!!
!! Part 1: Skip over title and program identification information.
!!
        READ (IO_UNIT%LMDI,*)
!!
!! Part 2: Skip over mesh data counters.
!!
        READ (IO_UNIT%LMDI,*)
!!
!! Part 3: Read data structures contained in mesh data input.
!!
        DO i = NUMQA-MSHQA+1,NUMQA
          READ (IO_UNIT%LMDI,*) QA_RECORD(i)
        ENDDO
        DO i = NUMG1-MSHG1+1,NUMG1
          READ (IO_UNIT%LMDI,*) GAUGE1D(i)%PAR
        ENDDO
        DO i = NUMG2-MSHG2+1,NUMG2
          READ (IO_UNIT%LMDI,*) GAUGE2D(i)%PAR
        ENDDO
        DO i = NUMG3-MSHG3+1,NUMG3
          READ (IO_UNIT%LMDI,*) GAUGE3D(i)%PAR
        ENDDO
        DO i = NUMMT-MSHMT+1,NUMMT
          READ (IO_UNIT%LMDI,*) MATERIAL(i)
        ENDDO
        DO i = NUMLU-MSHLU+1,NUMLU
          READ (IO_UNIT%LMDI,*) LAYERING(i)
        ENDDO
        DO i = NUMS2-MSHS2+1,NUMS2
          READ (IO_UNIT%LMDI,*) SECTION_2D(i)
        ENDDO
        DO i = NUMS1-MSHS1+1,NUMS1
          READ (IO_UNIT%LMDI,*) SECTION_1D(i)
        ENDDO
        DO i = NUMRB-MSHRB+1,NUMRB
          READ (IO_UNIT%LMDI,*) RIGID_BODY(i)
        ENDDO
        DO i = NUMRM-MSHRM+1,NUMRM
          READ (IO_UNIT%LMDI,*) RIGID_BODY_MASS(i)
        ENDDO
        DO i = NUMCM-MSHCM+1,NUMCM
          READ (IO_UNIT%LMDI,*) NODAL_POINT_MASS(i)
        ENDDO
        DO i = NUMDC-MSHDC+1,NUMDC
          READ (IO_UNIT%LMDI,*) DISPLACEMENT_BC(i)
        ENDDO
        DO i = NUMTC-MSHTC+1,NUMTC
          READ (IO_UNIT%LMDI,*) TIED_BC(i)
        ENDDO
        DO i = NUMSW-MSHSW+1,NUMSW
          READ (IO_UNIT%LMDI,*) SPOT_WELD(i)
        ENDDO
        DO i = NUMWC-MSHWC+1,NUMWC
          READ (IO_UNIT%LMDI,*) RIGID_WALL_BC(i)
        ENDDO
        DO i = NUMBF-MSHBF+1,NUMBF
          READ (IO_UNIT%LMDI,*) BODY_FORCE(i)
        ENDDO
        DO i = NUMPC-MSHPC+1,NUMPC
          READ (IO_UNIT%LMDI,*) PRESSURE_BC(i)
        ENDDO
        DO i = NUMFC-MSHFC+1,NUMFC
          READ (IO_UNIT%LMDI,*) FORCE_BC(i)
        ENDDO
        DO i = NUMSC-MSHSC+1,NUMSC
          READ (IO_UNIT%LMDI,*) SPRING_BC(i)
        ENDDO
        DO i = NUMVC-MSHVC+1,NUMVC
          READ (IO_UNIT%LMDI,*) DAMPER_BC(i)
        ENDDO
        DO i = NUMCC-MSHCC+1,NUMCC
          READ (IO_UNIT%LMDI,*) PERIODIC_BC(i)
        ENDDO
        DO i = NUMNR-MSHNR+1,NUMNR
          READ (IO_UNIT%LMDI,*) NONREFLECTING_BC(i)
        ENDDO
        DO i = NUMSI-MSHSI+1,NUMSI
          READ (IO_UNIT%LMDI,*) SLIDING_INTERFACE(i)
        ENDDO
        DO i = NUMTF-MSHTF+1,NUMTF
          READ (IO_UNIT%LMDI,*) TABULATED_FUNCTION(i)
        ENDDO
        DO i = NUMNS-MSHNS+1,NUMNS
          READ (IO_UNIT%LMDI,*) NODE_SET(i)
        ENDDO
        DO i = NUMES-MSHES+1,NUMES
          READ (IO_UNIT%LMDI,*) ELEMENT_SET(i)
        ENDDO
        DO i = NUMSS-MSHSS+1,NUMSS
          READ (IO_UNIT%LMDI,*) SEGMENT_SET(i)
        ENDDO
        DO i = NUMSG-MSHSG+1,NUMSG
          READ (IO_UNIT%LMDI,*) SEGMENT(i)%PAR
        ENDDO
        DO i = NUMNP-MSHNP+1,NUMNP
          READ (IO_UNIT%LMDI,*) NODE(i)%ID,MOTION(i)%Px,                         &
     &                      MOTION(i)%Py,MOTION(i)%Pz
        ENDDO
        DO i = NUMHX-MSHHX+1,NUMHX
          READ (IO_UNIT%LMDI,*) HEXAH(i)%PAR
        ENDDO
        DO i = NUMPX-MSHPX+1,NUMPX
          READ (IO_UNIT%LMDI,*) PENTA(i)%PAR
        ENDDO
        DO i = NUMTX-MSHTX+1,NUMTX
          READ (IO_UNIT%LMDI,*) TETRA(i)%PAR
        ENDDO
        DO i = NUMLS-MSHLS+1,NUMLS
          READ (IO_UNIT%LMDI,*) LSOLD(i)%PAR
        ENDDO
        DO i = NUMLX-MSHLX+1,NUMLX
          READ (IO_UNIT%LMDI,*) LSHEX(i)%PAR
        ENDDO
        DO i = NUMLM-MSHLM+1,NUMLM
          READ (IO_UNIT%LMDI,*) LSMBQ(i)%PAR
        ENDDO
        DO i = NUMM3-MSHM3+1,NUMM3
          READ (IO_UNIT%LMDI,*) MEMBT(i)%PAR
        ENDDO
        DO i = NUMM4-MSHM4+1,NUMM4
          READ (IO_UNIT%LMDI,*) MEMBQ(i)%PAR
        ENDDO
        DO i = NUMTR-MSHTR+1,NUMTR
          READ (IO_UNIT%LMDI,*) TRUSS(i)%PAR
        ENDDO
        DO i = NUMP3-MSHP3+1,NUMP3
          READ (IO_UNIT%LMDI,*) PLATT(i)%PAR
        ENDDO
        DO i = NUMP4-MSHP4+1,NUMP4
          READ (IO_UNIT%LMDI,*) PLATQ(i)%PAR
        ENDDO
        DO i = NUMBM-MSHBM+1,NUMBM
          READ (IO_UNIT%LMDI,*) BEAM(i)%PAR
        ENDDO
        DO i = NUMSP-MSHSP+1,NUMSP
          READ (IO_UNIT%LMDI,*) SPRING(i)%PAR
        ENDDO
        DO i = NUMDM-MSHDM+1,NUMDM
          READ (IO_UNIT%LMDI,*) DAMPER(i)%PAR
        ENDDO
        DO i = NUMIC-MSHIC+1,NUMIC
          READ (IO_UNIT%LMDI,*) VELOCITY_IC(i)
        ENDDO
        DO i = NUMNC-MSHNC+1,NUMNC
          READ (IO_UNIT%LMDI,*) CONSTRAINED_NODE(i)
        ENDDO
!!
!! Part 4: The following reads are deferred until the main program is ready to
!! build the nodal point, element and segment set member lists. The mesh data
!! file is closed after building nodal point, element and segment member lists.
!!
!!        DO i = NUMNE-MSHNE+1,NUMNE
!!          READ (IO_UNIT%LMDI,*) NNPSETS(i)
!!        ENDDO
!!        DO i = NUMEE-MSHEE+1,NUMEE
!!          READ (IO_UNIT%LMDI,*) NELSETS(i)
!!        ENDDO
!!        DO i = NUMSE-MSHSE+1,NUMSE
!!          READ (IO_UNIT%LMDI,*) NSGSETS(i)
!!        ENDDO
!!
!!        CLOSE (UNIT=IO_UNIT%LMDI, STATUS='KEEP')
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE WRITE_MESH_DATA
!!
!! Copyright (c) by KEY Associates, 10-FEB-1992 19:12:44
!!
!! Purpose: Write a mesh data file that combines the material from the
!! standard ASCII input file and any other previously read mesh data.
!!
      USE shared_common_data;
      USE indx_;           USE node_;           USE tabulated_function_;
      USE beam_;           USE coord_;          USE sliding_interface_;
      USE value_;          USE force_;          USE constrained_node_;
      USE hexah_;          USE penta_;          USE nonreflecting_bc_;
      USE tetra_;          USE lsold_;          USE nodal_point_mass_;
      USE membq_;          USE membt_;          USE rigid_body_mass_;
      USE truss_;          USE platq_;          USE state_variables_;
      USE platt_;          USE motion_;         USE enumerated_sets_;
      USE spring_;         USE damper_;         USE displacement_bc_;
      USE stress_;         USE segment_;        USE contact_surface_;
      USE tied_bc_;        USE results_;        USE relink_scratch_;
      USE gauge1d_;        USE gauge2d_;        USE rigid_wall_bc_;
      USE gauge3d_;        USE massprop_;       USE include_file_;
      USE material_;       USE layering_;       USE sliding_node_;
      USE force_bc_;       USE node_set_;       USE contact_node_;
      USE nrbc_data_;      USE spring_bc_;      USE periodic_bc_;
      USE damper_bc_;      USE spot_weld_;      USE pressure_bc_;
      USE qa_record_;      USE plate_pair_;     USE segment_set_;
      USE body_force_;     USE section_2d_;     USE element_set_;
      USE section_1d_;     USE rigid_body_;     USE plate_pair_;
      USE velocity_ic_;
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered WRITE_MESH_DATA.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open mesh data output file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LMDO,                                        &
     &          FILE   = 'fmamdo',                                             &
     &          STATUS = 'NEW',                                                &
     &          FORM   = 'FORMATTED',                                        &
     &          ERR    =  100                                                  &
     &          )
      IOERROR = .FALSE.
!!
!! Warning exit for failed OPEN operation.
!!
 100    IF (IOERROR) THEN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'WRITE_MESH_DATA.001.01'//                               &
     &          MSGL//'Unable To Execute OPEN On: '//'fmamdo'                  &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ELSE
!!
!! Inform user of write of mesh data file.
!!
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'WRITE_MESH_DATA.001.02'//                               &
     &          MSGL//'Mesh Data Title: '//JOB_ID_RECORD%CURRENT%TITLE//       &
     &          MSGL//'Mesh Data  Date: '//JOB_ID_RECORD%CURRENT%DATEE//       &
     &          MSGL//'Mesh Data  Time: '//JOB_ID_RECORD%CURRENT%TIMEE         &
     &          )
!!
!! Part 1: Write title and program identification information.
!!
        WRITE (IO_UNIT%LMDO,*)                                                   &
     &          JOB_ID_RECORD%CURRENT%TITLE,                                   &
     &          JOB_ID_RECORD%CURRENT%DATEE,                                   &
     &          JOB_ID_RECORD%CURRENT%TIMEE,                                   &
     &          JOB_ID_RECORD%PROGRAM,                                         &
     &          JOB_ID_RECORD%VERSION
!!
!! Part 2: Write counters to the mesh data file.
!!
        WRITE (IO_UNIT%LMDO,*)                                                   &
     &          NUMIF,NUMQA,NUMNP,NUMEL,NUMHX,NUMPX,NUMTX,NUMLS,NUMLX,         &
     &          NUMLM,NUMM4,NUMM3,NUMTR,NUMP4,NUMP3,NUMBM,NUMSP,NUMDM,         &
     &          NUMSG,NUMDC,NUMTC,NUMSW,NUMWC,NUMBF,NUMPC,NUMFC,NUMSC,         &
     &          NUMVC,NUMCC,NUMNR,NUMND,NUMFS,NUMIT,NUMSI,NUMSN,NUMCE,         &
     &          NUMCN,NUMCX,NUMNS,NUMNE,NUMES,NUMEE,NUMSS,NUMSE,NUMMT,         &
     &          NUMLU,NUMTF,NUMFP,NUMRF,NUMPV,NUMS1,NUMS2,NUMG1,NUMG2,         &
     &          NUMG3,NUMMP,NUMRM,NUMCM,NUMIC,NUMRB,NUMRT,NUMST,NUMAX,         &
     &          NUMPP,NUMNC,MAXRE,MXEID
!!
!! Part 3: Write all avaiable data structures.
!!
        DO i = 1,NUMQA
          WRITE (IO_UNIT%LMDO,*) QA_RECORD(i)
        ENDDO
        DO i = 1,NUMG1
          WRITE (IO_UNIT%LMDO,*) GAUGE1D(i)%PAR
        ENDDO
        DO i = 1,NUMG2
          WRITE (IO_UNIT%LMDO,*) GAUGE2D(i)%PAR
        ENDDO
        DO i = 1,NUMG3
          WRITE (IO_UNIT%LMDO,*) GAUGE3D(i)%PAR
        ENDDO
        DO i = 1,NUMMT
          WRITE (IO_UNIT%LMDO,*) MATERIAL(i)
        ENDDO
        DO i = 1,NUMLU
          WRITE (IO_UNIT%LMDO,*) LAYERING(i)
        ENDDO
        DO i = 1,NUMS2
          WRITE (IO_UNIT%LMDO,*) SECTION_2D(i)
        ENDDO
        DO i = 1,NUMS1
          WRITE (IO_UNIT%LMDO,*) SECTION_1D(i)
        ENDDO
        DO i = 1,NUMRB
          WRITE (IO_UNIT%LMDO,*) RIGID_BODY(i)
        ENDDO
        DO i = 1,NUMRM
          WRITE (IO_UNIT%LMDO,*) RIGID_BODY_MASS(i)
        ENDDO
        DO i = 1,NUMCM
          WRITE (IO_UNIT%LMDO,*) NODAL_POINT_MASS(i)
        ENDDO
        DO i = 1,NUMDC
          WRITE (IO_UNIT%LMDO,*) DISPLACEMENT_BC(i)
        ENDDO
        DO i = 1,NUMTC
          WRITE (IO_UNIT%LMDO,*) TIED_BC(i)
        ENDDO
        DO i = 1,NUMSW
          WRITE (IO_UNIT%LMDO,*) SPOT_WELD(i)
        ENDDO
        DO i = 1,NUMWC
          WRITE (IO_UNIT%LMDO,*) RIGID_WALL_BC(i)
        ENDDO
        DO i = 1,NUMBF
          WRITE (IO_UNIT%LMDO,*) BODY_FORCE(i)
        ENDDO
        DO i = 1,NUMPC
          WRITE (IO_UNIT%LMDO,*) PRESSURE_BC(i)
        ENDDO
        DO i = 1,NUMFC
          WRITE (IO_UNIT%LMDO,*) FORCE_BC(i)
        ENDDO
        DO i = 1,NUMSC
          WRITE (IO_UNIT%LMDO,*) SPRING_BC(i)
        ENDDO
        DO i = 1,NUMVC
          WRITE (IO_UNIT%LMDO,*) DAMPER_BC(i)
        ENDDO
        DO i = 1,NUMCC
          WRITE (IO_UNIT%LMDO,*) PERIODIC_BC(i)
        ENDDO
        DO i = 1,NUMNR
          WRITE (IO_UNIT%LMDO,*) NONREFLECTING_BC(i)
        ENDDO
        DO i = 1,NUMSI
          WRITE (IO_UNIT%LMDO,*) SLIDING_INTERFACE(i)
        ENDDO
        DO i = 1,NUMTF
          WRITE (IO_UNIT%LMDO,*) TABULATED_FUNCTION(i)
        ENDDO
        DO i = 1,NUMNS
          WRITE (IO_UNIT%LMDO,*) NODE_SET(i)
        ENDDO
        DO i = 1,NUMES
          WRITE (IO_UNIT%LMDO,*) ELEMENT_SET(i)
        ENDDO
        DO i = 1,NUMSS
          WRITE (IO_UNIT%LMDO,*) SEGMENT_SET(i)
        ENDDO
        DO i = 1,NUMSG
          WRITE (IO_UNIT%LMDO,*) SEGMENT(i)%PAR
        ENDDO
        DO i = 1,NUMNP
          WRITE (IO_UNIT%LMDO,*) NODE(i)%ID,MOTION(i)%Px,                        &
     &                       MOTION(i)%Py,MOTION(i)%Pz
        ENDDO
        DO i = 1,NUMHX
          WRITE (IO_UNIT%LMDO,*) HEXAH(i)%PAR
        ENDDO
        DO i = 1,NUMPX
          WRITE (IO_UNIT%LMDO,*) PENTA(i)%PAR
        ENDDO
        DO i = 1,NUMTX
          WRITE (IO_UNIT%LMDO,*) TETRA(i)%PAR
        ENDDO
        DO i = 1,NUMLS
          WRITE (IO_UNIT%LMDO,*) LSOLD(i)%PAR
        ENDDO
        DO i = 1,NUMLX
          WRITE (IO_UNIT%LMDO,*) LSHEX(i)%PAR
        ENDDO
        DO i = 1,NUMLM
          WRITE (IO_UNIT%LMDO,*) LSMBQ(i)%PAR
        ENDDO
        DO i = 1,NUMM3
          WRITE (IO_UNIT%LMDO,*) MEMBT(i)%PAR
        ENDDO
        DO i = 1,NUMM4
          WRITE (IO_UNIT%LMDO,*) MEMBQ(i)%PAR
        ENDDO
        DO i = 1,NUMTR
          WRITE (IO_UNIT%LMDO,*) TRUSS(i)%PAR
        ENDDO
        DO i = 1,NUMP3
          WRITE (IO_UNIT%LMDO,*) PLATT(i)%PAR
        ENDDO
        DO i = 1,NUMP4
          WRITE (IO_UNIT%LMDO,*) PLATQ(i)%PAR
        ENDDO
        DO i = 1,NUMBM
          WRITE (IO_UNIT%LMDO,*) BEAM(i)%PAR
        ENDDO
        DO i = 1,NUMSP
          WRITE (IO_UNIT%LMDO,*) SPRING(i)%PAR
        ENDDO
        DO i = 1,NUMDM
          WRITE (IO_UNIT%LMDO,*) DAMPER(i)%PAR
        ENDDO
        DO i = 1,NUMIC
          WRITE (IO_UNIT%LMDO,*) VELOCITY_IC(i)
        ENDDO
        DO i = 1,NUMNC
          WRITE (IO_UNIT%LMDO,*) CONSTRAINED_NODE(i)
        ENDDO
!!
!! Part 4: Write the nodal point, element and segment set member lists.
!!
        DO i = 1,NUMNE
          WRITE (IO_UNIT%LMDO,*) NNPSETS(i)
        ENDDO
        DO i = 1,NUMEE
          WRITE (IO_UNIT%LMDO,*) NELSETS(i)
        ENDDO
        DO i = 1,NUMSE
          WRITE (IO_UNIT%LMDO,*) NSGSETS(i)
        ENDDO
!!
        CLOSE (UNIT=IO_UNIT%LMDO, STATUS='KEEP')
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE SCAN_RESTART_DATA
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 10:11:16
!!
!! Purpose: Read restart data counters to add to counter values obtained
!! from scanning the ASCII input file(s).
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SCAN_RESTART_DATA.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open restart data input file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LRDI,                                        &
     &          FILE   = 'fmardi',                                             &
     &          STATUS = 'OLD',                                                &
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
     &          MSGL//'SCAN_RESTART_DATA.001.01'//                             &
     &          MSGL//'Unable To Execute OPEN On: '//'fmardi'                  &
     &          )
      ELSE
!!
!! Read title and program identification information.
!!
        READ (IO_UNIT%LRDI)                                                    &
     &          JOB_ID_RECORD%RESTART%TITLE,                                   &
     &          JOB_ID_RECORD%RESTART%DATEE,                                   &
     &          JOB_ID_RECORD%RESTART%TIMEE,                                   &
     &          JOB_ID_RECORD%RESTART%SEQUENCE_NUMBER
!!
        JOB_ID_RECORD%CURRENT%SEQUENCE_NUMBER =                                &
     &    JOB_ID_RECORD%RESTART%SEQUENCE_NUMBER
!!
!! Inform user of scan of restart data file.
!!
        WRITE (MSG1,'(I8)') JOB_ID_RECORD%RESTART%SEQUENCE_NUMBER
        MSG1(1:5) = '.rdo.'
        CALL USER_MESSAGE                                                      &
     &     (                                                                   &
     &     MSGL//'INFORM'//                                                    &
     &     MSGL//'SCAN_RESTART_DATA.001.02'//                                  &
     &     MSGL//'Restart Data Title: '//JOB_ID_RECORD%RESTART%TITLE//         &
     &     MSGL//'Restart Data  Date: '//JOB_ID_RECORD%RESTART%DATEE//         &
     &     MSGL//'Restart Data  Time: '//JOB_ID_RECORD%RESTART%TIMEE//         &
     &     MSGL//'Restart Sequence #: '//MSG1                                  &
     &     )
!!
!! Skip over time and energy balance data.
!!
        READ (IO_UNIT%LRDI)
!!
!! Read counters from old restart data file.
!!
        READ (IO_UNIT%LRDI)                                                    &
     &          NRSIF,NRSQA,NRSNP,NRSEL,NRSHX,NRSPX,NRSTX,NRSLS,NRSLX,         &
     &          NRSLM,NRSM4,NRSM3,NRSTR,NRSP4,NRSP3,NRSBM,NRSSP,NRSDM,         &
     &          NRSSG,NRSDC,NRSTC,NRSSW,NRSWC,NRSBF,NRSPC,NRSFC,NRSSC,         &
     &          NRSVC,NRSCC,NRSNR,NRSND,NRSFS,NRSIT,NRSSI,NRSSN,NRSCE,         &
     &          NRSCN,NRSCX,NRSNS,NRSNE,NRSES,NRSEE,NRSSS,NRSSE,NRSMT,         &
     &          NRSLU,NRSTF,NRSFP,NRSRF,NRSPV,NRSS1,NRSS2,NRSG1,NRSG2,         &
     &          NRSG3,NRSMP,NRSRM,NRSCM,NRSIC,NRSRB,NRSRT,NRSST,NRSAX,         &
     &          NRSPP,NRSNC,NRSRE,NRSID
!!
        NUMIF = NUMIF         ! INCLUDE, new included files counter.
        NUMQA = NUMQA + NRSQA ! QAREC, added QA records.
        NUMNP = NRSNP
        NUMEL = NRSEL
        NUMHX = NRSHX
        NUMPX = NRSPX
        NUMTX = NRSTX
        NUMLS = NRSLS
        NUMLX = NRSLX
        NUMLM = NRSLM
        NUMM4 = NRSM4
        NUMM3 = NRSM3
        NUMTR = NRSTR
        NUMP4 = NRSP4
        NUMP3 = NRSP3
        NUMBM = NRSBM
        NUMSP = NRSSP
        NUMDM = NRSDM
        NUMSG = NRSSG
        NUMDC = NRSDC
        NUMTC = NRSTC
        NUMSW = NRSSW
        NUMWC = NRSWC
        NUMBF = NRSBF
        NUMPC = NRSPC
        NUMFC = NRSFC
        NUMSC = NRSSC
        NUMVC = NRSVC
        NUMCC = NRSCC
        NUMNR = NRSNR
        NUMND = NRSND
        NUMFS = NRSFS
        NUMIT = NUMIT         ! INTERFACE, new sliding interface times
        NUMSI = NRSSI
        NUMSN = NRSSN
        NUMCE = NRSCE
        NUMCN = NRSCN
        NUMCX = NRSCX
        NUMNS = NRSNS
        NUMNE = NRSNE
        NUMES = NRSES
        NUMEE = NRSEE
        NUMSS = NRSSS
        NUMSE = NRSSE
        NUMMT = NRSMT
        NUMLU = NRSLU
        NUMTF = NRSTF
        NUMFP = NRSFP
        NUMRF = NUMRF         ! RESULTS, new results output requests.
        NUMPV = NUMPV         ! PARAM, new parameter values
        NUMS1 = NRSS1
        NUMS2 = NRSS2
        NUMG1 = NRSG1
        NUMG2 = NRSG2
        NUMG3 = NRSG3
        NUMMP = NUMMP         ! MASSPROP, new mass property output requests.
        NUMRM = NRSRM
        NUMCM = NRSCM
        NUMIC = NRSIC
        NUMRB = NRSRB
        NUMST = NRSST
        NUMAX = NRSAX
        NUMRT = NRSRT
        NUMPP = NRSPP
        NUMNC = NRSNC
        MXEID = NRSID
!!
        CLOSE (UNIT=IO_UNIT%LRDI, STATUS='KEEP')
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_RESTART_DATA
!!
!! Copyright (c) by KEY Associates, 19-FEB-1992 19:33:13
!!
!! Purpose: READ all of the data needed to restart the analysis.
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
      USE nrbc_data_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE plate_pair_
!!
      USE enumerated_sets_
      USE stress_
      USE state_variables_
      USE coord_
      USE indx_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      LOGICAL                                                                  &
     &          IOERROR
      CHARACTER                                                                &
     &          RECORD_BEING_READ*32    ! Text for read error message.
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered READ_RESTART_DATA.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Open restart data input file.
!!
      IOERROR = .TRUE.
      OPEN                                                                     &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LRDI,                                        &
     &          FILE   = 'fmardi',                                             &
     &          STATUS = 'OLD',                                                &
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
     &          MSGL//'READ_RESTART_DATA.001.01'//                             &
     &          MSGL//'Unable To Execute OPEN On: '//'fmardi'                  &
     &          )
      ELSE
!!
!! Skip over title and program identification information.
!!
        IRECORD = 0
        RECORD_BEING_READ = 'Title, Program ID'
        READ (IO_UNIT%LRDI,ERR=999)
!!
!! Read time and energy balance data, however, save TIMSIM%Stop in case it has
!! been reset by a STOP input record.
!!
        STOP_TIME = TIMSIM%Stop
        RECORD_BEING_READ = 'Common Blocks'
        READ (IO_UNIT%LRDI,ERR=999) TIME,ENERGY,COUNTER,REZONE,                &
     &    SUBPROCESS_TIME
        TIMSIM%Stop = STOP_TIME
!!
!! Skip over counters in restart data file.
!!
        RECORD_BEING_READ = 'Counters'
        READ (IO_UNIT%LRDI,ERR=999)
!!
!! Read all data structures needed to restart simulation.
!!
        RECORD_BEING_READ = 'QA_RECORD'
        DO i = NUMQA-NRSQA+1,NUMQA
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) QA_RECORD(i)
        ENDDO
        RECORD_BEING_READ = 'GAUGE1D'
        DO i = 1,NUMG1
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) GAUGE1D(i)
        ENDDO
        RECORD_BEING_READ = 'GAUGE2D'
        DO i = 1,NUMG2
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) GAUGE2D(i)
        ENDDO
        RECORD_BEING_READ = 'GAUGE3D'
        DO i = 1,NUMG3
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) GAUGE3D(i)
        ENDDO
        RECORD_BEING_READ = 'MATERIAL'
        DO i = 1,NUMMT
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) MATERIAL(i)
        ENDDO
        RECORD_BEING_READ = 'LAYERING'
        DO i = 1,NUMLU
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) LAYERING(i)
        ENDDO
        RECORD_BEING_READ = 'SECTION_2D'
        DO i = 1,NUMS2
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SECTION_2D(i)
        ENDDO
        RECORD_BEING_READ = 'SECTION_1D'
        DO i = 1,NUMS1
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SECTION_1D(i)
        ENDDO
        RECORD_BEING_READ = 'RIGID_BODY'
        DO i = 1,NUMRB
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) RIGID_BODY(i)
        ENDDO
        RECORD_BEING_READ = 'RIGID_BODY_MASS'
        DO i = 1,NUMRM
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) RIGID_BODY_MASS(i)
        ENDDO
        RECORD_BEING_READ = 'NODAL_POINT_MASS'
        DO i = 1,NUMCM
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NODAL_POINT_MASS(i)
        ENDDO
        RECORD_BEING_READ = 'DISPLACEMENT_BC'
        DO i = 1,NUMDC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) DISPLACEMENT_BC(i)
        ENDDO
        RECORD_BEING_READ = 'TIED_BC'
        DO i = 1,NUMTC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) TIED_BC(i)
        ENDDO
        RECORD_BEING_READ = 'SPOT_WELD'
        DO i = 1,NUMSW
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SPOT_WELD(i)
        ENDDO
        RECORD_BEING_READ = 'RIGID_WALL'
        DO i = 1,NUMWC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) RIGID_WALL_BC(i)
        ENDDO
        RECORD_BEING_READ = 'BODY_FORCE'
        DO i = 1,NUMBF
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) BODY_FORCE(i)
        ENDDO
        RECORD_BEING_READ = 'PRESSURE_BC'
        DO i = 1,NUMPC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PRESSURE_BC(i)
        ENDDO
        RECORD_BEING_READ = 'FORCE_BC'
        DO i = 1,NUMFC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) FORCE_BC(i)
        ENDDO
        RECORD_BEING_READ = 'SPRING_BC'
        DO i = 1,NUMSC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SPRING_BC(i)
        ENDDO
        RECORD_BEING_READ = 'DAMPER_BC'
        DO i = 1,NUMVC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) DAMPER_BC(i)
        ENDDO
        RECORD_BEING_READ = 'PERIODIC_BC'
        DO i = 1,NUMCC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PERIODIC_BC(i)
        ENDDO
        RECORD_BEING_READ = 'NONREFLECTING_BC'
        DO i = 1,NUMNR
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NONREFLECTING_BC(i)
        ENDDO
        RECORD_BEING_READ = 'SLIDING_INTERFACE'
        DO i = 1,NUMSI
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SLIDING_INTERFACE(i)
        ENDDO
        RECORD_BEING_READ = 'TABULATED_FUNCTION'
        DO i = 1,NUMTF
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) TABULATED_FUNCTION(i)
        ENDDO
        RECORD_BEING_READ = 'NODE_SET'
        DO i = 1,NUMNS
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NODE_SET(i)
        ENDDO
        RECORD_BEING_READ = 'ELEMENT_SET'
        DO i = 1,NUMES
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) ELEMENT_SET(i)
        ENDDO
        RECORD_BEING_READ = 'SEGMENT_SET'
        DO i = 1,NUMSS
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SEGMENT_SET(i)
        ENDDO
        RECORD_BEING_READ = 'SEGMENT'
        DO i = 1,NUMSG
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SEGMENT(i)
        ENDDO
        RECORD_BEING_READ = 'NODE,MOTION,FORCE'
        DO i = 1,NUMRT
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NODE(i),MOTION(i),FORCE(i)
        ENDDO
        RECORD_BEING_READ = 'HEXAH'
        DO i = 1,NUMHX
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) HEXAH(i)
        ENDDO
        RECORD_BEING_READ = 'PENTA'
        DO i = 1,NUMPX
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PENTA(i)
        ENDDO
        RECORD_BEING_READ = 'TETRA'
        DO i = 1,NUMTX
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) TETRA(i)
        ENDDO
        RECORD_BEING_READ = 'LSOLD'
        DO i = 1,NUMLS
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) LSOLD(i)
        ENDDO
        RECORD_BEING_READ = 'LSHEX'
        DO i = 1,NUMLX
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) LSHEX(i)
        ENDDO
        RECORD_BEING_READ = 'LSMBQ'
        DO i = 1,NUMLM
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) LSMBQ(i)
        ENDDO
        RECORD_BEING_READ = 'MEMBT'
        DO i = 1,NUMM3
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) MEMBT(i)
        ENDDO
        RECORD_BEING_READ = 'MEMBQ'
        DO i = 1,NUMM4
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) MEMBQ(i)
        ENDDO
        RECORD_BEING_READ = 'TRUSS'
        DO i = 1,NUMTR
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) TRUSS(i)
        ENDDO
        RECORD_BEING_READ = 'PLATT'
        DO i = 1,NUMP3
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PLATT(i)
        ENDDO
        RECORD_BEING_READ = 'PLATQ'
        DO i = 1,NUMP4
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PLATQ(i)
        ENDDO
        RECORD_BEING_READ = 'BEAM'
        DO i = 1,NUMBM
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) BEAM(i)
        ENDDO
        RECORD_BEING_READ = 'SPRING'
        DO i = 1,NUMSP
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) SPRING(i)
        ENDDO
        RECORD_BEING_READ = 'DAMPER'
        DO i = 1,NUMDM
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) DAMPER(i)
        ENDDO
        RECORD_BEING_READ = 'NNPSETS'
        DO i = 1,NUMNE
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NNPSETS(i)
        ENDDO
        RECORD_BEING_READ = 'NELSETS'
        DO i = 1,NUMEE
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NELSETS(i)
        ENDDO
        RECORD_BEING_READ = 'NSGSETS'
        DO i = 1,NUMSE
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) NSGSETS(i)
        ENDDO
        RECORD_BEING_READ = 'PLATE_PAIR'
        DO i = 1,NUMPP
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) PLATE_PAIR(i)
        ENDDO
        RECORD_BEING_READ = 'CONSTRAINED_NODE'
        DO i = 1,NUMNC
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) CONSTRAINED_NODE(i)
        ENDDO
!!
!! Read bulk write of plate stresses.
!!
        RECORD_BEING_READ = 'STRESS Array'
        DO i = 1,NUMST
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999) (STRESS(k,i),k=1,6)
        ENDDO
!!
!! Skip over element-by-element write of plate stresses.
!!
        RECORD_BEING_READ = 'Element-by-Element Stresses'
        DO i = 1,NUMST
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999)
        ENDDO
!!
!! Read bulk write of state variables.
!!
        IRECORD = 0
        RECORD_BEING_READ = 'STATE_VARIABLES Array'
        READ (IO_UNIT%LRDI,ERR=999) (STATE_VARIABLES(i), i=1,NUMAX)
!!
!! Skip over element-by-element write of state variables.
!!
        NUMBER_OF_ELEMENTS_WRITTEN = NUMHX+NUMPX+NUMTX+NUMLX+NUMLM+            &
     &    NUMM3+NUMM4+NUMTR+NUMP3+NUMP4+NUMBM+NUMSP+NUMDM+NUMSC+NUMVC

        RECORD_BEING_READ = 'Elem-by-Elem State Vars.'
        DO i = 1,NUMBER_OF_ELEMENTS_WRITTEN
          IRECORD = i
          READ (IO_UNIT%LRDI,ERR=999)
        ENDDO
!!
        IRECORD = 0
        RECORD_BEING_READ = 'NRBC_DATA'
        READ (IO_UNIT%LRDI,ERR=999) (NRBC_DATA(i),       i=1,NUMND)
        RECORD_BEING_READ = 'SLIDING_NODE'
        READ (IO_UNIT%LRDI,ERR=999) (SLIDING_NODE(i),    i=1,NUMSN)
        RECORD_BEING_READ = 'CONTACT_SURFACE'
        READ (IO_UNIT%LRDI,ERR=999) (CONTACT_SURFACE(i), i=1,NUMCE)
        RECORD_BEING_READ = 'CONTACT_NODE'
        READ (IO_UNIT%LRDI,ERR=999) (CONTACT_NODE(i),    i=1,NUMCN)
!!
        CLOSE (UNIT=IO_UNIT%LRDI, STATUS='KEEP')
!!
!! Inform user that a restart file has been read.
!!
        WRITE (MSGF,'(1PE12.4)') TIMSIM%Total
        WRITE (MSG1,'(I8)') TIMSIM%Step
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_RESTART_DATA.001.00'//                             &
     &          MSGL//'Restart File Read.'//                                   &
     &          MSGL//'Restarting at Simulation Time:'//MSGF//                 &
     &          MSGL//'Restarting at Simulation Step:'//MSG1                   &
     &          )
      ENDIF
!!
      RETURN
!!
!! Read error exit.
!!
 999    CALL READ_RESTART_DATA_ERROR_EXIT (RECORD_BEING_READ,IRECORD)
      END
!!_
      SUBROUTINE READ_RESTART_DATA_ERROR_EXIT(RECORD_BEING_READ,IRECORD)
!!
!! Copyright (c) by KEY Associates; 31-DEC-1995 11:55:43.00
!!
!! Purpose: Common read-error message processing for reading restart data.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          IRECORD
      CHARACTER                                                                &
     &          RECORD_BEING_READ*(*)
!!
      L = LEN (RECORD_BEING_READ)
      WRITE (MSG1,'(I8)') IRECORD
      CALL USER_MESSAGE                                                        &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_RESTART_DATA_ERROR_EXIT.001.00'//                  &
     &          MSGL//'Problem Reading Restart Record: '                       &
     &              //RECORD_BEING_READ(1:L)//                                 &
     &          MSGL//'Loop Index (If Non-Zero):'//MSG1//                      &
     &          MSGL//'Likely Cause: Too Much Data Requested By READ.'//       &
     &          MSGL//'Possibly A Mix Up In NRS** Values.'                     &
     &          )
!!
      RETURN
      END
!!_
      SUBROUTINE RESTART_CONSISTENCY_CHECK
!!
!! Copyright (c) by KEY Associates; 11-JAN-1996 17:46:59.00
!!
!! Purpose: Check data arrays modified by rezoning.
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
      USE nrbc_data_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE plate_pair_
!!
      USE enumerated_sets_
      USE stress_
      USE state_variables_
      USE coord_
      USE indx_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! First, clear an array to use for scratch space during consistency
!! checking.
!!
      ERROR%COUNT = 0

!$OMP PARALLEL DO
      DO n = 1,NUMRT
        NODE(n)%ICF = 0
      ENDDO
!!
!! Check that each node with rotational components points to a unique,
!! in-bounds location for the rotational degrees of freedom.
!!
      DO n = 1,NUMNP
        IF (NODE(n)%IRT .GT. 0) THEN
          IF (NODE(n)%IRT.LE.NUMNP .OR. NODE(n)%IRT.GT.NUMRT) THEN
            WRITE (MSG1,'(I8)') NODE(n)%ID
            WRITE (MSG2,'(I8)') NODE(n)%IRT
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'RESTART_CONSISTENCY_CHECK.001.01'//                     &
     &          MSGL//'Nodal Point With ID:'//MSG1//                           &
     &          MSGL//'References An Out-Of-Bounds'//                          &
     &          MSGL//'Rotational Location:'//MSG2                             &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE IF (NODE(NODE(n)%IRT)%ICF .NE. 0) THEN
            WRITE (MSG1,'(I8)') NODE(n)%ID
            WRITE (MSG2,'(I8)') NODE(NODE(n)%IRT)%ICF
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'RESTART_CONSISTENCY_CHECK.001.02'//                     &
     &          MSGL//'Nodal Point With ID:'//MSG1//                           &
     &          MSGL//'References A Rotational Location Used By A'//           &
     &          MSGL//'Nodal Point With ID:'//MSG2                             &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            NODE(NODE(n)%IRT)%ICF = NODE(n)%ID
          ENDIF
        ENDIF
      ENDDO
!!
!! Check that each 4-node plate element's translational dof's point
!! to the same place that *.IX(5:8) points.
!!
      DO n = 1,NUMP4
        DO i = 1,4
          IF (NODE(PLATQ(n)%PAR%IX(i))%IRT.NE.PLATQ(n)%PAR%IX(i+4)) THEN
            WRITE (MSG1,'(I8)') PLATQ(n)%PAR%EleID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'RESTART_CONSISTENCY_CHECK.002.00'//                     &
     &          MSGL//'4-Node Plate (P4EL) With ID:'//MSG1//                   &
     &          MSGL//'Has A Mismatch Between Translational'//                 &
     &          MSGL//'And Rotational References --'//                         &
     &          MSGL//'NODE(*.IX(i))%IRT .ne. *.IX(i+4)'                       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDDO
      ENDDO
!!
!! Use midside nodal constraint intialization routine to check on
!! location of midside nodes relative to bounding nodes.
!!
      CALL INITIALIZE_CONSTRAINED_NODES
!!
!! Check for error count not equal to zero.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'RESTART_CONSISTENCY_CHECK.003.00'//                     &
     &          MSGL//'Total Number of Inconsistencies Found:'//MSG1           &
     &          )
      ENDIF
!!
!! Clean up array used for consistency checking.
!!
      DO n = 1,NUMRT
        NODE(n)%ICF = 0
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE RESTART_INTERFACE_TIMES
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 15:06:31
!!
!! Purpose: Apply changes to sliding interface begining and ending times.
!!
      USE shared_common_data
      USE sliding_interface_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered RESTART_INTERFACE_TIMES.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Modify sliding interface *.Begin (calculation) and *.End (calculation) on
!! the basis of values read into INTERFACE_TIME%
!!
      DO i = 1,NUMIT
        Nsi = INTERFACE_TIME%SI(i)%ID
        SLIDING_INTERFACE(Nsi)%Begin = INTERFACE_TIME%SI(i)%Begin
        SLIDING_INTERFACE(Nsi)%End   = INTERFACE_TIME%SI(i)%End
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE READ_INCLUDED_FILES
!!
!! Copyright (c) by KEY Associates, 21-JAN-1992 20:40:37
!!
!! Purpose: Open and read each included file in turn.
!!
      USE shared_common_data
      USE include_file_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          IFINDEX         ! -/S Points to INCLUDE_FILE(IFINDEX)
      LOGICAL                                                                  &
     &          IOERROR
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered READ_INCLUDED_FILES.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Loop over all included files.
!!
      ERROR%COUNT = 0
      DO NIF = 1,NUMIF
      IFINDEX = NIF
!!
!! Open included file.
!!
        IOERROR = .TRUE.
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(NIF)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
        IOERROR = .FALSE.
 200    CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
        IF (IOERROR) THEN
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_INCLUDED_FILE.001.01'//                            &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(NIF)%Full_Name)                        &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of read of included input record file.
!!
        ELSE
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_INCLUDED_FILES.001.02'//                           &
     &          MSGL//'INCLUDE File Read: '                                    &
     &              //'fma3d.in'                                               &
     &          )
!!
!! Read and process included file's input records.
!!
          CALL READ_INPUT_RECORDS (IFINDEX)
!!
!! Disconnect included file.
!!
          CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
!! End of was-open-sucessful if-test.
!!
        ENDIF
!!
!! End of do-loop over included files.
!!
      ENDDO
!!
!! Print total number of input record processing errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_INCLUDED_FILES.004.00'//                           &
     &          MSGL//'Total Number Of Input Record Errors:'//MSG1//           &
     &          MSGL//'Execution Terminated By Program.'                       &
     &          )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_INPUT_RECORDS (IFINDEX)
!!
!! Copyright (c) by KEY Associates, 15-JUL-1990 10:42:13
!!
!! Purpose: Processes the input record file and load into the internal
!! data structures the numeric data defining the analysis to be conducted.
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: IFINDEX  ! Points to INCLUDE_FILE(IFINDEX)

      CHARACTER                                                                &
     &          TEXT*80,                                                       &
     &          INPERR*8,                                                      &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      INTEGER   :: ILCOUNT   ! S/S Input file line counter
      INTEGER   :: ILINDEX   ! S/S Start-of-Record line in input file

      LOGICAL   :: EOF
      LOGICAL   :: FIRST_CALL1 = .TRUE.
      LOGICAL   :: FIRST_CALL2 = .TRUE.

      EXTERNAL C_VALUE
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered READ_INPUT_RECORDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Read input records and process on the basis of keywords.
!!
      IERR = 0
      REWIND (IO_UNIT%LSDI)
      ILCOUNT = 0
      EOF = .FALSE.
      DO WHILE (.NOT.EOF)
        NVALS = MAXRE
        CALL GETIRV                                                            &
     &  (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
        KEY_WORD = C_VALUE(1)
!!
!!
        IF (INDEX(KEY_WORD,'TITLE'         ) .NE. 0) THEN
!!
!! Read "title" for job identification record.
!!
          CALL READ_JOB_IDENTIFICATION (NVALS,TEXT)
!!
        ELSE IF (INDEX(KEY_WORD,'QAREC'    ) .NE. 0) THEN
!!
!! Get user's analysis QA records.
!!
          CALL READ_QA_RECORD (NVALS,TEXT)
!!
        ELSE IF (INDEX(KEY_WORD,'CONTROL'  ) .NE. 0) THEN
!!
!! Get analysis execution control values.
!!
          CALL READ_CONTROL_VALUES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'REZONE'   ) .NE. 0) THEN
!!
!! Get analysis rezoning control values.
!!
          CALL READ_REZONE_CONTROL_VALUES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'INTERFACE') .NE. 0) THEN
!!
!! Get analysis interface begin/end time values.
!!
          CALL READ_INTERFACE_TIME_VALUES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'STOP'     ) .NE. 0) THEN
!!
!! Get analysis stop time.
!!
          CALL READ_ANALYSIS_STOP_TIME (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'PARAM'    ) .NE. 0) THEN
!!
!! Get adjusted algorithm parameter values.
!!
          CALL READ_PARAMETER_VALUES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'PRINT'    ) .NE. 0) THEN
!!
!! Get printed output parameter values.
!!
          CALL READ_PRINTED_OUTPUT_VALUES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'RESULTS'  ) .NE. 0) THEN
!!
!! Get Plotting_Database output specifications.
!!
          CALL READ_OUTPUT_REQUESTS (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'MASSPROP' ) .NE. 0) THEN
!!
!! Get mass property calculation request.
!!
          CALL READ_MASS_PROPERTY_REQUEST (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'GAUGE1D'  ) .NE. 0) THEN
!!
!! Get user specifications for 1-dimensional strain gauge.
!!
          CALL READ_STRAIN_GAUGE1D (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'GAUGE2D'  ) .NE. 0) THEN
!!
!! Get user specifications for 2-dimensional strain gauge.
!!
          CALL READ_STRAIN_GAUGE2D (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'GAUGE3D'  ) .NE. 0) THEN
!!
!! Get user specifications for 3-dimensional strain gauge.
!!
          CALL READ_STRAIN_GAUGE3D (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'MATERIAL' ) .NE. 0) THEN
!!
!! Get material property data.
!!
          CALL READ_MATERIAL_PROPERTIES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'LAYUP'    ) .NE. 0) THEN
!!
!! Get lay-up specifications for layered solid.
!!
          CALL READ_LAYERED_SOLID_LAYUP (FIRST_CALL1,NVALS)
          FIRST_CALL1 = .FALSE.
!!
        ELSE IF (INDEX(KEY_WORD,'PSECTION' ) .NE. 0) THEN
!!
!! Get plate section properties.
!!
          CALL READ_PLATE_SECTION_PROPERTIES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'BSECTION' ) .NE. 0) THEN
!!
!! Get beam and truss section properties.
!!
          CALL READ_BEAM_SECTION_PROPERTIES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'RBODY'    ) .NE. 0) THEN
!!
!! Get user specified rigid body.
!!
          CALL READ_RIGID_BODY (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'RBMASS'   ) .NE. 0) THEN
!!
!! Get user specified concentrated rigid body mass.
!!
          CALL READ_RB_CONCENTRATED_MASS (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'NPMASS'   ) .NE. 0) THEN
!!
!! Get user specified concentrated nodal point mass.
!!
          CALL READ_NP_CONCENTRATED_MASS (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'DISPBC'   ) .NE. 0) THEN
!!
!! Get kinematic boundary condition specification.
!!
          CALL READ_DISPLACEMENT_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'TIEDBC'   ) .NE. 0) THEN
!!
!! Get multipoint kinematic constraint specification.
!!
          CALL READ_TIED_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'SPOTWELD' ) .NE. 0) THEN
!!
!! Get spot weld constraint specification.
!!
          CALL READ_SPOT_WELD (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'WALLBC'   ) .NE. 0) THEN
!!
!! Get rigid wall boundary condition specification.
!!
          CALL READ_RIGID_WALL_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'BODYFORCE') .NE. 0) THEN
!!
!! Get body force specification.
!!
          CALL READ_BODY_FORCE (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'PRESSBC'  ) .NE. 0) THEN
!!
!! Get pressure boundary condition specification.
!!
          CALL READ_PRESSURE_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'FORCEBC'  ) .NE. 0) THEN
!!
!! Get concentrated nodal force boundary condition specification.
!!
          CALL READ_FORCE_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'SPRINGBC' ) .NE. 0) THEN
!!
!! Get "spring" boundary condition specification.
!!
          CALL READ_SPRING_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'DAMPERBC' ) .NE. 0) THEN
!!
!! Get "damper" boundary condition specification.
!!
          CALL READ_DAMPER_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'PERIODBC' ) .NE. 0) THEN
!!
!! Get periodic boundary condition specification.
!!
          CALL READ_PERIODIC_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'NRBC'     ) .NE. 0) THEN
!!
!! Get non-reflecting boundary condition specification.
!!
          CALL READ_NONREFLECTING_BC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'SLIDE'    ) .NE. 0) THEN
!!
!! Get sliding interface specifications.
!!
          CALL READ_SLIDING_INTERFACE (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'TABFTN'   ) .NE. 0) THEN
!!
!! Get tabulated function specification.
!!
          CALL READ_TABULATED_FUNCTION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'NPSET'    ) .NE. 0) THEN
!!
!! Log node set location in input file for later processing.
!!
          CALL LOG_NODE_SET_LOCATION (IFINDEX,ILINDEX)
!!
        ELSE IF (INDEX(KEY_WORD,'ELSET'    ) .NE. 0) THEN
!!
!! Log element set location in input file for later processing.
!!
          CALL LOG_ELEMENT_SET_LOCATION (IFINDEX,ILINDEX)
!!
        ELSE IF (INDEX(KEY_WORD,'SEGSET'   ) .NE. 0) THEN
!!
!! Log segment set location in input file for later processing.
!!
          CALL LOG_SEGMENT_SET_LOCATION (IFINDEX,ILINDEX)
!!
        ELSE IF (INDEX(KEY_WORD,'NPT'      ) .NE. 0) THEN
!!
!! Get nodal point coordinates.
!!
          CALL READ_NODAL_POINT_COORDINATES (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'HXEL'     ) .NE. 0) THEN
!!
!! Get hexahedron element definition.
!!
          CALL READ_HXEL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'PXEL'     ) .NE. 0) THEN
!!
!! Get pentahedron element definition.
!!
          CALL READ_PXEL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'TXEL'     ) .NE. 0) THEN
!!
!! Get tetrahedon element definition.
!!
          CALL READ_TXEL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'LSEL'     ) .NE. 0) THEN
!!
!! Get layered solid element definition.
!!
          CALL READ_LSEL_DEFINITION (FIRST_CALL2,NVALS)
          FIRST_CALL2 = .FALSE.
!!
        ELSE IF (INDEX(KEY_WORD,'M4EL'     ) .NE. 0) THEN
!!
!! Get 4-node membrane element definition.
!!
          CALL READ_M4EL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'M3EL'     ) .NE. 0) THEN
!!
!! Get 3-node membrane element definition.
!!
          CALL READ_M3EL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'TRUSS'     ) .NE. 0) THEN
!!
!! Get truss element definition.
!!
          CALL READ_TRUSS_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'P4EL'     ) .NE. 0) THEN
!!
!! Get 4-node plate element definition.
!!
          CALL READ_P4EL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'P3EL'     ) .NE. 0) THEN
!!
!! Get 3-node plate element definition.
!!
          CALL READ_P3EL_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'BEAM'     ) .NE. 0) THEN
!!
!! Get beam element definition.
!!
          CALL READ_BEAM_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'SPRING'   ) .NE. 0) THEN
!!
!! Get spring element definition.
!!
          CALL READ_SPRING_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'DAMPER'   ) .NE. 0) THEN
!!
!! Get damper element definition.
!!
          CALL READ_DAMPER_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'SEGMENT'  ) .NE. 0) THEN
!!
!! Get segment definition.
!!
          CALL READ_SEGMENT_DEFINITION (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'VELIC'    ) .NE. 0) THEN
!!
!! Get initial conditions on translational and rotational velocity.
!!
          CALL READ_VELOCITY_IC (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'NPCON1'   ) .NE. 0) THEN
!!
!! Get initial conditions on translational and rotational velocity.
!!
          CALL READ_NODAL_CONSTRAINT_1 (NVALS)
!!
        ELSE IF (INDEX(KEY_WORD,'INCLUDE'  ) .NE. 0) THEN
!!
!! Do nothing (and also avoid collecting unrecognized key_word error).
!!
        ELSE IF (INDEX(KEY_WORD,'ENDDATA'  ) .NE. 0) THEN
          EOF = .TRUE.
        ELSE
!!
!! Unrecognized keyword.
!!
          WRITE (MSG1,'(I8)') ILCOUNT
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_INPUT_RECORDS.001.01'//                            &
     &          MSGL//'Unrecognized Key Word In Input: '//KEY_WORD//           &
     &          MSGL//'Found While Reading Line Number:'//MSG1//               &
     &          MSGL//'Reading Will Continue to End Of Input.'                 &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
!!
        ENDIF
      ENDDO
!!
      ERROR%COUNT = ERROR%COUNT + IERR
!!
      RETURN
      END
!!_
      SUBROUTINE READ_CONTROL_VALUES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      DO i = 2,NVALS,2
        KEY_WORD = C_VALUE(i)
        SELECT CASE (KEY_WORD)
          CASE ('SPRINT')
            CONTROL%SPRINT = I_VALUE(i+1)
          CASE ('INPUT_ECHO')
            CONTROL%INECHO = I_VALUE(i+1)
          CASE ('WR_RESTART')
            CONTROL%WRSTAR = I_VALUE(i+1)
          CASE ('RD_RESTART')
            CONTROL%RDSTAR = I_VALUE(i+1)
          CASE ('WR_MESH')
            CONTROL%WRMESH = I_VALUE(i+1)
          CASE ('RD_MESH')
            CONTROL%RDMESH = I_VALUE(i+1)
          CASE ('DATA_CHECK')
            CONTROL%DCHECK = I_VALUE(i+1)
          CASE ('DT_TAB_FTN')
            CONTROL%DTTABF = I_VALUE(i+1)
          CASE ('SUBCYCLING')
            CONTROL%SUBCYC = I_VALUE(i+1)
          CASE ('POLAR_DECOMP')
            CONTROL%POLARD = I_VALUE(i+1)
          CASE ('MID_INTERVAL')
            CONTROL%MIDINT = I_VALUE(i+1)
          CASE ('BT_PLATE')
            CONTROL%BTPLTQ = I_VALUE(i+1)
          CASE ('REZONING')
            CONTROL%REZONE = I_VALUE(i+1)
          CASE ('RZ_RESTART')
            CONTROL%RZSTAR = I_VALUE(i+1)
!!
!! Unexpected execution control key word found.
!!
          CASE DEFAULT
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_CONTROL_VALUES'//                                  &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
        END SELECT
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE READ_REZONE_CONTROL_VALUES (NVALS)
!!
!! Copyright (c) by KEY Associates, 8-MAR-1991 21:27:27
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      DO i = 2,NVALS,2
        KEY_WORD = C_VALUE(i)
        SELECT CASE (KEY_WORD)
          CASE ('CHK_INTERVAL')
            REZONE%INTERVAL = I_VALUE(i+1)
          CASE ('MAX_FLAGS')
            REZONE%MAXCNT   = I_VALUE(i+1)
          CASE ('MAX_FLAG_CYC')
            REZONE%MAXFCY   = I_VALUE(i+1)
          CASE ('MAX_LEVELS')
            REZONE%MAXLEV   = I_VALUE(i+1)
          CASE ('WR_MESH')
            REZONE%WRMESH   = I_VALUE(i+1)
!!
!! Unexpected printed output key word found.
!!
          CASE DEFAULT
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_REZONE_CONTROL_VALUES.001.00'//                    &
     &          MSGL//'REZONE Input Record Contains'//                         &
     &          MSGL//'Unexpected Parameter Word: '//KEY_WORD                  &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
        END SELECT
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE READ_INTERFACE_TIME_VALUES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .GT. INTERFACE_TIME%Number_of_Interfaces) THEN
!!
!! Sliding interface control record storage exceeded.
!!
        WRITE (MSG1,'(I8)') INTERFACE_TIME%Number_of_Interfaces
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_INTERFACE_TIME_VALUES.001.01'//                    &
     &          MSGL//'More INTERFACE Records Read Than'//                     &
     &          MSGL//'Can Be Stored. Current Limit:'//MSG1                    &
     &          )
        ERROR%COUNT = ERROR%COUNT + 1
      ELSE IF (NCOUNT .LE. NUMIT) THEN
        DO i = 2,NVALS,2
          KEY_WORD = C_VALUE(i)
          SELECT CASE (KEY_WORD)
            CASE ('INTERFACE_ID')
              INTERFACE_TIME%SI(NCOUNT)%ID    = I_VALUE(i+1)
            CASE ('BEGIN')
              INTERFACE_TIME%SI(NCOUNT)%BEGIN = D_VALUE(i+1)
            CASE ('END')
              INTERFACE_TIME%SI(NCOUNT)%END   = D_VALUE(i+1)
!!
!! Unexpected interface time control control key word found.
!!
            CASE DEFAULT
              CALL USER_MESSAGE                                                &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'READ_INTERFACE_TIME_VALUES.001.02'//                  &
     &            MSGL//'Unexpected INTERFACE Keyword Found: '//KEY_WORD       &
     &            )
              ERROR%COUNT = ERROR%COUNT + 1
          END SELECT
        ENDDO
      ELSE
!!
!! More interface time records to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PARAMETER_VALUES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      DO i = 2,NVALS,2
        KEY_WORD = C_VALUE(i)
        SELECT CASE (KEY_WORD)
          CASE ('LINEAR_Q')
            PARAMVALUE%Blk1      = D_VALUE(i+1)
          CASE ('QUADRATIC_Q')
            PARAMVALUE%Blk2      = D_VALUE(i+1)
          CASE ('HG_VISCOSITY')
            PARAMVALUE%HGV       = D_VALUE(i+1)
          CASE ('HG_STIFFNESS')
            PARAMVALUE%HGK       = D_VALUE(i+1)
          CASE ('DT_RATIO')
            PARAMVALUE%DTratio   = D_VALUE(i+1)
          CASE ('DT_SCALE')
            PARAMVALUE%DTscale   = D_VALUE(i+1)
          CASE ('DT_GROW')
            PARAMVALUE%DTgrow    = D_VALUE(i+1)
          CASE ('STATUS')
            PARAMVALUE%STATUS    = D_VALUE(i+1)
          CASE ('ENG_BAL')
            PARAMVALUE%ENG_BAL   = D_VALUE(i+1)
          CASE ('SI_FACTOR')
            PARAMVALUE%Factor    = D_VALUE(i+1)
          CASE ('SI_CAPTURE')
            PARAMVALUE%Capture   = D_VALUE(i+1)
          CASE ('SI_BORDER')
            PARAMVALUE%Border    = D_VALUE(i+1)
          CASE ('SI_SORT_FREQ')
            PARAMVALUE%Sort_Freq = D_VALUE(i+1)
          CASE ('CYCLE_RATIO')
            PARAMVALUE%Cycle_R   = D_VALUE(i+1)
          CASE ('CYCLE_SPREAD')
            PARAMVALUE%Cycle_S   = D_VALUE(i+1)
          CASE ('CYCLE_FREQ')
            PARAMVALUE%Cycle_F   = D_VALUE(i+1)
          CASE ('MAT32_STEPS')
            PARAMVALUE%Max_Steps = D_VALUE(i+1)
          CASE ('NRBC_SCALE')
            PARAMVALUE%NRBC_Q    = D_VALUE(i+1)
          CASE ('LSOLD_EQ_TOL')
            PARAMVALUE%EQUI_TOL  = D_VALUE(i+1)
          CASE ('MAT22_ROTATE')
            PARAMVALUE%FIBROT    = D_VALUE(i+1)
          CASE ('BCS_MAX_ITER')
            PARAMVALUE%BCSMXIT   = D_VALUE(i+1)
          CASE ('BCS_GROW_LIM')
            PARAMVALUE%BCSGRLM   = D_VALUE(i+1)
          CASE ('BCS_RELAX')
            PARAMVALUE%BCSRLAX   = D_VALUE(i+1)
          CASE ('RZ_MAX_CPDOT')
            PARAMVALUE%CPDMAX    = D_VALUE(i+1)
          CASE ('RZ_MIN_CPDOT')
            PARAMVALUE%CPDMIN    = D_VALUE(i+1)
          CASE ('RZ_MAX_CPCHG')
            PARAMVALUE%CSPMAX    = D_VALUE(i+1)
          CASE ('RZ_MIN_CPCHG')
            PARAMVALUE%CSPMIN    = D_VALUE(i+1)
!!
!! Unexpected parameter key word found.
!!
          CASE DEFAULT
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_PARAMETER_VALUES'//                                &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
        END SELECT
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PRINTED_OUTPUT_VALUES (NVALS)
!!
!! Copyright (c) by KEY Associates, 8-MAR-1991 21:27:27
!!
      USE shared_common_data
      USE value_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      DO i = 2,NVALS,2
        KEY_WORD = C_VALUE(i)
        SELECT CASE (KEY_WORD)
          CASE ('BEGIN')
            PRINT%Begin  = D_VALUE(i+1)
          CASE ('END')
            PRINT%End    = D_VALUE(i+1)
          CASE ('DELTA')
            PRINT%Delta  = D_VALUE(i+1)
          CASE ('TIME')
            PRINT%Time   = D_VALUE(i+1)
          CASE DEFAULT
            IF (VALUE(i+1)%VTYP .EQ. 'C') THEN
              KEY_WORD = C_VALUE(i+1)
              IF (KEY_WORD .EQ. 'ALL') THEN
                IVAL = -1
              ELSE
!!
!! Unexpected printed output key word found.
!!
                CALL USER_MESSAGE                                              &
     &                  (                                                      &
     &                  MSGL//'WARN'//                                         &
     &                  MSGL//'READ_PRINTED_OUTPUT_VALUES.001.01'//            &
     &                  MSGL//'PRINT Input Record Contains'//                  &
     &                  MSGL//'Unexpected Keyword: '//KEY_WORD//               &
     &                  MSGL//'Keyword Expected: ALL'                          &
     &                  )
                ERROR%COUNT = ERROR%COUNT + 1
              ENDIF
            ELSE
              IVAL = I_VALUE(i+1)
            ENDIF
            KEY_WORD = C_VALUE(i)
            SELECT CASE (KEY_WORD)
              CASE ('NPT')
                PRINT%NODES  = IVAL
              CASE ('HXEL')
                PRINT%HEXAH  = IVAL
              CASE ('PXEL')
                PRINT%PENTA  = IVAL
              CASE ('TXEL')
                PRINT%TETRA  = IVAL
              CASE ('M3EL')
                PRINT%MEMBT  = IVAL
              CASE ('M4EL')
                PRINT%MEMBQ  = IVAL
              CASE ('TRUSS')
                PRINT%TRUSS  = IVAL
              CASE ('P3EL')
                PRINT%PLATT  = IVAL
              CASE ('P4EL')
                PRINT%PLATQ  = IVAL
              CASE ('BEAM')
                PRINT%BEAMS  = IVAL
              CASE ('SPRING')
                PRINT%SPRING = IVAL
              CASE ('DAMPER')
                PRINT%DAMPER = IVAL
!!
!! Unexpected printed output key word found.
!!
              CASE DEFAULT
                CALL USER_MESSAGE                                              &
     &                  (                                                      &
     &                  MSGL//'WARN'//                                         &
     &                  MSGL//'READ_PRINTED_OUTPUT_VALUES.001.02'//            &
     &                  MSGL//'PRINT Input Record Contains'//                  &
     &                  MSGL//'Unexpected Parameter Word: '//KEY_WORD          &
     &                  )
                ERROR%COUNT = ERROR%COUNT + 1
            END SELECT
        END SELECT
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE READ_OUTPUT_REQUESTS (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE results_
      USE output_
      USE value_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER(32)   :: C_VALUE
      CHARACTER(12)   :: KEY_WORD
      INTEGER, SAVE   :: NCOUNT = 0
      LOGICAL         :: FOUND
      LOGICAL, SAVE   :: FIRST = .TRUE.
      EXTERNAL        C_VALUE
!!
      IF (FIRST) THEN
        OUTPUT%NAME = (/'BEGIN ','END   ','DELTA ','TIME  ','MASSP ',          &
     &         'GAUGE ','NPT   ','HXEL  ','PXEL  ','TXEL  ','LSEL  ',          &
     &         'M3EL  ','M4EL  ','TRUSS ','P3EL  ','P4EL  ','BEAM  ',          &
     &         'SPRING','DAMPER','GLOBAL'/)
        FIRST = .FALSE.
      ENDIF
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMRF) THEN
        RESULTS(NCOUNT)%ResID = I_VALUE( 2)
        RESULTS(NCOUNT)%Title = C_VALUE( 3)
        RESULTS(NCOUNT)%OpfID = I_VALUE( 4)
        DO i = 5,NVALS,2
          KEY_WORD = C_VALUE(i)
          SELECT CASE (KEY_WORD)
            CASE ('BEGIN')
              RESULTS(NCOUNT)%Begin = D_VALUE(i+1)
            CASE ('END')
              RESULTS(NCOUNT)%End   = D_VALUE(i+1)
            CASE ('DELTA')
              RESULTS(NCOUNT)%Delta = D_VALUE(i+1)
            CASE ('TIME')
              RESULTS(NCOUNT)%Time  = D_VALUE(i+1)
            CASE DEFAULT
              IF (VALUE(i+1)%VTYP .EQ. 'C') THEN
                KEY_WORD = C_VALUE(i+1)
                IF (KEY_WORD .EQ. 'ALL') THEN
                  IVAL = -1
                ELSE
!!
!! Unexpected results output key word found.
!!
                  WRITE (MSG1,'(I8)') RESULTS(NCOUNT)%ResID
                  CALL USER_MESSAGE                                            &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'READ_OUTPUT_REQUESTS.001.01'//                        &
     &            MSGL//'RESULTS Input Record ID:'//MSG1//                     &
     &            MSGL//'Contains Unexpected Keyword: '//KEY_WORD//            &
     &            MSGL//'Keyword Expected: ALL'                                &
     &            )
                  ERROR%COUNT = ERROR%COUNT + 1
                ENDIF
              ELSE
                IVAL = I_VALUE(i+1)
              ENDIF
              KEY_WORD = C_VALUE(i)
              SELECT CASE (KEY_WORD)
                CASE ('MASSP')
                  RESULTS(NCOUNT)%MASSP  = IVAL
                CASE ('GAUGE')
                  RESULTS(NCOUNT)%GAUGE  = IVAL
                CASE ('NPT')
                  RESULTS(NCOUNT)%NODES  = IVAL
                CASE ('HXEL')
                  RESULTS(NCOUNT)%HEXAH  = IVAL
                CASE ('PXEL')
                  RESULTS(NCOUNT)%PENTA  = IVAL
                CASE ('TXEL')
                  RESULTS(NCOUNT)%TETRA  = IVAL
                CASE ('LSEL')
                  RESULTS(NCOUNT)%LSOLD  = IVAL
                CASE ('M3EL')
                  RESULTS(NCOUNT)%MEMBT  = IVAL
                CASE ('M4EL')
                  RESULTS(NCOUNT)%MEMBQ  = IVAL
                CASE ('TRUSS')
                  RESULTS(NCOUNT)%TRUSS  = IVAL
                CASE ('P3EL')
                  RESULTS(NCOUNT)%PLATT  = IVAL
                CASE ('P4EL')
                  RESULTS(NCOUNT)%PLATQ  = IVAL
                CASE ('BEAM')
                  RESULTS(NCOUNT)%BEAMS  = IVAL
                CASE ('SPRING')
                  RESULTS(NCOUNT)%SPRING = IVAL
                CASE ('DAMPER')
                  RESULTS(NCOUNT)%DAMPER = IVAL
                CASE ('GLOBAL')
                  RESULTS(NCOUNT)%GLOBAL = IVAL
!!
!! Unexpected key word found.
!!
                CASE DEFAULT
                  WRITE (MSG1,'(I8)') RESULTS(NCOUNT)%ResID
                  CALL USER_MESSAGE                                            &
     &        (                                                                &
     &        MSGL//'WARN'//                                                   &
     &        MSGL//'READ_OUTPUT_REQUESTS.001.02'//                            &
     &        MSGL//'RESULTS Input Record ID:'//MSG1//                         &
     &        MSGL//'Contains Unexpected Parameter Word: '//KEY_WORD           &
     &        )
                  ERROR%COUNT = ERROR%COUNT + 1
              END SELECT
          END SELECT
        ENDDO
      ELSE
!!
!! More results output requests to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_ANALYSIS_STOP_TIME (NVALS)
!!
!! Copyright (c) by KEY Associates, 26-FEB-1991 21:03:29
!!
      USE shared_common_data
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      TIMSIM%Stop = D_VALUE(2)
!!
      RETURN
      END
!!_
      SUBROUTINE READ_STRAIN_GAUGE1D (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE gauge1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMG1) THEN
        GAUGE1D(NCOUNT)%PAR%GauID    = I_VALUE(2)
        GAUGE1D(NCOUNT)%PAR%TBflg    = I_VALUE(3)
        IF (GAUGE1D(NCOUNT)%PAR%TBflg .GT. 0) THEN
          GAUGE1D(NCOUNT)%PAR%EleID  = I_VALUE(4)
          GAUGE1D(NCOUNT)%PAR%GauLoc = I_VALUE(5)
        ELSE
          GAUGE1D(NCOUNT)%PAR%IX(1)  = I_VALUE(4)
          GAUGE1D(NCOUNT)%PAR%IX(2)  = I_VALUE(5)
          GAUGE1D(NCOUNT)%PAR%IX(3)  = I_VALUE(6)
          GAUGE1D(NCOUNT)%PAR%IX(4)  = I_VALUE(7)
        ENDIF
      ELSE
!!
!! More 1-D strain gauge specifications to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_STRAIN_GAUGE2D (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE gauge2d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMG2) THEN
        GAUGE2D(NCOUNT)%PAR%GauID      = I_VALUE( 2)
        GAUGE2D(NCOUNT)%PAR%MPQTflg    = I_VALUE( 3)
        IF (GAUGE2D(NCOUNT)%PAR%MPQTflg .GT. 0) THEN
          GAUGE2D(NCOUNT)%PAR%EleID    = I_VALUE( 4)
          GAUGE2D(NCOUNT)%PAR%GauLoc   = I_VALUE( 5)
        ELSE
          GAUGE2D(NCOUNT)%PAR%IX(1)    = I_VALUE( 4)
          GAUGE2D(NCOUNT)%PAR%IX(2)    = I_VALUE( 5)
          GAUGE2D(NCOUNT)%PAR%IX(3)    = I_VALUE( 6)
          GAUGE2D(NCOUNT)%PAR%IX(4)    = I_VALUE( 7)
        ENDIF
      ELSE
!!
!! More 2-D strain gauge specifications to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_STRAIN_GAUGE3D (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE gauge3d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMG3) THEN
        GAUGE3D(NCOUNT)%PAR%GauID    = I_VALUE( 2)
        GAUGE3D(NCOUNT)%PAR%HPTflg   = I_VALUE( 3)
        IF (GAUGE3D(NCOUNT)%PAR%HPTflg .GT. 0) THEN
          GAUGE3D(NCOUNT)%PAR%EleID  = I_VALUE( 4)
        ELSE
          GAUGE3D(NCOUNT)%PAR%IX(1)  = I_VALUE( 4)
          GAUGE3D(NCOUNT)%PAR%IX(2)  = I_VALUE( 5)
          GAUGE3D(NCOUNT)%PAR%IX(3)  = I_VALUE( 6)
          GAUGE3D(NCOUNT)%PAR%IX(4)  = I_VALUE( 7)
          GAUGE3D(NCOUNT)%PAR%IX(5)  = I_VALUE( 8)
          GAUGE3D(NCOUNT)%PAR%IX(6)  = I_VALUE( 9)
          GAUGE3D(NCOUNT)%PAR%IX(7)  = I_VALUE(10)
          GAUGE3D(NCOUNT)%PAR%IX(8)  = I_VALUE(11)
        ENDIF
      ELSE
!!
!! More 3-D strain gauge specifications to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_MASS_PROPERTY_REQUEST (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE massprop_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMMP) THEN
        MASSPROP(NCOUNT)%MPID   = I_VALUE( 2)
        MASSPROP(NCOUNT)%SetID  = I_VALUE( 3)
        MASSPROP(NCOUNT)%Irot   = I_VALUE( 4)
        MASSPROP(NCOUNT)%Xnert  = D_VALUE( 5)
        MASSPROP(NCOUNT)%Ynert  = D_VALUE( 6)
        MASSPROP(NCOUNT)%Znert  = D_VALUE( 7)
        MASSPROP(NCOUNT)%Vxnert = D_VALUE( 8)
        MASSPROP(NCOUNT)%Vynert = D_VALUE( 9)
        MASSPROP(NCOUNT)%Vznert = D_VALUE(10)
      ELSE
!!
!! More mass property output requests to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_MATERIAL_PROPERTIES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE material_
      USE property_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      LOGICAL, SAVE :: FOUND, FIRST = .TRUE.

      EXTERNAL                                                                 &
     &          C_VALUE
!!
      IF (FIRST) THEN
        CALL INITIALIZE_PROPERTY_NAMES
        FIRST = .FALSE.
      ENDIF
!!
!! Define names and pointers in the PROPERTY structure.
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMMT) THEN
        MATERIAL(NCOUNT)%MatID = I_VALUE(2)
        MATERIAL(NCOUNT)%Label = C_VALUE(3)
        MATERIAL(NCOUNT)%Type  = I_VALUE(4)
        DO i = 5,NVALS,2
          FOUND = .FALSE.
          KEY_WORD = C_VALUE(i)
          DO n = 1,PROPERTY%Number_of_Entries
            IF (PROPERTY%NAME(n) .EQ. KEY_WORD) THEN
              m = PROPERTY%Location(n)
              MATERIAL(NCOUNT)%PVAL(m) = D_VALUE(i+1)
              FOUND = .TRUE.
            ENDIF
          ENDDO
          IF (.NOT.FOUND) THEN
!!
!! Unexpected key word found.
!!
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_MATERIAL_PROPERTIES'//                             &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ENDIF
        ENDDO
      ELSE
!!
!! More material property requests to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_PROPERTY_NAMES
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE property_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Define sets, names and pointers in structure PROPERTY.
!!
! 1. Mass density
      PROPERTY%Set( 1)        = 00
      PROPERTY%Name( 1)       = 'DENSITY'
      PROPERTY%Location( 1)   = 01

! 2. Linear bulk viscosity coefficient
      PROPERTY%Set( 2)        = 00
      PROPERTY%Name( 2)       = 'LINEAR_Q'
      PROPERTY%Location( 2)   = 02

! 3. Quadratic bulk viscosity coefficient
      PROPERTY%Set( 3)        = 00
      PROPERTY%Name( 3)       = 'QUADRATIC_Q'
      PROPERTY%Location( 3)   = 03

! 4. Hourglass viscosity coefficient
      PROPERTY%Set( 4)        = 00
      PROPERTY%Name( 4)       = 'HG_VISCOSITY'
      PROPERTY%Location( 4)   = 04

! 5. Hourglass stiffness coefficient
      PROPERTY%Set( 5)        = 00
      PROPERTY%Name( 5)       = 'HG_STIFFNESS'
      PROPERTY%Location( 5)   = 05

! MATERIAL_S1
!  6. Spring constant/tabulated function ID
      PROPERTY%Set( 6)        = 01
      PROPERTY%Name( 6)       = 'K1'
      PROPERTY%Location( 6)   = 06

! MATERIAL_S3
      PROPERTY%Set( 7)        = 02
      PROPERTY%Name( 7)       = 'K1'
      PROPERTY%Location( 7)   = 06
      PROPERTY%Set( 8)        = 02
      PROPERTY%Name( 8)       = 'YIELD'
      PROPERTY%Location( 8)   = 10
      PROPERTY%Set( 9)        = 02
      PROPERTY%Name( 9)       = 'HARD'
      PROPERTY%Location( 9)   = 11
      PROPERTY%Set(10)        = 02
      PROPERTY%Name(10)       = 'BETA'
      PROPERTY%Location(10)   = 13
      PROPERTY%Set(11)        = 02
      PROPERTY%Name(11)       = 'P_EXPONENT'
      PROPERTY%Location(11)   = 16
      PROPERTY%Set(12)        = 02
      PROPERTY%Name(12)       = 'D_SCALING'
      PROPERTY%Location(12)   = 17

! MATERIAL_S3
!  6. Spring constant/tabulated function ID
      PROPERTY%Set(13)        = 03
      PROPERTY%Name(13)       = 'K1'
      PROPERTY%Location(13)   = 06

! MATERIAL_D6
      PROPERTY%Set(14)        = 06
      PROPERTY%Name(14)       = 'D1'
      PROPERTY%Location(14)   = 06

! MATERIAL_D7
      PROPERTY%Set(15)        = 07
      PROPERTY%Name(15)       = 'D1'
      PROPERTY%Location(15)   = 06
      PROPERTY%Set(16)        = 07
      PROPERTY%Name(16)       = 'D2'
      PROPERTY%Location(16)   = 07

! MATERIAL_D8
      PROPERTY%Set(17)        = 08
      PROPERTY%Name(17)       = 'D1'
      PROPERTY%Location(17)   = 06

! MATERIAL_10,20,30,40,50 (Elastic)
! 6. Young's modulus
! 7. Poisson's ratio
! 8. Lamé parameter, Lambda             (derived)
! 9. Lamé parameter, Mu (Shear modulus) (derived)
      PROPERTY%Set(18)        = 10
      PROPERTY%Name(18)       = 'E'
      PROPERTY%Location(18)   = 06
      PROPERTY%Set(19)        = 10
      PROPERTY%Name(19)       = 'NU'
      PROPERTY%Location(19)   = 07

! MATERIAL_11,21,31,41,51 (Plastic)
!  6. Young's modulus
!  7. Poisson's ratio
!  8. Lamé parameter, Lambda            (derived)
!  9. Lamé parameter, Mu (Shear modulus)(derived)
! 10. Uniaxial-stress yield stress
! 11. Linear hardening modulus
! 12. Plastic modulus
! 13. Isotropic/kinematic hardening mix
! 14. Compression/tension cut-off stress
! 15. Stress type (comp./two-way/tension)=(-1/0/1)
! 16. Rate sensitivty exponent
! 17. Rate sensitivty scaling
      PROPERTY%Set(20)        = 11
      PROPERTY%Name(20)       = 'E'
      PROPERTY%Location(20)   = 06
      PROPERTY%Set(21)        = 11
      PROPERTY%Name(21)       = 'NU'
      PROPERTY%Location(21)   = 07
      PROPERTY%Set(22)        = 11
      PROPERTY%Name(22)       = 'YIELD'
      PROPERTY%Location(22)   = 10
      PROPERTY%Set(23)        = 11
      PROPERTY%Name(23)       = 'HARD'
      PROPERTY%Location(23)   = 11
      PROPERTY%Set(24)        = 11
      PROPERTY%Name(24)       = 'BETA'
      PROPERTY%Location(24)   = 13
      PROPERTY%Set(25)        = 11
      PROPERTY%Name(25)       = 'P_EXPONENT'
      PROPERTY%Location(25)   = 16
      PROPERTY%Set(26)        = 11
      PROPERTY%Name(26)       = 'D_SCALING'
      PROPERTY%Location(26)   = 17

! MATERIAL_22, (Bi-directional fabric)
!  6. Fiber A, linear    stiffness
!  7. Fiber A, quadratic stiffness
!  8. Fiber B, linear    stiffness
!  9. Fiber B, quadratic stiffness
! 10. Shear coefficient of friction
! 11. Bulk modulus in compression
! 12. Fully compacted aerial strain
! 13. Shear model (0,1=friction,elastic)
! 14. dG/dP, P = sqrt (Stress_A x Stress_B)
      PROPERTY%Set(27)        = 22
      PROPERTY%Name(27)       = 'FIBER_A_L'
      PROPERTY%Location(27)   = 06
      PROPERTY%Set(28)        = 22
      PROPERTY%Name(28)       = 'FIBER_A_Q'
      PROPERTY%Location(28)   = 07
      PROPERTY%Set(29)        = 22
      PROPERTY%Name(29)       = 'FIBER_B_L'
      PROPERTY%Location(29)   = 08
      PROPERTY%Set(30)        = 22
      PROPERTY%Name(30)       = 'FIBER_B_Q'
      PROPERTY%Location(30)   = 09
      PROPERTY%Set(31)        = 22
      PROPERTY%Name(31)       = 'FRICTION'
      PROPERTY%Location(31)   = 10
      PROPERTY%Set(32)        = 22
      PROPERTY%Name(32)       = 'FABRIC_BULK'
      PROPERTY%Location(32)   = 11
      PROPERTY%Set(33)        = 22
      PROPERTY%Name(33)       = 'F_COMPACTED'
      PROPERTY%Location(33)   = 12
      PROPERTY%Set(103)       = 22
      PROPERTY%Name(103)      = 'SHEAR_MODEL'
      PROPERTY%Location(103)  = 13
      PROPERTY%Set(104)       = 22
      PROPERTY%Name(104)      = 'FIBER_G_COEF'
      PROPERTY%Location(104)  = 14

! MATERIAL_25,35,45 (Orthotropic elastic)
!  6. Young's modulus, a-axis
!  7. Young's modulus, b-axis
!  8. Young's modulus, c-axis
!  9. Poisson's ratio, Daa=-Nuba(Sbb/Ebb)
! 10. Poisson's ratio, Daa=-Nuca(Scc/Ecc)
! 11. Poisson's ratio, Dbb=-Nucb(Scc/Ecc)
! 12. Shear modulus
! 13. Shear modulus
! 14. Shear modulus
! MATERIAL_25,35,45 (dervied values)
!  6. Modulus, a-axis
!  7. Modulus, b-axis
!  8. Modulus, c-axis
!  9. Coupling modulus
! 10. Coupling modulus
! 11. Coupling modulus
! 15. Uniaxial-strain stiffness
! 16. Planar-shear-strain stiffness
      PROPERTY%Set(34)        = 25
      PROPERTY%Name(34)       = 'EAA'
      PROPERTY%Location(34)   = 06
      PROPERTY%Set(35)        = 25
      PROPERTY%Name(35)       = 'EBB'
      PROPERTY%Location(35)   = 07
      PROPERTY%Set(36)        = 25
      PROPERTY%Name(36)       = 'ECC'
      PROPERTY%Location(36)   = 08
      PROPERTY%Set(37)        = 25
      PROPERTY%Name(37)       = 'NUBA'
      PROPERTY%Location(37)   = 09
      PROPERTY%Set(38)        = 25
      PROPERTY%Name(38)       = 'NUCA'
      PROPERTY%Location(38)   = 10
      PROPERTY%Set(39)        = 25
      PROPERTY%Name(39)       = 'NUCB'
      PROPERTY%Location(39)   = 11
      PROPERTY%Set(40)        = 25
      PROPERTY%Name(40)       = 'GAB'
      PROPERTY%Location(40)   = 12
      PROPERTY%Set(41)        = 25
      PROPERTY%Name(41)       = 'GAC'
      PROPERTY%Location(41)   = 13
      PROPERTY%Set(42)        = 25
      PROPERTY%Name(42)       = 'GBC'
      PROPERTY%Location(42)   = 14

! MATERIAL_32 (Smooth hardening plasticity)
!  6. Young's modulus   (see MATERIAL_31)
!  7. Poisson's ratio           "
!  8. Lamé parameter, Lambda    "
!  9. Lamé parameter, Mu        "
! 10. Yield stress              "
! 11. Kinematic hardening modulus
! 12. Isotropic hardening modulus
! 13. Hardening parameter p
! 14. Hardening parameter r
      PROPERTY%Set(43)        = 32
      PROPERTY%Name(43)       = 'E'
      PROPERTY%Location(43)   = 06
      PROPERTY%Set(44)        = 32
      PROPERTY%Name(44)       = 'NU'
      PROPERTY%Location(44)   = 07
      PROPERTY%Set(45)        = 32
      PROPERTY%Name(45)       = 'YIELD'
      PROPERTY%Location(45)   = 10
      PROPERTY%Set(46)        = 32
      PROPERTY%Name(46)       = 'E_KINEMATIC'
      PROPERTY%Location(46)   = 11
      PROPERTY%Set(47)        = 32
      PROPERTY%Name(47)       = 'E_ISOTROPIC'
      PROPERTY%Location(47)   = 12
      PROPERTY%Set(48)        = 32
      PROPERTY%Name(48)       = 'P_PARAMETER'
      PROPERTY%Location(48)   = 13
      PROPERTY%Set(49)        = 32
      PROPERTY%Name(49)       = 'R_PARAMETER'
      PROPERTY%Location(49)   = 14

! MATERIAL_33 (Soil/rock cap model)
!  9. Shear Modulus
! 10. Bulk modulus
! 11. Failure envelope constant
! 12. Failure envelope constant
! 13. Failure envelope constant
! 14. Cap envelope constant
! 15. Cap envelope constant
! 16. State function constant
! 17. State function constant
! 18. State function constant
! 19. Tensile cutoff, tension positive
! 20. Failure envelope intersection w/ J1
! 21. Tensile cutoff parameter, TCUT = 3*TMAX
! 22. Derived, Starting value for EL
      PROPERTY%Set(50)        = 33
      PROPERTY%Name(50)       = 'SHEAR'
      PROPERTY%Location(50)   = 09
      PROPERTY%Set(51)        = 33
      PROPERTY%Name(51)       = 'BULK'
      PROPERTY%Location(51)   = 10
      PROPERTY%Set(52)        = 33
      PROPERTY%Name(52)       = 'FAILURE_A'
      PROPERTY%Location(52)   = 11
      PROPERTY%Set(53)        = 33
      PROPERTY%Name(53)       = 'FAILURE_B'
      PROPERTY%Location(53)   = 12
      PROPERTY%Set(54)        = 33
      PROPERTY%Name(54)       = 'FAILURE_C'
      PROPERTY%Location(54)   = 13
      PROPERTY%Set(55)        = 33
      PROPERTY%Name(55)       = 'CAP_R'
      PROPERTY%Location(55)   = 14
      PROPERTY%Set(56)        = 33
      PROPERTY%Name(56)       = 'CAP_POSITION'
      PROPERTY%Location(56)   = 15
      PROPERTY%Set(57)        = 33
      PROPERTY%Name(57)       = 'VOLUME_D'
      PROPERTY%Location(57)   = 16
      PROPERTY%Set(58)        = 33
      PROPERTY%Name(58)       = 'VOLUME_W'
      PROPERTY%Location(58)   = 17
      PROPERTY%Set(59)        = 33
      PROPERTY%Name(59)       = 'LTYPE'
      PROPERTY%Location(59)   = 18
      PROPERTY%Set(60)        = 33
      PROPERTY%Name(60)       = 'TMAX'
      PROPERTY%Location(60)   = 19

! MATERIAL_36 (Linear viscoelastic)
!  6. Instantaneous bulk modulus
!  7. Long-term bulk modulus
!  8. Bulk exponential decay constant
!  9. Instantaneous shear modulus
! 10. Long-term shear modulus
! 11. Shear exponential decay constant
      PROPERTY%Set(61)        = 36
      PROPERTY%Name(61)       = 'KZERO'
      PROPERTY%Location(61)   = 06
      PROPERTY%Set(62)        = 36
      PROPERTY%Name(62)       = 'KINF'
      PROPERTY%Location(62)   = 07
      PROPERTY%Set(63)        = 36
      PROPERTY%Name(63)       = 'KBETA'
      PROPERTY%Location(63)   = 08
      PROPERTY%Set(64)        = 36
      PROPERTY%Name(64)       = 'GZERO'
      PROPERTY%Location(64)   = 09
      PROPERTY%Set(65)        = 36
      PROPERTY%Name(65)       = 'GINF'
      PROPERTY%Location(65)   = 10
      PROPERTY%Set(66)        = 36
      PROPERTY%Name(66)       = 'GBETA'
      PROPERTY%Location(66)   = 11

! MATERIAL_S5 (User defined spring)
      PROPERTY%Set(67)        = 05
      PROPERTY%Name(67)       = 'UK1'
      PROPERTY%Location(67)   = 06
      PROPERTY%Set(68)        = 05
      PROPERTY%Name(68)       = 'UK2'
      PROPERTY%Location(68)   = 07
      PROPERTY%Set(69)        = 05
      PROPERTY%Name(69)       = 'UK3'
      PROPERTY%Location(69)   = 08
      PROPERTY%Set(70)        = 05
      PROPERTY%Name(70)       = 'UK4'
      PROPERTY%Location(70)   = 09
      PROPERTY%Set(71)        = 05
      PROPERTY%Name(71)       = 'UK5'
      PROPERTY%Location(71)   = 10
      PROPERTY%Set(72)        = 05
      PROPERTY%Name(72)       = 'UK6'
      PROPERTY%Location(72)   = 11
      PROPERTY%Set(73)        = 05
      PROPERTY%Name(73)       = 'UK7'
      PROPERTY%Location(73)   = 12
      PROPERTY%Set(74)        = 05
      PROPERTY%Name(74)       = 'UK8'
      PROPERTY%Location(74)   = 13
      PROPERTY%Set(75)        = 05
      PROPERTY%Name(75)       = 'UK9'
      PROPERTY%Location(75)   = 14
      PROPERTY%Set(76)        = 05
      PROPERTY%Name(76)       = 'UK10'
      PROPERTY%Location(76)   = 15
      PROPERTY%Set(77)        = 05
      PROPERTY%Name(77)       = 'UK11'
      PROPERTY%Location(77)   = 16
      PROPERTY%Set(78)        = 05
      PROPERTY%Name(78)       = 'UK12'
      PROPERTY%Location(78)   = 17
      PROPERTY%Set(79)        = 05
      PROPERTY%Name(79)       = 'UK13'
      PROPERTY%Location(79)   = 18
      PROPERTY%Set(80)        = 05
      PROPERTY%Name(80)       = 'UK14'
      PROPERTY%Location(80)   = 19
      PROPERTY%Set(81)        = 05
      PROPERTY%Name(81)       = 'UK15'
      PROPERTY%Location(81)   = 20
      PROPERTY%Set(82)        = 05
      PROPERTY%Name(82)       = 'UK16'
      PROPERTY%Location(82)   = 21
      PROPERTY%Set(83)        = 05
      PROPERTY%Name(83)       = 'UK17'
      PROPERTY%Location(83)   = 22

! MATERIAL_D9 (User defined damper)
      PROPERTY%Set(84)        = 09
      PROPERTY%Name(84)       = 'UD1'
      PROPERTY%Location(84)   = 06
      PROPERTY%Set(85)        = 09
      PROPERTY%Name(85)       = 'UD2'
      PROPERTY%Location(85)   = 07
      PROPERTY%Set(86)        = 09
      PROPERTY%Name(86)       = 'UD3'
      PROPERTY%Location(86)   = 08
      PROPERTY%Set(87)        = 09
      PROPERTY%Name(87)       = 'UD4'
      PROPERTY%Location(87)   = 09
      PROPERTY%Set(88)        = 09
      PROPERTY%Name(88)       = 'UD5'
      PROPERTY%Location(88)   = 10
      PROPERTY%Set(89)        = 09
      PROPERTY%Name(89)       = 'UD6'
      PROPERTY%Location(89)   = 11
      PROPERTY%Set(90)        = 09
      PROPERTY%Name(90)       = 'UD7'
      PROPERTY%Location(90)   = 12
      PROPERTY%Set(91)        = 09
      PROPERTY%Name(91)       = 'UD8'
      PROPERTY%Location(91)   = 13
      PROPERTY%Set(92)        = 09
      PROPERTY%Name(92)       = 'UD9'
      PROPERTY%Location(92)   = 14
      PROPERTY%Set(93)        = 09
      PROPERTY%Name(93)       = 'UD10'
      PROPERTY%Location(93)   = 15
      PROPERTY%Set(94)        = 09
      PROPERTY%Name(94)       = 'UD11'
      PROPERTY%Location(94)   = 16
      PROPERTY%Set(95)        = 09
      PROPERTY%Name(95)       = 'UD12'
      PROPERTY%Location(95)   = 17
      PROPERTY%Set(96)        = 09
      PROPERTY%Name(96)       = 'UD13'
      PROPERTY%Location(96)   = 18
      PROPERTY%Set(97)        = 09
      PROPERTY%Name(97)       = 'UD14'
      PROPERTY%Location(97)   = 19
      PROPERTY%Set(98)        = 09
      PROPERTY%Name(98)       = 'UD15'
      PROPERTY%Location(98)   = 20
      PROPERTY%Set(99)        = 09
      PROPERTY%Name(99)       = 'UD16'
      PROPERTY%Location(99)   = 21
      PROPERTY%Set(100)       = 09
      PROPERTY%Name(100)      = 'UD17'
      PROPERTY%Location(100)  = 22

! MATERIAL_17,27,37 (Rubber: BULK/SHEAR only))
      PROPERTY%Set(101)       = 17
      PROPERTY%Name(101)      = 'BULK'
      PROPERTY%Location(101)  = 10
      PROPERTY%Set(102)       = 17
      PROPERTY%Name(102)      = 'TENSION_TEST'
      PROPERTY%Location(102)  = 09

! MATERIAL_38 (Soil and Crushable Foam)
!  6. Yield function constant
!  7. Yield function constant
!  8. Yield function constant
!  9. Loading shear modulus  (see Mat 33/17)
! 10. Unloading bulk modulus (see Mat 33/17)
! 11. Pressure = f(volume strain) ftn ID
! 12. (2 x SHEAR)               (derived)
! 13. (4 x SHEAR)/3             (derived)
! 14. Fracture pressure         (derived)
      PROPERTY%Set(105)       = 23
      PROPERTY%Name(105)      = 'A0'
      PROPERTY%Location(105)  =  6
      PROPERTY%Set(106)       = 23
      PROPERTY%Name(106)      = 'A1'
      PROPERTY%Location(106)  =  7
      PROPERTY%Set(107)       = 23
      PROPERTY%Name(107)      = 'A2'
      PROPERTY%Location(107)  =  8
      PROPERTY%Set(108)       = 23
      PROPERTY%Name(108)      = 'SHEAR'
      PROPERTY%Location(108)  =  9
      PROPERTY%Set(109)       = 23
      PROPERTY%Name(109)      = 'BULK'
      PROPERTY%Location(109)  = 10
      PROPERTY%Set(110)       = 23
      PROPERTY%Name(110)      = 'PVFTN'
      PROPERTY%Location(110)  = 11
!!
      PROPERTY%Number_of_Entries = NPNMV
!!
      RETURN
      END
!!_
      SUBROUTINE READ_LAYERED_SOLID_LAYUP (FIRST,NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE layering_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      LOGICAL, INTENT(IN) :: FIRST
      INTEGER, INTENT(IN) :: NVALS
!!
!! Local variables.
      INTEGER, SAVE :: NCOUNT
!!
      IF (FIRST) NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMLU) THEN
        LAYERING(NCOUNT)%LupID = I_VALUE(2)
        LAYERING(NCOUNT)%Isys  = I_VALUE(3)
        LMax = MXNLY ! LAYERING(NCOUNT)%Number_of_Layers
        Layer = 0
        i = 3
        DO WHILE (i+13 .LE. NVALS .AND. Layer+1 .LE. LMax)
          Layer = Layer + 1
          LAYERING(NCOUNT)%Number_of_Layers = Layer
          i = i + 1
          LAYERING(NCOUNT)%LayID(Layer) = I_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Ltype(Layer) = I_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%MatID(Layer) = I_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%H(1,Layer)   = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%H(2,Layer)   = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%H(3,Layer)   = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%H(4,Layer)   = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Ax(Layer)    = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Ay(Layer)    = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Az(Layer)    = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Bx(Layer)    = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%By(Layer)    = D_VALUE(i)
          i = i + 1
          LAYERING(NCOUNT)%Bz(Layer)    = D_VALUE(i)
        ENDDO
        IF ((NVALS-3)/13 .GT. LMax) THEN
          WRITE (MSG1,'(I8)') LAYERING(NCOUNT)%LupID
          WRITE (MSG2,'(I8)') LMax
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_LAYERED_SOLID_LAYUP.001.00'//                      &
     &          MSGL//'Layering (LAYUP) Input Record ID:'//MSG1//              &
     &          MSGL//'Specifies More Layers Than Can Be Stored:'//MSG2        &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
      ELSE
!!
!! More lay-up sequences to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_BEAM_SECTION_PROPERTIES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE section_1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument
      INTEGER, INTENT(IN) :: NVALS  ! I/- Number of input record values.
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER :: C_VALUE*32,KEY_WORD*12

      CHARACTER, PARAMETER :: PROFILE(8)*4 =  &
     &  (/'ROD ','TUBE','BAR ','BOX ','I   ','T   ','HAT ','    '/)
!!
!! PROFILE(1) = 'ROD',   Soild  ellipitical/circular cross section
!! PROFILE(2) = 'TUBE',  Hollow ellipitical/circular cross section
!! PROFILE(3) = 'BAR',   Soild  rectangular/square   cross section
!! PROFILE(4) = 'BOX',   Hollow rectangular/square   cross section
!! PROFILE(5) = 'I',     I-Beam cross section
!! PROFILE(6) = 'T',     T-Beam cross section
!! PROFILE(7) = 'HAT',   Hat cross section with effective skin cover
!! PROFILE(8) = ' '      (unused)
!!
      LOGICAL                                                                  &
     &          FOUND
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMS1) THEN
        SECTION_1D(NCOUNT)%SecID = I_VALUE( 2)
!!
!! Convert truss/beam cross section profile name into an integer value.
!!
        FOUND = .FALSE.
        KEY_WORD = C_VALUE(3)
        DO n = 1,5
          IF (PROFILE(n) .EQ. KEY_WORD) THEN
            SECTION_1D(NCOUNT)%Section = n
            FOUND = .TRUE.
          ENDIF
        ENDDO
        IF (.NOT.FOUND) THEN
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_BEAM_SECTION_PROPERTIES'//                         &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        SECTION_1D(NCOUNT)%Iprop     = I_VALUE( 4)
        SECTION_1D(NCOUNT)%NPLoc     = I_VALUE( 5)
        IF (SECTION_1D(NCOUNT)%Iprop .EQ. 0) THEN
          SECTION_1D(NCOUNT)%Width   = D_VALUE( 6)
          SECTION_1D(NCOUNT)%Height  = D_VALUE( 7)
          SECTION_1D(NCOUNT)%Twall   = D_VALUE( 8)
          SECTION_1D(NCOUNT)%Tflange = D_VALUE( 9)
          SECTION_1D(NCOUNT)%Yrefloc = D_VALUE(10)
          SECTION_1D(NCOUNT)%Zrefloc = D_VALUE(11)
        ELSE IF (SECTION_1D(NCOUNT)%Iprop .GT. 0) THEN
          SECTION_1D(NCOUNT)%Area    = D_VALUE( 6)
          SECTION_1D(NCOUNT)%By      = D_VALUE( 7)
          SECTION_1D(NCOUNT)%Bz      = D_VALUE( 8)
          SECTION_1D(NCOUNT)%Yrefloc = D_VALUE( 9)
          SECTION_1D(NCOUNT)%Zrefloc = D_VALUE(10)
        ENDIF
      ELSE
!!
!! More beam/truss section properties to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PLATE_SECTION_PROPERTIES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE section_2d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMS2) THEN
        SECTION_2D(NCOUNT)%SecID     = I_VALUE( 2)
        SECTION_2D(NCOUNT)%Ipts      = I_VALUE( 3)
        SECTION_2D(NCOUNT)%Irule     = I_VALUE( 4)
        SECTION_2D(NCOUNT)%Thickness = D_VALUE( 5)
        SECTION_2D(NCOUNT)%RefLoc    = D_VALUE( 6)
        SECTION_2D(NCOUNT)%Isys      = I_VALUE( 7)
        SECTION_2D(NCOUNT)%Ax        = D_VALUE( 8)
        SECTION_2D(NCOUNT)%Ay        = D_VALUE( 9)
        SECTION_2D(NCOUNT)%Az        = D_VALUE(10)
        SECTION_2D(NCOUNT)%Bx        = D_VALUE(11)
        SECTION_2D(NCOUNT)%By        = D_VALUE(12)
        SECTION_2D(NCOUNT)%Bz        = D_VALUE(13)
      ELSE
!!
!! More plate/membrane section properties to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_RIGID_BODY (NVALS)
!!
!! Copyright (c) by KEY Associates, 10-NOV-1990 13:56:24
!!
      USE shared_common_data
      USE rigid_body_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMRB) THEN
        RIGID_BODY(NCOUNT)%RBID  = I_VALUE( 2) ! Rigid body ID
        RIGID_BODY(NCOUNT)%ParID = I_VALUE( 3) ! RB Part ID
        RIGID_BODY(NCOUNT)%Prop  = I_VALUE( 4) ! Comp'd prop.
        RIGID_BODY(NCOUNT)%CMID  = I_VALUE( 5) ! Rigid mass ID
      ELSE
!!
!! More rigid bodies to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_RB_CONCENTRATED_MASS (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE rigid_body_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMRM) THEN
        RIGID_BODY_MASS(NCOUNT)%RMID     = I_VALUE( 2) ! Rigid mass ID
        RIGID_BODY_MASS(NCOUNT)%Mass     = D_VALUE( 3) ! Mass
        RIGID_BODY_MASS(NCOUNT)%B(1,1)   = D_VALUE( 4) ! Ixx
        RIGID_BODY_MASS(NCOUNT)%B(2,2)   = D_VALUE( 5) ! Iyy
        RIGID_BODY_MASS(NCOUNT)%B(3,3)   = D_VALUE( 6) ! Izz
        RIGID_BODY_MASS(NCOUNT)%B(1,2)   = D_VALUE( 7) ! Ixy
        RIGID_BODY_MASS(NCOUNT)%B(1,3)   = D_VALUE( 8) ! Ixz
        RIGID_BODY_MASS(NCOUNT)%B(2,3)   = D_VALUE( 9) ! Iyz
        RIGID_BODY_MASS(NCOUNT)%NPID     = I_VALUE(10) ! Nodal pt. ID
        RIGID_BODY_MASS(NCOUNT)%Pzero(1) = D_VALUE(11) ! Px of C.M.
        RIGID_BODY_MASS(NCOUNT)%Pzero(2) = D_VALUE(12) ! Py of C.M.
        RIGID_BODY_MASS(NCOUNT)%Pzero(3) = D_VALUE(13) ! Pz of C.M.
        RIGID_BODY_MASS(NCOUNT)%Vel(1)   = D_VALUE(14) ! Vx of C.M.
        RIGID_BODY_MASS(NCOUNT)%Vel(2)   = D_VALUE(15) ! Vy of C.M.
        RIGID_BODY_MASS(NCOUNT)%Vel(3)   = D_VALUE(16) ! Vz of C.M.
        RIGID_BODY_MASS(NCOUNT)%Omega(1) = D_VALUE(17) ! Ox of C.M.
        RIGID_BODY_MASS(NCOUNT)%Omega(2) = D_VALUE(18) ! Oy of C.M.
        RIGID_BODY_MASS(NCOUNT)%Omega(3) = D_VALUE(19) ! Oz of C.M.
      ELSE
!!
!! More rigid body substitute/added mass records to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_NP_CONCENTRATED_MASS (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE nodal_point_mass_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMCM) THEN
        NODAL_POINT_MASS(NCOUNT)%NMID   = I_VALUE( 2) ! Con'd mass ID
        NODAL_POINT_MASS(NCOUNT)%Mass   = D_VALUE( 3) ! Mass
        NODAL_POINT_MASS(NCOUNT)%B(1,1) = D_VALUE( 4) ! Ixx
        NODAL_POINT_MASS(NCOUNT)%B(2,2) = D_VALUE( 5) ! Iyy
        NODAL_POINT_MASS(NCOUNT)%B(3,3) = D_VALUE( 6) ! Izz
        NODAL_POINT_MASS(NCOUNT)%B(1,2) = D_VALUE( 7) ! Ixy
        NODAL_POINT_MASS(NCOUNT)%B(1,3) = D_VALUE( 8) ! Ixz
        NODAL_POINT_MASS(NCOUNT)%B(2,3) = D_VALUE( 9) ! Iyz
        NODAL_POINT_MASS(NCOUNT)%NPID   = I_VALUE(10) ! Nodal point ID
      ELSE
!!
!! More concentrated nodal point masses to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_DISPLACEMENT_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE displacement_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMDC) THEN
        DISPLACEMENT_BC(NCOUNT)%DBCID = I_VALUE( 2) ! Disp. BC ID
        DISPLACEMENT_BC(NCOUNT)%SetID = I_VALUE( 3) ! Nodal set ID
        DISPLACEMENT_BC(NCOUNT)%Code  = I_VALUE( 4) ! Constraint
        DISPLACEMENT_BC(NCOUNT)%HstID = I_VALUE( 5) ! History ftn.
        DISPLACEMENT_BC(NCOUNT)%Kavd  = I_VALUE( 6) ! (1/2/3=a/v/d)
        DISPLACEMENT_BC(NCOUNT)%Scale = D_VALUE( 7) ! Scale factor
        DISPLACEMENT_BC(NCOUNT)%Ax    = D_VALUE( 8) ! x-direction
        DISPLACEMENT_BC(NCOUNT)%Ay    = D_VALUE( 9) ! y-direction
        DISPLACEMENT_BC(NCOUNT)%Az    = D_VALUE(10) ! z-direction
      ELSE
!!
!! More displacement BC's to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_SPOT_WELD (NVALS)
!!
!! Copyright (c) by KEY Associates; 14-JUN-1994 21:44:42.00
!!
      USE shared_common_data
      USE spot_weld_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER :: C_VALUE*32,KEY_WORD*12
!!
!! LOCATION(1) = 'NODES',     Nodal points spot welded together.
!! LOCATION(2) = 'POSITION',  Position     spot welded together.
!!
      CHARACTER :: LOCATION(2)*12 = (/'NODES   ','POSITION'/)

      LOGICAL                                                                  &
     &          FOUND
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSW) THEN
        SPOT_WELD(NCOUNT)%SWID = I_VALUE( 2)   ! Spot weld ID
!!
!! Distinguish between spot welds connecting coincident nodal points and
!! spot welds located at a user specified coordinate position.
!!
        FOUND = .FALSE.
        KEY_WORD = C_VALUE(3)
        IF (KEY_WORD .EQ. LOCATION(1)) THEN
          FOUND = .TRUE.
          SPOT_WELD(NCOUNT)%PLACE   = LOCATION(1)         ! Nodal point based
          SPOT_WELD(NCOUNT)%NPID(1) = I_VALUE( 4) ! Nodal point ID 1
          SPOT_WELD(NCOUNT)%NPID(2) = I_VALUE( 5) ! Nodal point ID 2
          SPOT_WELD(NCOUNT)%NPID(3) = I_VALUE( 6) ! Nodal point ID 3
          SPOT_WELD(NCOUNT)%Fmax    = D_VALUE( 7) ! Failure force

        ELSE IF (KEY_WORD .EQ. LOCATION(2)) THEN
          FOUND = .TRUE.
          SPOT_WELD(NCOUNT)%PLACE   = LOCATION(2)         ! Position based
          SPOT_WELD(NCOUNT)%Px      = D_VALUE( 4) ! X-coordinate
          SPOT_WELD(NCOUNT)%Py      = D_VALUE( 5) ! Y-coordinate
          SPOT_WELD(NCOUNT)%Pz      = D_VALUE( 6) ! Z-coordinate
          SPOT_WELD(NCOUNT)%Fmax    = D_VALUE( 7) ! Failure force

        ENDIF
        IF (.NOT.FOUND) THEN
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_SPOT_WELD'//                                       &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
      ELSE
!!
!! More spot welds to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_TIED_BC (NVALS)
!!
!! Copyright (c) by KEY Associates;  6-JUL-1992 11:58:15.77
!!
      USE shared_common_data
      USE tied_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMTC) THEN
        TIED_BC(NCOUNT)%TBCID = I_VALUE( 2) ! Tied BC ID
        TIED_BC(NCOUNT)%SetID = I_VALUE( 3) ! Nodal set ID
        TIED_BC(NCOUNT)%Code  = I_VALUE( 4) ! Constraint
        TIED_BC(NCOUNT)%Fmax  = D_VALUE( 5) ! Max force
        TIED_BC(NCOUNT)%Ax    = D_VALUE( 6) ! x-direction
        TIED_BC(NCOUNT)%Ay    = D_VALUE( 7) ! y-direction
        TIED_BC(NCOUNT)%Az    = D_VALUE( 8) ! z-direction
      ELSE
!!
!! More tied BC's to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_RIGID_WALL_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE rigid_wall_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMWC) THEN
        RIGID_WALL_BC(NCOUNT)%RWID   = I_VALUE( 2) ! Rigid wall ID
        RIGID_WALL_BC(NCOUNT)%SetID  = I_VALUE( 3) ! Nodal set  ID
        RIGID_WALL_BC(NCOUNT)%P1(1)  = D_VALUE( 4) ! X of point P1
        RIGID_WALL_BC(NCOUNT)%P1(2)  = D_VALUE( 5) ! Y of point P1
        RIGID_WALL_BC(NCOUNT)%P1(3)  = D_VALUE( 6) ! Z of point P1
        RIGID_WALL_BC(NCOUNT)%P2(1)  = D_VALUE( 7) ! X of point P2
        RIGID_WALL_BC(NCOUNT)%P2(2)  = D_VALUE( 8) ! Y of point P2
        RIGID_WALL_BC(NCOUNT)%P2(3)  = D_VALUE( 9) ! Z of point P2
        RIGID_WALL_BC(NCOUNT)%P3(1)  = D_VALUE(10) ! X of point P3
        RIGID_WALL_BC(NCOUNT)%P3(2)  = D_VALUE(11) ! Y of point P3
        RIGID_WALL_BC(NCOUNT)%P3(3)  = D_VALUE(12) ! Z of point P3
        RIGID_WALL_BC(NCOUNT)%Width  = D_VALUE(13) ! 12-direction
        RIGID_WALL_BC(NCOUNT)%Height = D_VALUE(14) ! 13-direction
        RIGID_WALL_BC(NCOUNT)%CoF    = D_VALUE(15) ! Friction coef.
        RIGID_WALL_BC(NCOUNT)%Kode   = I_VALUE(16) ! rigid/inert/move
        IF (RIGID_WALL_BC(NCOUNT)%Kode .GT. 0) THEN
          RIGID_WALL_BC(NCOUNT)%Mass   = D_VALUE(17) ! mass at P1
          RIGID_WALL_BC(NCOUNT)%B(1,1) = D_VALUE(18) ! Ixx at P1
          RIGID_WALL_BC(NCOUNT)%B(2,2) = D_VALUE(19) ! Iyy at P1
          RIGID_WALL_BC(NCOUNT)%B(3,3) = D_VALUE(20) ! Izz at P1
          RIGID_WALL_BC(NCOUNT)%B(1,2) = D_VALUE(21) ! Ixy at P1
          RIGID_WALL_BC(NCOUNT)%B(1,3) = D_VALUE(22) ! Ixz at P1
          RIGID_WALL_BC(NCOUNT)%B(2,3) = D_VALUE(23) ! Iyz at P1
          RIGID_WALL_BC(NCOUNT)%Code   = I_VALUE(24) ! Disp. BC on P1
          IF (RIGID_WALL_BC(NCOUNT)%Code .GT. 0) THEN
            RIGID_WALL_BC(NCOUNT)%HstID  = I_VALUE(25) ! History ID
            RIGID_WALL_BC(NCOUNT)%Kavd   = I_VALUE(26) ! (1/2/3=a/v/d)
            RIGID_WALL_BC(NCOUNT)%Scale  = D_VALUE(27) ! scale factor
            RIGID_WALL_BC(NCOUNT)%Ax     = D_VALUE(28) ! BC x-direct.
            RIGID_WALL_BC(NCOUNT)%Ay     = D_VALUE(29) ! BC y-direct.
            RIGID_WALL_BC(NCOUNT)%Az     = D_VALUE(30) ! BC z-direct.
          ENDIF
        ENDIF
      ELSE
!!
!! More rigid wall conditions to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_BODY_FORCE (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE body_force_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMBF) THEN
        BODY_FORCE(NCOUNT)%BFID    = I_VALUE( 2) ! Body Force ID
        BODY_FORCE(NCOUNT)%SetID   = I_VALUE( 3) ! Nodal point Set ID
        BODY_FORCE(NCOUNT)%HstID   = I_VALUE( 4) ! History ID
        BODY_FORCE(NCOUNT)%Scale   = D_VALUE( 5) ! Scale factor
        BODY_FORCE(NCOUNT)%Gravity = D_VALUE( 6) ! Gravity constant
        BODY_FORCE(NCOUNT)%Ax      = D_VALUE( 7) ! X-direction
        BODY_FORCE(NCOUNT)%Ay      = D_VALUE( 8) ! Y-direction
        BODY_FORCE(NCOUNT)%Az      = D_VALUE( 9) ! Z-direction
        BODY_FORCE(NCOUNT)%Delay   = D_VALUE(10) ! Arrival time
      ELSE
!!
!! More body forces to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PRESSURE_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE pressure_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMPC) THEN
        PRESSURE_BC(NCOUNT)%PBCID = I_VALUE( 2) ! Pressure BC ID
        PRESSURE_BC(NCOUNT)%SetID = I_VALUE( 3) ! Seg Set ID
        PRESSURE_BC(NCOUNT)%HstID = I_VALUE( 4) ! History ID
        PRESSURE_BC(NCOUNT)%Scale = D_VALUE( 5) ! Scale factor
        PRESSURE_BC(NCOUNT)%PI    = D_VALUE( 6) ! Multiplier at I
        PRESSURE_BC(NCOUNT)%PJ    = D_VALUE( 7) ! Multiplier at J
        PRESSURE_BC(NCOUNT)%PK    = D_VALUE( 8) ! Multiplier at K
        PRESSURE_BC(NCOUNT)%PL    = D_VALUE( 9) ! Multiplier at L
        PRESSURE_BC(NCOUNT)%Delay = D_VALUE(10) ! Arrival time
      ELSE
!!
!! More pressure boundary conditions to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_FORCE_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE force_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMFC) THEN
        FORCE_BC(NCOUNT)%CFID   = I_VALUE( 2) ! Con. force ID
        FORCE_BC(NCOUNT)%SetID  = I_VALUE( 3) ! NP set ID
        FORCE_BC(NCOUNT)%HstID  = I_VALUE( 4) ! History ID
        FORCE_BC(NCOUNT)%Type   = I_VALUE( 5) ! 0/1=force/torque
        FORCE_BC(NCOUNT)%Follow = I_VALUE( 6) ! 0/1=no/yes
        FORCE_BC(NCOUNT)%Force  = D_VALUE( 7) ! Force/torque mag.
        FORCE_BC(NCOUNT)%Cx     = D_VALUE( 8) ! X-direction
        FORCE_BC(NCOUNT)%Cy     = D_VALUE( 9) ! Y-direction
        FORCE_BC(NCOUNT)%Cz     = D_VALUE(10) ! Z-direction
        FORCE_BC(NCOUNT)%Delay  = D_VALUE(11) ! Delay time
!!
      ELSE
!!
!! More concentrated force conditions to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_SPRING_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE spring_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSC) THEN
        SPRING_BC(NCOUNT)%SprID     = I_VALUE( 2) ! Spring ID
        SPRING_BC(NCOUNT)%SetID     = I_VALUE( 3) ! Nodal ID
        SPRING_BC(NCOUNT)%MatID     = I_VALUE( 4) ! Material ID
        SPRING_BC(NCOUNT)%Type      = I_VALUE( 5) ! 0/1=axial/tors'l
        SPRING_BC(NCOUNT)%Follow    = I_VALUE( 6) ! 0/1/2=U/V/Axis
        SPRING_BC(NCOUNT)%Axis(1)   = D_VALUE( 7) ! x-direct/axis
        SPRING_BC(NCOUNT)%Axis(2)   = D_VALUE( 8) ! y-direct/axis
        SPRING_BC(NCOUNT)%Axis(3)   = D_VALUE( 9) ! z-direct/axis
        SPRING_BC(NCOUNT)%RES%Force = D_VALUE(10) ! Initial force
      ELSE
!!
!! More spring restraints to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_DAMPER_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE damper_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMVC) THEN
        DAMPER_BC(NCOUNT)%DprID   = I_VALUE(2) ! Damper ID
        DAMPER_BC(NCOUNT)%SetID   = I_VALUE(3) ! Nodal ID
        DAMPER_BC(NCOUNT)%MatID   = I_VALUE(4) ! Material ID
        DAMPER_BC(NCOUNT)%Type    = I_VALUE(5) ! 0/1=axial/tors'l
        DAMPER_BC(NCOUNT)%Follow  = I_VALUE(6) ! 0/1/2=U/V/Axis
        DAMPER_BC(NCOUNT)%Axis(1) = D_VALUE(7) ! x-direct/axis
        DAMPER_BC(NCOUNT)%Axis(2) = D_VALUE(8) ! y-direct/axis
        DAMPER_BC(NCOUNT)%Axis(3) = D_VALUE(9) ! z-direct/axis
      ELSE
!!
!! More DAMPER restraints to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PERIODIC_BC (NVALS)
!!
!! Copyright (c) by KEY Associates; 17-APR-1994 17:01:40.26
!!
      USE shared_common_data
      USE periodic_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMCC) THEN
        PERIODIC_BC(NCOUNT)%PerID = I_VALUE(2)         ! Periodic BC ID
!!
!! Read side 1 definition.
!!
        KEY_WORD = C_VALUE(3)                          ! Set Type
        IF (INDEX(KEY_WORD,'SEGSET') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Typ1 = 0                 ! Segment
        ELSE IF (INDEX(KEY_WORD,'NPSET') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Typ1 = 1                 ! Node
        ELSE
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_PERIODIC_BC.001.01'//                              &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        PERIODIC_BC(NCOUNT)%S1ID = I_VALUE(4)          ! Set ID
!!
!! Read side 2 definition.
!!
        KEY_WORD = C_VALUE(5)                          ! Set Type
        IF (INDEX(KEY_WORD,'SEGSET') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Typ2 = 0                 ! Segment
        ELSE IF (INDEX(KEY_WORD,'NPSET') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Typ2 = 1                 ! Node
        ELSE
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_PERIODIC_BC.001.02'//                              &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        PERIODIC_BC(NCOUNT)%S2ID = I_VALUE(6)          ! Set ID
!!
!! Read form of periodic boundary condition.
!!
        KEY_WORD = C_VALUE(7)                          ! Geometry
        IF (INDEX(KEY_WORD,'LINEAR') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Type = 'LINEAR'          ! Linear
        ELSE IF (INDEX(KEY_WORD,'CYCLIC') .NE. 0) THEN
          PERIODIC_BC(NCOUNT)%Type = 'CYCLIC'          ! Cyclic
        ELSE
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_PERIODIC_BC.001.03'//                              &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
!! Read other parameters for sliding interface calculations.
!!
        PERIODIC_BC(NCOUNT)%Axis(1)   = D_VALUE( 8) ! x-direct/axis
        PERIODIC_BC(NCOUNT)%Axis(2)   = D_VALUE( 9) ! y-direct/axis
        PERIODIC_BC(NCOUNT)%Axis(3)   = D_VALUE(10) ! z-direct/axis
        PERIODIC_BC(NCOUNT)%Origin(1) = D_VALUE(11) ! x-position
        PERIODIC_BC(NCOUNT)%Origin(2) = D_VALUE(12) ! y-position
        PERIODIC_BC(NCOUNT)%Origin(3) = D_VALUE(13) ! z-position
        PERIODIC_BC(NCOUNT)%Theta     = D_VALUE(14) ! Angle, degrees
        PERIODIC_BC(NCOUNT)%Advance   = D_VALUE(15) ! Helixcal adv./deg.
      ELSE
!!
!! More PERIODBC restraints to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_NONREFLECTING_BC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE nonreflecting_bc_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMNR) THEN
        NONREFLECTING_BC(NCOUNT)%NRID  = I_VALUE(2) ! NRBC ID
        NONREFLECTING_BC(NCOUNT)%SetID = I_VALUE(3) ! Seg Set ID
      ELSE
!!
!! More non-reflecting boundary conditions to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_SLIDING_INTERFACE (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE sliding_interface_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER                                                                &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12
      EXTERNAL                                                                 &
     &          C_VALUE
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSI) THEN
        SLIDING_INTERFACE(NCOUNT)%SIID = I_VALUE(2)  ! Interface ID
!!
!! Read side 1 definition.
!!
        KEY_WORD = C_VALUE(3)                        ! Set Type
        IF (INDEX(KEY_WORD,'SEGSET') .NE. 0) THEN
          SLIDING_INTERFACE(NCOUNT)%Typ1 = 0                 ! Segment
        ELSE IF (INDEX(KEY_WORD,'NPSET') .NE. 0) THEN
          SLIDING_INTERFACE(NCOUNT)%Typ1 = 1                 ! Node
        ELSE
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_SLIDING_INTERFACE.001.01'//                        &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        SLIDING_INTERFACE(NCOUNT)%S1ID = I_VALUE( 4) ! Set ID
!!
!! Read side 2 definition.
!!
        KEY_WORD = C_VALUE(5)                        ! Set Type
        IF (INDEX(KEY_WORD,'SEGSET') .NE. 0) THEN
          SLIDING_INTERFACE(NCOUNT)%Typ2 = 0                 ! Segment
        ELSE IF (INDEX(KEY_WORD,'NPSET') .NE. 0) THEN
          SLIDING_INTERFACE(NCOUNT)%Typ2 = 1                 ! Node
        ELSE
!!
!! Unexpected key word found.
!!
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_SLIDING_INTERFACE.001.02'//                        &
     &          MSGL//'Unexpected Keyword Found: '//KEY_WORD                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
        SLIDING_INTERFACE(NCOUNT)%S2ID = I_VALUE( 6) ! Set ID
!!
!! Read other parameters for sliding interface calculations.
!!
        SLIDING_INTERFACE(NCOUNT)%Type      = I_VALUE( 7) ! n=0 only
        SLIDING_INTERFACE(NCOUNT)%Isym      = I_VALUE( 8) ! 0/1/2=sy/ms/sm
        SLIDING_INTERFACE(NCOUNT)%CoF       = D_VALUE( 9) ! Fric Coef.
        SLIDING_INTERFACE(NCOUNT)%Begin     = D_VALUE(10) ! Start time
        SLIDING_INTERFACE(NCOUNT)%End       = D_VALUE(11) ! End time
        SLIDING_INTERFACE(NCOUNT)%Factor    = D_VALUE(12) ! Relaxation
        SLIDING_INTERFACE(NCOUNT)%Capture   = D_VALUE(13) ! Thickness
        SLIDING_INTERFACE(NCOUNT)%Border    = D_VALUE(14) ! El border
!!        SLIDING_INTERFACE(NCOUNT)%Sort_Freq = I_VALUE(15) ! Resort Frq
      ELSE
!!
!! More sliding interfaces to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_TABULATED_FUNCTION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE tabulated_function_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMTF) THEN
        TABULATED_FUNCTION(NCOUNT)%TFID = I_VALUE(2) ! Tab ftn ID
        Nmax = TABULATED_FUNCTION(NCOUNT)%Number_of_Pairs
        n = 0
        i = 2
        DO WHILE (i+2 .LE. NVALS .AND. n+1 .LE. Nmax)
          n = n + 1
          TABULATED_FUNCTION(NCOUNT)%Number_of_Pairs = n
          i = i + 1
          TABULATED_FUNCTION(NCOUNT)%X(n) = D_VALUE(i) ! x-value
          i = i + 1
          TABULATED_FUNCTION(NCOUNT)%Y(n) = D_VALUE(i) ! y-value
        ENDDO
      ELSE
!!
!! More tabulated functions to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE LOG_NODE_SET_LOCATION (IFINDEX,ILINDEX)
!!
!! Copyright (c) by KEY Associates, 10-DEC-1991 20:29:39
!!
      USE shared_common_data
      USE node_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: IFINDEX  ! I/S Points to INCLUDE_FILE(IFINDEX)
      INTEGER, INTENT(IN) :: ILINDEX  ! S/S Start-of-Record line in input file
!!
!! Local variable.
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMNS) THEN
        NODE_SET(NCOUNT)%File_Number = IFINDEX
        NODE_SET(NCOUNT)%Line_Number = ILINDEX
      ELSE
!!
!! More node sets to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE LOG_ELEMENT_SET_LOCATION (IFINDEX,ILINDEX)
!!
!! Copyright (c) by KEY Associates, 10-DEC-1991 20:29:39
!!
      USE shared_common_data
      USE element_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: IFINDEX  ! I/S Points to INCLUDE_FILE(IFINDEX)
      INTEGER, INTENT(IN) :: ILINDEX  ! S/S Start-of-Record line in input file
!!
!! Local variable.
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMES) THEN
        ELEMENT_SET(NCOUNT)%File_Number = IFINDEX
        ELEMENT_SET(NCOUNT)%Line_Number = ILINDEX
      ELSE
!!
!! More element sets to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE LOG_SEGMENT_SET_LOCATION (IFINDEX,ILINDEX)
!!
!! Copyright (c) by KEY Associates, 10-DEC-1991 20:29:39
!!
      USE shared_common_data
      USE segment_set_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: IFINDEX  ! I/S Points to INCLUDE_FILE(IFINDEX)
      INTEGER, INTENT(IN) :: ILINDEX  ! S/S Start-of-Record line in input file
!!
!! Local variable.
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSS) THEN
        SEGMENT_SET(NCOUNT)%File_Number = IFINDEX
        SEGMENT_SET(NCOUNT)%Line_Number = ILINDEX
      ELSE
!!
!! More segment sets to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_AND_BUILD_NODE_SETS (ACTION)
!!
!! Copyright (c) by KEY Associates, 12-DEC-1991 14:03:15
!! Copyright (c) by KEY Associates; 11-NOV-1996 22:51:33
!!
!! Purpose: Read node set input records and build node set member lists.
!! 
      USE shared_common_data
      USE value_
      USE node_
      USE node_set_
      USE location_
      USE include_file_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      CHARACTER, INTENT(IN) :: ACTION*(*)  ! I/- First 'COUNT' then 'BUILD'
!!
      CHARACTER                                                                &
     &          TEXT*80,                                                       &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      INTEGER   :: ILCOUNT      ! S/S Input file line counter
      INTEGER   :: ILINDEX      ! S/- Start-of-Record line number in input file

      LOGICAL                                                                  &
     &          EOF,                                                           &
     &          FOUND,                                                         &
     &          IOERROR,                                                       &
     &          GENERATE,                                                      &
     &          Nbgn_STORED,                                                   &
     &          Nend_STORED
!!
!! Read mesh data if it exists.
!!
      IF (CONTROL%RDMESH .NE. 0 .AND. MSHNE .GT. 0) THEN
        IF (ACTION .EQ. 'BUILD') THEN
          DO i = 1,MSHNE
            READ (IO_UNIT%LMDI,*) NNPSETS(i)
          ENDDO
          LOCATION%Next_Node = MSHNE + 1
        ELSE
          LOCATION%Next_Node = MSHNE + 1
        ENDIF
      ELSE
        LOCATION%Next_Node = 1
      ENDIF
!!
!! Read NODE_SET records only and build NNPSETS as required. This module
!! expects to process an NPSET input record in one of the following forms:
!!      NPSET 100 Nodes_For_DISPBC_100 = ALL
!!      NPSET 100 Nodes_For_DISPBC_100 = 1 2 3 4 ...
!!      NPSET 100 Nodes_For_DISPBC_100 = 1 2 3 THRU 732 806 ...
!!      NPSET 100 Nodes_For_DISPBC_100 = THRU 123
!! Only those nodal point ID's found in the list of nodal ID's
!! {NODE(1:NUMNP)%ID} will be used.
!!
      IERR = 0
      IF_OPEN = 0
      ILCOUNT = 0
      EOF = .FALSE.
      ERROR%COUNT = 0
      DO 300 N = MSHNS+1,NUMNS
        IFN = NODE_SET(N)%File_Number
!!
!! Open the appropriate included file as dictated by the file number.
!!
        IF (IF_OPEN .NE. IFN) THEN
!!
!! Close old included file.
!!
          IF (IF_OPEN .NE. 0) CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(IFN)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200      CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_NODE_SETS.001.01'//                      &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of read of included input file for nodal point set record.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_AND_BUILD_NODE_SETS.001.02'//                      &
     &          MSGL//'INCLUDE File Pre-Read For  NPSET Records: '             &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            IF_OPEN = IFN
          ENDIF
        ENDIF
!!
!! Position input file to record for node set N.
!!
        REWIND (IO_UNIT%LSDI)
        ILINDEX = NODE_SET(N)%Line_Number
        DO i = 1,ILINDEX - 1
          READ (IO_UNIT%LSDI,100)
 100      FORMAT ()
        ENDDO
!!
!! Read record for node set N.
!!
        NVALS = MAXRE
        CALL GETIRV                                                            &
     &  (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
        KEY_WORD = C_VALUE(1)
!!
        IF (INDEX(KEY_WORD,'NPSET') .NE. 0) THEN
!!
!! Store set-ID and set-label.
!!
          NODE_SET(N)%SetID = I_VALUE(2)
          NODE_SET(N)%Label = C_VALUE(3)
          Ibgn = 4
!!
!! Check to insure that there is at least ONE entry enumerating
!! the nodal ID's comprising this set. (Must check all remaining
!! entries since interspersed null values are permitted.) We will
!! also find the first non-null entry in the node ID list.
!!
          i = Ibgn
          Jbgn = Ibgn
          FOUND = .FALSE.
          DO WHILE (.NOT.FOUND .AND. i.LE.NVALS)
            FOUND = FOUND .OR. (VALUE(i)%VTYP .NE. 'N')
            IF (FOUND) Jbgn = i
            i = i + 1
          ENDDO

          IF (.NOT.FOUND) THEN
            WRITE (MSG1,'(I8)') NODE_SET(N)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_NODE_SETS.003.01'//                      &
     &          MSGL//'NPSET (Nodal Point Set) Input Record ID:'//MSG1//       &
     &          MSGL//'Has A Null Nodal Point ID List.'                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
            GO TO 300
          ENDIF
!!
!! If the first entry in the ID list is a character constant, check
!! to see if it is one the two valid values: "ALL" or "THRU"
!!
          IF (VALUE(Jbgn)%VTYP .EQ. 'C') THEN
            KEY_WORD = C_VALUE(Jbgn)

            IF ((INDEX(KEY_WORD, 'ALL') .EQ. 0) .AND.                          &
     &            (INDEX(KEY_WORD,'THRU') .EQ. 0)) THEN

              WRITE (MSG1,'(I8)') NODE_SET(N)%SetID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_NODE_SETS.003.02'//                      &
     &          MSGL//'NPSET (Nodal Point Set) Input Record ID:'//MSG1//       &
     &          MSGL//'Has An Unrecognized Key Word:'//C_VALUE(Jbgn)//         &
     &          MSGL//'At The Start Of The Nodal Point ID List.'               &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
              GO TO 300
            ENDIF
          ELSE
            KEY_WORD = ' '
          ENDIF
!!
!! We can be moderately certain that the ID list is valid. If the next
!! entry is the key-word "ALL," define flag-word to be "ALL." We are
!! finished; no ID entries in NNPSETS are required for this case.
!!
          IF (INDEX(KEY_WORD,'ALL') .NE. 0) THEN
            NODE_SET(N)%Istart = 0
            NODE_SET(N)%Iend   = 0
            NODE_SET(N)%Flag   = 'ALL '
!!
!! There are integers and THRU's to process.
!!
          ELSE
            NODE_SET(N)%Istart = LOCATION%Next_Node
            Nbgn = 1
            GENERATE = .FALSE.
            Nbgn_STORED = .FALSE.
            DO i = Jbgn,NVALS
!!
!! Check for key-word "THRU"
!!
              IF (VALUE(i)%VTYP .EQ. 'C') THEN
                KEY_WORD = C_VALUE(i)
                IF (INDEX(KEY_WORD,'THRU') .NE. 0) THEN
                  GENERATE = .TRUE.
                ELSE
                  GENERATE = .FALSE.
                  WRITE (MSG1,'(I8)') NODE_SET(N)%SetID
                  CALL USER_MESSAGE                                            &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'READ_AND_BUILD_NODE_SETS.003.03'//                    &
     &            MSGL//'NPSET (Nodal Point Set) Input Record ID:'             &
     &                //MSG1//                                                 &
     &            MSGL//'Key-Word Expected: THRU'//                            &
     &            MSGL//'Key-Word Found: '//KEY_WORD//                         &
     &            MSGL//'No Generation Will Take Place.'                       &
     &            )
                  ERROR%COUNT = ERROR%COUNT + 1
                ENDIF
              ELSE IF (VALUE(i)%VTYP .NE. 'N') THEN
!!
!! Generation flag was NOT set on last trip through loop.
!!
                IF (.NOT.GENERATE) THEN
                  Next = I_VALUE(i)
                  Nbgn = Next
                  FOUND = .FALSE.
                  DO M = 1,NUMNP
                    FOUND = FOUND .OR. (Next .EQ. NODE(M)%ID)
                  ENDDO
                  IF (FOUND) THEN
                    IF (ACTION .EQ. 'BUILD') THEN
                      Nbgn_STORED = .TRUE.
                      NNPSETS(LOCATION%Next_Node) = Next
                      LOCATION%Next_Node = LOCATION%Next_Node + 1
                    ELSE
                      Nbgn_STORED = .TRUE.
                      LOCATION%Next_Node = LOCATION%Next_Node + 1
                    ENDIF
                  ELSE
                    Nbgn_STORED = .FALSE.
                  ENDIF
!!
!! Generation flag was set on last trip through loop.
!!
                ELSE
                  Nend = I_VALUE(i)
                  Nlow = MIN (Nbgn,Nend)
                  Nhgh = MAX (Nbgn,Nend)
                  IF (Nbgn_STORED) THEN
                    IF (Nbgn .EQ. Nlow) THEN
                      Nlow = Nlow + 1
                    ELSE
                      Nhgh = Nhgh - 1
                    ENDIF
                  ENDIF
                  Nend_STORED = .FALSE.
                  DO Next = Nlow,Nhgh
                    FOUND = .FALSE.
                    DO M = 1,NUMNP
                      FOUND = FOUND .OR. (Next .EQ. NODE(M)%ID)
                    ENDDO
                    IF (FOUND) THEN
                      IF (ACTION .EQ. 'BUILD') THEN
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        NNPSETS(LOCATION%Next_Node) = Next
                        LOCATION%Next_Node = LOCATION%Next_Node + 1
                      ELSE
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        LOCATION%Next_Node = LOCATION%Next_Node + 1
                      ENDIF
                    ENDIF
                  ENDDO
                  Nbgn = Nend  !  (An odd but acceptable situation)
                  Nbgn_STORED = Nend_STORED
                  GENERATE = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            NODE_SET(N)%Iend = LOCATION%Next_Node - 1
            NODE_SET(N)%Flag = 'E   '
          ENDIF
        ELSE
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_AND_BUILD_NODE_SETS.004.00'//                      &
     &          MSGL//'NPSET Input Record Expected.'//                         &
     &          MSGL//'Input Record Key-Word Found: '//KEY_WORD//              &
     &          MSGL//'Logic Error. Call KEY Associates.'                      &
     &          )
        ENDIF
!!
 300    ENDDO
!!
!! Define number of entries in nodal-point-set member list.
!!
      NUMNE = LOCATION%Next_Node - 1
!!
!! Close "last" included file.
!!
      IF (IF_OPEN .NE. 0) THEN
        CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
        IF_OPEN = 0
      ENDIF
!!
!! Print total number of nodal point set record errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &  (                                                                      &
     &  MSGL//'FATAL'//                                                        &
     &  MSGL//'READ_AND_BUILD_NODE_SETS.005.00'//                              &
     &  MSGL//'Total Number Of Open & Format Errors In Input:'//MSG1//         &
     &  MSGL//'Execution Terminated By Program.'                               &
     &  )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_AND_BUILD_ELEMENT_SETS (ACTION)
!!
!! Copyright (c) by KEY Associates, 12-DEC-1991 14:03:15
!!
!! Purpose: Read element set input records and build element set member lists.
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
      USE value_
      USE element_set_
      USE location_
      USE include_file_
      USE enumerated_sets_, ONLY: NELSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      CHARACTER, INTENT(IN) :: ACTION*(*)  ! I/- First 'COUNT' then 'BUILD'
!!
      CHARACTER                                                                &
     &          TEXT*80,                                                       &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      INTEGER   :: ILCOUNT      ! S/S Input file line counter
      INTEGER   :: ILINDEX      ! S/- Start-of-Record line in input file

      LOGICAL                                                                  &
     &          EOF,                                                           &
     &          FOUND,                                                         &
     &          IOERROR,                                                       &
     &          GENERATE,                                                      &
     &          Nbgn_STORED,                                                   &
     &          Nend_STORED
!!
!! Read mesh data if it exists.
!!
      IF (CONTROL%RDMESH .NE. 0 .AND. MSHEE .GT. 0) THEN
        IF (ACTION .EQ. 'BUILD') THEN
          DO i = 1,MSHEE
            READ (IO_UNIT%LMDI,*) NELSETS(i)
          ENDDO
          LOCATION%Next_Element = MSHEE + 1
        ELSE
          LOCATION%Next_Element = MSHEE + 1
        ENDIF
      ELSE
        LOCATION%Next_Element = 1
      ENDIF
!!
!! Read ELEMENT_SET records only and build NELSETS as required. This module
!! expects to process an ELSET input record in one of the following forms:
!!      ELSET 100 = ALL
!!      ELSET 100 = 1 2 3 4 ...
!!      ELSET 100 = 1 2 3 THRU 732 806 ...
!!      ELSET 100 = THRU 123 ...
!! When "THRU" occurs, only those nodal point ID's found in the elements will
!! be used.
!!
      IERR = 0
      IF_OPEN = 0
      ILCOUNT = 0
      EOF = .FALSE.
      ERROR%COUNT = 0
      DO 300 N = MSHES+1,NUMES
        IFN = ELEMENT_SET(N)%File_Number
!!
!! Open the the appropriate included file as dictated by the file number.
!!
        IF (IF_OPEN .NE. IFN) THEN
!!
!! Close old included file.
!!
          IF (IF_OPEN .NE. 0) CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(IFN)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200        CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_ELEMENT_SETS.001.01'//                   &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of pre-read of included input record file.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_AND_BUILD_ELEMENT_SETS.001.02'//                   &
     &          MSGL//'INCLUDE File Pre-Read For  ELSET Records: '             &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            IF_OPEN = IFN
          ENDIF
        ENDIF
!!
!! Position input file to element set record N.
!!
        REWIND (IO_UNIT%LSDI)
        ILINDEX = ELEMENT_SET(N)%Line_Number
        DO i = 1,ILINDEX - 1
          READ (IO_UNIT%LSDI,100)
 100      FORMAT ()
        ENDDO
!!
!! Read element set record N.
!!
        NVALS = MAXRE
        CALL GETIRV                                                            &
     &  (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
        KEY_WORD = C_VALUE(1)
!!
        IF (INDEX(KEY_WORD,'ELSET') .NE. 0) THEN
!!
!! Store set-ID and set-label.
!!
          ELEMENT_SET(N)%SetID = I_VALUE(2)
          ELEMENT_SET(N)%Label = C_VALUE(3)
          Ibgn = 4
!!
!! Check to insure that there is at least ONE entry enumerating
!! the element ID's comprising this set. (Must check all remaining
!! entries since interspersed null values are permitted.) We will
!! also find the first non-null entry in the element ID list.
!!
          i = Ibgn
          Jbgn = Ibgn
          FOUND = .FALSE.
          DO WHILE (.NOT.FOUND .AND. i.LE.NVALS)
            FOUND = FOUND .OR. (VALUE(i)%VTYP .NE. 'N')
            IF (FOUND) Jbgn = i
            i = i + 1
          ENDDO

          IF (.NOT.FOUND) THEN
            WRITE (MSG1,'(I8)') ELEMENT_SET(N)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_ELEMENT_SETS.003.01'//                   &
     &          MSGL//'ELSET (Element Set) Input Record ID:'//MSG1//           &
     &          MSGL//'Has A Null Element ID List.'                            &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
            GO TO 300
          ENDIF
!!
!! If the first entry in the ID list is a character constant, check
!! to see if it is one the two valid values: "ALL" or "THRU"
!!
          IF (VALUE(Jbgn)%VTYP .EQ. 'C') THEN
            KEY_WORD = C_VALUE(Jbgn)

            IF ((INDEX(KEY_WORD, 'ALL') .EQ. 0) .AND.                          &
     &            (INDEX(KEY_WORD,'THRU') .EQ. 0)) THEN

              WRITE (MSG1,'(I8)') ELEMENT_SET(N)%SetID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_ELEMENT_SETS.003.02'//                   &
     &          MSGL//'ELSET (Element Set) Input Record ID:'//MSG1//           &
     &          MSGL//'Has An Unrecognized Key Word:'//C_VALUE(Jbgn)//         &
     &          MSGL//'At The Start Of The Element ID List.'                   &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
              GO TO 300
            ENDIF
          ELSE
            KEY_WORD = ' '
          ENDIF
!!
!! We can be moderately certain that the ID list is valid. If the next
!! entry is the key-word "ALL," define flag-word to be "ALL." We are
!! finished; no ID entries in NELSETS are required for this case.
!!
          IF (INDEX(KEY_WORD,'ALL') .NE. 0) THEN
            ELEMENT_SET(N)%Istart = 0
            ELEMENT_SET(N)%Iend   = 0
            ELEMENT_SET(N)%Flag   = 'ALL '
!!
!! There are integers and THRU's to process.
!!
          ELSE
            ELEMENT_SET(N)%Istart = LOCATION%Next_Element
            Nbgn = 1
            GENERATE = .FALSE.
            Nbgn_STORED = .FALSE.
            DO i = Jbgn,NVALS
!!
!! Check for key-word "THRU"
!!
              IF (VALUE(i)%VTYP .EQ. 'C') THEN
                KEY_WORD = C_VALUE(i)
                IF (INDEX(KEY_WORD,'THRU') .NE. 0) THEN
                  GENERATE = .TRUE.
                ELSE
                  GENERATE = .FALSE.
                  WRITE (MSG1,'(I8)') ELEMENT_SET(N)%SetID
                  CALL USER_MESSAGE                                            &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'READ_AND_BUILD_ELEMENT_SETS.003.03'//                 &
     &            MSGL//'ELSET (Element Set) Input Record ID:'//MSG1//         &
     &            MSGL//'Key-Word Expected: THRU'//                            &
     &            MSGL//'Key-Word Found: '//KEY_WORD//                         &
     &            MSGL//'No Generation Will Take Place.'                       &
     &            )
                  ERROR%COUNT = ERROR%COUNT + 1
                ENDIF
              ELSE IF (VALUE(i)%VTYP .NE. 'N') THEN
!!
!! Case (a): Generation flag was NOT set on last trip through loop.
!!
                IF (.NOT.GENERATE) THEN
                  Next = I_VALUE(i)
                  Nbgn = Next
                  FOUND = .FALSE.
                  DO M = 1,NUMHX
                    FOUND = FOUND .OR. Next .EQ. HEXAH(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMPX
                    FOUND = FOUND .OR. Next .EQ. PENTA(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMTX
                    FOUND = FOUND .OR. Next .EQ. TETRA(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMLS
                    FOUND = FOUND .OR. Next .EQ. LSOLD(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMM3
                    FOUND = FOUND .OR. Next .EQ. MEMBT(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMM4
                    FOUND = FOUND .OR. Next .EQ. MEMBQ(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMTR
                    FOUND = FOUND .OR. Next .EQ. TRUSS(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMP3
                    FOUND = FOUND .OR. Next .EQ. PLATT(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMP4
                    FOUND = FOUND .OR. Next .EQ. PLATQ(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMBM
                    FOUND = FOUND .OR. Next .EQ. BEAM(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMSP
                    FOUND = FOUND .OR. Next .EQ. SPRING(M)%PAR%EleID
                  ENDDO
                  DO M = 1,NUMDM
                    FOUND = FOUND .OR. Next .EQ. DAMPER(M)%PAR%EleID
                  ENDDO
                  IF (FOUND) THEN
                    IF (ACTION .EQ. 'BUILD') THEN
                      Nbgn_STORED = .TRUE.
                      NELSETS(LOCATION%Next_Element) = Next
                      LOCATION%Next_Element =                                  &
     &                  LOCATION%Next_Element + 1
                    ELSE
                      Nbgn_STORED = .TRUE.
                      LOCATION%Next_Element =                                  &
     &                  LOCATION%Next_Element + 1
                    ENDIF
                  ELSE
                    Nbgn_STORED = .FALSE.
                  ENDIF
!!
!! Case (b): Generation flag was set on last trip through loop.
!!
                ELSE
                  Nend = I_VALUE(i)
                  Nlow = MIN (Nbgn,Nend)
                  Nhgh = MAX (Nbgn,Nend)
                  IF (Nbgn_STORED) THEN
                    IF (Nbgn .EQ. Nlow) THEN
                      Nlow = Nlow + 1
                    ELSE
                      Nhgh = Nhgh - 1
                    ENDIF
                  ENDIF
                  Nend_STORED = .FALSE.
                  DO Next = Nlow,Nhgh
                    FOUND = .FALSE.
                    DO M = 1,NUMHX
                      FOUND = FOUND .OR. Next .EQ. HEXAH(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMPX
                      FOUND = FOUND .OR. Next .EQ. PENTA(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMTX
                      FOUND = FOUND .OR. Next .EQ. TETRA(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMLS
                      FOUND = FOUND .OR. Next .EQ. LSOLD(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMM3
                      FOUND = FOUND .OR. Next .EQ. MEMBT(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMM4
                      FOUND = FOUND .OR. Next .EQ. MEMBQ(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMTR
                      FOUND = FOUND .OR. Next .EQ. TRUSS(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMP3
                      FOUND = FOUND .OR. Next .EQ. PLATT(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMP4
                      FOUND = FOUND .OR. Next .EQ. PLATQ(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMBM
                      FOUND = FOUND .OR. Next .EQ. BEAM(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMSP
                      FOUND = FOUND .OR. Next .EQ. SPRING(M)%PAR%EleID
                    ENDDO
                    DO M = 1,NUMDM
                      FOUND = FOUND .OR. Next .EQ. DAMPER(M)%PAR%EleID
                    ENDDO
                    IF (FOUND) THEN
                      IF (ACTION .EQ. 'BUILD') THEN
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        NELSETS(LOCATION%Next_Element) = Next
                        LOCATION%Next_Element = LOCATION%Next_Element+1
                      ELSE
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        LOCATION%Next_Element = LOCATION%Next_Element+1
                      ENDIF
                    ENDIF
                  ENDDO
                  Nbgn = Nend  !  (An odd but acceptable situation)
                  Nbgn_STORED = Nend_STORED
                  GENERATE = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            ELEMENT_SET(N)%Iend = LOCATION%Next_Element - 1
            ELEMENT_SET(N)%Flag = 'E   '
          ENDIF
        ELSE
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_AND_BUILD_ELEMENT_SETS.001.00'//                   &
     &          MSGL//'ELSET Input Record Expected.'//                         &
     &          MSGL//'Input Record Key-Word Found: '//KEY_WORD//              &
     &          MSGL//'Logic Error. Call KEY Associates.'                      &
     &          )
        ENDIF
!!
 300    ENDDO
!!
!! Define number of entries in element-set member list.
!!
      NUMEE = LOCATION%Next_Element - 1
!!
!! Close "last" included file.
!!
      IF (IF_OPEN .NE. 0) THEN
        CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
        IF_OPEN = 0
      ENDIF
!!
!! Print total number of element set record errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'READ_AND_BUILD_ELEMENT_SETS.004.00'//                         &
     &    MSGL//'Total Number Of Open & Format Errors In Input:'//MSG1//       &
     &    MSGL//'Execution Terminated By Program.'                             &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_AND_BUILD_SEGMENT_SETS (ACTION)
!!
!! Copyright (c) by KEY Associates, 24-JAN-1992 12:01:26
!! Copyright (c) by KEY Associates; 11-NOV-1996 22:51:33
!!
!! Purpose: Read segment set input records and build segment set member lists.
!!
      USE shared_common_data
      USE value_
      USE segment_
      USE segment_set_
      USE location_
      USE include_file_
      USE enumerated_sets_, ONLY: NSGSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      CHARACTER, INTENT(IN) :: ACTION*(*)  ! I/- First 'COUNT' then 'BUILD'
!!
      CHARACTER                                                                &
     &          TEXT*80,                                                       &
     &          C_VALUE*32,                                                    &
     &          KEY_WORD*12

      INTEGER   :: ILCOUNT      ! S/S Input file line counter
      INTEGER   :: ILINDEX      ! S/- Start-of-Record line number in input file

      LOGICAL                                                                  &
     &          EOF,                                                           &
     &          FOUND,                                                         &
     &          IOERROR,                                                       &
     &          GENERATE,                                                      &
     &          Nbgn_STORED,                                                   &
     &          Nend_STORED
!!
!! Read mesh data if it exists.
!!
      IF (CONTROL%RDMESH .NE. 0 .AND. MSHSE .GT. 0) THEN
        IF (ACTION .EQ. 'BUILD') THEN
          DO i = 1,MSHSE
            READ (IO_UNIT%LMDI,*) NSGSETS(i)
          ENDDO
          LOCATION%Next_Segment = MSHSE + 1
        ELSE
          LOCATION%Next_Segment = MSHSE + 1
        ENDIF
      ELSE
        LOCATION%Next_Segment = 1
      ENDIF
!!
!! Read SEGMENT_SET records only and build NSGSETS as required. This module
!! expects to process an SEGSET input record in one of the following forms:
!!      SEGSET 100 Segments_For_PRESSBC_100 = ALL
!!      SEGSET 100 Segments_For_PRESSBC_100 = 1 2 3 4 ...
!!      SEGSET 100 Segments_For_PRESSBC_100 = 1 2 3 THRU 732 806 ...
!!      SEGSET 100 Segments_For_PRESSBC_100 = THRU 123
!! Only those segment ID's found in the list of segment ID's
!! {SEGMENT(1:NUMSG)%PAR%SegID} will be used.
!!
      IERR = 0
      IF_OPEN = 0
      ILCOUNT = 0
      EOF = .FALSE.
      ERROR%COUNT = 0
      DO 300 N = MSHSS+1,NUMSS
        IFN = SEGMENT_SET(N)%File_Number
!!
!! Open the appropriate included file as dictated by the file number.
!!
        IF (IF_OPEN .NE. IFN) THEN
!!
!! Close old included file.
!!
          IF (IF_OPEN .NE. 0) CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
!!
!! Open the required included file.
!!
          IOERROR = .TRUE.
          OPEN                                                                 &
     &          (                                                              &
     &          UNIT   =  IO_UNIT%LSDI,                                        &
     &          FILE   =  INCLUDE_FILE(IFN)%Full_Name,                         &
     &          STATUS = 'OLD',                                                &
     &          FORM   = 'FORMATTED',                                          &
     &          ERR    =  200                                                  &
     &          )
          IOERROR = .FALSE.
 200        CONTINUE
!!
!! Warning exit for failed OPEN operation.
!!
          IF (IOERROR) THEN
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_SEGMENT_SETS.001.01'//                   &
     &          MSGL//'Unable To Execute OPEN On: '                            &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
!!
!! Inform user of read of included input file for segment set record.
!!
          ELSE
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'READ_AND_BUILD_SEGMENT_SETS.001.02'//                   &
     &          MSGL//'INCLUDE File Pre-Read For SEGSET Records: '             &
     &              //TRIM(INCLUDE_FILE(IFN)%Full_Name)                        &
     &          )
            IF_OPEN = IFN
          ENDIF
        ENDIF
!!
!! Position input file to record for segment set N.
!!
        REWIND (IO_UNIT%LSDI)
        ILINDEX = SEGMENT_SET(N)%Line_Number
        DO i = 1,ILINDEX - 1
          READ (IO_UNIT%LSDI,100)
 100      FORMAT ()
        ENDDO
!!
!! Read record for segment set N.
!!
        NVALS = MAXRE
        CALL GETIRV                                                            &
     &  (NVALS,IO_UNIT%LSDI,IO_UNIT%LELO,EOF,IERR,TEXT,ILINDEX,ILCOUNT)
        KEY_WORD = C_VALUE(1)
!!
        IF (INDEX(KEY_WORD,'SEGSET') .NE. 0) THEN
!!
!! Store set-ID and set-label.
!!
          SEGMENT_SET(N)%SetID = I_VALUE(2)
          SEGMENT_SET(N)%Label = C_VALUE(3)
          Ibgn = 4
!!
!! Check to insure that there is at least ONE entry enumerating
!! the segment ID's comprising this set. (Must check all remaining
!! entries since interspersed null values are permitted.) We will
!! also find the first non-null entry in the segment ID list.
!!
          i = Ibgn
          Jbgn = Ibgn
          FOUND = .FALSE.
          DO WHILE (.NOT.FOUND .AND. i.LE.NVALS)
            FOUND = FOUND .OR. (VALUE(i)%VTYP .NE. 'N')
            IF (FOUND) Jbgn = i
            i = i + 1
          ENDDO

          IF (.NOT.FOUND) THEN
            WRITE (MSG1,'(I8)') SEGMENT_SET(N)%SetID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_SEGMENT_SETS.003.01'//                   &
     &          MSGL//'SEGSET (Segment Set) Input Record ID:'//MSG1//          &
     &          MSGL//'Has A Null Segment ID List.'                            &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
            GO TO 300
          ENDIF
!!
!! If the first entry in the ID list is a character constant, check
!! to see if it is one the two valid values: "ALL" or "THRU"
!!
          IF (VALUE(Jbgn)%VTYP .EQ. 'C') THEN
            KEY_WORD = C_VALUE(Jbgn)

            IF ((INDEX(KEY_WORD, 'ALL') .EQ. 0) .AND.                          &
     &            (INDEX(KEY_WORD,'THRU') .EQ. 0)) THEN

              WRITE (MSG1,'(I8)') SEGMENT_SET(N)%SetID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'READ_AND_BUILD_SEGMENT_SETS.003.02'//                   &
     &          MSGL//'SEGSET (Segment Set) Input Record ID:'//MSG1//          &
     &          MSGL//'Has An Unrecognized Key Word:'//C_VALUE(Jbgn)//         &
     &          MSGL//'At The Start Of The Segment ID List.'                   &
     &          )
              ERROR%COUNT = ERROR%COUNT + 1
              GO TO 300
            ENDIF
          ELSE
            KEY_WORD = ' '
          ENDIF
!!
!! We can be moderately certain that the ID list is valid. If the next
!! entry is the key-word "ALL," define flag-word to be "ALL." We are
!! finished; no ID entries in NSGSETS are required for this case.
!!
          IF (INDEX(KEY_WORD,'ALL') .NE. 0) THEN
            SEGMENT_SET(N)%Istart = 0
            SEGMENT_SET(N)%Iend   = 0
            SEGMENT_SET(N)%Flag   = 'ALL '
!!
!! There are integers and THRU's to process.
!!
          ELSE
            SEGMENT_SET(N)%Istart = LOCATION%Next_Segment
            Nbgn = 1
            GENERATE = .FALSE.
            Nbgn_STORED = .FALSE.
            DO i = Jbgn,NVALS
!!
!! Check for key-word "THRU"
!!
              IF (VALUE(i)%VTYP .EQ. 'C') THEN
                KEY_WORD = C_VALUE(i)
                IF (INDEX(KEY_WORD,'THRU') .NE. 0) THEN
                  GENERATE = .TRUE.
                ELSE
                  GENERATE = .FALSE.
                  WRITE (MSG1,'(I8)') SEGMENT_SET(N)%SetID
                  CALL USER_MESSAGE                                            &
     &            (                                                            &
     &            MSGL//'WARN'//                                               &
     &            MSGL//'READ_AND_BUILD_SEGMENT_SETS.003.03'//                 &
     &            MSGL//'SEGSET (Segment Set) Input Record ID:'//MSG1//        &
     &            MSGL//'Key-Word Expected: THRU'//                            &
     &            MSGL//'Key-Word Found: '//KEY_WORD//                         &
     &            MSGL//'No Generation Will Take Place.'                       &
     &            )
                  ERROR%COUNT = ERROR%COUNT + 1
                ENDIF
              ELSE IF (VALUE(i)%VTYP .NE. 'N') THEN
!!
!! Generation flag was NOT set on last trip through loop.
!!
                IF (.NOT.GENERATE) THEN
                  Next = I_VALUE(i)
                  Nbgn = Next
                  FOUND = .FALSE.
                  DO M = 1,NUMSG
                  FOUND = FOUND .OR. (Next .EQ. SEGMENT(M)%PAR%SegID)
                  ENDDO
                  IF (FOUND) THEN
                    IF (ACTION .EQ. 'BUILD') THEN
                      Nbgn_STORED = .TRUE.
                      NSGSETS(LOCATION%Next_Segment) = Next
                      LOCATION%Next_Segment = LOCATION%Next_Segment+1
                    ELSE
                      Nbgn_STORED = .TRUE.
                      LOCATION%Next_Segment = LOCATION%Next_Segment+1
                    ENDIF
                  ELSE
                    Nbgn_STORED = .FALSE.
                  ENDIF
!!
!! Generation flag was set on last trip through loop.
!!
                ELSE
                  Nend = I_VALUE(i)
                  Nlow = MIN (Nbgn,Nend)
                  Nhgh = MAX (Nbgn,Nend)
                  IF (Nbgn_STORED) THEN
                    IF (Nbgn .EQ. Nlow) THEN
                      Nlow = Nlow + 1
                    ELSE
                      Nhgh = Nhgh - 1
                    ENDIF
                  ENDIF
                  Nend_STORED = .FALSE.
                  DO Next = Nlow,Nhgh
                    FOUND = .FALSE.
                    DO M = 1,NUMSG
                      FOUND = FOUND.OR.(Next.EQ.SEGMENT(M)%PAR%SegID)
                    ENDDO
                    IF (FOUND) THEN
                      IF (ACTION .EQ. 'BUILD') THEN
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        NSGSETS(LOCATION%Next_Segment) = Next
                        LOCATION%Next_Segment = LOCATION%Next_Segment+1
                      ELSE
                        Nend_STORED = Nend_STORED .OR. (Next .EQ. Nend)
                        LOCATION%Next_Segment = LOCATION%Next_Segment+1
                      ENDIF
                    ENDIF
                  ENDDO
                  Nbgn = Nend  !  (An odd but acceptable situation)
                  Nbgn_STORED = Nend_STORED
                  GENERATE = .FALSE.
                ENDIF
              ENDIF
            ENDDO
            SEGMENT_SET(N)%Iend = LOCATION%Next_Segment - 1
            SEGMENT_SET(N)%Flag = 'E   '
          ENDIF
        ELSE
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'READ_AND_BUILD_SEGMENT_SETS.004.00'//                   &
     &          MSGL//'SEGSET Input Record Expected.'//                        &
     &          MSGL//'Input Record Key-Word Found: '//KEY_WORD//              &
     &          MSGL//'Logic Error. Call KEY Associates.'                      &
     &          )
        ENDIF
!!
 300    ENDDO
!!
!! Define number of entries in segment-set member list.
!!
      NUMSE = LOCATION%Next_Segment - 1
!!
!! Close "last" included file.
!!
      IF (IF_OPEN .NE. 0) THEN
        CLOSE (UNIT=IO_UNIT%LSDI,STATUS='KEEP')
        IF_OPEN = 0
      ENDIF
!!
!! Print total number of segment set record errors detected.
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'READ_AND_BUILD_SEGMENT_SETS.005.00'//                         &
     &    MSGL//'Total Number Of Open & Format Errors In Input:'//MSG1//       &
     &    MSGL//'Execution Terminated By Program.'                             &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_NODAL_POINT_COORDINATES (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE node_
      USE motion_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMNP) THEN
        NODE(NCOUNT)%ID   = I_VALUE(2) ! Nodal Point ID
        MOTION(NCOUNT)%Px = D_VALUE(3) ! X-position
        MOTION(NCOUNT)%Py = D_VALUE(4) ! Y-position
        MOTION(NCOUNT)%Pz = D_VALUE(5) ! Z-position
      ELSE
!!
!! More nodal point records to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_HXEL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE hexah_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMHX) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        HEXAH(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        HEXAH(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        HEXAH(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        HEXAH(NCOUNT)%PAR%LupID = I_VALUE(5) ! Lay-up ID
        DO i = 1,8
          HEXAH(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More hexahedron elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_PXEL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE penta_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMPX) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        PENTA(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        PENTA(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        PENTA(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        PENTA(NCOUNT)%PAR%LupID = I_VALUE(5) ! Lay-up ID
        DO i = 1,6
          PENTA(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More pentahedron elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_TXEL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE tetra_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMTX) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        TETRA(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        TETRA(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        TETRA(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        TETRA(NCOUNT)%PAR%LupID = I_VALUE(5) ! Lay-up ID
        DO i = 1,4
          TETRA(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More tetrahedron elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_LSEL_DEFINITION (FIRST,NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE lsold_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      LOGICAL, INTENT(IN) :: FIRST
      INTEGER, INTENT(IN) :: NVALS
!!
!! Local variables.
      INTEGER, SAVE :: NCOUNT
!!
      IF (FIRST) NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMLS) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        LSOLD(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        LSOLD(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        LSOLD(NCOUNT)%PAR%LupID = I_VALUE(4) ! Lay-up ID
        DO i = 1,8
          LSOLD(NCOUNT)%PAR%IX(i) = I_VALUE(i+4) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More layered solid elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_M4EL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE membq_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMM4) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        MEMBQ(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        MEMBQ(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        MEMBQ(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        MEMBQ(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,4
          MEMBQ(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More 4-node membrane elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_M3EL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE membt_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMM3) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        MEMBT(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        MEMBT(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        MEMBT(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        MEMBT(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,3
          MEMBT(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More 3-node membrane elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_TRUSS_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE truss_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMTR) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        TRUSS(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        TRUSS(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        TRUSS(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        TRUSS(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,2
          TRUSS(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More 2-node axial force elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_P4EL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE platq_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMP4) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        PLATQ(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        PLATQ(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        PLATQ(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        PLATQ(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,4
          PLATQ(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More 4-node plate elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_P3EL_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE platt_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMP3) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        PLATT(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        PLATT(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        PLATT(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        PLATT(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,3
          PLATT(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! DOF Indicies
        ENDDO
      ELSE
!!
!! More 3-node plate elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_BEAM_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE beam_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMBM) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        BEAM(NCOUNT)%PAR%EleID = I_VALUE(2) ! Element ID
        BEAM(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        BEAM(NCOUNT)%PAR%MatID = I_VALUE(4) ! Material ID
        BEAM(NCOUNT)%PAR%SecID = I_VALUE(5) ! Section ID
        DO i = 1,2
          BEAM(NCOUNT)%PAR%IX(i) = I_VALUE(i+5) ! Nodal point ID's
        ENDDO
        DO i = 1,3
          BEAM(NCOUNT)%RES%Zaxis(i) = D_VALUE(i+7) ! Sec.Orientation
        ENDDO
      ELSE
!!
!! More 2-node beam elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_SPRING_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE spring_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSP) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        SPRING(NCOUNT)%PAR%EleID   = I_VALUE( 2) ! Element ID
        SPRING(NCOUNT)%PAR%ParID   = I_VALUE( 3) ! Part ID
        SPRING(NCOUNT)%PAR%MatID   = I_VALUE( 4) ! Material ID
        SPRING(NCOUNT)%PAR%Type    = I_VALUE( 5) ! 0/1=axial/torsional
        SPRING(NCOUNT)%PAR%IX(1)   = I_VALUE( 6) ! DOF Indicies
        SPRING(NCOUNT)%PAR%IX(2)   = I_VALUE( 7) ! DOF Indicies
        SPRING(NCOUNT)%PAR%Idir    = I_VALUE( 8) ! 0/1/2/3/4=i,j/x/...
        SPRING(NCOUNT)%PAR%Axis(1) = D_VALUE( 9) ! Ax, (Idir=4)
        SPRING(NCOUNT)%PAR%Axis(2) = D_VALUE(10) ! Ay, (Idir=4)
        SPRING(NCOUNT)%PAR%Axis(3) = D_VALUE(11) ! Az, (Idir=4)
      ELSE
!!
!! More 2-node spring elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_DAMPER_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE damper_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMDM) THEN
        MXEID = MAX (MXEID,I_VALUE(2))
        DAMPER(NCOUNT)%PAR%EleID   = I_VALUE( 2) ! Element ID
        DAMPER(NCOUNT)%PAR%ParID   = I_VALUE( 3) ! Part ID
        DAMPER(NCOUNT)%PAR%MatID   = I_VALUE( 4) ! Material ID
        DAMPER(NCOUNT)%PAR%Type    = I_VALUE( 5) ! 0/1=axial/torsional
        DAMPER(NCOUNT)%PAR%IX(1)   = I_VALUE( 6) ! DOF Indicies
        DAMPER(NCOUNT)%PAR%IX(2)   = I_VALUE( 7) ! DOF Indicies
        DAMPER(NCOUNT)%PAR%Idir    = I_VALUE( 8) ! 0/1/2/3/4=i,j/x/...
        DAMPER(NCOUNT)%PAR%Axis(1) = D_VALUE( 9) ! Ax, (Idir=4)
        DAMPER(NCOUNT)%PAR%Axis(2) = D_VALUE(10) ! Ay, (Idir=4)
        DAMPER(NCOUNT)%PAR%Axis(3) = D_VALUE(11) ! Az, (Idir=4)
      ELSE
!!
!! More 2-node damper elements to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_SEGMENT_DEFINITION (NVALS)
!!
!! Copyright (c) by KEY Associates, 16-DEC-1990 15:31:20
!!
      USE shared_common_data
      USE segment_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMSG) THEN
        SEGMENT(NCOUNT)%PAR%SegID = I_VALUE(2) ! Segment ID
        SEGMENT(NCOUNT)%PAR%ParID = I_VALUE(3) ! Part ID
        SEGMENT(NCOUNT)%PAR%IX(1) = I_VALUE(4) ! DOF Indicies
        SEGMENT(NCOUNT)%PAR%IX(2) = I_VALUE(5) ! DOF Indicies
        SEGMENT(NCOUNT)%PAR%IX(3) = I_VALUE(6) ! DOF Indicies
        SEGMENT(NCOUNT)%PAR%IX(4) = I_VALUE(7) ! DOF Indicies
      ELSE
!!
!! More segments to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_VELOCITY_IC (NVALS)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE velocity_ic_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMIC) THEN
        VELOCITY_IC(NCOUNT)%ICID = I_VALUE(2) ! Velocity IC ID
        VELOCITY_IC(NCOUNT)%NPID = I_VALUE(3) ! Nodal point ID
        VELOCITY_IC(NCOUNT)%Vx   = D_VALUE(4) ! Translational
        VELOCITY_IC(NCOUNT)%Vy   = D_VALUE(5) ! Translational
        VELOCITY_IC(NCOUNT)%Vz   = D_VALUE(6) ! Translational
        VELOCITY_IC(NCOUNT)%Ox   = D_VALUE(7) ! Rotational
        VELOCITY_IC(NCOUNT)%Oy   = D_VALUE(8) ! Rotational
        VELOCITY_IC(NCOUNT)%Oz   = D_VALUE(9) ! Rotational
      ELSE
!!
!! More velocity IC's to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_NODAL_CONSTRAINT_1 (NVALS)
!!
!! Copyright (c) by KEY Associates; 16-OCT-1995 21:05:45.00
!!
      USE shared_common_data
      USE constrained_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMNC) THEN
        CONSTRAINED_NODE(NCOUNT)%ID      = I_VALUE(2) ! Constraint ID
        CONSTRAINED_NODE(NCOUNT)%CNID    = I_VALUE(3) ! Constrained Node ID
        CONSTRAINED_NODE(NCOUNT)%NPID(1) = I_VALUE(4) ! 1st Node ID
        CONSTRAINED_NODE(NCOUNT)%NPID(2) = I_VALUE(5) ! 2nd Node ID
      ELSE
!!
!! More constrained node's to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_QA_RECORD (NVALS,TEXT)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE value_
      USE qa_record_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER, SAVE :: NCOUNT = 0
!!
      CHARACTER                                                                &
     &          TEXT*80
!!
      NCOUNT = NCOUNT + 1
      IF (NCOUNT .LE. NUMQA) THEN
!!
        QA_RECORD(NCOUNT)%LINE = ADJUSTL(TEXT(VALUE(2)%LOC:))
      ELSE
!!
!! More user QA records to read than counted.
!!
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE READ_JOB_IDENTIFICATION (NVALS,TEXT)
!!
!! Copyright (c) by KEY Associates, 17-JUL-1990 19:39:54
!!
      USE shared_common_data
      USE value_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          TEXT*80
!!
      JOB_ID_RECORD%CURRENT%TITLE = TEXT(VALUE(2)%LOC:)
      CALL LEFT_JUSTIFY                                                        &
     &          (                                                              &
     &          JOB_ID_RECORD%CURRENT%TITLE,                                   &
     &          JOB_ID_RECORD%TITLE_LENGTH                                     &
     &          )
!!
      RETURN
      END
!!_
      SUBROUTINE CONNECT_LAYERED_SOLIDS
!!
!! Copyright (c) by KEY Associates, 6-JAN-1991 09:31:39
!!
!! Purpose: Connect up layered solid elements with hexahedrons
!! and membranes which will represent the respective layers.
!!
      USE shared_common_data
      USE layering_
      USE lsold_
      USE hexah_, ONLY: LSHEX
      USE membq_, ONLY: LSMBQ
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          INTERNAL_ID,                                                   &
     &          LY_INTERNAL_ID
      EXTERNAL                                                                 &
     &          LY_INTERNAL_ID
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered CONNECT_LAYERED_SOLIDS.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Define auxillary hexa's and membranes for use with layered solids.
!!
      NXTHX = 0
      NXTM4 = 0
      DO n = 1,NUMLS
        INTERNAL_ID = LY_INTERNAL_ID (LSOLD(n)%PAR%LupID)
        DO i = 1,LAYERING(INTERNAL_ID)%Number_of_Layers
          IF (LAYERING(INTERNAL_ID)%Ltype(i) .EQ. 0) THEN
            NXTHX = NXTHX + 1
            MXEID = MXEID + 1
            LSOLD(n)%PAR%ID(i) = MXEID
            LSHEX(NXTHX)%PAR%EleID = MXEID
            LSHEX(NXTHX)%PAR%ParID = LSOLD(n)%PAR%ParID
            LSHEX(NXTHX)%PAR%MatID = LAYERING(INTERNAL_ID)%MatID(i)
          ELSE
            NXTM4 = NXTM4 + 1
            MXEID = MXEID + 1
            LSOLD(n)%PAR%ID(i) = MXEID
            LSMBQ(NXTM4)%PAR%EleID = MXEID
            LSMBQ(NXTM4)%PAR%ParID = LSOLD(n)%PAR%ParID
            LSMBQ(NXTM4)%PAR%MatID = LAYERING(INTERNAL_ID)%MatID(i)
            LSMBQ(NXTM4)%PAR%SecID = LAYERING(INTERNAL_ID)%Ltype(i)
          ENDIF
          DO k = 1,4
            LSOLD(n)%RES%H(k,i) = LAYERING(INTERNAL_ID)%H(k,i)
          ENDDO
        ENDDO
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_WORKING_STORAGE
!!
!! Copyright (c) by KEY Associates, 27-MAR-1992 16:30:09
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
      USE nrbc_data_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE plate_pair_
!!
      USE enumerated_sets_
      USE stress_
      USE state_variables_
      USE coord_
      USE indx_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                    &
     &          SegID,           &
     &          SetID,           &
     &          ParID,           &
     &          HXface(4,6),     & ! -/- Index array to obtain hexahedron  face
     &          PXface(4,5),     & ! -/- Index array to obtain pentahedron face
     &          TXface(3,4),     & ! -/- Index array to obtain tetrahedron face
     &          NEIGHBOR(10)       ! -/- Local work array used to build ICE
      REAL(KIND(0D0))            &
     &          Xce(4),          & ! -/- Contact element nodal coordinates     
     &          Yce(4),          & ! -/- Contact element nodal coordinates     
     &          Zce(4)             ! -/- Contact element nodal coordinates
      REAL(KIND(0D0))            &
     &          DTR,             & ! -/- Coverts degrees into radians          
     &          Theta,           & ! -/- Angle of rotation in periodic BC      
     &          CosThem1,        & ! -/- Cosine(Theta) minus 1.0               
     &          SinTheta,        & ! -/- Sine(Theta)                           
     &          C(3),            & ! -/- Vector in periodic BC                 
     &          RTX(3,3),        & ! -/- Rotation in periodic BC               
     &          PRX(3,3)           ! -/- Scratch in periodic BC
      LOGICAL                                                                  &
     &          MATCH,                                                         &
     &          FOUND,                                                         &
     &          COINCIDENT,                                                    &
     &          NEXT_EL_ID,                                                    &
     &          NEXT_NP_ID,                                                    &
     &          NEXT_SEG_ID
      DATA                                                                     &
     &          HXface /1,4,3,2,5,6,7,8,3,4,8,7,2,3,7,6,1,2,6,5,1,5,8,4/       &
     &          PXface /1,2,4,5,2,3,6,5,1,3,6,4,1,3,2,0,4,5,6,0/               &
     &          TXface /1,3,2,1,2,4,1,4,3,2,3,4/
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered INITIALIZE_WORKING_STORAGE.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
      PI = ATAN2 ( -1.0D+0 , 0.0D+0 )
      DTR = PI / 180.D+0
!!
      ERROR%COUNT = 0
!!
!! ### ZERO OUT STRESS AND STATE VARIABLE ARRAYS.
!!
      DO N = 1,NUMST
        STRESS(1,N) = 0.0
        STRESS(2,N) = 0.0
        STRESS(3,N) = 0.0
        STRESS(4,N) = 0.0
        STRESS(5,N) = 0.0
        STRESS(6,N) = 0.0
      ENDDO
      DO N = 1,NUMAX
        STATE_VARIABLES(N) = 0.0
      ENDDO
!!
!! ### DEFINE ENDING SIMULATION TIME.
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
!! ### DEFINE STATE VARIABLE STORAGE.
!! Initialize material data structures with state variable storage count.
!!
      DO N = 1,NUMMT
        MATERIAL(N)%Nsv = MAUX(MATERIAL(N)%Type)
      ENDDO
!!
!! ### DEFINE PROPERTY SET POINTER.
!! Initialize material data structures with property set pointer.
!!
      DO N = 1,NUMMT
        MATERIAL(N)%Set = MSET(MATERIAL(N)%Type)
      ENDDO
!!
!! ### DEFINE STATE VARIABLES ARRAY POINTERS.
!! Isv is a pointer containing the starting location, element by element, of
!! the material state variable storage. Note: shell elements contain several
!! integration stations through the thickness and beam elements use up to 16
!! points over the cross-section at which to evaluate the stress-strain model.
!!
      NV = 1
      NS = 1
!!
      DO N = 1,NUMHX
        HEXAH(N)%PAR%Isv = NV
        NV = NV + MATERIAL(HEXAH(N)%PAR%MatID)%Nsv
      ENDDO
      DO N = 1,NUMPX
        PENTA(N)%PAR%Isv = NV
        NV = NV + MATERIAL(PENTA(N)%PAR%MatID)%Nsv
      ENDDO
      DO N = 1,NUMTX
        TETRA(N)%PAR%Isv = NV
        NV = NV + MATERIAL(TETRA(N)%PAR%MatID)%Nsv
      ENDDO
!!
      DO N = 1,NUMLS
        LupID = LSOLD(N)%PAR%LupID
        DO i = 1,LAYERING(LupID)%Number_of_Layers
          IF (LAYERING(LupID)%Ltype(i) .EQ. 0) THEN
            LSHEX(LSOLD(N)%PAR%ID(i))%PAR%Isv = NV
          ELSE
            LSMBQ(LSOLD(N)%PAR%ID(i))%PAR%Isv = NV
          ENDIF
          NV = NV + MATERIAL(LAYERING(LupID)%MatID(i))%Nsv
        ENDDO
      ENDDO
!!
      DO N = 1,NUMM3
        MEMBT(N)%PAR%Isv = NV
        NV = NV + MATERIAL(MEMBT(N)%PAR%MatID)%Nsv
      ENDDO
      DO N = 1,NUMM4
        MEMBQ(N)%PAR%Isv = NV
        NV = NV + MATERIAL(MEMBQ(N)%PAR%MatID)%Nsv
      ENDDO
!!
      DO N = 1,NUMTR
        TRUSS(N)%PAR%Isv = NV
        NV = NV + MATERIAL(TRUSS(N)%PAR%MatID)%Nsv
      ENDDO
!!
      DO N = 1,NUMP3
        NPL = N
        Ipts = Ipts_PLATT(NPL)
        PLATT(N)%PAR%Isv = NV
        NV = NV + Ipts * MATERIAL(PLATT(N)%PAR%MatID)%Nsv
        PLATT(N)%PAR%Ist = NS
        NS = NS + Ipts
      ENDDO
      DO N = 1,NUMP4
        NPL = N
        Ipts = Ipts_PLATQ(NPL)
        PLATQ(N)%PAR%Isv = NV
        NV = NV + Ipts * MATERIAL(PLATQ(N)%PAR%MatID)%Nsv
        PLATQ(N)%PAR%Ist = NS
        NS = NS + Ipts
      ENDDO
!!
      DO N = 1,NUMBM
        BEAM(N)%PAR%Isv = NV
        Ipts = 16
        NV = NV + Ipts * MATERIAL(BEAM(N)%PAR%MatID)%Nsv
      ENDDO
!!
      DO N = 1,NUMSP
        SPRING(N)%PAR%Isv = NV
        NV = NV + MATERIAL(SPRING(N)%PAR%MatID)%Nsv
      ENDDO
      DO N = 1,NUMDM
        DAMPER(N)%PAR%Isv = NV
        NV = NV + MATERIAL(DAMPER(N)%PAR%MatID)%Nsv
      ENDDO
!!
      DO N = 1,NUMSC
        SPRING_BC(N)%Isv = NV
        NV = NV + MATERIAL(SPRING_BC(N)%MatID)%Nsv
      ENDDO
      DO N = 1,NUMVC
        DAMPER_BC(N)%Isv = NV
        NV = NV + MATERIAL(DAMPER_BC(N)%MatID)%Nsv
      ENDDO
!!
!! ### CONSTRUCT POINTERS FOR ROTATIONAL DEGREES OF FREEDOM.
!! To save storage and unneeded computations the data for rotational degrees
!! of freedom are stored in NODE(NUMNP+1:NUMRT), FORCE(NUMNP+1:NUMRT) and
!! MOTION(NUMNP+1:NUMRT). Up to this point NUMRT is only an estimate of the
!! number of locations needed. After the pointer NODE(*)%IRT is defined,
!! NUMRT is reset to the total number of locations actually needed in NODE,
!! MOTION and FORCE.
!!
      IF (NUMRT .GT. NUMNP) THEN
        DO N = 1,NUMNP
          NODE(N)%IRT = 0
        ENDDO
!!
!! Mark each node which occurs in a shell or beam element, or a type 1 spring
!! and damper element.
!!
        DO N = 1,NUMP3
          NODE(PLATT(N)%PAR%IX(1))%IRT = 1
          NODE(PLATT(N)%PAR%IX(2))%IRT = 1
          NODE(PLATT(N)%PAR%IX(3))%IRT = 1
        ENDDO
        DO N = 1,NUMP4
          NODE(PLATQ(N)%PAR%IX(1))%IRT = 1
          NODE(PLATQ(N)%PAR%IX(2))%IRT = 1
          NODE(PLATQ(N)%PAR%IX(3))%IRT = 1
          NODE(PLATQ(N)%PAR%IX(4))%IRT = 1
        ENDDO
        DO N = 1,NUMBM
          NODE(BEAM(N)%PAR%IX(1))%IRT = 1
          NODE(BEAM(N)%PAR%IX(2))%IRT = 1
        ENDDO
        DO N = 1,NUMSP
          IF (SPRING(N)%PAR%Type .EQ. 1) THEN
            NODE(SPRING(N)%PAR%IX(1))%IRT = 1
            NODE(SPRING(N)%PAR%IX(2))%IRT = 1
          ENDIF
        ENDDO
        DO N = 1,NUMDM
          IF (DAMPER(N)%PAR%Type .EQ. 1) THEN
            NODE(DAMPER(N)%PAR%IX(1))%IRT = 1
            NODE(DAMPER(N)%PAR%IX(2))%IRT = 1
          ENDIF
        ENDDO
!!
!! Count nodes with rotational degrees-of-freedom and construct pointer
!! to expanded nodal point arrays. Pointer (*%IRT) in expanded NODE array
!! points back to node with which it is associated; the node identifier
!! (*%ID) is assigned the same value as its "parent" (useful for error
!! message processing and doing a "changed" restart).
!!
        IRT = NUMNP
        DO N = 1,NUMNP
          IF (NODE(N)%IRT .EQ. 1) THEN
            IRT = IRT + 1
            IF (IRT .LE. NUMRT) THEN
              NODE(N)%IRT = IRT
              NODE(IRT)%IRT = N
              NODE(IRT)%ID = NODE(N)%ID
            ELSE
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.00'//                    &
     &          MSGL//'Rotational DOF''s Exceed Allowed Storage.'//            &
     &          MSGL//'Logic Error. Call KEY Associates.'                      &
     &          )
            ENDIF
          ENDIF
        ENDDO
        NUMRT = IRT
!!
!! Add rotational "nodes" to shell element definitions.
!!
        DO N = 1,NUMP3
          PLATT(N)%PAR%IX(4) = NODE(PLATT(N)%PAR%IX(1))%IRT
          PLATT(N)%PAR%IX(5) = NODE(PLATT(N)%PAR%IX(2))%IRT
          PLATT(N)%PAR%IX(6) = NODE(PLATT(N)%PAR%IX(3))%IRT
        ENDDO
        DO N = 1,NUMP4
          PLATQ(N)%PAR%IX(5) = NODE(PLATQ(N)%PAR%IX(1))%IRT
          PLATQ(N)%PAR%IX(6) = NODE(PLATQ(N)%PAR%IX(2))%IRT
          PLATQ(N)%PAR%IX(7) = NODE(PLATQ(N)%PAR%IX(3))%IRT
          PLATQ(N)%PAR%IX(8) = NODE(PLATQ(N)%PAR%IX(4))%IRT
        ENDDO
!!
!! Add rotational "nodes" to beam element definitions.
!!
        DO N = 1,NUMBM
          BEAM(N)%PAR%IX(3) = NODE(BEAM(N)%PAR%IX(1))%IRT
          BEAM(N)%PAR%IX(4) = NODE(BEAM(N)%PAR%IX(2))%IRT
        ENDDO
!!
!! Change rotational springs and dampers to point to rotational DOF's.
!!
        DO N = 1,NUMSP
          IF (SPRING(N)%PAR%Type .EQ. 1) THEN
            SPRING(N)%PAR%IX(1) = NODE(SPRING(N)%PAR%IX(1))%IRT
            SPRING(N)%PAR%IX(2) = NODE(SPRING(N)%PAR%IX(2))%IRT
          ENDIF
        ENDDO
        DO N = 1,NUMDM
          IF (DAMPER(N)%PAR%Type .EQ. 1) THEN
            DAMPER(N)%PAR%IX(1) = NODE(DAMPER(N)%PAR%IX(1))%IRT
            DAMPER(N)%PAR%IX(2) = NODE(DAMPER(N)%PAR%IX(2))%IRT
          ENDIF
        ENDDO
!!
      ENDIF
!!
!! ### INITIALIZE NODAL VELOCITIES WITH VELOCITY IC's.
!!
      IF (NUMIC .GT. 0) CALL INITIALIZE_MOTION
!!
!! ### TRUSS AND BEAM SECTION PROPERTIES.
!! Based on whether the user supplied truss/beam cross section dimensions
!! or supplied truss/beam cross section inertia properties compute the
!! complimentary quantities. Note that the use of section dimensions is
!! the preferred form of input. It is possible to specify "incompatible"
!! area and inertia data for which there does not exist positive or
!! realistic cross section dimensions.
!!
      DO N = 1,NUMS1
        NS1 = N
        CALL SECTION_PROPERTIES (NS1)
      ENDDO
!!
!! ### COMPLETE DATA SET FOR SPOT WELDS
!! Initialize spot weld data in order to minimize the computational
!! effort to impose the kinematic constraints a spot weld represents.
!!
      DO N = 1,NUMSW
        NSW = N
        CALL SPOT_WELD_INITIALIZATION (NSW)
      ENDDO
!!
!! ### PERIODIC BOUNDARY CONDITION DATA CHECK.
!! A periodic boundary condition allows geometry, loads and response
!! to repeat cyclically, either linearly or about an axis. By providing
!! an "advance," repeated helical response is obtained. The bounding
!! planes or surfaces defining the periodically repeated section need
!! be neither flat nor perpendicular to the direction of periodicity.
!! However, the bounding planes must self-similar, and, with the
!! current implementation, the meshing on the bounding planes must also
!! be identical. The nodal point sets pointed to by PERIODIC_BC(*)%S1ID
!! and PERIODIC_BC(*)%S2ID do not have to match pair-wise on input, but
!! must contain all the nodal point pairs needed to produce the periodic
!! motion.
!!
!! Loop over all periodic boundary conditions.
!!
      DO NCC = 1,NUMCC
!!
!! Do some simple checks on the node sets from face 1 and face 2.
!! First, Check to make certain both faces are "enumerated" sets.
!!
        SetID = PERIODIC_BC(NCC)%S1ID
        IF (NODE_SET(SetID)%Flag .NE. 'E') THEN
          WRITE (MSG1,'(I8)') PERIODIC_BC(NCC)%PerID
          WRITE (MSG2,'(I8)') NODE_SET(SetID)%SetID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.10'//                    &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'References An NPSET With ID:'//MSG2//                   &
     &          MSGL//'Which Is Not An Enumerated Node Set.'                   &
     &          )
        ENDIF

        SetID = PERIODIC_BC(NCC)%S2ID
        IF (NODE_SET(SetID)%Flag .NE. 'E') THEN
          WRITE (MSG1,'(I8)') PERIODIC_BC(NCC)%PerID
          WRITE (MSG2,'(I8)') NODE_SET(SetID)%SetID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.11'//                    &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'References An NPSET With ID:'//MSG2//                   &
     &          MSGL//'Which Is Not An Enumerated Node Set.'                   &
     &          )
        ENDIF
!!
!! Second, check to see that both node sets have the same number
!! of members.
!!
        SetID = PERIODIC_BC(NCC)%S1ID
        NN1 = NODE_SET(SetID)%Iend - NODE_SET(SetID)%Istart
        SetID = PERIODIC_BC(NCC)%S2ID
        NN2 = NODE_SET(SetID)%Iend - NODE_SET(SetID)%Istart

        IF (NN1 .NE. NN2) THEN
          WRITE (MSG1,'(I8)') PERIODIC_BC(NCC)%PerID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.12'//                    &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'References NPSET Records With'//                        &
     &          MSGL//'Different Numbers Of Nodal Points.'                     &
     &          )
        ENDIF

        Number_of_points = NN1 + 1
!!
!! Clear NODE(*)%IRB to use as a working array to reorder the enumeration
!! of surface 2 in NNPSETS into a pair-wise match with surface 1.
!!
        DO N = 1,NUMNP
          NODE(N)%IRB = 0
        ENDDO
!!
!! Temporarily use MOTION(*)%Ax,%Ay,%Az to hold the transformed positions
!! of surface 1. (If the input data is correct, the surfaces should
!! overlay.)
!!
!! Linear Periodicity --
!!
        IF (PERIODIC_BC(NCC)%TYPE .EQ. 'LINEAR') THEN
!!
!! Find "average" location of both surfaces.
!!
          X1ave = 0.0
          Y1ave = 0.0
          Z1ave = 0.0
          SetID = PERIODIC_BC(NCC)%S1ID
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
            X1ave = X1ave + MOTION(N)%Px
            Y1ave = Y1ave + MOTION(N)%Py
            Z1ave = Z1ave + MOTION(N)%Pz
          ENDDO

          X2ave = 0.0
          Y2ave = 0.0
          Z2ave = 0.0
          SetID = PERIODIC_BC(NCC)%S2ID
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
            X2ave = X2ave + MOTION(N)%Px
            Y2ave = Y2ave + MOTION(N)%Py
            Z2ave = Z2ave + MOTION(N)%Pz
          ENDDO
!!
!! Find translation distance from surface 1 to surface 2 in the
!! direction A = (Ax,Ay,Az). (Since the average position of face
!! 1 must map to the average position of face 2, the "direction"
!! of periodicity is the direction defined by Delta.
!!
          Delta_X = (X2ave - X1ave) / DBLE (Number_of_points)
          Delta_Y = (Y2ave - Y1ave) / DBLE (Number_of_points)
          Delta_Z = (Z2ave - Z1ave) / DBLE (Number_of_points)
          PERIODIC_BC(NCC)%Axis(1) = Delta_X
          PERIODIC_BC(NCC)%Axis(2) = Delta_Y
          PERIODIC_BC(NCC)%Axis(3) = Delta_Z
!!
!! Translate surface 1 by an amount Delta = (Delta_X,Delta_Y,Delta_Z)
!! so that surface 1 will coincide with surface 2.
!!
          SetID = PERIODIC_BC(NCC)%S1ID
          N = 0
          DO WHILE (NEXT_NP_ID(SetID,N))
            MOTION(N)%Ax = MOTION(N)%Px + Delta_X
            MOTION(N)%Ay = MOTION(N)%Py + Delta_Y
            MOTION(N)%Az = MOTION(N)%Pz + Delta_Z
          ENDDO
          Delta = ABS(Delta_X) + ABS(Delta_Y) + ABS(Delta_Z)
!!
!! Cyclical Periodicity --
!!
        ELSE IF (PERIODIC_BC(NCC)%TYPE .EQ. 'CYCLIC') THEN
!!
!! Retrieve angle (in degrees) between face 1 and face 2; do not convert
!! to radians
!!
          Theta = PERIODIC_BC(NCC)%Theta
!!
!! Retrieve advance in the event that this is Helical Periodicity.
!!
          Advance = Theta * PERIODIC_BC(NCC)%Advance
!!
!! Contruct orthogonal rotation R about the axis A = (Ax,Ay,Az)
!!
          Ax = PERIODIC_BC(NCC)%Axis(1)
          Ay = PERIODIC_BC(NCC)%Axis(2)
          Az = PERIODIC_BC(NCC)%Axis(3)
          Amag = SQRT (Ax*Ax + Ay*Ay + Az*Az)
          IF (Amag .GT. 1.0D-18) THEN
            Ax = Ax * (1.0 / Amag)
            Ay = Ay * (1.0 / Amag)
            Az = Az * (1.0 / Amag)
            PERIODIC_BC(NCC)%Axis(1) = Ax
            PERIODIC_BC(NCC)%Axis(2) = Ay
            PERIODIC_BC(NCC)%Axis(3) = Az
            C(1) = Ax
            C(2) = Ay
            C(3) = Az
          ELSE
            WRITE (MSG1,'(I8)') PERIODIC_BC(NCC)%PerID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.13'//                    &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'With Cyclical Periodicity Has A'//                      &
     &          MSGL//'Zero Length Axis Vector.'                               &
     &          )
          ENDIF
!!
!! With the angle of rotation, evaluate the trigonometric functions
!! SinTheta = Sin(Theta) and CosThem1 = Cos(Theta) - 1.
!!
          CosThem1 = COS ((DTR*Theta)) - 1.D+0
          SinTheta = SIN ((DTR*Theta))
!!
!! Construct rotation operator RTX. R = (Cos(phi)-1.0)*(I-CC) + Sin(phi)*(C x .)
!! The operator RTX is an "incremental" proper-orthogonal rotation. It is
!! constructed to be an "incremental operator" to avoid round-off.
!!
          DO n = 1,3
            DO m = 1,3
              RTX(m,n) = 0.D+0
              PRX(m,n) = -C(m) * C(n)
            ENDDO
          ENDDO
          DO n = 1,3
            PRX(n,n) = 1.D+0 + PRX(n,n)
          ENDDO
          RTX(1,2) = -(SinTheta*C(3))
          RTX(1,3) =  (SinTheta*C(2))
          RTX(2,1) =  (SinTheta*C(3))
          RTX(2,3) = -(SinTheta*C(1))
          RTX(3,1) = -(SinTheta*C(2))
          RTX(3,2) =  (SinTheta*C(1))
          DO n = 1,3
            DO m = 1,3
              RTX(m,n) = CosThem1*PRX(m,n) + RTX(m,n)
              PERIODIC_BC(NCC)%RTX(m,n) = RTX(m,n)
            ENDDO
          ENDDO
!!
!! Retrieve center of rotation.
!!
          Px = PERIODIC_BC(NCC)%Origin(1)
          Py = PERIODIC_BC(NCC)%Origin(2)
          Pz = PERIODIC_BC(NCC)%Origin(3)
!!
!! Rotate (and translate by an amount Advance) surface 1 so that
!! surface 1 will coincide with surface 2.
!!
          SetID = PERIODIC_BC(NCC)%S1ID
          N = 0
          Delta = 0.0
          DO WHILE (NEXT_NP_ID(SetID,N))
            Delta_X = RTX(1,1) * (MOTION(N)%Px - Px)                           &
     &              + RTX(1,2) * (MOTION(N)%Py - Py)                           &
     &              + RTX(1,3) * (MOTION(N)%Pz - Pz)
            Delta_Y = RTX(2,1) * (MOTION(N)%Px - Px)                           &
     &              + RTX(2,2) * (MOTION(N)%Py - Py)                           &
     &              + RTX(2,3) * (MOTION(N)%Pz - Pz)
            Delta_Z = RTX(3,1) * (MOTION(N)%Px - Px)                           &
     &              + RTX(3,2) * (MOTION(N)%Py - Py)                           &
     &              + RTX(3,3) * (MOTION(N)%Pz - Pz)
            MOTION(N)%Ax = MOTION(N)%Px + Delta_X + Advance * C(1)
            MOTION(N)%Ay = MOTION(N)%Py + Delta_Y + Advance * C(2)
            MOTION(N)%Az = MOTION(N)%Pz + Delta_Z + Advance * C(3)
            Delta = Delta + ABS(Delta_X) + ABS(Delta_Y) + ABS(Delta_Z)
          ENDDO
          Delta = Delta / DBLE (Number_of_points)
        ENDIF
!!
!! Taking surface 1 points one at a time, search surface 2 for a
!! match. Record location in NODE(N)%IRB for later reordering.
!! (Since NEXT_NP_ID(...) cannot be nested, NODE(*)%IRB will also
!! be used to construct the outer loop.)
!!
!! Move the outer-loop indicies out of NNPSETS into NODE(*)%IRB
!! (The sequence from NNPSETS must be perserved.)
!!
        N = 0
        SetID = PERIODIC_BC(NCC)%S1ID
        DO I = NODE_SET(SetID)%Istart,NODE_SET(SetID)%Iend
          N = N + 1
          NODE(N)%IRB = NNPSETS(I)
        ENDDO
!!
!! Now do nested loops looking for surface pairs.
!!
        SetID = PERIODIC_BC(NCC)%S2ID
        DO I = 1,Number_of_points
          N1 = NODE(I)%IRB
          Px = MOTION(N1)%Ax
          Py = MOTION(N1)%Ay
          Pz = MOTION(N1)%Az
          N2 = 0
          FOUND = .FALSE.
          DO WHILE (NEXT_NP_ID(SetID,N2) .AND. .NOT.FOUND)
            Distance = ABS(MOTION(N2)%Px - Px)                                 &
     &               + ABS(MOTION(N2)%Py - Py)                                 &
     &               + ABS(MOTION(N2)%Pz - Pz)
            IF (Distance .LE. 1.0D-5 * Delta) THEN
              FOUND = .TRUE.
              NODE(I)%IRB = N2
            ENDIF
          ENDDO
          IF (.NOT.FOUND) THEN
            WRITE (MSG1,'(I8)') PERIODIC_BC(NCC)%PerID
            WRITE (MSG2,'(I8)') NODE(N1)%ID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.001.13'//                    &
     &          MSGL//'PERIODBC Input Record ID:'//MSG1//                      &
     &          MSGL//'Cannot Find A Match On Surface 2 For'//                 &
     &          MSGL//'Nodal Point ID'//MSG2//' From Surface 1.'               &
     &          )
          ENDIF
        ENDDO
!!
!! Reorder the entries for surface 2 in NNPSETS using the
!! pairing constructed with NODE(*)%IRB
!!
        N = 0
        SetID = PERIODIC_BC(NCC)%S2ID
        DO I = NODE_SET(SetID)%Istart,NODE_SET(SetID)%Iend
          N = N + 1
          NNPSETS(I) = NODE(N)%IRB
        ENDDO
!!
!! Clean up "scratch" storage in MOTION and NODE.
!!
        DO N = 1,NUMNP
          NODE(N)%IRB  = 0
          MOTION(N)%Ax = 0.0
          MOTION(N)%Ay = 0.0
          MOTION(N)%Az = 0.0
        ENDDO
      ENDDO
!!
!! ### SLIDING INTERFACE DATA ARRAYS.
!! There are three data structures used for processing sliding interfaces:
!! SLIDING_NODE(1:NUMSN), CONTACT_SURFACE(1:NUMCE), CONTACT_NODE(1:NUMCN)%
!! The algorithm works by examining the sliding nodal points from the list
!! SLIDING_NODE to see if they have penetrated the element facets contained
!! in CONTACT_SURFACE which are currently opposing them. CONTACT_NODE is a
!! sublist of the nodes paraticipating in a sliding interface.
!!
      CALL SLIDING_INTERFACE_INIT
!!
!! ### LINKED LIST FOR RIGID BODY DOMAINS.
!! Rigid body domains are identified with Part ID's. All elements with the
!! same Part ID constitute a rigid body domain. NODE(*)%IRB is structured as
!! a nodal point linked list for each rigid body domain. The final result is
!! that RIGID_BODY(*)%FirstNP will contain the starting nodal point for each
!! rigid body domain. The domain linked lists are terminated by a minus one
!! in NODE(*)%IRB. First, loop over all element types and mark NODE(*)%IRB
!! for each nodal point occurring in a rigid body domain.
!!
      DO N = 1,NUMNP
        NODE(N)%IRB = 0
      ENDDO
!!
      DO M = 1,NUMRB
        ParID = RIGID_BODY(M)%ParID
        DO N = 1,NUMHX
          IF (ParID .EQ. HEXAH(N)%PAR%ParID) THEN
            DO i = 1,8
              NODE(HEXAH(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            HEXAH(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
        DO N = 1,NUMPX
          IF (ParID .EQ. PENTA(N)%PAR%ParID) THEN
            DO i = 1,6
              NODE(PENTA(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            PENTA(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
        DO N = 1,NUMTX
          IF (ParID .EQ. TETRA(N)%PAR%ParID) THEN
            DO i = 1,4
              NODE(TETRA(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            TETRA(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMLS
          IF (ParID .EQ. LSOLD(N)%PAR%ParID) THEN
            DO i = 1,8
              NODE(LSOLD(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            LSOLD(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMM4
          IF (ParID .EQ. MEMBQ(N)%PAR%ParID) THEN
            DO i = 1,4
              NODE(MEMBQ(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            MEMBQ(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
        DO N = 1,NUMM3
          IF (ParID .EQ. MEMBT(N)%PAR%ParID) THEN
            DO i = 1,3
              NODE(MEMBT(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            MEMBT(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMTR
          IF (ParID .EQ. TRUSS(N)%PAR%ParID) THEN
            DO i = 1,2
              NODE(TRUSS(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            TRUSS(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMP4
          IF (ParID .EQ. PLATQ(N)%PAR%ParID) THEN
            DO i = 1,4
              NODE(PLATQ(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            PLATQ(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
        DO N = 1,NUMP3
          IF (ParID .EQ. PLATT(N)%PAR%ParID) THEN
            DO i = 1,3
              NODE(PLATT(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            PLATT(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMBM
          IF (ParID .EQ. BEAM(N)%PAR%ParID) THEN
            DO i = 1,2
              NODE(BEAM(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            BEAM(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
!!
        DO N = 1,NUMSP
          IF (ParID .EQ. SPRING(N)%PAR%ParID) THEN
            DO i = 1,2
              NODE(SPRING(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            SPRING(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
        DO N = 1,NUMDM
          IF (ParID .EQ. DAMPER(N)%PAR%ParID) THEN
            DO i = 1,2
              NODE(DAMPER(N)%PAR%IX(i))%IRB = ParID
            ENDDO
            DAMPER(N)%PAR%ParID = -ParID
          ENDIF
        ENDDO
      ENDDO
!!
!! Scan NODE(*)%IRB in reverse to link successive locations.
!!
      DO M = 1,NUMRB
        RIGID_BODY(M)%FirstNP = -1
        DO N = NUMNP,1,-1
          IF (RIGID_BODY(M)%ParID .EQ. NODE(N)%IRB) THEN
            NODE(N)%IRB = RIGID_BODY(M)%FirstNP
            RIGID_BODY(M)%FirstNP = N
          ENDIF
        ENDDO
      ENDDO
!!
!! Determine if the mass and inertial properties of the rigid body have
!! been given explicitly by the user or are to be calculated from the
!! finite element discretization.
!!
      DO i = 1,NUMRB

        IF (RIGID_BODY(i)%Prop .NE. 0) THEN
!!
!! Use user-supplied mass and inertia.
!!
          ID = RIGID_BODY(i)%CMID
!!
!! Transfer mass and inertia from rigid-body-mass record to rigid-body record.
!!
          RIGID_BODY(i)%Mass = RIGID_BODY_MASS(ID)%Mass
          DO k = 1,9
            RIGID_BODY(i)%B(k,1) = RIGID_BODY_MASS(ID)%B(k,1)
          ENDDO
!!
!! Define center of mass location.
!!
          IF (RIGID_BODY_MASS(ID)%NPID .EQ. 0) THEN
            RIGID_BODY(i)%Px = RIGID_BODY_MASS(ID)%Pzero(1)
            RIGID_BODY(i)%Py = RIGID_BODY_MASS(ID)%Pzero(2)
            RIGID_BODY(i)%Pz = RIGID_BODY_MASS(ID)%Pzero(3)
            RIGID_BODY(i)%Vx = RIGID_BODY_MASS(ID)%Vel(1)
            RIGID_BODY(i)%Vy = RIGID_BODY_MASS(ID)%Vel(2)
            RIGID_BODY(i)%Vz = RIGID_BODY_MASS(ID)%Vel(3)
            RIGID_BODY(i)%Ox = RIGID_BODY_MASS(ID)%Omega(1)
            RIGID_BODY(i)%Oy = RIGID_BODY_MASS(ID)%Omega(2)
            RIGID_BODY(i)%Oz = RIGID_BODY_MASS(ID)%Omega(3)
          ELSE
            NPID = RIGID_BODY_MASS(ID)%NPID
            RIGID_BODY(i)%Px = MOTION(NPID)%Px
            RIGID_BODY(i)%Py = MOTION(NPID)%Py
            RIGID_BODY(i)%Pz = MOTION(NPID)%Pz
            RIGID_BODY(i)%Vx = MOTION(NPID)%Vx
            RIGID_BODY(i)%Vy = MOTION(NPID)%Vy
            RIGID_BODY(i)%Vz = MOTION(NPID)%Vz
            IF (NODE(NPID)%IRT .NE. 0) THEN
              NPID = NODE(NPID)%IRT
              RIGID_BODY(i)%Ox = MOTION(NPID)%Vx
              RIGID_BODY(i)%Oy = MOTION(NPID)%Vy
              RIGID_BODY(i)%Oz = MOTION(NPID)%Vz
            ELSE
              RIGID_BODY(i)%Ox = 0.0
              RIGID_BODY(i)%Oy = 0.0
              RIGID_BODY(i)%Oz = 0.0
            ENDIF
          ENDIF
          WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
          WRITE (MSG2,'(I8)') RIGID_BODY_MASS(ID)%RMID
          CALL USER_MESSAGE                                                    &
     &      (                                                                  &
     &      MSGL//'INFORM'//                                                   &
     &      MSGL//'INITIALIZE_WORKING_STORAGE.003.02'//                        &
     &      MSGL//'Total Mass And Inertia For Rigid Body ID:'//MSG1//          &
     &      MSGL//'Defined By Data From  Rigid Body Mass ID:'//MSG2            &
     &      )
        ELSE
          WRITE (MSG1,'(I8)') RIGID_BODY(i)%RBID
          CALL USER_MESSAGE                                                    &
     &      (                                                                  &
     &      MSGL//'INFORM'//                                                   &
     &      MSGL//'INITIALIZE_WORKING_STORAGE.003.01'//                        &
     &      MSGL//'Total Mass And Inertia For Rigid Body ID:'//MSG1//          &
     &      MSGL//'Will Be Computed Based on Nodal Masses.'                    &
     &      )
        ENDIF
      ENDDO
!!
!! ### LINK NRBC BOUNDARY SEGMENTS TO SOLID ELEMENTS.
!! Load segment ID's into auxillary data for nonreflecting BC's.
!!
      NXND = 0
      DO NNR = 1,NUMNR
        SetID = NONREFLECTING_BC(NNR)%SetID
        N = 0
        DO WHILE (NEXT_SEG_ID(SetID,N))
          NXND = NXND + 1
          NRBC_DATA(NXND)%SegID = N
!!
!! Until element match is found use .Mel to store parent nonreflecting BC ID.
!!
          NRBC_DATA(NXND)%Mel = NONREFLECTING_BC(NNR)%NRID
        ENDDO
      ENDDO
!!
!! Find solid element with which each nonreflecting boundary boundary
!! segment is associated.
!!
      DO NND = 1,NUMND
        FOUND = .FALSE.
        SegID = NRBC_DATA(NND)%SegID
!!
!! Scan hexahedrons for a match.
!!
        IF (SEGMENT(SegID)%PAR%IX(4) .NE. 0) THEN
          N = 1
          DO WHILE (.NOT.FOUND .AND. N .LE. NUMHX)
            DO i = 1,6
              IHX1 = HEXAH(N)%PAR%IX(HXface(1,i))
              IHX3 = HEXAH(N)%PAR%IX(HXface(3,i))
              MATCH = .FALSE.
              IF (IHX1 .EQ. SEGMENT(SegID)%PAR%IX(1)) THEN
                MATCH = IHX3 .EQ. SEGMENT(SegID)%PAR%IX(3)
              ELSE IF (IHX1 .EQ. SEGMENT(SegID)%PAR%IX(2)) THEN
                MATCH = IHX3 .EQ. SEGMENT(SegID)%PAR%IX(4)
              ELSE IF (IHX1 .EQ. SEGMENT(SegID)%PAR%IX(3)) THEN
                MATCH = IHX3 .EQ. SEGMENT(SegID)%PAR%IX(1)
              ELSE IF (IHX1 .EQ. SEGMENT(SegID)%PAR%IX(4)) THEN
                MATCH = IHX3 .EQ. SEGMENT(SegID)%PAR%IX(2)
              ENDIF
              IF (MATCH) THEN
                FOUND = .TRUE.
                NRBC_DATA(NND)%Mel = N
                NRBC_DATA(NND)%Type = 0
              ENDIF
            ENDDO
            N = N + 1
          ENDDO
        ENDIF
!!
!! Scan pentahedrons for a match.
!!
        IF (.NOT.FOUND .AND. SEGMENT(SegID)%PAR%IX(4) .NE. 0) THEN
          N = 1
          DO WHILE (.NOT.FOUND .AND. N .LE. NUMPX)
            DO i = 1,3
              IPX1 = PENTA(N)%PAR%IX(PXface(1,i))
              IPX3 = PENTA(N)%PAR%IX(PXface(3,i))
              MATCH = .FALSE.
              IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(1)) THEN
                MATCH = IPX3 .EQ. SEGMENT(SegID)%PAR%IX(3)
              ELSE IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(2)) THEN
                MATCH = IPX3 .EQ. SEGMENT(SegID)%PAR%IX(4)
              ELSE IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(3)) THEN
                MATCH = IPX3 .EQ. SEGMENT(SegID)%PAR%IX(1)
              ELSE IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(4)) THEN
                MATCH = IPX3 .EQ. SEGMENT(SegID)%PAR%IX(2)
              ENDIF
              IF (MATCH) THEN
                FOUND = .TRUE.
                NRBC_DATA(NND)%Mel = N
                NRBC_DATA(NND)%Type = 1
              ENDIF
            ENDDO
            N = N + 1
          ENDDO
        ELSE
          N = 1
          DO WHILE (.NOT.FOUND .AND. N .LE. NUMPX)
            DO i = 4,5
              IPX1 = PENTA(N)%PAR%IX(PXface(1,i))
              IPX2 = PENTA(N)%PAR%IX(PXface(2,i))
              IPX3 = PENTA(N)%PAR%IX(PXface(3,i))
              MATCH = .FALSE.
              IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(1)) THEN
                MATCH = IPX2 .EQ. SEGMENT(SegID)%PAR%IX(2)                     &
     &              .AND. IPX3 .EQ. SEGMENT(SegID)%PAR%IX(3)
              ELSE IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(2)) THEN
                MATCH = IPX2 .EQ. SEGMENT(SegID)%PAR%IX(3)                     &
     &              .AND. IPX3 .EQ. SEGMENT(SegID)%PAR%IX(1)
              ELSE IF (IPX1 .EQ. SEGMENT(SegID)%PAR%IX(3)) THEN
                MATCH = IPX2 .EQ. SEGMENT(SegID)%PAR%IX(1)                     &
     &              .AND. IPX3 .EQ. SEGMENT(SegID)%PAR%IX(2)
              ENDIF
              IF (MATCH) THEN
                FOUND = .TRUE.
                NRBC_DATA(NND)%Mel = N
                NRBC_DATA(NND)%Type = 1
              ENDIF
            ENDDO
            N = N + 1
          ENDDO
        ENDIF
!!
!! Scan tetrahedrons for a match.
!!
        IF (.NOT.FOUND .AND. SEGMENT(SegID)%PAR%IX(4) .EQ. 0) THEN
          N = 1
          DO WHILE (.NOT.FOUND .AND. N .LE. NUMTX)
            DO i = 1,4
              ITX1 = TETRA(N)%PAR%IX(TXface(1,i))
              ITX2 = TETRA(N)%PAR%IX(TXface(2,i))
              ITX3 = TETRA(N)%PAR%IX(TXface(3,i))
              MATCH = .FALSE.
              IF (ITX1 .EQ. SEGMENT(SegID)%PAR%IX(1)) THEN
                MATCH = ITX2 .EQ. SEGMENT(SegID)%PAR%IX(2)                     &
     &              .AND. ITX3 .EQ. SEGMENT(SegID)%PAR%IX(3)
              ELSE IF (ITX1 .EQ. SEGMENT(SegID)%PAR%IX(2)) THEN
                MATCH = ITX2 .EQ. SEGMENT(SegID)%PAR%IX(3)                     &
     &              .AND. ITX3 .EQ. SEGMENT(SegID)%PAR%IX(1)
              ELSE IF (ITX1 .EQ. SEGMENT(SegID)%PAR%IX(3)) THEN
                MATCH = ITX2 .EQ. SEGMENT(SegID)%PAR%IX(1)                     &
     &              .AND. ITX3 .EQ. SEGMENT(SegID)%PAR%IX(2)
              ENDIF
              IF (MATCH) THEN
                FOUND = .TRUE.
                NRBC_DATA(NND)%Mel = N
                NRBC_DATA(NND)%Type = 2
              ENDIF
            ENDDO
            N = N + 1
          ENDDO
        ENDIF
!!
!! Report any segments without a corresponding element.
!!
        IF (.NOT.FOUND) THEN
          WRITE (MSG1,'(I8)') NRBC_DATA(NND)%Mel
          WRITE (MSG2,'(I8)') SEGMENT(NRBC_DATA(NND)%SegID)%PAR%SegID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'INITIALIZE_WORKING_STORAGE.003.01'//                    &
     &          MSGL//'NRBC Input Record ID:'//MSG1//                          &
     &          MSGL//'Which References Segment ID:'//MSG2//                   &
     &          MSGL//'Does Not Match The Face Any Solid Element.'             &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ENDIF
!!
      ENDDO
!!
!! ### CONSTRAINED NODES.
!! Locate mid-side nodes wrt bounding nodes.
!!
      IF (NUMNC .GT. 0) CALL INITIALIZE_CONSTRAINED_NODES
!!
!! ### REZONE DATA.
!! If rezoning has been requested, initialize data need for rezoning.
!!
      IF (CONTROL%REZONE .NE. 0) CALL INITIALIZE_REZONE_DATA
!!
!! Record time spent in scanning input records.
!!
!SPEC_CPU2000      CALL TIMER (1)
!!
      IF (ERROR%COUNT .GT. 0) THEN
        WRITE (MSG1,'(I8)') ERROR%COUNT
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'INITIALIZE_WORKING_STORAGE.003.02'//                          &
     &    MSGL//'Total Number Of Initilization Errors In Input:'//MSG1//       &
     &    MSGL//'Execution Terminated By Program.'                             &
     &    )
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_CONSTRAINED_NODES
!!
!! Copyright (c) by KEY Associates; 21-OCT-1995 16:10:20.00
!!
!! Purpose: Check location of mid-side node and compute weight.
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE constrained_node_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered INITIALIZE_CONSTRAINED_NODES.'
      WRITE (IO_UNIT%LELO,*) ' '

      DO i = 1,NUMNC
        IMS = CONSTRAINED_NODE(i)%CNID
        ID1 = CONSTRAINED_NODE(i)%NPID(1)
        ID2 = CONSTRAINED_NODE(i)%NPID(2)
        Px0 = MOTION(IMS)%Px + MOTION(IMS)%Ux
        Py0 = MOTION(IMS)%Py + MOTION(IMS)%Uy
        Pz0 = MOTION(IMS)%Pz + MOTION(IMS)%Uz
        Px1 = MOTION(ID1)%Px + MOTION(ID1)%Ux
        Py1 = MOTION(ID1)%Py + MOTION(ID1)%Uy
        Pz1 = MOTION(ID1)%Pz + MOTION(ID1)%Uz
        Px2 = MOTION(ID2)%Px + MOTION(ID2)%Ux
        Py2 = MOTION(ID2)%Py + MOTION(ID2)%Uy
        Pz2 = MOTION(ID2)%Pz + MOTION(ID2)%Uz
!!
!! First, check to see if the mid-side node is close enough to the
!! line connecting the bounding nodes.
!!
        Distance_01 = SQRT ((Px1-Px0)**2 + (Py1-Py0)**2 + (Pz1-Pz0)**2)
        Distance_02 = SQRT ((Px2-Px0)**2 + (Py2-Py0)**2 + (Pz2-Pz0)**2)
        Distance_12 = SQRT ((Px1-Px2)**2 + (Py1-Py2)**2 + (Pz1-Pz2)**2)

        Delta = ABS(Distance_12 - Distance_01 - Distance_02)

        IF (Delta .GT. 0.1*Distance_12) THEN
          WRITE (MSG1,'(I8)') CONSTRAINED_NODE(i)%ID
          WRITE (MSG2,'(I8)') NODE(CONSTRAINED_NODE(i)%CNID)%ID
          CALL USER_MESSAGE                                                    &
     & (                                                                       &
     & MSGL//'WARN'//                                                          &
     & MSGL//'INITIALIZE_CONSTRAINED_NODES.01.001'//                           &
     & MSGL//'Nodal Point Constraint (NPCON1) Input Record ID:'//MSG1//        &
     & MSGL//'References A Mid-Side Node With ID:'//MSG2//                     &
     & MSGL//'That Lies Too Far Off The Line'                                  &
     &     //'Connecting The Bounding Nodes.'                                  &
     & )
          ERROR%COUNT = ERROR%COUNT + 1
!!
!! Put mid-side node exactly on line between bounding nodes by constructing
!! interpolation weights.
!!
        ELSE
          Weight = Distance_02 / (Distance_01 + Distance_02)

          CONSTRAINED_NODE(i)%Weight = Weight
          Wg1 = Weight
          Wg2 = 1.0 - Weight

          MOTION(IMS)%Px = Wg1 * MOTION(ID1)%Px + Wg2 * MOTION(ID2)%Px
          MOTION(IMS)%Py = Wg1 * MOTION(ID1)%Py + Wg2 * MOTION(ID2)%Py
          MOTION(IMS)%Pz = Wg1 * MOTION(ID1)%Pz + Wg2 * MOTION(ID2)%Pz

          MOTION(IMS)%Ux = Wg1 * MOTION(ID1)%Ux + Wg2 * MOTION(ID2)%Ux
          MOTION(IMS)%Uy = Wg1 * MOTION(ID1)%Uy + Wg2 * MOTION(ID2)%Uy
          MOTION(IMS)%Uz = Wg1 * MOTION(ID1)%Uz + Wg2 * MOTION(ID2)%Uz
!!
!! Fix up any problems with initial conditions on velocity. (Assumes
!! bounding nodes are correct.)
!!
          MOTION(IMS)%Vx = Wg1 * MOTION(ID1)%Vx + Wg2 * MOTION(ID2)%Vx
          MOTION(IMS)%Vy = Wg1 * MOTION(ID1)%Vy + Wg2 * MOTION(ID2)%Vy
          MOTION(IMS)%Vz = Wg1 * MOTION(ID1)%Vz + Wg2 * MOTION(ID2)%Vz

          IF (NODE(IMS)%IRT .GT. 0) THEN
            IMS = NODE(IMS)%IRT
            ID1 = NODE(ID1)%IRT
            ID2 = NODE(ID2)%IRT

            MOTION(IMS)%Vx = Wg1 * MOTION(ID1)%Vx + Wg2 * MOTION(ID2)%Vx
            MOTION(IMS)%Vy = Wg1 * MOTION(ID1)%Vy + Wg2 * MOTION(ID2)%Vy
            MOTION(IMS)%Vz = Wg1 * MOTION(ID1)%Vz + Wg2 * MOTION(ID2)%Vz
          ENDIF

        ENDIF

      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_REZONE_DATA
!!
!! Copyright (c) by KEY Associates; 15-OCT-1995 21:03:54.00
!!
!! Purpose: Initialize rezone maximum ID variables and initial angle
!! between plate pairs.
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
      USE nrbc_data_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE plate_pair_
!!
      USE enumerated_sets_
      USE stress_
      USE state_variables_
      USE coord_
      USE indx_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Find maximum ID's in use.
!!
      IDXNP = 0
      DO i = 1,NUMNP
        IDXNP = MAX (IDXNP, NODE(i)%ID)
      ENDDO
!!
      IDXDC = 0
      DO i = 1,NUMDC
        IDXDC = MAX (IDXDC, DISPLACEMENT_BC(i)%DBCID)
      ENDDO
!!
      IDXNC = 0
      DO i = 1,NUMNC
        IDXNC = MAX (IDXNC, CONSTRAINED_NODE(i)%ID)
      ENDDO
!!
      IDXPP = 0
      DO i = 1,NUMPP
        IDXPP = MAX (IDXPP, PLATE_PAIR(i)%ID)
      ENDDO
!!
!! Construct plate element normal vectors for use in computing the initial
!! angle between adjacent plates.
!!
      DO i = 1,NUMP4
!!
!! Construct diagonal vectors.
!!
        P1 = MOTION(PLATQ(i)%PAR%IX(1))%Px
        P2 = MOTION(PLATQ(i)%PAR%IX(2))%Px
        P3 = MOTION(PLATQ(i)%PAR%IX(3))%Px
        P4 = MOTION(PLATQ(i)%PAR%IX(4))%Px
        Ax = P3 - P1
        Bx = P4 - P2
        P1 = MOTION(PLATQ(i)%PAR%IX(1))%Py
        P2 = MOTION(PLATQ(i)%PAR%IX(2))%Py
        P3 = MOTION(PLATQ(i)%PAR%IX(3))%Py
        P4 = MOTION(PLATQ(i)%PAR%IX(4))%Py
        Ay = P3 - P1
        By = P4 - P2
        P1 = MOTION(PLATQ(i)%PAR%IX(1))%Pz
        P2 = MOTION(PLATQ(i)%PAR%IX(2))%Pz
        P3 = MOTION(PLATQ(i)%PAR%IX(3))%Pz
        P4 = MOTION(PLATQ(i)%PAR%IX(4))%Pz
        Az = P3 - P1
        Bz = P4 - P2
!!
!! Plate normal based on cross product of diagonals and inverse of the magnitude
!!
        Ex = Ay*Bz - Az*By
        Ey = Az*Bx - Ax*Bz
        Ez = Ax*By - Ay*Bx
        Emi = 1.0 / SQRT (Ex*Ex + Ey*Ey + Ez*Ez)
!!
        PLATQ(i)%ADP%Ax = Ex * Emi
        PLATQ(i)%ADP%Ay = Ey * Emi
        PLATQ(i)%ADP%Az = Ez * Emi
!!
!! Assign initial rezoning level to be 1.
!!
        PLATQ(i)%ADP%Level = 1
!!
      ENDDO
!!
!! Compute the initial value of the cosine of the angle between adjacent
!! plate normals
!!
      DO N = 1,NUMPP
!!
!! Paired element ID's
!!
        ID1 = PLATE_PAIR(N)%IDS1%ID
        ID2 = PLATE_PAIR(N)%IDS2%ID
!!
!! Retrieve paired element normal vectors.
!!
        Ex = PLATQ(ID1)%ADP%Ax
        Ey = PLATQ(ID1)%ADP%Ay
        Ez = PLATQ(ID1)%ADP%Az
        Fx = PLATQ(ID2)%ADP%Ax
        Fy = PLATQ(ID2)%ADP%Ay
        Fz = PLATQ(ID2)%ADP%Az
!!
!! Compute the cosine of the angle between the element normals.
!!
        PLATE_PAIR(N)%CSPINI = Ex*Fx + Ey*Fy + Ez*Fz
!!
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE BUILD_PLATE_PAIRS (ACTION)
!!
!! Copyright (c) by KEY Associates; 16-OCT-1995 20:00:04.00
!!
!! Purpose: Build table of plate pairs. That is, identify all interfaces
!! between adjacent quadrilateral plates. This table is used to rapidly
!! evaluate rezoning criteria based on the behavior of one plate with
!! respect to another plate.
!!
      USE shared_common_data
      USE platq_
      USE plate_pair_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Argument.
      CHARACTER, INTENT(IN) :: ACTION*(*)  ! I/- First 'COUNT' then 'BUILD'
!!
      INTEGER, DIMENSION(5) :: IX      ! -/- Local temporary storage
      INTEGER, DIMENSION(5) :: JX      ! -/- Local temporary storage
      INTEGER               :: KOUNT=1 ! -/- Local integer flag
      INTEGER               :: BUILD=0 ! -/- Local integer flag
      INTEGER               :: IACTION ! -/- Local integer for ACTION
!!
!! Determine ACTION (count or build) being requested.
!!
      IF (INDEX(ACTION,'COUNT') .NE. 0) THEN
        IACTION = KOUNT
      ELSE IF (INDEX(ACTION,'BUILD') .NE. 0) THEN
        IACTION = BUILD
      ELSE
        L = LEN (ACTION)
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'BUILD_PLATE_PAIRS.002.01'//                             &
     &          MSGL//'Unknown ACTION Requested:'//ACTION(1:L)//               &
     &          MSGL//'Kown Actions: COUNT, BUILD'                             &
     &          )
      ENDIF
!!
!! Examine each side in each quadrilateral and look for a matching side in
!! another quadrilateral.
!!
      NUMPP = 0
      IDXPP = 0
      DO I = 1,NUMP4
!!
!! Examine each side in turn.
!!
!!      L---------K     Position of sides in
!!      |    4    |     relation to element
!!      |         |     Definition:
!!      |1       3|
!!      |         |     NP's    = (I,J,K,L)
!!      |    2    |     Sides   = (1,2,3,4)
!!      I---------J     IX(1:5) = (L,I,J,K,L)
!!
        IX(1) = PLATQ(I)%PAR%IX(4)
        IX(2) = PLATQ(I)%PAR%IX(1)
        IX(3) = PLATQ(I)%PAR%IX(2)
        IX(4) = PLATQ(I)%PAR%IX(3)
        IX(5) = PLATQ(I)%PAR%IX(4)
!!
!! Loop over all four sides.
!!
        DO IS = 1,4
          Imin = MIN (IX(IS),IX(IS+1))
          Imax = MAX (IX(IS),IX(IS+1))
!!
!! Look at all "future" quadrilaterals. (By not stopping at the first
!! match, T-intersections, X-intersections et cetera will be found.)
!!
          DO J = I+1,NUMP4
            JX(1) = PLATQ(J)%PAR%IX(4)
            JX(2) = PLATQ(J)%PAR%IX(1)
            JX(3) = PLATQ(J)%PAR%IX(2)
            JX(4) = PLATQ(J)%PAR%IX(3)
            JX(5) = PLATQ(J)%PAR%IX(4)
!!
!! Look at all four sides of candidate neighbor quadrilateral.
!!
            DO JS = 1,4
              Jmin = MIN (JX(JS),JX(JS+1))
              Jmax = MAX (JX(JS),JX(JS+1))
!!
              IF (Imin .EQ. Jmin) THEN
                IF (Imax .EQ. Jmax) THEN
!!
!! A matching side has been found.
!!
                  IF (IACTION .EQ. KOUNT) THEN
                    NUMPP = NUMPP + 1
                  ELSE ! IF (IACTION .EQ. BUILD) THEN
                    NUMPP = NUMPP + 1
                    IDXPP = IDXPP + 1
                    PLATE_PAIR(NUMPP)%ID = IDXPP
                    PLATE_PAIR(NUMPP)%IDS1%ID = I
                    PLATE_PAIR(NUMPP)%IDS1%IS = IS
                    PLATE_PAIR(NUMPP)%IDS2%ID = J
                    PLATE_PAIR(NUMPP)%IDS2%IS = JS
                    PLATE_PAIR(NUMPP)%CSPINI = 0.0
                    PLATE_PAIR(NUMPP)%CosPhi = 0.0
                    PLATE_PAIR(NUMPP)%CosDot = 0.0
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE SLIDING_INTERFACE_INIT
!!
!! Copyright (c) by KEY Associates; 29-JAN-1995 19:07:18.00
!!
!! Purpose: This module initializes the data structures required to
!! implement the sliding interface algorithms. The sliding interface
!! data is a "separate" set of nodes and elements that have data passed
!! to them from the "underlying" parent finite elements of the model.
!! The sliding interface algorithms then looks for and process any
!! contacts taking place amongst its members. When the processing is
!! finished, the results are passed back to the nodes and elements
!! making up the whole model.
!!
!! There are three data structures used for processing sliding interfaces:
!! SLIDING_NODE(1:NUMSN), CONTACT_SURFACE(1:NUMCE), CONTACT_NODE(1:NUMCN)%
!! The algorithm works by examining the sliding nodal points from the list
!! SLIDING_NODE to see if they have penetrated the element facets contained
!! in CONTACT_SURFACE which are currently opposing them. CONTACT_NODE is a
!! sublist of the entire simulation node list that is paraticipating in a
!! sliding interface.
!!
      USE shared_common_data
      USE sliding_interface_
      USE segment_
      USE node_
      USE motion_
      USE sliding_node_
      USE contact_surface_
      USE contact_node_
      USE node_set_
      USE element_set_
      USE segment_set_
      USE enumerated_sets_, ONLY: NNPSETS,NSGSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                    &
     &          SetID,           &
     &          NEIGHBOR(10)       ! -/- Local work array used to build ICE
      REAL(KIND(0D0))            &
     &          Xce(4),          & ! -/- Contact element nodal coordinates     
     &          Yce(4),          & ! -/- Contact element nodal coordinates     
     &          Zce(4)             ! -/- Contact element nodal coordinates
      LOGICAL                    &
     &          FOUND,           &
     &          NEXT_NP_ID,      & ! External function for unpacking node sets
     &          NEXT_SEG_ID        ! External function for unpacking segment set
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered SLIDING_INTERFACE_INIT.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!!      Isym  = 0/1/2 = symmetric/master-slave/slave-master
!!      Itype = 0/1   = segment-set/nodal-set
!!
!! Set End (calculation), Factor, Capture, Border and Sort_Frequency
!! to default values in the event zero's were input.
!!
      DO Nsi = 1,NUMSI
        IF (SLIDING_INTERFACE(Nsi)%End .EQ. 0.0) THEN
          SLIDING_INTERFACE(Nsi)%End = HUGE( SLIDING_INTERFACE(Nsi)%End)
        ENDIF
        IF (SLIDING_INTERFACE(Nsi)%Factor .EQ. 0.0) THEN
          SLIDING_INTERFACE(Nsi)%Factor = PARAMVALUE%Factor
        ENDIF
        IF (SLIDING_INTERFACE(Nsi)%Capture .EQ. 0.0) THEN
          SLIDING_INTERFACE(Nsi)%Capture = PARAMVALUE%Capture
        ENDIF
        IF (SLIDING_INTERFACE(Nsi)%Border .EQ. 0.0) THEN
          SLIDING_INTERFACE(Nsi)%Border = PARAMVALUE%Border
        ENDIF
        IF (SLIDING_INTERFACE(Nsi)%Sort_Freq .EQ. 0.0) THEN
          SLIDING_INTERFACE(Nsi)%Sort_Freq = NINT (PARAMVALUE%Sort_Freq)
        ENDIF
      ENDDO
!!
!! For single-surface interfaces, force interface to be slave-master.
!! (This is the logical place for this operation, however, it must be
!! done in GET_AUXILIARY_STORAGE_LENGTH where the values of .Type and
!! .Isym are first interogated.)
!!
!!      DO Nsi = 1,NUMSI
!!        IF (SLIDING_INTERFACE(Nsi)%Type .EQ. 1) THEN
!!          SLIDING_INTERFACE(Nsi)%Isym = 2
!!        ENDIF
!!      ENDDO
!!
!! Modify sliding interface *.Begin (calculation) and *.End (calculation) on
!! the basis of values read into INTERFACE_TIME%
!!
      DO i = 1,NUMIT
        Nsi = INTERFACE_TIME%SI(i)%ID
        SLIDING_INTERFACE(Nsi)%Begin = INTERFACE_TIME%SI(i)%Begin
        SLIDING_INTERFACE(Nsi)%End   = INTERFACE_TIME%SI(i)%End
      ENDDO
!!
!! Initialize to zero the sliding-node and contact-element start and
!! end pointers.
!!
      DO i = 1,NUMSI
        SLIDING_INTERFACE(i)%ISN1bgn = 0
        SLIDING_INTERFACE(i)%ISN1end = 0
        SLIDING_INTERFACE(i)%ICE1bgn = 0
        SLIDING_INTERFACE(i)%ICE1end = 0
        SLIDING_INTERFACE(i)%ISN2bgn = 0
        SLIDING_INTERFACE(i)%ISN2end = 0
        SLIDING_INTERFACE(i)%ICE2bgn = 0
        SLIDING_INTERFACE(i)%ICE2end = 0
      ENDDO
!!
!! ### BUILD SLIDING NODAL POINT LIST: SLIDING_NODE.
!! Record beginning and ending locations for each sliding interface.
!!
      NXSN = 1
      DO Nsi = 1,NUMSI
        Isym = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Clear NODE(*)%IRB to use as a working array to construct sliding interface
!! pointers. (We can only do this because NODE(*)%IRB has not yet been defined
!! by rigid body initialization.)
!!
        DO N = 1,NUMNP
          NODE(N)%IRB = 0
        ENDDO
!!
!! Collect sliding nodal points from Side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
          WRITE (MSG1,'(I8)') SLIDING_INTERFACE(Nsi)%SIID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SLIDING_INTERFACE_INIT.001.01'//                        &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Invalid Use Of Node Set On Side 1'//                    &
     &          MSGL//'With Symmetric Or Master-Slave Type.'                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                IF (SEGMENT(N)%PAR%IX(4) .EQ. 0) THEN
                  NODE(SEGMENT(N)%PAR%IX(1))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(2))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(3))%IRB = 1
                ELSE
                  NODE(SEGMENT(N)%PAR%IX(1))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(2))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(3))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(4))%IRB = 1
                ENDIF
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Count the number of entries in NODE(*)%IRB and transfer data to SLIDING_
!! NODE(*)%Nsn
!!
        SLIDING_INTERFACE(Nsi)%ISN1bgn = NXSN
        DO N = 1,NUMNP
          IF (NODE(N)%IRB .NE. 0) THEN
            SLIDING_NODE(NXSN)%Nsn = N
            NXSN = NXSN + 1
          ENDIF
        ENDDO
        SLIDING_INTERFACE(Nsi)%ISN1end = NXSN - 1
!!
!! Clear NODE(*)%IRB to use as a working array to construct sliding interface
!! pointers.
!!
        DO N = 1,NUMNP
          NODE(N)%IRB = 0
        ENDDO
!!
!! Collect sliding nodal points from Side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
          WRITE (MSG1,'(I8)') SLIDING_INTERFACE(Nsi)%SIID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SLIDING_INTERFACE_INIT.001.02'//                        &
     &          MSGL//'SLIDE Input Record ID:'//MSG1//                         &
     &          MSGL//'Invalid Use Of Node Set On Side 2'//                    &
     &          MSGL//'With Symmetric Or Slave-Master Type.'                   &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                IF (SEGMENT(N)%PAR%IX(4) .EQ. 0) THEN
                  NODE(SEGMENT(N)%PAR%IX(1))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(2))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(3))%IRB = 1
                ELSE
                  NODE(SEGMENT(N)%PAR%IX(1))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(2))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(3))%IRB = 1
                  NODE(SEGMENT(N)%PAR%IX(4))%IRB = 1
                ENDIF
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Count the number of entries in NODE(*)%IRB and transfer data to SLIDING_
!! NODE(*)%Nsn
!!
        SLIDING_INTERFACE(Nsi)%ISN2bgn = NXSN
        DO N = 1,NUMNP
          IF (NODE(N)%IRB .NE. 0) THEN
            SLIDING_NODE(NXSN)%Nsn = N
            NXSN = NXSN + 1
          ENDIF
        ENDDO
        SLIDING_INTERFACE(Nsi)%ISN2end = NXSN - 1
!!
      ENDDO
!!
!! Consistentcy check. The total number of sliding nodes NXSN-1
!! used in the sliding interface should equal the number of
!! sliding nodes NUMSN computed earlier.
!!
      IF (NUMSN .NE. NXSN-1) THEN
        WRITE (MSG1,'(I8)') NUMSN
        WRITE (MSG2,'(I8)') NXSN-1
        CALL USER_MESSAGE                                                      &
     &  (                                                                      &
     &  MSGL//'FATAL'//                                                        &
     &  MSGL//'SLIDING_INTERFACE_INIT.002.01'//                                &
     &  MSGL//'Logic Error. Original Count Of Sliding Nodes'//                 &
     &  MSGL//'NUMSN Does Not Equal Count Obtained During'//                   &
     &  MSGL//'Initialization. NUMSN:'//MSG1//'  NXSN-1:'//MSG2                &
     &  )
      ENDIF
!!
!! ### BUILD CONTACT ELEMENT LIST: CONTACT_SURFACE.
!! Record beginning and ending locations for each sliding interface.
!!
      NXCE = 1
      DO Nsi = 1,NUMSI
        Isym = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Collect contact elements from Side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        SLIDING_INTERFACE(Nsi)%ICE1bgn = NXCE
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
!!
!! This error condition processed above in building SLIDING_NODE%
!!
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,4
                  CONTACT_SURFACE(NXCE)%IX(k) = SEGMENT(N)%PAR%IX(k)
                ENDDO
                NXCE = NXCE + 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        SLIDING_INTERFACE(Nsi)%ICE1end = NXCE - 1
!!
!! Collect contact elements from Side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        SLIDING_INTERFACE(Nsi)%ICE2bgn = NXCE
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
!!
!! This error condition processed above in building SLIDING_NODE%
!!
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,4
                  CONTACT_SURFACE(NXCE)%IX(k) = SEGMENT(N)%PAR%IX(k)
                ENDDO
                NXCE = NXCE + 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        SLIDING_INTERFACE(Nsi)%ICE2end = NXCE - 1
!!
      ENDDO
!!
!! Consistentcy check. The total number of contact elements NXCE-1
!! used in the sliding interface should equal the number of
!! contact elements NUMCE computed earlier.
!!
      IF (NUMCE .NE. NXCE-1) THEN
        WRITE (MSG1,'(I8)') NUMCE
        WRITE (MSG2,'(I8)') NXCE-1
        CALL USER_MESSAGE                                                      &
     &    (                                                                    &
     &    MSGL//'FATAL'//                                                      &
     &    MSGL//'SLIDING_INTERFACE_INIT.002.02'//                              &
     &    MSGL//'Logic Error. Original Count Of Contact Elements'//            &
     &    MSGL//'NUMCE Does Not Equal Count Obtained During'//                 &
     &    MSGL//'Initialization. NUMCE:'//MSG1//'  NXCE-1:'//MSG2              &
     &    )
      ENDIF
!!
!! ### BUILD LIST OF ALL NODES USED IN SLIDING INTERFACES: CONTACT_NODE
!! Redefine all nodal point ID's (global pointers) in SLIDING_NODE(*)%Nsn
!! and CONTACT_SURFACE(*)%IX(1:4) to pointers into CONTACT_NODES(*)
!!
!! Clear NODE(*)%IRB to use as a working array to construct sliding interface
!! pointers.
!!
      DO N = 1,NUMNP
        NODE(N)%IRB = 0
      ENDDO
!!
!! Mark every use of a nodal point in a sliding interface. Note that the native
!! triangles have a zero in *%IX(4).
!!
      DO N = 1,NUMSN
        NODE(SLIDING_NODE(N)%Nsn)%IRB = 1
      ENDDO
      DO M = 1,NUMCE
        IF (CONTACT_SURFACE(M)%IX(4) .EQ. 0) THEN
          NODE(CONTACT_SURFACE(M)%IX(1))%IRB = 1
          NODE(CONTACT_SURFACE(M)%IX(2))%IRB = 1
          NODE(CONTACT_SURFACE(M)%IX(3))%IRB = 1
        ELSE
          NODE(CONTACT_SURFACE(M)%IX(1))%IRB = 1
          NODE(CONTACT_SURFACE(M)%IX(2))%IRB = 1
          NODE(CONTACT_SURFACE(M)%IX(3))%IRB = 1
          NODE(CONTACT_SURFACE(M)%IX(4))%IRB = 1
        ENDIF
      ENDDO
!!
!! Define global nodal point ID's in CONTACT_NODE(*)%NPID
!!
      NXCN = 0
      DO N = 1,NUMNP
        IF (NODE(N)%IRB .NE. 0) THEN
          NXCN = NXCN + 1
          CONTACT_NODE(NXCN)%NPID = N
        ENDIF
      ENDDO
!!
!! Consistentcy check. The number of model nodal points NXCN used in
!! the sliding interface should equal the number of contact nodes
!! NUMCT computed earlier.
!!
      IF (NUMCN .NE. NXCN) THEN
        WRITE (MSG1,'(I8)') NUMCN
        WRITE (MSG2,'(I8)') NXCN
        CALL USER_MESSAGE                                                      &
     &          (                                                              &
     &          MSGL//'FATAL'//                                                &
     &          MSGL//'SLIDING_INTERFACE_INIT.002.03'//                        &
     &          MSGL//'Logic Error. Original Count Of Contact Nodes'//         &
     &          MSGL//'NUMCN Does Not Equal Count Obtained During'//           &
     &          MSGL//'Initialization. NUMCN:'//MSG1//'   NXCN:'//MSG2         &
     &          )
      ENDIF
!!
!! Redefine global nodal point ID's in SLIDING_NODE and CONTACT_SURFACE
!!
      DO N = 1,NUMSN
        M = 1
        DO WHILE                                                               &
     &          (                                                              &
     &          SLIDING_NODE(N)%Nsn .NE. CONTACT_NODE(M)%NPID                  &
     &          .AND.                                                          &
     &          M .LT. NUMCN                                                   &
     &          )
          M = M + 1
        ENDDO
        SLIDING_NODE(N)%Nsn = M
      ENDDO
!!
      DO M = 1,NUMCE
        DO k = 1,3+MIN(1,CONTACT_SURFACE(M)%IX(4))
          N = 1
          DO WHILE                                                             &
     &          (                                                              &
     &          CONTACT_SURFACE(M)%IX(k) .NE. CONTACT_NODE(N)%NPID             &
     &          .AND.                                                          &
     &          N .LT. NUMCN                                                   &
     &          )
            N = N + 1
          ENDDO
          CONTACT_SURFACE(M)%IX(k) = N
        ENDDO
      ENDDO
!!
!! ### FIND AND RECORD CONTACT_SURFACE NEIGHBORING ELEMENTS.
!! Complete CONTACT_SURFACE by locating elements adjacent to contact elements.
!!
      DO N = 1,NUMCE
        IF (CONTACT_SURFACE(N)%IX(1) .EQ. 0) THEN
          CONTACT_SURFACE(N)%NX(1) = -N
          CONTACT_SURFACE(N)%NX(2) = -N
          CONTACT_SURFACE(N)%NX(3) = -N
        ELSE
          CONTACT_SURFACE(N)%NX(1) = -N
          CONTACT_SURFACE(N)%NX(2) = -N
          CONTACT_SURFACE(N)%NX(3) = -N
          CONTACT_SURFACE(N)%NX(4) = -N
        ENDIF
      ENDDO
!!
      DO Nsi = 1,NUMSI
        ICE1bgn = SLIDING_INTERFACE(Nsi)%ICE1bgn
        ICE2end = SLIDING_INTERFACE(Nsi)%ICE2end
        DO N = ICE1bgn,ICE2end
          i = 3+MIN(1,CONTACT_SURFACE(N)%IX(4))
          DO j = 1,3+MIN(1,CONTACT_SURFACE(N)%IX(4))
            I1 = CONTACT_SURFACE(N)%IX(i)
            I2 = CONTACT_SURFACE(N)%IX(j)
!!
!! Find all elements that have I1 as a nodal point.
!!
            Inext = 0
            DO M = ICE1bgn,ICE2end
              DO K = 1,3+MIN(1,CONTACT_SURFACE(M)%IX(4))
                IF (I1 .EQ. CONTACT_SURFACE(M)%IX(K)) THEN
                  Inext = Inext + 1
                  NEIGHBOR(Inext) = M
                ENDIF
              ENDDO
            ENDDO
!!
!! Check all neighbors for occurance of I2.
!!
            M = 0
            FOUND = .FALSE.
            DO WHILE (.NOT.FOUND .AND. M .LT. Inext)
              M = M + 1
              K = 0
              Kmax = 3+MIN(1,CONTACT_SURFACE(NEIGHBOR(M))%IX(4))
              DO WHILE (.NOT.FOUND .AND. K .LT. Kmax)
                K = K + 1
                IF (I2 .EQ. CONTACT_SURFACE(NEIGHBOR(M))%IX(K)) THEN
                  IF (NEIGHBOR(M) .NE. N) THEN
                    CONTACT_SURFACE(N)%NX(i) = NEIGHBOR(M)
                    FOUND = .TRUE.
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            i = j
          ENDDO
        ENDDO
      ENDDO
!!
!! ### FIND INITIAL SLIDING NODE-CONTACT ELEMENT PAIRINGS
!! Find contact elements opposite sliding nodal points.
!!
      DO 100 Nsi = 1,NUMSI
!!
!! Skip over this sliding interface if it is a single-surface interface.
!!
        IF (SLIDING_INTERFACE(Nsi)%Type .EQ. 1) GO TO 100
!!
!! For sliding nodes on Side 1, look for contact elements on Side 2.
!!
        ISN1bgn = SLIDING_INTERFACE(Nsi)%ISN1bgn
        ISN1end = SLIDING_INTERFACE(Nsi)%ISN1end
!!
!! Loop over all sliding nodal points.
!!
        DO Nsn = ISN1bgn,ISN1end
          Xsn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Px
          Ysn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Py
          Zsn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Pz
          Dmin = HUGE( Dmin)
          SLIDING_NODE(Nsn)%Mce = 0
          ICE2bgn = SLIDING_INTERFACE(Nsi)%ICE2bgn
          ICE2end = SLIDING_INTERFACE(Nsi)%ICE2end
!!
!! Loop over all contact elements.
!!
          DO Mce = ICE2bgn,ICE2end
!!
!! Gather contact-element Mce nodal point coordinates.
!!
            DO k = 1,3+MIN(1,CONTACT_SURFACE(Mce)%IX(4))
       Xce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Px
       Yce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Py
       Zce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Pz
            ENDDO
!!
!! Compute contact-element Mce center coordinates.
!!
            IF (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0) THEN
              Xave = 0.3333 * (Xce(1) + Xce(2) + Xce(3))
              Yave = 0.3333 * (Yce(1) + Yce(2) + Yce(3))
              Zave = 0.3333 * (Zce(1) + Zce(2) + Zce(3))
              RR1 = (Xce(1)-Xave)**2+(Yce(1)-Yave)**2+(Zce(1)-Zave)**2
              RR2 = (Xce(2)-Xave)**2+(Yce(2)-Yave)**2+(Zce(2)-Zave)**2
              RR3 = (Xce(3)-Xave)**2+(Yce(3)-Yave)**2+(Zce(3)-Zave)**2
              Rsquared = MAX (RR1,RR2,RR3)
            ELSE
              Xave = 0.25 * (Xce(1) + Xce(2) + Xce(3) + Xce(4))
              Yave = 0.25 * (Yce(1) + Yce(2) + Yce(3) + Yce(4))
              Zave = 0.25 * (Zce(1) + Zce(2) + Zce(3) + Zce(4))
              RR1 = (Xce(1)-Xave)**2+(Yce(1)-Yave)**2+(Zce(1)-Zave)**2
              RR2 = (Xce(2)-Xave)**2+(Yce(2)-Yave)**2+(Zce(2)-Zave)**2
              RR3 = (Xce(3)-Xave)**2+(Yce(3)-Yave)**2+(Zce(3)-Zave)**2
              RR4 = (Xce(4)-Xave)**2+(Yce(4)-Yave)**2+(Zce(4)-Zave)**2
              Rsquared = MAX (RR1,RR2,RR3,RR4)
            ENDIF
!!
!! Compute distance between contact-element Mce and sliding nodal point Nsn.
!!
            Dsquared = (Xave-Xsn)**2 + (Yave-Ysn)**2 + (Zave-Zsn)**2
!!
!! Check to see if this contact element is closer than last contact element.
!!
            IF (Dsquared .LT. Rsquared) THEN
              IF (Dsquared .LT. Dmin) THEN
                Dmin = Dsquared
                SLIDING_NODE(Nsn)%Mce = Mce
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!!
!! For sliding nodes on Side 2, look for contact elements on Side 1.
!!
        ISN2bgn = SLIDING_INTERFACE(Nsi)%ISN2bgn
        ISN2end = SLIDING_INTERFACE(Nsi)%ISN2end
!!
!! Loop over all sliding nodal points.
!!
        DO Nsn = ISN2bgn,ISN2end
          Xsn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Px
          Ysn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Py
          Zsn = MOTION(CONTACT_NODE(SLIDING_NODE(Nsn)%Nsn)%NPID)%Pz
          Dmin = HUGE( Dmin )
          SLIDING_NODE(Nsn)%Mce = 0
          ICE1bgn = SLIDING_INTERFACE(Nsi)%ICE1bgn
          ICE1end = SLIDING_INTERFACE(Nsi)%ICE1end
!!
!! Loop over all contact elements.
!!
          DO Mce = ICE1bgn,ICE1end
!!
!! Gather contact-element Mce nodal point coordinates.
!!
            DO k = 1,3+MIN(1,CONTACT_SURFACE(Mce)%IX(4))
       Xce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Px
       Yce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Py
       Zce(k) = MOTION(CONTACT_NODE(CONTACT_SURFACE(Mce)%IX(k))%NPID)%Pz
            ENDDO
!!
!! Compute contact-element Mce center coordinates.
!!
            IF (CONTACT_SURFACE(Mce)%IX(4) .EQ. 0) THEN
              Xave = 0.3333 * (Xce(1) + Xce(2) + Xce(3))
              Yave = 0.3333 * (Yce(1) + Yce(2) + Yce(3))
              Zave = 0.3333 * (Zce(1) + Zce(2) + Zce(3))
              RR1 = (Xce(1)-Xave)**2+(Yce(1)-Yave)**2+(Zce(1)-Zave)**2
              RR2 = (Xce(2)-Xave)**2+(Yce(2)-Yave)**2+(Zce(2)-Zave)**2
              RR3 = (Xce(3)-Xave)**2+(Yce(3)-Yave)**2+(Zce(3)-Zave)**2
              Rsquared = MAX (RR1,RR2,RR3)
            ELSE
              Xave = 0.25 * (Xce(1) + Xce(2) + Xce(3) + Xce(4))
              Yave = 0.25 * (Yce(1) + Yce(2) + Yce(3) + Yce(4))
              Zave = 0.25 * (Zce(1) + Zce(2) + Zce(3) + Zce(4))
              RR1 = (Xce(1)-Xave)**2+(Yce(1)-Yave)**2+(Zce(1)-Zave)**2
              RR2 = (Xce(2)-Xave)**2+(Yce(2)-Yave)**2+(Zce(2)-Zave)**2
              RR3 = (Xce(3)-Xave)**2+(Yce(3)-Yave)**2+(Zce(3)-Zave)**2
              RR4 = (Xce(4)-Xave)**2+(Yce(4)-Yave)**2+(Zce(4)-Zave)**2
              Rsquared = MAX (RR1,RR2,RR3,RR4)
            ENDIF
!!
!! Compute distance between contact-element Mce and sliding nodal point Nsn.
!!
            Dsquared = (Xave-Xsn)**2 + (Yave-Ysn)**2 + (Zave-Zsn)**2
!!
!! Check to see if sliding nodal point Nsn is in the neighborhood of the
!! contact-element Mce
!!
            IF (Dsquared .LT. Rsquared) THEN
              IF (Dsquared .LT. Dmin) THEN
                Dmin = Dsquared
                SLIDING_NODE(Nsn)%Mce = Mce
              ENDIF
            ENDIF
          ENDDO
        ENDDO
 100    ENDDO

      RETURN
      END
!!_
      SUBROUTINE SPOT_WELD_INITIALIZATION (NSW)
!!
!! Copyright (c) by KEY Associates;  9-JUL-1994 10:45:12.00
!!
!! Purpose: Initialize spot weld data in order to minimize the com-
!! putational effort to impose the kinematic constraints a spot weld
!! represents. Spot welds can be used to tie together the motion of
!! one, two or three pieces of sheet metal at a point. (A spot weld
!! "fastening together one piece of sheet metal" is an oxymoron and
!! a degenerate case of questionable utility. However, the coding is
!! designed to support this case.)
!!
!! Spot welds maybe specified in terms of "coincident" nodal points
!! on adjacent pieces of sheet metal fastened together, or may be
!! specified in terms of a geometric location interior to one or more
!! elements on adjacent pieces of sheet metal fastened together.
!! The initilization and implementation for these two cases are
!! different.
!!
      USE shared_common_data
      USE platt_
      USE platq_
      USE node_
      USE motion_
      USE spot_weld_

      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      INTEGER, PARAMETER :: LSRT = 12 ! Length of sorted-element-distance array

      INTEGER                    &
     &          NSW,             & ! I/- Current spot weld on which to work    &
     &          NP(3),           & ! -/- Scratch array for spot welded NP's    &
     &          EL(3),           & ! -/- Scratch array for spot welded EL's    &
     &          NELEM(LSRT),     & ! -/- Sort array for closest el's, element n&
     &          IETYP(LSRT)        ! -/- Sort array for closest el's, 0/1=tri/qu
      REAL(KIND(0D0))                                                          &
     &          Xi1EL(LSRT),     & ! -/- Scratch array for EL iso-coordinates  &
     &          Xi2EL(LSRT),     & ! -/- Scratch array for EL iso-coordinates  &
     &          DSQRD(LSRT)        ! -/- Sort array for closest el's, d**2
      LOGICAL                                                                  &
     &          COINCIDENT,                                                    &
     &          INSIDE_ELEMENT
!!
!! Initialize failure flag to false (not failed).
!!
      SPOT_WELD(NSW)%FAILED = .FALSE.
!!
!! CASE 1: "POSITION"
!! Spot welds located by position must have the "position" converted
!! to the elements containing the spot weld and the isoparametric
!! position within each element identified. Both quadrilateral and
!! triangular elements are supported.
!!
      IF (SPOT_WELD(NSW)%PLACE .EQ. 'POSITION') THEN
!!
!! Construct distance between spot weld position and the center of each
!! quadrilateral and triangular element. First, retrieve spot weld position.
!!
        Px = SPOT_WELD(NSW)%Px
        Py = SPOT_WELD(NSW)%Py
        Pz = SPOT_WELD(NSW)%Pz
!!
!! Initialize arrays for search for the LSRT-closest elements. (IETYP(*)
!! equal to "10" is a value other than "0" or "1")
!!
        DO n = 1,LSRT
          DSQRD(n) = HUGE( DSQRD(n) )
          NELEM(n) = 0
          IETYP(n) = 10
        ENDDO
!!
!! Loop over all triangular elements (IETYP(*) = 0).
!!
        DO i = 1,NUMP3
          NP1 = PLATT(i)%PAR%IX(1)
          NP2 = PLATT(i)%PAR%IX(2)
          NP3 = PLATT(i)%PAR%IX(3)

          Xctr = 0.3333 * (MOTION(NP1)%Px+MOTION(NP2)%Px+MOTION(NP3)%Px)
          Yctr = 0.3333 * (MOTION(NP1)%Py+MOTION(NP2)%Py+MOTION(NP3)%Py)
          Zctr = 0.3333 * (MOTION(NP1)%Pz+MOTION(NP2)%Pz+MOTION(NP3)%Pz)

          Distance_Squared = (Xctr-Px)**2 + (Yctr-Py)**2 + (Zctr-Pz)**2

          n = 1
          DO WHILE (Distance_Squared .GT. DSQRD(n) .AND. n .LE. LSRT)
            n = n + 1
          ENDDO
          IF (n .LT. LSRT) THEN
            DO k = LSRT-1,n,-1
              DSQRD(k+1) = DSQRD(k)
              NELEM(k+1) = NELEM(k)
              IETYP(k+1) = IETYP(k)
            ENDDO
          ENDIF
          IF (n .LE. LSRT) THEN
            DSQRD(n) = Distance_Squared
            NELEM(n) = i
            IETYP(n) = 0
          ENDIF

        ENDDO
!!
!! Loop over all quadrilateral elements (IETYP(*) = 1).
!!
        DO i = 1,NUMP4
          NP1 = PLATQ(i)%PAR%IX(1)
          NP2 = PLATQ(i)%PAR%IX(2)
          NP3 = PLATQ(i)%PAR%IX(3)
          NP4 = PLATQ(i)%PAR%IX(4)

          Xctr = 0.25 * (MOTION(NP1)%Px + MOTION(NP2)%Px +                     &
     &                   MOTION(NP3)%Px + MOTION(NP4)%Px )
          Yctr = 0.25 * (MOTION(NP1)%Py + MOTION(NP2)%Py +                     &
     &                   MOTION(NP3)%Py + MOTION(NP4)%Py )
          Zctr = 0.25 * (MOTION(NP1)%Pz + MOTION(NP2)%Pz +                     &
     &                   MOTION(NP3)%Pz + MOTION(NP4)%Pz )

          Distance_Squared = (Xctr-Px)**2 + (Yctr-Py)**2 + (Zctr-Pz)**2

          n = 1
          DO WHILE (Distance_Squared .GT. DSQRD(n) .AND. n .LE. LSRT)
            n = n + 1
          ENDDO
          IF (n .LT. LSRT) THEN
            DO k = LSRT-1,n,-1
              DSQRD(k+1) = DSQRD(k)
              NELEM(k+1) = NELEM(k)
              IETYP(k+1) = IETYP(k)
            ENDDO
          ENDIF
          IF (n .LE. LSRT) THEN
            DSQRD(n) = Distance_Squared
            NELEM(n) = i
            IETYP(n) = 1
          ENDIF

        ENDDO
!!
!! Heading for informational output on search for isoparameteric locations
!! of spot weld location.
!!
        WRITE (IO_UNIT%LELO,'(/A/A,A)')                                        &
     &    ' *** INFORMATION *** SPOT_WELD_INITIALIZATION.001.01',              &
     &    ' *** INFORMATION *** ELEMENT  ID   ,Iterate,=====Xi1=====,',        &
     &    '=====Xi2=====,==d(RR)/dXi1=,==d(RR)/dXi2='
!!
!! Find the two (or three) elements that contain the spot weld position.
!! The in-plane isoparametric position must lie between -1 and +1, and
!! the distance between the "position" and the shell middle surface
!! must be a "minimum."
!!
        DO 300 n = 1,LSRT
          IEL = NELEM(n)
!!
!! Find Xi1 and Xi2 that minimizes r**2 = (Px - PxI*NI(Xi1,Xi2))**2 +
!! (Py - PyI*NI(Xi1,Xi2))**2 + (Pz - PzI*NI(Xi1,Xi2))**2
!!
!!         Xi2 ^
!!             |
!!          +1 * 3
!!             |\
!!             | \
!!             |  \
!!             | # \<---- Spot weld located at (Px,Py,Pz)
!!             |    \
!!             |1    \ 2
!!           0 *------*---> Xi1
!!             0     +1
!!
          IF (IETYP(n) .EQ. 0) THEN

            I1 = PLATT(IEL)%PAR%IX(1)
            I2 = PLATT(IEL)%PAR%IX(2)
            I3 = PLATT(IEL)%PAR%IX(3)

            X1 = MOTION(I1)%Px
            X2 = MOTION(I2)%Px
            X3 = MOTION(I3)%Px

            Y1 = MOTION(I1)%Py
            Y2 = MOTION(I2)%Py
            Y3 = MOTION(I3)%Py

            Z1 = MOTION(I1)%Pz
            Z2 = MOTION(I2)%Pz
            Z3 = MOTION(I3)%Pz
!!
!! Initial guess at isoparametric coordinates Xi1 and Xi2
!!
            Xi1 = 0.0
            Xi2 = 0.0
!!
!! Evaluate linear isoparameteric shape functions.
!!
            A1 = 1.0 - Xi1 - Xi2
            A2 = Xi1
            A3 = Xi2

            Qx = (A1*X1 + A2*X2 + A3*X3 - Px)
            Qy = (A1*Y1 + A2*Y2 + A3*Y3 - Py)
            Qz = (A1*Z1 + A2*Z2 + A3*Z3 - Pz)

            Qx1 = X2 - X1
            Qy1 = Y2 - Y1
            Qz1 = Z2 - Z1

            Qx2 = X3 - X1
            Qy2 = Y3 - Y1
            Qz2 = Z3 - Z1
!!
!! Compute slope of r**2 at current point (vanishes at point of
!! closest approach).
!!
            RR1 = 2.0 * (Qx*Qx1 + Qy*Qy1 + Qz*Qz1)
            RR2 = 2.0 * (Qx*Qx2 + Qy*Qy2 + Qz*Qz2)
!!
!! Compute Newton-Raphson gradient
!!
            RR11 = 2.0 * (Qx1*Qx1 + Qy1*Qy1 + Qz1*Qz1)
            RR12 = 2.0 * (Qx1*Qx2 + Qy1*Qy2 + Qz1*Qz2)
            RR22 = 2.0 * (Qx2*Qx2 + Qy2*Qy2 + Qz2*Qz2)
!!
!! Compute N-R iterative improvement. (In this case the gradient is
!! constant and the solution is obtained in one shot.)
!!
            DET = RR11 * RR22 - RR12 * RR12
            C11 = +RR22 / DET
            C22 = +RR11 / DET
            C12 = -RR12 / DET
            dXi1 = -(C11 * RR1 + C12 * RR2)
            dXi2 = -(C12 * RR1 + C22 * RR2)
            Xi1 = Xi1 + dXi1
            Xi2 = Xi2 + dXi2
!!
!! Compute squared distance between spot weld location and point
!! in element.
!!
            A1 = 1.0 - Xi1 - Xi2
            A2 = Xi1
            A3 = Xi2
            Qx = A1*X1 + A2*X2 + A3*X3 - Px
            Qy = A1*Y1 + A2*Y2 + A3*Y3 - Py
            Qz = A1*Z1 + A2*Z2 + A3*Z3 - Pz
            RR = Qx*Qx + Qy*Qy + Qz*Qz
            RR1 = 2.0 * (Qx*Qx1 + Qy*Qy1 + Qz*Qz1)
            RR2 = 2.0 * (Qx*Qx2 + Qy*Qy2 + Qz*Qz2)

            WRITE (IO_UNIT%LELO,'(21X,A,2I6,3X,4(1PE14.7))')                   &
     &          ' P3EL  ',PLATT(IEL)%PAR%EleID,1,Xi1,Xi2,RR1,RR2
!!
!! Modify LSRT-closest element data.
!!
            DSQRD(n) = RR
            Xi1EL(n) = Xi1
            Xi2EL(n) = Xi2
!!
!! QUADRILATERAL ELEMENT (IETYP(*)=1).
!! Find Xi1 and Xi2 that minimizes r**2 = (Px - PxI*NI(Xi1,Xi2))**2 +
!! (Py - PyI*NI(Xi1,Xi2))**2 + (Pz - PzI*NI(Xi1,Xi2))**2
!!
!!                  Xi2
!!            4     +1      3
!!             *-----|-----*
!!             |     |     |
!!             |     |     |
!!         -1 -------0------- +1 Xi1
!!             |     |     |
!!             |     |  #  |<---- Spot weld located at (Px,Py,Pz)
!!             *-----|-----*
!!            1     -1      2
!!
          ELSE IF (IETYP(n) .EQ. 1) THEN

            I1 = PLATQ(IEL)%PAR%IX(1)
            I2 = PLATQ(IEL)%PAR%IX(2)
            I3 = PLATQ(IEL)%PAR%IX(3)
            I4 = PLATQ(IEL)%PAR%IX(4)

            X1 = MOTION(I1)%Px
            X2 = MOTION(I2)%Px
            X3 = MOTION(I3)%Px
            X4 = MOTION(I4)%Px

            Y1 = MOTION(I1)%Py
            Y2 = MOTION(I2)%Py
            Y3 = MOTION(I3)%Py
            Y4 = MOTION(I4)%Py

            Z1 = MOTION(I1)%Pz
            Z2 = MOTION(I2)%Pz
            Z3 = MOTION(I3)%Pz
            Z4 = MOTION(I4)%Pz
!!
!! Initial guess at isoparametric coordinates Xi1 and Xi2
!!
            Xi1 = 0.0
            Xi2 = 0.0

            Icount = 0
 100          Icount = Icount + 1
!!
!! Evaluate bilinear isoparameteric shape functions.
!!
            A1 = 0.25*(1.0 - Xi1)*(1.0 - Xi2)
            A2 = 0.25*(1.0 + Xi1)*(1.0 - Xi2)
            A3 = 0.25*(1.0 + Xi1)*(1.0 + Xi2)
            A4 = 0.25*(1.0 - Xi1)*(1.0 + Xi2)
!!
!! Evaluate Xi1 derivative of isoparameteric shape functions
!!
            A11 = -0.25*(1.0 - Xi2)
            A21 = +0.25*(1.0 - Xi2)
            A31 = +0.25*(1.0 + Xi2)
            A41 = -0.25*(1.0 + Xi2)
!!
!! Evaluate Xi2 derivative of isoparameteric shape functions
!!
            A12 = -0.25*(1.0 - Xi1)
            A22 = -0.25*(1.0 + Xi1)
            A32 = +0.25*(1.0 + Xi1)
            A42 = +0.25*(1.0 - Xi1)
!!
!! Evaluate Xi1Xi2 derivative of isoparameteric shape functions
!!
            A112 = +0.25
            A212 = -0.25
            A312 = +0.25
            A412 = -0.25

            Qx = (A1*X1 + A2*X2 + A3*X3 + A4*X4 - Px)
            Qy = (A1*Y1 + A2*Y2 + A3*Y3 + A4*Y4 - Py)
            Qz = (A1*Z1 + A2*Z2 + A3*Z3 + A4*Z4 - Pz)

            Qx1 = (A11*X1 + A21*X2 + A31*X3 + A41*X4)
            Qy1 = (A11*Y1 + A21*Y2 + A31*Y3 + A41*Y4)
            Qz1 = (A11*Z1 + A21*Z2 + A31*Z3 + A41*Z4)

            Qx2 = (A12*X1 + A22*X2 + A32*X3 + A42*X4)
            Qy2 = (A12*Y1 + A22*Y2 + A32*Y3 + A42*Y4)
            Qz2 = (A12*Z1 + A22*Z2 + A32*Z3 + A42*Z4)

            Qx12 = (A112*X1 + A212*X2 + A312*X3 + A412*X4)
            Qy12 = (A112*Y1 + A212*Y2 + A312*Y3 + A412*Y4)
            Qz12 = (A112*Z1 + A212*Z2 + A312*Z3 + A412*Z4)
!!
!! Compute slope of r**2 at current point (vanishes at point of
!! closest approach).
!!
            RR1 = 2.0 * (Qx*Qx1 + Qy*Qy1 + Qz*Qz1)
            RR2 = 2.0 * (Qx*Qx2 + Qy*Qy2 + Qz*Qz2)
!!
!! Check for convergence.
!!
            IF (ABS (RR1) .LT. 1.0D-4 .AND. ABS (RR2) .LT. 1.0D-4) THEN
!!
!! Compute squared distance between spot weld location and point
!! in element.
!!
              RR = Qx*Qx + Qy*Qy + Qz*Qz
!!
!! Modify LSRT-closest element data.
!!
              DSQRD(n) = RR
              Xi1EL(n) = Xi1
              Xi2EL(n) = Xi2

              WRITE (IO_UNIT%LELO,'(21X,A,2I6,3X,4(1PE14.7))')                 &
     &          ' P4EL  ',PLATQ(IEL)%PAR%EleID,Icount,Xi1,Xi2,RR1,RR2
!!
!! Iterate to improve values of Xi1 and Xi2
!!
            ELSE IF (Icount .LE. 20) THEN
!!
!! Compute Newton-Raphson gradient
!!
              RR11 = 2.0 * (Qx1*Qx1 + Qy1*Qy1 + Qz1*Qz1)
              RR12 = 2.0 * (Qx1*Qx2 + Qy1*Qy2 + Qz1*Qz2 +                      &
     &                      Qx*Qx12 + Qy*Qy12 + Qz*Qz12)
              RR22 = 2.0 * (Qx2*Qx2 + Qy2*Qy2 + Qz2*Qz2)
!!
!! Compute N-R iterative improvement.
!!
              DET = RR11 * RR22 - RR12 * RR12
              C11 = +RR22 / DET
              C22 = +RR11 / DET
              C12 = -RR12 / DET
              dXi1 = -(C11 * RR1 + C12 * RR2)
              dXi2 = -(C12 * RR1 + C22 * RR2)
              Xi1 = Xi1 + dXi1
              Xi2 = Xi2 + dXi2

              WRITE (IO_UNIT%LELO,'(21X,A,2I6,3X,4(1PE14.7))')                 &
     &          ' P4EL  ',PLATQ(IEL)%PAR%EleID,Icount,Xi1,Xi2,RR1,RR2

              GO TO 100
!!
!! Too many N-R iterations.
!!
            ELSE
              ERROR%COUNT = ERROR%COUNT + 1
              WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
              WRITE (MSG2,'(I8)') PLATQ(IEL)%PAR%EleID
              CALL USER_MESSAGE                                                &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SPOT_WELD_INITIALIZATION.002.01'//                      &
     &          MSGL//'SPOTWELD Input Record ID:'//MSG1//                      &
     &          MSGL//'Location WRT Element ID: '//MSG2//                      &
     &          MSGL//'Iso-Coordinate Iterations Exceed 20.'                   &
     &          )
            ENDIF

          ENDIF
 300      ENDDO
!!
!! Select elements that are spot_welded together.
!!
        Icount = 0
        DO i = 1,3
          Nmin = 0
          Dmin = HUGE( Dmin )
          DO n = 1,LSRT
            IF (IETYP(n) .EQ. 0) THEN
              INSIDE_ELEMENT =                                                 &
     &          Xi1EL(n) .GE. 0.0 .AND. Xi2EL(n) .GE. 0.0
              INSIDE_ELEMENT =                                                 &
     &          INSIDE_ELEMENT .AND. Xi1EL(n)+Xi2EL(n).LE.1.0
            ELSE IF (IETYP(n) .EQ. 1) THEN
              INSIDE_ELEMENT =                                                 &
     &          ABS(Xi1EL(n)).LE.1.0 .AND. ABS(Xi2EL(n)).LE.1.0
            ENDIF
            IF (INSIDE_ELEMENT) THEN
              IF (DSQRD(n) .LT. Dmin) THEN
              Dmin = DSQRD(n)
              Nmin = n
              ENDIF
            ENDIF
          ENDDO
          IF (Nmin .GT. 0) THEN
            Icount = Icount + 1
            SPOT_WELD(NSW)%EleID(i) = Nmin
            SPOT_WELD(NSW)%Type(i)  = IETYP(Nmin)
            SPOT_WELD(NSW)%Xi1(i)   = Xi1EL(Nmin)
            SPOT_WELD(NSW)%Xi2(i)   = Xi2EL(Nmin)
!!
!! Jigger isoparametric coordinate so that element is skipped in
!! succeeding searches.
!!
            Xi1EL(Nmin) = 2.0
          ENDIF
        ENDDO
!!
!! Report the elements spot welded together.
!!
        IF (Icount .GT. 0) THEN
          WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
          DO i = 1,Icount
            IEL = SPOT_WELD(NSW)%EleID(i)
            IF (SPOT_WELD(NSW)%Type(i) .EQ. 0) THEN
              WRITE (MSG2,'(I8)') PLATT(IEL)%PAR%EleID
            ELSE
              WRITE (MSG2,'(I8)') PLATQ(IEL)%PAR%EleID
            ENDIF
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'SPOT_WELD_INITIALIZATION.001.02'//                      &
     &          MSGL//'SPOTWELD Input Record ID:'//MSG1//                      &
     &          MSGL//'Spot Welded Element ID:  '//MSG2                        &
     &          )
          ENDDO
!!
!! Report no-finds.
!!
        ELSE
          ERROR%COUNT = ERROR%COUNT + 1
          WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SPOT_WELD_INITIALIZATION.002.02'//                      &
     &          MSGL//'SPOTWELD Input Record ID:'//MSG1//                      &
     &          MSGL//'No Elements Were Found At User''s Location.'            &
     &          )
        ENDIF
!!
!! CASE 2: "NODES"
!!
      ELSE IF (SPOT_WELD(NSW)%PLACE .EQ. 'NODES') THEN

        NP1 = SPOT_WELD(NSW)%NPID(1)
        NP2 = SPOT_WELD(NSW)%NPID(2)
        NP3 = SPOT_WELD(NSW)%NPID(3)
!!
!! First, compute the distance between spot welded nodal points.
!!
        D12 = 0.0
        D13 = 0.0
        D23 = 0.0

        IF (NP1.GT.0 .AND. NP2.GT.0)                                           &
     &          D12 = (MOTION(NP1)%Px - MOTION(NP2)%Px)**2                     &
     &              + (MOTION(NP1)%Py - MOTION(NP2)%Py)**2                     &
     &              + (MOTION(NP1)%Pz - MOTION(NP2)%Pz)**2

        IF (NP1.GT.0 .AND. NP3.GT.0)                                           &
     &          D13 = (MOTION(NP1)%Px - MOTION(NP3)%Px)**2                     &
     &              + (MOTION(NP1)%Py - MOTION(NP3)%Py)**2                     &
     &              + (MOTION(NP1)%Pz - MOTION(NP3)%Pz)**2

        IF (NP2.GT.0 .AND. NP3.GT.0)                                           &
     &          D23 = (MOTION(NP2)%Px - MOTION(NP3)%Px)**2                     &
     &              + (MOTION(NP2)%Py - MOTION(NP3)%Py)**2                     &
     &              + (MOTION(NP2)%Pz - MOTION(NP3)%Pz)**2
!!
!! Second, check to see if the spot welded nodal points are coincident.
!!
        COINCIDENT = (D12 .LT. 1.0D-3)                                         &
     &           .AND. (D13 .LT. 1.0D-3)                                       &
     &           .AND. (D23 .LT. 1.0D-3)

        SPOT_WELD(NSW)%Coincident = COINCIDENT
!!
!! Count the number of nodes spot welded together. If there are three
!! non-coincident nodes spot welded together (rather than just two nodes
!! spot welded together), make the intermediate node NP2, and place the
!! intermediate node on a line between NP1 and NP3 (disturbing the spacing
!! as little as possible).
!!
        IF                                                                     &
     &          (                                                              &
     &          .NOT.COINCIDENT                                                &
     &          .AND.                                                          &
     &          (NP1.GT.0 .AND. NP2.GT.0 .AND. NP3.GT.0)                       &
     &          )                                                              &
     &    THEN
!!
!! Re-order nodal points.
!!
          Dmax = MAX (D12, D13, D23)
          IF (Dmax .EQ. D12) THEN
            SPOT_WELD(NSW)%NPID(2) = NP3
            SPOT_WELD(NSW)%NPID(3) = NP2
          ELSE IF (Dmax .EQ. D23) THEN
            SPOT_WELD(NSW)%NPID(1) = NP2
            SPOT_WELD(NSW)%NPID(2) = NP1
          ENDIF
!!
!! Reposition intermediate nodal point.
!!
          NP1 = SPOT_WELD(NSW)%NPID(1)
          NP2 = SPOT_WELD(NSW)%NPID(2)
          NP3 = SPOT_WELD(NSW)%NPID(3)

          D12 = (MOTION(NP1)%Px - MOTION(NP2)%Px)**2                           &
     &        + (MOTION(NP1)%Py - MOTION(NP2)%Py)**2                           &
     &        + (MOTION(NP1)%Pz - MOTION(NP2)%Pz)**2
          D23 = (MOTION(NP2)%Px - MOTION(NP3)%Px)**2                           &
     &        + (MOTION(NP2)%Py - MOTION(NP3)%Py)**2                           &
     &        + (MOTION(NP2)%Pz - MOTION(NP3)%Pz)**2

          A12 = SQRT (D12)
          A23 = SQRT (D23)
          W1 = A12 / (A12+A23)
          W3 = A23 / (A12+A23)
          MOTION(NP2)%Px = W3 * MOTION(NP1)%Px + W1 * MOTION(NP3)%Px
          MOTION(NP2)%Py = W3 * MOTION(NP1)%Py + W1 * MOTION(NP3)%Py
          MOTION(NP2)%Pz = W3 * MOTION(NP1)%Pz + W1 * MOTION(NP3)%Pz
!!
!! Since the spot welded nodal points are "coincident," co-locate them.
!!
        ELSE IF (COINCIDENT) THEN

          Ncount = 0
          Px = 0.0
          Py = 0.0
          Pz = 0.0
          IF (NP1 .GT. 0) THEN
            Ncount = Ncount + 1
            Px = Px + MOTION(NP1)%Px
            Py = Py + MOTION(NP1)%Py
            Pz = Pz + MOTION(NP1)%Pz
          ENDIF
          IF (NP2 .GT. 0) THEN
            Ncount = Ncount + 1
            Px = Px + MOTION(NP2)%Px
            Py = Py + MOTION(NP2)%Py
            Pz = Pz + MOTION(NP2)%Pz
          ENDIF
          IF (NP3 .GT. 0) THEN
            Ncount = Ncount + 1
            Px = Px + MOTION(NP3)%Px
            Py = Py + MOTION(NP3)%Py
            Pz = Pz + MOTION(NP3)%Pz
          ENDIF
          IF (Ncount .EQ. 0) THEN
            ERROR%COUNT = ERROR%COUNT + 1
            WRITE (MSG1,'(I8)') SPOT_WELD(NSW)%SWID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SPOT_WELD_INITIALIZATION.002.03'//                      &
     &          MSGL//'SPOTWELD Input Record ID:'//MSG1//                      &
     &          MSGL//'Has Zero Nodal Points Specified.'                       &
     &          )
          ELSE
            Px = Px / DBLE (Ncount)
            Py = Py / DBLE (Ncount)
            Pz = Pz / DBLE (Ncount)
            IF (NP1 .GT. 0) THEN
              MOTION(NP1)%Px = Px
              MOTION(NP1)%Py = Py
              MOTION(NP1)%Pz = Pz
            ENDIF
            IF (NP2 .GT. 0) THEN
              MOTION(NP2)%Px = Px
              MOTION(NP2)%Py = Py
              MOTION(NP2)%Pz = Pz
            ENDIF
            IF (NP3 .GT. 0) THEN
              MOTION(NP3)%Px = Px
              MOTION(NP3)%Py = Py
              MOTION(NP3)%Pz = Pz
            ENDIF
          ENDIF
        ENDIF
!!
!! Compute "dumb bell" center of mass and inertia.
!!
        IF (.NOT.COINCIDENT) THEN

          NP(1) = SPOT_WELD(NSW)%NPID(1)
          NP(2) = SPOT_WELD(NSW)%NPID(2)
          NP(3) = SPOT_WELD(NSW)%NPID(3)

          Px = 0.0
          Py = 0.0
          Pz = 0.0
          Qm = 0.0
          DO i = 1,3
            IF (NP(i) .GT. 0) THEN
              NTR = NP(i)
              Px = Px + NODE(NTR)%Mass * MOTION(NTR)%Px
              Py = Py + NODE(NTR)%Mass * MOTION(NTR)%Py
              Pz = Pz + NODE(NTR)%Mass * MOTION(NTR)%Pz
              Qm = Qm + NODE(NTR)%Mass
            ENDIF
          ENDDO
          Xcm = Px / Qm
          Ycm = Py / Qm
          Zcm = Pz / Qm
!!
!! Store "dumb bell" center of mass.
!!
          SPOT_WELD(NSW)%Xcm = Xcm
          SPOT_WELD(NSW)%Ycm = Ycm
          SPOT_WELD(NSW)%Zcm = Zcm
!!
!! Compute "dumb bell" inertia
!!
          Brr = 0.0
          DO i = 1,3
            IF (NP(i) .GT. 0) THEN
              NTR = NP(i)
              RSQ = (MOTION(NTR)%Px - Xcm)**2                                  &
     &            + (MOTION(NTR)%Py - Ycm)**2                                  &
     &            + (MOTION(NTR)%Pz - Zcm)**2
              Brr = Brr + RSQ * NODE(NTR)%Mass
            ENDIF
          ENDDO
          SPOT_WELD(NSW)%Brr = Brr
        ENDIF

      ENDIF

      RETURN
      END
!!_
      SUBROUTINE SECTION_PROPERTIES (NS1)
!!
!! Copyright (c) by KEY Associates; 19-MAY-1991 10:15:48
!!
!! Purpose: Calculate the undeformed area of truss and beam geometries,
!! and calculate the area monents of inertia for the cross section.
!!
      USE shared_common_data
      USE section_1d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NS1  ! Current section on which to work
!!
!! Locall variables.
      REAL(KIND(0D0))                                                          &
     &          AESS,AEHS,ARSS,ARHS,AIBS,     & ! Area of section              &
     &          BESr,BEHr,BRSr,BRHr,BIBr,     & ! I(r**2)da inertia            &
     &          BESy,BEHy,BRSy,BRHy,BIBy,     & ! I(y**2)da inertia            &
     &          BESz,BEHz,BRSz,BRHz,BIBz        ! I(z**2)da inertia
!!
      CHARACTER, PARAMETER :: &
     &  GEOMETRY(3)*12 = (/' Elliptical ',' Rectangular',' I-Section  '/)
!!
      REAL(KIND(0D0)), PARAMETER :: &
     &  DZ=0.0D+0, D1=1.0D+0, D2=2.0D+0, D3=3.0D+0, D4=4.0D+0, &
     &  D5=5.0D+0, D6=6.0D+0, D7=7.0D+0, D8=8.0D+0, D9=9.0D+0
!!
!! Pi over 4.0
!!
      REAL(KIND(0D0)), PARAMETER :: &
     &  QPi = (3.1415926535897932384626433D+0/4.0D+0)
!!
!! Convergence tolerance when back calculating section dimensions.
!!
      REAL(KIND(0D0)), PARAMETER :: Epsilon = 1.0D-6
!!
!! Define function statements for various truss and beam cross sections.
!!
!! 1a. Elliptical solid section of width w and heigth h. (0.0625 = 1/16)
!!
      AESS(w,h) = QPi*w*h
      BESr(w,h) = AESS(w,h) * (w*w+h*h) * 0.0625D+0
      BESy(w,h) = AESS(w,h) * (w*w)     * 0.0625D+0
      BESz(w,h) = AESS(w,h) * (h*h)     * 0.0625D+0
!!
!! 1b. Elliptical hollow section of width w, heigth h, and wall t.
!!
      AEHS(w,h,t) = AESS(w,h) - AESS(w-t-t,h-t-t)
      BEHr(w,h,t) = BESr(w,h) - BESr(w-t-t,h-t-t)
      BEHy(w,h,t) = BESy(w,h) - BESy(w-t-t,h-t-t)
      BEHz(w,h,t) = BESz(w,h) - BESz(w-t-t,h-t-t)
!!
!! 2a. Solid rectangle of width w and heigth h. (0.083333333= 1/12)
!!
      ARSS(w,h) = w*h
      BRSr(w,h) = ARSS(w,h) * (w*w+h*h) * 0.0833333333D+0
      BRSy(w,h) = ARSS(w,h) * (w*w)     * 0.0833333333D+0
      BRSz(w,h) = ARSS(w,h) * (h*h)     * 0.0833333333D+0
!!
!! 2b. Hollow rectangle of width w, heigth h, and wall t.
!!
      ARHS(w,h,t) = ARSS(w,h) - ARSS(w-t-t,h-t-t)
      BRHr(w,h,t) = BRSr(w,h) - BRSr(w-t-t,h-t-t)
      BRHy(w,h,t) = BRSy(w,h) - BRSy(w-t-t,h-t-t)
      BRHz(w,h,t) = BRSz(w,h) - BRSz(w-t-t,h-t-t)
!!
!! 3. I-section of width w, heigth h, web tw and flange tf.
!!
      AIBS(w,h,tw,tf) = ARSS(w,h) - ARSS(w,h-tf-tf) + ARSS(tw,h-tf-tf)
      BIBr(w,h,tw,tf) = BRSr(w,h) - BRSr(w,h-tf-tf) + BRSr(tw,h-tf-tf)
      BIBy(w,h,tw,tf) = BRSy(w,h) - BRSy(w,h-tf-tf) + BRSy(tw,h-tf-tf)
      BIBz(w,h,tw,tf) = BRSz(w,h) - BRSz(w,h-tf-tf) + BRSz(tw,h-tf-tf)
!!
!! Newton-Raphson gradient operator for both elliptical and rectangular
!! cross sections.
!!
      DADW(t) = t+t
      DADH(t) = t+t
      DADT(w,h,t) = D2*(w+h-t-t-t-t)
      DYDW(w,h,t) = D3*(w*w*h-(w-t-t)*(w-t-t)*(h-t-t))
      DYDH(w,h,t) = w*w*w-(w-t-t)*(w-t-t)*(w-t-t)
      DYDT(w,h,t) = D6*(w-t-t)*(w-t-t)*(h-t-t)+D2*(w-t-t)*(w-t-t)*(w-t-t)
      DZDW(w,h,t) = h*h*h-(h-t-t)*(h-t-t)*(h-t-t)
      DZDH(w,h,t) = D3*(h*h*w-(h-t-t)*(h-t-t)*(w-t-t))
      DZDT(w,h,t) = D6*(h-t-t)*(h-t-t)*(w-t-t)+D2*(h-t-t)*(h-t-t)*(h-t-t)
!!
!! Retrieve cross-section parameters and compute area and inertia properties.
!!
      IF (SECTION_1D(NS1)%Iprop .EQ. 0) THEN
        Isec    = SECTION_1D(NS1)%Section
        Width   = SECTION_1D(NS1)%Width
        Height  = SECTION_1D(NS1)%Height
        Twall   = SECTION_1D(NS1)%Twall
        Tflange = SECTION_1D(NS1)%Tflange
        IF (Isec .EQ. 1) THEN                 ! ROD, solid ellipse/circle
          Area = AESS(Width,Height)
          Br   = BESr(Width,Height)
          By   = BESy(Width,Height)
          Bz   = BESz(Width,Height)
        ELSE IF (Isec .EQ. 2) THEN            ! TUBE, hollow ellipse/circle
          Area = AEHS(Width,Height,Twall)
          Br   = BEHr(Width,Height,Twall)
          By   = BEHy(Width,Height,Twall)
          Bz   = BEHz(Width,Height,Twall)
        ELSE IF (Isec .EQ. 3) THEN            ! BAR, solid rectangle/square
          Area = ARSS(Width,Height)
          Br   = BRSr(Width,Height)
          By   = BRSy(Width,Height)
          Bz   = BRSz(Width,Height)
        ELSE  IF (Isec .EQ. 4) THEN           ! BOX, hollow rectangle/square
          Area = ARHS(Width,Height,Twall)
          Br   = BRHr(Width,Height,Twall)
          By   = BRHy(Width,Height,Twall)
          Bz   = BRHz(Width,Height,Twall)
        ELSE IF (Isec .EQ. 5) THEN            ! I, I-beam cross section
          Area = AIBS(Width,Height,Twall,Tflange)
          Br   = BIBr(Width,Height,Twall,Tflange)
          By   = BIBy(Width,Height,Twall,Tflange)
          Bz   = BIBz(Width,Height,Twall,Tflange)
        ENDIF
!!
!! Check to make certain calculated area of truss or beam is non-zero.
!!
        IF (Area .EQ. DZ) THEN
          WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SECTION_PROPERTIES.001.00'//                            &
     &          MSGL//'1-D (BSECTION) Section ID:'//MSG1//                     &
     &          MSGL//'Section Geometry:'//GEOMETRY(Isec)//                    &
     &          MSGL//'Has A Zero Cross-Sectional Area.'                       &
     &          )
          ERROR%COUNT = ERROR%COUNT + 1
        ELSE
          SECTION_1D(NS1)%Area = Area
          SECTION_1D(NS1)%Br   = Br
          SECTION_1D(NS1)%By   = By
          SECTION_1D(NS1)%Bz   = Bz
        ENDIF
!!
!! Retrieve cross-section area and inertia properties and compute the section
!! dimensions.
!!
      ELSE IF (SECTION_1D(NS1)%Iprop .NE. 0) THEN
        Isec = SECTION_1D(NS1)%Section
        Area = SECTION_1D(NS1)%Area
        By   = SECTION_1D(NS1)%By
        Bz   = SECTION_1D(NS1)%Bz
        BCSMXIT = PARAMVALUE%BCSMXIT
        BCSRLAX = PARAMVALUE%BCSRLAX
        GRLMmax = D1 + PARAMVALUE%BCSGRLM
        GRLMmin = D1 - PARAMVALUE%BCSGRLM
!!
!! Define the torsional inertia in terms of the other two cross section
!! inertia's. (This expression is an identity or invariant relation that
!! is required to be true.)
!!
        Br = By + Bz
        SECTION_1D(NS1)%Br = Br
!!
!! ROD, solid ellipse/circle.
!!
        IF (Isec .EQ. 1) THEN
          w = D4 * SQRT (By/Area)
          h = D4 * SQRT (Bz/Area)
          IF (ABS(Area-AESS(w,h)) .GT. Epsilon*Area) THEN
            WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SECTION_PROPERTIES.002.01'//                            &
     &          MSGL//'BSECTION Input Record ID:'//MSG1//                      &
     &          MSGL//'A ROD Width/Height-pair Does Not Exist'//               &
     &          MSGL//'That Reproduces The Input Cross Section Moments.'       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            SECTION_1D(NS1)%Width  = w
            SECTION_1D(NS1)%Height = h
            SECTION_1D(NS1)%Twall  = DZ
          ENDIF
!!
!! TUBE, hollow ellipse/circle.
!!
        ELSE IF (Isec .EQ. 2) THEN
!!
!! Iterate for cross section dimensions that reproduce input section properties.
!!
          d = SQRT(Area/0.36D+0)
          Ratio = (By/Bz)**(0.33333333333D+0)
          w = d * Ratio
          h = d / Ratio
          t = d * 0.1D+0
          Icount = 0
          Imaxit = NINT (BCSMXIT)
          C1 = AEHS(w,h,t) - Area
          C2 = BEHy(w,h,t) - By
          C3 = BEHz(w,h,t) - Bz
          Residual = ABS(C1) + ABS(C2) + ABS(C3)

          WRITE (IO_UNIT%LELO,'(/2X,A,I8,3X,A/1X,A,A)')                        &
     &          'Beam Section (BSECTION) ID:',SECTION_1D(NS1)%SecID,           &
     &          '(TUBE -- Hollow Ellipse/Circle)',                             &
     &          'Iteration Count,===Residual==,====Width====,',                &
     &          '====Height===,====Twall===='

          DO WHILE (Residual .GT. Epsilon*Area .AND. Icount .LE. Imaxit)

            WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                       &
     &        Icount,Residual,w,h,t

            A11 = QPi*DADW(t)
            A12 = QPi*DADH(t)
            A13 = QPi*DADT(w,h,t)
            A21 = QPi*DYDW(w,h,t)*0.0625D+0
            A22 = QPi*DYDH(w,h,t)*0.0625D+0
            A23 = QPi*DYDT(w,h,t)*0.0625D+0
            A31 = QPi*DZDW(w,h,t)*0.0625D+0
            A32 = QPi*DZDH(w,h,t)*0.0625D+0
            A33 = QPi*DZDT(w,h,t)*0.0625D+0
            DET = A11*A22*A33+A12*A23*A31+A13*A21*A32                          &
     &          - A31*A22*A13-A32*A23*A11-A33*A21*A12
            B11 = (A22*A33-A32*A23)*(D1/DET)
            B12 = (A32*A13-A12*A33)*(D1/DET)
            B13 = (A12*A23-A13*A22)*(D1/DET)
            B21 = (A31*A23-A21*A33)*(D1/DET)
            B22 = (A11*A33-A31*A13)*(D1/DET)
            B23 = (A21*A13-A11*A23)*(D1/DET)
            B31 = (A21*A32-A31*A22)*(D1/DET)
            B32 = (A31*A12-A11*A32)*(D1/DET)
            B33 = (A22*A11-A12*A21)*(D1/DET)
            x = w - BCSRLAX * (B11*C1 + B12*C2 + B13*C3)
            y = h - BCSRLAX * (B21*C1 + B22*C2 + B23*C3)
            z = t - BCSRLAX * (B31*C1 + B32*C2 + B33*C3)
            w = MAX (GRLMmin*w, MIN (GRLMmax*w, x))
            h = MAX (GRLMmin*h, MIN (GRLMmax*h, y))
            t = MAX (GRLMmin*t, MIN (GRLMmax*t, z))
            Icount = Icount + 1
            C1 = AEHS(w,h,t) - Area
            C2 = BEHy(w,h,t) - By
            C3 = BEHz(w,h,t) - Bz
            Residual = ABS(C1) + ABS(C2) + ABS(C3)
          ENDDO

          WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                         &
     &      Icount,Residual,w,h,t

          IF (ABS(Area-AEHS(w,h,t)) .GT. Epsilon*Area) THEN
            WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
            CALL USER_MESSAGE                                                  &
     &      (                                                                  &
     &      MSGL//'WARN'//                                                     &
     &      MSGL//'SECTION_PROPERTIES.002.02'//                                &
     &      MSGL//'BSECTION Input Record ID:'//MSG1//                          &
     &      MSGL//'A TUBE Width/Height/Twall-triple Does Not Exist'//          &
     &      MSGL//'That Reproduces The Input Cross Section Moments.'           &
     &      )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            SECTION_1D(NS1)%Width  = w
            SECTION_1D(NS1)%Height = h
            SECTION_1D(NS1)%Twall  = t
          ENDIF
!!
!! BAR, solid rectangle/square.
!!
        ELSE IF (Isec .EQ. 3) THEN
          w = 3.464101615D+0 * SQRT (By/Area)
          h = 3.464101615D+0 * SQRT (Bz/Area)
          IF (ABS(Area-ARSS(w,h)) .GT. Epsilon*Area) THEN
            WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SECTION_PROPERTIES.002.03'//                            &
     &          MSGL//'BSECTION Input Record ID:'//MSG1//                      &
     &          MSGL//'A BAR Width/Height-pair Does Not Exist'//               &
     &          MSGL//'That Reproduces The Input Cross Section Moments.'       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            SECTION_1D(NS1)%Width  = w
            SECTION_1D(NS1)%Height = h
            SECTION_1D(NS1)%Twall  = DZ
          ENDIF
!!
!! BOX, hollow rectangle/square.
!!
        ELSE IF (Isec .EQ. 4) THEN
!!
!! Iterate for cross section dimensions that reproduce input section properties.
!!
          d = SQRT(Area/0.36D+0)
          Ratio = (By/Bz)**(0.33333333333D+0)
          w = d * Ratio
          h = d / Ratio
          t = d * 0.1D+0
          Icount = 0
          Imaxit = NINT (BCSMXIT)
          C1 = ARHS(w,h,t) - Area
          C2 = BRHy(w,h,t) - By
          C3 = BRHz(w,h,t) - Bz
          Residual = ABS(C1) + ABS(C2) + ABS(C3)

          WRITE (IO_UNIT%LELO,'(/2X,A,I8,3X,A/1X,A,A)')                        &
     &          'Beam Section (BSECTION) ID:',SECTION_1D(NS1)%SecID,           &
     &          '(BOX -- Hollow Rectangle/Square)',                            &
     &          'Iteration Count,===Residual==,====Width====,',                &
     &          '====Height===,====Twall===='

          DO WHILE (Residual .GT. Epsilon*Area .AND. Icount .LE. Imaxit)

            WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                       &
     &        Icount,Residual,w,h,t

            A11 = DADW(t)
            A12 = DADH(t)
            A13 = DADT(w,h,t)
            A21 = DYDW(w,h,t)*0.0833333333D+0
            A22 = DYDH(w,h,t)*0.0833333333D+0
            A23 = DYDT(w,h,t)*0.0833333333D+0
            A31 = DZDW(w,h,t)*0.0833333333D+0
            A32 = DZDH(w,h,t)*0.0833333333D+0
            A33 = DZDT(w,h,t)*0.0833333333D+0
            DET = A11*A22*A33+A12*A23*A31+A13*A21*A32                          &
     &          - A31*A22*A13-A32*A23*A11-A33*A21*A12
            B11 = (A22*A33-A32*A23)*(D1/DET)
            B12 = (A32*A13-A12*A33)*(D1/DET)
            B13 = (A12*A23-A13*A22)*(D1/DET)
            B21 = (A31*A23-A21*A33)*(D1/DET)
            B22 = (A11*A33-A31*A13)*(D1/DET)
            B23 = (A21*A13-A11*A23)*(D1/DET)
            B31 = (A21*A32-A31*A22)*(D1/DET)
            B32 = (A31*A12-A11*A32)*(D1/DET)
            B33 = (A22*A11-A12*A21)*(D1/DET)
            x = w - BCSRLAX * (B11*C1 + B12*C2 + B13*C3)
            y = h - BCSRLAX * (B21*C1 + B22*C2 + B23*C3)
            z = t - BCSRLAX * (B31*C1 + B32*C2 + B33*C3)
            w = MAX (GRLMmin*w, MIN (GRLMmax*w, x))
            h = MAX (GRLMmin*h, MIN (GRLMmax*h, y))
            t = MAX (GRLMmin*t, MIN (GRLMmax*t, z))
            Icount = Icount + 1
            C1 = ARHS(w,h,t) - Area
            C2 = BRHy(w,h,t) - By
            C3 = BRHz(w,h,t) - Bz
            Residual = ABS(C1) + ABS(C2) + ABS(C3)
          ENDDO

          WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                         &
     &      Icount,Residual,w,h,t

          IF (ABS(Area-ARHS(w,h,t)) .GT. Epsilon*Area) THEN
            WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
            CALL USER_MESSAGE                                                  &
     &          (                                                              &
     &          MSGL//'WARN'//                                                 &
     &          MSGL//'SECTION_PROPERTIES.002.04'//                            &
     &          MSGL//'BSECTION Input Record ID:'//MSG1//                      &
     &          MSGL//'A BOX Width/Height/Twall-triple Does Not Exist'//       &
     &          MSGL//'That Reproduces The Input Cross Section Moments.'       &
     &          )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            SECTION_1D(NS1)%Width  = w
            SECTION_1D(NS1)%Height = h
            SECTION_1D(NS1)%Twall  = t
          ENDIF
!!
!! I-beam cross section
!!
        ELSE IF (Isec .EQ. 5) THEN
!!
!! Iterate for cross section dimensions that reproduce input section properties.
!!
          w = 3.464101615D+0 * SQRT (By/Area)
          h = SQRT (Bz/Area)
          t = Area / (w + w + h)
          Icount = 0
          Imaxit = NINT (BCSMXIT)
          C1 = AIBS(w,h,t,t) - Area
          C2 = BIBy(w,h,t,t) - By
          C3 = BIBz(w,h,t,t) - Bz
          Residual = ABS(C1) + ABS(C2) + ABS(C3)

          WRITE (IO_UNIT%LELO,'(/2X,A,I8,3X,A/1X,A,A)')                        &
     &          'Beam Section (BSECTION) ID:',SECTION_1D(NS1)%SecID,           &
     &          '(I-Beam -- T-Flange .EQS. T-Web)',                            &
     &          'Iteration Count,===Residual==,====Width====,',                &
     &          '====Height===,====T-Web===='

          DO WHILE (Residual .GT. Epsilon*Area .AND. Icount .LE. Imaxit)

            WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                       &
     &        Icount,Residual,w,h,t

            A11 = t+t
            A12 = t
            A13 = w+w+h-t-t-t-t
            A21 = 0.0833333333D+0*(D6*w*w*t)
            A22 = 0.0833333333D+0*(t*t*t)
            A23 = 0.0833333333D+0*(D2*(w**3)+(h+h+h-D8)*(t**3))
            A31 = 0.0833333333D+0*((h**3)-((h-t-t)**3))
            A32 = 0.2500000000D+0*(w*h*h-(w-t)*(h-t-t)*(h-t-t))
            A33 = 0.0833333333D+0*(h-t-t)*(h-t-t)*(h+D6*w-D8*t)
            DET = A11*A22*A33+A12*A23*A31+A13*A21*A32                          &
     &          - A31*A22*A13-A32*A23*A11-A33*A21*A12
            B11 = (A22*A33-A32*A23)*(D1/DET)
            B12 = (A32*A13-A12*A33)*(D1/DET)
            B13 = (A12*A23-A13*A22)*(D1/DET)
            B21 = (A31*A23-A21*A33)*(D1/DET)
            B22 = (A11*A33-A31*A13)*(D1/DET)
            B23 = (A21*A13-A11*A23)*(D1/DET)
            B31 = (A21*A32-A31*A22)*(D1/DET)
            B32 = (A31*A12-A11*A32)*(D1/DET)
            B33 = (A22*A11-A12*A21)*(D1/DET)
            x = w - BCSRLAX * (B11*C1 + B12*C2 + B13*C3)
            y = h - BCSRLAX * (B21*C1 + B22*C2 + B23*C3)
            z = t - BCSRLAX * (B31*C1 + B32*C2 + B33*C3)
            w = MAX (GRLMmin*w, MIN (GRLMmax*w, x))
            h = MAX (GRLMmin*h, MIN (GRLMmax*h, y))
            t = MAX (GRLMmin*t, MIN (GRLMmax*t, z))
            Icount = Icount + 1
            C1 = AIBS(w,h,t,t) - Area
            C2 = BIBy(w,h,t,t) - By
            C3 = BIBz(w,h,t,t) - Bz
            Residual = ABS(C1) + ABS(C2) + ABS(C3)
          ENDDO

          WRITE (IO_UNIT%LELO,'(4X,I8,4X,4(1PE14.4))')                         &
     &      Icount,Residual,w,h,t

          IF (ABS(Area-AIBS(w,h,t,t)) .GT. Epsilon*Area) THEN
            WRITE (MSG1,'(I8)') SECTION_1D(NS1)%SecID
            CALL USER_MESSAGE                                                  &
     &      (                                                                  &
     &      MSGL//'WARN'//                                                     &
     &      MSGL//'SECTION_PROPERTIES.002.05'//                                &
     &      MSGL//'BSECTION Input Record ID:'//MSG1//                          &
     &      MSGL//'An I-Beam Width/Height/Twall/'                              &
     &          //'Tflange-tuple Does Not Exist'//                             &
     &      MSGL//'That Reproduces The Input Cross Section Moments.'           &
     &      )
            ERROR%COUNT = ERROR%COUNT + 1
          ELSE
            SECTION_1D(NS1)%Width   = w
            SECTION_1D(NS1)%Height  = h
            SECTION_1D(NS1)%Twall   = t
            SECTION_1D(NS1)%Tflange = t
          ENDIF
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE INITIALIZE_MOTION
!!
!! Copyright (c) by KEY Associates, 26-FEB-1991 20:35:32
!!
!! Purpose: For those nodal points for which a non-zero velocity initial
!! condition is specified move data from the temporary array VELOCITY_IC
!! to the permanant array MOTION%
!!
      USE shared_common_data
      USE node_
      USE motion_
      USE velocity_ic_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered INITIALIZE_MOTION.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Initialize translational velocity.
!!
      DO n = 1,NUMIC
        NP = VELOCITY_IC(n)%NPID
        MOTION(NP)%Vx = VELOCITY_IC(n)%Vx
        MOTION(NP)%Vy = VELOCITY_IC(n)%Vy
        MOTION(NP)%Vz = VELOCITY_IC(n)%Vz
      ENDDO
!!
!! Initialize rotational velocity.
!!
      IF (NUMRT .GT. NUMNP) THEN
        DO n = 1,NUMIC
          NP = NODE(VELOCITY_IC(n)%NPID)%IRT
          IF (NP .GT. 0) THEN
            MOTION(NP)%Vx = VELOCITY_IC(n)%Ox
            MOTION(NP)%Vy = VELOCITY_IC(n)%Oy
            MOTION(NP)%Vz = VELOCITY_IC(n)%Oz
          ENDIF
        ENDDO
      ENDIF
!!
      RETURN
      END
!!_
      SUBROUTINE GET_AUXILIARY_STORAGE_LENGTH
!!
!! Copyright (c) by KEY Associates, 10-NOV-1990 10:03:44
!!
!! Purpose: Compute additional storage needed to conduct calculation:
!!      1. Stress components for through-the-thickness integration of shell
!!         elements to obtain membrane and bending stress resultants.
!!      2. State variable storage locations for stress-strain models.
!!      3. Sliding interface working arrays.
!!      4. Nonreflecting boundary condition working array.
!!
      USE shared_common_data
      USE material_
      USE layering_
      USE section_2d_
      USE node_
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
      USE nonreflecting_bc_
      USE sliding_interface_
      USE tabulated_function_
      USE segment_
      USE node_set_
      USE element_set_
      USE segment_set_
      USE enumerated_sets_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                                                                  &
     &          Ipts,                                                          &
     &          LupID,                                                         &
     &          SecID,                                                         &
     &          SetID
      LOGICAL                                                                  &
     &          NEXT_NP_ID,                                                    &
     &          NEXT_EL_ID,                                                    &
     &          NEXT_SEG_ID
!!
      WRITE (IO_UNIT%LELO,*) ' '
      WRITE (IO_UNIT%LELO,*) ' Entered GET_AUXILIARY_STORAGE_LENGTH.'
      WRITE (IO_UNIT%LELO,*) ' '
!!
!! Compute length of auxiliary storage. Auxiliary storage is used for the
!! multiple stress states in plates and the state variables needed by the
!! respective constitutive models.
!!
      NUMST = 0
      NUMAX = 0
!!
      DO N = 1,NUMHX
        NUMAX = NUMAX + MAUX(MATERIAL(HEXAH(N)%PAR%MatID)%Type)
      ENDDO
      DO N = 1,NUMPX
        NUMAX = NUMAX + MAUX(MATERIAL(PENTA(N)%PAR%MatID)%Type)
      ENDDO
      DO N = 1,NUMTX
        NUMAX = NUMAX + MAUX(MATERIAL(TETRA(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMLS
        LupID = LSOLD(N)%PAR%LupID
        DO i = 1,LAYERING(LupID)%Number_of_Layers
          NUMAX = NUMAX + MAUX(MATERIAL(LAYERING(LupID)%MatID(i))%Type)
        ENDDO
      ENDDO
!!
      DO N = 1,NUMM4
        NUMAX = NUMAX + MAUX(MATERIAL(MEMBQ(N)%PAR%MatID)%Type)
      ENDDO
      DO N = 1,NUMM3
        NUMAX = NUMAX + MAUX(MATERIAL(MEMBT(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMTR
        NUMAX = NUMAX + MAUX(MATERIAL(TRUSS(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMP3
        NPL = N
        Ipts = Ipts_PLATT(NPL)
        NUMST = NUMST + Ipts
        NUMAX = NUMAX + Ipts*MAUX(MATERIAL(PLATT(N)%PAR%MatID)%Type)
      ENDDO
      DO N = 1,NUMP4
        NPL = N
        Ipts = Ipts_PLATQ(NPL)
        NUMST = NUMST + Ipts
        NUMAX = NUMAX + Ipts*MAUX(MATERIAL(PLATQ(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMBM
        NUMAX = NUMAX + 16*MAUX(MATERIAL(BEAM(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMSP
        NUMAX = NUMAX + MAUX(MATERIAL(SPRING(N)%PAR%MatID)%Type)
      ENDDO
      DO N = 1,NUMDM
        NUMAX = NUMAX + MAUX(MATERIAL(DAMPER(N)%PAR%MatID)%Type)
      ENDDO
!!
      DO N = 1,NUMSC
        NUMAX = NUMAX + MAUX(MATERIAL(SPRING_BC(N)%MatID)%Type)
      ENDDO
      DO N = 1,NUMVC
        NUMAX = NUMAX + MAUX(MATERIAL(DAMPER_BC(N)%MatID)%Type)
      ENDDO
!!
!! COMPUTE LENGTH OF DATA ARRAYS NEEDED TO IMPLEMENT THE SLIDING INTERFACE
!! CALCULATIONS.
!!
!! For single-surface interfaces, force interface to be slave-master.
!!
      DO Nsi = 1,NUMSI
        IF (SLIDING_INTERFACE(Nsi)%Type .EQ. 1) THEN
          SLIDING_INTERFACE(Nsi)%Isym = 2
        ENDIF
      ENDDO
!!
!! The sliding interface calculations need 10 arrays the length of the longest
!! set of sliding nodal points found from all of the sliding interfaces present.
!!
      NUMCX = 0
!!
!! Count the number of data records needed for the list of sliding nodal points
!! SLIDING_NODE%
!!
      NUMSN = 0

      DO Nsi = 1,NUMSI
        Isym = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Clear NODE(*)%IRB to use as a working array to count the number of
!! sliding nodal points in each side of a sliding interface.
!!
        DO N = 1,NUMNP
          NODE(N)%IRB = 0
        ENDDO
!!
!! Collect sliding nodal points from Side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
!!
!! Excluded case; node sets not allowed with symmetric and master-slave.
!!
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Count the number of entries in NODE(*)%IRB for Side 1 of sliding interface.
!!
        KOUNT = 0
        DO N = 1,NUMNP
          KOUNT = KOUNT + NODE(N)%IRB
        ENDDO
        NUMSN = NUMSN + KOUNT
        NUMCX = MAX (NUMCX, KOUNT)
!!
!! Clear NODE(*)%IRB to use as a working array to count the number of
!! sliding nodal points in each side of a sliding interface.
!!
        DO N = 1,NUMNP
          NODE(N)%IRB = 0
        ENDDO
!!
!! Collect sliding nodal points from Side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
!!
!! Excluded case; node sets not allowed with symmetric and slave-master.
!!
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Count the number of entries in NODE(*)%IRB for Side 2 of sliding interface.
!!
        KOUNT = 0
        DO N = 1,NUMNP
          KOUNT = KOUNT + NODE(N)%IRB
        ENDDO
        NUMSN = NUMSN + KOUNT
        NUMCX = MAX (NUMCX, KOUNT)

      ENDDO
!!
!! Count the number of data records needed for the list of contact elements
!! CONTACT_SURFACE%
!!
!!      Itype = 0/1   = segment-set/nodal-set
!!      Isym  = 0/1/2 = symmetric/master-slave/slave-master
!!
      NUMCE = 0
      DO Nsi = 1,NUMSI
        Isym = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Count contact elements from side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                NUMCE = NUMCE + 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Collect contact elements from side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                NUMCE = NUMCE + 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! Count the number of data records needed for list of nodal points used in
!! the construction of the sliding interface CONTACT_NODE% (Note: NODE(*)%IRB
!! is used as a scratch array to identify nodes used in one or more sliding
!! interfaces.)
!!
!!      Itype = 0/1   = segment-set/nodal-set
!!      Isym  = 0/1/2 = symmetric/master-slave/slave-master
!!
!! Clear NODE(1:NUMNP)%IRB.
!!
      DO i = 1,NUMNP
        NODE(i)%IRB = 0
      ENDDO
!!
!! Mark each node occurring as a sliding node.
!!
      DO Nsi = 1,NUMSI
        Isym  = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Mark sliding nodal points from Side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
!!
!! Excluded case; node sets not allowed with symmetric and master-slave.
!!
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Mark sliding nodal points from Side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
!!
!! Excluded case; node sets not allowed with symmetric and slave-master.
!!
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ELSE IF (Itype .EQ. 1) THEN
              N = 0
              DO WHILE (NEXT_NP_ID(SetID,N))
                NODE(N)%IRB = 1
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! Mark nodes used in contact elements.
!!
      DO Nsi = 1,NUMSI
        Isym  = SLIDING_INTERFACE(Nsi)%Isym
!!
!! Mark nodes in contact elements from Side 1 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ1
        SetID = SLIDING_INTERFACE(Nsi)%S1ID
        IF (Isym .LE. 1 .AND. Itype .EQ. 1) THEN
        ELSE
          IF (Isym .LT. 2) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!!
!! Mark nodes in contact elements from Side 2 of sliding interface.
!!
        Itype = SLIDING_INTERFACE(Nsi)%Typ2
        SetID = SLIDING_INTERFACE(Nsi)%S2ID
        IF (Isym .NE. 1 .AND. Itype .EQ. 1) THEN
        ELSE
          IF (Isym .NE. 1) THEN
            IF (Itype .EQ. 0) THEN
              N = 0
              DO WHILE (NEXT_SEG_ID(SetID,N))
                DO k = 1,3+MIN(1,SEGMENT(N)%PAR%IX(4))
                  NODE(SEGMENT(N)%PAR%IX(k))%IRB = 1
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!!
!! Count the number of marked nodes to determine the number of CONTACT_NODE
!! date records needed.
!!
      NUMCN = 0
      DO i = 1,NUMNP
        NUMCN = NUMCN + NODE(i)%IRB
      ENDDO
!!
!! Count number of segments used to define nonreflecting BC's.
!!
      NXND = 0
      DO NNR = 1,NUMNR
        SetID = NONREFLECTING_BC(NNR)%SetID
        N = 0
        DO WHILE (NEXT_SEG_ID(SetID,N))
          NXND = NXND + 1
        ENDDO
      ENDDO
      NUMND = NXND
!!
      RETURN
      END
