      SUBROUTINE SUBCYCLING_PARTITION (SUBCYCLING)
!!
!! Copyright (c) by S W Key, 13-NOV-1991 08:48:00
!!
!! Purpose: Set up the node and element groupings for subcycling.
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
      USE segment_
      USE node_
      USE rigid_body_
      USE contact_node_
      USE spring_bc_
      USE damper_bc_
      USE periodic_bc_
      USE displacement_bc_
      USE tied_bc_
      USE rigid_wall_bc_
      USE nonreflecting_bc_
      USE spot_weld_
      USE constrained_node_
      USE node_set_
      USE enumerated_sets_, ONLY: NNPSETS
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      LOGICAL, INTENT(IN) :: SUBCYCLING
!!
!! Local variables.
      CHARACTER(32) :: CGROUP(8)*32 ! Partition table for printing

      INTEGER       :: ISImin,ISImax
      INTEGER       :: IX1          ! Local scratch variable
      INTEGER       :: IX2          ! Local scratch variable
      INTEGER, SAVE :: NCR          ! Ratio of one subcycling group to the next
      INTEGER, SAVE :: KCS          ! Number of times to loop on defining *.ISI
      INTEGER, SAVE :: KPF          ! Frequncy of mandatory partitioning
      INTEGER       :: SetID
      INTEGER, SAVE :: IPCONT = 0   ! Partition counter
      INTEGER, SAVE :: IPPRNT = 1   ! Partition info print control counter
      INTEGER       :: Number_per_Group(8) ! Number of elements/group

      REAL(KIND(0D0)), SAVE :: LGNI

      LOGICAL       :: QUIESCENT
      LOGICAL       :: NEXT_NP_ID
      LOGICAL       :: NEXT_SEG_ID
      LOGICAL       :: REPARTITION
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! Obtain partition ratio constants NCR and LGNI, interface spreading loop-
!! limit KCS, and frequency for mandatory repartitioning KPF.
!!
      IF (FIRST) THEN
        IF (SUBCYCLING) THEN
          NCR = NINT (MIN (  4.0D+0, MAX (2.0D+0,PARAMVALUE%Cycle_R)))
          LGNI = ONE / LOG (DBLE (NCR))
          KCS = NINT (MIN (  5.0D+0, MAX (ONE,PARAMVALUE%Cycle_S)))
          KPF = NINT (MIN (100.0D+0, MAX (ONE,PARAMVALUE%Cycle_F)))
        ELSE
          NCR = 1
        ENDIF
        FIRST = .FALSE.
      ENDIF
!!
!! If the cycle ratio NCR has been set to 1, no subcycling will occur.
!! The following block of code is a condensed version of this module for
!! the case of no subcycling. That is, the case where every nodal point
!! and element is integrated with the same time step. At the bottom of
!! the if-test a return is executed.
!!
      IF (NCR .EQ. 1) THEN
        TIMSIM%Cymax = 1
!$OMP PARALLEL DO
        DO N = 1,NUMRT
          NODE(N)%ISI = 1
          NODE(N)%DTlast = TIMSIM%DTlast
          NODE(N)%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMHX
          HEXAH(N)%RES%ISI = 1
          HEXAH(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMPX
          PENTA(N)%RES%ISI = 1
          PENTA(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMTX
          TETRA(N)%RES%ISI = 1
          TETRA(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMLS
          LSOLD(N)%RES%ISI = 1
          LSOLD(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMM3
          MEMBT(N)%RES%ISI = 1
          MEMBT(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMM4
          MEMBQ(N)%RES%ISI = 1
          MEMBQ(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMTR
          TRUSS(N)%RES%ISI = 1
          TRUSS(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMP3
          PLATT(N)%RES%ISI = 1
          PLATT(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
!$OMP PARALLEL DO
        DO N = 1,NUMP4
          PLATQ(N)%RES%ISI = 1
          PLATQ(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMBM
          BEAM(N)%RES%ISI = 1
          BEAM(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMSP
          SPRING(N)%RES%ISI = 1
          SPRING(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        DO N = 1,NUMDM
          DAMPER(N)%RES%ISI = 1
          DAMPER(N)%RES%DTnext = TIMSIM%DTnext
        ENDDO
        RETURN
      ENDIF
!!
!! Before embarking on the straight forward but arduous task of creating
!! a subcycling partition of the nodes and elements, check to see if the
!! existing partition will suffice. Note: The frequecy parameter KPF may
!! force a repartition.
!!

      ITSTEP = TIMSIM%Step
      REPARTITION = (MOD(ITSTEP,KPF).EQ.0) .OR. REZONE%REPARTITION
!!
      IF (.NOT.REPARTITION) THEN
        DO N = 1,NUMHX
          REPARTITION = REPARTITION .OR.                                       &
     &    (HEXAH(N)%RES%DTelt .LT. DBLE(HEXAH(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
        DO N = 1,NUMPX
          REPARTITION = REPARTITION .OR.                                       &
     &    (PENTA(N)%RES%DTelt .LT. DBLE(PENTA(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
        DO N = 1,NUMTX
          REPARTITION = REPARTITION .OR.                                       &
     &    (TETRA(N)%RES%DTelt .LT. DBLE(TETRA(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMLS
          REPARTITION = REPARTITION .OR.                                       &
     &    (LSOLD(N)%RES%DTelt .LT. DBLE(LSOLD(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMM3
          REPARTITION = REPARTITION .OR.                                       &
     &    (MEMBT(N)%RES%DTelt .LT. DBLE(MEMBT(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
        DO N = 1,NUMM4
          REPARTITION = REPARTITION .OR.                                       &
     &    (MEMBQ(N)%RES%DTelt .LT. DBLE(MEMBQ(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMTR
          REPARTITION = REPARTITION .OR.                                       &
     &    (TRUSS(N)%RES%DTelt .LT. DBLE(TRUSS(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMP3
          REPARTITION = REPARTITION .OR.                                       &
     &    (PLATT(N)%RES%DTelt .LT. DBLE(PLATT(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
        DO N = 1,NUMP4
          REPARTITION = REPARTITION .OR.                                       &
     &    (PLATQ(N)%RES%DTelt .LT. DBLE(PLATQ(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMBM
          REPARTITION = REPARTITION .OR.                                       &
     &    (BEAM(N)%RES%DTelt .LT. DBLE(BEAM(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
!!
        DO N = 1,NUMSP
          REPARTITION = REPARTITION .OR.                                       &
     &    (SPRING(N)%RES%DTelt .LT. DBLE(SPRING(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
        DO N = 1,NUMDM
        REPARTITION = REPARTITION .OR.                                         &
     &    (DAMPER(N)%RES%DTelt .LT. DBLE(DAMPER(N)%RES%ISI)*TIMSIM%DTnext)
        ENDDO
      ENDIF
!!
!! Check to see if a repartition is required.
!!
      IF (REPARTITION) THEN
!!
!! In the event that a rezoning caused the repartion, set the flag to false
!! and wait for the next rezoning to set the flag back to true.
!!
        REZONE%REPARTITION = .FALSE.
!!
!! Scan each element type and find integer multiple of TIMSIM%DTmin (The
!! effect of using TIMSIM%DTmin in place of TIMSIM%DTnext is to scale all
!! time steps down. Thus, the more conservative TIMSIM%DTnext is, the more
!! conservative is the entire subcycling process.) The following module
!! is a table look-up for ISI = NCR**INT(LOG(DTelt/TIMSIM%DTmin)*LGNI)
!!
        CALL ELEMENT_PARTITION ( NCR,ISImax,Number_per_Group )
!!
        TIMSIM%Cymax = ISImax
!!
!! Loop over node/element *.ISI definitions. Succeeding loops spread or
!! diffuse the interface further and further into slower elements, that
!! is, elements with larger critical time steps.
!!
        DO K = 1,KCS
!!
!! Assign each nodal point the minimum subcycling index ISI of the elements
!! attached to it.
!!
          DO N = 1,NUMRT
            NODE(N)%ISI = ISImax
          ENDDO
!!
          DO N = 1,NUMHX
            DO i = 1,8
              NODE(HEXAH(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(HEXAH(N)%PAR%IX(i))%ISI, HEXAH(N)%RES%ISI)
            ENDDO
          ENDDO
          DO N = 1,NUMPX
            DO i = 1,6
              NODE(PENTA(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(PENTA(N)%PAR%IX(i))%ISI, PENTA(N)%RES%ISI)
            ENDDO
          ENDDO
          DO N = 1,NUMTX
            DO i = 1,4
              NODE(TETRA(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(TETRA(N)%PAR%IX(i))%ISI, TETRA(N)%RES%ISI)
            ENDDO
          ENDDO
!!
          DO N = 1,NUMLS
            DO i = 1,8
              NODE(LSOLD(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(LSOLD(N)%PAR%IX(i))%ISI, LSOLD(N)%RES%ISI)
            ENDDO
          ENDDO
!!
          DO N = 1,NUMM3
            DO i = 1,3
              NODE(MEMBT(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(MEMBT(N)%PAR%IX(i))%ISI, MEMBT(N)%RES%ISI)
            ENDDO
          ENDDO
          DO N = 1,NUMM4
            DO i = 1,4
              NODE(MEMBQ(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(MEMBQ(N)%PAR%IX(i))%ISI, MEMBQ(N)%RES%ISI)
            ENDDO
          ENDDO
!!
          DO N = 1,NUMTR
            DO i = 1,2
              NODE(TRUSS(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(TRUSS(N)%PAR%IX(i))%ISI, TRUSS(N)%RES%ISI)
            ENDDO
          ENDDO
!!
          DO N = 1,NUMP3
            DO i = 1,3
              NODE(PLATT(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(PLATT(N)%PAR%IX(i))%ISI, PLATT(N)%RES%ISI)
            ENDDO
          ENDDO
          DO N = 1,NUMP4
            DO i = 1,4
              NODE(PLATQ(N)%PAR%IX(i))%ISI =                                   &
     &          MIN (NODE(PLATQ(N)%PAR%IX(i))%ISI, PLATQ(N)%RES%ISI)
            ENDDO
          ENDDO
!!
          DO N = 1,NUMBM
            DO i = 1,2
              NODE(BEAM(N)%PAR%IX(i))%ISI =                                    &
     &          MIN (NODE(BEAM(N)%PAR%IX(i))%ISI, BEAM(N)%RES%ISI)
            ENDDO
          ENDDO
!!
!! Note: The rotational springs and dampers (Type = 1) reference the
!! rotational degress of freedom.
!!
          DO N = 1,NUMSP
            IF (SPRING(N)%PAR%Type .EQ. 0) THEN
              IX1 = SPRING(N)%PAR%IX(1)
              IX2 = SPRING(N)%PAR%IX(2)
            ELSE
              IX1 = NODE(SPRING(N)%PAR%IX(1))%IRT
              IX2 = NODE(SPRING(N)%PAR%IX(2))%IRT
            ENDIF
            NODE(IX1)%ISI = MIN (NODE(IX1)%ISI, SPRING(N)%RES%ISI)
            NODE(IX2)%ISI = MIN (NODE(IX2)%ISI, SPRING(N)%RES%ISI)
          ENDDO
          DO N = 1,NUMDM
            IF (DAMPER(N)%PAR%Type .EQ. 0) THEN
              IX1 = DAMPER(N)%PAR%IX(1)
              IX2 = DAMPER(N)%PAR%IX(2)
            ELSE
              IX1 = NODE(DAMPER(N)%PAR%IX(1))%IRT
              IX2 = NODE(DAMPER(N)%PAR%IX(2))%IRT
            ENDIF
            NODE(IX1)%ISI = MIN (NODE(IX1)%ISI, DAMPER(N)%RES%ISI)
            NODE(IX2)%ISI = MIN (NODE(IX2)%ISI, DAMPER(N)%RES%ISI)
          ENDDO
!!
!! Put all nodes that are contained in a rigid body in the lowest partition.
!!
          DO N = 1,NUMRB
            M = RIGID_BODY(N)%FirstNP
            DO WHILE (M .GT. 0)
              NODE(M)%ISI = 1
              M = NODE(M)%IRB
            ENDDO
          ENDDO
!!
!! Put all nodes that are on a contact surface in the lowest partition.
!!
          DO N = 1,NUMCN
            NODE(CONTACT_NODE(N)%NPID)%ISI = 1
          ENDDO
!!
!! Put all nodes that have displacement BC's in the lowest partition.
!! (The imposition of displacement BC's (IMPOSE_DISPLACEMENT_BC) should
!! now be able to handle subcycling.)
!!
!!!         DO NDC = 1,NUMDC
!!!           SetID = DISPLACEMENT_BC(NDC)%SetID
!!!           N = 0
!!!           DO WHILE (NEXT_NP_ID(SetID,N))
!!!             NODE(N)%ISI = 1
!!!           ENDDO
!!!         ENDDO
!!
!! Put all nodes that have tied BC's in the lowest partition.
!!
          DO NTC = 1,NUMTC
            SetID = TIED_BC(NTC)%SetID
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              NODE(N)%ISI = 1
            ENDDO
          ENDDO
!!
!! Put all nodes that have spring or damper BC's in the lowest partition.
!!
          DO NSC = 1,NUMSC
            SetID = SPRING_BC(NSC)%SetID
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              NODE(N)%ISI = 1
            ENDDO
          ENDDO
          DO NVC = 1,NUMVC
            SetID = DAMPER_BC(NVC)%SetID
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              NODE(N)%ISI = 1
            ENDDO
          ENDDO
!!
!! Put all nodes that participate in a rigid wall BC in the lowest partition.
!!
          DO NWC = 1,NUMWC
            SetID = RIGID_WALL_BC(NWC)%SetID
            N = 0
            DO WHILE (NEXT_NP_ID(SetID,N))
              NODE(N)%ISI = 1
            ENDDO
          ENDDO
!!
!! Put all nodes that participate in a nonreflecting boundary condition
!! in the lowest partition.
!!
          DO NNR = 1,NUMNR
            SetID = NONREFLECTING_BC(NNR)%SetID
            N = 0
            DO WHILE (NEXT_SEG_ID(SetID,N))
              DO i = 1,4
                IF (SEGMENT(N)%PAR%IX(i) .GT. 0) THEN
                  NODE(SEGMENT(N)%PAR%IX(i))%ISI = 1
                ENDIF
              ENDDO
            ENDDO
          ENDDO
!!
!! Insure that nodal pairs occurring in periodic boundary conditions
!! both have the same (and smallest between them) partition.
!!
          DO NCC = 1,NUMCC
            SetID = PERIODIC_BC(NCC)%S1ID
            N1 = NODE_SET(SetID)%Istart - 1
            SetID = PERIODIC_BC(NCC)%S2ID
            N2 = NODE_SET(SetID)%Istart - 1
            DO i = NODE_SET(SetID)%Istart,NODE_SET(SetID)%Iend
              N1 = N1 + 1
              N2 = N2 + 1
              NP1 = NNPSETS(N1)
              NP2 = NNPSETS(N2)
              ISImin = MIN (NODE(NP1)%ISI,NODE(NP2)%ISI)
              NODE(NP1)%ISI = ISImin
              NODE(NP2)%ISI = ISImin
            ENDDO
          ENDDO
!!
!! Insure that nodal pairs (or triples) occurring in spot weld
!! connections all have the same (and smallest between them) partition.
!!
          DO NSW = 1,NUMSW
            IF (SPOT_WELD(NSW)%PLACE .EQ. 'POSITION') THEN
              ISImin = ISImax
              DO i = 1,3
                IF (SPOT_WELD(NSW)%EleID(i) .GT. 0) THEN
                  IF (SPOT_WELD(NSW)%Type(i) .EQ. 0) THEN
                    DO n = 1,3
                      ISImin = MIN (ISImin,                                    &
     &              NODE(PLATT(SPOT_WELD(NSW)%EleID(i))%PAR%IX(n))%ISI)
                    ENDDO
                  ELSE IF (SPOT_WELD(NSW)%Type(i) .EQ. 1) THEN
                    DO n = 1,4
                      ISImin = MIN (ISImin,                                    &
     &              NODE(PLATQ(SPOT_WELD(NSW)%EleID(i))%PAR%IX(n))%ISI)
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
              DO i = 1,3
                IF (SPOT_WELD(NSW)%EleID(i) .GT. 0) THEN
                  IF (SPOT_WELD(NSW)%Type(i) .EQ. 0) THEN
                    DO n = 1,3
                      NODE(PLATT(SPOT_WELD(NSW)                                &
     &                  %EleID(i))%PAR%IX(n))%ISI = ISImin
                    ENDDO
                  ELSE IF (SPOT_WELD(NSW)%Type(i) .EQ. 1) THEN
                    DO n = 1,4
                      NODE(PLATQ(SPOT_WELD(NSW)                                &
     &                  %EleID(i))%PAR%IX(n))%ISI = ISImin
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
            ELSE IF (SPOT_WELD(NSW)%PLACE .EQ. 'NODES') THEN
              ISImin = ISImax
              DO i = 1,3
                IF (SPOT_WELD(NSW)%NPID(i) .GT. 0) THEN
                  ISImin = MIN (ISImin,NODE(SPOT_WELD(NSW)%NPID(i))%ISI)
                ENDIF
              ENDDO
              DO i = 1,3
                IF (SPOT_WELD(NSW)%NPID(i) .GT. 0) THEN
                  NODE(SPOT_WELD(NSW)%NPID(i))%ISI = ISImin
                ENDIF
              ENDDO
            ENDIF
          ENDDO
!!
!! Insure that nodal-triples occurring in constrained node boundary
!! conditions all have the same (and smallest between them) partition.
!! Because it is possible for two constraind nodes to have a common
!! end point, the following do-loop must be repeated until there are
!! no more changes.
!!
          QUIESCENT = .FALSE.
          DO WHILE (.NOT.QUIESCENT)
            QUIESCENT = .TRUE.
            DO NNC = 1,NUMNC
              IMS = CONSTRAINED_NODE(NNC)%CNID
              ID1 = CONSTRAINED_NODE(NNC)%NPID(1)
              ID2 = CONSTRAINED_NODE(NNC)%NPID(2)
              ISImin = MIN (NODE(IMS)%ISI,NODE(ID1)%ISI,NODE(ID2)%ISI)
              ISIsup = MAX (NODE(IMS)%ISI,NODE(ID1)%ISI,NODE(ID2)%ISI)
              NODE(IMS)%ISI = ISImin
              NODE(ID1)%ISI = ISImin
              NODE(ID2)%ISI = ISImin
              QUIESCENT = QUIESCENT .AND. (ISImin .EQ. ISIsup)
            ENDDO
          ENDDO
!!
!! Reset each element subcycling index to the minimum found among the nodal
!! points defining the element.
!!
          DO N = 1,NUMHX
            ISImin = ISImax
            DO i = 1,8
              ISImin = MIN (ISImin,NODE(HEXAH(N)%PAR%IX(i))%ISI)
            ENDDO
            HEXAH(N)%RES%ISI = ISImin
          ENDDO
          DO N = 1,NUMPX
            ISImin = ISImax
            DO i = 1,6
              ISImin = MIN (ISImin,NODE(PENTA(N)%PAR%IX(i))%ISI)
            ENDDO
            PENTA(N)%RES%ISI = ISImin
          ENDDO
          DO N = 1,NUMTX
            ISImin = ISImax
            DO i = 1,4
              ISImin = MIN (ISImin,NODE(TETRA(N)%PAR%IX(i))%ISI)
            ENDDO
            TETRA(N)%RES%ISI = ISImin
          ENDDO
!!
          DO N = 1,NUMLS
            ISImin = ISImax
            DO i = 1,8
              ISImin = MIN (ISImin,NODE(LSOLD(N)%PAR%IX(i))%ISI)
            ENDDO
            LSOLD(N)%RES%ISI = ISImin
          ENDDO
!!
          DO N = 1,NUMM3
            ISImin = ISImax
            DO i = 1,3
              ISImin = MIN (ISImin,NODE(MEMBT(N)%PAR%IX(i))%ISI)
            ENDDO
            MEMBT(N)%RES%ISI = ISImin
          ENDDO
          DO N = 1,NUMM4
            ISImin = ISImax
            DO i = 1,4
              ISImin = MIN (ISImin,NODE(MEMBQ(N)%PAR%IX(i))%ISI)
            ENDDO
            MEMBQ(N)%RES%ISI = ISImin
          ENDDO
!!
          DO N = 1,NUMTR
            TRUSS(N)%RES%ISI =                                                 &
     &        MIN (NODE(TRUSS(N)%PAR%IX(1))%ISI,                               &
     &             NODE(TRUSS(N)%PAR%IX(2))%ISI)
          ENDDO
!!
          DO N = 1,NUMP3
            ISImin = ISImax
            DO i = 1,3
              ISImin = MIN (ISImin,NODE(PLATT(N)%PAR%IX(i))%ISI)
            ENDDO
            PLATT(N)%RES%ISI = ISImin
          ENDDO
          DO N = 1,NUMP4
            ISImin = ISImax
            DO i = 1,4
              ISImin = MIN (ISImin,NODE(PLATQ(N)%PAR%IX(i))%ISI)
            ENDDO
            PLATQ(N)%RES%ISI = ISImin
          ENDDO
!!
          DO N = 1,NUMBM
            BEAM(N)%RES%ISI =                                                  &
     &        MIN (NODE(BEAM(N)%PAR%IX(1))%ISI,                                &
     &             NODE(BEAM(N)%PAR%IX(2))%ISI)
          ENDDO
!!
          DO N = 1,NUMSP
            IF (SPRING(N)%PAR%Type .EQ. 0) THEN
              IX1 = SPRING(N)%PAR%IX(1)
              IX2 = SPRING(N)%PAR%IX(2)
            ELSE
              IX1 = NODE(SPRING(N)%PAR%IX(1))%IRT
              IX2 = NODE(SPRING(N)%PAR%IX(2))%IRT
            ENDIF
            SPRING(N)%RES%ISI = MIN (NODE(IX1)%ISI,NODE(IX2)%ISI)
          ENDDO
          DO N = 1,NUMDM
            IF (DAMPER(N)%PAR%Type .EQ. 0) THEN
              IX1 = DAMPER(N)%PAR%IX(1)
              IX2 = DAMPER(N)%PAR%IX(2)
            ELSE
              IX1 = NODE(DAMPER(N)%PAR%IX(1))%IRT
              IX2 = NODE(DAMPER(N)%PAR%IX(2))%IRT
            ENDIF
            DAMPER(N)%RES%ISI = MIN (NODE(IX1)%ISI,NODE(IX2)%ISI)
          ENDDO
!!
!! End of KCS subcycling interface spreading do-loop.
!!
        ENDDO
!!
!! Update integration subcyling indicies (*.ISI) for rotational degrees of
!! freedom to match integration subcyling indicies for translational degrees
!! of freedom. Since only the finite elements with rotational degress of
!! freedom "know" where the translational degrees of freedom are stored,
!! the element pointers are used. (The storage of rotational degrees of freedom
!! are easily accessed in Version 7 of the code, but the storage of rotational
!! degrees of freedom in Version 8 with adaptivity is not easily accessed except
!! through the elements with rotational degrees of freedom due to the dynamic
!! addition and deletion of elements.)
!!
        DO N = 1,NUMP3
          DO i = 1,3
            NODE(PLATT(N)%PAR%IX(i+3))%ISI =                                   &
     &        NODE(PLATT(N)%PAR%IX(i))%ISI
          ENDDO
        ENDDO
        DO N = 1,NUMP4
          DO i = 1,4
            NODE(PLATQ(N)%PAR%IX(i+4))%ISI =                                   &
     &        NODE(PLATQ(N)%PAR%IX(i))%ISI
          ENDDO
        ENDDO
        DO N = 1,NUMBM
          NODE(BEAM(N)%PAR%IX(3))%ISI =                                        &
     &      NODE(BEAM(N)%PAR%IX(1))%ISI
          NODE(BEAM(N)%PAR%IX(4))%ISI =                                        &
     &      NODE(BEAM(N)%PAR%IX(2))%ISI
        ENDDO
        DO N = 1,NUMSP
          IF (SPRING(N)%PAR%Type .EQ. 1) THEN
            IX1 = SPRING(N)%PAR%IX(1)
            IX2 = SPRING(N)%PAR%IX(2)
            NODE(IX1)%ISI = NODE(NODE(IX1)%IRT)%ISI
            NODE(IX2)%ISI = NODE(NODE(IX2)%IRT)%ISI
          ENDIF
        ENDDO
        DO N = 1,NUMDM
          IF (DAMPER(N)%PAR%Type .EQ. 1) THEN
            IX1 = DAMPER(N)%PAR%IX(1)
            IX2 = DAMPER(N)%PAR%IX(2)
            NODE(IX1)%ISI = NODE(NODE(IX1)%IRT)%ISI
            NODE(IX2)%ISI = NODE(NODE(IX2)%IRT)%ISI
          ENDIF
        ENDDO
!!
!! Inform user of subcycling partition.
!!
        IPCONT = IPCONT + 1
        IF (IPCONT .GE. IPPRNT) THEN
          IPPRNT = IPPRNT + IPPRNT
          DO i = 1,8
            CGROUP(i) = ' '
            CGROUP(i)(10:10) = '/'
            CGROUP(i)(21:21) = '/'
            WRITE (CGROUP(i)( 4: 4),'(I1)') i
            WRITE (CGROUP(i)(13:18),'(I6)') NCR**(i-1)
            WRITE (CGROUP(i)(24:29),'(I6)') Number_per_Group(i)
          ENDDO

          WRITE (MSG1,'(     I8)') TIMSIM%Step
          WRITE (MSGF,'(1PE12.4)') TIMSIM%Total
          WRITE (MSG2,'(     I8)') IPCONT
          CALL USER_MESSAGE                                                    &
     &          (                                                              &
     &          MSGL//'INFORM'//                                               &
     &          MSGL//'PARTITION.001.00'//                                     &
     &          MSGL//'Time Step Number:   '//MSG1//                           &
     &          MSGL//'Integration Time:   '//MSGF//                           &
     &          MSGL//'Partitioning Number:'//MSG2//                           &
     &          MSGL//' Group # /  Cycles  / Elements'//                       &
     &          MSGL//CGROUP(1)//MSGL//CGROUP(2)//                             &
     &          MSGL//CGROUP(3)//MSGL//CGROUP(4)//                             &
     &          MSGL//CGROUP(5)//MSGL//CGROUP(6)//                             &
     &          MSGL//CGROUP(7)//MSGL//CGROUP(8)                               &
     &          )
        ENDIF
!!
!! End of check for need-to-repartition nodes and elements if-test.
!!
      ENDIF
!!
!! Define DTlast and DTnext for all nodes and elements.
!!
      DO N = 1,NUMRT
        NODE(N)%DTlast = NODE(N)%DTnext
        NODE(N)%DTnext = DBLE (NODE(N)%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMHX
        HEXAH(N)%RES%DTnext = DBLE (HEXAH(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
      DO N = 1,NUMPX
        PENTA(N)%RES%DTnext = DBLE (PENTA(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
      DO N = 1,NUMTX
        TETRA(N)%RES%DTnext = DBLE (TETRA(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMLS
        LSOLD(N)%RES%DTnext = DBLE (LSOLD(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMM3
        MEMBT(N)%RES%DTnext = DBLE (MEMBT(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
      DO N = 1,NUMM4
        MEMBQ(N)%RES%DTnext = DBLE (MEMBQ(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMTR
        TRUSS(N)%RES%DTnext = DBLE (TRUSS(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMP3
        PLATT(N)%RES%DTnext = DBLE (PLATT(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
      DO N = 1,NUMP4
        PLATQ(N)%RES%DTnext = DBLE (PLATQ(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMBM
        BEAM(N)%RES%DTnext = DBLE (BEAM(N)%RES%ISI) * TIMSIM%DTnext
      ENDDO
!!
      DO N = 1,NUMSP
        SPRING(N)%RES%DTnext = DBLE (SPRING(N)%RES%ISI)*TIMSIM%DTnext
      ENDDO
      DO N = 1,NUMDM
        DAMPER(N)%RES%DTnext = DBLE (DAMPER(N)%RES%ISI)*TIMSIM%DTnext
      ENDDO
!!
      RETURN
      END
!!_
      SUBROUTINE ELEMENT_PARTITION ( NCR,ISImax,Number_per_Group )
!!
!! Copyright (c) by S W Key, 23-JAN-1992 13:31:40
!!
!! Purpose: Set up the element groupings for subcycling.
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
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)  :: NCR     ! Ratio between subcycling groups.
      INTEGER, INTENT(OUT) :: ISImax  ! Maximum subcycling index found.
      INTEGER, INTENT(OUT) :: Number_per_Group(8) ! Elements/group
!!
!! Local variables.
      INTEGER, SAVE :: INTmax = 32767 ! Maximum 2-byte integer (2**15 - 1).
      INTEGER       :: IGroup         ! Group number assigned to ISI
!!
!! External functions
      INTEGER       :: IS_INDEX  ! Function that returns appropriate IS index.
      EXTERNAL         IS_INDEX
!!
!! Invert the minimum time step.
!!
      DTinv = ONE / TIMSIM%DTmin
!!
!! Largest permitted element time step DTemx. (This value is needed to skip
!! over "massless" elements with time steps equal to 1.0D+37)
!!
      DTemx = TIMSIM%DTmin * DBLE (INTmax)
!!
!! Divide the element time step by the minimum time step; convert to an integer.
!!
      DO i = 1,NUMHX
        IF (HEXAH(i)%RES%DTelt .LT. DTemx) THEN
          HEXAH(i)%RES%ISI = INT (DTinv * HEXAH(i)%RES%DTelt)
        ELSE
          HEXAH(i)%RES%ISI = 1
        ENDIF
      ENDDO
      DO i = 1,NUMPX
        IF (PENTA(i)%RES%DTelt .LT. DTemx) THEN
          PENTA(i)%RES%ISI = INT (DTinv * PENTA(i)%RES%DTelt)
        ELSE
          PENTA(i)%RES%ISI = 1
        ENDIF
      ENDDO
      DO i = 1,NUMTX
        IF (TETRA(i)%RES%DTelt .LT. DTemx) THEN
          TETRA(i)%RES%ISI = INT (DTinv * TETRA(i)%RES%DTelt)
        ELSE
          TETRA(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMLS
        IF (LSOLD(i)%RES%DTelt .LT. DTemx) THEN
          LSOLD(i)%RES%ISI = INT (DTinv * LSOLD(i)%RES%DTelt)
        ELSE
          LSOLD(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMM3
        IF (MEMBT(i)%RES%DTelt .LT. DTemx) THEN
          MEMBT(i)%RES%ISI = INT (DTinv * MEMBT(i)%RES%DTelt)
        ELSE
          MEMBT(i)%RES%ISI = 1
        ENDIF
      ENDDO
      DO i = 1,NUMM4
        IF (MEMBQ(i)%RES%DTelt .LT. DTemx) THEN
          MEMBQ(i)%RES%ISI = INT (DTinv * MEMBQ(i)%RES%DTelt)
        ELSE
          MEMBQ(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMTR
        IF (TRUSS(i)%RES%DTelt .LT. DTemx) THEN
          TRUSS(i)%RES%ISI = INT (DTinv * TRUSS(i)%RES%DTelt)
        ELSE
          TRUSS(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMP3
        IF (PLATT(i)%RES%DTelt .LT. DTemx) THEN
          PLATT(i)%RES%ISI = INT (DTinv * PLATT(i)%RES%DTelt)
        ELSE
          PLATT(i)%RES%ISI = 1
        ENDIF
      ENDDO
      DO i = 1,NUMP4
        IF (PLATQ(i)%RES%DTelt .LT. DTemx) THEN
          PLATQ(i)%RES%ISI = INT (DTinv * PLATQ(i)%RES%DTelt)
        ELSE
          PLATQ(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMBM
        IF (BEAM(i)%RES%DTelt .LT. DTemx) THEN
          BEAM(i)%RES%ISI = INT (DTinv * BEAM(i)%RES%DTelt)
        ELSE
          BEAM(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
      DO i = 1,NUMSP
        IF (SPRING(i)%RES%DTelt .LT. DTemx) THEN
          SPRING(i)%RES%ISI = INT (DTinv * SPRING(i)%RES%DTelt)
        ELSE
          SPRING(i)%RES%ISI = 1
        ENDIF
      ENDDO
      DO i = 1,NUMDM
        IF (DAMPER(i)%RES%DTelt .LT. DTemx) THEN
          DAMPER(i)%RES%ISI = INT (DTinv * DAMPER(i)%RES%DTelt)
        ELSE
          DAMPER(i)%RES%ISI = 1
        ENDIF
      ENDDO
!!
!! Clear counter in order to record the number of elements in each group.
!!
      DO i = 1,8
        Number_per_Group(i) = 0
      ENDDO
!!
!! Select integration group with highest multiple less than *.ISI%
!!
      ISImax = 0
      DO i = 1,NUMHX
        HEXAH(i)%RES%ISI = IS_INDEX (NCR,HEXAH(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,HEXAH(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
      DO i = 1,NUMPX
        PENTA(i)%RES%ISI = IS_INDEX (NCR,PENTA(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,PENTA(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
      DO i = 1,NUMTX
        TETRA(i)%RES%ISI = IS_INDEX (NCR,TETRA(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,TETRA(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMLS
        LSOLD(i)%RES%ISI = IS_INDEX (NCR,LSOLD(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,LSOLD(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMM3
        MEMBT(i)%RES%ISI = IS_INDEX (NCR,MEMBT(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,MEMBT(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
      DO i = 1,NUMM4
        MEMBQ(i)%RES%ISI = IS_INDEX (NCR,MEMBQ(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,MEMBQ(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMTR
        TRUSS(i)%RES%ISI = IS_INDEX (NCR,TRUSS(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,TRUSS(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMP3
        PLATT(i)%RES%ISI = IS_INDEX (NCR,PLATT(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,PLATT(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
      DO i = 1,NUMP4
        PLATQ(i)%RES%ISI = IS_INDEX (NCR,PLATQ(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,PLATQ(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMBM
        BEAM(i)%RES%ISI = IS_INDEX (NCR,BEAM(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,BEAM(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      DO i = 1,NUMSP
        SPRING(i)%RES%ISI = IS_INDEX (NCR,SPRING(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,SPRING(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
      DO i = 1,NUMDM
        DAMPER(i)%RES%ISI = IS_INDEX (NCR,DAMPER(i)%RES%ISI,IGroup)
        ISImax = MAX (ISImax,DAMPER(i)%RES%ISI)
        Number_per_Group(IGroup) = Number_per_Group(IGroup) + 1
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION IS_INDEX ( NCR,ISI,IGroup )
!!
!! Copyright (c) by S W Key, 23-JAN-1992 12:33:40
!!
!! Purpose: Select integration subcycling index from ISI_GROUPS that is
!! the largest value less than the input value ISI.
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN)  :: NCR     ! Cycle ratio between groups
      INTEGER, INTENT(IN)  :: ISI     ! Element time step / minimum step
      INTEGER, INTENT(OUT) :: IGroup  ! Group number assigned to ISI
!!
!! Local variables.
      INTEGER, SAVE :: ISI_GROUPS(1:8,2:4) ! Available subcycling groups
      LOGICAL, SAVE :: FIRST = .TRUE.
!!
!! The available subcycling groups are given by ISI_GROUPS(m,N) = N**(m-1)
!!
      IF (FIRST) THEN
        ISI_GROUPS(1:8,2) = (/1,2,4,8,16,32,64,128/)
        ISI_GROUPS(1:8,3) = (/1,3,9,27,81,243,729,2187/)
        ISI_GROUPS(1:8,4) = (/1,4,16,64,256,1024,4096,16384/)
        FIRST = .FALSE.
      ENDIF
!!
!! Scan ISI_GROUPS(1:8,NCR) to find the maximum permissable IS index.
!!
      IGroup = 1
      IS_INDEX = 1
      DO I = 2,8
        IF (ISI .LT. ISI_GROUPS(I,NCR)) RETURN
        IGroup = I
        IS_INDEX = ISI_GROUPS(I,NCR)
      ENDDO
!!
      RETURN
      END
