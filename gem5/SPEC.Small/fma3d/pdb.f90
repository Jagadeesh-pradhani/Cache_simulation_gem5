      SUBROUTINE WRITE_TO_PLOTTING_DATABASE
!!
!! Copyright (c) by KEY Associates, 16-MAR-1991 14:33:15
!!
!! Purpose: Write results to binary file suitable for G'N'P postprocessing.
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
      INTEGER                                                                  &
     &          IXSW(4),                                                       &
     &          MatID_DATA              ! External function to return MatID
      REAL(KIND(0D0))                 &
     &          STRESS_DATA,          & ! External function to return stress
     &          BULK_STRAIN,          & ! External function to rtn bulk strain
     &          INTERNAL_ENERGY,      & ! External function to rtn int energy
     &          PRESSURE,             & ! External function to rtn pressure
     &          EFFECTIVE_STRESS,     & ! External function to rtn eff. stress
     &          INTEGRATION_BIN,      & ! External function to rtn integr'n bin
     &          MATERIAL_DIRECTION      ! External function to rtn fiber angles
!!
!! Variable names of data written to database.
!!
      CHARACTER                                                                &
     &          NAMEXY( 3)*10,                                                 &
     &          NAMENP(14)*10,                                                 &
     &          NAMEEL(13)*10,                                                 &
     &          NAMEGL( 7)*10,                                                 &
     &          NAMEMP(11)*10,                                                 &
     &          GAUGE_LABEL*10,                                                &
     &          MASSPROP_LABEL*10
      DATA                                                                     &
     &          NAMEXY /'X','Y','Z'/,                                          &
     &          NAMENP /'Displ_X','Displ_Y','Displ_Z',                         &
     &                  'Vel_X'  ,'Vel_Y'  ,'Vel_Z'  ,                         &
     &                  'Accel_X','Accel_Y','Accel_Z',                         &
     &                  'Force_X','Force_Y','Force_Z',                         &
     &                  'Group','Contact' /,                                   &
     &          NAMEEL /'Sig_X'     ,'Sig_Y'  ,'Sig_Z'   ,                     &
     &                  'Sig_XY'    ,'Sig_XZ' ,'Sig_YZ'  ,                     &
     &                  'Vol_Strain','Int_Eng','Pressure',                     &
     &                  'Sig_Eff'   ,'Bin'    ,'Angle_1' ,                     &
     &                  'Angle_2'   /,                                         &
     &          NAMEGL /'KE','IE','EE','KE+IE','Bal','Dt','Step'/,             &
     &          NAMEMP /'Mass','Xmom','Ymom','Zmom','PxCM','PyCM',             &
     &                  'PzCM','VxCM','VyCM','VzCM','KEng'/
!!
      LOGICAL, SAVE :: EXSTAT = .FALSE.
      LOGICAL, SAVE :: FIRST = .TRUE.
      LOGICAL       :: MATERIAL_TYPE_22
      INTEGER, SAVE :: NPASS = 1
      REAL(KIND(0D0)), PARAMETER :: RTD = 57.29577951D0
!!
!! Names of arrays written to the Plotting_Database.
!!
!!        COORD       NODAL       ELEMENT       GLOBAL
!!
!!          X         Displ_X      Sig_X          KE
!!          Y         Displ_Y      Sig_Y          IE
!!          Z         Displ_Z      Sig_Z          EE
!!                      Vel_X      Sig_XY         KE+IE
!!                      Vel_Y      Sig_XZ         Bal
!!                      Vel_Z      Sig_YZ         Dt
!!                    Accel_X      Vol_Strain     Step
!!                    Accel_Y      Int_E          SG#_1234/5
!!                    Accel_Z      Pressure       ...
!!                    Force_X      Sig_Eff        MP___6xxxx
!!                    Force_Y      Bin            ...
!!                    Force_Z      Angle_1
!!                      Group      Angle_2
!!                    Contact
!!
      IF (FIRST) THEN
        MATERIAL_TYPE_22 = .FALSE.
        DO N = 1,NUMMT
          MATERIAL_TYPE_22 =                                                   &
     &      MATERIAL_TYPE_22 .OR. (MATERIAL(N)%Type .EQ. 22)
        ENDDO
        NCOORD = 3
        NUMXL  = NUMHX + NUMPX + NUMTX      & ! Solid elements
     &         + NUMLS                      & ! Layered solid elements
     &         + 2 * (NUMM4 + NUMM3)        & ! Membrane elements
     &         + 0 * NUMTR                  & ! Truss elements (axial force o
     &         + 2 * (NUMP4 + NUMP3)        & ! Plate elements
     &         + 0 * NUMBM                  & ! Beam elements
     &         + NUMSP + NUMDM              & ! Spring and Damper elements
     &         + NUMSG                        ! Segments
        NVARNP = 14
        NVAREL = 11
        IF (MATERIAL_TYPE_22) NVAREL = 13
        NGLOBL = 7
        NVARGA = 2*NUMG1 + 3*NUMG2 + 6*NUMG3
        NVARMP = 11*NUMMP
        NVARGL = NGLOBL + NVARGA + NVARMP
        NUMIX  = 9
        NUMBT  = NUMBM + NUMTR  !  Inserted into the old PDB definition.
        IBLKNP = 0
        IBLKEL = 0
        IPACK  = 1
        FIRST = .FALSE.
      ENDIF
!!
!! Open and close Plotting_Database on each entry to allow simultaneous
!! monitoring and plotting as the calculations proceed. If this is the start
!! of an "unchanged" restart (NPASS equals 1 and CONTROL%RDSTAR equals 1),
!! then extend Plotting_Database if it already exists. Otherwise, create a
!! new Plotting_Database.
!!
      IF (NPASS .EQ. 1 .AND. CONTROL%RDSTAR .NE. 1) THEN
        EXSTAT = .FALSE.
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =     IO_UNIT%LPDB,                                     &
     &          FILE   =    'fmapdb',                                          &
     &          STATUS =    'NEW',                                             &
     &          FORM   =    'UNFORMATTED'                                      &
     &          )
!!
      ELSE
        INQUIRE (FILE = 'fmapdb', EXIST = EXSTAT)
        OPEN                                                                   &
     &          (                                                              &
     &          UNIT   =     IO_UNIT%LPDB,                                     &
     &          FILE   =    'fmapdb',                                          &
     &          STATUS =    'UNKNOWN',                                         &
     &          FORM   =    'UNFORMATTED',                                     &
     &          POSITION =    'APPEND'                                         &
     &          )
      ENDIF
!!
!! I. HEADER DATA.
!! Initialize Plotting_Database with size parameters, variable names,
!! initial nodal coordinates, element definitions, and element material
!! numbers.
!!
      IF (NPASS .EQ. 1 .AND. .NOT.EXSTAT) THEN
        NPASS = 0
        REWIND (IO_UNIT%LPDB)
!!
        WRITE  (IO_UNIT%LPDB)                  &
     &          JOB_ID_RECORD%CURRENT%TITLE,   & ! Job title
     &          JOB_ID_RECORD%PROGRAM,         & ! Program name
     &          JOB_ID_RECORD%VERSION,         & ! Program version
     &          JOB_ID_RECORD%CURRENT%DATEE,   & ! Date on which job started
     &          JOB_ID_RECORD%CURRENT%TIMEE      ! Time at which job started
!!
!! Write file size parameters.
!!
        WRITE (IO_UNIT%LPDB)  &
     &          NCOORD,       & ! Number of coordinates
     &          NUMNP,        & ! Number of nodal points
     &          NUMXL,        & ! Number of "elements"
     &          NUMIX,        & ! Number of nodes/element
     &          NUMMT,        & ! Number of materials
     &          NVARNP,       & ! Number of variables/node
     &          NVAREL,       & ! Number of variables/element
     &          NVARGL,       & ! Number of global variables
     &          NUMBT,        & ! Number of beams and trusses
     &          IBLKNP,       & ! Number of nodes/block
     &          IBLKEL,       & ! Number of elements/block
     &          IPACK           ! Number of reals/word
!!
!! Write variable names.
!!
        WRITE (IO_UNIT%LPDB) (NAMEXY(i),i = 1,NCOORD)
        WRITE (IO_UNIT%LPDB) (NAMENP(i),i = 1,NVARNP)
        WRITE (IO_UNIT%LPDB) (NAMEEL(i),i = 1,NVAREL)
        WRITE (IO_UNIT%LPDB) (NAMEGL(i),i = 1,NGLOBL),                         &
     &    (GAUGE_LABEL(i), i = 1,NVARGA),                                      &
     &    (MASSPROP_LABEL(NAMEMP,i), i = 1,NVARMP)
!!
!! Write original nodal point positions.
!!
        WRITE (IO_UNIT%LPDB) (MOTION(i)%Px,   i = 1,NUMNP),                    &
     &                       (MOTION(i)%Py,   i = 1,NUMNP),                    &
     &                       (MOTION(i)%Pz,   i = 1,NUMNP),                    &
     &                       (NODE(i)%ID,     i = 1,NUMNP)
!!
!! Write beam/truss records.
!!
        WRITE (IO_UNIT%LPDB)                                                   &
     &          (                                                              &
     &          BEAM(n)%PAR%EleID,                                             &
     &          BEAM(n)%PAR%MatID,                                             &
     &          BEAM(n)%PAR%IX(1),                                             &
     &          BEAM(n)%PAR%IX(2),                                             &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Section,                         &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%NPLoc,                           &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Width,                           &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Height,                          &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Twall,                           &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Tflange,                         &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Tcover,                          &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Yrefloc,                         &
     &          SECTION_1D(BEAM(n)%PAR%SecID)%Zrefloc,                         &
     &          BEAM(n)%RES%Zaxis(1),                                          &
     &          BEAM(n)%RES%Zaxis(2),                                          &
     &          BEAM(n)%RES%Zaxis(3),                                          &
     &          n = 1,NUMBM                                                    &
     &          ),                                                             &
     &          (                                                              &
     &          TRUSS(m)%PAR%EleID,                                            &
     &          TRUSS(m)%PAR%MatID,                                            &
     &          TRUSS(m)%PAR%IX(1),                                            &
     &          TRUSS(m)%PAR%IX(2),                                            &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Section,                        &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%NPLoc,                          &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Width,                          &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Height,                         &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Twall,                          &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Tflange,                        &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Tcover,                         &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Yrefloc,                        &
     &          SECTION_1D(TRUSS(m)%PAR%SecID)%Zrefloc,                        &
     &          0.0,                                                           &
     &          0.0,                                                           &
     &          0.0,                                                           &
     &          m = 1,NUMTR                                                    &
     &          )
!!
!! Write element connectivity.
!! Element-type flag values:
!!      =  -1, 8-node hexahedron    (HXEL/HEXAH)
!!      =  -2, 8-node pentahedron   (PXEL/PENTA)
!!      =  -3, 8-node tetrahedron   (TXEL/TETRA)
!!      =  -4, 6-node pentahedron   (PXEL/PENTA)
!!      =  -5, 4-node tetrahedron   (TXEL/TETRA)
!!      =  -6, 8-node layered solid (LSEL/LSOLD)
!!      =  -7, 3-node membrane      (M3EL/MEMBT)
!!      =  -8, 4-node membrane      (M4EL/MEMBQ)
!!      =  -9, 3-node plate         (P3EL/PLATT)
!!      = -10, 4-node plate         (P4EL/PLATQ)
!!      = -11, 2-node spring        (SPRING/SPRING)
!!      = -12, 2-node damper        (DAMPER/DAMPER)
!!      = -13, 3-node segment       (SEGMENT/SEGMENT)
!!      = -14, 4-node segment       (SEGMENT/SEGMENT)
!!
        NTUPLE = NUMIX - 1

        DO N = 1,NUMHX
          WRITE (IO_UNIT%LPDB) -1,(HEXAH(N)%PAR%IX(i),i = 1,NTUPLE)
        ENDDO
        DO N = 1,NUMPX
          WRITE (IO_UNIT%LPDB)                                                 &
     &      -4,(PENTA(N)%PAR%IX(i),i = 1,6),(0,k=1,NTUPLE-6)
        ENDDO
        DO N = 1,NUMTX
          WRITE (IO_UNIT%LPDB)                                                 &
     &      -5,(TETRA(N)%PAR%IX(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMLS
          WRITE (IO_UNIT%LPDB) -6,(LSOLD(N)%PAR%IX(i),i = 1,NTUPLE)
        ENDDO
        DO N = 1,NUMM3
          IXSW(1) = MEMBT(N)%PAR%IX(1)
          IXSW(2) = MEMBT(N)%PAR%IX(2)
          IXSW(3) = MEMBT(N)%PAR%IX(3)
          IXSW(4) = MEMBT(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -7,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
          IXSW(1) = MEMBT(N)%PAR%IX(2)
          IXSW(2) = MEMBT(N)%PAR%IX(1)
          IXSW(3) = MEMBT(N)%PAR%IX(3)
          IXSW(4) = MEMBT(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -7,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMM4
          IXSW(1) = MEMBQ(N)%PAR%IX(1)
          IXSW(2) = MEMBQ(N)%PAR%IX(2)
          IXSW(3) = MEMBQ(N)%PAR%IX(3)
          IXSW(4) = MEMBQ(N)%PAR%IX(4)
          WRITE (IO_UNIT%LPDB) -8,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
          IXSW(1) = MEMBQ(N)%PAR%IX(2)
          IXSW(2) = MEMBQ(N)%PAR%IX(1)
          IXSW(3) = MEMBQ(N)%PAR%IX(4)
          IXSW(4) = MEMBQ(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -8,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
!!!       DO N = 1,NUMTR
!!!         IXSW(1) = TRUSS(N)%PAR%IX(1)
!!!         IXSW(2) = TRUSS(N)%PAR%IX(2)
!!!         IXSW(3) = TRUSS(N)%PAR%IX(1)
!!!         IXSW(4) = TRUSS(N)%PAR%IX(2)
!!!         WRITE (IO_UNIT%LPDB) (IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
!!!       ENDDO
        DO N = 1,NUMP3
          IXSW(1) = PLATT(N)%PAR%IX(1)
          IXSW(2) = PLATT(N)%PAR%IX(2)
          IXSW(3) = PLATT(N)%PAR%IX(3)
          IXSW(4) = PLATT(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -9,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
          IXSW(1) = PLATT(N)%PAR%IX(2)
          IXSW(2) = PLATT(N)%PAR%IX(1)
          IXSW(3) = PLATT(N)%PAR%IX(3)
          IXSW(4) = PLATT(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -9,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMP4
          IXSW(1) = PLATQ(N)%PAR%IX(1)
          IXSW(2) = PLATQ(N)%PAR%IX(2)
          IXSW(3) = PLATQ(N)%PAR%IX(3)
          IXSW(4) = PLATQ(N)%PAR%IX(4)
          WRITE (IO_UNIT%LPDB) -10,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
          IXSW(1) = PLATQ(N)%PAR%IX(2)
          IXSW(2) = PLATQ(N)%PAR%IX(1)
          IXSW(3) = PLATQ(N)%PAR%IX(4)
          IXSW(4) = PLATQ(N)%PAR%IX(3)
          WRITE (IO_UNIT%LPDB) -10,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMSP
          IXSW(1) = SPRING(N)%PAR%IX(1)
          IXSW(2) = SPRING(N)%PAR%IX(2)
          IXSW(3) = SPRING(N)%PAR%IX(1)
          IXSW(4) = SPRING(N)%PAR%IX(2)
          WRITE (IO_UNIT%LPDB) -11,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMDM
          IXSW(1) = DAMPER(N)%PAR%IX(1)
          IXSW(2) = DAMPER(N)%PAR%IX(2)
          IXSW(3) = DAMPER(N)%PAR%IX(1)
          IXSW(4) = DAMPER(N)%PAR%IX(2)
          WRITE (IO_UNIT%LPDB) -12,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
        ENDDO
        DO N = 1,NUMSG
          IXSW(1) = SEGMENT(N)%PAR%IX(1)
          IXSW(2) = SEGMENT(N)%PAR%IX(2)
          IXSW(3) = SEGMENT(N)%PAR%IX(3)
          IXSW(4) = SEGMENT(N)%PAR%IX(4)
          IF (IXSW(4) .LE. 0) THEN
            WRITE (IO_UNIT%LPDB) -13,(IXSW(i),i = 1,3),(0,k=1,NTUPLE-3)
          ELSE
            WRITE (IO_UNIT%LPDB) -14,(IXSW(i),i = 1,4),(0,k=1,NTUPLE-4)
          ENDIF
        ENDDO
!!
!! Write material property numbers.
!!
        IF (NUMMT .GT. 1) THEN
          WRITE (IO_UNIT%LPDB)(MatID_DATA(n),n = 1,NUMXL)
        ENDIF
      ENDIF
!!
!! II. CURRENT CONFIGURATION DATA
!! Write current time to Plotting_Database.
!!
      WRITE (IO_UNIT%LPDB) TIMSIM%Total
!!
!! Write nodal point based data to Plotting_Database.
!!
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Ux,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Uy,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Uz,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Vx,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Vy,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Vz,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Ax,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Ay,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (MOTION(i)%Az,       i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (FORCE(i)%Xext,      i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (FORCE(i)%Yext,      i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (FORCE(i)%Zext,      i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (DBLE(NODE(i)%ISI),  i = 1,NUMNP)
      WRITE (IO_UNIT%LPDB) (DBLE(NODE(i)%ICF),  i = 1,NUMNP)
!!
!! Write element based data to Plotting_Database.
!!
      DO i = 1,6
        WRITE (IO_UNIT%LPDB) (STRESS_DATA(n,i), n = 1,NUMXL)
      ENDDO
      WRITE (IO_UNIT%LPDB) (BULK_STRAIN(n),      n = 1,NUMXL)
      WRITE (IO_UNIT%LPDB) (INTERNAL_ENERGY(n),  n = 1,NUMXL)
      WRITE (IO_UNIT%LPDB) (PRESSURE(n),         n = 1,NUMXL)
      WRITE (IO_UNIT%LPDB) (EFFECTIVE_STRESS(n), n = 1,NUMXL)
      WRITE (IO_UNIT%LPDB) (INTEGRATION_BIN(n),  n = 1,NUMXL)
      IF (MATERIAL_TYPE_22) THEN
        WRITE (IO_UNIT%LPDB) (MATERIAL_DIRECTION(n,1), n = 1,NUMXL)
        WRITE (IO_UNIT%LPDB) (MATERIAL_DIRECTION(n,2), n = 1,NUMXL)
      ENDIF
!!
!! Write global data to Plotting_Database.
!!
      ENGK = 0.0
      DO N = 1,NUMRT
        ENGK = ENGK + NODE(N)%Mass *                                           &
     &          (                                                              &
     &          MOTION(N)%Vx * MOTION(N)%Vx +                                  &
     &          MOTION(N)%Vy * MOTION(N)%Vy +                                  &
     &          MOTION(N)%Vz * MOTION(N)%Vz                                    &
     &          )
      ENDDO
      ENERGY%Kinetic = 0.5 * ENGK
!!
      IF (ENERGY%External .EQ. 0.0) THEN
        ENERGY%Balance = 0.0
      ELSE
        ENERGY%Balance = 100.0 * (ENERGY%External-ENERGY%Kinetic-              &
     &                              ENERGY%Internal)/ENERGY%External
      ENDIF
!!
      WRITE (IO_UNIT%LPDB)                                                     &
     &        ENERGY%Kinetic,                                                  &
     &        ENERGY%Internal,                                                 &
     &        ENERGY%External,                                                 &
     &        ENERGY%Kinetic+ENERGY%Internal,                                  &
     &        ENERGY%Balance,                                                  &
     &        TIMSIM%DTlast,                                                   &
     &        DBLE(TIMSIM%Step),                                               &
     &        (GAUGE1D(i)%RES%Strain,GAUGE1D(i)%RES%Torsion,i=1,NUMG1),        &
     &        ((GAUGE2D(i)%RES%Strain(k),k=1,3),i=1,NUMG2),                    &
     &        ((GAUGE3D(i)%RES%Strain(k),k=1,6),i=1,NUMG3),                    &
     &        (                                                                &
     &        MASSPROP(i)%Mass,                                                &
     &        MASSPROP(i)%Xmv,                                                 &
     &        MASSPROP(i)%Ymv,                                                 &
     &        MASSPROP(i)%Zmv,                                                 &
     &        MASSPROP(i)%Xcm,                                                 &
     &        MASSPROP(i)%Ycm,                                                 &
     &        MASSPROP(i)%Zcm,                                                 &
     &        MASSPROP(i)%Vxcm,                                                &
     &        MASSPROP(i)%Vycm,                                                &
     &        MASSPROP(i)%Vzcm,                                                &
     &        MASSPROP(i)%KE,                                                  &
     &        i=1,NUMMP                                                        &
     &        )
!!
      CLOSE (UNIT=IO_UNIT%LPDB, STATUS='KEEP')
!!
      RETURN
      END
!!_
      CHARACTER*10 FUNCTION GAUGE_LABEL (Istrain)
!!
!! Copyright (c) by KEY Associates, 5-DEC-1991 22:07:37
!!
!! Purpose: Construct strain gauge labels for Plotting_Database. The label
!! will have the form "SG#_2107/4" where "2107" is the gauge ID and "4" is
!! the strain component.
!!
      USE shared_common_data
      USE gauge1d_
      USE gauge2d_
      USE gauge3d_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          BLANK_LABEL*10
      DATA                                                                     &
     &          BLANK_LABEL  /'SG#_____/_'/
      INTEGER                                                                  &
     &          Istrain
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NG1,NG2,NG3
!!
      IF (FIRST) THEN
        NG1 = 2*NUMG1
        NG2 = NG1 + 3*NUMG2
        NG3 = NG2 + 6*NUMG3
        FIRST = .FALSE.
      ENDIF

      GAUGE_LABEL = BLANK_LABEL
      IF (Istrain .LE. NG1) THEN
        N = (Istrain + 1) / 2
        M = Istrain - 2*(N-1)
        WRITE (GAUGE_LABEL(10:10),'(I1)') M
        WRITE (GAUGE_LABEL(4:8),'(I5)') GAUGE1D(N)%PAR%GauID
      ELSE IF (Istrain .LE. NG2) THEN
        N = (Istrain - NG1 + 2) / 3
        M = (Istrain - NG1) - 3*(N-1)
        WRITE (GAUGE_LABEL(10:10),'(I1)') M
        WRITE (GAUGE_LABEL(4:8),'(I5)') GAUGE2D(N)%PAR%GauID
      ELSE IF (Istrain .LE. NG3) THEN
        N = (Istrain - NG2 + 5) / 6
        M = (Istrain - NG2) - 6*(N-1)
        WRITE (GAUGE_LABEL(10:10),'(I1)') M
        WRITE (GAUGE_LABEL(4:8),'(I5)') GAUGE3D(N)%PAR%GauID
      ENDIF
!!
!! Replace blanks with proper characters.
!!
      DO i = 1,10
        IF (GAUGE_LABEL(i:i) .EQ. ' ') THEN
          GAUGE_LABEL(i:i) = BLANK_LABEL(i:i)
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      CHARACTER*10 FUNCTION MASSPROP_LABEL (NAMEMP,Iname)
!!
!! Copyright (c) by KEY Associates, 8-DEC-1991 17:59:29
!!
!! Purpose: Construct mass property labels for Plotting_Database. The label
!! will have the form "MP_2/MV" where "2" is the mass property ID and "MV"
!! is the mass property result.
!!
      USE shared_common_data
      USE massprop_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      CHARACTER                                                                &
     &          NAMEMP(1)*10,                                                  &
     &          BLANK_LABEL*10
      DATA                                                                     &
     &          BLANK_LABEL /'MP________'/
      INTEGER                                                                  &
     &          Iname
!!
      MASSPROP_LABEL = BLANK_LABEL
      N = (Iname + 10) / 11
      M = Iname - 11*(N-1)
      MASSPROP_LABEL(7:10) = NAMEMP(M)(1:4)
      WRITE (MASSPROP_LABEL(1:6),'(I6)') MASSPROP(N)%MPID
!!
!! Replace blanks with proper characters.
!!
      DO i = 1,10
        IF (MASSPROP_LABEL(i:i) .EQ. ' ') THEN
          MASSPROP_LABEL(i:i) = BLANK_LABEL(i:i)
        ENDIF
      ENDDO
!!
      RETURN
      END
!!_
      INTEGER FUNCTION MatID_DATA ( NEL )
!!
!! Copyright (c) by KEY Associates, 12-APR-1991 18:59:44
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
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL  ! I/- N-th element
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM,NSG
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        NSG = NDM + NUMSG
        FIRST = .FALSE.
      ENDIF
!!
      MatID_DATA = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        MatID_DATA = HEXAH(N)%PAR%MatID
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        MatID_DATA = PENTA(N)%PAR%MatID
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        MatID_DATA = TETRA(N)%PAR%MatID
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        MatID_DATA = LSOLD(N)%PAR%LupID
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        MatID_DATA = MEMBT(N)%PAR%MatID
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        MatID_DATA = MEMBQ(N)%PAR%MatID
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        MatID_DATA = TRUSS(N)%PAR%MatID
      ELSE IF (NEL .LE. NP3) THEN
        N = (NEL - NTR + 1) / 2
        MatID_DATA = PLATT(N)%PAR%MatID
      ELSE IF (NEL .LE. NP4) THEN
        N = (NEL - NP3 + 1) / 2
        MatID_DATA = PLATQ(N)%PAR%MatID
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NP4
        MatID_DATA = SPRING(N)%PAR%MatID
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
        MatID_DATA = DAMPER(N)%PAR%MatID
      ELSE IF (NEL .LE. NSG) THEN
        N = NEL - NDM
        MatID_DATA = SEGMENT(N)%PAR%ParID
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION STRESS_DATA ( NEL,I )
!!
!! Copyright (c) by KEY Associates, 8-APR-1991 20:33:30
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
      USE stress_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
      INTEGER                 &
     &          I,            & ! I/- Stress component requested
     &          NEL,          & ! I/- N-th element
     &          Ist,          & ! -/- Pointer into STRESS(I,Ist)
     &          Ipts            ! -/- Number of integration points in shell

      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      STRESS_DATA = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        STRESS_DATA = HEXAH(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        STRESS_DATA = PENTA(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        STRESS_DATA = TETRA(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        STRESS_DATA = 0.0
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        IF (I .LE. 2) THEN
          STRESS_DATA = MEMBT(N)%RES%Stress(I)
        ELSE IF (I .EQ. 4) THEN
          STRESS_DATA = MEMBT(N)%RES%Stress(3)
        ENDIF
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        IF (I .LE. 2) THEN
          STRESS_DATA = MEMBQ(N)%RES%Stress(I)
        ELSE IF (I .EQ. 4) THEN
          STRESS_DATA = MEMBQ(N)%RES%Stress(3)
        ENDIF
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        IF (I .EQ. 1) THEN
          STRESS_DATA = TRUSS(N)%RES%Stress
        ELSE
          STRESS_DATA = 0.0
        ENDIF
      ELSE IF (NEL .LE. NP3) THEN
        M = NEL - NTR + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATT(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATT(N)
          Ist = PLATT(N)%PAR%Ist + Ipts - 1
        ENDIF
        STRESS_DATA = STRESS(I,Ist)
      ELSE IF (NEL .LE. NP4) THEN
        M = NEL - NP3 + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATQ(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATQ(N)
          Ist = PLATQ(N)%PAR%Ist + Ipts - 1
        ENDIF
        STRESS_DATA = STRESS(I,Ist)
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NP4
        IF (I .EQ. 1) THEN
          STRESS_DATA = SPRING(N)%RES%Force
        ELSE
          STRESS_DATA = 0.0
        ENDIF
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
        IF (I .EQ. 1) THEN
          STRESS_DATA = DAMPER(N)%RES%Force
        ELSE
          STRESS_DATA = 0.0
        ENDIF
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION BULK_STRAIN ( NEL )
!!
!! Copyright (c) by KEY Associates, 12-APR-1991 18:59:44
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
      INTEGER, INTENT(IN) :: NEL  ! I/- N-th element
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
      REAL(KIND(0D0))  X
      LOG_POS(X) = LOG(MAX(1.0D-4,X))
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      BULK_STRAIN = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        BULK_STRAIN = LOG_POS(HEXAH(N)%RES%Volume/HEXAH(N)%PAR%Volume)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        BULK_STRAIN = LOG_POS(PENTA(N)%RES%Volume/PENTA(N)%PAR%Volume)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        BULK_STRAIN = LOG_POS(TETRA(N)%RES%Volume/TETRA(N)%PAR%Volume)
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        BULK_STRAIN = LOG_POS(LSOLD(N)%RES%Volume/LSOLD(N)%PAR%Volume)
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        BULK_STRAIN = LOG_POS(MEMBT(N)%RES%Area/MEMBT(N)%PAR%Area)
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        BULK_STRAIN = LOG_POS(MEMBQ(N)%RES%Area/MEMBQ(N)%PAR%Area)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        BULK_STRAIN = LOG_POS(TRUSS(N)%RES%Length/TRUSS(N)%PAR%Length)
      ELSE IF (NEL .LE. NP3) THEN
        N = (NEL - NTR + 1) / 2
        BULK_STRAIN = LOG_POS(PLATT(N)%RES%Area/PLATT(N)%PAR%Area)
      ELSE IF (NEL .LE. NP4) THEN
        N = (NEL - NP3 + 1) / 2
        BULK_STRAIN = LOG_POS(PLATQ(N)%RES%Area/PLATQ(N)%PAR%Area)
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NP4
        BULK_STRAIN = SPRING(N)%RES%Delta
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
        BULK_STRAIN = DAMPER(N)%RES%Delta
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION INTERNAL_ENERGY ( NEL )
!!
!! Copyright (c) by KEY Associates, 12-APR-1991 18:59:44
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
      INTEGER, INTENT(IN) :: NEL   ! I/- N-th element
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      INTERNAL_ENERGY = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        INTERNAL_ENERGY = HEXAH(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        INTERNAL_ENERGY = PENTA(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        INTERNAL_ENERGY = TETRA(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        INTERNAL_ENERGY = LSOLD(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        INTERNAL_ENERGY = MEMBT(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        INTERNAL_ENERGY = MEMBQ(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        INTERNAL_ENERGY = TRUSS(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NP3) THEN
        N = (NEL - NTR + 1) / 2
        INTERNAL_ENERGY = PLATT(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NP4) THEN
        N = (NEL - NP3 + 1) / 2
        INTERNAL_ENERGY = PLATQ(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NP4
        INTERNAL_ENERGY = SPRING(N)%RES%Int_Eng
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
        INTERNAL_ENERGY = DAMPER(N)%RES%Int_Eng
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION PRESSURE ( NEL )
!!
!! Copyright (c) by KEY Associates, 8-APR-1991 20:33:30
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
      USE stress_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL  ! I/- N-th element
!!
!! Local Variables.
      INTEGER :: Ist              ! -/- Pointer into STRESS(1:6,Ist)
      INTEGER :: Ipts             ! -/- Number of integration points in shell
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      PRESSURE = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          HEXAH(N)%RES%Stress(1) +                                       &
     &          HEXAH(N)%RES%Stress(2) +                                       &
     &          HEXAH(N)%RES%Stress(3)                                         &
     &          )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          PENTA(N)%RES%Stress(1) +                                       &
     &          PENTA(N)%RES%Stress(2) +                                       &
     &          PENTA(N)%RES%Stress(3)                                         &
     &          )
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          TETRA(N)%RES%Stress(1) +                                       &
     &          TETRA(N)%RES%Stress(2) +                                       &
     &          TETRA(N)%RES%Stress(3)                                         &
     &          )
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        PRESSURE = 0.0
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress(1) +                                       &
     &          MEMBT(N)%RES%Stress(2)                                         &
     &          )
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress(1) +                                       &
     &          MEMBQ(N)%RES%Stress(2)                                         &
     &          )
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress                                            &
     &          )
      ELSE IF (NEL .LE. NP3) THEN
        M = NEL - NTR + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATT(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATT(N)
          Ist = PLATT(N)%PAR%Ist + Ipts - 1
        ENDIF
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          STRESS(1,Ist) +                                                &
     &          STRESS(2,Ist)                                                  &
     &          )
      ELSE IF (NEL .LE. NP4) THEN
        M = NEL - NP3 + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATQ(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATQ(N)
          Ist = PLATQ(N)%PAR%Ist + Ipts - 1
        ENDIF
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          STRESS(1,Ist) +                                                &
     &          STRESS(2,Ist)                                                  &
     &          )
      ELSE IF (NEL .LE. NSP) THEN
      ELSE IF (NEL .LE. NDM) THEN
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION EFFECTIVE_STRESS ( NEL )
!!
!! Copyright (c) by KEY Associates, 8-APR-1991 20:33:30
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
      USE stress_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL  ! I/- N-th element
!!
!! Local variables.
      INTEGER :: Ist              ! -/- Pointer into STRESS(1:6,Ist)
      INTEGER :: Ipts             ! -/- Number of integration points in shell
      REAL(KIND(0D0)) :: PRESSURE ! -/- Used to get deviatoric stress components
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
!! Effective stress calculation using diviatoric stress components.
!!
      EFF_STR (S1,S2,S3,S4) = SQRT (1.5*(S1*S1+S2*S2+S3*S3+S4+S4))
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      EFFECTIVE_STRESS = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          HEXAH(N)%RES%Stress(1) +                                       &
     &          HEXAH(N)%RES%Stress(2) +                                       &
     &          HEXAH(N)%RES%Stress(3)                                         &
     &          )
        S1 = HEXAH(N)%RES%Stress(1) - PRESSURE
        S2 = HEXAH(N)%RES%Stress(2) - PRESSURE
        S3 = HEXAH(N)%RES%Stress(3) - PRESSURE
        S4 = HEXAH(N)%RES%Stress(4) * HEXAH(N)%RES%Stress(4) +                 &
     &       HEXAH(N)%RES%Stress(5) * HEXAH(N)%RES%Stress(5) +                 &
     &       HEXAH(N)%RES%Stress(6) * HEXAH(N)%RES%Stress(6)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          PENTA(N)%RES%Stress(1) +                                       &
     &          PENTA(N)%RES%Stress(2) +                                       &
     &          PENTA(N)%RES%Stress(3)                                         &
     &          )
        S1 = PENTA(N)%RES%Stress(1) - PRESSURE
        S2 = PENTA(N)%RES%Stress(2) - PRESSURE
        S3 = PENTA(N)%RES%Stress(3) - PRESSURE
        S4 = PENTA(N)%RES%Stress(4) * PENTA(N)%RES%Stress(4) +                 &
     &       PENTA(N)%RES%Stress(5) * PENTA(N)%RES%Stress(5) +                 &
     &       PENTA(N)%RES%Stress(6) * PENTA(N)%RES%Stress(6)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          TETRA(N)%RES%Stress(1) +                                       &
     &          TETRA(N)%RES%Stress(2) +                                       &
     &          TETRA(N)%RES%Stress(3)                                         &
     &          )
        S1 = TETRA(N)%RES%Stress(1) - PRESSURE
        S2 = TETRA(N)%RES%Stress(2) - PRESSURE
        S3 = TETRA(N)%RES%Stress(3) - PRESSURE
        S4 = TETRA(N)%RES%Stress(4) * TETRA(N)%RES%Stress(4) +                 &
     &       TETRA(N)%RES%Stress(5) * TETRA(N)%RES%Stress(5) +                 &
     &       TETRA(N)%RES%Stress(6) * TETRA(N)%RES%Stress(6)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        EFFECTIVE_STRESS = 0.0
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          MEMBT(N)%RES%Stress(1) +                                       &
     &          MEMBT(N)%RES%Stress(2)                                         &
     &          )
        S1 = MEMBT(N)%RES%Stress(1) - PRESSURE
        S2 = MEMBT(N)%RES%Stress(2) - PRESSURE
        S3 =                        - PRESSURE
        S4 = MEMBT(N)%RES%Stress(3) * MEMBT(N)%RES%Stress(3)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          MEMBQ(N)%RES%Stress(1) +                                       &
     &          MEMBQ(N)%RES%Stress(2)                                         &
     &          )
        S1 = MEMBQ(N)%RES%Stress(1) - PRESSURE
        S2 = MEMBQ(N)%RES%Stress(2) - PRESSURE
        S3 =                        - PRESSURE
        S4 = MEMBQ(N)%RES%Stress(3) * MEMBQ(N)%RES%Stress(3)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          TRUSS(N)%RES%Stress                                            &
     &          )
        S1 = TRUSS(N)%RES%Stress - PRESSURE
        S2 =                     - PRESSURE
        S3 =                     - PRESSURE
        S4 = 0.0
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NP3) THEN
        M = NEL - NTR + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATT(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATT(N)
          Ist = PLATT(N)%PAR%Ist + Ipts - 1
        ENDIF
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          STRESS(1,Ist) +                                                &
     &          STRESS(2,Ist)                                                  &
     &          )
        S1 = STRESS(1,Ist) - PRESSURE
        S2 = STRESS(2,Ist) - PRESSURE
        S3 =               - PRESSURE
        S4 = STRESS(4,Ist) * STRESS(4,Ist)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NP4) THEN
        M = NEL - NP3 + 1
        N = M / 2
        IF (MOD(M,2) .EQ. 0) THEN
          Ist = PLATQ(N)%PAR%Ist
        ELSE
          Ipts = Ipts_PLATQ(N)
          Ist = PLATQ(N)%PAR%Ist + Ipts - 1
        ENDIF
        PRESSURE = 0.3333333333 *                                              &
     &          (                                                              &
     &          STRESS(1,Ist) +                                                &
     &          STRESS(2,Ist)                                                  &
     &          )
        S1 = STRESS(1,Ist) - PRESSURE
        S2 = STRESS(2,Ist) - PRESSURE
        S3 =               - PRESSURE
        S4 = STRESS(4,Ist) * STRESS(4,Ist)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NSP) THEN
      ELSE IF (NEL .LE. NDM) THEN
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION INTEGRATION_BIN ( NEL )
!!
!! Copyright (c) by KEY Associates, 11-NOV-1991 13:20:34
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
      INTEGER, INTENT(IN) :: NEL  ! I/- N-th element
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        FIRST = .FALSE.
      ENDIF
!!
      INTEGRATION_BIN = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        INTEGRATION_BIN = HEXAH(N)%RES%ISI
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        INTEGRATION_BIN = PENTA(N)%RES%ISI
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
        INTEGRATION_BIN = TETRA(N)%RES%ISI
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
        INTEGRATION_BIN = LSOLD(N)%RES%ISI
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        INTEGRATION_BIN = MEMBT(N)%RES%ISI
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        INTEGRATION_BIN = MEMBQ(N)%RES%ISI
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
        INTEGRATION_BIN = TRUSS(N)%RES%ISI
      ELSE IF (NEL .LE. NP3) THEN
        N = (NEL - NTR + 1) / 2
        INTEGRATION_BIN = PLATT(N)%RES%ISI
      ELSE IF (NEL .LE. NP4) THEN
        N = (NEL - NP3 + 1) / 2
        INTEGRATION_BIN = PLATQ(N)%RES%ISI
      ELSE IF (NEL .LE. NBM) THEN
        N = NEL - NP4
        INTEGRATION_BIN = BEAM(N)%RES%ISI
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NBM
        INTEGRATION_BIN = SPRING(N)%RES%ISI
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
        INTEGRATION_BIN = DAMPER(N)%RES%ISI
      ENDIF
!!
      RETURN
      END
!!_
      REAL(KIND(0D0))  FUNCTION MATERIAL_DIRECTION ( NEL,Iangle )
!!
!! Copyright (c) by KEY Associates, 2-DEC-1991 20:26:34
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
      USE material_
      USE state_variables_
!!
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)
!!
!! Arguments.
      INTEGER, INTENT(IN) :: NEL     ! I/- N-th element
      INTEGER, INTENT(IN) :: IAngle  ! I/- I-th material direction
!!
!! Local variables.
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: NHX,NPX,NTX,NLS,NM3,NM4,NTR,NP3,NP4,NBM,NSP,NDM
      REAL(KIND(0D0)), SAVE :: Sign_A
      REAL(KIND(0D0)), PARAMETER :: RTD = 57.2957795131D0
      REAL(KIND(0D0)) :: Ar,As,Br,Bs
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NTX = NPX + NUMTX
        NLS = NTX + NUMLS
        NM3 = NLS + 2*NUMM3
        NM4 = NM3 + 2*NUMM4
        NTR = NM4 + 0 ! NUMTR
        NP3 = NTR + 2*NUMP3
        NP4 = NP3 + 2*NUMP4
        NBM = NP4 + 0 ! NUMBM
        NSP = NBM + NUMSP
        NDM = NSP + NUMDM
        Sign_A = -ONE
        FIRST = .FALSE.
      ENDIF
!!
      MATERIAL_DIRECTION = 720.0D0
      IF (NEL .LE. NHX) THEN
        N = NEL
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPX
      ELSE IF (NEL .LE. NLS) THEN
        N = NEL - NTX
      ELSE IF (NEL .LE. NM3) THEN
        N = (NEL - NLS + 1) / 2
        Sign_A = -Sign_A
        IF (MATERIAL(MEMBT(N)%PAR%MatID)%Type .EQ. 22) THEN
          Isv = MEMBT(N)%PAR%Isv
          IF (IAngle .EQ. 1) THEN
            Ar = STATE_VARIABLES(Isv+0)
            As = STATE_VARIABLES(Isv+1)
            MATERIAL_DIRECTION = Sign_A * SIGN (RTD*ACOS(Ar), As)
          ELSE IF (IAngle .EQ. 2) THEN
            Br = STATE_VARIABLES(Isv+2)
            Bs = STATE_VARIABLES(Isv+3)
            MATERIAL_DIRECTION = Sign_A * SIGN (RTD*ACOS(Br), Bs)
          ENDIF
        ENDIF
      ELSE IF (NEL .LE. NM4) THEN
        N = (NEL - NM3 + 1) / 2
        Sign_A = -Sign_A
        IF (MATERIAL(MEMBQ(N)%PAR%MatID)%Type .EQ. 22) THEN
          Isv = MEMBQ(N)%PAR%Isv
          IF (IAngle .EQ. 1) THEN
            Ar = STATE_VARIABLES(Isv+0)
            As = STATE_VARIABLES(Isv+1)
            MATERIAL_DIRECTION = Sign_A * SIGN (RTD*ACOS(Ar), As)
          ELSE IF (IAngle .EQ. 2) THEN
            Br = STATE_VARIABLES(Isv+2)
            Bs = STATE_VARIABLES(Isv+3)
            MATERIAL_DIRECTION = Sign_A * SIGN (RTD*ACOS(Br), Bs)
          ENDIF
        ENDIF
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NM4
      ELSE IF (NEL .LE. NP3) THEN
        N = (NEL - NTR + 1) / 2
      ELSE IF (NEL .LE. NP4) THEN
        N = (NEL - NP3 + 1) / 2
      ELSE IF (NEL .LE. NSP) THEN
        N = NEL - NP4
      ELSE IF (NEL .LE. NDM) THEN
        N = NEL - NSP
      ENDIF
!!
      RETURN
      END
