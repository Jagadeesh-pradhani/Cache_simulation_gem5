!!
!!       Common Blocks And Structures Used Throughout FMA-3D.
!!
!! Copyright (c) by KEY Associates;  9-SEP-1997 22:37:39.00
!!
      MODULE shared_common_data
!!
!! Initially, counters for the ASCII input records to be read and stored.
!! When the simulation storage has finally been determined, these counters
!! will indicate the total number of records to read and store. Some of
!! the counters are for temporary storage, for example, NUMIC which is the
!! number of initial condition records to read. Some of the counters are
!! just for the heck of it, for example, NUMPV which is the number of
!! "PARAM" input records, an input record that only requires a fixed amount
!! of storage irrespective of the number of records in the input file.
!!

!!    IMPLICIT LOGICAL (A-Z)
      IMPLICIT REAL(KIND(0D0)) (A-H,O-Z)

      REAL(KIND(0D0)), PARAMETER :: ONE = 1.0D+0

      COMMON /BLOCK_03_NUM/ NUMIF,NUMQA,NUMNP,NUMEL,NUMHX,NUMPX,NUMTX,NU&
     &MLS,NUMLX, NUMLM,NUMM4,NUMM3,NUMTR,NUMP4,NUMP3,NUMBM,NUMSP,NUMDM, &
     &NUMSG,NUMDC,NUMTC,NUMSW,NUMWC,NUMBF,NUMPC,NUMFC,NUMSC, NUMVC,NUMCC&
     &,NUMNR,NUMND,NUMFS,NUMIT,NUMSI,NUMSN,NUMCE, NUMCN,NUMCX,NUMNS,NUMN&
     &E,NUMES,NUMEE,NUMSS,NUMSE,NUMMT, NUMLU,NUMTF,NUMFP,NUMRF,NUMPV,NUM&
     &S1,NUMS2,NUMG1,NUMG2, NUMG3,NUMMP,NUMRM,NUMCM,NUMIC,NUMRB,NUMRT,NU&
     &MST,NUMAX, NUMPP,NUMNC,MAXRE,MXEID,MAUX(0:60),MSET(0:60)

      INTEGER               NUMIF,NUMQA,NUMNP,NUMEL,NUMHX,NUMPX,NUMTX,NU&
     &MLS,NUMLX, NUMLM,NUMM4,NUMM3,NUMTR,NUMP4,NUMP3,NUMBM,NUMSP,NUMDM, &
     &NUMSG,NUMDC,NUMTC,NUMSW,NUMWC,NUMBF,NUMPC,NUMFC,NUMSC, NUMVC,NUMCC&
     &,NUMNR,NUMND,NUMFS,NUMIT,NUMSI,NUMSN,NUMCE, NUMCN,NUMCX,NUMNS,NUMN&
     &E,NUMES,NUMEE,NUMSS,NUMSE,NUMMT, NUMLU,NUMTF,NUMFP,NUMRF,NUMPV,NUM&
     &S1,NUMS2,NUMG1,NUMG2, NUMG3,NUMMP,NUMRM,NUMCM,NUMIC,NUMRB,NUMRT,NU&
     &MST,NUMAX, NUMPP,NUMNC,MAXRE,MXEID,MAUX,MSET
!!
!! Counters for the records contained in the mesh data input. For ease of
!! management, the mesh counters exactly parallel the NUMxx counters. Not
!! all MSHxx counters have corresponding mesh data input.
!!
      COMMON /BLOCK_03_MSH/ MSHIF,MSHQA,MSHNP,MSHEL,MSHHX,MSHPX,MSHTX,MS&
     &HLS,MSHLX, MSHLM,MSHM4,MSHM3,MSHTR,MSHP4,MSHP3,MSHBM,MSHSP,MSHDM, &
     &MSHSG,MSHDC,MSHTC,MSHSW,MSHWC,MSHBF,MSHPC,MSHFC,MSHSC, MSHVC,MSHCC&
     &,MSHNR,MSHND,MSHFS,MSHIT,MSHSI,MSHSN,MSHCE, MSHCN,MSHCX,MSHNS,MSHN&
     &E,MSHES,MSHEE,MSHSS,MSHSE,MSHMT, MSHLU,MSHTF,MSHFP,MSHRF,MSHPV,MSH&
     &S1,MSHS2,MSHG1,MSHG2, MSHG3,MSHMP,MSHRM,MSHCM,MSHIC,MSHRB,MSHRT,MS&
     &HST,MSHAX, MSHPP,MSHNC,MSHRE,MSHID

      INTEGER               MSHIF,MSHQA,MSHNP,MSHEL,MSHHX,MSHPX,MSHTX,MS&
     &HLS,MSHLX, MSHLM,MSHM4,MSHM3,MSHTR,MSHP4,MSHP3,MSHBM,MSHSP,MSHDM, &
     &MSHSG,MSHDC,MSHTC,MSHSW,MSHWC,MSHBF,MSHPC,MSHFC,MSHSC, MSHVC,MSHCC&
     &,MSHNR,MSHND,MSHFS,MSHIT,MSHSI,MSHSN,MSHCE, MSHCN,MSHCX,MSHNS,MSHN&
     &E,MSHES,MSHEE,MSHSS,MSHSE,MSHMT, MSHLU,MSHTF,MSHFP,MSHRF,MSHPV,MSH&
     &S1,MSHS2,MSHG1,MSHG2, MSHG3,MSHMP,MSHRM,MSHCM,MSHIC,MSHRB,MSHRT,MS&
     &HST,MSHAX, MSHPP,MSHNC,MSHRE,MSHID
!!
!! Counters for the records contained in the restart input. For ease of
!! management, the mesh counters exactly parallel the NUMxx counters. Not
!! all NRSxx counters have corresponding restart data input.
!!
      COMMON /BLOCK_03_NRS/ NRSIF,NRSQA,NRSNP,NRSEL,NRSHX,NRSPX,NRSTX,NR&
     &SLS,NRSLX, NRSLM,NRSM4,NRSM3,NRSTR,NRSP4,NRSP3,NRSBM,NRSSP,NRSDM, &
     &NRSSG,NRSDC,NRSTC,NRSSW,NRSWC,NRSBF,NRSPC,NRSFC,NRSSC, NRSVC,NRSCC&
     &,NRSNR,NRSND,NRSFS,NRSIT,NRSSI,NRSSN,NRSCE, NRSCN,NRSCX,NRSNS,NRSN&
     &E,NRSES,NRSEE,NRSSS,NRSSE,NRSMT, NRSLU,NRSTF,NRSFP,NRSRF,NRSPV,NRS&
     &S1,NRSS2,NRSG1,NRSG2, NRSG3,NRSMP,NRSRM,NRSCM,NRSIC,NRSRB,NRSRT,NR&
     &SST,NRSAX, NRSPP,NRSNC,NRSRE,NRSID

      INTEGER               NRSIF,NRSQA,NRSNP,NRSEL,NRSHX,NRSPX,NRSTX,NR&
     &SLS,NRSLX, NRSLM,NRSM4,NRSM3,NRSTR,NRSP4,NRSP3,NRSBM,NRSSP,NRSDM, &
     &NRSSG,NRSDC,NRSTC,NRSSW,NRSWC,NRSBF,NRSPC,NRSFC,NRSSC, NRSVC,NRSCC&
     &,NRSNR,NRSND,NRSFS,NRSIT,NRSSI,NRSSN,NRSCE, NRSCN,NRSCX,NRSNS,NRSN&
     &E,NRSES,NRSEE,NRSSS,NRSSE,NRSMT, NRSLU,NRSTF,NRSFP,NRSRF,NRSPV,NRS&
     &S1,NRSS2,NRSG1,NRSG2, NRSG3,NRSMP,NRSRM,NRSCM,NRSIC,NRSRB,NRSRT,NR&
     &SST,NRSAX, NRSPP,NRSNC,NRSRE,NRSID
!!
!! Counters for the maximum ID values in use for each record. For ease of
!! management, the maximum ID counters exactly parallel the NUMxx counters.
!! Not all IDXxx counters have corresponding record ID's to which they can
!! be related.
!!
      INTEGER               IDXIF,IDXQA,IDXNP,IDXEL,IDXHX,IDXPX,IDXTX,ID&
     &XLS,IDXLX, IDXLM,IDXM4,IDXM3,IDXTR,IDXP4,IDXP3,IDXBM,IDXSP,IDXDM, &
     &IDXSG,IDXDC,IDXTC,IDXSW,IDXWC,IDXBF,IDXPC,IDXFC,IDXSC, IDXVC,IDXCC&
     &,IDXNR,IDXND,IDXFS,IDXIT,IDXSI,IDXSN,IDXCE, IDXCN,IDXCX,IDXNS,IDXN&
     &E,IDXES,IDXEE,IDXSS,IDXSE,IDXMT, IDXLU,IDXTF,IDXFP,IDXRF,IDXPV,IDX&
     &S1,IDXS2,IDXG1,IDXG2, IDXG3,IDXMP,IDXRM,IDXCM,IDXIC,IDXRB,IDXRT,ID&
     &XST,IDXAX, IDXPP,IDXNC,IDXRE,IDXID

      COMMON /BLOCK_03_IDX/ IDXIF,IDXQA,IDXNP,IDXEL,IDXHX,IDXPX,IDXTX,ID&
     &XLS,IDXLX, IDXLM,IDXM4,IDXM3,IDXTR,IDXP4,IDXP3,IDXBM,IDXSP,IDXDM, &
     &IDXSG,IDXDC,IDXTC,IDXSW,IDXWC,IDXBF,IDXPC,IDXFC,IDXSC, IDXVC,IDXCC&
     &,IDXNR,IDXND,IDXFS,IDXIT,IDXSI,IDXSN,IDXCE, IDXCN,IDXCX,IDXNS,IDXN&
     &E,IDXES,IDXEE,IDXSS,IDXSE,IDXMT, IDXLU,IDXTF,IDXFP,IDXRF,IDXPV,IDX&
     &S1,IDXS2,IDXG1,IDXG2, IDXG3,IDXMP,IDXRM,IDXCM,IDXIC,IDXRB,IDXRT,ID&
     &XST,IDXAX, IDXPP,IDXNC,IDXRE,IDXID
!!
!! To the extent practical, all non-local variables are contained within
!! a derived data type.
!!
!! 1. USER's JOB_ID_RECORD: Current and restart file title, date, and time.
!!
      INTEGER, PARAMETER :: NPNTL = 80  ! Number of Characters in title

        TYPE :: CUR_type
          CHARACTER(NPNTL)    TITLE            ! Current job-title
          CHARACTER(12)       DATEE            ! Current date-of-execution
          CHARACTER(8)        TIMEE            ! Current time-of-execution
          INTEGER             SEQUENCE_NUMBER  ! Current restart sequence number
        END TYPE

        TYPE :: MSH_type
          CHARACTER(NPNTL)    TITLE            ! Mesh data title
          CHARACTER(12)       DATEE            ! Mesh data date-of-creation
          CHARACTER(8)        TIMEE            ! Mesh data time-of-creation
        END TYPE

        TYPE :: RST_type
          CHARACTER(NPNTL)    TITLE            ! Restart file job-title
          CHARACTER(12)       DATEE            ! Restart file date-of-execution
          CHARACTER(8)        TIMEE            ! Restart file time-of-execution
          INTEGER             SEQUENCE_NUMBER  ! Restart file sequence number
        END TYPE

      TYPE :: JOBIDR
        INTEGER       TITLE_LENGTH
        CHARACTER(8)  PROGRAM                 ! Program name
        CHARACTER(8)  VERSION                 ! Program version
        CHARACTER(32) RELEASE                 ! Program release date
        TYPE (CUR_type) :: CURRENT
        TYPE (MSH_type) :: MESH
        TYPE (RST_type) :: RESTART
      END TYPE
!!
!! 2. ANALYSIS CONTROL VARIABLES: Flags controling execution.
!!
      INTEGER, PARAMETER :: NPNCV = 14  ! Number of named entries

      TYPE :: EXECON
        INTEGER :: Number_of_Entries
        INTEGER          SPRINT  ! Skip unessential calculations, default: no
        INTEGER          INECHO  ! Write input echo, default: no
        INTEGER          WRSTAR  ! Write Restart, default: no
        INTEGER          RDSTAR  ! Read  Restart, default: no
        INTEGER          WRMESH  ! Write Mesh, default: no
        INTEGER          RDMESH  ! Read  Mesh, default: no
        INTEGER          DCHECK  ! Data check flag, default: no
        INTEGER          DTTABF  ! Time step tabulated ftn, default: none
        INTEGER          SUBCYC  ! Time integration Subcycling, default: no
        INTEGER          POLARD  ! Polar decomposition rotation, default: no
        INTEGER          MIDINT  ! Mid-interval gradient evaluation, default: no
        INTEGER          BTPLTQ  ! Belytscko-Tsay 4-node plate, default: no
        INTEGER          REZONE  ! Rezoning/remeshing processor: no
        INTEGER          RZSTAR  ! Rezoning/remeshing restart: no
        CHARACTER(12), DIMENSION(NPNCV) :: NAME
      END TYPE
!!
!! 3. SLIDING INTERFACE CONTROL: Override for sliding interface execution.
!!
      INTEGER, PARAMETER :: NPNIV =  3  ! Number of named entries
      INTEGER, PARAMETER :: NPNIT = 24  ! Number of NPNIV-tuples

      TYPE :: SI_type
        INTEGER          ID          ! Sliding interface ID
        REAL(KIND(0D0))  BEGIN       ! Start calculations
        REAL(KIND(0D0))  END         ! End calculations
      END TYPE

      TYPE :: INTERF
        INTEGER    :: Number_of_Interfaces
        INTEGER    :: Number_of_Entries
        CHARACTER(12),  DIMENSION(NPNIV) :: NAME
        TYPE (SI_type), DIMENSION(NPNIT) :: SI
      END TYPE
!!
!! 4. NUMERICAL PROCEDURE CONSTANTS: User adjustable coefficients.
!!
      INTEGER, PARAMETER :: NPNPA = 27  ! Number of named entries

      TYPE :: PARAMS
        REAL(KIND(0D0))  Blk1            ! Linear bulk viscosity coefficeint
        REAL(KIND(0D0))  Blk2            ! Quad.  bulk viscosity coefficeint
        REAL(KIND(0D0))  HGV             ! Hourglass viscosity
        REAL(KIND(0D0))  HGK             ! Hourglass stiffness
        REAL(KIND(0D0))  DTratio         ! Limit ratio, DTmin(t)/Dtmin(0)
        REAL(KIND(0D0))  DTscale         ! Scaling factor for critical Dt
        REAL(KIND(0D0))  DTgrow          ! Step to step growth allowed in DTnext
        REAL(KIND(0D0))  STATUS          ! Status request checking interval
        REAL(KIND(0D0))  ENG_BAL         ! Energy balance reporting interval
        REAL(KIND(0D0))  Factor          ! Sliding Interface relaxation factor
        REAL(KIND(0D0))  Capture         ! Sliding Interface capture distance
        REAL(KIND(0D0))  Border          ! Contact element "exterior" border
        REAL(KIND(0D0))  Sort_Freq       ! Sliding Interface sorting frequency
        REAL(KIND(0D0))  Cycle_R         ! Subcycling group ratio
        REAL(KIND(0D0))  Cycle_S         ! Subcycling interface spreading
        REAL(KIND(0D0))  Cycle_F         ! Subcycling repartitioning frequency
        REAL(KIND(0D0))  Max_Steps       ! Material 32 substep upper limit
        REAL(KIND(0D0))  NRBC_Q          ! Nonreflecting BC viscous scaling
        REAL(KIND(0D0))  EQUI_TOL        ! Layered solid equilibrium tolerance
        REAL(KIND(0D0))  FIBROT          ! Material 22 fiber rotation form
        REAL(KIND(0D0))  BCSMXIT         ! Beam cross section maximum iterations
        REAL(KIND(0D0))  BCSGRLM         ! Beam cross section iterate growth lim
        REAL(KIND(0D0))  BCSRLAX         ! Beam cross section under relaxtion
        REAL(KIND(0D0))  CPDMAX          ! Maximum rate of change for Cos(Phi)
        REAL(KIND(0D0))  CPDMIN          ! Minimum rate of change for Cos(Phi)
        REAL(KIND(0D0))  CSPMAX          ! Maximum change for Cos(Phi)
        REAL(KIND(0D0))  CSPMIN          ! Minimum change for Cos(Phi)
        INTEGER :: Number_of_Entries
        CHARACTER(12), DIMENSION(NPNPA) :: NAME
      END TYPE
!!
!! 5. SOUND SPEED: Element longitudinal and transverse sound speeds
!!
      TYPE :: SSPDAT
        REAL(KIND(0D0))  Density ! Current material density
        REAL(KIND(0D0))  RCL2    ! Rho * (Longitudinal_wave_speed)**2
        REAL(KIND(0D0))  RCS2    ! Rho * (Transverse_shear_wave_speed)**2
      END TYPE
!!
!! 6. SIMULATION TIME: Simulation time and associated critical time step info
!!
      INTEGER, PARAMETER :: NPNDT = 10  ! Number of minimum DT's retained

      TYPE :: TIMESIMULATION
        REAL(KIND(0D0))  Stop            ! End of simulation
        REAL(KIND(0D0))  Total           ! Total simulation time (master clock)
        INTEGER          Step            ! Time step counter
        INTEGER          Cycle           ! Subcycling cycle counter
        INTEGER          Cymax           ! Maximum number of subcycle steps
        INTEGER          DTcntrl         ! Time step control function
        REAL(KIND(0D0))  DTnext          ! Next time step (master time step)
        REAL(KIND(0D0))  DTlast          ! Last time step
        REAL(KIND(0D0))  DTmin           ! Minimum over all elements
        REAL(KIND(0D0))  DTmax           ! Maximum over all elements
        REAL(KIND(0D0))  DTlim           ! Smallest permissable time step
        REAL(KIND(0D0))  DTsav           ! Last element minimum time step
        CHARACTER(32)    Source          ! Source of minimum time step
        REAL(KIND(0D0))  DTHxx           ! Hexahedron  DTMax, (8 nodes)
        REAL(KIND(0D0))  DTPnx           ! Pentahedron DTMax, (6 nodes)
        REAL(KIND(0D0))  DTTtx           ! Tetrahedron DTMax, (4 nodes)
        REAL(KIND(0D0))  DTLSx           ! Lay. Solid  DTMax, (8 nodes)
        REAL(KIND(0D0))  DTM3x           ! Membrane  DTMax,   (3 nodes)
        REAL(KIND(0D0))  DTM4x           ! Membrane  DTMax,   (4 nodes)
        REAL(KIND(0D0))  DTTrx           ! Truss     DTMax,   (2 nodes)
        REAL(KIND(0D0))  DTP3x           ! Plate     DTMax,   (3 nodes)
        REAL(KIND(0D0))  DTP4x           ! Plate     DTMax,   (4 nodes)
        REAL(KIND(0D0))  DTBmx           ! Beam      DTMax,   (2 nodes)
        REAL(KIND(0D0))  DTSpx           ! Spring    DTMax,   (2 nodes)
        REAL(KIND(0D0))  DTDmx           ! Damper    DTMax,   (2 nodes)
        REAL(KIND(0D0))  DTSBx           ! Spring BC DTMax,   (1 node )
        REAL(KIND(0D0))  DTDBx           ! Damper BC DTMax,   (1 node )
        REAL(KIND(0D0))  DTHex(NPNDT)    ! Hexahedron  DTCrit, (8 nodes)
        REAL(KIND(0D0))  DTPen(NPNDT)    ! Pentahedron DTCrit, (6 nodes)
        REAL(KIND(0D0))  DTTet(NPNDT)    ! Tetrahedron DTCrit, (4 nodes)
        REAL(KIND(0D0))  DTLYS(NPNDT)    ! Lay. Solid  DTCrit, (8 nodes)
        REAL(KIND(0D0))  DTMb3(NPNDT)    ! Membrane  DTCrit,   (3 nodes)
        REAL(KIND(0D0))  DTMb4(NPNDT)    ! Membrane  DTCrit,   (4 nodes)
        REAL(KIND(0D0))  DTTru(NPNDT)    ! Truss     DTCrit,   (2 nodes)
        REAL(KIND(0D0))  DTPl3(NPNDT)    ! Plate     DTCrit,   (3 nodes)
        REAL(KIND(0D0))  DTPl4(NPNDT)    ! Plate     DTCrit,   (4 nodes)
        REAL(KIND(0D0))  DTBms(NPNDT)    ! Beam      DTCrit,   (2 nodes)
        REAL(KIND(0D0))  DTSpr(NPNDT)    ! Spring    DTCrit,   (2 nodes)
        REAL(KIND(0D0))  DTDmp(NPNDT)    ! Damper    DTCrit,   (2 nodes)
        REAL(KIND(0D0))  DTSBC(NPNDT)    ! Spring BC DTCrit,   (1 node )
        REAL(KIND(0D0))  DTDBC(NPNDT)    ! Damper BC DTCrit,   (1 node )
        INTEGER          Ncrit           ! Element ID w/ critical time step
        INTEGER          Hexah(NPNDT)    ! Hexahedron  IDs, (8 nodes)
        INTEGER          Penta(NPNDT)    ! Pentahedron IDs, (6 nodes)
        INTEGER          Tetra(NPNDT)    ! Tetrahedron IDs, (4 nodes)
        INTEGER          LSold(NPNDT)    ! Lay. Solid  IDs, (8 nodes)
        INTEGER          Memb3(NPNDT)    ! Membrane IDs,    (3 nodes)
        INTEGER          Memb4(NPNDT)    ! Membrane IDs,    (4 nodes)
        INTEGER          Truss(NPNDT)    ! Truss  IDs,      (2 nodes)
        INTEGER          Plat3(NPNDT)    ! Plate  IDs,      (3 nodes)
        INTEGER          Plat4(NPNDT)    ! Plate  IDs,      (4 nodes)
        INTEGER          Beams(NPNDT)    ! Beam   IDs,      (2 nodes)
        INTEGER          Spring(NPNDT)   ! Spring IDs,      (2 nodes)
        INTEGER          Damper(NPNDT)   ! Damper IDs,      (2 nodes)
        INTEGER          SprBC(NPNDT)    ! Spring BC IDs,   (1 node )
        INTEGER          DmpBC(NPNDT)    ! Damper BC IDs,   (1 node )
      END TYPE
!!
!! 7. ENERGY BALANCE: Energy balance terms.
!!
      TYPE :: ENGDAT
        REAL(KIND(0D0))  External  ! External energy
        REAL(KIND(0D0))  Internal  ! Internal energy, material
        REAL(KIND(0D0))  Sliding   ! Internal energy, tangential friction
        REAL(KIND(0D0))  Contact   ! Internal energy, normal contact
        REAL(KIND(0D0))  Hourglass ! Internal energy, hourglass
        REAL(KIND(0D0))  Bulk_Vis  ! Internal energy, bulk viscosity
        REAL(KIND(0D0))  Kinetic   ! Kinetic energy
        REAL(KIND(0D0))  Eng_Sum   ! Kinetic + Internal
        REAL(KIND(0D0))  Balance   ! Energy balance
      END TYPE
!!
!! 8. IO_UNIT NUMBERS: Integers assigned to I/O units. The use of terse,
!! lower-case, incomprehensible file names such as "fmasdi" is for
!! compatability with the mono-syllable argot of unix-speak and the
!! unix file system directory structure and file naming conventions.
!!
      TYPE :: IOUNTS

        ! Input files, logical unit numbers
        INTEGER      LSDI ! fmasdi, simulation data input      (ASCII)
        INTEGER      LSRI ! fmasri, status request input       (ASCII)
        INTEGER      LMDI ! fmamdi, mesh data input            (binary)
        INTEGER      LRDI ! fmardi, restart data input         (binary)

        ! Output files, logical unit numbers
        INTEGER      LELO ! fmaelo, execution log output       (ASCII)
        INTEGER      LSDO ! fmasdo, simulation data output     (ASCII)
        INTEGER      LCRO ! fmacro, computed results output    (ASCII)
        INTEGER      LSRO ! fmasro, status request output      (ASCII)
        INTEGER      LMDO ! fmamdo, mesh data output           (binary)
        INTEGER      LRDO ! fmardo, restart data output        (binary)
        INTEGER      LPDB ! fmapdb, plotting database          (binary)

        ! Rezone restart file, logical unit number
        INTEGER      LRZD ! FMARZD, rezone restart data out/in (binary)

      END TYPE
!!
!! 9. PRINTED OUTPUT REQUESTS:
!!
      INTEGER, PARAMETER :: NPNPF = 16  ! Number of named entries

      TYPE :: PRINTF
        INTEGER :: Number_of_Entries
        CHARACTER(6), DIMENSION(NPNPF) :: NAME
        REAL(KIND(0D0))  Begin   ! Begin writing to printer file
        REAL(KIND(0D0))  End     ! End writing to printer file
        REAL(KIND(0D0))  Delta   ! Time increment between writes
        REAL(KIND(0D0))  Time    ! Next time to write
        INTEGER          NODES   ! Flag (-1/n=ALL/Nodal Point SetID)
        INTEGER          HEXAH   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          PENTA   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          TETRA   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          MEMBT   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          MEMBQ   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          TRUSS   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          PLATT   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          PLATQ   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          BEAMS   ! Flag (-1/n=ALL/Element SetID)
        INTEGER          SPRING  ! Flag (-1/n=ALL/Element SetID)
        INTEGER          DAMPER  ! Flag (-1/n=ALL/Element SetID)
      END TYPE
!!
!! 10. OPERATION COUNTERS: Counters for number of times each element is executed
!!
      TYPE :: ECOUNT
        INTEGER          HEXAH
        INTEGER          PENTA
        INTEGER          TETRA
        INTEGER          LSOLD
        INTEGER          MEMBT
        INTEGER          MEMBQ
        INTEGER          TRUSS
        INTEGER          PLATT
        INTEGER          PLATQ
        INTEGER          BEAM
        INTEGER          SPRING
        INTEGER          DAMPER
        INTEGER          GAUGE1D
        INTEGER          GAUGE2D
        INTEGER          GAUGE3D
      END TYPE
!!
!! 11. INPUT RECORD EORRORS: Counter for number of errors on input records.
!!
      TYPE :: IOERRS
        INTEGER          COUNT   ! Count of input record errors.
      END TYPE
!!
!! 12. REZONING PROCEDURE CONSTANTS: User adjustable coefficients.
!!
      INTEGER, PARAMETER :: NPNRZ = 13  ! Number of named entries

      TYPE :: RZMESH
        INTEGER    :: Number_of_Entries
        CHARACTER(12), DIMENSION(NPNRZ) :: NAME
        INTEGER          INTERVAL    ! Master time step check interval
        INTEGER          TMSTEP      ! Time step at which first flag set
        INTEGER          FCOUNT      ! Flag counter
        INTEGER          FLGCYC      ! Flag-cycles counter
        INTEGER          MAXCNT      ! Maximum flags before rezone
        INTEGER          MAXFCY      ! Maximum flag-cycles
        INTEGER          MAXLEV      ! Maximum rezone levels
        INTEGER          WRMESH      ! Write "rezone" mesh
        LOGICAL          FLAG_IS_SET ! Rezoning control flag
        LOGICAL          FIRST       ! Flag for starting new statistics
        LOGICAL          REPARTITION ! Flag to force a repartition
        LOGICAL          WRDATA      ! Flag for writing rezone data
        LOGICAL          FAILED      ! Flag to indicate array overflow
      END TYPE
!!
!! Module TIMER is used to accumulate subprocess time increments in the
!! array SUBPROCESS_TIME. SUBPROCESS_TIME is carried here to make it
!! available for restart data reads and writes.
!!
      REAL(KIND(0D0))  SUBPROCESS_TIME(100)
!!
!! The character variable MSGL must be used to start each line of a message
!! written by USER_MESSAGE. The character variables MSG1, MSG2 and MSGF are
!! internal buffers for formatting data for use in constructing messages.
!!
      CHARACTER MSGL*1, MSG1*8, MSG2*8, MSGF*12

      TYPE (JOBIDR) :: JOB_ID_RECORD
      TYPE (EXECON) :: CONTROL
      TYPE (INTERF) :: INTERFACE_TIME
      TYPE (PARAMS) :: PARAMVALUE
      TYPE (SSPDAT) :: SOUND_SPEED
      TYPE (TIMESIMULATION) :: TIMSIM
      TYPE (ENGDAT) :: ENERGY
      TYPE (IOUNTS) :: IO_UNIT
      TYPE (PRINTF) :: PRINT
      TYPE (ECOUNT) :: COUNTER
      TYPE (IOERRS) :: ERROR
      TYPE (RZMESH) :: REZONE

      END MODULE shared_common_data
