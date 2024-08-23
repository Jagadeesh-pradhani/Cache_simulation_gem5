!!
!! 14. SECTION PROPERTY DATA: Beam and truss section structures.
!!
      MODULE section_1d_                                                
                                                                        
      TYPE :: section_1d_type                                           
        INTEGER          SecID    ! Section ID (User defined value retained)
        INTEGER          Section  ! Beam cross section number (1,2,3,...)  
        INTEGER          Iprop    ! Geometric properties (0/1=dimension/inertia)
        INTEGER          NPLoc    ! Nodal point location (0,1,2,3,4,5)     
        REAL(KIND(0D0))  Width    ! Cross section width            (Iprop=0)
        REAL(KIND(0D0))  Height   ! Cross section height           (Iprop=0)
        REAL(KIND(0D0))  Twall    ! Wall or web thickness          (Iprop=0)
        REAL(KIND(0D0))  Tflange  ! Flange thickness (I-beam)      (Iprop=0)
        REAL(KIND(0D0))  Tcover   ! Cover thickness (hat-section)  (Iprop=0)
        REAL(KIND(0D0))  Area     ! Cross section area             (Iprop=1)
        REAL(KIND(0D0))  Br       ! I(r**2)da inertia about x-axis         
        REAL(KIND(0D0))  By       ! I(y**2)da inertia about z-axis (Iprop=1)
        REAL(KIND(0D0))  Bz       ! I(z**2)da inertia about y-axis (Iprop=1)
        REAL(KIND(0D0))  Yrefloc  ! Reference axis location bcsys  (NPLoc=5)
        REAL(KIND(0D0))  Zrefloc  ! Reference axis location bcsys  (NPLoc=5)
      END TYPE                                                          
                                                                        
      TYPE (section_1d_type), DIMENSION(:), ALLOCATABLE :: SECTION_1D   
                                                                        
      END MODULE section_1d_
