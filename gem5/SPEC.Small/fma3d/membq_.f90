!!
!! 47. QUADRILATERAL MEMBRANE: Element parameter and results structures.
!!
      MODULE membq_                                                     
                                                                        
      TYPE :: PAR_membq                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_2D)
        INTEGER          IX(4)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Area     ! Initial area                           
      END TYPE                                                          
                                                                        
      TYPE :: RES_membq                                                 
        REAL(KIND(0D0))  Area     ! Current area                           
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Beta     ! Angle in polar decomp of deformation grad.
        REAL(KIND(0D0))  Stress(3)! In-plane stress (plane stress)         
        REAL(KIND(0D0))  Pr       ! Hourglass control forces, r-component  
        REAL(KIND(0D0))  Ps       ! Hourglass control forces, s-component  
        REAL(KIND(0D0))  Pt       ! Hourglass control forces, t-component  
        REAL(KIND(0D0))  Xint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: membq_type                                                
        TYPE (PAR_membq) :: PAR                                         
        TYPE (RES_membq) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (membq_type), DIMENSION(:), ALLOCATABLE :: MEMBQ             
!!
!! LAYERED SOLID QUADRILATERAL: Element parameter and results structures.
!! Note that this structure must "match" the MEMBQ structure since it is
!! passed to MEMBQ subroutines for processing!!
!!
      TYPE (membq_type), DIMENSION(:), ALLOCATABLE :: LSMBQ             
                                                                        
      END MODULE membq_
