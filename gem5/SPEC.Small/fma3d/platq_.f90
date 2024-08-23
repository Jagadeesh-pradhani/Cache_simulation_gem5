!!
!! 51. QUADRILATERAL PLATE: Element parameter and results structures.
!!
      MODULE platq_                                                     
                                                                        
      TYPE :: PAR_platq                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_2D)
        INTEGER          IX(8)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          Ist      ! Starting location for stress           
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Area     ! Initial area                           
      END TYPE                                                          
                                                                        
      TYPE :: RES_platq                                                 
        REAL(KIND(0D0))  Area     ! Current area                           
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Shear(2) ! Transeverse shear                      
        REAL(KIND(0D0))  Pr(2)    ! Hourglass control forces, r-components 
        REAL(KIND(0D0))  Ps(2)    ! Hourglass control forces, s-components 
        REAL(KIND(0D0))  Pt(2)    ! Hourglass control forces, t-components 
        REAL(KIND(0D0))  Xint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: ADP_platq                                                 
        LOGICAL          REFINE   ! Flag to indicate rezoning is needed    
        INTEGER          LEVEL    ! Rezone level                           
        REAL(KIND(0D0))  Ax,Ay,Az ! Element normal vector                  
      END TYPE                                                          
                                                                        
      TYPE :: platq_type                                                
        TYPE (PAR_platq) :: PAR                                         
        TYPE (RES_platq) :: RES                                         
        TYPE (ADP_platq) :: ADP                                         
      END TYPE                                                          
                                                                        
      TYPE (platq_type), DIMENSION(:), ALLOCATABLE :: PLATQ             
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (platq_type), DIMENSION(:), ALLOCATABLE :: RZ_PLATQ          
                                                                        
      END MODULE platq_
