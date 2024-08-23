!!
!! 44. TETRAHEDRON: Element parameter and results structures.
!!
      MODULE tetra_                                                     
                                                                        
      TYPE :: PAR_tetra                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          LupID    ! Lay-up ID    (Reset: points to LAYERING)
        INTEGER          IX(4)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Volume   ! Initial volume                         
      END TYPE                                                          
                                                                        
      TYPE :: RES_tetra                                                 
        REAL(KIND(0D0))  Volume   ! Current volume                         
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Stress(6)! Stress                                 
        REAL(KIND(0D0))  Xint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: tetra_type                                                
        TYPE (PAR_tetra) :: PAR                                         
        TYPE (RES_tetra) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (tetra_type), DIMENSION(:), ALLOCATABLE :: TETRA             
                                                                        
      END MODULE tetra_
