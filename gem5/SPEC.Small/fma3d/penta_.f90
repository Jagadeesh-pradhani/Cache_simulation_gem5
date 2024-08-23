!!
!! 43. PENTAHEDRON: Element parameter and results structures.
!!
      MODULE penta_                                                     
                                                                        
      TYPE :: PAR_penta                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          LupID    ! Lay-up ID    (Reset: points to LAYERING)
        INTEGER          IX(6)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Volume   ! Initial volume                         
      END TYPE                                                          
                                                                        
      TYPE :: RES_penta                                                 
        REAL(KIND(0D0))  Volume   ! Current volume                         
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Stress(6)! Stress                                 
        REAL(KIND(0D0))  Px(4)    ! Hourglass control forces, x-components 
        REAL(KIND(0D0))  Py(4)    ! Hourglass control forces, y-components 
        REAL(KIND(0D0))  Pz(4)    ! Hourglass control forces, z-components 
        REAL(KIND(0D0))  Xint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: penta_type                                                
        TYPE (PAR_penta) :: PAR                                         
        TYPE (RES_penta) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (penta_type), DIMENSION(:), ALLOCATABLE :: PENTA             
                                                                        
      END MODULE penta_
