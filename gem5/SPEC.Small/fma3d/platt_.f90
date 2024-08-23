!!
!! 52. TRIANGULAR PLATE: Element parameter and results structures.
!!
      MODULE platt_                                                     
                                                                        
      TYPE :: PAR_platt                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_2D)
        INTEGER          IX(6)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Ist      ! Starting location for stress           
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Area     ! Initial area                           
      END TYPE                                                          
                                                                        
      TYPE :: RES_platt                                                 
        REAL(KIND(0D0))  Area     ! Current area                           
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Xint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(6)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: platt_type                                                
        TYPE (PAR_platt) :: PAR                                         
        TYPE (RES_platt) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (platt_type), DIMENSION(:), ALLOCATABLE :: PLATT             
                                                                        
      END MODULE platt_
