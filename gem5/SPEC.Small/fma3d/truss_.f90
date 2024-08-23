!!
!! 50. AXIAL FORCE TRUSS: Element parameter and results structures.
!!
      MODULE truss_                                                     
                                                                        
      TYPE :: PAR_truss                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_1D)
        INTEGER          IX(2)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Length   ! Initial length                         
      END TYPE                                                          
                                                                        
      TYPE :: RES_truss                                                 
        REAL(KIND(0D0))  Length   ! Current length                         
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Stress   ! Axial stress (uniaxial stress)         
        REAL(KIND(0D0))  Xint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: truss_type                                                
        TYPE (PAR_truss) :: PAR                                         
        TYPE (RES_truss) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (truss_type), DIMENSION(:), ALLOCATABLE :: TRUSS             
                                                                        
      END MODULE truss_
