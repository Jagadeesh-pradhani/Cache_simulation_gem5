!!
!! 53. BEAM: Element parameters and results.
!!
      MODULE beam_                                                      
                                                                        
      TYPE :: PAR_beam                                                  
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_1D)
        INTEGER          IX(4)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Length   ! Initial length                         
      END TYPE                                                          
                                                                        
      TYPE :: RES_beam                                                  
        REAL(KIND(0D0))  Length   ! Current length                         
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Zaxis(3) ! Cross section orientation              
        REAL(KIND(0D0))  Axial(16)! Axial stress                           
        REAL(KIND(0D0))  Shear(16)! Shear stress                           
        REAL(KIND(0D0))  Trs      ! Transverse shear stress                
        REAL(KIND(0D0))  Trt      ! Transverse shear stress                
        REAL(KIND(0D0))  Xint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(4)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: beam_type                                                 
        TYPE (PAR_beam) :: PAR                                          
        TYPE (RES_beam) :: RES                                          
      END TYPE                                                          
                                                                        
      TYPE (beam_type), DIMENSION(:), ALLOCATABLE :: BEAM               
                                                                        
      END MODULE beam_
