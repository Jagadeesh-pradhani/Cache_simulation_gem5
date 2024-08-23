!!
!! 56. AXIAL/TORSIONAL DAMPER: Element parameters and results.
!!
      MODULE damper_                                                    
                                                                        
      TYPE :: PAR_damper                                                
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          Type     ! Axial or torsional (0/1=axial/torsional)
        INTEGER          IX(2)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Idir     ! Direction of action (0/1/2/3/4=i,j/x/y/z/a)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Axis(3)  ! Action direction, (Idir=4)             
      END TYPE                                                          
                                                                        
      TYPE :: RES_damper                                                
        REAL(KIND(0D0))  Delta    ! Relative displacement/twist            
        REAL(KIND(0D0))  Int_Eng  ! Internal energy                        
        REAL(KIND(0D0))  Force    ! Axial force/torque                     
        REAL(KIND(0D0))  Axis(3)  ! Action direction                       
        REAL(KIND(0D0))  Xint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(2)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: damper_type                                               
        TYPE (PAR_damper) :: PAR                                        
        TYPE (RES_damper) :: RES                                        
      END TYPE                                                          
                                                                        
      TYPE (damper_type), DIMENSION(:), ALLOCATABLE :: DAMPER           
                                                                        
      END MODULE damper_
