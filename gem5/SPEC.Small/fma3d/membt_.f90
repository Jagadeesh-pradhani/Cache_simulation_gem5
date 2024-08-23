!!
!! 49. TRIANGULAR MEMBRANE: Element parameter and results structures.
!!
      MODULE membt_                                                     
                                                                        
      TYPE :: PAR_membt                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          SecID    ! Section  ID  (Reset: points to SECTION_2D)
        INTEGER          IX(3)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Area     ! Initial area                           
      END TYPE                                                          
                                                                        
      TYPE :: RES_membt                                                 
        REAL(KIND(0D0))  Area     ! Current area                           
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Beta     ! Angle in polar decomp of deformation grad.
        REAL(KIND(0D0))  Stress(3)! In-plane stress (plane stress)         
        REAL(KIND(0D0))  Xint(3)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(3)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(3)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: membt_type                                                
        TYPE (PAR_membt) :: PAR                                         
        TYPE (RES_membt) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (membt_type), DIMENSION(:), ALLOCATABLE :: MEMBT             
                                                                        
      END MODULE membt_
