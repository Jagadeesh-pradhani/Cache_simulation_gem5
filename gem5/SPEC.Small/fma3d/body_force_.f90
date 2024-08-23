!!
!! 22. BODY FORCE
!!
      MODULE body_force_                                                
                                                                        
      TYPE :: body_force_type                                           
        INTEGER          BFID    ! Pressure BC ID                          
        INTEGER          SetID   ! Nodal point set ID (-n/n=-NPID/Set ID)  
        INTEGER          HstID   ! History ID (tabulated function ID)      
        REAL(KIND(0D0))  Scale   ! History function scale factor           
        REAL(KIND(0D0))  Gravity ! Gravitational constant                  
        REAL(KIND(0D0))  Ax      ! Direction in which gravity acts, x-component
        REAL(KIND(0D0))  Ay      ! Direction in which gravity acts, y-component
        REAL(KIND(0D0))  Az      ! Direction in which gravity acts, z-component
        REAL(KIND(0D0))  Delay   ! Body force arrival time                 
      END TYPE                                                          
                                                                        
      TYPE (body_force_type), DIMENSION(:), ALLOCATABLE :: BODY_FORCE   
                                                                        
      END MODULE body_force_
