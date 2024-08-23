!!
!! 24. CONCENTRATED FORCE BC:
!!
      MODULE force_bc_                                                  
                                                                        
      TYPE :: force_bc_type                                             
        INTEGER          CFID    ! Concentrated force ID                   
        INTEGER          SetID   ! Node set ID (-n/0/n=-NPID/all/Set ID)   
        INTEGER          HstID   ! History ID (force/torque = mag*f(t))    
        INTEGER          Type    ! Force or torque (0/1=force/torque)      
        INTEGER          Follow  ! Flag for follower force (0/1=no/yes)    
        REAL(KIND(0D0))  Force   ! Concentrated force/torque magnitude     
        REAL(KIND(0D0))  Cx      ! X-component of unit direction vector    
        REAL(KIND(0D0))  Cy      ! Y-component of unit direction vector    
        REAL(KIND(0D0))  Cz      ! Z-component of unit direction vector    
        REAL(KIND(0D0))  Delay   ! Delay time to use with tabulated function
      END TYPE                                                          
                                                                        
      TYPE (force_bc_type), DIMENSION(:), ALLOCATABLE :: FORCE_BC       
                                                                        
      END MODULE force_bc_
