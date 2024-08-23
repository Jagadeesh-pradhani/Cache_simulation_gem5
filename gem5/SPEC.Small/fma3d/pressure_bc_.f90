!!
!! 23. PRESSURE BC
!!
      MODULE pressure_bc_                                               
                                                                        
      TYPE :: pressure_bc_type                                          
        INTEGER          PBCID   ! Pressure BC ID                          
        INTEGER          SetID   ! Segment set ID (-n/0/n=-SGID/all/Set ID)
        INTEGER          HstID   ! History ID (tabulated function ID)      
        REAL(KIND(0D0))  Scale   ! History function scale factor           
        REAL(KIND(0D0))  PI      ! Pressure multiplier at node I/ on Seg Set
        REAL(KIND(0D0))  PJ      ! Pressure multiplier at node J           
        REAL(KIND(0D0))  PK      ! Pressure multiplier at node K           
        REAL(KIND(0D0))  PL      ! Pressure multiplier at node L           
        REAL(KIND(0D0))  Delay   ! Pressure arrival time                   
      END TYPE                                                          
                                                                        
      TYPE (pressure_bc_type), DIMENSION(:), ALLOCATABLE :: PRESSURE_BC 
                                                                        
      END MODULE pressure_bc_
