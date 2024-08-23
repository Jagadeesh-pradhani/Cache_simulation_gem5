!!
!! 18. DISPLACEMENT BC
!!
      MODULE displacement_bc_                                           
                                                                        
      TYPE :: displacement_bc_type                                      
        INTEGER          DBCID   ! Displacement BC ID                      
        INTEGER          SetID   ! Node set ID (-n/0/n=-NPID/all/Set ID)   
        INTEGER          Code    ! Constraint code                         
        INTEGER          HstID   ! History ID (tabulated function ID)      
        INTEGER          Kavd    ! Constrained kinematic variable (1/2/3=a/v/d)
        REAL(KIND(0D0))  Scale   ! History function scale factor           
        REAL(KIND(0D0))  Ax      ! BC direction, x-component (Code=4,40,   
        REAL(KIND(0D0))  Ay      ! BC direction, y-component  44,8,80,88)  
        REAL(KIND(0D0))  Az      ! BC direction, z-component               
      END TYPE                                                          
                                                                        
      TYPE (displacement_bc_type), &
     &  DIMENSION(:), ALLOCATABLE :: DISPLACEMENT_BC                                                           
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (displacement_bc_type), &
     &  DIMENSION(:), ALLOCATABLE :: RZ_DISPLACEMENT_BC                                                        
                                                                        
      END MODULE displacement_bc_
