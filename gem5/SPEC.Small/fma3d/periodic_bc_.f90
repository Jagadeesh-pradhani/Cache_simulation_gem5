!!
!! 27. PERIODIC BC: Periodic boundary condition, "repeated surfaces"
!!
      MODULE periodic_bc_                                               
                                                                        
      TYPE :: periodic_bc_type                                          
        INTEGER          PerID       ! Periodic ID (User defined value retained)
        INTEGER          Typ1        ! Set Type (0/1=segment/node)         
        INTEGER          S1ID        ! Node set ID (-n/n=-NPID/Set ID)     
        INTEGER          Typ2        ! Set Type (0/1=segment/node)         
        INTEGER          S2ID        ! Node set ID (-n/n=-NPID/Set ID)     
        CHARACTER(8)     Type        ! LINEAR or CYCLIC repeat             
        REAL(KIND(0D0))  Axis(3)     ! Rotation axis orientation           
        REAL(KIND(0D0))  Origin(3)   ! Rotation axis location              
        REAL(KIND(0D0))  Theta       ! Rotation angle, degrees             
        REAL(KIND(0D0))  Advance     ! Helical advance per degree          
        REAL(KIND(0D0))  RTX(3,3)    ! Orthogonal rotation operator, derived
      END TYPE                                                          
                                                                        
      TYPE (periodic_bc_type), DIMENSION(:), ALLOCATABLE :: PERIODIC_BC 
                                                                        
      END MODULE periodic_bc_
