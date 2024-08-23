!!
!! 31. SLIDING INTERFACE
!!
      MODULE contact_node_                                              
                                                                        
      TYPE :: contact_node_type                                         
        INTEGER          NPID    ! Nodal point ID of contact surface node  
        REAL(KIND(0D0))  Mass    ! Mass of contact surface node            
        REAL(KIND(0D0))  Ax      ! Acceleration component                  
        REAL(KIND(0D0))  Ay      ! Acceleration component                  
        REAL(KIND(0D0))  Az      ! Acceleration component                  
        REAL(KIND(0D0))  Vx      ! Velocity component                      
        REAL(KIND(0D0))  Vy      ! Velocity component                      
        REAL(KIND(0D0))  Vz      ! Velocity component                      
        REAL(KIND(0D0))  Px      ! Position component                      
        REAL(KIND(0D0))  Py      ! Position component                      
        REAL(KIND(0D0))  Pz      ! Position component                      
      END TYPE                                                          
                                                                        
      TYPE (contact_node_type), &
     &  DIMENSION(:), ALLOCATABLE :: CONTACT_NODE                                                                 
                                                                        
      TYPE (contact_node_type), DIMENSION(4) :: TRIANGLE_NODE           
                                                                        
      END MODULE contact_node_
