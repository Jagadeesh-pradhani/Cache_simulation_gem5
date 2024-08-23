!!
!! 58. VELOCITY IC: Initial condition data structure.
!!
      MODULE velocity_ic_                                               
                                                                        
      TYPE :: velocity_ic_type                                          
        INTEGER          ICID    ! Nodal point ID (User defined value retained)
        INTEGER          NPID    ! Nodal point ID (Reset: points to internal ID)
        REAL(KIND(0D0))  Vx      ! X-component of translational velocity   
        REAL(KIND(0D0))  Vy      ! Y-component of translational velocity   
        REAL(KIND(0D0))  Vz      ! Z-component of translational velocity   
        REAL(KIND(0D0))  Ox      ! X-component of rotational velocity      
        REAL(KIND(0D0))  Oy      ! Y-component of rotational velocity      
        REAL(KIND(0D0))  Oz      ! Z-component of rotational velocity      
      END TYPE                                                          
                                                                        
      TYPE (velocity_ic_type), DIMENSION(:), ALLOCATABLE :: VELOCITY_IC 
                                                                        
      END MODULE velocity_ic_
