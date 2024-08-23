!!
!! 33. SLIDING INTERFACE
!!
      MODULE contact_surface_                                           
                                                                        
      TYPE :: contact_surface_type                                      
        INTEGER          IX(4)   ! Contact segment definition              
        INTEGER          NX(4)   ! Neighboring contact segment pointer     
        REAL(KIND(0D0))  An(3)   ! Contact unit segment normal vector      
        REAL(KIND(0D0))  Area    ! Contact segment area                    
        REAL(KIND(0D0))  Rsqd    ! Radius squared of the enclosing sphere  
        REAL(KIND(0D0))  Xave    ! X-position of segment center            
        REAL(KIND(0D0))  Yave    ! Y-position of segment center            
        REAL(KIND(0D0))  Zave    ! Z-position of segment center            
      END TYPE                                                          
                                                                        
      TYPE (contact_surface_type), &
     &  DIMENSION(:), ALLOCATABLE :: CONTACT_SURFACE                                                           
                                                                        
      TYPE (contact_surface_type) :: TRIANGLE_SURFACE                   
                                                                        
      END MODULE contact_surface_
