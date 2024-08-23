!!
!! 20. SPOT WELD
!!
      MODULE spot_weld_                                                 
                                                                        
      TYPE :: spot_weld_type                                            
        INTEGER          SWID       ! Spot weld ID                         
        CHARACTER(8)     PLACE      ! Keyword "NODE" or "POSITION"         
        INTEGER          NPID(3)    ! Node point ID's                      
        LOGICAL          Coincident ! Flag for co-located nodes            
        REAL(KIND(0D0))  Xcm        ! X-coordinate of center of mass       
        REAL(KIND(0D0))  Ycm        ! Y-coordinate of center of mass       
        REAL(KIND(0D0))  Zcm        ! Z-coordinate of center of mass       
        REAL(KIND(0D0))  Brr        ! Dumb bell inertia                    
        REAL(KIND(0D0))  Px         ! X-component of spot weld position    
        REAL(KIND(0D0))  Py         ! Y-component of spot weld position    
        REAL(KIND(0D0))  Pz         ! Z-component of spot weld position    
        REAL(KIND(0D0))  Fmax       ! Fail constraint force, 0 => no-fail  
        REAL(KIND(0D0))  FORCE      ! Current constraint force b/ breaking 
        LOGICAL          FAILED     ! Flag for failed spot weld            
        INTEGER          EleID(3)   ! Element ID's                         
        INTEGER          Type(3)    ! 0/1 = triangle/quadrilateral         
        REAL(KIND(0D0))  Xi1(3)     ! Isoparametric element location       
        REAL(KIND(0D0))  Xi2(3)     ! Isoparametric element location       
      END TYPE                                                          
                                                                        
      TYPE (spot_weld_type),  DIMENSION(:), ALLOCATABLE :: SPOT_WELD    
                                                                        
      END MODULE spot_weld_
