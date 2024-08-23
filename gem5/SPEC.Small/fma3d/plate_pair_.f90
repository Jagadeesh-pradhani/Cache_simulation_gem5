!!
!! 59. REZONING: Plate pairs that share a common boundary.
!!
      MODULE plate_pair_                                                
                                                                        
      TYPE :: IBIDIS                                                    
       INTEGER           ID              ! ID of element                   
       INTEGER           IS              ! Side number of element          
      END TYPE                                                          
                                                                        
      TYPE :: plate_pair_type                                           
        INTEGER          :: ID           ! Internal ID/pointer to next location
        TYPE (IBIDIS)    :: IDS1                                           
        TYPE (IBIDIS)    :: IDS2                                           
        REAL(KIND(0D0))  :: CSPINI       ! Cosine of initial angle between plates
        REAL(KIND(0D0))  :: CosPhi       ! Cosine of angle between plates  
        REAL(KIND(0D0))  :: CosDot       ! Rate at which Cos(Phi) is changing
      END TYPE                                                          
                                                                        
      TYPE (plate_pair_type), DIMENSION(:), ALLOCATABLE :: PLATE_PAIR   
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (plate_pair_type), DIMENSION(:), ALLOCATABLE :: RZ_PLATE_PAIR
                                                                        
      END MODULE plate_pair_
