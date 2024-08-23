!!
!! 40. NODAL POINT DATA/MOTION: Nodal point data structures.
!!
      MODULE motion_                                                    
                                                                        
      TYPE :: motion_type                                               
        REAL(KIND(0D0))  Px      ! Initial x-position                      
        REAL(KIND(0D0))  Py      ! Initial y-position                      
        REAL(KIND(0D0))  Pz      ! Initial z-position                      
        REAL(KIND(0D0))  Ux      ! X displacement                          
        REAL(KIND(0D0))  Uy      ! Y displacement                          
        REAL(KIND(0D0))  Uz      ! Z displacement                          
        REAL(KIND(0D0))  Vx      ! X velocity                              
        REAL(KIND(0D0))  Vy      ! Y velocity                              
        REAL(KIND(0D0))  Vz      ! Z velocity                              
        REAL(KIND(0D0))  Ax      ! X acceleration                          
        REAL(KIND(0D0))  Ay      ! Y acceleration                          
        REAL(KIND(0D0))  Az      ! Z acceleration                          
      END TYPE                                                          
                                                                        
      TYPE (motion_type), DIMENSION(:), ALLOCATABLE :: MOTION           
                                                                        
      TYPE (motion_type), DIMENSION(8)              :: EMOTION          
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (motion_type), DIMENSION(:), ALLOCATABLE :: RZ_MOTION        
                                                                        
      END MODULE motion_
