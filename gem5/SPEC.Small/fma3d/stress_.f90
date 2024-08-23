!!
!! 62. Concatenated shell element thru-the-thickness stress components.
!!
      MODULE stress_                                                    
                                                                        
      REAL(KIND(0D0)), DIMENSION(:,:), ALLOCATABLE :: STRESS                       
!!
!! Array used for subdividing the current mesh.
!!
      REAL(KIND(0D0)), DIMENSION(:,:), ALLOCATABLE :: RZ_STRESS                    
                                                                        
      END MODULE stress_
