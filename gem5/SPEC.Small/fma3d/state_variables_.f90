!!
!! 63. Concatenated material state variable storage.
!!
      MODULE state_variables_                                           
                                                                        
      REAL(KIND(0D0)), DIMENSION(:), ALLOCATABLE :: STATE_VARIABLES                
!!
!! Array used for subdividing the current mesh.
!!
      REAL(KIND(0D0)), DIMENSION(:), ALLOCATABLE :: RZ_STATE_VARIABLES             
                                                                        
      END MODULE state_variables_
