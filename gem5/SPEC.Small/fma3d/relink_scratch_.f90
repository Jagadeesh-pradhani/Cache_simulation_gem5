!!
!! 66. Relinking scratch arrays for node, element and segment ID sorting.
!!
      MODULE relink_scratch_                                            
                                                                        
      INTEGER, DIMENSION(:), ALLOCATABLE :: INPV, JNPV                  
      INTEGER, DIMENSION(:), ALLOCATABLE :: IELV, JELV                  
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISGV, JSGV                  
                                                                        
      END MODULE relink_scratch_
