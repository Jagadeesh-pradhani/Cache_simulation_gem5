!!
!! 38. LOCATION: Location pointers for building the lists containing set
!! entity ID's.
!!
      MODULE location_                                                  
                                                                        
      TYPE :: location_type                                             
        INTEGER    :: Next_Node     ! = 1                                  
        INTEGER    :: Next_Element  ! = 1                                  
        INTEGER    :: Next_Segment  ! = 1                                  
      END TYPE                                                          
                                                                        
      TYPE (location_type) :: LOCATION                                  
                                                                        
      END MODULE location_
