!!
!! 37. SEGMENT SET: set structure.
!!
      MODULE segment_set_                                               
                                                                        
      TYPE :: segment_set_type                                          
        INTEGER           File_Number    ! On pass 1 the file no. is noted.
        INTEGER           Line_Number    ! On pass 1 the line no. is noted.
        INTEGER           SetID          ! Set ID                          
        INTEGER           Istart         ! Starting location               
        INTEGER           Iend           ! Ending location                 
        CHARACTER(4)      Flag           ! ALL = All nodes, elements or segments.
        CHARACTER(32)     Label          ! User's text ID/description      
      END TYPE                                                          
                                                                        
      TYPE (segment_set_type), DIMENSION(:), ALLOCATABLE :: SEGMENT_SET 
                                                                        
      END MODULE segment_set_
