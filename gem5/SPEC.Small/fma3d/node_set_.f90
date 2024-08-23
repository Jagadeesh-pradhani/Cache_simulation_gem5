!!
!! 35. NODE SET: set structure.
!!
      MODULE node_set_                                                  
                                                                        
      TYPE :: node_set_type                                             
        INTEGER           File_Number    ! On pass 1 the file no. is noted.
        INTEGER           Line_Number    ! On pass 1 the line no. is noted.
        INTEGER           SetID          ! Set ID                          
        INTEGER           Istart         ! Starting location               
        INTEGER           Iend           ! Ending location                 
        CHARACTER(4)      Flag           ! ALL = All nodes, elements or segments.
        CHARACTER(32)     Label          ! User's text ID/description      
      END TYPE                                                          
                                                                        
      TYPE (node_set_type), DIMENSION(:), ALLOCATABLE :: NODE_SET       
                                                                        
      END MODULE node_set_
