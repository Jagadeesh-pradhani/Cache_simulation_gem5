!!
!! 57. SEGMENT: Segment parameters
!!
      MODULE segment_                                                   
                                                                        
      TYPE :: PAR_segment                                               
        INTEGER          SegID    ! Segment ID   (User defined value retained)
        INTEGER          ParID    ! Defining ID  (Reset: points to internal ID)
        INTEGER          IX(4)    ! NP Indicies  (Reset: internal node nums)
      END TYPE                                                          
                                                                        
      TYPE :: segment_type                                              
        TYPE (PAR_segment) :: PAR                                       
      END TYPE                                                          
                                                                        
      TYPE (segment_type), DIMENSION(:), ALLOCATABLE :: SEGMENT         
                                                                        
      END MODULE segment_
