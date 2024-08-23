!!
!! 28. NONREFLECTING BC
!!
      MODULE nonreflecting_bc_                                          
                                                                        
      TYPE :: nonreflecting_bc_type                                     
        INTEGER          NRID    ! Nonreflecting boundary condition ID     
        INTEGER          SetID   ! Segment set ID (-n/n=-SegID/SetID)      
        INTEGER          INRbgn  ! Starting location in NRBC_DATA          
        INTEGER          INRend  ! Ending location in NRBC_DATA            
      END TYPE                                                          
                                                                        
      TYPE (nonreflecting_bc_type), &
     &  DIMENSION(:), ALLOCATABLE :: NONREFLECTING_BC                                                         
                                                                        
      END MODULE nonreflecting_bc_
