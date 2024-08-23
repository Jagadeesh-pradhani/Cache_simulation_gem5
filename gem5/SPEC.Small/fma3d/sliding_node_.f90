!!
!! 32. SLIDING INTERFACE
!!
      MODULE sliding_node_                                              
                                                                        
      TYPE :: sliding_node_type                                         
        INTEGER          Nsn     ! Sliding node ID, points to CONTACT_NODE(*)
        INTEGER          Mce     ! Contact el't ID, points to CONTACT_SURFACE(*)
        REAL(KIND(0D0))  P(4)    ! Current contact element nodal weighting 
        REAL(KIND(0D0))  Force   ! Force to place node on contact el't at n+1
        REAL(KIND(0D0))  Sign    ! Indicates contact came from + or - direction
        REAL(KIND(0D0))  Tcontact! Time at which contact occurred.         
        REAL(KIND(0D0))  Xcontact! Place in Mce where contact occurred.    
        REAL(KIND(0D0))  Ycontact! Place in Mce where contact occurred.    
        REAL(KIND(0D0))  Zcontact! Place in Mce where contact occurred.    
        REAL(KIND(0D0))  Xdepth  ! Vector from contact point to position.  
        REAL(KIND(0D0))  Ydepth  ! Vector from contact point to position.  
        REAL(KIND(0D0))  Zdepth  ! Vector from contact point to position.  
        INTEGER          INTEXT  ! Contact point wrt interior (1/0=int/ext)
      END TYPE                                                          
                                                                        
      TYPE (sliding_node_type), &
     &  DIMENSION(:), ALLOCATABLE :: SLIDING_NODE                                                                 
                                                                        
      TYPE (sliding_node_type) :: TRIANGLE_SLINODE                      
                                                                        
      END MODULE sliding_node_
