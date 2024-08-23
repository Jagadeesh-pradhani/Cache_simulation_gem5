!!
!! 39. NODAL POINT DATA/COORDINATES: Nodal point data structures.
!!
      MODULE node_                                                      
                                                                        
      TYPE :: node_type                                                 
        INTEGER          ID      ! Nodal Point ID (User supplied value)    
        INTEGER          ISI     ! Subcycling index                        
        INTEGER          IRB     ! Linked list pointer for rigid bodies    
        INTEGER          IRT     ! Pointer for rotational DOF's            
        REAL(KIND(0D0))  Time    ! Nodal point time                        
        REAL(KIND(0D0))  DTlast  ! Nodal point time step                   
        REAL(KIND(0D0))  DTnext  ! Nodal point time step                   
        REAL(KIND(0D0))  Mass    ! Nodal point mass                        
        REAL(KIND(0D0))  Minv    ! Inverse nodal point mass                
        LOGICAL          On      ! Used to expidite nodal loop if-tests    
        INTEGER          ICF     ! Contact flag, 0/1/2=none/as-a-np/as-an-el
      END TYPE                                                          
                                                                        
      TYPE (node_type), DIMENSION(:), ALLOCATABLE :: NODE               
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (node_type), DIMENSION(:), ALLOCATABLE :: RZ_NODE            
                                                                        
      END MODULE node_
