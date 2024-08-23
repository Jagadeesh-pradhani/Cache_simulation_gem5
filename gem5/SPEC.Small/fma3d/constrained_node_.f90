!!
!! 60. REZONING: Constraints for midside nodes.
!!
      MODULE constrained_node_                                          
                                                                        
      TYPE :: constrained_node_type                                     
        INTEGER          ID              ! Internal ID/pointer to next location
        INTEGER          CNID            ! Constrained nodal point, ID     
        INTEGER          NPID(2)         ! Controlling nodal points, ID's  
        REAL(KIND(0D0))  TMass(2)        ! Effective mass at controlling nodes
        REAL(KIND(0D0))  RMass(2)        ! Effective mass at controlling nodes
        REAL(KIND(0D0))  Weight          ! Vcnid = V1*Weight + (1-Weight)*V2
      END TYPE                                                          
                                                                        
      TYPE (constrained_node_type), &
     &  DIMENSION(:), ALLOCATABLE :: CONSTRAINED_NODE                                                         
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (constrained_node_type), &
     &  DIMENSION(:), ALLOCATABLE :: RZ_CONSTRAINED_NODE                                                      
                                                                        
      END MODULE constrained_node_
