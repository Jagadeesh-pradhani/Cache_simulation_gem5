!!
!! 30. SLIDING INTERFACE
!!
      MODULE sliding_interface_                                         
                                                                        
      TYPE :: sliding_interface_type                                    
        INTEGER          SIID     ! Sliding interface ID                   
        INTEGER          Typ1     ! Type (0/1=segment-set/node-set)        
        INTEGER          S1ID     ! Side 1 ID                              
        INTEGER          Typ2     ! Type (0/1=segment-set/node-set)        
        INTEGER          S2ID     ! Side 2 ID (0=single-surface interface) 
        INTEGER          Type     ! Type (n=0,1,2,3,...) [Type 0 only 7/21/90]
        INTEGER          Isym     ! Flag, 0/1/2=symmetric/master-slave/vice versa
        REAL(KIND(0D0))  CoF      ! Coefficient of friction                
        REAL(KIND(0D0))  Begin    ! Begin calculations                     
        REAL(KIND(0D0))  End      ! End calculations                       
        REAL(KIND(0D0))  Factor   ! Impact force relaxation factor         
        REAL(KIND(0D0))  Capture  ! Capture distance for penetration test  
        REAL(KIND(0D0))  Border   ! Contact element "exterior" border      
        INTEGER          Sort_Freq! Sorting frequency for initiating contact
        INTEGER          ISN1bgn  ! Starting location in SLIDING_NODE, side 1
        INTEGER          ISN1end  ! Ending   location in SLIDING_NODE, side 1
        INTEGER          ISN2bgn  ! Starting location in SLIDING_NODE, side 2
        INTEGER          ISN2end  ! Ending   location in SLIDING_NODE, side 2
        INTEGER          ICE1bgn  ! Starting location in CONTACT_SURFACE, side 1
        INTEGER          ICE1end  ! Ending   location in CONTACT_SURFACE, side 1
        INTEGER          ICE2bgn  ! Starting location in CONTACT_SURFACE, side 2
        INTEGER          ICE2end  ! Ending   location in CONTACT_SURFACE, side 2
      END TYPE                                                          
                                                                        
      TYPE (sliding_interface_type), &
     & DIMENSION(:), ALLOCATABLE :: SLIDING_INTERFACE                                                       
                                                                        
      END MODULE sliding_interface_
