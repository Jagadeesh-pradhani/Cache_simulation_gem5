!!
!! 29. NONREFLECTING BC DATA
!!
      MODULE nrbc_data_                                                 
                                                                        
      TYPE :: nrbc_data_type                                            
        INTEGER          SegID   ! Boundary segment ID                     
        INTEGER          Mel     ! Solid element w/w segment is associated 
        INTEGER          Type    ! Element type, 0/1/2=hexa/penta/tetra    
        REAL(KIND(0D0))  RCL     ! Density x longitudinal_wave_speed       
        REAL(KIND(0D0))  RCS     ! Density x transverse_shear_wave_speed   
      END TYPE                                                          
                                                                        
      TYPE (nrbc_data_type),  DIMENSION(:), ALLOCATABLE :: NRBC_DATA    
                                                                        
      END MODULE nrbc_data_
