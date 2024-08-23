!!
!! 11. MATERIAL PROPERTY DATA NAMES: Used with input/output operations.
!!
      MODULE property_                                                  
                                                                        
      INTEGER   , PARAMETER :: NPNMV = 110       ! Number of variable names
                                                                        
      TYPE :: property_type                                             
        INTEGER                         :: Number_of_Entries ! = NPNMV  
        INTEGER,       DIMENSION(NPNMV) :: Set                          
        CHARACTER(12), DIMENSION(NPNMV) :: Name                         
        INTEGER,       DIMENSION(NPNMV) :: Location                     
      END TYPE                                                          
                                                                        
      TYPE (property_type) :: PROPERTY                                  
                                                                        
      END MODULE property_
