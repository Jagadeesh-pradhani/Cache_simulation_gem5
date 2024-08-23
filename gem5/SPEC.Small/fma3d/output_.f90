!!
!! 5. OUTPUT DATA NAMES: Used with input/output operations.
!!
      MODULE output_                                                    
                                                                        
      INTEGER   , PARAMETER :: MXNRN = 20 ! Number of named entries in structure
                                                                        
      TYPE :: output_type                                               
        INTEGER    :: Number_of_Entries   !   = MXNRN                      
        CHARACTER(6), DIMENSION(MXNRN) :: NAME                          
      END TYPE                                                          
                                                                        
      TYPE (output_type) :: OUTPUT                                      
                                                                        
      END MODULE output_
