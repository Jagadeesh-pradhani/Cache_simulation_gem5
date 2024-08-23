!!
!! 3. QA RECORD: Captured text lines from users input file.
!!
      MODULE  qa_record_                                                
                                                                        
      INTEGER   , PARAMETER :: MXNQA = 76 ! Number of Characters in QA record.
                                                                        
      TYPE :: qa_record_type                                            
        INTEGER               LSIZE                                     
        CHARACTER(MXNQA)      LINE     ! User QA record                 
      END TYPE                                                          
                                                                        
      TYPE (qa_record_type), DIMENSION(:), ALLOCATABLE :: QA_RECORD     
                                                                        
      END MODULE qa_record_
