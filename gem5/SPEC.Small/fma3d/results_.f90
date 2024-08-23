!!
!! 4. RESULTS OUTPUT (PLOTTING_DATABASE): Results file control structure.
!!
      MODULE  results_                                                  
                                                                        
      INTEGER   , PARAMETER :: MXNRS = 20 ! Number of named entries in structure
                                                                        
      TYPE :: results_type                                              
        CHARACTER(12)    Title      ! User's results file title            
        INTEGER          ResID      ! Results file ID                      
        INTEGER          OpfID      ! Output file ID                       
        INTEGER          Number_of_Entries                                 
        REAL(KIND(0D0))  Begin      ! Begin writing to printer file        
        REAL(KIND(0D0))  End        ! End writing to printer file          
        REAL(KIND(0D0))  Delta      ! Time increment between writes        
        REAL(KIND(0D0))  Time       ! Next time to write                   
        INTEGER          MASSP      ! Flag (0/n=no/yes)                    
        INTEGER          GAUGE      ! Flag (0/n=no/yes)                    
        INTEGER          NODES      ! Flag (-1/n=ALL/Nodal Point SetID)    
        INTEGER          HEXAH      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          PENTA      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          TETRA      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          LSOLD      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          MEMBT      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          MEMBQ      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          TRUSS      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          PLATT      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          PLATQ      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          BEAMS      ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          SPRING     ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          DAMPER     ! Flag (-1/n=ALL/Element SetID)        
        INTEGER          GLOBAL     ! Flag (0/n=no/yes)                    
      END TYPE                                                          
                                                                        
      TYPE (results_type), DIMENSION(:), ALLOCATABLE :: RESULTS         
                                                                        
      END MODULE results_
