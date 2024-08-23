!!
!! 34. TABULATED FUNCTION
!!
      MODULE tabulated_function_                                        
                                                                        
      INTEGER   , PARAMETER :: MXNTF = 1260      ! Maximum number of x,y-pairs.
                                                                        
      TYPE :: tabulated_function_type                                   
        INTEGER          TFID            ! Tabulated function ID           
        INTEGER          Number_of_Pairs ! Number of x,y-pairs read.       
        INTEGER          Last_Location   ! Interval used last access       
        REAL(KIND(0D0))  X(MXNTF)        ! Column of x-values              
        REAL(KIND(0D0))  Y(MXNTF)        ! Column of y-values              
        REAL(KIND(0D0))  SLOPE(MXNTF)    ! Slope dY/dX                     
      END TYPE                                                          
                                                                        
      TYPE (tabulated_function_type), &
     &  DIMENSION(:), ALLOCATABLE :: TABULATED_FUNCTION                                                     
                                                                        
      END MODULE tabulated_function_
