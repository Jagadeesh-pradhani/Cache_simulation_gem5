!!
!! 2. INPUT RECORD TOKENS: Character, real, and integer "stacking."
!!
      MODULE  value_                                                    
                                                                        
      TYPE :: value_type                                                
        INTEGER                  LOC     ! Location in input record line   
        CHARACTER(4)             VTYP    ! Type of value found, 'N/C/I/R/D/U'
        CHARACTER(32)            CVAL    ! Character value                 
        INTEGER                  IVAL    ! Integer value (I*4)             
        REAL(KIND(0D0))          RVAL    ! Real value (R*8)                
        REAL(KIND(0D0))          DVAL    ! Real value (R*8)                
      END TYPE                                                          
                                                                        
      TYPE (value_type), DIMENSION(:), ALLOCATABLE :: VALUE             
                                                                        
      END MODULE value_
