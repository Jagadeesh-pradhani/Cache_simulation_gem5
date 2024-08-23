!!
!! 12. SECTION PROPERTY DATA: Layered solid section structures.
!!
      MODULE layering_                                                  
                                                                        
      INTEGER   , PARAMETER :: MXNLY = 32 ! Maximum number of layers in structure
                                                                        
      TYPE :: layering_type                                             
        INTEGER          LupID        ! Lay-up ID (User defined value retained)
        INTEGER          Number_of_Layers                                  
                                      ! Number of layers in lay-up read.   
        INTEGER          Isys         ! Coordinate sytem, (0/n=global/special)
        INTEGER          LayID(MXNLY) ! Layer ID (User defined value retained)
        INTEGER          Ltype(MXNLY) ! Layer type (0/n=solid/membr. section ID)
        INTEGER          MatID(MXNLY) ! Material ID (must match materials avail.)
        REAL(KIND(0D0))  H(4,MXNLY)   ! Corner thicknesses, defines sub-elements
        REAL(KIND(0D0))  Ax(MXNLY)    ! Fiber orientation vector A, x-component
        REAL(KIND(0D0))  Ay(MXNLY)    ! Fiber orientation vector A, y-component
        REAL(KIND(0D0))  Az(MXNLY)    ! Fiber orientation vector A, z-component
        REAL(KIND(0D0))  Bx(MXNLY)    ! Fiber orientation vector B, x-component
        REAL(KIND(0D0))  By(MXNLY)    ! Fiber orientation vector B, y-component
        REAL(KIND(0D0))  Bz(MXNLY)    ! Fiber orientation vector B, z-component
      END TYPE                                                          
                                                                        
      TYPE (layering_type), DIMENSION(:), ALLOCATABLE :: LAYERING       
                                                                        
      END MODULE layering_
