!!
!! 13. SECTION PROPERTY DATA: Plate and membrane section structures.
!!
      MODULE section_2d_                                                
                                                                        
      TYPE :: section_2d_type                                           
        INTEGER          SecID    ! Section ID (User defined value retained)
        INTEGER          Ipts     ! Number of thickness integration points 
        INTEGER          Irule    ! Ftn: zeta coord's & weights f/ integration
        REAL(KIND(0D0))  Thickness! Thickness                              
        REAL(KIND(0D0))  RefLoc   ! Reference surface location (middle = 0.0)
        INTEGER          Isys     ! Coordinate sytem, (0/n=global/special) 
        REAL(KIND(0D0))  Ax       ! Fiber orientation vector A, x-component
        REAL(KIND(0D0))  Ay       ! Fiber orientation vector A, y-component
        REAL(KIND(0D0))  Az       ! Fiber orientation vector A, z-component
        REAL(KIND(0D0))  Bx       ! Fiber orientation vector B, x-component
        REAL(KIND(0D0))  By       ! Fiber orientation vector B, y-component
        REAL(KIND(0D0))  Bz       ! Fiber orientation vector B, z-component
      END TYPE                                                          
                                                                        
      TYPE (section_2d_type), DIMENSION(:), ALLOCATABLE :: SECTION_2D   
                                                                        
      END MODULE section_2d_
