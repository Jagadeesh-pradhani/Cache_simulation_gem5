!!
!! 41. NODAL POINT DATA/FORCE: Nodal point data struct
!!
      MODULE force_                                                     
                                                                        
      TYPE :: force_type                                                
        REAL(KIND(0D0))  Xext    ! X direction external force              
        REAL(KIND(0D0))  Yext    ! Y direction external force              
        REAL(KIND(0D0))  Zext    ! Z direction external force              
        REAL(KIND(0D0))  Xint    ! X direction internal force              
        REAL(KIND(0D0))  Yint    ! Y direction internal force              
        REAL(KIND(0D0))  Zint    ! Z direction internal force              
      END TYPE                                                          
                                                                        
      TYPE (force_type), DIMENSION(:), ALLOCATABLE :: FORCE             
!!
!! Array used for subdividing the current mesh.
!!
      TYPE (force_type), DIMENSION(:), ALLOCATABLE :: RZ_FORCE          
                                                                        
      END MODULE force_
