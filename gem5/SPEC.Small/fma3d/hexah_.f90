!!
!! 42. HEXAHEDRON: Element parameter and results structures.
!!
      MODULE hexah_                                                     
                                                                        
      TYPE :: PAR_hexah                                                 
        INTEGER          EleID    ! Element ID   (User defined value retained)
        INTEGER          ParID    ! Part ID      (User defined value retained)
        INTEGER          MatID    ! Material ID  (Reset: points to MATERIAL)
        INTEGER          LupID    ! Lay-up ID    (Reset: points to LAYERING)
        INTEGER          IX(8)    ! NP Indicies  (Reset: internal node nums)
        INTEGER          Isv      ! Starting location for state variables  
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Volume   ! Initial volume                         
      END TYPE                                                          
                                                                        
      TYPE :: RES_hexah                                                 
        REAL(KIND(0D0))  Volume   ! Current volume                         
        REAL(KIND(0D0))  Int_Eng  ! Internal energy density                
        REAL(KIND(0D0))  Stress(6)! Stress                                 
        REAL(KIND(0D0))  Px(4)    ! Hourglass control forces, x-components 
        REAL(KIND(0D0))  Py(4)    ! Hourglass control forces, y-components 
        REAL(KIND(0D0))  Pz(4)    ! Hourglass control forces, z-components 
        REAL(KIND(0D0))  Xint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Yint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Zint(8)  ! Forces from stress divergence          
        REAL(KIND(0D0))  Time     ! Element local integration time         
        REAL(KIND(0D0))  DTelt    ! Element local critical time step       
        REAL(KIND(0D0))  DTnext   ! Element local integration time step    
        INTEGER          ISI      ! Subcycling index                       
      END TYPE                                                          
                                                                        
      TYPE :: hexah_type                                                
        TYPE (PAR_hexah) :: PAR                                         
        TYPE (RES_hexah) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (hexah_type), DIMENSION(:), ALLOCATABLE :: HEXAH             
!!
!! LAYERED SOLID HEXAH: Element parameter and results structures.
!! NOTE: This structure must match that used for HEXAH since these
!! LSOLD hexahedrons are "fed" to the HEXAH processing routines!!
!!
      TYPE (hexah_type), DIMENSION(:), ALLOCATABLE :: LSHEX             
                                                                        
      END MODULE hexah_
