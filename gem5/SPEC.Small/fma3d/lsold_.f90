!!
!! 45. LAYERED SOLID: element parameter and results structures.
!!
      MODULE lsold_                                                     
                                                                        
      INTEGER   , PARAMETER :: MXLSL = 32 ! Maximum number of layers in structure
                                                                        
      TYPE :: PAR_lsold                                                 
        INTEGER          EleID     ! Element ID   (User defined value retained)
        INTEGER          ParID     ! Part ID      (User defined value retained)
        INTEGER          LupID     ! Lay-up ID    (Reset: points to LAYERING)
        INTEGER          IX(8)     ! NP Indicies  (Reset: internal node nums)
        INTEGER          ID(MXLSL) ! Element ID's used for lay-up calcs    
        INTEGER          IGR      ! Parallel-safe accumulation group number
        REAL(KIND(0D0))  Volume    ! Initial volume                        
      END TYPE                                                          
                                                                        
      TYPE :: RES_lsold                                                 
        REAL(KIND(0D0))  Volume    ! Current volume                        
        REAL(KIND(0D0))  Int_Eng   ! Internal energy density               
        REAL(KIND(0D0))  H(4,MXLSL)! Layer thicknesses at each corner      
        REAL(KIND(0D0))  Xint(8)   ! Forces from stress divergence         
        REAL(KIND(0D0))  Yint(8)   ! Forces from stress divergence         
        REAL(KIND(0D0))  Zint(8)   ! Forces from stress divergence         
        REAL(KIND(0D0))  Time      ! Element local integration time        
        REAL(KIND(0D0))  DTelt     ! Element local critical time step      
        REAL(KIND(0D0))  DTnext    ! Element local integration time step   
        INTEGER          ISI       ! Subcycling index                      
      END TYPE                                                          
                                                                        
      TYPE :: lsold_type                                                
        TYPE (PAR_lsold) :: PAR                                         
        TYPE (RES_lsold) :: RES                                         
      END TYPE                                                          
                                                                        
      TYPE (lsold_type), DIMENSION(:), ALLOCATABLE :: LSOLD             
                                                                        
      END MODULE lsold_
