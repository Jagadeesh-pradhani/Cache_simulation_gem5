!!
!! 8. 2-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      MODULE gauge2d_                                                   
                                                                        
      TYPE :: PAR_gauge2d                                               
        INTEGER          GauID     ! Gauge ID (User defined value retained)
        INTEGER          MPQTflg   ! Membrane/Plate/Quad/Triangle flag     
        INTEGER          EleID     ! Element ID (Reset: points to MEMB*,PLAT*)
        INTEGER          GauLoc    ! Plate surface (-1,0,+1; fiber strain) 
        INTEGER          IX(8)     ! Nodal DOF's defining strain gauge.    
        INTEGER          NUMIX     ! Actual number of nodes specified      
        REAL(KIND(0D0))  Thickness ! From element section data for PLATT/Q 
        REAL(KIND(0D0))  Refloc    ! From element section data for PLATT/Q 
      END TYPE                                                          
                                                                        
      TYPE :: RES_gauge2d                                               
        REAL(KIND(0D0))  Beta      ! Angle, 2-D defor. gradient polar decomp
        REAL(KIND(0D0))  Strain(3) ! Logrithmic strain, LOGe(L/Lo)         
      END TYPE                                                          
                                                                        
      TYPE :: gauge2d_type                                              
        TYPE (PAR_gauge2d) :: PAR                                       
        TYPE (RES_gauge2d) :: RES                                       
      END TYPE                                                          
                                                                        
      TYPE (gauge2d_type), DIMENSION(:), ALLOCATABLE :: GAUGE2D         
                                                                        
      END MODULE gauge2d_
