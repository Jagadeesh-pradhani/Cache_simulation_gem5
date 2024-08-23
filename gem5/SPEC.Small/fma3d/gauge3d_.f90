!!
!! 9. 3-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      MODULE gauge3d_                                                   
                                                                        
      TYPE :: PAR_gauge3d                                               
        INTEGER          GauID     ! Gauge ID (User defined value retained)
        INTEGER          HPTflg    ! Hexa/Penta/Tetra flag (0/1/2/3=no/H/P/T)
        INTEGER          EleID     ! Element ID (Reset: pts. to HEXA,PENTA,TETRA)
        INTEGER          IX(8)     ! Nodal points defining strain.         
        INTEGER          NUMIX     ! Actual number of nodes specified      
      END TYPE                                                          
                                                                        
      TYPE :: RES_gauge3d                                               
        REAL(KIND(0D0))  Strain(6) ! Logrithmic strain, LOGe(L/Lo)         
      END TYPE                                                          
                                                                        
      TYPE :: gauge3d_type                                              
        TYPE (PAR_gauge3d) :: PAR                                       
        TYPE (RES_gauge3d) :: RES                                       
      END TYPE                                                          
                                                                        
      TYPE (gauge3d_type), DIMENSION(:), ALLOCATABLE :: GAUGE3D         
                                                                        
      END MODULE gauge3d_
