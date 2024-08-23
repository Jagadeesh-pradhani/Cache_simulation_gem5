!!
!! 7. 1-D STRAIN GAUGE DEFINITION: Gauge definition and results structures.
!!
      MODULE gauge1d_                                                   
                                                                        
      TYPE :: PAR_gauge1d                                               
        INTEGER          GauID   ! Gauge ID (User defined value retained)  
        INTEGER          TBflg   ! Truss/Beam Flag (0/1/2=no/truss/beam)   
        INTEGER          EleID   ! Element ID (Reset: points to TRUSS/BEAMS)
        INTEGER          GauLoc  ! Cross section location (fiber strain)   
        INTEGER          IX(4)   ! Nodal points defining gauge section     
        INTEGER          NUMIX   ! Actual number of nodes specified        
      END TYPE                                                          
                                                                        
      TYPE :: RES_gauge1d                                               
        REAL(KIND(0D0))  Strain  ! Logrithmic strain, LOGe(L/Lo)           
        REAL(KIND(0D0))  Torsion ! Logrithmic strain, LOGe(THETA/THETAo)   
      END TYPE                                                          
                                                                        
      TYPE :: gauge1d_type                                              
        TYPE (PAR_gauge1d) :: PAR                                       
        TYPE (RES_gauge1d) :: RES                                       
      END TYPE                                                          
                                                                        
      TYPE (gauge1d_type), DIMENSION(:), ALLOCATABLE :: GAUGE1D         
                                                                        
      END MODULE gauge1d_
