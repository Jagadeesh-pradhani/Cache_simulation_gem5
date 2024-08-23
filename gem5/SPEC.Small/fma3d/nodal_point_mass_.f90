!!
!! 17. CONCENTRATED NODAL MASS: Parameter and results structures.
!!
      MODULE nodal_point_mass_                                          
                                                                        
      TYPE :: nodal_point_mass_type                                     
        INTEGER          NMID     ! Nodal mass ID (User defined value retained)
        INTEGER          NPID     ! Nodal point at which mass is concentrated
        INTEGER          Prop     ! Computed mass properties not_used/used (0/1)
        REAL(KIND(0D0))  Pzero(3) ! Initial X,Y,Z position of C.M.         
        REAL(KIND(0D0))  Disp(3)  ! Displacement of C.M.                   
        REAL(KIND(0D0))  Vel(3)   ! Velocity of C.M.                       
        REAL(KIND(0D0))  Accel(3) ! Acceleration of C.M.                   
        REAL(KIND(0D0))  Omega(3) ! Angular velocity of C.M.               
        REAL(KIND(0D0))  Theta(3) ! Angular velocity integral              
        REAL(KIND(0D0))  Mass     ! Concentrated mass                      
        REAL(KIND(0D0))  B(3,3)   ! Concentrated inertia                   
        REAL(KIND(0D0))  Force(3) ! Forces acting on C.M.                  
        REAL(KIND(0D0))  Torque(3)! Torques acting on C.M.                 
      END TYPE                                                          
                                                                        
      TYPE (nodal_point_mass_type), &
     &  DIMENSION(:), ALLOCATABLE :: NODAL_POINT_MASS                                                         
                                                                        
      END MODULE nodal_point_mass_
