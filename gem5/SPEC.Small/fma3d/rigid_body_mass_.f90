!!
!! 16. RIGID BODY MASS: Parameter and results structures.
!!
      MODULE rigid_body_mass_                                           
                                                                        
      TYPE :: rigid_body_mass_type                                      
        INTEGER          RMID     ! Rigid mass ID (User defined value retained)
        INTEGER          NPID     ! Nodal point at which mass is concentrated
        REAL(KIND(0D0))  Pzero(3) ! Initial X,Y,Z position of C.M.         
        REAL(KIND(0D0))  Disp(3)  ! Displacement of C.M.                   
        REAL(KIND(0D0))  Vel(3)   ! Velocity of C.M.                       
        REAL(KIND(0D0))  Accel(3) ! Acceleration of C.M.                   
        REAL(KIND(0D0))  Omega(3) ! Angular velocity of C.M.               
        REAL(KIND(0D0))  Mass     ! Concentrated mass                      
        REAL(KIND(0D0))  B(3,3)   ! Concentrated inertia                   
        REAL(KIND(0D0))  Force(3) ! Forces acting on C.M.                  
        REAL(KIND(0D0))  Torque(3)! Torques acting on C.M.                 
      END TYPE                                                          
                                                                        
      TYPE (rigid_body_mass_type), &
     &  DIMENSION(:), ALLOCATABLE :: RIGID_BODY_MASS                                                           
                                                                        
      END MODULE rigid_body_mass_
