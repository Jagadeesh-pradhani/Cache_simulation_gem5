!!
!! 21. RIGID WALL BC
!!
      MODULE rigid_wall_bc_                                             
                                                                        
      TYPE :: rigid_wall_bc_type                                        
        INTEGER          RWID     ! Rigid wall ID                          
        INTEGER          SetID    ! Node set ID (-n/0/n=-NPID/all/Set ID)  
        REAL(KIND(0D0))  P1(3)    ! X,y,z of point on wall (center of mass)
        REAL(KIND(0D0))  P2(3)    ! X,y,z of point on wall (width direction)
        REAL(KIND(0D0))  P3(3)    ! X,y,z of point on wall (height direction)
        REAL(KIND(0D0))  Cn(3)    ! Wall normal vector, Cn = -(P2-P1)x(P3-P1)
        REAL(KIND(0D0))  Width    ! Width in the 12-direction (-W/2 to W/2)
                                  ! (Width = 0, infinite width)            
        REAL(KIND(0D0))  Height   ! Height in the 13-direction (-H/2 to H/2)
                                  ! (Height =0, infinite height)           
        REAL(KIND(0D0))  CoF      ! Coefficient of friction                
        INTEGER          Kode     ! Constraint code (0/1/2=rigid/inertial/moving)
        REAL(KIND(0D0))  Mass     ! Concentrated mass at P1                
        REAL(KIND(0D0))  B(3,3)   ! Concentrated inertia tensor at P1      
        REAL(KIND(0D0))  C(3,3)   ! Inverse of concentrated inertia tensor at P1
        REAL(KIND(0D0))  Disp(3)  ! Displacement of P1                     
        REAL(KIND(0D0))  Vel(3)   ! Velocity of P1                         
        REAL(KIND(0D0))  Accel(3) ! Acceleration of P1                     
        REAL(KIND(0D0))  Omega(3) ! Angular velocity of P1                 
        REAL(KIND(0D0))  Force(3) ! Forces acting on P1                    
        REAL(KIND(0D0))  Torque(3)! Torques acting on P1                   
        INTEGER          Code     ! Constraint code on P1 (see Displacement BC)
        INTEGER          HstID    ! History ID (tabulated function ID)     
        INTEGER          Kavd     ! Constrained kinematic variable (1/2/3=a/v/d)
        REAL(KIND(0D0))  Scale    ! History function scale factor          
        REAL(KIND(0D0))  Ax       ! BC direction, x-component (Code=4,40,  
        REAL(KIND(0D0))  Ay       ! BC direction, y-component  44,8,80,88) 
        REAL(KIND(0D0))  Az       ! BC direction, z-component              
      END TYPE                                                          
                                                                        
      TYPE (rigid_wall_bc_type), &
     &  DIMENSION(:), ALLOCATABLE :: RIGID_WALL_BC                                                               
                                                                        
      END MODULE rigid_wall_bc_
