!!
!! 15. RIGID BODY: Parameter and results structures.
!!
      MODULE rigid_body_                                                
                                                                        
      TYPE :: rigid_body_type                                           
        INTEGER          RBID     ! Rigid body ID (User defined value retained)
        INTEGER          ParID    ! Part ID, Defines rigid body domain     
        INTEGER          Prop     ! Computed mass properties not_used/used (0/1)
        INTEGER          CMID     ! Pts. to user spec'd added/subst'd masses
        INTEGER          FirstNP  ! Pts. to 1st node in NODE(*).IRB defining RB
        INTEGER          Code(2)  ! Constraint code on CM (see Displacement BC)
        INTEGER          HstID(2) ! History ID (tabulated function ID)     
        INTEGER          Kavd(2)  ! Constrained kinematic variable (1/2/3=a/v/d)
        REAL(KIND(0D0))  Scale(2) ! History function scale factor          
        REAL(KIND(0D0))  Bx(2)    ! BC direction, x-component (Code=4,40,  
        REAL(KIND(0D0))  By(2)    ! BC direction, y-component  44,8,80,88) 
        REAL(KIND(0D0))  Bz(2)    ! BC direction, z-component              
        REAL(KIND(0D0))  Px       ! Initial X position of C.M.             
        REAL(KIND(0D0))  Py       ! Initial Y position of C.M.             
        REAL(KIND(0D0))  Pz       ! Initial Z position of C.M.             
        REAL(KIND(0D0))  Ux       ! X displacement of C.M.                 
        REAL(KIND(0D0))  Uy       ! Y displacement of C.M.                 
        REAL(KIND(0D0))  Uz       ! Z displacement of C.M.                 
        REAL(KIND(0D0))  Vx       ! X velocity of C.M.                     
        REAL(KIND(0D0))  Vy       ! Y velocity of C.M.                     
        REAL(KIND(0D0))  Vz       ! Z velocity of C.M.                     
        REAL(KIND(0D0))  Ax       ! X acceleration of C.M.                 
        REAL(KIND(0D0))  Ay       ! Y acceleration of C.M.                 
        REAL(KIND(0D0))  Az       ! Z acceleration of C.M.                 
        REAL(KIND(0D0))  Ox       ! X angular velocity of C.M.             
        REAL(KIND(0D0))  Oy       ! Y angular velocity of C.M.             
        REAL(KIND(0D0))  Oz       ! Z angular velocity of C.M.             
        REAL(KIND(0D0))  Mass     ! Concentrated mass                      
        REAL(KIND(0D0))  B(3,3)   ! Concentrated inertia                   
        REAL(KIND(0D0))  Force(3) ! Forces acting on C.M.                  
        REAL(KIND(0D0))  Torque(3)! Torques acting on C.M.                 
      END TYPE                                                          
                                                                        
      TYPE (rigid_body_type), DIMENSION(:), ALLOCATABLE :: RIGID_BODY   
                                                                        
      END MODULE rigid_body_
