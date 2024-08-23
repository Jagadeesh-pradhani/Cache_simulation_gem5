!!
!! 25. SPRING BC: Nodal point axial/torsional spring restraint.
!!
      MODULE spring_bc_                                                 
                                                                        
      TYPE :: RES_spring_bc                                             
        REAL(KIND(0D0))  Delta   ! Axial displacement/rotation             
        REAL(KIND(0D0))  Int_Eng ! Internal energy                         
        REAL(KIND(0D0))  Force   ! Axial force/torque                      
        REAL(KIND(0D0))  DTelt   ! Element local critical time step        
      END TYPE                                                          
                                                                        
      TYPE :: spring_bc_type                                            
        INTEGER          SprID   ! Spring ID (User defined value retained) 
        INTEGER          SetID   ! Node set ID (-n/0/n=-NPID/all/Set ID)   
        INTEGER          MatID   ! Material ID  (Reset: points to MATERIAL)
        INTEGER          Type    ! Axial or torsional (0/1=axial/torsional)
        INTEGER          Follow  ! Flag for follower spring (0/1=U/Axis)   
        INTEGER          Isv     ! Starting location for state variables   
        REAL(KIND(0D0))  Axis(3) ! Direction/hinge-axis orientation        
        TYPE (RES_spring_bc) :: RES                                     
      END TYPE                                                          
                                                                        
      TYPE (spring_bc_type), DIMENSION(:), ALLOCATABLE :: SPRING_BC     
                                                                        
      END MODULE spring_bc_
