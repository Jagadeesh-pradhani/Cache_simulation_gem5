!!
!! 26. DAMPER BC: Nodal point axial/torsional viscous restraint.
!!
      MODULE damper_bc_                                                 
                                                                        
      TYPE :: RES_damper_bc                                             
        REAL(KIND(0D0))  Delta   ! Axial velocity/rotation rate            
        REAL(KIND(0D0))  Int_Eng ! Internal energy                         
        REAL(KIND(0D0))  Force   ! Axial force/torque                      
        REAL(KIND(0D0))  DTelt   ! Element local critical time step        
      END TYPE                                                          
                                                                        
      TYPE :: damper_bc_type                                            
        INTEGER          DprID   ! Damper ID (User defined value retained) 
        INTEGER          SetID   ! Node set ID (-n/0/n=-NPID/all/Set ID)   
        INTEGER          MatID   ! Material ID  (Reset: points to MATERIAL)
        INTEGER          Type    ! Axial or torsional (0/1=axial/torsional)
        INTEGER          Follow  ! Flag for follower Damper (0/1/2=U/V/Axis)
        INTEGER          Isv     ! Starting location for state variables   
        REAL(KIND(0D0))  Axis(3) ! Direction/hinge-axis orientation        
        TYPE (RES_damper_bc) :: RES                                     
      END TYPE                                                          
                                                                        
      TYPE (damper_bc_type), DIMENSION(:), ALLOCATABLE :: DAMPER_BC     
                                                                        
      END MODULE damper_bc_
