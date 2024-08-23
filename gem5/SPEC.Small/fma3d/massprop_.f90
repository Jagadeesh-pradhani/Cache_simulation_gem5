!!
!! 6. RESULTS OUTPUT/MASS PROPERTY: Mass property calculation request structure.
!!
      MODULE massprop_                                                  
                                                                        
      TYPE :: massprop_type                                             
        INTEGER          MPID    ! Mass property request ID                
        INTEGER          SetID   ! Nodal point set ID                      
        REAL(KIND(0D0))  Mass    ! Total mass for node set                 
        REAL(KIND(0D0))  Xmv     ! X-momentum                              
        REAL(KIND(0D0))  Ymv     ! Y-momentum                              
        REAL(KIND(0D0))  Zmv     ! Z-momentum                              
        REAL(KIND(0D0))  Xcm     ! X-coordinate for center of mass         
        REAL(KIND(0D0))  Ycm     ! Y-coordinate for center of mass         
        REAL(KIND(0D0))  Zcm     ! Z-coordinate for center of mass         
        REAL(KIND(0D0))  Vxcm    ! X-velocity component for center of mass 
        REAL(KIND(0D0))  Vycm    ! Y-velocity component for center of mass 
        REAL(KIND(0D0))  Vzcm    ! Z-velocity component for center of mass 
        REAL(KIND(0D0))  KE      ! Translational kinetic energy            
        INTEGER          Irot    ! Rotation (0/1/2/3=no/wrt cm/wrt 0/wrt pt)
        REAL(KIND(0D0))  Xnert   ! X-coordinate for Irot = 3               
        REAL(KIND(0D0))  Ynert   ! Y-coordinate for Irot = 3               
        REAL(KIND(0D0))  Znert   ! Z-coordinate for Irot = 3               
        REAL(KIND(0D0))  Vxnert  ! X-velocity   for Irot = 3               
        REAL(KIND(0D0))  Vynert  ! Y-velocity   for Irot = 3               
        REAL(KIND(0D0))  Vznert  ! Z-velocity   for Irot = 3               
        REAL(KIND(0D0))  B(6)    ! Inertia, B(1:6)=(Ixx,Iyy,Izz,Ixy,Ixz,Iyz)
        REAL(KIND(0D0))  Oxmv    ! X-component of angular momemtum         
        REAL(KIND(0D0))  Oymv    ! Y-component of angular momemtum         
        REAL(KIND(0D0))  Ozmv    ! Z-component of angular momemtum         
        REAL(KIND(0D0))  Omega   ! Angular velocity rads/sec               
        REAL(KIND(0D0))  Ax      ! Axis of rotation, x-component           
        REAL(KIND(0D0))  Ay      ! Axis of rotation, y-component           
        REAL(KIND(0D0))  Az      ! Axis of rotation, z-component           
      END TYPE                                                          
                                                                        
      TYPE (massprop_type), DIMENSION(:), ALLOCATABLE :: MASSPROP       
                                                                        
      END MODULE massprop_
